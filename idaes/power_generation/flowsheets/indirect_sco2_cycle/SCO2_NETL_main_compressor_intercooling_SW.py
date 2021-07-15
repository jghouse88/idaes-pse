"""
SCO2 baslien cycle from the NETL report

--case baseline620 - Turbine inlet temperature 893.15 K (620 C).
--case basleine760 - Turbine inlet temperature 1033.15 K (760 C).

"""

import argparse
from pyomo.environ import (ConcreteModel,
                           Block,
                           Var,
                           Param,
                           Constraint,
                           SolverFactory,
                           TransformationFactory, TerminationCondition,
                           value, Expression, minimize, units)
from pyomo.network import Arc, SequentialDecomposition

# Import IDAES libraries
from idaes.core import FlowsheetBlock, UnitModelBlockData
from idaes.generic_models.unit_models import (Mixer, MomentumMixingType,
                                              PressureChanger, Heater,
                                              Separator, HeatExchanger)
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
# from idaes.unit_models.heat_exchanger import (delta_temperature_lmtd_callback,
#                                               HeatExchangerFlowPattern)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
import sys
sys.path.append('../using_surrogate_properties')
from SCO2_properties_surrogate import SCO2ParameterBlock
from idaes.generic_models.properties.swco2 import SWCO2ParameterBlock, StateVars, htpx

import idaes.logger as idaeslog

parser = argparse.ArgumentParser(description='Select case to simulate')
parser.add_argument("--case",
                    choices=["baseline620", "baseline760"],
                    default="baseline620",
                    help="Case to simulate. baseline620 is TIT 620 C and "
                    "basleine760 is TIT 760 C.")
parser.add_argument("--deltaTmin",
                    action="store_true",
                    help="Flag to activate deltaTmin constraint for "
                    "recuperators.")
args = parser.parse_args()

if args.case == "baseline620":
    flow_mol = 112.29e3  # mol/s
    boiler_inlet_temperature = 671.15  # K
    turbine_inlet_temperature = 893.15  # K
    split_frac_bypass_comp = 0.1625
    split_frac_fg_cooler = 0.0604
else:
    flow_mol = 89.116e3  # mol/s
    boiler_inlet_temperature = 786.15  # K
    turbine_inlet_temperature = 1033.15  # K
    split_frac_bypass_comp = 0.1566
    split_frac_fg_cooler = 0.0669

# Approach tempertaure for recuperators
deltaTmin_HTR = 5.6  # (or 10 F according to report)
deltaTmin_LTR = 5.6  # (or 10 F according to report)

# Setup solver and options
solver = SolverFactory('ipopt')
outlvl = idaeslog.INFO
tee = True

# Set up concrete model
m = ConcreteModel()

# Create a flowsheet block
m.fs = FlowsheetBlock(default={'dynamic': False})

# Create the properties param block
# m.fs.properties = SCO2ParameterBlock()
m.fs.properties = SWCO2ParameterBlock()

# Add unit models to the flowsheet
m.fs.boiler = Heater(default={'dynamic': False,
                              'property_package': m.fs.properties,
                              'has_pressure_change': True})

m.fs.turbine = PressureChanger(
    default={'dynamic': False,
             'property_package': m.fs.properties,
             'compressor': False,
             'thermodynamic_assumption': ThermodynamicAssumption.isentropic})

m.fs.HTR_pseudo_shell = Heater(default={'dynamic': False,
                                        'property_package': m.fs.properties,
                                        'has_pressure_change': True})

m.fs.HTR_pseudo_tube = Heater(default={'dynamic': False,
                                       'property_package': m.fs.properties,
                                       'has_pressure_change': True})

m.fs.LTR_pseudo_shell = Heater(default={'dynamic': False,
                                        'property_package': m.fs.properties,
                                        'has_pressure_change': True})

m.fs.LTR_pseudo_tube = Heater(default={'dynamic': False,
                                       'property_package': m.fs.properties,
                                       'has_pressure_change': True})

m.fs.splitter_1 = Separator(default={"property_package": m.fs.properties,
                                     "outlet_list": ["bypass", "to_cooler"]})

m.fs.co2_cooler = Heater(default={'dynamic': False,
                                  'property_package': m.fs.properties,
                                  'has_pressure_change': True})

m.fs.main_compressor_stage_1 = PressureChanger(
    default={'dynamic': False,
             'property_package': m.fs.properties,
             'compressor': True,
             'thermodynamic_assumption': ThermodynamicAssumption.isentropic})

m.fs.inter_cooler = Heater(default={'dynamic': False,
                                    'property_package': m.fs.properties,
                                    'has_pressure_change': True})

m.fs.main_compressor_stage_2 = PressureChanger(
    default={'dynamic': False,
             'property_package': m.fs.properties,
             'compressor': True,
             'thermodynamic_assumption': ThermodynamicAssumption.isentropic})

m.fs.bypass_compressor = PressureChanger(
    default={'dynamic': False,
             'property_package': m.fs.properties,
             'compressor': True,
             'thermodynamic_assumption': ThermodynamicAssumption.isentropic})

m.fs.splitter_2 = Separator(default={"property_package": m.fs.properties,
                                     "ideal_separation": False,
                                     "outlet_list": ["to_FG_cooler",
                                                     "to_LTR"]})

m.fs.FG_cooler = Heater(default={'dynamic': False,
                                 'property_package': m.fs.properties,
                                 'has_pressure_change': True})

m.fs.mixer = Mixer(default={"property_package": m.fs.properties,
                            # "momentum_mixing_type":
                            # MomentumMixingType.equality,
                            "inlet_list": ["FG_out", "LTR_out", "bypass"]})

# Not using for now; Instead using the common approach of two heater blocks
# and connect the heat streams.

# m.fs.LTR = HeatExchanger(
#     default={'dynamic': False,
#              'delta_temperature_callback': delta_temperature_lmtd_callback,
#              'flow_pattern': HeatExchangerFlowPattern.countercurrent,
#              'shell': {'property_package': m.fs.properties},
#              'tube': {'property_package': m.fs.properties}})

# # Connect the flowsheet
m.fs.s01 = Arc(source=m.fs.boiler.outlet,
               destination=m.fs.turbine.inlet)
m.fs.s02 = Arc(source=m.fs.turbine.outlet,
               destination=m.fs.HTR_pseudo_shell.inlet)
m.fs.s03 = Arc(source=m.fs.HTR_pseudo_shell.outlet,
               destination=m.fs.LTR_pseudo_shell.inlet)
m.fs.s04 = Arc(source=m.fs.LTR_pseudo_shell.outlet,
               destination=m.fs.splitter_1.inlet)
m.fs.s05 = Arc(source=m.fs.splitter_1.to_cooler,
               destination=m.fs.co2_cooler.inlet)
m.fs.s06 = Arc(source=m.fs.splitter_1.bypass,
               destination=m.fs.bypass_compressor.inlet)
m.fs.s07 = Arc(source=m.fs.co2_cooler.outlet,
               destination=m.fs.main_compressor_stage_1.inlet)
m.fs.s08 = Arc(source=m.fs.main_compressor_stage_1.outlet,
               destination=m.fs.inter_cooler.inlet)
m.fs.s09 = Arc(source=m.fs.inter_cooler.outlet,
               destination=m.fs.main_compressor_stage_2.inlet)
m.fs.s10 = Arc(source=m.fs.bypass_compressor.outlet,
               destination=m.fs.mixer.bypass)
m.fs.s11 = Arc(source=m.fs.main_compressor_stage_2.outlet,
               destination=m.fs.splitter_2.inlet)
m.fs.s12 = Arc(source=m.fs.splitter_2.to_FG_cooler,
               destination=m.fs.FG_cooler.inlet)
m.fs.s13 = Arc(source=m.fs.splitter_2.to_LTR,
               destination=m.fs.LTR_pseudo_tube.inlet)
m.fs.s14 = Arc(source=m.fs.LTR_pseudo_tube.outlet,
               destination=m.fs.mixer.LTR_out)
m.fs.s15 = Arc(source=m.fs.FG_cooler.outlet,
               destination=m.fs.mixer.FG_out)
m.fs.s16 = Arc(source=m.fs.mixer.outlet,
               destination=m.fs.HTR_pseudo_tube.inlet)

m.fs.boiler.inlet.flow_mol.fix(flow_mol)
# m.fs.boiler.inlet.temperature.fix(685.15)
m.fs.boiler.inlet.pressure.fix(34.51e6)
m.fs.boiler.inlet.enth_mol.fix(htpx(T=boiler_inlet_temperature*units.K,
                                    P=m.fs.boiler.inlet.pressure[0].value*units.Pa))

# m.fs.boiler.outlet.temperature.fix(893.15)  # Turbine inlet T = 620 C
m.fs.boiler.deltaP.fix(-0.21e6)
m.fs.boiler.outlet.enth_mol.fix(htpx(
    T=turbine_inlet_temperature*units.K,
    P=value(m.fs.boiler.inlet.pressure[0]+m.fs.boiler.deltaP[0])*units.Pa))

m.fs.boiler.initialize(outlvl=outlvl)

propagate_state(m.fs.s01)

m.fs.turbine.ratioP.fix(1/3.68)
m.fs.turbine.efficiency_isentropic.fix(0.927)
m.fs.turbine.initialize(outlvl=outlvl)

propagate_state(m.fs.s02)

# m.fs.HTR_pseudo_shell.outlet.temperature.fix(489.15)
m.fs.HTR_pseudo_shell.deltaP.fix(-0.07e6)
m.fs.HTR_pseudo_shell.outlet.enth_mol.fix(htpx(
    T=446.15*units.K,
    P=value(m.fs.HTR_pseudo_shell.inlet.pressure[0] +
            m.fs.HTR_pseudo_shell.deltaP[0])*units.Pa))

m.fs.HTR_pseudo_shell.initialize(outlvl=outlvl)

propagate_state(m.fs.s03)

# m.fs.LTR_pseudo_shell.outlet.temperature.fix(354.15)
m.fs.LTR_pseudo_shell.deltaP.fix(-0.07e6)
m.fs.LTR_pseudo_shell.outlet.enth_mol.fix(htpx(
    T=330.15*units.K,
    P=value(m.fs.LTR_pseudo_shell.inlet.pressure[0] +
            m.fs.LTR_pseudo_shell.deltaP[0])*units.Pa))
m.fs.LTR_pseudo_shell.initialize(outlvl=outlvl)

propagate_state(m.fs.s04)
m.fs.splitter_1.split_fraction[0, "bypass"].fix(split_frac_bypass_comp)
m.fs.splitter_1.initialize(outlvl=outlvl)

propagate_state(m.fs.s05)
# m.fs.co2_cooler.outlet.temperature.fix(308.15)
m.fs.co2_cooler.deltaP.fix(-0.07e6)
m.fs.co2_cooler.outlet.enth_mol.fix(htpx(
    T=308.15*units.K,
    P=value(m.fs.co2_cooler.inlet.pressure[0]+m.fs.co2_cooler.deltaP[0])
    * units.Pa))
m.fs.co2_cooler.initialize(outlvl=outlvl)

propagate_state(m.fs.s06)
m.fs.bypass_compressor.efficiency_isentropic.fix(0.85)
m.fs.bypass_compressor.ratioP.fix(3.75)
m.fs.bypass_compressor.initialize(outlvl=outlvl)

propagate_state(m.fs.s07)
m.fs.main_compressor_stage_1.efficiency_isentropic.fix(0.85)
m.fs.main_compressor_stage_1.ratioP.fix(1.95)
m.fs.main_compressor_stage_1.initialize(outlvl=outlvl)

propagate_state(m.fs.s08)
m.fs.inter_cooler.deltaP.fix(-0.07e6)
m.fs.inter_cooler.outlet.enth_mol.fix(htpx(
    T=308.15*units.K,
    P=value(m.fs.inter_cooler.inlet.pressure[0]+m.fs.inter_cooler.deltaP[0])
    * units.Pa))
m.fs.inter_cooler.initialize(outlvl=outlvl)

propagate_state(m.fs.s09)
m.fs.main_compressor_stage_2.efficiency_isentropic.fix(0.85)
m.fs.main_compressor_stage_2.ratioP.fix(1.95)
m.fs.main_compressor_stage_2.initialize(outlvl=idaeslog.INFO_HIGH)

propagate_state(m.fs.s11)
m.fs.splitter_2.split_fraction[0, "to_FG_cooler"].fix(split_frac_fg_cooler)
m.fs.splitter_2.initialize(outlvl=outlvl)

propagate_state(m.fs.s12)
# m.fs.FG_cooler.outlet.temperature.fix(483.15)
m.fs.FG_cooler.deltaP.fix(-0.06e6)
m.fs.FG_cooler.outlet.enth_mol.fix(htpx(
    T=440.15*units.K,
    P=value(m.fs.FG_cooler.inlet.pressure[0]+m.fs.FG_cooler.deltaP[0])
    * units.Pa))
m.fs.FG_cooler.initialize(outlvl=outlvl)


propagate_state(m.fs.s13)

m.fs.LTR_pseudo_tube.deltaP.fix(-0.07e6)
m.fs.LTR_pseudo_tube.heat_duty[0].\
    fix(-value(m.fs.LTR_pseudo_shell.heat_duty[0]))
m.fs.LTR_pseudo_tube.initialize(outlvl=outlvl)

# Add constraint heats of the LTR_pseudo shell and tube
m.fs.LTR_pseudo_tube.heat_duty[0].unfix()
m.fs.c1 = Constraint(expr=m.fs.LTR_pseudo_shell.heat_duty[0] ==
                     -m.fs.LTR_pseudo_tube.heat_duty[0])

propagate_state(m.fs.s10)
propagate_state(m.fs.s14)
propagate_state(m.fs.s15)

m.fs.mixer.initialize(outlvl=outlvl)

propagate_state(m.fs.s16)

m.fs.HTR_pseudo_tube.heat_duty[0].\
    fix(-value(m.fs.HTR_pseudo_shell.heat_duty[0]))
m.fs.HTR_pseudo_tube.deltaP.fix(-0.07e6)
m.fs.HTR_pseudo_tube.initialize(outlvl=outlvl)

m.fs.HTR_pseudo_tube.heat_duty[0].unfix()
m.fs.c2 = Constraint(expr=m.fs.HTR_pseudo_shell.heat_duty[0] ==
                     -m.fs.HTR_pseudo_tube.heat_duty[0])

TransformationFactory("network.expand_arcs").apply_to(m.fs)

if args.deltaTmin:
    # HTR
    m.fs.HTR_pseudo_shell.outlet.enth_mol.unfix()
    m.fs.c3 = Constraint(
        expr=m.fs.HTR_pseudo_shell.control_volume.properties_out[0].
        temperature == m.fs.HTR_pseudo_tube.control_volume.properties_in[0].
        temperature + deltaTmin_HTR)

    # LTR
    m.fs.LTR_pseudo_shell.outlet.enth_mol.unfix()
    m.fs.c4 = Constraint(
        expr=m.fs.LTR_pseudo_shell.control_volume.properties_out[0].
        temperature == m.fs.LTR_pseudo_tube.control_volume.properties_in[0].
        temperature + deltaTmin_LTR)
if degrees_of_freedom(m) == 0:
    solver.solve(m, tee=tee)
else:
    raise Exception("Degrees of freedom is not 0.")


# Print reports
for i in m.fs.component_objects(Block):
    if isinstance(i, UnitModelBlockData):
        i.report()

m.gross_cycle_power_output = \
    Expression(expr=(-m.fs.turbine.work_mechanical[0] -
               m.fs.main_compressor_stage_1.work_mechanical[0] -
               m.fs.main_compressor_stage_2.work_mechanical[0] -
               m.fs.bypass_compressor.work_mechanical[0]))

# account for generator loss = 1.5% of gross power output
m.net_cycle_power_output = Expression(expr=0.985*m.gross_cycle_power_output)

m.total_cycle_power_input = Expression(
    expr=(m.fs.boiler.heat_duty[0]
          + m.fs.FG_cooler.heat_duty[0]))

m.cycle_efficiency = Expression(
    expr=m.net_cycle_power_output/m.total_cycle_power_input*100)

# Expression to compute recovered duty in recuperators
m.recuperator_duty = Expression(
    expr=(m.fs.HTR_pseudo_tube.heat_duty[0] +
          m.fs.LTR_pseudo_tube.heat_duty[0]))

print('Net power = ', value(m.net_cycle_power_output)*1e-6, ' MW')
print('Total power input = ', value(m.total_cycle_power_input)*1e-6, ' MW')
print('cycle efficiency = ',
      value(m.cycle_efficiency), ' %')
print('Total recuperator heat duty = ',
      value(m.recuperator_duty)*1e-6, ' MW')
