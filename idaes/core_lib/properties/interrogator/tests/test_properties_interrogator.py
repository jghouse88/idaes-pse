#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for Property Interrogator Tool

@author: alee
"""
import pytest

from pyomo.environ import ConcreteModel

from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.unit_models import Flash, HeatExchanger1D
from idaes.unit_models.pressure_changer import (PressureChanger,
                                                ThermodynamicAssumption)
from idaes.unit_models.separator import Separator, SplittingType
from idaes.core_lib.properties.interrogator import PropertyInterrogatorBlock


def test_interrogator_parameter_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    # Check that parameter block has expected attributes
    assert isinstance(m.fs.params.required_properties, dict)
    assert len(m.fs.params.required_properties) == 0


def test_interrogator_state_block_methods():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.props = m.fs.params.state_block_class(
            [0],
            default={"parameters": m.fs.params})

    # Check get_term methods return an unindexed dummy var
    assert m.fs.props[0].get_material_flow_terms("Liq", "A") is \
        m.fs.props[0]._dummy_var
    assert m.fs.props[0].get_enthalpy_flow_terms("Liq") is \
        m.fs.props[0]._dummy_var
    assert m.fs.props[0].get_material_density_terms("Liq", "A") is \
        m.fs.props[0]._dummy_var
    assert m.fs.props[0].get_energy_density_terms("Liq") is \
        m.fs.props[0]._dummy_var

    # Check that get_term calls were logged correctly
    assert m.fs.params.required_properties == {
            "material flow terms": ["fs.props"],
            "material density terms": ["fs.props"],
            "enthalpy flow terms": ["fs.props"],
            "energy density terms": ["fs.props"]}


def test_interrogator_state_block_unindexed_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.props = m.fs.params.state_block_class(
            [0],
            default={"parameters": m.fs.params})

    # Check get_term methods return an unindexed dummy var
    assert m.fs.props[0].prop_unindexed is \
        m.fs.props[0]._dummy_var

    # Call again to make sure duplicates are skipped in required_properties
    assert m.fs.props[0].prop_unindexed is \
        m.fs.props[0]._dummy_var

    # Check that get_term calls were logged correctly
    assert m.fs.params.required_properties == {
            "prop_unindexed": ["fs.props"]}


def test_interrogator_state_block_phase_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.props = m.fs.params.state_block_class(
            [0],
            default={"parameters": m.fs.params})

    # Check get_term methods return an unindexed dummy var
    assert m.fs.props[0].prop_phase["Liq"] is \
        m.fs.props[0]._dummy_var_phase["Liq"]
    assert m.fs.props[0].prop_phase["Vap"] is \
        m.fs.props[0]._dummy_var_phase["Vap"]

    # Check that get_term calls were logged correctly
    assert m.fs.params.required_properties == {
            "prop_phase": ["fs.props"]}


def test_interrogator_state_block_comp_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.props = m.fs.params.state_block_class(
            [0],
            default={"parameters": m.fs.params})

    # Check get_term methods return an unindexed dummy var
    assert m.fs.props[0].prop_comp["A"] is \
        m.fs.props[0]._dummy_var_comp["A"]
    assert m.fs.props[0].prop_comp["B"] is \
        m.fs.props[0]._dummy_var_comp["B"]

    # Check that get_term calls were logged correctly
    assert m.fs.params.required_properties == {
            "prop_comp": ["fs.props"]}


def test_interrogator_state_block_phase_comp_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.props = m.fs.params.state_block_class(
            [0],
            default={"parameters": m.fs.params})

    # Check get_term methods return an unindexed dummy var
    assert m.fs.props[0].prop_phase_comp["Liq", "A"] is \
        m.fs.props[0]._dummy_var_phase_comp["Liq", "A"]
    assert m.fs.props[0].prop_phase_comp["Vap", "B"] is \
        m.fs.props[0]._dummy_var_phase_comp["Vap", "B"]

    # Check that get_term calls were logged correctly
    assert m.fs.params.required_properties == {
            "prop_phase_comp": ["fs.props"]}


def test_interrogator_report_method():
    # Display method should return an TypeError
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.props = m.fs.params.state_block_class(
            [0],
            default={"parameters": m.fs.params})

    with pytest.raises(TypeError,
                       match="Models constructed using the Property "
                       "Interrogator package should not be used for report "
                       "methods."):
        m.fs.props.report()


def test_interrogator_initialize_method():
    # Initialize method should return an TypeError
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.props = m.fs.params.state_block_class(
            [0],
            default={"parameters": m.fs.params})

    with pytest.raises(TypeError,
                       match="Models constructed using the Property "
                       "Interrogator package cannot be used to solve a "
                       "flowsheet. Please rebuild your flowsheet using a "
                       "valid property package."):
        m.fs.props.initialize()


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.P01 = PressureChanger(default={
            "property_package": m.fs.params,
            "thermodynamic_assumption": ThermodynamicAssumption.isentropic})

    m.fs.HX02 = HeatExchanger1D(default={
                "shell_side": {"property_package": m.fs.params},
                "tube_side": {"property_package": m.fs.params}})

    m.fs.F03 = Flash(default={"property_package": m.fs.params})

    return m


def test_interrogate_flowsheet(model):
    assert model.fs.params.required_properties == {
            "material flow terms": ["fs.P01", "fs.HX02", "fs.F03"],
            "enthalpy flow terms": ["fs.P01", "fs.HX02", "fs.F03"],
            "material density terms": ["fs.P01", "fs.HX02"],
            "energy density terms": ["fs.P01", "fs.HX02"],
            "pressure": ["fs.P01", "fs.HX02", "fs.F03"],
            "temperature": ["fs.HX02"],
            "entr_mol": ["fs.P01"]}


def test_list_required_properties(model):
    prop_list = model.fs.params.list_required_properties()

    assert prop_list == ["material density terms",
                         "material flow terms",
                         "enthalpy flow terms",
                         "energy density terms",
                         "pressure",
                         "entr_mol",
                         "temperature"]


def test_list_models_requiring_property(model):
    for k in model.fs.params.required_properties.keys():
        model_list = model.fs.params.list_models_requiring_property(k)

        if k == "entr_mol":
            assert model_list == ["fs.P01"]
        elif k == "temperature":
            assert model_list == ["fs.HX02"]
        elif k in ["material density terms", "energy density terms"]:
            assert model_list == ["fs.P01", "fs.HX02"]
        else:
            assert model_list == ["fs.P01", "fs.HX02", "fs.F03"]


def test_list_properties_required_by_model_by_name(model):
    prop_list = model.fs.params.list_properties_required_by_model("fs.P01")
    assert prop_list == ["material density terms",
                         "material flow terms",
                         "enthalpy flow terms",
                         "energy density terms",
                         "pressure",
                         "entr_mol"]


def test_list_properties_required_by_model_by_object(model):
    prop_list = model.fs.params.list_properties_required_by_model(
            model.fs.P01)
    assert prop_list == ["material density terms",
                         "material flow terms",
                         "enthalpy flow terms",
                         "energy density terms",
                         "pressure",
                         "entr_mol"]


def test_list_properties_required_by_model_invalid_model(model):
    with pytest.raises(ValueError):
        model.fs.params.list_properties_required_by_model("foo")


def test_print_required_properties(model, capsys):
    model.fs.params.print_required_properties()

    captured = capsys.readouterr()
    assert captured.out == """
==========================================================================
Property Interrogator Summary

The Flowsheet requires the following properties (times required):

    material density terms                                               2
    material flow terms                                                  3
    enthalpy flow terms                                                  3
    energy density terms                                                 2
    pressure                                                             3
    entr_mol                                                             1
    temperature                                                          1

Note: User constraints may require additional properties which are not
reported here.
"""


def test_print_models_requiring_property(model, capsys):
    model.fs.params.print_models_requiring_property("pressure")

    captured = capsys.readouterr()
    assert captured.out == """
The following models in the Flowsheet require pressure:
    fs.P01
    fs.HX02
    fs.F03
"""


def test_print_properties_reqruied_by_model(model, capsys):
    model.fs.params.print_properties_required_by_model("fs.P01")

    captured = capsys.readouterr()
    assert captured.out == """
The following properties are required by model fs.P01:
    material density terms
    material flow terms
    enthalpy flow terms
    energy density terms
    pressure
    entr_mol
"""


# =============================================================================
# Extra checks to make sure Interrogator works with all Separator options
# All that really matters here is that these run, but do asserts anyway
def test_Separator_1():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.S01 = Separator(default={
                "property_package": m.fs.params,
                "material_balance_type": MaterialBalanceType.componentPhase,
                "split_basis": SplittingType.totalFlow,
                "outlet_list": ["a", "B", "c"],
                "ideal_separation": False,
                "has_phase_equilibrium": False})

    assert len(m.fs.params.required_properties) == 3
    for k in m.fs.params.required_properties.keys():
        assert k in ["material flow terms",
                     "pressure",
                     "temperature"]


def test_Separator_2():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.S01 = Separator(default={
                "property_package": m.fs.params,
                "material_balance_type": MaterialBalanceType.componentPhase,
                "split_basis": SplittingType.phaseFlow,
                "outlet_list": ["a", "B", "c"],
                "ideal_separation": False,
                "has_phase_equilibrium": False})

    assert len(m.fs.params.required_properties) == 3
    for k in m.fs.params.required_properties.keys():
        assert k in ["material flow terms",
                     "pressure",
                     "temperature"]


def test_Separator_3():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.S01 = Separator(default={
                "property_package": m.fs.params,
                "material_balance_type": MaterialBalanceType.componentPhase,
                "split_basis": SplittingType.componentFlow,
                "outlet_list": ["a", "B", "c"],
                "ideal_separation": False,
                "has_phase_equilibrium": False})

    assert len(m.fs.params.required_properties) == 3
    for k in m.fs.params.required_properties.keys():
        assert k in ["material flow terms",
                     "pressure",
                     "temperature"]


def test_ideal_Separator_1():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.S01 = Separator(default={
                "property_package": m.fs.params,
                "num_outlets": 4,
                "ideal_separation": True,
                "split_basis": SplittingType.phaseComponentFlow,
                "ideal_split_map": {("Vap", "A"): "outlet_1",
                                    ("Vap", "B"): "outlet_2",
                                    ("Liq", "A"): "outlet_3",
                                    ("Liq", "B"): "outlet_4"}})

    # Ideal Separator should require no property calls
    # Only dummy state variables are required
    assert len(m.fs.params.required_properties) == 0


def test_ideal_Separator_2():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.S01 = Separator(default={
                "property_package": m.fs.params,
                "num_outlets": 2,
                "ideal_separation": True,
                "split_basis": SplittingType.phaseFlow,
                "ideal_split_map": {("Vap"): "outlet_1",
                                    ("Liq"): "outlet_2"}})

    # Ideal Separator should require no property calls
    # Only dummy state variables are required
    assert len(m.fs.params.required_properties) == 0


def test_ideal_Separator_3():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = PropertyInterrogatorBlock()

    m.fs.S01 = Separator(default={
                "property_package": m.fs.params,
                "num_outlets": 2,
                "ideal_separation": True,
                "split_basis": SplittingType.componentFlow,
                "ideal_split_map": {("A"): "outlet_1",
                                    ("B"): "outlet_2"}})

    # Ideal Separator should require no property calls
    # Only dummy state variables are required
    assert len(m.fs.params.required_properties) == 0