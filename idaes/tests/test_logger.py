import pytest
import idaes.logger as idaeslog
import logging
import pyomo.environ as pyo

__author__ = "John Eslick"

def test_get_idaes_logger(caplog):
    caplog.set_level(logging.DEBUG)
    log = idaeslog.getLogger("My Test Logger 1")
    log.setLevel(logging.DEBUG)
    assert log.name == "idaes.My Test Logger 1"
    log.flowsheet("Hello!")
    log.unit("Hello!")
    log.info("Hello!")
    log.unit_high("Hello!")
    log.cv("Hello!")
    log.prop("Hello!")
    assert caplog.records[0].levelname == "FLOWSHEET"
    assert caplog.records[1].levelname == "UNIT"
    assert caplog.records[2].levelname == "INFO"
    assert caplog.records[3].levelname == "UNIT_HIGH"
    assert caplog.records[4].levelname == "CV"
    assert caplog.records[5].levelname == "PROP"

    log = idaeslog.getLogger("idaes.My Test Logger 2")
    assert log.name == "idaes.My Test Logger 2"

def test_get_model_logger(caplog):
    log = idaeslog.getModelLogger("My Model 1")
    assert isinstance(log, logging.Logger)
    assert log.name == "idaes.model.My Model 1"
    caplog.set_level(idaeslog.UNIT)
    log.flowsheet("Hello! from flowsheet")
    log.unit("Hello! from unit")
    log.info("Hello! from info")
    log.unit_high("Hello! from unit high")
    log.cv("Hello! from cv")
    log.prop("Hello! from prop")
    for record in caplog.records:
        assert record.message in ["Hello! from flowsheet", "Hello! from unit"]
    log = idaeslog.getModelLogger("idaes.My Model 2")
    assert log.name == "idaes.model.My Model 2"

def test_get_init_logger():
    log = idaeslog.getInitLogger("My Init 1")
    assert log.name == "idaes.init.My Init 1"

def test_solver_condition():
    # test the results that can be tested without a solver
    assert idaeslog.condition(None) == "Error, no result"
    assert idaeslog.condition("something else") == "something else"

@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
def test_solver_condition2():
    solver = pyo.SolverFactory('ipopt')
    model = pyo.ConcreteModel("Solver Result Test Model")
    model.x = pyo.Var([1,2])
    model.y = pyo.Var(initialize=5)
    model.x.fix(2)
    model.y.unfix()
    model.c = pyo.Constraint(expr=model.x[1] + model.x[2]==model.y)
    res = solver.solve(model)
    assert idaeslog.condition(res).startswith("optimal") #better solve
    model.c2 = pyo.Constraint(expr=model.x[1]==model.y)
    res = solver.solve(model)
    assert idaeslog.condition(res).startswith("other") #too few degrees of freedom

@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
def test_solver_log(caplog):
    solver = pyo.SolverFactory('ipopt')
    model = pyo.ConcreteModel("Solver Result Test Model")
    model.x = pyo.Var([1,2])
    model.y = pyo.Var(initialize=5)
    model.x.fix(2)
    model.y.unfix()
    model.c = pyo.Constraint(expr=model.x[1] + model.x[2]==model.y)

    log = idaeslog.getLogger("solver")
    caplog.set_level(idaeslog.DEBUG)
    log.setLevel(idaeslog.DEBUG)

    idaeslog.solver_capture_on()
    with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
        res = solver.solve(model, tee=True)
    assert(not slc.thread.is_alive()) # make sure logging thread is down
    s = ""
    for record in caplog.records:
        s += record.message
    assert "Optimal" in s

    # test that an excpetion still results in tread terminating right
    try:
        with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
            res = solver.solve(modelf, tee=True)
    except NameError:
        pass # expect name error
    assert(not slc.thread.is_alive()) # make sure logging thread is down
