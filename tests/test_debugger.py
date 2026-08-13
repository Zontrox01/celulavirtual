"""
tests/test_debugger.py

Validación de engine/debugger.py.

Incluye el criterio de éxito de la sección 11 del whitepaper: el modo
paso a paso no debe alterar la distribución estadística de resultados
respecto al modo run() normal, ya que el depurador solo debe cambiar
la forma de observar la simulación, no su física.
"""

import pytest
from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager
from engine.ssa import GillespieSSA, SimulationEnded
from engine.debugger import Debugger, DebuggerEndedError


def build_system(seed, a=20, extra_reaction=False):
    sm = SpeciesManager()
    sm.add_species("A", a)
    sm.add_species("B", 0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    if extra_reaction:
        sm.add_species("ATP", 100)
        sm.add_species("C", 0)
        rm.add_reaction(Reaction("consume_ATP", {"ATP": 1}, {"C": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=seed)
    return sm, rm, ssa


# ---------------------------------------------------------------------
# Criterio de validación central: paso a paso == run() normal
# ---------------------------------------------------------------------

def test_debugger_matches_direct_run_with_same_seed():
    """
    Recorrer la simulación con el depurador (step a step) debe producir
    exactamente la misma secuencia de eventos que ssa.run() con el
    mismo seed — el depurador no debe introducir ninguna diferencia
    estadística, solo control de ejecución.
    """
    _, _, ssa_direct = build_system(seed=42)
    direct_events = ssa_direct.run()

    _, _, ssa_debugged = build_system(seed=42)
    debugger = Debugger(ssa_debugged)
    debugged_events = []
    while not debugger.has_ended:
        batch = debugger.step(n=1)
        if not batch:
            break
        debugged_events.extend(batch)

    assert [(e.reaction_id, e.time) for e in direct_events] == [
        (e.reaction_id, e.time) for e in debugged_events
    ]


# ---------------------------------------------------------------------
# Ejecución paso a paso
# ---------------------------------------------------------------------

def test_step_default_advances_one_reaction():
    sm, rm, ssa = build_system(seed=1)
    debugger = Debugger(ssa)

    events = debugger.step()
    assert len(events) == 1
    assert sm.get_count("A") == 19


def test_step_n_advances_multiple_reactions():
    sm, rm, ssa = build_system(seed=2)
    debugger = Debugger(ssa)

    events = debugger.step(n=7)
    assert len(events) == 7
    assert sm.get_count("A") == 13


def test_step_raises_after_simulation_ended():
    sm, rm, ssa = build_system(seed=3, a=2)
    debugger = Debugger(ssa)

    debugger.step(n=2)  # agota las 2 moléculas de A
    assert debugger.has_ended is False  # todavía no se ha intentado un paso de más

    empty = debugger.step()  # este intento choca con PropensityExhausted internamente
    assert empty == []
    assert debugger.has_ended is True

    with pytest.raises(DebuggerEndedError):
        debugger.step()


def test_max_time_pause_does_not_mark_debugger_as_permanently_ended():
    """
    Regresión: si el GillespieSSA subyacente tiene max_time y lo
    alcanza, el depurador NO debe marcarse como terminado de forma
    permanente (a diferencia de agotar la propensidad total).
    """
    sm, rm, _ = build_system(seed=4, a=10_000)
    ssa = GillespieSSA(sm, rm, seed=4, max_time=0.001)
    debugger = Debugger(ssa)

    debugger.run_until_breakpoint(max_steps=100_000)
    assert debugger.has_ended is False

    ssa.max_time = 0.002
    more_events = debugger.step(n=5)
    assert len(more_events) > 0, "debería poder seguir avanzando tras subir max_time"


# ---------------------------------------------------------------------
# Breakpoints
# ---------------------------------------------------------------------

def test_step_stops_early_when_breakpoint_triggers():
    sm, rm, ssa = build_system(seed=4)
    debugger = Debugger(ssa)
    debugger.set_breakpoint("pocas_A", lambda s: s["A"] <= 15)

    events = debugger.step(n=100)  # pide 100, pero debe pararse mucho antes

    assert sm.get_count("A") == 15
    assert debugger.last_triggered_breakpoint == "pocas_A"
    assert len(events) < 100


def test_run_until_breakpoint_stops_at_condition():
    sm, rm, ssa = build_system(seed=5)
    debugger = Debugger(ssa)
    debugger.set_breakpoint("mitad", lambda s: s["A"] <= 10)

    debugger.run_until_breakpoint()

    assert sm.get_count("A") == 10
    assert debugger.last_triggered_breakpoint == "mitad"


def test_run_until_breakpoint_without_breakpoints_respects_max_steps():
    sm, rm, ssa = build_system(seed=6)
    debugger = Debugger(ssa)  # sin breakpoints definidos

    events = debugger.run_until_breakpoint(max_steps=5)

    assert len(events) == 5
    assert debugger.last_triggered_breakpoint is None


def test_remove_breakpoint():
    sm, rm, ssa = build_system(seed=7)
    debugger = Debugger(ssa)
    debugger.set_breakpoint("nunca_debe_pararse_aqui", lambda s: s["A"] <= 19)
    debugger.remove_breakpoint("nunca_debe_pararse_aqui")

    events = debugger.step(n=5)
    assert len(events) == 5  # no se paró, porque el breakpoint fue eliminado


def test_clear_breakpoints():
    sm, rm, ssa = build_system(seed=8)
    debugger = Debugger(ssa)
    debugger.set_breakpoint("bp1", lambda s: s["A"] <= 19)
    debugger.set_breakpoint("bp2", lambda s: s["A"] <= 18)
    debugger.clear_breakpoints()

    events = debugger.step(n=5)
    assert len(events) == 5


# ---------------------------------------------------------------------
# Inspección
# ---------------------------------------------------------------------

def test_current_event_and_current_state():
    sm, rm, ssa = build_system(seed=9)
    debugger = Debugger(ssa)

    assert debugger.current_event is None  # todavía no se ha dado ningún paso

    debugger.step()
    assert debugger.current_event is not None
    assert debugger.current_event.reaction_id == "A_to_B"
    assert debugger.current_state() == {"A": 19, "B": 1}


def test_pending_propensities_includes_all_candidate_reactions():
    sm, rm, ssa = build_system(seed=10, extra_reaction=True)
    debugger = Debugger(ssa)

    props = debugger.pending_propensities()
    assert set(props.keys()) == {"A_to_B", "consume_ATP"}
    assert props["A_to_B"] == 20.0
    assert props["consume_ATP"] == 100.0


def test_history_accumulates_across_multiple_step_calls():
    sm, rm, ssa = build_system(seed=11)
    debugger = Debugger(ssa)

    debugger.step(n=3)
    debugger.step(n=2)

    assert len(debugger.history()) == 5


def test_repr_contains_key_information():
    sm, rm, ssa = build_system(seed=12)
    debugger = Debugger(ssa)
    debugger.set_breakpoint("bp1", lambda s: False)
    debugger.step(n=2)

    text = repr(debugger)
    assert "eventos=2" in text
    assert "bp1" in text
