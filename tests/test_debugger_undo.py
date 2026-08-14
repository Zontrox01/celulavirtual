"""
tests/test_debugger_undo.py

Validación de la integración de retroceso (undo) en engine/debugger.py
(whitepaper, secciones 4.4.b y 5.3). Complementa tests/test_debugger.py
y tests/test_snapshot.py.
"""

import pytest
from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager
from engine.ssa import GillespieSSA
from engine.debugger import Debugger, DebuggerEndedError


def build_system(seed=1, a=200):
    sm = SpeciesManager()
    sm.add_species("A", a)
    sm.add_species("B", 0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    return sm, GillespieSSA(sm, rm, seed=seed)


def test_snapshots_disabled_by_default():
    _, ssa = build_system()
    debugger = Debugger(ssa)
    assert debugger.snapshots_enabled is False


def test_undo_without_snapshots_enabled_raises():
    _, ssa = build_system()
    debugger = Debugger(ssa)
    debugger.step(n=5)

    with pytest.raises(DebuggerEndedError):
        debugger.undo_to(2)


def test_undo_restores_species_state():
    sm, ssa = build_system(seed=7)
    debugger = Debugger(ssa, enable_snapshots=True, snapshot_interval=10)

    debugger.step(n=30)
    state_at_30 = sm.snapshot()

    debugger.undo_to(15)
    state_at_15 = sm.snapshot()

    assert state_at_15 != state_at_30
    assert ssa.step_count == 15


def test_undo_truncates_debugger_history():
    sm, ssa = build_system(seed=7)
    debugger = Debugger(ssa, enable_snapshots=True, snapshot_interval=10)

    debugger.step(n=30)
    assert len(debugger.history()) == 30

    debugger.undo_to(15)
    assert len(debugger.history()) == 15
    assert all(e.step <= 15 for e in debugger.history())


def test_undo_then_continue_reproduces_the_same_original_future():
    """
    El test decisivo: retroceder y avanzar de nuevo debe reproducir
    exactamente lo que ya había pasado la primera vez, evento a evento.
    """
    sm, ssa = build_system(seed=42, a=500)
    debugger = Debugger(ssa, enable_snapshots=True, snapshot_interval=10)

    debugger.step(n=40)
    original_events = debugger.history()
    original_tail = original_events[20:40]

    debugger.undo_to(20)
    replayed = []
    while len(replayed) < 20:
        batch = debugger.step(n=1)
        replayed.extend(batch)

    assert [(e.reaction_id, e.time) for e in original_tail] == [
        (e.reaction_id, e.time) for e in replayed
    ]


def test_undo_resets_ended_flag_when_going_back_before_exhaustion():
    sm, ssa = build_system(seed=1, a=5)
    debugger = Debugger(ssa, enable_snapshots=True, snapshot_interval=1)

    debugger.step(n=5)  # agota las 5 moléculas de A
    debugger.step()  # intento fallido -> has_ended True
    assert debugger.has_ended is True

    debugger.undo_to(3)
    assert debugger.has_ended is False

    events = debugger.step(n=2)  # debería poder seguir avanzando
    assert len(events) == 2


def test_available_snapshot_steps_empty_when_disabled():
    _, ssa = build_system()
    debugger = Debugger(ssa)
    assert debugger.available_snapshot_steps() == []


def test_available_snapshot_steps_reflects_interval():
    _, ssa = build_system()
    debugger = Debugger(ssa, enable_snapshots=True, snapshot_interval=5)
    debugger.step(n=12)
    assert debugger.available_snapshot_steps() == [0, 5, 10]


# ---------------------------------------------------------------------
# Regresión: nada de lo ya validado en test_debugger.py debe romperse
# ---------------------------------------------------------------------

def test_regression_debugger_without_snapshots_behaves_as_before():
    sm, ssa = build_system(seed=1)
    debugger = Debugger(ssa)  # sin enable_snapshots, comportamiento por defecto

    events = debugger.step(n=10)
    assert len(events) == 10
    assert debugger.current_event.reaction_id == "A_to_B"
    assert debugger.snapshots_enabled is False
