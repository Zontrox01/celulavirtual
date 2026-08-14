"""
tests/test_cell_debug_undo.py

Validación del cierre de la limitación anotada en el whitepaper
(sección 8, Fase 5.3): cell.debug().undo_to() debe sincronizar también
el histórico de eventos de la propia Cell, no solo el del Debugger.
"""

import pytest
from cell import Cell
from engine.reactions import Reaction


def build_simple_cell(seed=1, a=200):
    cell = Cell(seed=seed)
    cell.add_species("A", a).add_species("B", 0)
    cell.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    return cell


def test_debug_accepts_enable_snapshots_backward_compatible_default():
    """Regresión: cell.debug() sin argumentos debe seguir funcionando como antes."""
    cell = build_simple_cell(seed=8)
    debugger = cell.debug()
    events = debugger.step(n=5)
    assert len(events) == 5
    assert debugger.snapshots_enabled is False


def test_undo_via_cell_debug_truncates_cell_events():
    cell = build_simple_cell(seed=7)
    debugger = cell.debug(enable_snapshots=True, snapshot_interval=10)

    debugger.step(n=30)
    assert len(cell.get_events()) == 30

    debugger.undo_to(15)

    assert len(cell.get_events()) == 15
    assert all(e.step <= 15 for e in cell.get_events())
    assert cell.get_events() == debugger.history()


def test_cell_time_reflects_the_rewound_state():
    cell = build_simple_cell(seed=7)
    debugger = cell.debug(enable_snapshots=True, snapshot_interval=10)

    debugger.step(n=30)
    time_at_30 = cell.time

    debugger.undo_to(15)

    assert cell.time < time_at_30
    assert cell.time == debugger.time


def test_cell_get_trajectory_is_consistent_after_undo():
    cell = build_simple_cell(seed=7)
    debugger = cell.debug(enable_snapshots=True, snapshot_interval=10)

    debugger.step(n=30)
    debugger.undo_to(15)

    df = cell.get_trajectory()
    assert len(df) == 16  # fila inicial (t=0) + 15 eventos
    assert (df["A"] + df["B"] == cell.species.get_count("A") + cell.species.get_count("B")).all()
    # el estado final de la trayectoria debe coincidir con el estado real actual
    assert df.iloc[-1]["A"] == cell.species.get_count("A")
    assert df.iloc[-1]["B"] == cell.species.get_count("B")


def test_cell_run_after_debug_undo_continues_from_the_rewound_point():
    """
    Tras retroceder vía el depurador, seguir con cell.run() normal debe
    continuar desde el punto reanudado, sin duplicar ni perder eventos.
    """
    cell = build_simple_cell(seed=7)
    debugger = cell.debug(enable_snapshots=True, snapshot_interval=10)

    debugger.step(n=30)
    debugger.undo_to(15)

    new_events = cell.run(max_steps=10)

    assert len(new_events) == 10
    assert len(cell.get_events()) == 25  # 15 recortados + 10 nuevos
    # los pasos deben ser consecutivos, sin huecos ni duplicados
    steps = [e.step for e in cell.get_events()]
    assert steps == list(range(1, 26))


def test_multiple_undo_and_redo_cycles_stay_consistent():
    cell = build_simple_cell(seed=9, a=500)
    debugger = cell.debug(enable_snapshots=True, snapshot_interval=5)

    debugger.step(n=50)
    debugger.undo_to(30)
    assert len(cell.get_events()) == 30

    debugger.step(n=10)  # avanza de nuevo desde 30 -> 40 (reproduce lo mismo que ya había pasado)
    assert len(cell.get_events()) == 40

    debugger.undo_to(10)
    assert len(cell.get_events()) == 10
    assert cell.get_events() == debugger.history()


def test_undo_without_snapshots_raises_and_does_not_touch_cell_events():
    cell = build_simple_cell(seed=1)
    debugger = cell.debug()  # sin enable_snapshots

    debugger.step(n=5)
    events_before = cell.get_events()

    from engine.debugger import DebuggerEndedError
    with pytest.raises(DebuggerEndedError):
        debugger.undo_to(2)

    assert cell.get_events() == events_before  # nada debe haber cambiado
