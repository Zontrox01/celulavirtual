"""
tests/test_cell.py

Validación de cell.py.
"""

import pytest
from cell import Cell
from engine.reactions import Reaction


def build_simple_cell(seed=1):
    cell = Cell(seed=seed)
    cell.add_species("A", 20).add_species("B", 0)
    cell.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    return cell


def test_add_species_and_reaction_are_chainable():
    cell = Cell()
    result = cell.add_species("A", 10).add_species("B", 0)
    assert result is cell
    assert cell.species.get_count("A") == 10
    assert cell.species.get_count("B") == 0


def test_run_executes_until_exhaustion_and_conserves_species():
    cell = build_simple_cell(seed=1)
    events = cell.run()

    assert len(events) == 20
    assert cell.species.get_count("A") == 0
    assert cell.species.get_count("B") == 20
    assert cell.time > 0.0


def test_run_respects_max_steps():
    cell = build_simple_cell(seed=2)
    events = cell.run(max_steps=5)

    assert len(events) == 5
    assert cell.species.get_count("A") == 15
    assert cell.species.get_count("B") == 5


def test_run_can_be_called_multiple_times_to_continue():
    cell = build_simple_cell(seed=3)
    first_batch = cell.run(max_steps=5)
    second_batch = cell.run(max_steps=5)

    assert len(first_batch) == 5
    assert len(second_batch) == 5
    assert cell.species.get_count("A") == 10
    assert len(cell.get_events()) == 10


def test_run_respects_max_time():
    cell = Cell(seed=4)
    cell.add_species("A", 10_000).add_species("B", 0)
    cell.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))

    events = cell.run(max_time=0.001)
    assert all(e.time <= 0.001 for e in events)
    assert cell.time <= 0.001


def test_debug_returns_a_working_debugger():
    cell = build_simple_cell(seed=8)
    debugger = cell.debug()

    events = debugger.step(n=5)
    assert len(events) == 5
    assert cell.species.get_count("A") == 15


def test_debug_events_are_synced_into_cell_history():
    cell = build_simple_cell(seed=9)
    debugger = cell.debug()
    debugger.step(n=4)

    assert len(cell.get_events()) == 4
    assert cell.get_events() == debugger.history()


def test_run_and_debug_can_alternate_without_losing_events():
    cell = build_simple_cell(seed=10)
    cell.run(max_steps=3)
    debugger = cell.debug()
    debugger.step(n=3)
    cell.run(max_steps=3)

    assert len(cell.get_events()) == 9
    assert cell.species.get_count("A") == 11  # 20 - 9

    df = cell.get_trajectory()
    assert len(df) == 10  # fila inicial + 9 eventos
    assert (df["A"] + df["B"] == 20).all()


def test_get_trajectory_raises_before_any_run():
    cell = build_simple_cell()
    with pytest.raises(RuntimeError):
        cell.get_trajectory()


def test_get_trajectory_reflects_initial_state_and_final_state():
    cell = build_simple_cell(seed=7)
    cell.run()  # corre hasta agotar A (20 -> 0)

    df = cell.get_trajectory()

    assert list(df.columns) == ["step", "time", "reaction_id", "A", "B"]
    assert len(df) == 21  # fila inicial (t=0) + 20 eventos

    first_row = df.iloc[0]
    assert first_row["step"] == 0
    assert first_row["time"] == 0.0
    assert first_row["A"] == 20
    assert first_row["B"] == 0

    last_row = df.iloc[-1]
    assert last_row["A"] == 0
    assert last_row["B"] == 20

    # conservación en todas las filas
    assert (df["A"] + df["B"] == 20).all()
    # el tiempo debe ser no decreciente
    assert (df["time"].diff().dropna() > 0).all()


def test_get_events_returns_a_copy_not_internal_reference():
    cell = build_simple_cell(seed=5)
    cell.run(max_steps=3)

    events = cell.get_events()
    events.append("algo_falso")

    assert len(cell.get_events()) == 3  # el histórico interno no debe verse afectado


def test_repr_contains_key_information():
    cell = build_simple_cell(seed=6)
    cell.run(max_steps=2)
    text = repr(cell)

    assert "especies=2" in text
    assert "reacciones=1" in text
    assert "eventos=2" in text
