"""
tests/test_trajectory.py

Validación de engine/trajectory.py.
"""

import pandas as pd
import pytest
from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager
from engine.ssa import GillespieSSA
from engine.trajectory import build_trajectory, TrajectoryError


def build_system(seed, initial_a=10):
    sm = SpeciesManager()
    sm.add_species("A", initial_a)
    sm.add_species("B", 0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=seed)
    return sm, rm, ssa


def test_trajectory_has_initial_row_plus_one_per_event():
    initial_state = {"A": 10, "B": 0}
    sm, rm, ssa = build_system(seed=1)
    events = ssa.run()

    df = build_trajectory(initial_state, events, rm)

    assert len(df) == len(events) + 1
    assert df.iloc[0]["step"] == 0
    assert df.iloc[0]["time"] == 0.0
    assert pd.isna(df.iloc[0]["reaction_id"])


def test_trajectory_columns_are_step_time_reaction_and_species_sorted():
    initial_state = {"B": 0, "A": 10}  # orden de inserción invertido a propósito
    sm, rm, ssa = build_system(seed=2)
    events = ssa.run(max_steps=3)

    df = build_trajectory(initial_state, events, rm)

    assert list(df.columns) == ["step", "time", "reaction_id", "A", "B"]


def test_trajectory_conserves_total_across_all_rows():
    initial_state = {"A": 10, "B": 0}
    sm, rm, ssa = build_system(seed=3)
    events = ssa.run()

    df = build_trajectory(initial_state, events, rm)
    assert (df["A"] + df["B"] == 10).all()


def test_trajectory_matches_live_species_manager_at_the_end():
    initial_state = {"A": 10, "B": 0}
    sm, rm, ssa = build_system(seed=4)
    events = ssa.run()

    df = build_trajectory(initial_state, events, rm)

    assert df.iloc[-1]["A"] == sm.get_count("A")
    assert df.iloc[-1]["B"] == sm.get_count("B")


def test_trajectory_unknown_reaction_raises():
    from engine.ssa import SSAEvent

    initial_state = {"A": 5, "B": 0}
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))

    fake_event = SSAEvent(
        step=1, time=0.5, reaction_id="reaccion_inexistente",
        propensities={}, total_propensity=0.0,
    )

    with pytest.raises(TrajectoryError):
        build_trajectory(initial_state, [fake_event], rm)


def test_trajectory_species_missing_from_initial_state_raises():
    from engine.ssa import SSAEvent

    initial_state = {"A": 5}  # falta "B", que la reacción produce
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))

    fake_event = SSAEvent(
        step=1, time=0.5, reaction_id="A_to_B",
        propensities={}, total_propensity=0.0,
    )

    with pytest.raises(TrajectoryError):
        build_trajectory(initial_state, [fake_event], rm)


def test_trajectory_with_no_events_returns_only_initial_row():
    initial_state = {"A": 5, "B": 0}
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))

    df = build_trajectory(initial_state, [], rm)
    assert len(df) == 1
    assert df.iloc[0]["A"] == 5
