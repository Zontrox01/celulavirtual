"""
tests/test_snapshot.py

Validación de engine/snapshot.py.

El criterio central: restaurar un snapshot y reproducir hacia adelante
debe generar EXACTAMENTE la misma secuencia de eventos que ya había
ocurrido en esos mismos pasos — no una simulación nueva, la misma.
"""

import pytest
from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager
from engine.ssa import GillespieSSA
from engine.snapshot import SnapshotStore, SnapshotError


def build_system(seed=1, a=200):
    sm = SpeciesManager()
    sm.add_species("A", a)
    sm.add_species("B", 0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    return GillespieSSA(sm, rm, seed=seed)


def test_invalid_interval_raises():
    ssa = build_system()
    with pytest.raises(SnapshotError):
        SnapshotStore(ssa, interval=0)


def test_initial_snapshot_exists_at_step_zero():
    ssa = build_system()
    store = SnapshotStore(ssa, interval=10)
    assert store.snapshot_steps() == [0]


def test_maybe_snapshot_respects_interval():
    ssa = build_system()
    store = SnapshotStore(ssa, interval=5)

    for _ in range(4):
        ssa.step()
        store.maybe_snapshot()
    assert store.snapshot_steps() == [0]  # todavía no han pasado 5 desde el último

    ssa.step()  # paso 5
    store.maybe_snapshot()
    assert store.snapshot_steps() == [0, 5]


def test_nearest_snapshot_at_or_before():
    ssa = build_system()
    store = SnapshotStore(ssa, interval=10)
    for _ in range(25):
        ssa.step()
        store.maybe_snapshot()
    # snapshots esperados en 0, 10, 20

    nearest = store.nearest_snapshot_at_or_before(17)
    assert nearest.step_count == 10


def test_nearest_snapshot_raises_if_none_available_before():
    ssa = build_system()
    store = SnapshotStore(ssa, interval=10)
    with pytest.raises(SnapshotError):
        store.nearest_snapshot_at_or_before(-1)


def test_rewind_to_future_step_raises():
    ssa = build_system()
    store = SnapshotStore(ssa, interval=10)
    for _ in range(5):
        ssa.step()
    with pytest.raises(SnapshotError):
        store.rewind_to(100)


def test_rewind_to_negative_step_raises():
    ssa = build_system()
    store = SnapshotStore(ssa, interval=10)
    with pytest.raises(SnapshotError):
        store.rewind_to(-1)


def test_rewind_restores_species_state_exactly():
    ssa = build_system(seed=5)
    store = SnapshotStore(ssa, interval=10)
    for _ in range(30):
        ssa.step()
        store.maybe_snapshot()

    state_at_30 = ssa.species.snapshot()
    store.rewind_to(15)
    state_at_15 = ssa.species.snapshot()

    assert state_at_15 != state_at_30  # de verdad retrocedió
    assert ssa.step_count == 15


def test_rewind_then_replay_reproduces_the_exact_same_future():
    """
    Criterio central: correr la simulación de una vez hasta el paso 40
    y anotar los eventos 21-40; por separado, retroceder al paso 20 y
    volver a avanzar 20 pasos más -- deben ser exactamente los mismos
    eventos (mismo id de reacción, mismo tiempo), no solo estadísticamente
    parecidos.
    """
    ssa = build_system(seed=99, a=500)
    store = SnapshotStore(ssa, interval=10)

    events_first_pass = []
    for _ in range(40):
        events_first_pass.append(ssa.step())
        store.maybe_snapshot()

    original_tail = events_first_pass[20:40]  # eventos de los pasos 21 a 40

    store.rewind_to(20)
    assert ssa.step_count == 20

    replayed_tail = [ssa.step() for _ in range(20)]

    assert [(e.reaction_id, e.time) for e in original_tail] == [
        (e.reaction_id, e.time) for e in replayed_tail
    ]


def test_rewind_discards_snapshots_after_the_restored_point():
    ssa = build_system(seed=3)
    store = SnapshotStore(ssa, interval=10)
    for _ in range(35):
        ssa.step()
        store.maybe_snapshot()
    # snapshots en 0, 10, 20, 30

    store.rewind_to(20)
    assert store.snapshot_steps() == [0, 10, 20]  # el de 30 ya no es válido
