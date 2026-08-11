"""
tests/test_ssa.py

Validación de engine/ssa.py.

Dos tipos de test, siguiendo el criterio de éxito de la Fase 0
(whitepaper, sección 11):
  1. Mecánicos: la simulación aplica las reacciones correctamente,
     conserva especies, respeta max_time, termina cuando debe.
  2. Estadísticos: la distribución de tiempos de espera y de elección
     de reacción converge al comportamiento teórico esperado del
     método directo de Gillespie (no solo "no crashea").
"""

import math
import statistics

import pytest

from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager
from engine.ssa import GillespieSSA, SimulationEnded, PropensityExhausted


def make_species(**counts):
    sm = SpeciesManager()
    for species_id, count in counts.items():
        sm.add_species(species_id, count)
    return sm


# ---------------------------------------------------------------------
# Tests mecánicos
# ---------------------------------------------------------------------

def test_single_step_applies_the_only_reaction():
    sm = make_species(A=10, B=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=1)

    event = ssa.step()

    assert event.reaction_id == "A_to_B"
    assert event.step == 1
    assert event.time > 0.0
    assert sm.get_count("A") == 9
    assert sm.get_count("B") == 1


def test_conservation_across_full_run():
    # A -> B, cantidad total (A+B) debe mantenerse constante en todo momento
    sm = make_species(A=50, B=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=7)

    for _ in ssa:
        assert sm.get_count("A") + sm.get_count("B") == 50

    assert sm.get_count("A") == 0
    assert sm.get_count("B") == 50


def test_time_is_monotonically_increasing():
    sm = make_species(A=30, B=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=3)

    last_time = -1.0
    for event in ssa:
        assert event.time > last_time
        last_time = event.time


def test_simulation_ends_when_propensity_reaches_zero():
    sm = make_species(A=5, B=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=2)

    events = ssa.run()
    assert len(events) == 5
    assert sm.get_count("A") == 0

    with pytest.raises(PropensityExhausted):
        ssa.step()


def test_max_time_stops_simulation():
    sm = make_species(A=10_000, B=0)  # propensidad alta -> pasos muy frecuentes
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=5, max_time=0.001)

    events = ssa.run()
    assert all(e.time <= 0.001 for e in events)
    assert ssa.time <= 0.001


def test_max_time_is_a_resumable_pause_not_a_permanent_end():
    """
    Regresión: alcanzar max_time NO debe marcar la simulación como
    terminada de forma permanente. Subir max_time y volver a llamar a
    run()/step() debe poder continuar generando eventos nuevos.
    """
    sm = make_species(A=10_000, B=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=5, max_time=0.001)

    first_batch = ssa.run()
    assert len(first_batch) > 0

    ssa.max_time = 0.002
    second_batch = ssa.run()
    assert len(second_batch) > 0, "la simulación debería poder continuar tras subir max_time"


def test_exhaustion_remains_permanent_unlike_max_time():
    """
    A diferencia de max_time, agotar la propensidad total (sin
    reactivos suficientes para ninguna reacción) sí debe ser
    permanente: PropensityExhausted en cada intento posterior.
    """
    sm = make_species(A=2, B=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=1)

    ssa.run()  # agota las 2 moléculas de A
    with pytest.raises(PropensityExhausted):
        ssa.step()
    with pytest.raises(PropensityExhausted):
        ssa.step()  # sigue siendo permanente en llamadas sucesivas


def test_run_respects_max_steps():
    sm = make_species(A=100, B=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=9)

    events = ssa.run(max_steps=10)
    assert len(events) == 10
    assert sm.get_count("A") == 90


def test_reproducibility_with_same_seed():
    def run_once():
        sm = make_species(A=20, B=0)
        rm = ReactionManager()
        rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
        ssa = GillespieSSA(sm, rm, seed=123)
        return [(e.reaction_id, e.time) for e in ssa]

    assert run_once() == run_once()


def test_event_records_all_candidate_propensities():
    sm = make_species(A=5, B=5, C=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_C", {"A": 1}, {"C": 1}, rate_constant=1.0))
    rm.add_reaction(Reaction("B_to_C", {"B": 1}, {"C": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=4)

    event = ssa.step()
    # el evento debe registrar la propensidad de AMBAS reacciones candidatas,
    # no solo la que disparó (necesario para depurador y modo forense)
    assert set(event.propensities.keys()) == {"A_to_C", "B_to_C"}


# ---------------------------------------------------------------------
# Tests estadísticos
# ---------------------------------------------------------------------

def test_waiting_time_matches_exponential_distribution():
    """
    Con una propensidad total constante (reactivo abundante, un solo
    paso por simulación), el tiempo de espera debe distribuirse
    exponencialmente con media 1/lambda. Se comprueba con muchas
    repeticiones independientes y una tolerancia estadística generosa.
    """
    total_propensity = 50.0  # A=50, k=1.0 -> propensidad = 50
    n_samples = 4000
    taus = []

    for i in range(n_samples):
        sm = make_species(A=50, B=0)
        rm = ReactionManager()
        rm.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
        ssa = GillespieSSA(sm, rm, seed=i)
        event = ssa.step()
        taus.append(event.time)

    expected_mean = 1.0 / total_propensity
    sample_mean = statistics.mean(taus)
    sample_std = statistics.stdev(taus)

    # error estándar de la media, para una tolerancia basada en la propia varianza muestral
    standard_error = sample_std / math.sqrt(n_samples)
    tolerance = 5 * standard_error

    assert abs(sample_mean - expected_mean) < tolerance, (
        f"media muestral {sample_mean:.6f} fuera de tolerancia "
        f"respecto a la media teórica {expected_mean:.6f} (tol={tolerance:.6f})"
    )


def test_reaction_choice_proportional_to_propensity():
    """
    Con dos reacciones de propensidades distintas (3:1), la frecuencia
    empírica de elección de cada una debe converger a esa proporción.
    """
    n_samples = 4000
    counts = {"fast": 0, "slow": 0}

    for i in range(n_samples):
        sm = make_species(A=10, B=10, C=0)
        rm = ReactionManager()
        # "fast" tiene el triple de constante que "slow" -> propensidad 3:1
        rm.add_reaction(Reaction("fast", {"A": 1}, {"C": 1}, rate_constant=3.0))
        rm.add_reaction(Reaction("slow", {"B": 1}, {"C": 1}, rate_constant=1.0))
        ssa = GillespieSSA(sm, rm, seed=i)
        event = ssa.step()
        counts[event.reaction_id] += 1

    fraction_fast = counts["fast"] / n_samples
    expected_fraction = 3.0 / 4.0  # 3 / (3+1)

    # tolerancia basada en el error estándar de una proporción binomial
    standard_error = math.sqrt(expected_fraction * (1 - expected_fraction) / n_samples)
    tolerance = 5 * standard_error

    assert abs(fraction_fast - expected_fraction) < tolerance, (
        f"fracción empírica de 'fast' = {fraction_fast:.4f}, "
        f"esperada ~{expected_fraction:.4f} (tol={tolerance:.4f})"
    )