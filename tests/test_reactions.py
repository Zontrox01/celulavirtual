"""
tests/test_reactions.py

Validación de engine/reactions.py.
"""

import pytest
from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager, InvalidReactionError


def make_species(**counts):
    sm = SpeciesManager()
    for species_id, count in counts.items():
        sm.add_species(species_id, count)
    return sm


def test_unimolecular_propensity():
    # A -> B, propensidad = k * count(A)
    sm = make_species(A=10, B=0)
    r = Reaction("r1", reactants={"A": 1}, products={"B": 1}, rate_constant=2.0)
    assert r.propensity(sm) == 20.0


def test_bimolecular_different_species_propensity():
    # A + B -> C, propensidad = k * count(A) * count(B)
    sm = make_species(A=5, B=4, C=0)
    r = Reaction("r2", reactants={"A": 1, "B": 1}, products={"C": 1}, rate_constant=1.5)
    assert r.propensity(sm) == pytest.approx(1.5 * 5 * 4)


def test_bimolecular_same_species_propensity():
    # 2A -> B, propensidad = k * C(count(A), 2) = k * count*(count-1)/2
    sm = make_species(A=6, B=0)
    r = Reaction("r3", reactants={"A": 2}, products={"B": 1}, rate_constant=1.0)
    assert r.propensity(sm) == pytest.approx(6 * 5 / 2)


def test_propensity_zero_if_not_enough_reactant():
    sm = make_species(A=1, B=0)
    r = Reaction("r4", reactants={"A": 2}, products={"B": 1}, rate_constant=10.0)
    assert r.propensity(sm) == 0.0


def test_propensity_unknown_reactant_species_raises():
    sm = make_species(A=5)
    r = Reaction("r5", reactants={"NOEXISTE": 1}, products={"A": 1}, rate_constant=1.0)
    with pytest.raises(InvalidReactionError):
        r.propensity(sm)


def test_apply_consumes_reactants_and_produces_products():
    sm = make_species(A=5, B=4, C=0)
    r = Reaction("r6", reactants={"A": 1, "B": 1}, products={"C": 2}, rate_constant=1.0)
    r.apply(sm)
    assert sm.get_count("A") == 4
    assert sm.get_count("B") == 3
    assert sm.get_count("C") == 2


def test_negative_rate_constant_raises():
    with pytest.raises(InvalidReactionError):
        Reaction("r7", reactants={"A": 1}, products={"B": 1}, rate_constant=-1.0)


def test_zero_or_negative_coefficient_raises():
    with pytest.raises(InvalidReactionError):
        Reaction("r8", reactants={"A": 0}, products={"B": 1}, rate_constant=1.0)


def test_empty_reaction_raises():
    with pytest.raises(InvalidReactionError):
        Reaction("r9", reactants={}, products={}, rate_constant=1.0)


def test_involved_species():
    r = Reaction("r10", reactants={"A": 1, "B": 1}, products={"C": 1}, rate_constant=1.0)
    assert r.involved_species() == {"A", "B", "C"}


def test_reaction_manager_add_and_get():
    rm = ReactionManager()
    r = Reaction("r1", reactants={"A": 1}, products={"B": 1}, rate_constant=1.0)
    rm.add_reaction(r)
    assert rm.get_reaction("r1") is r
    assert len(rm) == 1


def test_reaction_manager_duplicate_id_raises():
    rm = ReactionManager()
    r1 = Reaction("r1", reactants={"A": 1}, products={"B": 1}, rate_constant=1.0)
    r2 = Reaction("r1", reactants={"B": 1}, products={"A": 1}, rate_constant=1.0)
    rm.add_reaction(r1)
    with pytest.raises(ValueError):
        rm.add_reaction(r2)


def test_reaction_manager_total_propensity():
    sm = make_species(A=10, B=5, C=0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("r1", reactants={"A": 1}, products={"C": 1}, rate_constant=1.0))  # 10
    rm.add_reaction(Reaction("r2", reactants={"B": 1}, products={"C": 1}, rate_constant=2.0))  # 10
    assert rm.total_propensity(sm) == pytest.approx(20.0)
    props = rm.propensities(sm)
    assert props["r1"] == pytest.approx(10.0)
    assert props["r2"] == pytest.approx(10.0)
