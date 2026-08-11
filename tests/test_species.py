"""
tests/test_species.py

Validación de engine/species.py.
"""

import pytest
from engine.species import SpeciesManager, UnknownSpeciesError, NegativeCountError


def test_add_and_get_species():
    sm = SpeciesManager()
    sm.add_species("ATP", initial_count=100)
    assert sm.get_count("ATP") == 100
    assert len(sm) == 1


def test_add_species_default_zero():
    sm = SpeciesManager()
    sm.add_species("mRNA_geneA")
    assert sm.get_count("mRNA_geneA") == 0


def test_add_duplicate_species_raises():
    sm = SpeciesManager()
    sm.add_species("ATP", 10)
    with pytest.raises(ValueError):
        sm.add_species("ATP", 5)


def test_add_species_negative_initial_raises():
    sm = SpeciesManager()
    with pytest.raises(NegativeCountError):
        sm.add_species("ATP", -1)


def test_get_unknown_species_raises():
    sm = SpeciesManager()
    with pytest.raises(UnknownSpeciesError):
        sm.get_count("no_existe")


def test_change_count_increments_and_decrements():
    sm = SpeciesManager()
    sm.add_species("ATP", 100)
    sm.change_count("ATP", -30)
    assert sm.get_count("ATP") == 70
    sm.change_count("ATP", 5)
    assert sm.get_count("ATP") == 75


def test_change_count_below_zero_raises_and_does_not_mutate():
    sm = SpeciesManager()
    sm.add_species("ATP", 10)
    with pytest.raises(NegativeCountError):
        sm.change_count("ATP", -20)
    # el estado no debe haber cambiado tras el intento fallido
    assert sm.get_count("ATP") == 10


def test_change_count_unknown_species_raises():
    sm = SpeciesManager()
    with pytest.raises(UnknownSpeciesError):
        sm.change_count("no_existe", 1)


def test_set_count_overwrites():
    sm = SpeciesManager()
    sm.add_species("ATP", 10)
    sm.set_count("ATP", 999)
    assert sm.get_count("ATP") == 999


def test_has_species():
    sm = SpeciesManager()
    sm.add_species("ATP", 10)
    assert sm.has_species("ATP") is True
    assert sm.has_species("GTP") is False


def test_species_ids():
    sm = SpeciesManager()
    sm.add_species("ATP", 10)
    sm.add_species("GTP", 5)
    assert set(sm.species_ids()) == {"ATP", "GTP"}


def test_snapshot_and_restore():
    sm = SpeciesManager()
    sm.add_species("ATP", 100)
    sm.add_species("GTP", 50)

    snap = sm.snapshot()
    sm.change_count("ATP", -100)
    assert sm.get_count("ATP") == 0

    sm.restore(snap)
    assert sm.get_count("ATP") == 100
    assert sm.get_count("GTP") == 50


def test_snapshot_is_independent_copy():
    sm = SpeciesManager()
    sm.add_species("ATP", 100)
    snap = sm.snapshot()
    sm.change_count("ATP", -10)
    # la snapshot tomada antes no debe verse afectada por el cambio posterior
    assert snap["ATP"] == 100
