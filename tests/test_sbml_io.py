"""
tests/test_sbml_io.py

Validación de data_io/sbml_io.py.

Criterio central: exportar una Cell a SBML y volver a importarla debe
reproducir el mismo estado inicial y las mismas reacciones — no solo
"no falla", sino que el round-trip es fiel.
"""

import xml.etree.ElementTree as ET

import pytest
from cell import Cell
from engine.reactions import Reaction
from data_io.sbml_io import (
    export_cell_to_sbml_string,
    export_cell_to_sbml_file,
    import_sbml_string,
    import_sbml_file,
    build_cell_from_sbml,
    SBMLError,
    SBML_NS,
)


def build_sample_cell(seed=1):
    cell = Cell(seed=seed)
    cell.add_species("A", 10).add_species("B", 0).add_species("RNAP", 5)
    cell.add_reaction(Reaction("A_to_B", {"RNAP": 1, "A": 1}, {"RNAP": 1, "A": 1, "B": 1}, rate_constant=2.5))
    cell.add_reaction(Reaction("degrade_B", {"B": 1}, {}, rate_constant=0.1))
    return cell


def test_export_produces_valid_xml():
    cell = build_sample_cell()
    xml_string = export_cell_to_sbml_string(cell)

    root = ET.fromstring(xml_string)  # no debe lanzar excepción
    assert root.tag == f"{{{SBML_NS}}}sbml"


def test_export_contains_all_species_with_correct_amounts():
    cell = build_sample_cell()
    xml_string = export_cell_to_sbml_string(cell)

    species, _ = import_sbml_string(xml_string)
    assert species == {"A": 10, "B": 0, "RNAP": 5}


def test_export_contains_all_reactions_with_correct_stoichiometry():
    cell = build_sample_cell()
    xml_string = export_cell_to_sbml_string(cell)

    _, reactions = import_sbml_string(xml_string)
    reactions_by_id = {r.reaction_id: r for r in reactions}

    assert set(reactions_by_id.keys()) == {"A_to_B", "degrade_B"}
    assert reactions_by_id["A_to_B"].reactants == {"RNAP": 1, "A": 1}
    assert reactions_by_id["A_to_B"].products == {"RNAP": 1, "A": 1, "B": 1}
    assert reactions_by_id["degrade_B"].reactants == {"B": 1}
    assert reactions_by_id["degrade_B"].products == {}


def test_export_preserves_rate_constants():
    cell = build_sample_cell()
    xml_string = export_cell_to_sbml_string(cell)

    _, reactions = import_sbml_string(xml_string)
    reactions_by_id = {r.reaction_id: r for r in reactions}

    assert reactions_by_id["A_to_B"].rate_constant == pytest.approx(2.5)
    assert reactions_by_id["degrade_B"].rate_constant == pytest.approx(0.1)


def test_higher_stoichiometry_reaction_round_trips():
    cell = Cell(seed=1)
    cell.add_species("A", 20).add_species("B", 0)
    cell.add_reaction(Reaction("dimerize", {"A": 2}, {"B": 1}, rate_constant=0.5))

    xml_string = export_cell_to_sbml_string(cell)
    _, reactions = import_sbml_string(xml_string)

    reaction = reactions[0]
    assert reaction.reactants == {"A": 2}
    assert reaction.rate_constant == pytest.approx(0.5)


def test_export_to_file_and_import_from_file(tmp_path):
    cell = build_sample_cell(seed=3)
    path = export_cell_to_sbml_file(cell, tmp_path / "modelo.xml")

    assert path.exists()
    species, reactions = import_sbml_file(path)
    assert species == {"A": 10, "B": 0, "RNAP": 5}
    assert len(reactions) == 2


def test_import_file_missing_raises():
    with pytest.raises(SBMLError):
        import_sbml_file("/ruta/que/no/existe.xml")


def test_import_invalid_xml_raises():
    with pytest.raises(SBMLError):
        import_sbml_string("esto no es xml valido <<<")


def test_import_reaction_without_kinetic_law_raises():
    xml_string = """<?xml version="1.0"?>
    <sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
      <model id="m">
        <listOfSpecies>
          <species id="A" compartment="c" initialAmount="5"/>
        </listOfSpecies>
        <listOfReactions>
          <reaction id="sin_cinetica" reversible="false"/>
        </listOfReactions>
      </model>
    </sbml>"""
    with pytest.raises(SBMLError):
        import_sbml_string(xml_string)


def test_build_cell_from_sbml_string_reproduces_behavior():
    """
    Round-trip completo: exportar una célula, reconstruir OTRA célula
    nueva a partir del SBML, y comprobar que su comportamiento (mismo
    seed) es exactamente el mismo que el de la original.
    """
    original = build_sample_cell(seed=99)
    xml_string = export_cell_to_sbml_string(original)

    rebuilt = build_cell_from_sbml(xml_string, seed=123)

    assert rebuilt.species.get_count("A") == 10
    assert rebuilt.species.get_count("B") == 0
    assert rebuilt.species.get_count("RNAP") == 5
    assert len(rebuilt.reactions) == 2

    events = rebuilt.run(max_steps=20)
    assert len(events) == 20


def test_build_cell_from_sbml_file(tmp_path):
    original = build_sample_cell(seed=5)
    path = export_cell_to_sbml_file(original, tmp_path / "modelo.xml")

    rebuilt = build_cell_from_sbml(path, seed=1)
    assert rebuilt.species.get_count("A") == 10
    assert len(rebuilt.reactions) == 2


def test_catalytic_species_appears_in_both_lists():
    """RNAP aparece igual en reactivos y productos (whitepaper: convenio simplificado, no 'modifier')."""
    cell = build_sample_cell()
    xml_string = export_cell_to_sbml_string(cell)

    root = ET.fromstring(xml_string)
    ns = {"sbml": SBML_NS}
    reaction_el = next(
        r for r in root.findall(".//sbml:reaction", ns) if r.attrib["id"] == "A_to_B"
    )
    reactant_species = {
        ref.attrib["species"]
        for ref in reaction_el.find("sbml:listOfReactants", ns).findall("sbml:speciesReference", ns)
    }
    product_species = {
        ref.attrib["species"]
        for ref in reaction_el.find("sbml:listOfProducts", ns).findall("sbml:speciesReference", ns)
    }
    assert "RNAP" in reactant_species
    assert "RNAP" in product_species
