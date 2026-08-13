"""
tests/test_regulation.py

Validación de biology/regulation.py.

Incluye el criterio de la Fase 2 del roadmap (whitepaper, sección 8):
"dinámica de un represor afectando otro gen" — se comprueba de forma
comparativa (con represión activa vs. sin ella), no solo que el código
no falle.
"""

import pytest
from cell import Cell
from data_io.genome_loader import GenomeData, GeneAnnotation
from biology.genome import GenomeModule
from biology.translation import TranslationModule
from biology.regulation import RegulationModule, RegulationModuleError


def make_genome_with_target_and_repressor():
    return GenomeData(
        sequence_id="toy",
        sequence="A" * 300,
        genes=[
            GeneAnnotation("target", 1, 90, "+", promoter_id="strong"),
            GeneAnnotation("repGene", 100, 190, "+", promoter_id="strong"),
        ],
    )


def build_cell_with_target_and_repressor(seed=1):
    cell = Cell(seed=seed)
    genome_module = GenomeModule(make_genome_with_target_and_repressor())
    installed_genes = genome_module.install(cell, rnap_initial_count=10)

    translation_module = TranslationModule(installed_genes)
    translated_genes = translation_module.install(cell, ribosome_initial_count=10)

    installed_by_id = {g.gene_id: g for g in installed_genes}
    translated_by_id = {g.gene_id: g for g in translated_genes}
    return cell, installed_by_id, translated_by_id


def test_add_repression_creates_repressed_species():
    cell, installed, translated = build_cell_with_target_and_repressor()
    regulation_module = RegulationModule()

    link = regulation_module.add_repression(cell, installed["target"], translated["repGene"])

    assert cell.species.has_species(link.repressed_gene_species_id)
    assert cell.species.get_count(link.repressed_gene_species_id) == 0


def test_binding_and_unbinding_reactions_have_correct_stoichiometry():
    cell, installed, translated = build_cell_with_target_and_repressor()
    regulation_module = RegulationModule()

    link = regulation_module.add_repression(cell, installed["target"], translated["repGene"])

    binding = cell.reactions.get_reaction(link.binding_reaction_id)
    unbinding = cell.reactions.get_reaction(link.unbinding_reaction_id)

    assert binding.reactants == {"gene_target": 1, "Prot_repGene": 1}
    assert binding.products == {link.repressed_gene_species_id: 1}
    assert unbinding.reactants == {link.repressed_gene_species_id: 1}
    assert unbinding.products == {"gene_target": 1, "Prot_repGene": 1}


def test_missing_target_gene_species_raises():
    cell, installed, translated = build_cell_with_target_and_repressor()
    regulation_module = RegulationModule()

    from biology.genome import InstalledGene
    fake_target = InstalledGene("fantasma", "gene_fantasma", "mRNA_fantasma", "transcribe_fantasma")

    with pytest.raises(RegulationModuleError):
        regulation_module.add_repression(cell, fake_target, translated["repGene"])


def test_missing_repressor_species_raises():
    cell, installed, translated = build_cell_with_target_and_repressor()
    regulation_module = RegulationModule()

    from biology.translation import TranslatedGene
    fake_repressor = TranslatedGene("fantasma", "mRNA_fantasma", "Prot_fantasma", "translate_fantasma")

    with pytest.raises(RegulationModuleError):
        regulation_module.add_repression(cell, installed["target"], fake_repressor)


def test_duplicate_repression_raises():
    cell, installed, translated = build_cell_with_target_and_repressor()
    regulation_module = RegulationModule()
    regulation_module.add_repression(cell, installed["target"], translated["repGene"])

    with pytest.raises(RegulationModuleError):
        regulation_module.add_repression(cell, installed["target"], translated["repGene"])


def test_get_link_before_install_raises():
    regulation_module = RegulationModule()
    with pytest.raises(RegulationModuleError):
        regulation_module.get_link("target", "repGene")


def test_get_link_after_install_returns_it():
    cell, installed, translated = build_cell_with_target_and_repressor()
    regulation_module = RegulationModule()
    link = regulation_module.add_repression(cell, installed["target"], translated["repGene"])

    assert regulation_module.get_link("target", "repGene") is link


def test_gene_copy_toggles_between_free_and_repressed():
    """
    Al ejecutar solo la reacción de unión, la copia libre del gen debe
    desaparecer (pasar a 0) y aparecer como reprimida; al ejecutar la
    de disociación, debe volver.
    """
    cell, installed, translated = build_cell_with_target_and_repressor()
    regulation_module = RegulationModule()
    link = regulation_module.add_repression(cell, installed["target"], translated["repGene"])

    # todavía no se ha corrido la simulación, así que no hay represor
    # traducido: se le da manualmente una copia para poder forzar la unión
    cell.species.change_count("Prot_repGene", 1)

    # forzar manualmente la reacción de unión (sin pasar por el motor SSA)
    binding = cell.reactions.get_reaction(link.binding_reaction_id)
    binding.apply(cell.species)

    assert cell.species.get_count("gene_target") == 0
    assert cell.species.get_count(link.repressed_gene_species_id) == 1

    unbinding = cell.reactions.get_reaction(link.unbinding_reaction_id)
    unbinding.apply(cell.species)

    assert cell.species.get_count("gene_target") == 1
    assert cell.species.get_count(link.repressed_gene_species_id) == 0


# ---------------------------------------------------------------------
# Criterio central de la Fase 2: represión activa reduce la expresión
# ---------------------------------------------------------------------

def test_active_repression_reduces_target_mrna_compared_to_no_repression():
    """
    Con represión fuerte (unión rápida, disociación lenta) y una buena
    cantidad inicial de represor ya traducido, la producción de ARNm
    del gen diana debe ser sustancialmente menor que sin represión.
    """
    # --- Escenario A: sin represión ---
    cell_a, installed_a, translated_a = build_cell_with_target_and_repressor(seed=100)
    cell_a.run(max_time=200.0)
    mrna_without_repression = cell_a.species.get_count("mRNA_target")

    # --- Escenario B: con represión fuerte ---
    cell_b, installed_b, translated_b = build_cell_with_target_and_repressor(seed=100)
    regulation_module = RegulationModule()
    regulation_module.add_repression(
        cell_b,
        installed_b["target"],
        translated_b["repGene"],
        binding_rate=50.0,    # unión muy rápida
        unbinding_rate=0.01,  # disociación muy lenta -> el gen pasa la mayor parte reprimido
    )
    cell_b.run(max_time=200.0)
    mrna_with_repression = cell_b.species.get_count("mRNA_target")

    assert mrna_with_repression < mrna_without_repression / 2, (
        f"con represión ({mrna_with_repression}) debería ser bastante menor que "
        f"sin represión ({mrna_without_repression})"
    )


def test_repression_link_reaction_ids_actually_fire_during_simulation():
    cell, installed, translated = build_cell_with_target_and_repressor(seed=5)
    regulation_module = RegulationModule()
    link = regulation_module.add_repression(
        cell, installed["target"], translated["repGene"], binding_rate=5.0, unbinding_rate=5.0
    )

    cell.run(max_time=50.0)

    used_reactions = {e.reaction_id for e in cell.get_events()}
    assert link.binding_reaction_id in used_reactions
    assert link.unbinding_reaction_id in used_reactions
