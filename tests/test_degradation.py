"""
tests/test_degradation.py

Validación de biology/degradation.py.
"""

import math
import pytest
from cell import Cell
from data_io.genome_loader import GenomeData, GeneAnnotation
from biology.genome import GenomeModule
from biology.translation import TranslationModule
from biology.degradation import DegradationModule, DegradationModuleError


def make_toy_genome_data():
    return GenomeData(
        sequence_id="toy",
        sequence="A" * 300,
        genes=[
            GeneAnnotation("ribo1", 1, 90, "+", promoter_id="strong"),
            GeneAnnotation("metA", 100, 190, "+", promoter_id="medium"),
        ],
    )


def build_cell_with_genome_and_translation(seed=1):
    cell = Cell(seed=seed)
    genome_module = GenomeModule(make_toy_genome_data())
    installed_genes = genome_module.install(cell)
    translation_module = TranslationModule(installed_genes)
    translated_genes = translation_module.install(cell)
    return cell, installed_genes, translated_genes


def test_install_mrna_only_when_no_translated_genes_given():
    cell = Cell(seed=1)
    genome_module = GenomeModule(make_toy_genome_data())
    installed_genes = genome_module.install(cell)

    degradation_module = DegradationModule(installed_genes)  # sin translated_genes
    records = degradation_module.install(cell)

    assert len(records) == 2  # solo los 2 ARNm, ninguna proteína
    assert all(r.species_kind == "mRNA" for r in records)


def test_install_mrna_and_protein_when_translated_genes_given():
    cell, installed_genes, translated_genes = build_cell_with_genome_and_translation()

    degradation_module = DegradationModule(installed_genes, translated_genes)
    records = degradation_module.install(cell)

    assert len(records) == 4  # 2 ARNm + 2 proteínas
    kinds = {r.species_kind for r in records}
    assert kinds == {"mRNA", "protein"}


def test_rate_constant_derived_correctly_from_half_life():
    cell, installed_genes, translated_genes = build_cell_with_genome_and_translation()
    degradation_module = DegradationModule(installed_genes, translated_genes)

    degradation_module.install(cell, mrna_half_life=100.0, protein_half_life=1000.0)

    mrna_record = degradation_module.get_degradation("mRNA_ribo1")
    protein_record = degradation_module.get_degradation("Prot_ribo1")

    mrna_reaction = cell.reactions.get_reaction(mrna_record.reaction_id)
    protein_reaction = cell.reactions.get_reaction(protein_record.reaction_id)

    assert mrna_reaction.rate_constant == pytest.approx(math.log(2) / 100.0)
    assert protein_reaction.rate_constant == pytest.approx(math.log(2) / 1000.0)


def test_per_gene_half_life_overrides_default():
    cell, installed_genes, translated_genes = build_cell_with_genome_and_translation()
    degradation_module = DegradationModule(installed_genes, translated_genes)

    degradation_module.install(
        cell,
        mrna_half_life=300.0,
        mrna_half_lives={"ribo1": 30.0},  # ARNm de ribo1 se degrada mucho más rápido
    )

    ribo1_reaction = cell.reactions.get_reaction(
        degradation_module.get_degradation("mRNA_ribo1").reaction_id
    )
    metA_reaction = cell.reactions.get_reaction(
        degradation_module.get_degradation("mRNA_metA").reaction_id
    )

    assert ribo1_reaction.rate_constant == pytest.approx(math.log(2) / 30.0)
    assert metA_reaction.rate_constant == pytest.approx(math.log(2) / 300.0)  # usa el default


def test_zero_or_negative_half_life_raises():
    cell, installed_genes, _ = build_cell_with_genome_and_translation()
    degradation_module = DegradationModule(installed_genes)

    with pytest.raises(DegradationModuleError):
        degradation_module.install(cell, mrna_half_life=0.0)


def test_reaction_has_no_products_species_is_simply_removed():
    cell, installed_genes, _ = build_cell_with_genome_and_translation()
    degradation_module = DegradationModule(installed_genes)
    degradation_module.install(cell)

    reaction = cell.reactions.get_reaction(
        degradation_module.get_degradation("mRNA_ribo1").reaction_id
    )
    assert reaction.reactants == {"mRNA_ribo1": 1}
    assert reaction.products == {}


def test_install_raises_if_mrna_species_missing():
    cell = Cell(seed=1)  # célula vacía
    from biology.genome import InstalledGene

    fake_installed = [
        InstalledGene("fantasma", "gene_fantasma", "mRNA_fantasma", "transcribe_fantasma")
    ]
    degradation_module = DegradationModule(fake_installed)

    with pytest.raises(DegradationModuleError):
        degradation_module.install(cell)


def test_degradation_module_rejects_empty_installed_genes():
    with pytest.raises(DegradationModuleError):
        DegradationModule([])


def test_get_degradation_before_install_raises():
    cell, installed_genes, _ = build_cell_with_genome_and_translation()
    degradation_module = DegradationModule(installed_genes)
    with pytest.raises(DegradationModuleError):
        degradation_module.get_degradation("mRNA_ribo1")


# ---------------------------------------------------------------------
# Pipeline completo: transcripción + traducción + degradación
# ---------------------------------------------------------------------

def test_full_pipeline_reaches_a_bounded_steady_state():
    """
    Con producción Y degradación activas, las cantidades de ARNm y
    proteína deben estabilizarse en vez de crecer sin límite. Se
    comprueba dejando correr la simulación durante un tramo de tiempo
    adicional (usando la reanudación real tras max_time) y verificando
    que las cantidades no se disparan proporcionalmente al tiempo
    extra — señal de que se ha alcanzado un equilibrio dinámico entre
    producción y degradación, no que el sistema sigue creciendo.
    """
    cell = Cell(seed=7)
    genome_module = GenomeModule(make_toy_genome_data())
    installed_genes = genome_module.install(cell, rnap_initial_count=10)

    translation_module = TranslationModule(installed_genes)
    translated_genes = translation_module.install(cell, ribosome_initial_count=10)

    degradation_module = DegradationModule(installed_genes, translated_genes)
    degradation_module.install(cell, mrna_half_life=5.0, protein_half_life=20.0)

    # 4 de producción (2 transcripción + 2 traducción) + 4 de degradación
    assert len(cell.reactions) == 8

    # tramo inicial (varias vidas medias de proteína, para llegar a régimen estacionario)
    cell.run(max_time=150.0)
    mrna_mid = cell.species.get_count("mRNA_ribo1")
    protein_mid = cell.species.get_count("Prot_ribo1")
    assert mrna_mid > 0
    assert protein_mid > 0

    # duplicar el tiempo simulado total (150 -> 300)
    cell.run(max_time=300.0)
    mrna_late = cell.species.get_count("mRNA_ribo1")
    protein_late = cell.species.get_count("Prot_ribo1")

    # si no hubiera degradación, duplicar el tiempo dispararía las
    # cantidades de forma aproximadamente proporcional; en equilibrio
    # deben quedarse en el mismo orden de magnitud
    assert mrna_late < mrna_mid * 3
    assert protein_late < protein_mid * 3
    assert mrna_late > 0
    assert protein_late > 0

    # y debe haber ocurrido degradación de verdad, no solo producción
    degradation_reaction_ids = {r.reaction_id for r in degradation_module._records.values()}
    used_reactions = {e.reaction_id for e in cell.get_events()}
    assert degradation_reaction_ids & used_reactions