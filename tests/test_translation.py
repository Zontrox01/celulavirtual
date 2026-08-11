"""
tests/test_translation.py

Validación de biology/translation.py.
"""

import pytest
from cell import Cell
from data_io.genome_loader import GenomeData, GeneAnnotation
from biology.genome import GenomeModule, InstalledGene
from biology.translation import TranslationModule, TranslationModuleError


def make_toy_genome_data():
    return GenomeData(
        sequence_id="toy",
        sequence="A" * 300,
        genes=[
            GeneAnnotation("ribo1", 1, 90, "+", promoter_id="strong"),
            GeneAnnotation("metA", 100, 190, "+", promoter_id="medium"),
        ],
    )


def build_cell_with_genome_installed(seed=1):
    cell = Cell(seed=seed)
    genome_module = GenomeModule(make_toy_genome_data())
    installed_genes = genome_module.install(cell)
    return cell, installed_genes


def test_install_creates_protein_species_for_each_gene():
    cell, installed_genes = build_cell_with_genome_installed()
    translation_module = TranslationModule(installed_genes)

    translation_module.install(cell)

    assert cell.species.has_species("Prot_ribo1")
    assert cell.species.get_count("Prot_ribo1") == 0
    assert cell.species.has_species("Prot_metA")
    assert cell.species.get_count("Prot_metA") == 0


def test_install_creates_ribosome_species_with_default_count():
    cell, installed_genes = build_cell_with_genome_installed()
    translation_module = TranslationModule(installed_genes)

    translation_module.install(cell)

    assert cell.species.has_species("Ribosome")
    assert cell.species.get_count("Ribosome") == 20


def test_install_does_not_duplicate_ribosome_if_already_present():
    cell, installed_genes = build_cell_with_genome_installed()
    cell.add_species("Ribosome", 777)
    translation_module = TranslationModule(installed_genes)

    translation_module.install(cell)

    assert cell.species.get_count("Ribosome") == 777


def test_install_creates_one_translation_reaction_per_gene():
    cell, installed_genes = build_cell_with_genome_installed()
    translation_module = TranslationModule(installed_genes)

    translated = translation_module.install(cell)

    assert len(translated) == 2
    for record in translated:
        reaction = cell.reactions.get_reaction(record.translation_reaction_id)
        assert reaction.reactants == {"Ribosome": 1, record.mrna_species_id: 1}
        assert reaction.products == {
            "Ribosome": 1,
            record.mrna_species_id: 1,
            record.protein_species_id: 1,
        }


def test_rbs_strengths_scale_the_translation_rate():
    cell, installed_genes = build_cell_with_genome_installed()
    translation_module = TranslationModule(installed_genes)

    translation_module.install(
        cell, base_translation_rate=2.0, rbs_strengths={"ribo1": 5.0}
    )

    ribo1 = translation_module.get_translated("ribo1")
    metA = translation_module.get_translated("metA")

    r1 = cell.reactions.get_reaction(ribo1.translation_reaction_id)
    r2 = cell.reactions.get_reaction(metA.translation_reaction_id)

    assert r1.rate_constant == pytest.approx(2.0 * 5.0)
    assert r2.rate_constant == pytest.approx(2.0 * 1.0)  # metA sin entrada -> multiplicador 1.0


def test_translation_module_rejects_empty_installed_genes():
    with pytest.raises(TranslationModuleError):
        TranslationModule([])


def test_install_raises_if_mrna_species_missing():
    cell = Cell(seed=1)  # célula vacía, sin GenomeModule instalado
    fake_installed = [
        InstalledGene(
            gene_id="fantasma",
            gene_species_id="gene_fantasma",
            mrna_species_id="mRNA_fantasma",
            transcription_reaction_id="transcribe_fantasma",
        )
    ]
    translation_module = TranslationModule(fake_installed)

    with pytest.raises(TranslationModuleError):
        translation_module.install(cell)


def test_installing_translation_twice_raises_name_conflict():
    cell, installed_genes = build_cell_with_genome_installed()
    translation_module = TranslationModule(installed_genes)
    translation_module.install(cell)

    with pytest.raises(TranslationModuleError):
        translation_module.install(cell)


def test_get_translated_before_install_raises():
    _, installed_genes = build_cell_with_genome_installed()
    translation_module = TranslationModule(installed_genes)
    with pytest.raises(TranslationModuleError):
        translation_module.get_translated("ribo1")


# ---------------------------------------------------------------------
# Flujo completo: ADN -> ARNm -> proteína sobre la misma célula
# ---------------------------------------------------------------------

def test_full_transcription_and_translation_pipeline():
    cell = Cell(seed=42)
    genome_module = GenomeModule(make_toy_genome_data())
    installed_genes = genome_module.install(cell, rnap_initial_count=10)

    translation_module = TranslationModule(installed_genes)
    translation_module.install(cell, ribosome_initial_count=10)

    # 4 reacciones activas: 2 de transcripción + 2 de traducción
    assert len(cell.reactions) == 4

    events = cell.run(max_steps=200)
    assert len(events) == 200

    # las copias de gen nunca se consumen (catalizador)
    assert cell.species.get_count("gene_ribo1") == 1
    assert cell.species.get_count("gene_metA") == 1

    # tiene que haberse producido tanto ARNm como proteína
    assert cell.species.get_count("mRNA_ribo1") > 0
    assert cell.species.get_count("Prot_ribo1") > 0
    assert cell.species.get_count("mRNA_metA") > 0
    assert cell.species.get_count("Prot_metA") > 0

    # cada evento ejecutado corresponde a exactamente una de las 4 reacciones
    reaction_ids_used = {e.reaction_id for e in events}
    assert reaction_ids_used <= {
        "transcribe_ribo1", "transcribe_metA",
        "translate_ribo1", "translate_metA",
    }
