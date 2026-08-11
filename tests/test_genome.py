"""
tests/test_genome.py

Validación de biology/genome.py.
"""

import pytest
from cell import Cell
from data_io.genome_loader import GenomeData, GeneAnnotation, load_genome
from biology.genome import GenomeModule, GenomeModuleError


def make_toy_genome_data():
    """GenomeData construido directamente en memoria (sin pasar por FASTA/CSV)."""
    return GenomeData(
        sequence_id="toy",
        sequence="A" * 300,  # contenido irrelevante para este módulo, solo importa la anotación
        genes=[
            GeneAnnotation("ribo1", 1, 90, "+", promoter_id="strong"),
            GeneAnnotation("ribo2", 100, 190, "+", promoter_id="strong"),
            GeneAnnotation("metA", 200, 290, "+", promoter_id="medium"),
        ],
    )


def test_install_creates_gene_and_mrna_species_for_each_gene():
    cell = Cell(seed=1)
    module = GenomeModule(make_toy_genome_data())

    module.install(cell)

    for gene_id in ("ribo1", "ribo2", "metA"):
        assert cell.species.has_species(f"gene_{gene_id}")
        assert cell.species.get_count(f"gene_{gene_id}") == 1
        assert cell.species.has_species(f"mRNA_{gene_id}")
        assert cell.species.get_count(f"mRNA_{gene_id}") == 0


def test_install_creates_rnap_species_with_default_count():
    cell = Cell(seed=1)
    module = GenomeModule(make_toy_genome_data())

    module.install(cell)

    assert cell.species.has_species("RNAP")
    assert cell.species.get_count("RNAP") == 30


def test_install_does_not_duplicate_rnap_if_already_present():
    cell = Cell(seed=1)
    cell.add_species("RNAP", 999)
    module = GenomeModule(make_toy_genome_data())

    module.install(cell)

    assert cell.species.get_count("RNAP") == 999  # no se sobrescribe


def test_install_creates_one_transcription_reaction_per_gene():
    cell = Cell(seed=1)
    module = GenomeModule(make_toy_genome_data())

    installed = module.install(cell)

    assert len(installed) == 3
    assert len(cell.reactions) == 3
    for record in installed:
        reaction = cell.reactions.get_reaction(record.transcription_reaction_id)
        assert reaction.reactants == {"RNAP": 1, record.gene_species_id: 1}
        assert reaction.products == {
            "RNAP": 1,
            record.gene_species_id: 1,
            record.mrna_species_id: 1,
        }


def test_promoter_strengths_scale_the_transcription_rate():
    cell = Cell(seed=1)
    module = GenomeModule(make_toy_genome_data())

    module.install(
        cell,
        base_transcription_rate=2.0,
        promoter_strengths={"strong": 3.0, "medium": 1.0},
    )

    ribo1 = module.get_installed("ribo1")
    metA = module.get_installed("metA")

    ribo1_reaction = cell.reactions.get_reaction(ribo1.transcription_reaction_id)
    metA_reaction = cell.reactions.get_reaction(metA.transcription_reaction_id)

    assert ribo1_reaction.rate_constant == pytest.approx(2.0 * 3.0)
    assert metA_reaction.rate_constant == pytest.approx(2.0 * 1.0)


def test_unlisted_promoter_defaults_to_base_rate():
    cell = Cell(seed=1)
    genome = GenomeData(
        sequence_id="toy",
        sequence="A" * 50,
        genes=[GeneAnnotation("sin_promotor_conocido", 1, 30, "+")],
    )
    module = GenomeModule(genome)

    module.install(cell, base_transcription_rate=5.0, promoter_strengths={"otro": 10.0})

    record = module.get_installed("sin_promotor_conocido")
    reaction = cell.reactions.get_reaction(record.transcription_reaction_id)
    assert reaction.rate_constant == pytest.approx(5.0)


def test_transcription_produces_mrna_and_conserves_gene_copy():
    cell = Cell(seed=2)
    module = GenomeModule(make_toy_genome_data())
    module.install(cell, rnap_initial_count=5)

    cell.run(max_steps=20)

    # el gen nunca se consume: siempre debe quedar exactamente 1 copia
    assert cell.species.get_count("gene_ribo1") == 1
    assert cell.species.get_count("gene_ribo2") == 1
    assert cell.species.get_count("gene_metA") == 1
    # y debe haberse producido algo de ARNm en total
    total_mrna = (
        cell.species.get_count("mRNA_ribo1")
        + cell.species.get_count("mRNA_ribo2")
        + cell.species.get_count("mRNA_metA")
    )
    assert total_mrna == 20  # cada paso ejecutado es, por diseño, una transcripción


def test_installing_the_same_genome_twice_raises_name_conflict():
    cell = Cell(seed=1)
    module = GenomeModule(make_toy_genome_data())
    module.install(cell)

    with pytest.raises(GenomeModuleError):
        module.install(cell)


def test_get_installed_before_install_raises():
    module = GenomeModule(make_toy_genome_data())
    with pytest.raises(GenomeModuleError):
        module.get_installed("ribo1")


# ---------------------------------------------------------------------
# Integración de extremo a extremo con el genoma de juguete real en disco
# ---------------------------------------------------------------------

def test_end_to_end_with_real_toy_genome_files():
    genome = load_genome(
        "examples/genomas/toy_genome.fasta",
        "examples/genomas/toy_genome_annotation.csv",
    )
    assert set(genome.gene_ids()) == {"ribo1", "ribo2", "metA", "reg1"}

    cell = Cell(seed=3)
    module = GenomeModule(genome)
    installed = module.install(
        cell,
        promoter_strengths={
            "strong_promoter": 3.0,
            "medium_promoter": 1.5,
            "weak_promoter": 0.5,
        },
    )

    assert len(installed) == 4
    assert len(cell.reactions) == 4

    events = cell.run(max_steps=50)
    assert len(events) == 50
    # las 4 copias de gen deben seguir presentes: la transcripción no las consume
    for gene_id in ("ribo1", "ribo2", "metA", "reg1"):
        assert cell.species.get_count(f"gene_{gene_id}") == 1
