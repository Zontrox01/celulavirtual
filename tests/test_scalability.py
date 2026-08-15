"""
tests/test_scalability.py

Fase 6 del roadmap (whitepaper, sección 8): prueba de escalabilidad.
Genera un genoma sintético de decenas/cientos de genes, lo carga con
el mismo GenomeModule/TranslationModule/DegradationModule usados para
el genoma de juguete de 4 genes — SIN NINGÚN CAMBIO DE CÓDIGO — y
perfila el rendimiento resultante (whitepaper, sección 6.3: "el
límite real es rendimiento, no arquitectura").
"""

import time

import pytest
from cell import Cell
from data_io.genome_loader import load_genome
from biology.genome import GenomeModule
from biology.translation import TranslationModule
from biology.degradation import DegradationModule
from tools.generate_synthetic_genome import generate_synthetic_genome


def build_cell_from_synthetic_genome(n_genes, tmp_dir, seed=1):
    fasta_path, annotation_path = generate_synthetic_genome(n_genes, tmp_dir, seed=seed)
    genome = load_genome(fasta_path, annotation_path)

    cell = Cell(seed=seed)
    genome_module = GenomeModule(genome)
    installed_genes = genome_module.install(cell, rnap_initial_count=30)

    translation_module = TranslationModule(installed_genes)
    translated_genes = translation_module.install(cell, ribosome_initial_count=30)

    degradation_module = DegradationModule(installed_genes, translated_genes)
    degradation_module.install(cell, mrna_half_life=300.0, protein_half_life=3600.0)

    return cell, genome


def test_generate_synthetic_genome_produces_valid_files(tmp_path):
    fasta_path, annotation_path = generate_synthetic_genome(10, tmp_path, seed=1)
    assert fasta_path.exists()
    assert annotation_path.exists()

    genome = load_genome(fasta_path, annotation_path)
    assert len(genome.gene_ids()) == 10


def test_synthetic_genome_has_mix_of_strands():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        fasta_path, annotation_path = generate_synthetic_genome(10, tmp_dir, seed=1)
        genome = load_genome(fasta_path, annotation_path)
        strands = {gene.strand for gene in genome.genes}
        assert strands == {"+", "-"}


def test_invalid_gene_count_raises():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        with pytest.raises(ValueError):
            generate_synthetic_genome(0, tmp_dir)


def test_no_code_changes_needed_to_scale_from_4_to_100_genes(tmp_path):
    """
    Criterio central de la Fase 6: el mismo GenomeModule/TranslationModule/
    DegradationModule que ya usamos para el genoma de juguete de 4 genes
    debe funcionar igual con 100, sin ninguna modificación de código.
    """
    cell, genome = build_cell_from_synthetic_genome(100, tmp_path, seed=42)

    assert len(genome.gene_ids()) == 100
    # RNAP + Ribosome + 100 genes + 100 mRNA + 100 proteínas
    assert len(cell.species) == 2 + 100 + 100 + 100
    # 100 transcripción + 100 traducción + 100 degradación de ARNm + 100 degradación de proteína
    assert len(cell.reactions) == 400

    events = cell.run(max_steps=50)
    assert len(events) == 50


def test_promoter_strengths_still_work_at_scale(tmp_path):
    """El mecanismo de promoter_strengths (sección 7) debe seguir funcionando igual a 50 genes."""
    fasta_path, annotation_path = generate_synthetic_genome(50, tmp_path, seed=3)
    genome = load_genome(fasta_path, annotation_path)

    cell = Cell(seed=3)
    genome_module = GenomeModule(genome)
    installed = genome_module.install(
        cell,
        promoter_strengths={"strong_promoter": 5.0, "medium_promoter": 2.0, "weak_promoter": 0.5},
    )

    strong_gene = next(g for g in genome.genes if g.promoter_id == "strong_promoter")
    weak_gene = next(g for g in genome.genes if g.promoter_id == "weak_promoter")

    strong_installed = next(g for g in installed if g.gene_id == strong_gene.gene_id)
    weak_installed = next(g for g in installed if g.gene_id == weak_gene.gene_id)

    strong_rate = cell.reactions.get_reaction(strong_installed.transcription_reaction_id).rate_constant
    weak_rate = cell.reactions.get_reaction(weak_installed.transcription_reaction_id).rate_constant

    assert strong_rate > weak_rate


def test_scalability_profiling_across_genome_sizes(tmp_path):
    """
    Perfila el tiempo por paso de simulación a distintos tamaños de
    genoma. No es un test de regresión de rendimiento estricto (el
    hardware varía) — imprime los resultados para inspección directa
    en la máquina de quien lo ejecute.
    """
    sizes = [10, 50, 100]
    steps_per_size = 500
    results = {}

    for n_genes in sizes:
        cell, _ = build_cell_from_synthetic_genome(n_genes, tmp_path / f"n{n_genes}", seed=1)

        start = time.perf_counter()
        cell.run(max_steps=steps_per_size)
        elapsed = time.perf_counter() - start

        results[n_genes] = elapsed
        print(
            f"   -> {n_genes:4d} genes | {len(cell.reactions):4d} reacciones | "
            f"{steps_per_size} pasos en {elapsed:.4f}s | "
            f"{elapsed / steps_per_size * 1000:.4f} ms/paso"
        )

    # sanity check mínimo: más genes no debería ser sistemáticamente MÁS
    # RÁPIDO (una tolerancia amplia, por el ruido de medición en un solo run)
    assert results[100] >= results[10] * 0.3
