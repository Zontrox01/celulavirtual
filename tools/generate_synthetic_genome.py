"""
tools/generate_synthetic_genome.py

Genera un genoma sintético de N genes (FASTA + anotación CSV), del
mismo tipo que examples/genomas/toy_genome.fasta pero a escala
configurable. Se usa para la Fase 6 del roadmap (whitepaper, sección
8): "prueba de escalabilidad — cargar un genoma anotado real/grande y
confirmar que el GenomeModule genera las reacciones correctamente sin
cambios de código; perfilar el rendimiento resultante."

No es un módulo de la librería (no vive en engine/ ni en biology/) —
es una herramienta de generación de datos de prueba, coherente con la
sección 6.2 del whitepaper: el genoma SIEMPRE se carga desde archivo,
nunca se codifica a mano; este script produce esos archivos.

No estaba prevista en la arquitectura original del whitepaper (que no
contemplaba un directorio `tools/`); se documenta el añadido en
FILES.md siguiendo la misma regla aplicada a `biology/regulation.py`.
"""

from __future__ import annotations
import csv
import random
from pathlib import Path
from typing import Tuple, Union

DEFAULT_PROMOTER_IDS = ("strong_promoter", "medium_promoter", "weak_promoter")


def _random_orf(rng: random.Random, n_codons: int) -> str:
    """ORF válido: empieza en ATG, termina en TAA, longitud múltiplo de 3."""
    body = "".join(rng.choice("ACGT") for _ in range(n_codons * 3 - 6))
    return "ATG" + body + "TAA"


def generate_synthetic_genome(
    n_genes: int,
    output_dir: Union[str, Path],
    prefix: str = "synthetic",
    seed: int = 0,
    min_codons: int = 20,
    max_codons: int = 60,
    spacer_length: int = 15,
) -> Tuple[Path, Path]:
    """
    Genera un genoma sintético de `n_genes` genes con ORFs válidos,
    intercalados con espaciadores aleatorios, y su anotación CSV
    (gene_id, start, end, strand, promoter_id). La mitad de los genes
    van en hebra +, la otra mitad en hebra -, para ejercitar ambos
    casos del cargador (data_io/genome_loader.py).

    Devuelve (ruta_fasta, ruta_anotacion_csv).
    """
    if n_genes <= 0:
        raise ValueError(f"n_genes debe ser positivo (recibido {n_genes}).")

    rng = random.Random(seed)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = output_dir / f"{prefix}_{n_genes}genes.fasta"
    annotation_path = output_dir / f"{prefix}_{n_genes}genes_annotation.csv"

    sequence = ""
    genes = []  # (gene_id, start, end, strand, promoter_id)

    for i in range(n_genes):
        gene_id = f"gene{i:04d}"
        n_codons = rng.randint(min_codons, max_codons)
        orf = _random_orf(rng, n_codons)
        strand = "+" if i % 2 == 0 else "-"
        promoter_id = DEFAULT_PROMOTER_IDS[i % len(DEFAULT_PROMOTER_IDS)]

        sequence += "".join(rng.choice("ACGT") for _ in range(spacer_length))
        start = len(sequence) + 1
        sequence += orf
        end = len(sequence)

        genes.append((gene_id, start, end, strand, promoter_id))

    sequence += "".join(rng.choice("ACGT") for _ in range(spacer_length))

    with open(fasta_path, "w") as f:
        f.write(f">{prefix}_{n_genes}genes genoma sintético para pruebas de escalabilidad\n")
        for i in range(0, len(sequence), 70):
            f.write(sequence[i : i + 70] + "\n")

    with open(annotation_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id", "start", "end", "strand", "promoter_id"])
        for gene_id, start, end, strand, promoter_id in genes:
            writer.writerow([gene_id, start, end, strand, promoter_id])

    return fasta_path, annotation_path
