"""
data_io/genome_loader.py

Carga de genoma desde archivos externos (FASTA + anotación CSV).

Responsabilidad única: leer la secuencia de ADN (vía BioPython) y su
anotación de genes (posición, hebra, promotor) desde archivos, y
exponerlos como una estructura de datos simple (GenomeData) que
biology/genome.py (Fase 1) usará para generar las reacciones de
transcripción. El genoma NUNCA se codifica a mano en Python — este
módulo es el único punto de entrada de datos genómicos externos
(whitepaper, sección 6.2).
"""

from __future__ import annotations
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


class GenomeLoadError(ValueError):
    """Se lanza si el FASTA o la anotación están mal formados o son inconsistentes entre sí."""


@dataclass(frozen=True)
class GeneAnnotation:
    """
    Anotación de un gen sobre la secuencia del genoma.

    Coordenadas 1-based, inclusivas en ambos extremos (convención
    habitual en anotaciones biológicas tipo GFF3).
    """

    gene_id: str
    start: int
    end: int
    strand: str  # '+' o '-'
    promoter_id: Optional[str] = None

    def __post_init__(self) -> None:
        if self.strand not in ("+", "-"):
            raise GenomeLoadError(
                f"Gen '{self.gene_id}': la hebra debe ser '+' o '-' (recibido '{self.strand}')."
            )
        if self.start < 1 or self.end < self.start:
            raise GenomeLoadError(
                f"Gen '{self.gene_id}': coordenadas inválidas (start={self.start}, end={self.end})."
            )

    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class GenomeData:
    """Genoma cargado: secuencia completa + lista de genes anotados."""

    sequence_id: str
    sequence: str
    genes: List[GeneAnnotation] = field(default_factory=list)

    def get_gene_sequence(self, gene_id: str) -> str:
        """
        Devuelve la secuencia del gen, ya orientada según su hebra
        (complemento reverso real vía BioPython si strand == '-').
        """
        gene = self._find_gene(gene_id)
        # coordenadas 1-based inclusivas -> slicing de Python 0-based, exclusivo al final
        raw = self.sequence[gene.start - 1 : gene.end]
        if gene.strand == "-":
            return str(Seq(raw).reverse_complement())
        return raw

    def _find_gene(self, gene_id: str) -> GeneAnnotation:
        for gene in self.genes:
            if gene.gene_id == gene_id:
                return gene
        raise GenomeLoadError(f"Gen '{gene_id}' no encontrado en la anotación.")

    def gene_ids(self) -> List[str]:
        return [g.gene_id for g in self.genes]


def _load_sequence(fasta_path: Path) -> Tuple[str, str]:
    """
    Lee el FASTA y devuelve (id_de_secuencia, secuencia). El genoma de
    juguete (whitepaper, sección 6.1) se representa como una única
    secuencia (un cromosoma/plásmido) — se exige que el FASTA
    contenga exactamente un registro.
    """
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) == 0:
        raise GenomeLoadError(f"El archivo FASTA '{fasta_path}' no contiene ninguna secuencia.")
    if len(records) > 1:
        raise GenomeLoadError(
            f"El archivo FASTA '{fasta_path}' contiene {len(records)} secuencias; "
            f"por ahora solo se admite una (genoma de una única molécula)."
        )
    record = records[0]
    return record.id, str(record.seq).upper()


def _load_annotation(annotation_path: Path) -> List[GeneAnnotation]:
    """
    Lee la anotación de genes desde un CSV con columnas obligatorias
    gene_id, start, end, strand, y una columna opcional promoter_id.
    """
    genes: List[GeneAnnotation] = []
    seen_ids = set()

    with open(annotation_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required_columns = {"gene_id", "start", "end", "strand"}
        missing = required_columns - set(reader.fieldnames or [])
        if missing:
            raise GenomeLoadError(
                f"Faltan columnas obligatorias en '{annotation_path}': {sorted(missing)}"
            )

        for row_num, row in enumerate(reader, start=2):  # la fila 1 es la cabecera
            gene_id = (row.get("gene_id") or "").strip()
            if not gene_id:
                raise GenomeLoadError(f"Fila {row_num}: gene_id vacío.")
            if gene_id in seen_ids:
                raise GenomeLoadError(f"Fila {row_num}: gene_id '{gene_id}' duplicado.")
            seen_ids.add(gene_id)

            try:
                start = int(row["start"])
                end = int(row["end"])
            except (ValueError, KeyError):
                raise GenomeLoadError(
                    f"Fila {row_num} (gen '{gene_id}'): start/end deben ser enteros."
                ) from None

            strand = (row.get("strand") or "").strip()
            promoter_id = (row.get("promoter_id") or "").strip() or None

            genes.append(GeneAnnotation(gene_id, start, end, strand, promoter_id))

    return genes


def load_genome(fasta_path: str, annotation_path: str) -> GenomeData:
    """
    Punto de entrada público: carga un genoma completo (secuencia +
    anotación) desde dos archivos externos.

    Valida que las coordenadas de cada gen quepan dentro de la
    secuencia cargada — los errores de anotación se detectan aquí,
    no silenciosamente más adelante en biology/genome.py.
    """
    fasta_path_p = Path(fasta_path)
    annotation_path_p = Path(annotation_path)

    if not fasta_path_p.exists():
        raise GenomeLoadError(f"No se encuentra el archivo FASTA: {fasta_path_p}")
    if not annotation_path_p.exists():
        raise GenomeLoadError(f"No se encuentra el archivo de anotación: {annotation_path_p}")

    sequence_id, sequence = _load_sequence(fasta_path_p)
    genes = _load_annotation(annotation_path_p)

    for gene in genes:
        if gene.end > len(sequence):
            raise GenomeLoadError(
                f"Gen '{gene.gene_id}': end={gene.end} excede la longitud de la "
                f"secuencia ({len(sequence)})."
            )

    return GenomeData(sequence_id=sequence_id, sequence=sequence, genes=genes)
