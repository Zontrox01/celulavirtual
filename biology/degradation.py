"""
biology/degradation.py

DegradationModule: tercer módulo biológico (Capa 2), construido sobre
genes ya transcritos (GenomeModule) y, opcionalmente, ya traducidos
(TranslationModule). Registra una reacción de degradación de un solo
paso por cada ARNm y por cada proteína: especie -> nada.

Sin este módulo, ARNm y proteínas solo se acumulan (whitepaper,
sección 5): GenomeModule y TranslationModule generan producción pero
ninguna vía de desaparición. Es la pieza que faltaba para que el
sistema alcance un equilibrio dinámico en vez de crecer sin límite.

La tasa de degradación se deriva de la vida media (half-life), la
magnitud que se suele reportar en literatura biológica, mediante la
relación estándar de desintegración de primer orden:
    k = ln(2) / vida_media
en vez de pedir directamente una constante de velocidad poco intuitiva.

No conoce nada de transcripción, traducción ni metabolismo — solo
necesita saber qué especies de ARNm/proteína existen ya en la célula.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from cell import Cell
from engine.reactions import Reaction
from biology.genome import InstalledGene
from biology.translation import TranslatedGene

DEFAULT_MRNA_HALF_LIFE = 300.0      # segundos (~5 min, orden de magnitud típico en bacterias)
DEFAULT_PROTEIN_HALF_LIFE = 3600.0  # segundos (~1 h, las proteínas suelen durar más que el ARNm)


class DegradationModuleError(ValueError):
    """Se lanza si la degradación no puede instalarse (especie inexistente, conflictos, vida media inválida)."""


@dataclass(frozen=True)
class DegradationRecord:
    """Referencia a la reacción de degradación generada para una especie concreta."""

    gene_id: str
    species_id: str
    species_kind: str  # "mRNA" o "protein"
    reaction_id: str
    half_life: float


def _rate_from_half_life(half_life: float, label: str) -> float:
    if half_life <= 0:
        raise DegradationModuleError(f"{label}: la vida media debe ser positiva (recibida {half_life}).")
    return math.log(2) / half_life


class DegradationModule:
    """
    Instala reacciones de degradación de ARNm (siempre) y de proteínas
    (si se proporcionan genes traducidos) sobre una célula ya montada
    con GenomeModule y, opcionalmente, TranslationModule.
    """

    def __init__(
        self,
        installed_genes: Sequence[InstalledGene],
        translated_genes: Sequence[TranslatedGene] = (),
    ) -> None:
        if not installed_genes:
            raise DegradationModuleError(
                "No se han pasado genes instalados: DegradationModule necesita el "
                "resultado de GenomeModule.install() para saber qué ARNm existen."
            )
        self.installed_genes = list(installed_genes)
        self.translated_genes = list(translated_genes)
        self._records: Dict[str, DegradationRecord] = {}

    def install(
        self,
        cell: Cell,
        mrna_half_life: float = DEFAULT_MRNA_HALF_LIFE,
        protein_half_life: float = DEFAULT_PROTEIN_HALF_LIFE,
        mrna_half_lives: Optional[Dict[str, float]] = None,
        protein_half_lives: Optional[Dict[str, float]] = None,
    ) -> List[DegradationRecord]:
        """
        Registra en `cell` una reacción de degradación por cada ARNm
        (siempre) y por cada proteína (si se pasaron genes traducidos
        al construir el módulo).

        `mrna_half_lives` / `protein_half_lives`: vida media específica
        por gen, indexada por gene_id. Los genes sin entrada usan el
        valor por defecto (`mrna_half_life` / `protein_half_life`).
        """
        mrna_half_lives = mrna_half_lives or {}
        protein_half_lives = protein_half_lives or {}

        records: List[DegradationRecord] = []

        for gene in self.installed_genes:
            if not cell.species.has_species(gene.mrna_species_id):
                raise DegradationModuleError(
                    f"El gen '{gene.gene_id}' referencia la especie de ARNm "
                    f"'{gene.mrna_species_id}', que no existe en la célula."
                )

            reaction_id = f"degrade_{gene.mrna_species_id}"
            if reaction_id in self._records:
                raise DegradationModuleError(
                    f"Conflicto: la degradación de '{gene.mrna_species_id}' ya está instalada."
                )

            half_life = mrna_half_lives.get(gene.gene_id, mrna_half_life)
            rate = _rate_from_half_life(half_life, f"ARNm de '{gene.gene_id}'")

            cell.add_reaction(
                Reaction(reaction_id, reactants={gene.mrna_species_id: 1}, products={}, rate_constant=rate)
            )

            record = DegradationRecord(
                gene_id=gene.gene_id,
                species_id=gene.mrna_species_id,
                species_kind="mRNA",
                reaction_id=reaction_id,
                half_life=half_life,
            )
            records.append(record)
            self._records[reaction_id] = record

        for gene in self.translated_genes:
            if not cell.species.has_species(gene.protein_species_id):
                raise DegradationModuleError(
                    f"El gen '{gene.gene_id}' referencia la especie de proteína "
                    f"'{gene.protein_species_id}', que no existe en la célula."
                )

            reaction_id = f"degrade_{gene.protein_species_id}"
            if reaction_id in self._records:
                raise DegradationModuleError(
                    f"Conflicto: la degradación de '{gene.protein_species_id}' ya está instalada."
                )

            half_life = protein_half_lives.get(gene.gene_id, protein_half_life)
            rate = _rate_from_half_life(half_life, f"proteína de '{gene.gene_id}'")

            cell.add_reaction(
                Reaction(reaction_id, reactants={gene.protein_species_id: 1}, products={}, rate_constant=rate)
            )

            record = DegradationRecord(
                gene_id=gene.gene_id,
                species_id=gene.protein_species_id,
                species_kind="protein",
                reaction_id=reaction_id,
                half_life=half_life,
            )
            records.append(record)
            self._records[reaction_id] = record

        return records

    def get_degradation(self, species_id: str) -> DegradationRecord:
        reaction_id = f"degrade_{species_id}"
        try:
            return self._records[reaction_id]
        except KeyError:
            raise DegradationModuleError(
                f"No hay degradación instalada para la especie '{species_id}'."
            ) from None
