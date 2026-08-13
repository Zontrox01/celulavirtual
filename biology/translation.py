"""
biology/translation.py

TranslationModule: segundo módulo biológico (Capa 2), construido sobre
genes ya transcritos por GenomeModule (biology/genome.py). Por cada
gen, registra una especie de proteína y una reacción de traducción de
un único paso (Ribosoma + ARNm -> Ribosoma + ARNm + Proteína) — misma
simplificación deliberada que en GenomeModule (whitepaper, sección 5,
nota de implementación de la Fase 1): sin estados intermedios de
elongación codón a codón, que quedan como refinamiento futuro sin
romper esta interfaz.

No conoce nada de transcripción, metabolismo ni depuración — solo
necesita saber qué especies de ARNm existen ya en la célula, a través
de los InstalledGene que devuelve GenomeModule.install().
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from cell import Cell
from engine.reactions import Reaction
from biology.genome import InstalledGene

DEFAULT_RIBOSOME_SPECIES_ID = "Ribosome"


class TranslationModuleError(ValueError):
    """Se lanza si la traducción no puede instalarse (ARNm inexistente, conflictos de nombres, etc.)."""


@dataclass(frozen=True)
class TranslatedGene:
    """Referencia a las especies/reacción generadas para la traducción de un gen concreto."""

    gene_id: str
    mrna_species_id: str
    protein_species_id: str
    translation_reaction_id: str


class TranslationModule:
    """
    Instala una reacción de traducción por cada gen ya transcrito.

    El ribosoma se trata como maquinaria preexistente de la célula
    (igual que la ARN-polimerasa en GenomeModule): su cantidad inicial
    se indica al instalar, no se deriva del genoma.
    """

    def __init__(self, installed_genes: Sequence[InstalledGene]) -> None:
        if not installed_genes:
            raise TranslationModuleError(
                "No se han pasado genes instalados: TranslationModule necesita el "
                "resultado de GenomeModule.install() para saber qué ARNm existen."
            )
        self.installed_genes = list(installed_genes)
        self._translated: Dict[str, TranslatedGene] = {}

    def install(
        self,
        cell: Cell,
        ribosome_species_id: str = DEFAULT_RIBOSOME_SPECIES_ID,
        ribosome_initial_count: int = 20,
        base_translation_rate: float = 1.0,
        rbs_strengths: Optional[Dict[str, float]] = None,
        atp_species_id: Optional[str] = None,
        atp_cost_per_translation: int = 0,
    ) -> List[TranslatedGene]:
        """
        Registra en `cell` las especies y reacciones de traducción de
        todos los genes ya transcritos.

        `rbs_strengths`: multiplicador de tasa por gen, indexado por
        gene_id. Los genes sin entrada usan multiplicador 1.0.

        `atp_species_id` / `atp_cost_per_translation`: acoplamiento
        energético opcional (whitepaper, Fase 3), mismo mecanismo que
        en GenomeModule.install(). Por defecto no hay coste de ATP.
        """
        rbs_strengths = rbs_strengths or {}
        couple_to_atp = atp_species_id is not None and atp_cost_per_translation > 0

        if couple_to_atp and not cell.species.has_species(atp_species_id):
            raise TranslationModuleError(
                f"Se pidió acoplar la traducción a ATP ('{atp_species_id}'), pero esa "
                f"especie no existe en la célula (¿se instaló primero MetabolismModule?)."
            )

        if not cell.species.has_species(ribosome_species_id):
            cell.add_species(ribosome_species_id, ribosome_initial_count)

        installed: List[TranslatedGene] = []
        for gene in self.installed_genes:
            if not cell.species.has_species(gene.mrna_species_id):
                raise TranslationModuleError(
                    f"El gen '{gene.gene_id}' referencia la especie de ARNm "
                    f"'{gene.mrna_species_id}', que no existe en la célula "
                    f"(¿se instaló primero GenomeModule sobre esta misma Cell?)."
                )

            protein_species_id = f"Prot_{gene.gene_id}"
            reaction_id = f"translate_{gene.gene_id}"

            if cell.species.has_species(protein_species_id):
                raise TranslationModuleError(
                    f"Conflicto al instalar la traducción de '{gene.gene_id}': "
                    f"'{protein_species_id}' ya existe en la célula "
                    f"(¿se ha instalado esta traducción dos veces?)."
                )

            cell.add_species(protein_species_id, 0)

            rate = base_translation_rate * rbs_strengths.get(gene.gene_id, 1.0)

            reactants = {ribosome_species_id: 1, gene.mrna_species_id: 1}
            products = {
                ribosome_species_id: 1,
                gene.mrna_species_id: 1,
                protein_species_id: 1,
            }
            if couple_to_atp:
                reactants[atp_species_id] = atp_cost_per_translation

            cell.add_reaction(
                Reaction(
                    reaction_id,
                    reactants=reactants,
                    products=products,
                    rate_constant=rate,
                )
            )

            record = TranslatedGene(
                gene_id=gene.gene_id,
                mrna_species_id=gene.mrna_species_id,
                protein_species_id=protein_species_id,
                translation_reaction_id=reaction_id,
            )
            installed.append(record)
            self._translated[gene.gene_id] = record

        return installed

    def get_translated(self, gene_id: str) -> TranslatedGene:
        try:
            return self._translated[gene_id]
        except KeyError:
            raise TranslationModuleError(
                f"El gen '{gene_id}' no se ha traducido todavía (¿llamaste a install()?)."
            ) from None
