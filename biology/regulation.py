"""
biology/regulation.py

RegulationModule: cuarto módulo biológico (Capa 2). Añade represión
transcripcional simple sobre un gen ya instalado (GenomeModule),
mediada por una proteína represora ya traducida (TranslationModule).

Modela la represión como un interruptor molecular de dos estados: el
gen libre (gene_X) puede transcribirse (ya lo hace vía la reacción de
transcripción que instaló GenomeModule); el represor puede unirse a él
formando un complejo reprimido (gene_X_repressed_by_Y), que no tiene
ninguna reacción de transcripción asociada. Mientras el gen está en
ese estado, su copia libre pasa a 0 — bloqueando la reacción de
transcripción ya existente de forma natural, sin tener que modificarla:

    Represor + gen_libre  <-- bind -->  gen_reprimido   (unbind)

No conoce nada de metabolismo ni depuración — solo añade esta pareja
de reacciones de unión/disociación sobre especies que ya existen en la
célula (instaladas por GenomeModule y TranslationModule).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict

from cell import Cell
from engine.reactions import Reaction
from biology.genome import InstalledGene
from biology.translation import TranslatedGene

DEFAULT_BINDING_RATE = 1.0
DEFAULT_UNBINDING_RATE = 1.0


class RegulationModuleError(ValueError):
    """Se lanza si la represión no puede instalarse (especies inexistentes, conflictos, etc.)."""


@dataclass(frozen=True)
class RepressionLink:
    """Referencia a las especies/reacciones generadas para un enlace de represión concreto."""

    target_gene_id: str
    repressor_gene_id: str
    free_gene_species_id: str
    repressed_gene_species_id: str
    repressor_species_id: str
    binding_reaction_id: str
    unbinding_reaction_id: str


class RegulationModule:
    """
    Instala represiones transcripcionales: un represor (proteína ya
    traducida) puede unirse a un gen diana (ya transcrito), formando un
    complejo reprimido que bloquea su transcripción mientras dura la
    unión — reversible, estocástica, sin pasos intermedios.

    Un mismo módulo puede instalar varios enlaces de represión (un
    represor sobre varios genes, o varios represores sobre genes
    distintos), para construir redes de regulación más complejas.
    """

    def __init__(self) -> None:
        self._links: Dict[str, RepressionLink] = {}

    def add_repression(
        self,
        cell: Cell,
        target_gene: InstalledGene,
        repressor: TranslatedGene,
        binding_rate: float = DEFAULT_BINDING_RATE,
        unbinding_rate: float = DEFAULT_UNBINDING_RATE,
    ) -> RepressionLink:
        """
        Instala en `cell` la represión de `target_gene` por la proteína
        de `repressor`, con las tasas de unión/disociación indicadas.
        """
        if not cell.species.has_species(target_gene.gene_species_id):
            raise RegulationModuleError(
                f"El gen diana '{target_gene.gene_id}' referencia la especie "
                f"'{target_gene.gene_species_id}', que no existe en la célula."
            )
        if not cell.species.has_species(repressor.protein_species_id):
            raise RegulationModuleError(
                f"El represor '{repressor.gene_id}' referencia la especie "
                f"'{repressor.protein_species_id}', que no existe en la célula."
            )

        repressed_species_id = f"{target_gene.gene_species_id}_repressed_by_{repressor.gene_id}"
        binding_reaction_id = f"bind_{repressor.gene_id}_to_{target_gene.gene_id}"
        unbinding_reaction_id = f"unbind_{repressor.gene_id}_from_{target_gene.gene_id}"

        if binding_reaction_id in self._links:
            raise RegulationModuleError(
                f"Conflicto: la represión de '{target_gene.gene_id}' por "
                f"'{repressor.gene_id}' ya está instalada."
            )

        cell.add_species(repressed_species_id, 0)

        cell.add_reaction(
            Reaction(
                binding_reaction_id,
                reactants={target_gene.gene_species_id: 1, repressor.protein_species_id: 1},
                products={repressed_species_id: 1},
                rate_constant=binding_rate,
            )
        )
        cell.add_reaction(
            Reaction(
                unbinding_reaction_id,
                reactants={repressed_species_id: 1},
                products={target_gene.gene_species_id: 1, repressor.protein_species_id: 1},
                rate_constant=unbinding_rate,
            )
        )

        link = RepressionLink(
            target_gene_id=target_gene.gene_id,
            repressor_gene_id=repressor.gene_id,
            free_gene_species_id=target_gene.gene_species_id,
            repressed_gene_species_id=repressed_species_id,
            repressor_species_id=repressor.protein_species_id,
            binding_reaction_id=binding_reaction_id,
            unbinding_reaction_id=unbinding_reaction_id,
        )
        self._links[binding_reaction_id] = link
        return link

    def get_link(self, target_gene_id: str, repressor_gene_id: str) -> RepressionLink:
        key = f"bind_{repressor_gene_id}_to_{target_gene_id}"
        try:
            return self._links[key]
        except KeyError:
            raise RegulationModuleError(
                f"No hay represión instalada de '{target_gene_id}' por '{repressor_gene_id}'."
            ) from None
