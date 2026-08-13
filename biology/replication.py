"""
biology/replication.py

ReplicationModule: sexto módulo biológico (Capa 2) y cierre del ciclo
celular central del MVP. Dispara la división de una célula cuando se
cumple una condición configurable (p. ej. masa proteica acumulada),
repartiendo estocásticamente (binomial, p=0.5) cada especie molecular
entre dos células hijas — tal como ocurre en la segregación real de
moléculas entre células hijas al dividirse.

Simplificación deliberada: no se modela la replicación del ADN como
una reacción bioquímica explícita (síntesis de nuevas hebras, consumo
de dNTPs); se asume que, en el instante de la división, cada especie
presente se reparte igual que cualquier otra. Modelar la replicación
del ADN paso a paso queda como refinamiento futuro, análogo a la
simplificación ya aplicada en GenomeModule/TranslationModule.

Aviso importante: con reparto binomial puro y una única copia de gen
antes de dividirse, una de las dos hijas puede quedarse sin ninguna
copia de ese gen (nunca podría volver a transcribirlo). Por eso se
incluye `replicate_gene_copies()`: duplicar las copias de gen justo
antes de dividir, representando que la replicación del ADN ya ha
ocurrido, sin modelar el mecanismo paso a paso.

No conoce nada de metabolismo, regulación ni depuración — solo decide
"¿toca dividirse?" y, si toca, reparte el estado.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import numpy as np

from cell import Cell

DivisionCondition = Callable[[Dict[str, int]], bool]


class ReplicationModuleError(ValueError):
    """Se lanza si la división no puede ejecutarse (condición inválida, etc.)."""


def total_count_condition(species_ids: Iterable[str], threshold: int) -> DivisionCondition:
    """
    Fábrica de condiciones de división: dispara cuando la suma de las
    cantidades de `species_ids` alcanza `threshold`. Ejemplo típico:
    dividir cuando la masa proteica total supera cierto umbral.
    """
    species_ids = tuple(species_ids)
    if threshold <= 0:
        raise ReplicationModuleError(f"threshold debe ser positivo (recibido {threshold}).")

    def condition(state: Dict[str, int]) -> bool:
        return sum(state.get(species_id, 0) for species_id in species_ids) >= threshold

    return condition


@dataclass(frozen=True)
class DivisionEvent:
    """Registro de una división ya ocurrida, para inspección posterior."""

    parent_time: float
    parent_event_count: int
    daughter_a_state: Dict[str, int]
    daughter_b_state: Dict[str, int]


class ReplicationModule:
    """
    Comprueba una condición de división sobre el estado de una célula
    y, cuando se cumple, produce dos células hijas con reparto
    binomial de todas las especies presentes.
    """

    def __init__(self, condition: DivisionCondition, seed: Optional[int] = None) -> None:
        self.condition = condition
        self._rng = np.random.default_rng(seed)
        self.history: List[DivisionEvent] = []

    def should_divide(self, cell: Cell) -> bool:
        return self.condition(cell.species.snapshot())

    def check_and_divide(self, cell: Cell) -> Optional[Tuple[Cell, Cell]]:
        """Si la condición se cumple, divide y devuelve las dos hijas; si no, devuelve None."""
        if not self.should_divide(cell):
            return None
        return self.divide(cell)

    def divide(self, cell: Cell) -> Tuple[Cell, Cell]:
        """
        Crea dos células hijas repartiendo binomialmente (p=0.5) cada
        especie del estado actual de `cell`, sin comprobar la condición
        de división (usar check_and_divide() para respetarla). Las
        hijas heredan las mismas reacciones (mismos objetos `Reaction`,
        inmutables) que la célula madre: dividir no altera la red de
        reacciones, solo el estado molecular. Son `Cell` nuevas: su
        reloj arranca en 0, como corresponde al nacimiento de una
        célula.
        """
        state = cell.species.snapshot()

        daughter_a_state: Dict[str, int] = {}
        daughter_b_state: Dict[str, int] = {}
        for species_id, count in state.items():
            n_a = int(self._rng.binomial(count, 0.5)) if count > 0 else 0
            daughter_a_state[species_id] = n_a
            daughter_b_state[species_id] = count - n_a

        daughter_a = self._build_daughter(daughter_a_state, cell)
        daughter_b = self._build_daughter(daughter_b_state, cell)

        self.history.append(
            DivisionEvent(
                parent_time=cell.time,
                parent_event_count=len(cell.get_events()),
                daughter_a_state=daughter_a_state,
                daughter_b_state=daughter_b_state,
            )
        )

        return daughter_a, daughter_b

    @staticmethod
    def replicate_gene_copies(cell: Cell, gene_species_ids: Iterable[str]) -> None:
        """
        Duplica el número de copias de las especies de gen indicadas,
        representando que la replicación del ADN ya ha ocurrido justo
        antes de dividir (ver aviso en el docstring del módulo). Llamar
        a esto antes de divide()/check_and_divide() si se quiere que
        ambas hijas tengan una probabilidad razonable de heredar una
        copia funcional de cada gen.
        """
        for species_id in gene_species_ids:
            current = cell.species.get_count(species_id)
            cell.species.change_count(species_id, current)  # de N a 2N

    def _build_daughter(self, state: Dict[str, int], parent: Cell) -> Cell:
        daughter_seed = int(self._rng.integers(0, 2**31 - 1))
        daughter = Cell(seed=daughter_seed)
        for species_id, count in state.items():
            daughter.add_species(species_id, count)
        for reaction in parent.reactions.reactions():
            daughter.add_reaction(reaction)
        return daughter
