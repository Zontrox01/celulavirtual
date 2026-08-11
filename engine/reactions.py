"""
engine/reactions.py

Gestor de reacciones del motor de simulación.

Responsabilidad única: definir reacciones (reactivos, productos, tasa)
y calcular su propensidad dado el estado actual de las especies
(engine/species.py). No decide cuál dispara ni cuándo — eso es
responsabilidad de engine/ssa.py (whitepaper_celula_virtual.md, sección 5.4).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List
import math

from engine.species import SpeciesManager, UnknownSpeciesError


class InvalidReactionError(ValueError):
    """Se lanza si una reacción está mal definida o referencia una especie inexistente."""


@dataclass(frozen=True)
class Reaction:
    """
    Reacción química elemental de cinética de acción de masas.

    reactants / products: {species_id: coeficiente_estequiométrico}
    rate_constant: constante cinética (k) de la reacción.

    La propensidad se calcula como:
        a = k * producto_sobre_reactivos( C(count_i, n_i) )
    donde C(count_i, n_i) es "count_i elige n_i" (combinaciones) — la
    forma estándar de propensidad de Gillespie para reacciones
    elementales (uni- o bimoleculares; no se contemplan órdenes
    superiores por ahora).
    """

    reaction_id: str
    reactants: Dict[str, int]
    products: Dict[str, int]
    rate_constant: float

    def __post_init__(self) -> None:
        if self.rate_constant < 0:
            raise InvalidReactionError(
                f"'{self.reaction_id}': la constante de velocidad no puede ser negativa "
                f"({self.rate_constant})."
            )
        for coeff_dict, label in ((self.reactants, "reactivos"), (self.products, "productos")):
            for species_id, coeff in coeff_dict.items():
                if coeff <= 0:
                    raise InvalidReactionError(
                        f"'{self.reaction_id}': coeficiente estequiométrico de "
                        f"'{species_id}' en {label} debe ser positivo (recibido {coeff})."
                    )
        if not self.reactants and not self.products:
            raise InvalidReactionError(
                f"'{self.reaction_id}': una reacción necesita al menos un reactivo o un producto."
            )

    def propensity(self, species: SpeciesManager) -> float:
        """
        Calcula la propensidad actual de esta reacción dado el estado de `species`.

        Devuelve 0.0 si no hay suficientes moléculas de algún reactivo
        para que la reacción pueda ocurrir, en vez de lanzar una excepción:
        propensidad nula es un resultado válido y esperado, no un error.
        """
        combinatorial_factor = 1.0
        for species_id, coeff in self.reactants.items():
            try:
                count = species.get_count(species_id)
            except UnknownSpeciesError:
                raise InvalidReactionError(
                    f"'{self.reaction_id}': el reactivo '{species_id}' no está registrado "
                    f"en el SpeciesManager."
                ) from None
            if count < coeff:
                return 0.0
            combinatorial_factor *= math.comb(count, coeff)
        return self.rate_constant * combinatorial_factor

    def apply(self, species: SpeciesManager) -> None:
        """
        Aplica el efecto estequiométrico de esta reacción: consume los
        reactivos y produce los productos. No comprueba la propensidad —
        quien llama (engine/ssa.py) es responsable de haber verificado
        que la reacción puede dispararse antes de invocar esto.
        """
        for species_id, coeff in self.reactants.items():
            species.change_count(species_id, -coeff)
        for species_id, coeff in self.products.items():
            species.change_count(species_id, coeff)

    def involved_species(self) -> set:
        """Conjunto de todas las especies implicadas (reactivos y productos)."""
        return set(self.reactants) | set(self.products)


class ReactionManager:
    """Catálogo de reacciones activas en el sistema."""

    def __init__(self) -> None:
        self._reactions: Dict[str, Reaction] = {}

    def add_reaction(self, reaction: Reaction) -> None:
        if reaction.reaction_id in self._reactions:
            raise ValueError(f"La reacción '{reaction.reaction_id}' ya está registrada.")
        self._reactions[reaction.reaction_id] = reaction

    def get_reaction(self, reaction_id: str) -> Reaction:
        try:
            return self._reactions[reaction_id]
        except KeyError:
            raise KeyError(f"Reacción '{reaction_id}' no registrada.") from None

    def reactions(self) -> List[Reaction]:
        """Todas las reacciones, en orden de inserción (relevante para reproducibilidad del SSA)."""
        return list(self._reactions.values())

    def propensities(self, species: SpeciesManager) -> Dict[str, float]:
        """Propensidad actual de cada reacción, indexada por reaction_id."""
        return {r.reaction_id: r.propensity(species) for r in self._reactions.values()}

    def total_propensity(self, species: SpeciesManager) -> float:
        return sum(self.propensities(species).values())

    def __len__(self) -> int:
        return len(self._reactions)

    def __repr__(self) -> str:
        return f"ReactionManager({list(self._reactions.keys())!r})"
