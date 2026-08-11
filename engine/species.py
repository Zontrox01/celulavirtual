"""
engine/species.py

Gestor de especies moleculares del motor de simulación.

Responsabilidad única: mantener el catálogo de especies (moléculas)
presentes en el sistema y sus cantidades actuales (números de copias).
No conoce reacciones ni cinética — esa lógica vive en engine/reactions.py
y engine/ssa.py respectivamente (whitepaper_celula_virtual.md, sección 5.4).
"""

from __future__ import annotations
from typing import Dict, Iterator


class UnknownSpeciesError(KeyError):
    """Se lanza al operar sobre una especie no registrada."""


class NegativeCountError(ValueError):
    """Se lanza si una operación dejaría la cantidad de una especie en negativo."""


class SpeciesManager:
    """
    Catálogo de especies moleculares y sus cantidades actuales.

    Cada especie se identifica por un ID único (str), siguiendo el
    esquema de identificación de la sección 7.2 del whitepaper
    (p.ej. "ATP", "mRNA_geneA", "Prot_geneA").
    """

    def __init__(self) -> None:
        self._counts: Dict[str, int] = {}

    def add_species(self, species_id: str, initial_count: int = 0) -> None:
        """Registra una nueva especie en el sistema."""
        if species_id in self._counts:
            raise ValueError(f"La especie '{species_id}' ya está registrada.")
        if initial_count < 0:
            raise NegativeCountError(
                f"No se puede crear '{species_id}' con cantidad negativa ({initial_count})."
            )
        self._counts[species_id] = initial_count

    def get_count(self, species_id: str) -> int:
        """Devuelve la cantidad actual de una especie."""
        try:
            return self._counts[species_id]
        except KeyError:
            raise UnknownSpeciesError(species_id) from None

    def set_count(self, species_id: str, count: int) -> None:
        """Fija directamente la cantidad de una especie (uso en tests / snapshots)."""
        if species_id not in self._counts:
            raise UnknownSpeciesError(species_id)
        if count < 0:
            raise NegativeCountError(f"'{species_id}' no puede tener cantidad negativa ({count}).")
        self._counts[species_id] = count

    def change_count(self, species_id: str, delta: int) -> None:
        """
        Aplica un cambio (positivo o negativo) a la cantidad de una especie.

        Es la operación que usará el bucle SSA (engine/ssa.py, Fase 0)
        al aplicar el efecto de una reacción disparada.
        """
        if species_id not in self._counts:
            raise UnknownSpeciesError(species_id)
        new_count = self._counts[species_id] + delta
        if new_count < 0:
            raise NegativeCountError(
                f"La reacción dejaría '{species_id}' en {new_count} (negativo)."
            )
        self._counts[species_id] = new_count

    def has_species(self, species_id: str) -> bool:
        return species_id in self._counts

    def species_ids(self) -> Iterator[str]:
        return iter(self._counts.keys())

    def snapshot(self) -> Dict[str, int]:
        """
        Copia inmutable del estado actual.

        La usará engine/debugger.py (Fase 0.5) para inspección,
        y engine/snapshot.py (Fase 5.3) como base del retroceso (undo).
        """
        return dict(self._counts)

    def restore(self, state: Dict[str, int]) -> None:
        """Restaura el estado completo desde una snapshot previa."""
        self._counts = dict(state)

    def __len__(self) -> int:
        return len(self._counts)

    def __repr__(self) -> str:
        return f"SpeciesManager({self._counts!r})"
