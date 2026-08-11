"""
engine/trajectory.py

Registro de trayectorias temporales.

Responsabilidad única: convertir el histórico de eventos (SSAEvent)
generado por engine/ssa.py en una trayectoria tiempo x especie x
cantidad, apta para análisis y graficado (whitepaper, sección 5).

No decide cuándo dispara una reacción — solo reconstruye, a partir de
un estado inicial y de la secuencia de reacciones que dispararon, cómo
evolucionó cada especie a lo largo del tiempo. Es determinista y no
depende del estado "en vivo" de ningún SpeciesManager: recibe el
estado inicial explícitamente y reaplica la estequiometría de cada
reacción registrada en `reactions`.
"""

from __future__ import annotations
from typing import Dict, List, Sequence

import pandas as pd

from engine.reactions import ReactionManager
from engine.ssa import SSAEvent


class TrajectoryError(ValueError):
    """Se lanza si los eventos no son reconstruibles con el estado/reacciones dados."""


def build_trajectory(
    initial_state: Dict[str, int],
    events: Sequence[SSAEvent],
    reactions: ReactionManager,
) -> pd.DataFrame:
    """
    Reconstruye la trayectoria completa a partir de un estado inicial
    y la secuencia de eventos ya ocurridos.

    Devuelve un DataFrame con una fila por evento (más la fila inicial,
    en tiempo 0), columnas 'step', 'time', 'reaction_id' (vacía en la
    fila inicial) y una columna por cada especie de `initial_state`.
    """
    species_ids = sorted(initial_state.keys())
    current = dict(initial_state)

    rows: List[dict] = [
        {"step": 0, "time": 0.0, "reaction_id": None, **current}
    ]

    for event in events:
        try:
            reaction = reactions.get_reaction(event.reaction_id)
        except KeyError:
            raise TrajectoryError(
                f"El evento del paso {event.step} referencia la reacción "
                f"'{event.reaction_id}', que no está registrada en `reactions`."
            ) from None

        for species_id, coeff in reaction.reactants.items():
            if species_id not in current:
                raise TrajectoryError(
                    f"La reacción '{event.reaction_id}' involucra la especie "
                    f"'{species_id}', ausente de initial_state."
                )
            current[species_id] -= coeff
        for species_id, coeff in reaction.products.items():
            if species_id not in current:
                raise TrajectoryError(
                    f"La reacción '{event.reaction_id}' involucra la especie "
                    f"'{species_id}', ausente de initial_state."
                )
            current[species_id] += coeff

        rows.append(
            {"step": event.step, "time": event.time, "reaction_id": event.reaction_id, **current}
        )

    df = pd.DataFrame(rows)
    ordered_columns = ["step", "time", "reaction_id"] + species_ids
    return df[ordered_columns]
