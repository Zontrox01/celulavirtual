"""
engine/ssa.py

Bucle de simulación estocástica exacta (algoritmo de Gillespie, método directo).

Responsabilidad única: dado un SpeciesManager y un ReactionManager,
generar la secuencia de eventos (qué reacción dispara, en qué instante)
que sigue la cinética estocástica exacta. Está escrito como iterador
para poder pausarse paso a paso — es el punto de enganche que usará
engine/debugger.py en la Fase 0.5 (whitepaper, secciones 4.2 y 5.4).

No conoce nada de depuración, breakpoints ni biología — solo el
algoritmo de Gillespie puro sobre las abstracciones de species.py
y reactions.py.
"""

from __future__ import annotations
import math
import random
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional

from engine.species import SpeciesManager
from engine.reactions import ReactionManager


class SimulationEnded(Exception):
    """
    Se lanza cuando step() no puede completar un paso más en esta
    llamada. Puede ser una pausa reanudable (se alcanzó max_time: basta
    con subir max_time y volver a llamar a step()/run()) o el fin
    permanente de la simulación (ver PropensityExhausted, que hereda
    de esta clase para que capturar SimulationEnded siga cubriendo
    ambos casos si no interesa la distinción).
    """


class PropensityExhausted(SimulationEnded):
    """
    Se lanza cuando la propensidad total es 0: ninguna reacción puede
    ocurrir ya, bajo ninguna circunstancia. A diferencia de alcanzar
    max_time, esto sí es un fin permanente e irreversible.
    """


@dataclass(frozen=True)
class SSAEvent:
    """
    Registro inmutable de un paso de la simulación.

    Es la unidad mínima que se acumulará después en el log de eventos
    usado por engine/trajectory.py, engine/debugger.py y
    engine/forensics.py (whitepaper, secciones 4.3 y 4.5). Guarda las
    propensidades de TODAS las reacciones candidatas, no solo la que
    disparó, porque esa información es la base de la visualización del
    espacio de propensidades (sección 5.1) y del modo forense.
    """

    step: int
    time: float
    reaction_id: str
    propensities: Dict[str, float]
    total_propensity: float


class GillespieSSA:
    """
    Motor de simulación estocástica exacta (Gillespie, 1977 — método directo).

    Uso básico (recorre toda la simulación hasta que se agota):
        ssa = GillespieSSA(species, reactions, seed=42)
        for event in ssa:
            ...

    Uso paso a paso (lo que reutilizará el depurador):
        event = ssa.step()
    """

    def __init__(
        self,
        species: SpeciesManager,
        reactions: ReactionManager,
        seed: Optional[int] = None,
        max_time: Optional[float] = None,
    ) -> None:
        self.species = species
        self.reactions = reactions
        self.max_time = max_time
        self.time: float = 0.0
        self.step_count: int = 0
        self._rng = random.Random(seed)
        self._ended = False

    def _draw_time_and_reaction(
        self, propensities: Dict[str, float], total: float
    ) -> tuple[float, str]:
        """
        Sortea el tiempo de espera (exponencial de parámetro `total`) y
        la reacción que dispara (proporcional a su propensidad),
        siguiendo el método directo de Gillespie.
        """
        r1 = self._rng.random()
        while r1 <= 0.0:  # evita log(1/0); random.random() puede devolver 0.0
            r1 = self._rng.random()
        tau = (1.0 / total) * math.log(1.0 / r1)

        r2 = self._rng.random() * total
        cumulative = 0.0
        chosen_id: Optional[str] = None
        for reaction_id, a in propensities.items():
            cumulative += a
            if r2 <= cumulative:
                chosen_id = reaction_id
                break
        if chosen_id is None:
            # protección ante errores de redondeo en el borde superior
            chosen_id = next(reversed(list(propensities)))

        return tau, chosen_id

    def step(self) -> SSAEvent:
        """
        Ejecuta exactamente un paso del algoritmo: calcula propensidades,
        sortea tiempo y reacción, la aplica, avanza el reloj, y devuelve
        el evento resultante.

        Lanza SimulationEnded si ninguna reacción tiene propensidad
        positiva, o si el siguiente evento superaría max_time.
        """
        if self._ended:
            raise PropensityExhausted("La simulación ya ha terminado (agotada permanentemente).")

        propensities = self.reactions.propensities(self.species)
        total = sum(propensities.values())

        if total <= 0.0:
            self._ended = True
            raise PropensityExhausted(
                "Propensidad total = 0: ninguna reacción puede ocurrir ya."
            )

        tau, chosen_id = self._draw_time_and_reaction(propensities, total)
        new_time = self.time + tau

        if self.max_time is not None and new_time > self.max_time:
            # Pausa reanudable, NO fin permanente: no se marca self._ended.
            # Subir max_time y volver a llamar a step()/run() debe poder
            # continuar la misma simulación desde aquí.
            raise SimulationEnded(
                f"Se alcanzó max_time ({self.max_time}) antes del siguiente evento."
            )

        reaction = self.reactions.get_reaction(chosen_id)
        reaction.apply(self.species)

        self.time = new_time
        self.step_count += 1

        return SSAEvent(
            step=self.step_count,
            time=self.time,
            reaction_id=chosen_id,
            propensities=propensities,
            total_propensity=total,
        )

    def __iter__(self) -> Iterator[SSAEvent]:
        return self

    def __next__(self) -> SSAEvent:
        try:
            return self.step()
        except SimulationEnded:
            raise StopIteration from None

    def run(self, max_steps: Optional[int] = None) -> List[SSAEvent]:
        """
        Atajo de conveniencia: ejecuta hasta max_steps pasos (o hasta que
        la simulación termine sola), devolviendo la lista de eventos.
        Para control fino paso a paso, usar step() directamente
        (es lo que hará engine/debugger.py en la Fase 0.5).
        """
        events: List[SSAEvent] = []
        for i, event in enumerate(self):
            events.append(event)
            if max_steps is not None and (i + 1) >= max_steps:
                break
        return events
