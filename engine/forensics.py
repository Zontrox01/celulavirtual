"""
engine/forensics.py

Análisis causal retrospectivo ("modo forense", whitepaper sección 4.5).

Responsabilidad única: dado un log de eventos ya ocurridos (lista de
SSAEvent, como el que acumula Cell/Debugger) y la red de reacciones
que los generó, responder preguntas retrospectivas sobre qué pasó y
por qué — sin ejecutar nada nuevo, solo consultando y correlacionando
el historial ya registrado.

Contraste con el depurador paso a paso (engine/debugger.py):

    Depurador paso a paso:          ¿Qué va a pasar?     (hacia adelante)
    Análisis causal retrospectivo:  ¿Por qué pasó?        (hacia atrás)
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from engine.reactions import ReactionManager
from engine.ssa import SSAEvent


@dataclass(frozen=True)
class SpeciesContribution:
    """Cuánto contribuyó una reacción concreta al cambio neto de una especie en una ventana."""

    reaction_id: str
    times_fired: int
    net_change: int  # positivo si esa reacción produjo más de lo que consumió, negativo si al revés


@dataclass(frozen=True)
class ForensicReport:
    """Resultado de un análisis retrospectivo sobre una especie y una ventana temporal."""

    species_id: str
    window_start: float
    window_end: float
    net_change: int
    contributions: List[SpeciesContribution]  # ordenadas por |net_change| descendente

    def top_consumer(self) -> Optional[SpeciesContribution]:
        """La reacción que más contribuyó a REDUCIR la especie en esta ventana, si alguna."""
        consumers = [c for c in self.contributions if c.net_change < 0]
        return min(consumers, key=lambda c: c.net_change) if consumers else None

    def top_producer(self) -> Optional[SpeciesContribution]:
        """La reacción que más contribuyó a AUMENTAR la especie en esta ventana, si alguna."""
        producers = [c for c in self.contributions if c.net_change > 0]
        return max(producers, key=lambda c: c.net_change) if producers else None


class ForensicsAnalyzer:
    """
    Analiza retrospectivamente un log de eventos ya ocurridos, sin
    ejecutar ninguna simulación nueva — solo lee `events` y consulta
    la estequiometría de cada reacción en `reactions`.
    """

    def __init__(self, events: Sequence[SSAEvent], reactions: ReactionManager) -> None:
        self.events = list(events)
        self.reactions = reactions

    def events_in_window(self, start: float, end: float) -> List[SSAEvent]:
        """Eventos cuyo instante cae dentro de [start, end]."""
        return [e for e in self.events if start <= e.time <= end]

    def last_event_before(self, t: float) -> Optional[SSAEvent]:
        """El último evento ocurrido en o antes del instante t, si alguno."""
        candidates = [e for e in self.events if e.time <= t]
        return max(candidates, key=lambda e: e.time) if candidates else None

    def reactions_fired_count(self, start: float, end: float) -> Dict[str, int]:
        """Cuántas veces disparó cada reacción dentro de la ventana [start, end]."""
        counts: Dict[str, int] = {}
        for event in self.events_in_window(start, end):
            counts[event.reaction_id] = counts.get(event.reaction_id, 0) + 1
        return counts

    def analyze_species(self, species_id: str, window_start: float, window_end: float) -> ForensicReport:
        """
        Reconstruye qué reacciones contribuyeron al cambio de
        `species_id` durante [window_start, window_end], y en qué
        medida cada una — la operación central del modo forense.

        Las reacciones catalíticas (que consumen y producen la misma
        cantidad de la especie, p. ej. un represor que se libera y se
        vuelve a unir) tienen cambio neto 0 y no aparecen en el
        informe: no "explican" ningún cambio, aunque hayan disparado.
        Eventos que referencian una reacción ausente de `reactions`
        se ignoran (log de otra red de reacciones), en vez de fallar.
        """
        window_events = self.events_in_window(window_start, window_end)

        times_fired: Dict[str, int] = {}
        net_change_per_reaction: Dict[str, int] = {}

        for event in window_events:
            try:
                reaction = self.reactions.get_reaction(event.reaction_id)
            except KeyError:
                continue

            delta = reaction.products.get(species_id, 0) - reaction.reactants.get(species_id, 0)
            if delta == 0:
                continue

            times_fired[event.reaction_id] = times_fired.get(event.reaction_id, 0) + 1
            net_change_per_reaction[event.reaction_id] = (
                net_change_per_reaction.get(event.reaction_id, 0) + delta
            )

        contributions = [
            SpeciesContribution(
                reaction_id=reaction_id,
                times_fired=times_fired[reaction_id],
                net_change=net_change_per_reaction[reaction_id],
            )
            for reaction_id in net_change_per_reaction
        ]
        contributions.sort(key=lambda c: abs(c.net_change), reverse=True)

        return ForensicReport(
            species_id=species_id,
            window_start=window_start,
            window_end=window_end,
            net_change=sum(c.net_change for c in contributions),
            contributions=contributions,
        )
