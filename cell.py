"""
cell.py

Orquestador `Cell` (Capa 3): compone SpeciesManager, ReactionManager y
GillespieSSA bajo una API pública única, tal como se definió en el
whitepaper (whitepaper_celula_virtual.md, sección 5.3).

En esta fase (Fase 0) expone solo ejecución normal (run). La ejecución
depurada (cell.debug()) se activará en la Fase 0.5, una vez exista
engine/debugger.py — el método ya está aquí como punto de entrada
reservado, para no tener que cambiar la interfaz pública más adelante.
"""

from __future__ import annotations
from typing import List, Optional

from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager
from engine.ssa import GillespieSSA, SSAEvent, SimulationEnded
from engine.trajectory import build_trajectory
from engine.debugger import Debugger


class Cell:
    """
    Célula virtual: contenedor de especies y reacciones, con capacidad
    de simular su evolución temporal mediante el motor SSA propio.

    Por ahora es un contenedor "vacío" en el sentido biológico — no
    carga genoma ni tiene módulos biológicos todavía (eso llega en la
    Fase 1 con biology/genome.py y el resto de módulos). Esta clase es
    el punto de composición al que esos módulos se irán enganchando.
    """

    def __init__(self, seed: Optional[int] = None) -> None:
        self.species = SpeciesManager()
        self.reactions = ReactionManager()
        self._seed = seed
        self._ssa: Optional[GillespieSSA] = None
        self._events: List[SSAEvent] = []
        self._initial_state: Optional[dict] = None

    # ------------------------------------------------------------------
    # Construcción del sistema
    # ------------------------------------------------------------------

    def add_species(self, species_id: str, initial_count: int = 0) -> "Cell":
        """Registra una especie molecular. Devuelve self para encadenar llamadas."""
        self.species.add_species(species_id, initial_count)
        return self

    def add_reaction(self, reaction: Reaction) -> "Cell":
        """Registra una reacción. Devuelve self para encadenar llamadas."""
        self.reactions.add_reaction(reaction)
        return self

    # ------------------------------------------------------------------
    # Ejecución
    # ------------------------------------------------------------------

    def _ensure_ssa(self) -> GillespieSSA:
        """
        Crea el motor SSA la primera vez que hace falta, sobre las
        especies/reacciones ya registradas. Creación perezosa para
        poder seguir añadiendo especies/reacciones hasta justo antes
        de ejecutar la primera vez.
        """
        if self._ssa is None:
            self._initial_state = self.species.snapshot()
            self._ssa = GillespieSSA(self.species, self.reactions, seed=self._seed)
        return self._ssa

    def run(
        self,
        max_steps: Optional[int] = None,
        max_time: Optional[float] = None,
    ) -> List[SSAEvent]:
        """
        Ejecuta la simulación en modo normal (sin depuración) hasta
        max_steps pasos, o hasta max_time de tiempo simulado, o hasta
        que la simulación se agote por sí sola — lo que ocurra antes.

        Los eventos generados se acumulan en el histórico de la célula,
        accesible después vía get_events(). Se puede llamar a run()
        varias veces seguidas para seguir avanzando la misma simulación.
        """
        ssa = self._ensure_ssa()
        if max_time is not None:
            ssa.max_time = max_time

        new_events: List[SSAEvent] = []
        try:
            for event in ssa:
                new_events.append(event)
                if max_steps is not None and len(new_events) >= max_steps:
                    break
        except SimulationEnded:
            pass

        self._events.extend(new_events)
        return new_events

    def debug(self) -> Debugger:
        """
        Punto de entrada a la ejecución depurada: paso a paso,
        breakpoints, inspección de propensidades (whitepaper, sección 4).

        El Debugger devuelto envuelve el mismo motor SSA que usa run(),
        y se le pasa un callback para que cada evento ejecutado ahí se
        añada también al histórico de la propia Cell — así run() y
        debug() pueden alternarse sin perder eventos, y get_trajectory()
        ve la simulación completa independientemente de por qué modo
        se haya avanzado cada paso.
        """
        ssa = self._ensure_ssa()
        return Debugger(ssa, on_event=self._events.append)

    # ------------------------------------------------------------------
    # Consulta de resultados
    # ------------------------------------------------------------------

    @property
    def time(self) -> float:
        """Tiempo simulado transcurrido hasta ahora."""
        return self._ssa.time if self._ssa is not None else 0.0

    def get_events(self) -> List[SSAEvent]:
        """Histórico completo de eventos ocurridos hasta ahora."""
        return list(self._events)

    def get_trajectory(self):
        """
        Trayectoria tiempo x especie x cantidad, apta para análisis y
        graficado (pandas.DataFrame), reconstruida a partir del estado
        inicial capturado antes del primer run() y del histórico de
        eventos ocurridos desde entonces (engine/trajectory.py).
        """
        if self._initial_state is None:
            raise RuntimeError(
                "Todavía no se ha ejecutado ninguna simulación (llama a run() primero)."
            )
        return build_trajectory(self._initial_state, self._events, self.reactions)

    def __repr__(self) -> str:
        return (
            f"Cell(especies={len(self.species)}, reacciones={len(self.reactions)}, "
            f"tiempo={self.time:.4f}, eventos={len(self._events)})"
        )
