"""
engine/debugger.py

Depurador interactivo del motor de simulación (Fase 0.5).

Responsabilidad única: envolver un GillespieSSA para dar control de
ejecución paso a paso, breakpoints condicionales sobre el estado
molecular, e inspección del estado y de las propensidades de todas
las reacciones candidatas en cada pausa (whitepaper, sección 4).

No altera la física del motor SSA subyacente: reutiliza directamente
GillespieSSA.step(), por lo que la distribución estadística de
resultados es idéntica a la del modo run() normal (criterio de
validación de la sección 11 del whitepaper).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional

from engine.ssa import GillespieSSA, SSAEvent, SimulationEnded, PropensityExhausted


BreakpointCondition = Callable[[Dict[str, int]], bool]


class DebuggerEndedError(Exception):
    """Se lanza al intentar avanzar un depurador cuya simulación ya ha terminado."""


@dataclass(frozen=True)
class Breakpoint:
    """Condición de parada, evaluada sobre una snapshot del estado de especies tras cada paso."""

    breakpoint_id: str
    condition: BreakpointCondition
    description: str = ""

    def is_triggered(self, state: Dict[str, int]) -> bool:
        return self.condition(state)


class Debugger:
    """
    Envoltorio de depuración sobre un GillespieSSA ya existente.

    Uso:
        debugger = Debugger(ssa)
        debugger.set_breakpoint("atp_bajo", lambda s: s.get("ATP", 0) < 50)
        debugger.step()                     # avanza una reacción
        debugger.step(n=10)                 # avanza hasta 10 reacciones (o hasta breakpoint)
        debugger.run_until_breakpoint()     # avanza hasta que se cumpla algún breakpoint
        debugger.current_event              # último evento disparado
        debugger.pending_propensities()     # propensidades de TODAS las reacciones candidatas

    `on_event`: callback opcional invocado tras cada evento ejecutado
    (justo después de añadirlo al histórico propio del depurador). Lo
    usa cell.py para mantener sincronizado su propio histórico de
    eventos cuando la ejecución pasa por el modo depurado en vez de
    por run(), de forma que get_trajectory() vea también esos eventos.
    """

    def __init__(
        self,
        ssa: GillespieSSA,
        on_event: Optional[Callable[[SSAEvent], None]] = None,
    ) -> None:
        self._ssa = ssa
        self._history: List[SSAEvent] = []
        self._breakpoints: Dict[str, Breakpoint] = {}
        self._last_triggered_breakpoint: Optional[str] = None
        self._ended = False
        self._on_event = on_event

    # ------------------------------------------------------------------
    # Breakpoints
    # ------------------------------------------------------------------

    def set_breakpoint(
        self, breakpoint_id: str, condition: BreakpointCondition, description: str = ""
    ) -> None:
        """Registra un breakpoint. Se evalúa después de cada paso, nunca antes."""
        self._breakpoints[breakpoint_id] = Breakpoint(breakpoint_id, condition, description)

    def remove_breakpoint(self, breakpoint_id: str) -> None:
        self._breakpoints.pop(breakpoint_id, None)

    def clear_breakpoints(self) -> None:
        self._breakpoints.clear()

    def _check_breakpoints(self) -> Optional[str]:
        """ID del primer breakpoint que se cumple con el estado actual, o None."""
        state = self._ssa.species.snapshot()
        for bp_id, bp in self._breakpoints.items():
            if bp.is_triggered(state):
                return bp_id
        return None

    # ------------------------------------------------------------------
    # Ejecución paso a paso
    # ------------------------------------------------------------------

    def step(self, n: int = 1) -> List[SSAEvent]:
        """
        Avanza hasta n reacciones, o hasta que se cumpla algún breakpoint,
        lo que ocurra antes. Devuelve los eventos ejecutados en esta
        llamada (puede ser más corta que n si se disparó un breakpoint
        o la simulación terminó).
        """
        if self._ended:
            raise DebuggerEndedError("La simulación ya ha terminado.")

        executed: List[SSAEvent] = []
        self._last_triggered_breakpoint = None

        for _ in range(n):
            try:
                event = self._ssa.step()
            except PropensityExhausted:
                self._ended = True
                break
            except SimulationEnded:
                # Pausa por max_time: no es un fin permanente. No se marca
                # self._ended, para que subir max_time y seguir llamando a
                # step() pueda reanudar la misma simulación.
                break

            executed.append(event)
            self._history.append(event)
            if self._on_event is not None:
                self._on_event(event)

            triggered = self._check_breakpoints()
            if triggered is not None:
                self._last_triggered_breakpoint = triggered
                break

        return executed

    def run_until_breakpoint(self, max_steps: Optional[int] = None) -> List[SSAEvent]:
        """
        Avanza reacción a reacción hasta que se dispare algún breakpoint,
        la simulación termine, o se alcance max_steps (protección frente
        a bucles infinitos si no hay breakpoints definidos).
        """
        executed: List[SSAEvent] = []
        steps_taken = 0
        while True:
            if max_steps is not None and steps_taken >= max_steps:
                break
            batch = self.step(n=1)
            executed.extend(batch)
            steps_taken += len(batch)
            if not batch:
                break  # la simulación terminó sin ejecutar nada más
            if self._last_triggered_breakpoint is not None:
                break
        return executed

    # ------------------------------------------------------------------
    # Inspección
    # ------------------------------------------------------------------

    @property
    def current_event(self) -> Optional[SSAEvent]:
        """Último evento ejecutado, o None si todavía no se ha dado ningún paso."""
        return self._history[-1] if self._history else None

    @property
    def last_triggered_breakpoint(self) -> Optional[str]:
        """ID del breakpoint que detuvo la última llamada a step()/run_until_breakpoint(), si alguno."""
        return self._last_triggered_breakpoint

    @property
    def time(self) -> float:
        return self._ssa.time

    @property
    def has_ended(self) -> bool:
        return self._ended

    def current_state(self) -> Dict[str, int]:
        """Snapshot del estado molecular actual."""
        return self._ssa.species.snapshot()

    def pending_propensities(self) -> Dict[str, float]:
        """
        Propensidad actual de TODAS las reacciones candidatas, no solo
        la que disparó el último evento. Es la base de la futura
        visualización del espacio de propensidades (sección 5.1) y del
        análisis causal retrospectivo (sección 4.5).
        """
        return self._ssa.reactions.propensities(self._ssa.species)

    def history(self) -> List[SSAEvent]:
        """Histórico completo de eventos ejecutados a través de este depurador."""
        return list(self._history)

    def __repr__(self) -> str:
        return (
            f"Debugger(tiempo={self.time:.4f}, eventos={len(self._history)}, "
            f"breakpoints={list(self._breakpoints.keys())}, terminado={self._ended})"
        )
