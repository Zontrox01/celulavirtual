"""
engine/snapshot.py

Snapshots de estado + retroceso (undo) del motor SSA (whitepaper,
secciones 4.4.b y 5.3).

En vez de intentar invertir matemáticamente cada reacción (frágil:
habría que deshacer también el sorteo aleatorio que la eligió), este
módulo guarda snapshots periódicos del estado COMPLETO (especies +
estado del generador aleatorio + tiempo + contador de pasos) y
reconstruye cualquier punto intermedio reproduciendo hacia adelante
desde el snapshot más cercano. Como el generador aleatorio se restaura
bit a bit, reproducir hacia adelante genera EXACTAMENTE la misma
secuencia de eventos que ya había ocurrido — no es una simulación
nueva, es una repetición determinista de la que ya pasó.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List

from engine.ssa import GillespieSSA


class SnapshotError(ValueError):
    """Se lanza si no se puede tomar o restaurar un snapshot, o el retroceso pedido es inválido."""


@dataclass(frozen=True)
class Snapshot:
    """Estado completo y reproducible del motor SSA en un instante dado."""

    step_count: int
    state: dict  # lo que devuelve GillespieSSA.get_full_state()


class SnapshotStore:
    """
    Guarda snapshots periódicos de un GillespieSSA y permite retroceder
    a cualquier paso ya ocurrido, reproduciendo desde el snapshot más
    cercano anterior o igual a ese paso.
    """

    def __init__(self, ssa: GillespieSSA, interval: int = 50) -> None:
        if interval <= 0:
            raise SnapshotError(f"interval debe ser positivo (recibido {interval}).")
        self._ssa = ssa
        self.interval = interval
        self._snapshots: List[Snapshot] = []
        self._take_snapshot()  # snapshot inicial, siempre disponible como suelo

    def _take_snapshot(self) -> None:
        self._snapshots.append(
            Snapshot(step_count=self._ssa.step_count, state=self._ssa.get_full_state())
        )

    def maybe_snapshot(self) -> bool:
        """
        Llamar después de cada step() del motor: toma un snapshot nuevo
        si han pasado `interval` pasos desde el último. Devuelve True
        si se tomó un snapshot nuevo.
        """
        last = self._snapshots[-1]
        if self._ssa.step_count - last.step_count >= self.interval:
            self._take_snapshot()
            return True
        return False

    def nearest_snapshot_at_or_before(self, step_count: int) -> Snapshot:
        candidates = [s for s in self._snapshots if s.step_count <= step_count]
        if not candidates:
            raise SnapshotError(
                f"No hay ningún snapshot disponible en o antes del paso {step_count} "
                f"(el más antiguo es el paso {self._snapshots[0].step_count})."
            )
        return max(candidates, key=lambda s: s.step_count)

    def rewind_to(self, step_count: int) -> int:
        """
        Retrocede el GillespieSSA hasta exactamente `step_count`:
        restaura el snapshot más cercano anterior o igual, y reproduce
        hacia adelante paso a paso (de forma determinista, mismo
        generador aleatorio) hasta alcanzar el paso pedido.

        Devuelve cuántos pasos se reprodujeron desde el snapshot
        restaurado.
        """
        if step_count < 0:
            raise SnapshotError(f"step_count no puede ser negativo (recibido {step_count}).")
        if step_count > self._ssa.step_count:
            raise SnapshotError(
                f"No se puede retroceder a un paso futuro ({step_count}) "
                f"cuando la simulación va por el paso {self._ssa.step_count}."
            )

        snapshot = self.nearest_snapshot_at_or_before(step_count)
        self._ssa.restore_full_state(snapshot.state)

        # descartar snapshots posteriores al punto de restauración: si a
        # partir de aquí se avanza de nuevo, la trayectoria futura podría
        # divergir de la que ya se había explorado
        self._snapshots = [s for s in self._snapshots if s.step_count <= snapshot.step_count]

        replayed = 0
        while self._ssa.step_count < step_count:
            self._ssa.step()
            replayed += 1

        return replayed

    def snapshot_steps(self) -> List[int]:
        """Los pasos en los que hay snapshots disponibles (para inspección)."""
        return [s.step_count for s in self._snapshots]
