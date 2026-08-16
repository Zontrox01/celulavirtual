"""
gui/simulation_worker.py

Ejecuta la simulación en un QThread aparte. Sin esto, correr miles de
pasos (o max_time largos) congelaría la ventana entera, porque Qt
procesa eventos de interfaz en el mismo hilo donde se ejecuta el
código de la aplicación a menos que se delegue explícitamente a otro
hilo.

No contiene lógica de negocio propia: delega en gui/app_state.py
(run_steps/run_time), y solo se encarga de la mecánica de hilos y de
emitir señales que la interfaz pueda escuchar seguro entre hilos.
"""

from __future__ import annotations
from typing import Optional

from PySide6.QtCore import QThread, Signal

from gui.app_state import AppState


class SimulationWorker(QThread):
    """
    Ejecuta run_steps(n) o run_time(seconds) sobre un AppState en un
    hilo aparte. Emite:
      - finished_run(int): número de eventos ejecutados, al terminar.
      - error(str): mensaje de error, si algo falla durante la ejecución.
    """

    finished_run = Signal(int)
    error = Signal(str)

    def __init__(
        self,
        app_state: AppState,
        steps: Optional[int] = None,
        seconds: Optional[float] = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        if (steps is None) == (seconds is None):
            raise ValueError("Debe indicarse exactamente uno de: steps o seconds.")
        self._app_state = app_state
        self._steps = steps
        self._seconds = seconds

    def run(self) -> None:  # se ejecuta en el hilo secundario, no en el de la interfaz
        try:
            if self._steps is not None:
                events = self._app_state.run_steps(self._steps)
            else:
                events = self._app_state.run_time(self._seconds)
            self.finished_run.emit(len(events))
        except Exception as exc:  # noqa: BLE001 — cualquier fallo debe llegar a la interfaz, no perderse en el hilo
            self.error.emit(str(exc))
