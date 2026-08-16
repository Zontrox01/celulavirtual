"""
gui/execution_panel.py

Panel de ejecución (pestaña "Ejecución"): correr la simulación en
segundo plano (gui/simulation_worker.py), ver el estado molecular
actual en una tabla, la trayectoria en una gráfica (matplotlib
embebido), y disparar una división celular (biology/replication.py)
cuando se cumple una condición.
"""

from __future__ import annotations
from typing import List, Optional

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QSpinBox, QDoubleSpinBox, QTableWidget,
    QTableWidgetItem, QMessageBox, QListWidget, QAbstractItemView,
)

from gui.app_state import AppState, AppStateError
from gui.simulation_worker import SimulationWorker


class ExecutionPanel(QWidget):
    """Emite cell_changed cuando el estado de la célula activa cambia (tras run o división)."""

    cell_changed = Signal()

    def __init__(self, app_state: AppState, parent=None) -> None:
        super().__init__(parent)
        self._app_state = app_state
        self._worker: Optional[SimulationWorker] = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QHBoxLayout(self)

        left_column = QVBoxLayout()

        # --- Controles de ejecución ---
        run_group = QGroupBox("Ejecutar")
        run_layout = QFormLayout(run_group)

        self._steps_spin = QSpinBox()
        self._steps_spin.setRange(1, 1_000_000)
        self._steps_spin.setValue(200)
        run_steps_button = QPushButton("Correr N pasos")
        run_steps_button.setObjectName("primary")
        run_steps_button.clicked.connect(self._run_steps)
        steps_row = QHBoxLayout()
        steps_row.addWidget(self._steps_spin)
        steps_row.addWidget(run_steps_button)
        run_layout.addRow("Pasos:", steps_row)

        self._seconds_spin = QDoubleSpinBox()
        self._seconds_spin.setRange(0.001, 1_000_000.0)
        self._seconds_spin.setValue(50.0)
        run_time_button = QPushButton("Correr N segundos")
        run_time_button.clicked.connect(self._run_time)
        time_row = QHBoxLayout()
        time_row.addWidget(self._seconds_spin)
        time_row.addWidget(run_time_button)
        run_layout.addRow("Tiempo simulado:", time_row)

        self._time_label = QLabel("Tiempo simulado: 0.0000")
        self._events_label = QLabel("Eventos totales: 0")
        run_layout.addRow(self._time_label)
        run_layout.addRow(self._events_label)

        left_column.addWidget(run_group)

        # --- Especies en vivo ---
        species_group = QGroupBox("Especies")
        species_layout = QVBoxLayout(species_group)
        self._species_table = QTableWidget(0, 2)
        self._species_table.setHorizontalHeaderLabels(["Especie", "Cantidad"])
        self._species_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        species_layout.addWidget(self._species_table)
        left_column.addWidget(species_group)

        # --- División celular ---
        division_group = QGroupBox("División celular")
        division_layout = QVBoxLayout(division_group)

        self._division_species_list = QListWidget()
        self._division_species_list.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        division_layout.addWidget(QLabel("Especies que suman la condición de división:"))
        division_layout.addWidget(self._division_species_list)

        threshold_row = QHBoxLayout()
        self._threshold_spin = QSpinBox()
        self._threshold_spin.setRange(1, 10_000_000)
        self._threshold_spin.setValue(30)
        threshold_row.addWidget(QLabel("Umbral:"))
        threshold_row.addWidget(self._threshold_spin)
        division_layout.addLayout(threshold_row)

        divide_button = QPushButton("Dividir célula (si se cumple la condición)")
        divide_button.setObjectName("danger")
        divide_button.clicked.connect(self._divide_cell)
        division_layout.addWidget(divide_button)

        self._division_status_label = QLabel("")
        division_layout.addWidget(self._division_status_label)

        left_column.addWidget(division_group)
        layout.addLayout(left_column, stretch=1)

        # --- Gráfica de trayectoria ---
        trajectory_group = QGroupBox("Trayectoria")
        trajectory_layout = QVBoxLayout(trajectory_group)
        self._figure = Figure(figsize=(6, 5))
        self._canvas = FigureCanvasQTAgg(self._figure)
        trajectory_layout.addWidget(self._canvas)
        refresh_plot_button = QPushButton("Actualizar gráfica")
        refresh_plot_button.clicked.connect(self.refresh)
        trajectory_layout.addWidget(refresh_plot_button)
        layout.addWidget(trajectory_group, stretch=2)

    # ------------------------------------------------------------------
    # Ejecución
    # ------------------------------------------------------------------

    def _run_steps(self) -> None:
        self._start_worker(steps=self._steps_spin.value())

    def _run_time(self) -> None:
        self._start_worker(seconds=self._seconds_spin.value())

    def _start_worker(self, steps: Optional[int] = None, seconds: Optional[float] = None) -> None:
        if self._app_state.cell is None:
            QMessageBox.warning(self, "Sin célula", "Carga e instala un genoma primero, en la pestaña Genoma.")
            return
        if self._worker is not None and self._worker.isRunning():
            QMessageBox.information(self, "Ya en marcha", "Ya hay una ejecución en curso.")
            return

        self._worker = SimulationWorker(self._app_state, steps=steps, seconds=seconds)
        self._worker.finished_run.connect(self._on_run_finished)
        self._worker.error.connect(self._on_run_error)
        self._worker.start()

    def _on_run_finished(self, n_events: int) -> None:
        self.refresh()
        self.cell_changed.emit()

    def _on_run_error(self, message: str) -> None:
        QMessageBox.critical(self, "Error durante la ejecución", message)

    # ------------------------------------------------------------------
    # Refresco de vistas
    # ------------------------------------------------------------------

    def refresh(self) -> None:
        """Actualiza tabla de especies, etiquetas y gráfica a partir del estado actual."""
        if self._app_state.cell is None:
            return

        cell = self._app_state.cell
        snapshot = cell.species.snapshot()

        self._species_table.setRowCount(len(snapshot))
        for row, (species_id, count) in enumerate(sorted(snapshot.items())):
            self._species_table.setItem(row, 0, QTableWidgetItem(species_id))
            self._species_table.setItem(row, 1, QTableWidgetItem(str(count)))

        self._division_species_list.clear()
        for species_id in sorted(snapshot.keys()):
            self._division_species_list.addItem(species_id)

        self._time_label.setText(f"Tiempo simulado: {cell.time:.4f}")
        self._events_label.setText(f"Eventos totales: {len(cell.get_events())}")

        self._refresh_trajectory_plot()

    def _refresh_trajectory_plot(self) -> None:
        try:
            df = self._app_state.trajectory()
        except AppStateError:
            return

        self._figure.clear()
        ax = self._figure.add_subplot(111)

        species_columns = [c for c in df.columns if c not in ("step", "time", "reaction_id")]
        for species_id in species_columns:
            ax.step(df["time"], df[species_id], where="post", label=species_id)

        ax.set_xlabel("Tiempo simulado")
        ax.set_ylabel("Cantidad")
        ax.set_title("Trayectoria de especies")
        if len(species_columns) <= 12:  # con demasiadas especies, la leyenda satura el gráfico
            ax.legend(fontsize="small", loc="upper left")
        self._figure.tight_layout()
        self._canvas.draw()

    # ------------------------------------------------------------------
    # División celular
    # ------------------------------------------------------------------

    def _divide_cell(self) -> None:
        selected_items = self._division_species_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "Sin especies", "Selecciona al menos una especie para la condición de división.")
            return

        threshold_species_ids: List[str] = [item.text() for item in selected_items]
        gene_species_ids = [s for s in self._app_state.species_snapshot().keys() if s.startswith("gene_")]

        try:
            daughter_a, daughter_b = self._app_state.divide(
                threshold_species_ids, self._threshold_spin.value(), gene_species_ids
            )
        except AppStateError as exc:
            self._division_status_label.setText(f"⚠️ {exc}")
            return

        self._app_state.adopt_daughter(daughter_a)
        self._division_status_label.setText(
            f"✅ División realizada (generación {self._app_state.generation}). "
            f"Continuando con una de las dos hijas; la otra queda descartada de esta vista "
            f"(no hay soporte de multi-célula en esta interfaz todavía)."
        )
        self.refresh()
        self.cell_changed.emit()
