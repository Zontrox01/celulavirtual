"""
gui/debugger_panel.py

Panel del depurador (pestaña "Depurador"): el diferenciador del
proyecto, expuesto en la interfaz. Ejecución paso a paso, breakpoints
condicionales, visualización del espacio de propensidades
(engine/propensity_view.py), y retroceso (engine/snapshot.py vía
Debugger.undo_to()).

El Debugger se crea una vez (al entrar en la pestaña o al pulsar
"Nuevo depurador") y se reutiliza mientras la célula no cambie de
identidad (p. ej. tras una división, hay que crear uno nuevo).
"""

from __future__ import annotations
from typing import Optional

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QPushButton, QLabel, QSpinBox, QDoubleSpinBox, QComboBox,
    QLineEdit, QListWidget, QMessageBox, QCheckBox,
)

from engine.debugger import DebuggerEndedError
from engine.propensity_view import prepare_propensity_bars
from gui.app_state import AppState
from gui.breakpoint_builder import build_breakpoint_condition, describe_breakpoint, OPERATORS


class DebuggerPanel(QWidget):
    """Emite cell_changed cuando avanzar/retroceder cambia el estado de la célula activa."""

    cell_changed = Signal()

    def __init__(self, app_state: AppState, parent=None) -> None:
        super().__init__(parent)
        self._app_state = app_state
        self._debugger = None  # type: Optional["Debugger"]
        self._breakpoint_conditions: dict = {}  # descripción -> callable
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QHBoxLayout(self)
        left_column = QVBoxLayout()

        # --- Depurador ---
        session_group = QGroupBox("Sesión de depuración")
        session_layout = QVBoxLayout(session_group)

        new_debugger_row = QHBoxLayout()
        self._enable_snapshots_check = QCheckBox("Habilitar retroceso (snapshots)")
        self._enable_snapshots_check.setChecked(True)
        self._snapshot_interval_spin = QSpinBox()
        self._snapshot_interval_spin.setRange(1, 10_000)
        self._snapshot_interval_spin.setValue(20)
        new_debugger_button = QPushButton("Nuevo depurador sobre la célula actual")
        new_debugger_button.setObjectName("primary")
        new_debugger_button.clicked.connect(self._new_debugger)
        new_debugger_row.addWidget(self._enable_snapshots_check)
        new_debugger_row.addWidget(QLabel("Intervalo:"))
        new_debugger_row.addWidget(self._snapshot_interval_spin)
        session_layout.addLayout(new_debugger_row)
        session_layout.addWidget(new_debugger_button)

        self._session_status_label = QLabel("Sin sesión de depuración activa.")
        session_layout.addWidget(self._session_status_label)

        left_column.addWidget(session_group)

        # --- Paso a paso ---
        step_group = QGroupBox("Paso a paso")
        step_layout = QVBoxLayout(step_group)

        step_row = QHBoxLayout()
        self._step_n_spin = QSpinBox()
        self._step_n_spin.setRange(1, 100_000)
        self._step_n_spin.setValue(1)
        step_button = QPushButton("Avanzar N pasos")
        step_button.clicked.connect(self._step)
        step_row.addWidget(self._step_n_spin)
        step_row.addWidget(step_button)
        step_layout.addLayout(step_row)

        run_until_bp_button = QPushButton("Avanzar hasta el próximo breakpoint")
        run_until_bp_button.clicked.connect(self._run_until_breakpoint)
        step_layout.addWidget(run_until_bp_button)

        self._current_event_label = QLabel("Último evento: (ninguno todavía)")
        self._current_event_label.setWordWrap(True)
        step_layout.addWidget(self._current_event_label)

        left_column.addWidget(step_group)

        # --- Breakpoints ---
        breakpoints_group = QGroupBox("Breakpoints")
        breakpoints_layout = QVBoxLayout(breakpoints_group)

        bp_form = QHBoxLayout()
        self._bp_species_edit = QLineEdit()
        self._bp_species_edit.setPlaceholderText("especie, p. ej. ATP")
        self._bp_operator_combo = QComboBox()
        self._bp_operator_combo.addItems(sorted(OPERATORS.keys()))
        self._bp_threshold_spin = QDoubleSpinBox()
        self._bp_threshold_spin.setRange(-1_000_000_000.0, 1_000_000_000.0)
        add_bp_button = QPushButton("Añadir")
        add_bp_button.clicked.connect(self._add_breakpoint)
        bp_form.addWidget(self._bp_species_edit)
        bp_form.addWidget(self._bp_operator_combo)
        bp_form.addWidget(self._bp_threshold_spin)
        bp_form.addWidget(add_bp_button)
        breakpoints_layout.addLayout(bp_form)

        self._breakpoints_list = QListWidget()
        breakpoints_layout.addWidget(self._breakpoints_list)

        remove_bp_button = QPushButton("Quitar seleccionado")
        remove_bp_button.clicked.connect(self._remove_selected_breakpoint)
        breakpoints_layout.addWidget(remove_bp_button)

        left_column.addWidget(breakpoints_group)

        # --- Retroceso ---
        undo_group = QGroupBox("Retroceso (undo)")
        undo_layout = QHBoxLayout(undo_group)
        self._undo_step_spin = QSpinBox()
        self._undo_step_spin.setRange(0, 10_000_000)
        undo_button = QPushButton("Retroceder al paso")
        undo_button.clicked.connect(self._undo_to)
        undo_layout.addWidget(self._undo_step_spin)
        undo_layout.addWidget(undo_button)
        left_column.addWidget(undo_group)

        layout.addLayout(left_column, stretch=1)

        # --- Espacio de propensidades ---
        propensity_group = QGroupBox("Espacio de propensidades")
        propensity_layout = QVBoxLayout(propensity_group)
        self._figure = Figure(figsize=(5, 5))
        self._canvas = FigureCanvasQTAgg(self._figure)
        propensity_layout.addWidget(self._canvas)
        layout.addWidget(propensity_group, stretch=1)

    # ------------------------------------------------------------------
    # Sesión
    # ------------------------------------------------------------------

    def _new_debugger(self) -> None:
        if self._app_state.cell is None:
            QMessageBox.warning(self, "Sin célula", "Carga e instala un genoma primero.")
            return

        self._debugger = self._app_state.new_debugger(
            enable_snapshots=self._enable_snapshots_check.isChecked(),
            snapshot_interval=self._snapshot_interval_spin.value(),
        )
        self._breakpoint_conditions.clear()
        self._breakpoints_list.clear()
        self._session_status_label.setText(
            f"Sesión activa (snapshots {'activados' if self._debugger.snapshots_enabled else 'desactivados'})."
        )
        self._refresh_propensity_plot()

    def _require_debugger(self):
        if self._debugger is None:
            QMessageBox.warning(self, "Sin sesión", "Crea un depurador nuevo primero.")
            return None
        return self._debugger

    # ------------------------------------------------------------------
    # Paso a paso
    # ------------------------------------------------------------------

    def _step(self) -> None:
        debugger = self._require_debugger()
        if debugger is None:
            return
        try:
            events = debugger.step(n=self._step_n_spin.value())
        except DebuggerEndedError:
            QMessageBox.information(self, "Simulación terminada", "No quedan más reacciones posibles.")
            return

        if events:
            last = events[-1]
            self._current_event_label.setText(
                f"Último evento: paso {last.step}, t={last.time:.4f}, reacción '{last.reaction_id}'"
            )
        if debugger.last_triggered_breakpoint is not None:
            self._current_event_label.setText(
                self._current_event_label.text() + f"  ⏸ breakpoint '{debugger.last_triggered_breakpoint}'"
            )

        self._refresh_propensity_plot()
        self.cell_changed.emit()

    def _run_until_breakpoint(self) -> None:
        debugger = self._require_debugger()
        if debugger is None:
            return
        if not self._breakpoint_conditions:
            QMessageBox.information(self, "Sin breakpoints", "Añade al menos un breakpoint primero.")
            return

        events = debugger.run_until_breakpoint(max_steps=200_000)
        if events:
            last = events[-1]
            triggered = debugger.last_triggered_breakpoint or "(fin de la simulación)"
            self._current_event_label.setText(
                f"Detenido en el paso {last.step}, t={last.time:.4f} — breakpoint: {triggered}"
            )
        self._refresh_propensity_plot()
        self.cell_changed.emit()

    # ------------------------------------------------------------------
    # Breakpoints
    # ------------------------------------------------------------------

    def _add_breakpoint(self) -> None:
        debugger = self._require_debugger()
        if debugger is None:
            return

        species_id = self._bp_species_edit.text().strip()
        operator = self._bp_operator_combo.currentText()
        threshold = self._bp_threshold_spin.value()

        try:
            condition = build_breakpoint_condition(species_id, operator, threshold)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Breakpoint inválido", str(exc))
            return

        description = describe_breakpoint(species_id, operator, threshold)
        debugger.set_breakpoint(description, condition)
        self._breakpoint_conditions[description] = condition
        self._breakpoints_list.addItem(description)

    def _remove_selected_breakpoint(self) -> None:
        debugger = self._require_debugger()
        if debugger is None:
            return
        for item in self._breakpoints_list.selectedItems():
            description = item.text()
            debugger.remove_breakpoint(description)
            self._breakpoint_conditions.pop(description, None)
            self._breakpoints_list.takeItem(self._breakpoints_list.row(item))

    # ------------------------------------------------------------------
    # Retroceso
    # ------------------------------------------------------------------

    def _undo_to(self) -> None:
        debugger = self._require_debugger()
        if debugger is None:
            return
        try:
            debugger.undo_to(self._undo_step_spin.value())
        except DebuggerEndedError as exc:
            QMessageBox.warning(self, "Retroceso no disponible", str(exc))
            return
        except Exception as exc:  # noqa: BLE001 — p. ej. SnapshotError si el paso pedido no existe
            QMessageBox.warning(self, "No se pudo retroceder", str(exc))
            return

        self._current_event_label.setText(f"Retrocedido al paso {self._undo_step_spin.value()}.")
        self._refresh_propensity_plot()
        self.cell_changed.emit()

    # ------------------------------------------------------------------
    # Visualización de propensidades
    # ------------------------------------------------------------------

    def _refresh_propensity_plot(self) -> None:
        debugger = self._debugger
        if debugger is None:
            return

        propensities = debugger.pending_propensities()
        bars = prepare_propensity_bars(propensities)

        self._figure.clear()
        ax = self._figure.add_subplot(111)

        if bars:
            labels = [b.reaction_id for b in bars]
            values = [b.propensity for b in bars]
            colors = [
                "#58a6ff" if debugger.current_event and b.reaction_id == debugger.current_event.reaction_id
                else "#4a4a52"
                for b in bars
            ]
            ax.barh(labels, values, color=colors)
            ax.invert_yaxis()  # la de mayor propensidad arriba
            ax.set_xlabel("Propensidad")
        else:
            ax.text(0.5, 0.5, "(sin reacciones activas)", ha="center", va="center")

        ax.set_title("Espacio de propensidades")
        self._figure.tight_layout()
        self._canvas.draw()
