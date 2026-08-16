"""
gui/forensics_panel.py

Panel de análisis forense (pestaña "Forense"): interfaz sobre
engine/forensics.py. El usuario elige una especie y una ventana
temporal, y se muestra qué reacciones contribuyeron más a su cambio
neto en ese intervalo — el "¿por qué pasó?" en contraste con el
"¿qué va a pasar?" del depurador paso a paso.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QLineEdit, QDoubleSpinBox, QTableWidget,
    QTableWidgetItem, QMessageBox, QAbstractItemView,
)

from gui.app_state import AppState, AppStateError


class ForensicsPanel(QWidget):
    def __init__(self, app_state: AppState, parent=None) -> None:
        super().__init__(parent)
        self._app_state = app_state
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)

        form_group = QGroupBox("Análisis causal retrospectivo")
        form_layout = QFormLayout(form_group)

        self._species_edit = QLineEdit()
        self._species_edit.setPlaceholderText("p. ej. ATP, mRNA_target")
        form_layout.addRow("Especie:", self._species_edit)

        self._window_start_spin = QDoubleSpinBox()
        self._window_start_spin.setRange(0.0, 1_000_000_000.0)
        self._window_start_spin.setValue(0.0)
        form_layout.addRow("Desde (tiempo):", self._window_start_spin)

        self._window_end_spin = QDoubleSpinBox()
        self._window_end_spin.setRange(0.0, 1_000_000_000.0)
        self._window_end_spin.setValue(1_000_000.0)
        form_layout.addRow("Hasta (tiempo):", self._window_end_spin)

        use_full_range_button = QPushButton("Usar todo el rango simulado")
        use_full_range_button.clicked.connect(self._use_full_range)
        form_layout.addRow(use_full_range_button)

        analyze_button = QPushButton("Analizar")
        analyze_button.setObjectName("primary")
        analyze_button.clicked.connect(self._analyze)
        form_layout.addRow(analyze_button)

        layout.addWidget(form_group)

        self._summary_label = QLabel("")
        self._summary_label.setWordWrap(True)
        layout.addWidget(self._summary_label)

        self._table = QTableWidget(0, 3)
        self._table.setHorizontalHeaderLabels(["Reacción", "Veces disparada", "Cambio neto"])
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self._table)

    def _use_full_range(self) -> None:
        if self._app_state.cell is None:
            return
        self._window_start_spin.setValue(0.0)
        self._window_end_spin.setValue(self._app_state.cell.time)

    def _analyze(self) -> None:
        species_id = self._species_edit.text().strip()
        if not species_id:
            QMessageBox.warning(self, "Falta la especie", "Escribe el id de la especie a analizar.")
            return

        try:
            report = self._app_state.analyze_species(
                species_id, self._window_start_spin.value(), self._window_end_spin.value()
            )
        except AppStateError as exc:
            QMessageBox.warning(self, "No se pudo analizar", str(exc))
            return

        producer = report.top_producer()
        consumer = report.top_consumer()
        summary_parts = [f"Cambio neto de '{species_id}' en la ventana: {report.net_change:+d}."]
        if producer is not None:
            summary_parts.append(
                f"Mayor productor: '{producer.reaction_id}' ({producer.net_change:+d}, "
                f"disparó {producer.times_fired} veces)."
            )
        if consumer is not None:
            summary_parts.append(
                f"Mayor consumidor: '{consumer.reaction_id}' ({consumer.net_change:+d}, "
                f"disparó {consumer.times_fired} veces)."
            )
        if producer is None and consumer is None:
            summary_parts.append("Ninguna reacción de la ventana afectó a esta especie.")
        self._summary_label.setText("  ".join(summary_parts))

        self._table.setRowCount(len(report.contributions))
        for row, contribution in enumerate(report.contributions):
            self._table.setItem(row, 0, QTableWidgetItem(contribution.reaction_id))
            self._table.setItem(row, 1, QTableWidgetItem(str(contribution.times_fired)))
            self._table.setItem(row, 2, QTableWidgetItem(f"{contribution.net_change:+d}"))
