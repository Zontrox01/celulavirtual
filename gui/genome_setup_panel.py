"""
gui/genome_setup_panel.py

Panel de carga e instalación del genoma (pestaña "Genoma"). Deja que
el usuario elija los archivos FASTA + anotación, configure las
cantidades iniciales de maquinaria (RNAP, ribosomas), las vidas medias
de degradación, y opcionalmente el metabolismo — y llama a
gui/app_state.py para construir la célula. No contiene lógica de
negocio propia.
"""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QLineEdit, QSpinBox, QDoubleSpinBox,
    QComboBox, QFileDialog, QMessageBox,
)

from gui.app_state import AppState, AppStateError, GenomeSetupOptions


class GenomeSetupPanel(QWidget):
    """Emite genome_installed cuando la célula se ha construido con éxito."""

    genome_installed = Signal()

    def __init__(self, app_state: AppState, parent=None) -> None:
        super().__init__(parent)
        self._app_state = app_state
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)

        # --- Archivos ---
        files_group = QGroupBox("Archivos del genoma")
        files_layout = QFormLayout(files_group)

        self._fasta_edit = QLineEdit()
        fasta_row = QHBoxLayout()
        fasta_row.addWidget(self._fasta_edit)
        fasta_button = QPushButton("Examinar…")
        fasta_button.clicked.connect(self._browse_fasta)
        fasta_row.addWidget(fasta_button)
        files_layout.addRow("FASTA:", fasta_row)

        self._annotation_edit = QLineEdit()
        annotation_row = QHBoxLayout()
        annotation_row.addWidget(self._annotation_edit)
        annotation_button = QPushButton("Examinar…")
        annotation_button.clicked.connect(self._browse_annotation)
        annotation_row.addWidget(annotation_button)
        files_layout.addRow("Anotación (CSV):", annotation_row)

        layout.addWidget(files_group)

        # --- Maquinaria y degradación ---
        machinery_group = QGroupBox("Maquinaria y degradación")
        machinery_layout = QFormLayout(machinery_group)

        self._seed_spin = QSpinBox()
        self._seed_spin.setRange(0, 2_147_483_647)
        self._seed_spin.setValue(1)
        machinery_layout.addRow("Semilla aleatoria:", self._seed_spin)

        self._rnap_spin = QSpinBox()
        self._rnap_spin.setRange(1, 100_000)
        self._rnap_spin.setValue(30)
        machinery_layout.addRow("RNAP inicial:", self._rnap_spin)

        self._ribosome_spin = QSpinBox()
        self._ribosome_spin.setRange(1, 100_000)
        self._ribosome_spin.setValue(20)
        machinery_layout.addRow("Ribosomas inicial:", self._ribosome_spin)

        self._mrna_half_life_spin = QDoubleSpinBox()
        self._mrna_half_life_spin.setRange(0.1, 1_000_000.0)
        self._mrna_half_life_spin.setValue(300.0)
        machinery_layout.addRow("Vida media ARNm (s):", self._mrna_half_life_spin)

        self._protein_half_life_spin = QDoubleSpinBox()
        self._protein_half_life_spin.setRange(0.1, 1_000_000.0)
        self._protein_half_life_spin.setValue(3600.0)
        machinery_layout.addRow("Vida media proteína (s):", self._protein_half_life_spin)

        layout.addWidget(machinery_group)

        # --- Metabolismo ---
        metabolism_group = QGroupBox("Metabolismo (opcional)")
        metabolism_layout = QFormLayout(metabolism_group)

        self._metabolism_combo = QComboBox()
        self._metabolism_combo.addItem("Sin metabolismo", "none")
        self._metabolism_combo.addItem("De fondo, con acoplamiento a ATP", "background_with_atp_coupling")
        self._metabolism_combo.addItem("Catalizado por una enzima (gen concreto)", "enzyme_catalyzed")
        self._metabolism_combo.currentIndexChanged.connect(self._on_metabolism_mode_changed)
        metabolism_layout.addRow("Modo:", self._metabolism_combo)

        self._enzyme_gene_edit = QLineEdit()
        self._enzyme_gene_edit.setPlaceholderText("id del gen, p. ej. metA")
        self._enzyme_gene_edit.setEnabled(False)
        metabolism_layout.addRow("Gen de la enzima:", self._enzyme_gene_edit)

        self._atp_cost_transcription_spin = QSpinBox()
        self._atp_cost_transcription_spin.setRange(0, 1000)
        self._atp_cost_transcription_spin.setValue(2)
        self._atp_cost_transcription_spin.setEnabled(False)
        metabolism_layout.addRow("Coste ATP / transcripción:", self._atp_cost_transcription_spin)

        self._atp_cost_translation_spin = QSpinBox()
        self._atp_cost_translation_spin.setRange(0, 1000)
        self._atp_cost_translation_spin.setValue(4)
        self._atp_cost_translation_spin.setEnabled(False)
        metabolism_layout.addRow("Coste ATP / traducción:", self._atp_cost_translation_spin)

        layout.addWidget(metabolism_group)

        # --- Instalar ---
        install_button = QPushButton("Cargar e instalar genoma")
        install_button.setObjectName("primary")
        install_button.clicked.connect(self._install_genome)
        layout.addWidget(install_button)

        self._status_label = QLabel("")
        layout.addWidget(self._status_label)

        layout.addStretch(1)

    def _on_metabolism_mode_changed(self) -> None:
        mode = self._metabolism_combo.currentData()
        self._enzyme_gene_edit.setEnabled(mode == "enzyme_catalyzed")
        self._atp_cost_transcription_spin.setEnabled(mode == "background_with_atp_coupling")
        self._atp_cost_translation_spin.setEnabled(mode == "background_with_atp_coupling")

    def _browse_fasta(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Elegir archivo FASTA", "", "FASTA (*.fasta *.fa);;Todos (*)")
        if path:
            self._fasta_edit.setText(path)

    def _browse_annotation(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Elegir anotación", "", "CSV (*.csv);;Todos (*)")
        if path:
            self._annotation_edit.setText(path)

    def _install_genome(self) -> None:
        if not self._fasta_edit.text() or not self._annotation_edit.text():
            QMessageBox.warning(self, "Faltan archivos", "Elige el FASTA y la anotación antes de instalar.")
            return

        options = GenomeSetupOptions(
            fasta_path=self._fasta_edit.text(),
            annotation_path=self._annotation_edit.text(),
            seed=self._seed_spin.value(),
            rnap_initial_count=self._rnap_spin.value(),
            ribosome_initial_count=self._ribosome_spin.value(),
            mrna_half_life=self._mrna_half_life_spin.value(),
            protein_half_life=self._protein_half_life_spin.value(),
            metabolism_mode=self._metabolism_combo.currentData(),
            enzyme_gene_id=self._enzyme_gene_edit.text() or None,
            atp_cost_per_transcription=self._atp_cost_transcription_spin.value(),
            atp_cost_per_translation=self._atp_cost_translation_spin.value(),
        )

        try:
            self._app_state.load_and_install_genome(options)
        except AppStateError as exc:
            QMessageBox.critical(self, "Error al instalar el genoma", str(exc))
            return
        except Exception as exc:  # noqa: BLE001 — cualquier otro fallo (archivo mal formado, etc.) también se muestra
            QMessageBox.critical(self, "Error inesperado", str(exc))
            return

        n_genes = len(self._app_state.installed_genes)
        self._status_label.setText(
            f"✅ Genoma instalado: {n_genes} genes, {len(self._app_state.cell.reactions)} reacciones."
        )
        self.genome_installed.emit()
