"""
gui/main_window.py

Ventana principal de la aplicación: una sola ventana grande con
pestañas para cada área funcional (Genoma, Ejecución, Depurador,
Forense), más un menú "Archivo" para exportar/importar SBML — tal
como se pidió: toda la funcionalidad accesible de una vez, en vez de
repartida en varias ventanas.

No contiene lógica de negocio propia: crea el AppState compartido y
lo pasa a cada panel; cada panel es responsable de su propia
interacción con él.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QMainWindow, QTabWidget, QFileDialog, QMessageBox,
)

from gui.app_state import AppState, AppStateError
from gui.genome_setup_panel import GenomeSetupPanel
from gui.execution_panel import ExecutionPanel
from gui.debugger_panel import DebuggerPanel
from gui.forensics_panel import ForensicsPanel


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Célula Virtual")
        self.resize(1200, 800)

        self._app_state = AppState()

        self._tabs = QTabWidget()
        self.setCentralWidget(self._tabs)

        self._genome_panel = GenomeSetupPanel(self._app_state)
        self._execution_panel = ExecutionPanel(self._app_state)
        self._debugger_panel = DebuggerPanel(self._app_state)
        self._forensics_panel = ForensicsPanel(self._app_state)

        self._tabs.addTab(self._genome_panel, "Genoma")
        self._tabs.addTab(self._execution_panel, "Ejecución")
        self._tabs.addTab(self._debugger_panel, "Depurador")
        self._tabs.addTab(self._forensics_panel, "Forense")

        # Cuando el genoma se instala, o cuando run()/step()/undo()/división
        # cambian el estado, se refresca la tabla y la gráfica de la pestaña
        # de ejecución — sin esto, el usuario tendría que cambiar de pestaña
        # a mano para ver los resultados de lo que hizo en otra.
        self._genome_panel.genome_installed.connect(self._execution_panel.refresh)
        self._debugger_panel.cell_changed.connect(self._execution_panel.refresh)
        self._execution_panel.cell_changed.connect(self._execution_panel.refresh)

        self._build_menu()

    def _build_menu(self) -> None:
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("&Archivo")

        export_action = file_menu.addAction("Exportar a SBML…")
        export_action.triggered.connect(self._export_sbml)

        import_action = file_menu.addAction("Importar desde SBML…")
        import_action.triggered.connect(self._import_sbml)

    def _export_sbml(self) -> None:
        if self._app_state.cell is None:
            QMessageBox.warning(self, "Sin célula", "Carga e instala un genoma primero.")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Exportar a SBML", "modelo.xml", "SBML (*.xml)")
        if not path:
            return

        try:
            self._app_state.export_sbml(path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Error al exportar", str(exc))
            return

        QMessageBox.information(self, "Exportado", f"Modelo exportado a:\n{path}")

    def _import_sbml(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Importar desde SBML", "", "SBML (*.xml)")
        if not path:
            return

        try:
            self._app_state.import_sbml(path)
        except AppStateError as exc:
            QMessageBox.critical(self, "Error al importar", str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Error inesperado", str(exc))
            return

        QMessageBox.information(
            self,
            "Importado",
            "Modelo importado. Nota: el SBML no conserva a qué módulo biológico "
            "pertenecía cada reacción, así que la pestaña Genoma no reflejará "
            "esta célula — pero Ejecución, Depurador y Forense sí funcionan con ella.",
        )
        self._execution_panel.refresh()
