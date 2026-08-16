"""
gui/theme.py

Tema oscuro de la aplicación, inspirado en la paleta Catppuccin Mocha
(fondo azulado muy oscuro, acentos de color saturado por función:
azul para acciones principales, verde para "avanzar", rosa/rojo para
detener o deshacer). Separado en su propio archivo para que el resto
de la GUI no se llene de valores de color hardcodeados.

Uso de las variantes de botón (vía objectName, ver *_panel.py):
  - sin objectName / "default": botón neutro
  - "primary": acción principal de la pestaña (instalar, correr, analizar)
  - "success": acción de avance/progreso (paso a paso)
  - "danger": acción de detener o deshacer

No modifica ningún archivo fuera de gui/: es puramente estético.
"""

from __future__ import annotations

from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication

# --- Paleta base ---
BG = "#1e1e2e"          # fondo de la ventana
SURFACE = "#232334"      # fondo de paneles/tarjetas
SURFACE_ALT = "#292a3a"  # fondo alterno (filas de tabla)
INPUT_BG = "#313244"     # fondo de campos de entrada
BORDER = "#45475a"       # bordes
BORDER_LIGHT = "#585b70"  # bordes en hover
TEXT = "#cdd6f4"          # texto principal
TEXT_MUTED = "#7f849c"    # texto secundario/deshabilitado

BLUE = "#89b4fa"     # acento principal
BLUE_HOVER = "#74c7ec"
GREEN = "#a6e3a1"    # avance / éxito
GREEN_HOVER = "#94e2d5"
PINK = "#f38ba8"     # detener / peligro
PINK_HOVER = "#eba0ac"
YELLOW = "#f9e2af"   # aviso
PURPLE = "#cba6f7"   # decorativo / resaltes secundarios


def dark_palette() -> QPalette:
    """Paleta de colores oscura, aplicada a nivel de QApplication."""
    palette = QPalette()

    palette.setColor(QPalette.ColorRole.Window, QColor(BG))
    palette.setColor(QPalette.ColorRole.WindowText, QColor(TEXT))
    palette.setColor(QPalette.ColorRole.Base, QColor(INPUT_BG))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(SURFACE_ALT))
    palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(SURFACE))
    palette.setColor(QPalette.ColorRole.ToolTipText, QColor(TEXT))
    palette.setColor(QPalette.ColorRole.Text, QColor(TEXT))
    palette.setColor(QPalette.ColorRole.Button, QColor(SURFACE))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor(TEXT))
    palette.setColor(QPalette.ColorRole.BrightText, QColor(PINK))
    palette.setColor(QPalette.ColorRole.Link, QColor(BLUE))
    palette.setColor(QPalette.ColorRole.Highlight, QColor(BLUE))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor(BG))

    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText, QColor(TEXT_MUTED))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, QColor(TEXT_MUTED))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, QColor(TEXT_MUTED))

    return palette


DARK_STYLESHEET = f"""
QMainWindow, QWidget {{
    background-color: {BG};
    color: {TEXT};
}}

QToolTip {{
    color: {TEXT};
    background-color: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 4px;
    padding: 4px 6px;
}}

QLabel {{
    color: {TEXT};
}}

QGroupBox {{
    background-color: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 8px;
    margin-top: 14px;
    padding: 12px 10px 10px 10px;
    font-weight: 600;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: 12px;
    top: -2px;
    padding: 0 6px;
    color: {BLUE};
}}

QLineEdit, QTextEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background-color: {INPUT_BG};
    color: {TEXT};
    border: 1px solid {BORDER};
    border-radius: 5px;
    padding: 5px 8px;
    selection-background-color: {BLUE};
    selection-color: {BG};
}}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
    border: 1px solid {BLUE};
}}
QComboBox::drop-down {{
    border: none;
    width: 22px;
}}

QPushButton {{
    background-color: {INPUT_BG};
    color: {TEXT};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 7px 16px;
    font-weight: 600;
}}
QPushButton:hover {{
    background-color: {BORDER_LIGHT};
}}
QPushButton:pressed {{
    background-color: {SURFACE};
}}
QPushButton:disabled {{
    color: {TEXT_MUTED};
    background-color: {SURFACE};
    border: 1px solid {SURFACE_ALT};
}}

/* Variantes semánticas: acción principal / avance / detener-deshacer */
QPushButton#primary {{
    background-color: {BLUE};
    color: {BG};
    border: 1px solid {BLUE};
}}
QPushButton#primary:hover {{
    background-color: {BLUE_HOVER};
}}
QPushButton#success {{
    background-color: {GREEN};
    color: {BG};
    border: 1px solid {GREEN};
}}
QPushButton#success:hover {{
    background-color: {GREEN_HOVER};
}}
QPushButton#danger {{
    background-color: {PINK};
    color: {BG};
    border: 1px solid {PINK};
}}
QPushButton#danger:hover {{
    background-color: {PINK_HOVER};
}}

QTableWidget, QListWidget {{
    background-color: {INPUT_BG};
    alternate-background-color: {SURFACE_ALT};
    gridline-color: {BORDER};
    border: 1px solid {BORDER};
    border-radius: 6px;
}}
QListWidget::item, QTableWidget::item {{
    padding: 4px;
}}
QListWidget::item:selected, QTableWidget::item:selected {{
    background-color: {BLUE};
    color: {BG};
}}
QHeaderView::section {{
    background-color: {SURFACE};
    color: {TEXT};
    padding: 6px;
    border: none;
    border-bottom: 1px solid {BORDER};
}}

QTabWidget::pane {{
    border: 1px solid {BORDER};
    border-radius: 8px;
    background-color: {SURFACE};
    top: -1px;
}}
QTabBar::tab {{
    background: {INPUT_BG};
    color: {TEXT_MUTED};
    padding: 9px 20px;
    border: 1px solid {BORDER};
    border-bottom: none;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    margin-right: 2px;
}}
QTabBar::tab:selected {{
    background: {SURFACE};
    color: {BLUE};
    font-weight: 600;
}}
QTabBar::tab:hover:!selected {{
    color: {TEXT};
}}

QProgressBar {{
    background-color: {INPUT_BG};
    border: 1px solid {BORDER};
    border-radius: 6px;
    text-align: center;
    color: {TEXT};
}}
QProgressBar::chunk {{
    background-color: {BLUE};
    border-radius: 5px;
}}

QMenuBar {{
    background-color: {SURFACE};
    color: {TEXT};
}}
QMenuBar::item:selected {{
    background-color: {INPUT_BG};
}}
QMenu {{
    background-color: {SURFACE};
    color: {TEXT};
    border: 1px solid {BORDER};
}}
QMenu::item:selected {{
    background-color: {BLUE};
    color: {BG};
}}

QScrollBar:vertical {{
    background: {BG};
    width: 12px;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {BORDER_LIGHT};
    border-radius: 5px;
    min-height: 24px;
}}
QScrollBar::handle:vertical:hover {{
    background: {BLUE};
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}

QCheckBox {{
    color: {TEXT};
    spacing: 8px;
}}
"""


def apply_dark_theme(app: QApplication) -> None:
    """Aplica el tema oscuro (paleta + hoja de estilos) a toda la aplicación."""
    app.setStyle("Fusion")
    app.setPalette(dark_palette())
    app.setStyleSheet(DARK_STYLESHEET)
