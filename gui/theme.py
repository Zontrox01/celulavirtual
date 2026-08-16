"""
gui/theme.py

Tema oscuro de la aplicación. Separado en su propio archivo para que
main_window.py no se llene de valores de color hardcodeados.
"""

from __future__ import annotations

from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication


def dark_palette() -> QPalette:
    """Paleta de colores oscura, aplicada a nivel de QApplication."""
    palette = QPalette()

    background = QColor(30, 30, 34)
    surface = QColor(40, 40, 46)
    text = QColor(220, 220, 225)
    disabled_text = QColor(120, 120, 126)
    accent = QColor(88, 166, 255)  # azul, para botones/enlaces resaltados
    highlight = QColor(61, 90, 128)

    palette.setColor(QPalette.ColorRole.Window, background)
    palette.setColor(QPalette.ColorRole.WindowText, text)
    palette.setColor(QPalette.ColorRole.Base, surface)
    palette.setColor(QPalette.ColorRole.AlternateBase, background)
    palette.setColor(QPalette.ColorRole.ToolTipBase, surface)
    palette.setColor(QPalette.ColorRole.ToolTipText, text)
    palette.setColor(QPalette.ColorRole.Text, text)
    palette.setColor(QPalette.ColorRole.Button, surface)
    palette.setColor(QPalette.ColorRole.ButtonText, text)
    palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 90, 90))
    palette.setColor(QPalette.ColorRole.Link, accent)
    palette.setColor(QPalette.ColorRole.Highlight, highlight)
    palette.setColor(QPalette.ColorRole.HighlightedText, text)

    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText, disabled_text)
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, disabled_text)
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, disabled_text)

    return palette


DARK_STYLESHEET = """
QToolTip {
    color: #dcdce1;
    background-color: #28282e;
    border: 1px solid #4a4a52;
    padding: 4px;
}
QGroupBox {
    border: 1px solid #4a4a52;
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 8px;
    font-weight: 600;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
}
QPushButton {
    background-color: #3a3a42;
    border: 1px solid #55555f;
    border-radius: 4px;
    padding: 6px 14px;
}
QPushButton:hover {
    background-color: #46464f;
}
QPushButton:pressed {
    background-color: #2c2c33;
}
QPushButton:disabled {
    color: #78787e;
    background-color: #2c2c33;
}
QPushButton#primary {
    background-color: #2f5fa8;
    border: 1px solid #3d6bb5;
}
QPushButton#primary:hover {
    background-color: #3a6bb8;
}
QTableWidget, QListWidget {
    background-color: #26262b;
    alternate-background-color: #2b2b31;
    gridline-color: #3c3c44;
    border: 1px solid #4a4a52;
}
QHeaderView::section {
    background-color: #33333a;
    padding: 4px;
    border: none;
}
QTabWidget::pane {
    border: 1px solid #4a4a52;
    border-radius: 4px;
}
QTabBar::tab {
    background: #2c2c33;
    padding: 8px 16px;
    border: 1px solid #4a4a52;
    border-bottom: none;
}
QTabBar::tab:selected {
    background: #3a3a42;
}
QProgressBar {
    border: 1px solid #4a4a52;
    border-radius: 4px;
    text-align: center;
}
QProgressBar::chunk {
    background-color: #2f5fa8;
}
"""


def apply_dark_theme(app: QApplication) -> None:
    """Aplica el tema oscuro (paleta + hoja de estilos) a toda la aplicación."""
    app.setStyle("Fusion")
    app.setPalette(dark_palette())
    app.setStyleSheet(DARK_STYLESHEET)
