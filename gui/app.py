"""
gui/app.py

Punto de entrada de la aplicación GUI. Ejecutar con:

    python -m gui.app

desde la raíz del proyecto (o `python gui/app.py` si el directorio
raíz ya está en PYTHONPATH).
"""

from __future__ import annotations
import sys

from PySide6.QtWidgets import QApplication

from gui.theme import apply_dark_theme
from gui.main_window import MainWindow


def main() -> int:
    app = QApplication(sys.argv)
    apply_dark_theme(app)

    window = MainWindow()
    window.show()

    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
