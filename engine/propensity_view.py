"""
engine/propensity_view.py

Visualización en tiempo real del espacio de propensidades (whitepaper,
secciones 4.4.a y 5.1). Convierte las propensidades de todas las
reacciones candidatas — que ya expone `Debugger.pending_propensities()`
— en una representación visual, para poder VER qué reacción "va
ganando" en cada instante, no solo leerlo en una tabla de números.

La lógica de preparación de datos (ordenar, calcular pesos relativos)
está deliberadamente separada de la renderización con `rich`, para
poder testearla sin necesitar una terminal real ni esa librería
instalada. Si `rich` no está disponible, se cae a una versión en texto
plano con el mismo contenido informativo.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass(frozen=True)
class PropensityBar:
    """Una fila del gráfico: una reacción, su propensidad, y su peso relativo."""

    reaction_id: str
    propensity: float
    fraction: float  # propensity / total, en [0, 1]


def prepare_propensity_bars(propensities: Dict[str, float]) -> List[PropensityBar]:
    """
    Ordena las reacciones por propensidad descendente y calcula el
    peso relativo de cada una. Si la propensidad total es 0 (ninguna
    reacción puede ocurrir), todas las fracciones son 0 también, en
    vez de dividir por cero.
    """
    total = sum(propensities.values())
    bars = [
        PropensityBar(
            reaction_id=reaction_id,
            propensity=propensity,
            fraction=(propensity / total) if total > 0 else 0.0,
        )
        for reaction_id, propensity in propensities.items()
    ]
    bars.sort(key=lambda b: b.propensity, reverse=True)
    return bars


def render_propensity_bars_as_text(
    bars: List[PropensityBar],
    bar_width: int = 30,
    highlight_reaction_id: Optional[str] = None,
) -> str:
    """
    Renderiza los datos ya preparados como texto plano (sin `rich`).
    Sirve tanto de alternativa cuando `rich` no está instalado como
    de formato fácil de testear sin depender de una librería externa.
    """
    if not bars:
        return "(sin reacciones)"

    max_id_len = max(len(b.reaction_id) for b in bars)
    lines = []
    for bar in bars:
        filled = int(round(bar.fraction * bar_width))
        bar_str = "#" * filled + "-" * (bar_width - filled)
        marker = " <-- disparó" if bar.reaction_id == highlight_reaction_id else ""
        lines.append(
            f"{bar.reaction_id.ljust(max_id_len)} | {bar_str} | "
            f"{bar.propensity:8.3f} ({bar.fraction * 100:5.1f}%){marker}"
        )
    return "\n".join(lines)


def render_propensity_bars_rich(
    bars: List[PropensityBar],
    highlight_reaction_id: Optional[str] = None,
    title: str = "Espacio de propensidades",
) -> None:
    """
    Imprime en terminal una tabla enriquecida (colores, barra gráfica)
    usando `rich`. Requiere `pip install rich` (whitepaper, sección
    7.1, "necesarias a partir de la Fase 5"). Lanza ImportError si no
    está instalado — `show_propensities()` captura ese caso y cae a
    texto plano.
    """
    from rich.console import Console
    from rich.table import Table
    from rich.bar import Bar

    console = Console()
    table = Table(title=title)
    table.add_column("Reacción")
    table.add_column("Propensidad", justify="right")
    table.add_column("Peso relativo")

    for bar in bars:
        style = "bold green" if bar.reaction_id == highlight_reaction_id else None
        label = bar.reaction_id + (" ◀" if bar.reaction_id == highlight_reaction_id else "")
        table.add_row(
            label,
            f"{bar.propensity:.3f}",
            Bar(size=1.0, begin=0, end=bar.fraction),
            style=style,
        )

    console.print(table)


def show_propensities(
    propensities: Dict[str, float],
    highlight_reaction_id: Optional[str] = None,
    use_rich: bool = True,
) -> str:
    """
    Punto de entrada de conveniencia: prepara los datos y los muestra,
    usando `rich` si está disponible (y use_rich=True), o cayendo a
    texto plano si no. Devuelve siempre el texto plano equivalente
    (se haya mostrado también con rich o no), útil para logs y tests.
    """
    bars = prepare_propensity_bars(propensities)
    text = render_propensity_bars_as_text(bars, highlight_reaction_id=highlight_reaction_id)

    if use_rich:
        try:
            render_propensity_bars_rich(bars, highlight_reaction_id)
            return text
        except ImportError:
            pass

    print(text)
    return text
