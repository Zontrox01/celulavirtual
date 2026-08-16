"""
gui/breakpoint_builder.py

Construye condiciones de breakpoint (Callable[[dict], bool], el tipo
que espera Debugger.set_breakpoint()) a partir de valores simples
elegidos en la interfaz — nunca evaluando texto arbitrario (`eval`),
aunque el código se ejecute localmente en la máquina del propio
usuario: es una precaución barata y evita cualquier ambigüedad futura
si este mecanismo se reutiliza en otro contexto.

Sin dependencia de Qt: los widgets de gui/debugger_panel.py solo leen
sus valores (especie, operador, umbral) y se los pasan a
build_breakpoint_condition().
"""

from __future__ import annotations
from typing import Callable, Dict

ComparisonOperator = str  # uno de OPERATORS.keys()

OPERATORS: Dict[ComparisonOperator, Callable[[float, float], bool]] = {
    "<": lambda a, b: a < b,
    "<=": lambda a, b: a <= b,
    ">": lambda a, b: a > b,
    ">=": lambda a, b: a >= b,
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
}


class BreakpointBuilderError(ValueError):
    """Se lanza si los parámetros para construir un breakpoint son inválidos."""


def build_breakpoint_condition(
    species_id: str, operator: ComparisonOperator, threshold: float
) -> Callable[[Dict[str, int]], bool]:
    """
    Construye una condición `state -> bool` que compara
    `state[species_id]` contra `threshold` usando `operator`.

    Especies ausentes del estado se tratan como 0 (igual que hace
    engine/forensics.py), en vez de lanzar KeyError en mitad de una
    simulación — un breakpoint sobre una especie que todavía no existe
    simplemente no se dispara hasta que aparezca.
    """
    if operator not in OPERATORS:
        raise BreakpointBuilderError(
            f"Operador '{operator}' no soportado. Usa uno de: {sorted(OPERATORS.keys())}."
        )
    if not species_id:
        raise BreakpointBuilderError("species_id no puede estar vacío.")

    compare = OPERATORS[operator]

    def condition(state: Dict[str, int]) -> bool:
        return compare(state.get(species_id, 0), threshold)

    return condition


def describe_breakpoint(species_id: str, operator: ComparisonOperator, threshold: float) -> str:
    """Descripción legible de un breakpoint, para mostrar en la lista de la GUI."""
    return f"{species_id} {operator} {threshold:g}"
