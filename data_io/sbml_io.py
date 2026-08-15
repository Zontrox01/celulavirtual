"""
data_io/sbml_io.py

Import/export de modelos en formato SBML (whitepaper, sección 7).

Escrito directamente sobre `xml.etree.ElementTree` (librería estándar)
en vez de sobre `libsbml`: genera SBML Level 3 Version 1 válido, sin
depender de una librería adicional. `libsbml` sigue siendo la opción
recomendada en el whitepaper si hace falta VALIDACIÓN estricta contra
el esquema oficial o importar modelos SBML arbitrarios de terceros —
este módulo cubre el caso de uso principal (exportar nuestros propios
modelos para inspeccionarlos en herramientas externas, y reimportarlos
tal cual) sin esa dependencia extra.

Simplificaciones deliberadas:
  - La cinética se exporta como acción de masas CONTINUA estándar de
    SBML (tasa = k * Π [especie_i]^coeficiente_i), la convención
    habitual en biología de sistemas. Para reacciones con coeficiente
    estequiométrico > 1 sobre la misma especie (p. ej. 2A -> B), esto
    difiere ligeramente de la propensidad discreta combinatoria exacta
    que usa nuestro motor SSA (engine/reactions.py) — es una
    aproximación continua estándar, no un error, pero vale la pena
    saber que no es una copia exacta bit a bit de la cinética interna.
  - Las especies catalíticas (que aparecen igual en reactivos y
    productos, como RNAP o un ribosoma) se exportan apareciendo en
    ambas listas, en vez de usar el concepto SBML de "modifier" — más
    simple y con round-trip exacto para nuestro propio formato, a
    costa de no ser el modismo SBML "canónico" para catalizadores.
  - El import está pensado para reabrir SBML generado por ESTE mismo
    módulo (round-trip fiel), no para parsear cualquier SBML arbitrario
    de terceros — eso exigiría un parser de MathML general, fuera de
    alcance del MVP.
"""

from __future__ import annotations
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Union

from engine.reactions import Reaction

SBML_NS = "http://www.sbml.org/sbml/level3/version1/core"
MATHML_NS = "http://www.w3.org/1998/Math/MathML"
DEFAULT_COMPARTMENT_ID = "citoplasma"

ET.register_namespace("", SBML_NS)
ET.register_namespace("math", MATHML_NS)


class SBMLError(ValueError):
    """Se lanza ante cualquier problema al generar o interpretar SBML."""


# ---------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------

def _build_kinetic_law_math(reaction: Reaction) -> ET.Element:
    """Construye el MathML de acción de masas: k * especie1^coef1 * especie2^coef2 * ..."""
    math_el = ET.Element(f"{{{MATHML_NS}}}math")

    factors: List[ET.Element] = []
    k_ci = ET.Element(f"{{{MATHML_NS}}}ci")
    k_ci.text = "k"
    factors.append(k_ci)

    for species_id, coeff in reaction.reactants.items():
        species_ci = ET.Element(f"{{{MATHML_NS}}}ci")
        species_ci.text = species_id
        if coeff == 1:
            factors.append(species_ci)
        else:
            power = ET.Element(f"{{{MATHML_NS}}}apply")
            ET.SubElement(power, f"{{{MATHML_NS}}}power")
            power.append(species_ci)
            exponent = ET.SubElement(power, f"{{{MATHML_NS}}}cn")
            exponent.text = str(coeff)
            factors.append(power)

    if len(factors) == 1:
        math_el.append(factors[0])
        return math_el

    apply_el = ET.SubElement(math_el, f"{{{MATHML_NS}}}apply")
    ET.SubElement(apply_el, f"{{{MATHML_NS}}}times")
    for factor in factors:
        apply_el.append(factor)

    return math_el


def _build_reaction_element(reaction: Reaction) -> ET.Element:
    reaction_el = ET.Element(
        "reaction",
        {"id": reaction.reaction_id, "reversible": "false", "fast": "false"},
    )

    if reaction.reactants:
        reactants_el = ET.SubElement(reaction_el, "listOfReactants")
        for species_id, coeff in reaction.reactants.items():
            ET.SubElement(
                reactants_el,
                "speciesReference",
                {"species": species_id, "stoichiometry": str(coeff), "constant": "true"},
            )

    if reaction.products:
        products_el = ET.SubElement(reaction_el, "listOfProducts")
        for species_id, coeff in reaction.products.items():
            ET.SubElement(
                products_el,
                "speciesReference",
                {"species": species_id, "stoichiometry": str(coeff), "constant": "true"},
            )

    kinetic_law_el = ET.SubElement(reaction_el, "kineticLaw")
    kinetic_law_el.append(_build_kinetic_law_math(reaction))
    params_el = ET.SubElement(kinetic_law_el, "listOfLocalParameters")
    ET.SubElement(params_el, "localParameter", {"id": "k", "value": repr(reaction.rate_constant)})

    return reaction_el


def export_cell_to_sbml_string(
    cell,
    model_id: str = "celula_virtual",
    model_name: str = "Célula Virtual",
) -> str:
    """
    Genera el SBML (como string) que representa el estado actual de
    `cell`: todas sus especies (con su cantidad actual como amount
    inicial) y todas sus reacciones (con su cinética de acción de
    masas). No incluye histórico de eventos ni trayectorias — solo la
    red de reacciones y el estado molecular, que es lo que SBML sabe
    representar.
    """
    sbml_el = ET.Element(
        "sbml",
        {"xmlns": SBML_NS, "level": "3", "version": "1"},
    )
    model_el = ET.SubElement(sbml_el, "model", {"id": model_id, "name": model_name})

    compartments_el = ET.SubElement(model_el, "listOfCompartments")
    ET.SubElement(
        compartments_el,
        "compartment",
        {"id": DEFAULT_COMPARTMENT_ID, "constant": "true", "size": "1"},
    )

    species_el = ET.SubElement(model_el, "listOfSpecies")
    for species_id in cell.species.species_ids():
        count = cell.species.get_count(species_id)
        ET.SubElement(
            species_el,
            "species",
            {
                "id": species_id,
                "compartment": DEFAULT_COMPARTMENT_ID,
                "initialAmount": str(count),
                "hasOnlySubstanceUnits": "true",
                "boundaryCondition": "false",
                "constant": "false",
            },
        )

    reactions_el = ET.SubElement(model_el, "listOfReactions")
    for reaction in cell.reactions.reactions():
        reactions_el.append(_build_reaction_element(reaction))

    return ET.tostring(sbml_el, encoding="unicode", xml_declaration=False)


def export_cell_to_sbml_file(
    cell,
    path: Union[str, Path],
    model_id: str = "celula_virtual",
    model_name: str = "Célula Virtual",
) -> Path:
    """Igual que export_cell_to_sbml_string(), pero escribe el resultado en `path`."""
    path = Path(path)
    xml_string = export_cell_to_sbml_string(cell, model_id=model_id, model_name=model_name)
    path.write_text('<?xml version="1.0" encoding="UTF-8"?>\n' + xml_string, encoding="utf-8")
    return path


# ---------------------------------------------------------------------
# Import (round-trip fiel para SBML exportado por este mismo módulo)
# ---------------------------------------------------------------------

def _parse_reaction_element(reaction_el: ET.Element) -> Reaction:
    ns = {"sbml": SBML_NS}
    reaction_id = reaction_el.attrib["id"]

    reactants: Dict[str, int] = {}
    reactants_el = reaction_el.find("sbml:listOfReactants", ns)
    if reactants_el is not None:
        for ref in reactants_el.findall("sbml:speciesReference", ns):
            reactants[ref.attrib["species"]] = int(float(ref.attrib["stoichiometry"]))

    products: Dict[str, int] = {}
    products_el = reaction_el.find("sbml:listOfProducts", ns)
    if products_el is not None:
        for ref in products_el.findall("sbml:speciesReference", ns):
            products[ref.attrib["species"]] = int(float(ref.attrib["stoichiometry"]))

    kinetic_law_el = reaction_el.find("sbml:kineticLaw", ns)
    if kinetic_law_el is None:
        raise SBMLError(f"La reacción '{reaction_id}' no tiene kineticLaw: no se puede recuperar la tasa.")

    params_el = kinetic_law_el.find("sbml:listOfLocalParameters", ns)
    if params_el is None:
        raise SBMLError(f"La reacción '{reaction_id}' no tiene listOfLocalParameters con 'k'.")

    rate_constant = None
    for param in params_el.findall("sbml:localParameter", ns):
        if param.attrib.get("id") == "k":
            rate_constant = float(param.attrib["value"])
            break
    if rate_constant is None:
        raise SBMLError(f"La reacción '{reaction_id}' no define el parámetro local 'k'.")

    return Reaction(reaction_id, reactants=reactants, products=products, rate_constant=rate_constant)


def import_sbml_string(xml_string: str) -> Tuple[Dict[str, int], List[Reaction]]:
    """
    Recupera (especies_iniciales, reacciones) a partir de un SBML
    generado por export_cell_to_sbml_string()/export_cell_to_sbml_file()
    de este mismo módulo. No es un parser de SBML arbitrario de
    terceros: asume el convenio de exportación de arriba (cinética de
    acción de masas con un único parámetro local 'k').
    """
    try:
        root = ET.fromstring(xml_string)
    except ET.ParseError as e:
        raise SBMLError(f"XML inválido: {e}") from None

    ns = {"sbml": SBML_NS}
    model_el = root.find("sbml:model", ns)
    if model_el is None:
        raise SBMLError("El documento SBML no contiene un elemento <model>.")

    species: Dict[str, int] = {}
    species_list_el = model_el.find("sbml:listOfSpecies", ns)
    if species_list_el is not None:
        for species_el in species_list_el.findall("sbml:species", ns):
            species_id = species_el.attrib["id"]
            amount = int(float(species_el.attrib.get("initialAmount", "0")))
            species[species_id] = amount

    reactions: List[Reaction] = []
    reactions_list_el = model_el.find("sbml:listOfReactions", ns)
    if reactions_list_el is not None:
        for reaction_el in reactions_list_el.findall("sbml:reaction", ns):
            reactions.append(_parse_reaction_element(reaction_el))

    return species, reactions


def import_sbml_file(path: Union[str, Path]) -> Tuple[Dict[str, int], List[Reaction]]:
    """Igual que import_sbml_string(), leyendo el XML desde `path`."""
    path = Path(path)
    if not path.exists():
        raise SBMLError(f"No se encuentra el archivo SBML: {path}")
    return import_sbml_string(path.read_text(encoding="utf-8"))


def build_cell_from_sbml(path_or_string: Union[str, Path], seed=None):
    """
    Reconstruye una Cell nueva a partir de un SBML exportado por este
    módulo: crea las especies con sus cantidades y registra todas las
    reacciones. Import tardío de Cell para evitar un ciclo de imports
    (cell.py no depende de data_io, pero data_io sí puede depender de
    cell.py solo en esta función de conveniencia).
    """
    from cell import Cell

    if isinstance(path_or_string, Path) or (
        isinstance(path_or_string, str) and Path(path_or_string).exists() and "<sbml" not in path_or_string
    ):
        species, reactions = import_sbml_file(path_or_string)
    else:
        species, reactions = import_sbml_string(path_or_string)

    cell = Cell(seed=seed)
    for species_id, count in species.items():
        cell.add_species(species_id, count)
    for reaction in reactions:
        cell.add_reaction(reaction)

    return cell
