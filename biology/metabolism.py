"""
biology/metabolism.py

MetabolismModule: quinto módulo biológico (Capa 2). Provee un
metabolismo energético mínimo: importación de glucosa desde el
entorno (reacción de orden cero) y una vía glicolítica simplificada
que la convierte en ATP, opcionalmente catalizada por una enzima ya
traducida (p. ej. el producto de un gen metabólico instalado por
GenomeModule + TranslationModule).

Es el punto de acoplamiento energético del modelo (whitepaper, sección
5): las reacciones de transcripción/traducción pueden configurarse
(vía GenomeModule.install(atp_species_id=...) y
TranslationModule.install(atp_species_id=...)) para consumir el ATP
que este módulo produce, de forma que la expresión génica se
ralentiza si la energía escasea — cerrando el ciclo entre metabolismo
y expresión génica que motivó este proyecto.

Simplificación deliberada: la glicólisis real tiene ~10 pasos
enzimáticos; aquí se colapsa en un único paso estocástico
(Glucosa [+ Enzima] -> [Enzima +] N x ATP), mismo criterio de
simplificación que en GenomeModule/TranslationModule. El rendimiento
neto real de la glicólisis es de 2 ATP por glucosa; se usa por
defecto pero es configurable.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

from cell import Cell
from engine.reactions import Reaction

DEFAULT_GLUCOSE_SPECIES_ID = "Glucose"
DEFAULT_ATP_SPECIES_ID = "ATP"
DEFAULT_ATP_YIELD_PER_GLUCOSE = 2
DEFAULT_GLUCOSE_IMPORT_RATE = 1.0
DEFAULT_GLYCOLYSIS_RATE = 0.1


class MetabolismModuleError(ValueError):
    """Se lanza si el metabolismo no puede instalarse (enzima inexistente, conflictos, parámetros inválidos)."""


@dataclass(frozen=True)
class MetabolismInstallation:
    """Referencia a las especies/reacciones instaladas por el metabolismo."""

    glucose_species_id: str
    atp_species_id: str
    import_reaction_id: str
    glycolysis_reaction_id: str
    enzyme_species_id: Optional[str]


class MetabolismModule:
    """
    Instala un metabolismo energético mínimo: importación continua de
    glucosa desde el entorno (reacción de orden cero, sin reactivos —
    representa el suministro externo, no modelado como especie) más
    una reacción de glicólisis simplificada que la convierte en ATP.
    """

    def __init__(self) -> None:
        self._installation: Optional[MetabolismInstallation] = None

    def install(
        self,
        cell: Cell,
        enzyme_species_id: Optional[str] = None,
        glucose_species_id: str = DEFAULT_GLUCOSE_SPECIES_ID,
        atp_species_id: str = DEFAULT_ATP_SPECIES_ID,
        glucose_initial_count: int = 0,
        atp_initial_count: int = 100,
        glucose_import_rate: float = DEFAULT_GLUCOSE_IMPORT_RATE,
        glycolysis_rate: float = DEFAULT_GLYCOLYSIS_RATE,
        atp_yield_per_glucose: int = DEFAULT_ATP_YIELD_PER_GLUCOSE,
    ) -> MetabolismInstallation:
        """
        Instala en `cell` la importación de glucosa y la glicólisis.

        `enzyme_species_id`: especie de proteína que cataliza la
        glicólisis (p. ej. "Prot_metA", el producto de un gen
        metabólico ya traducido). Si se indica, debe existir ya en la
        célula, y la glicólisis la requiere como catalizador (no se
        consume: aparece también como producto). Si se omite, la
        glicólisis ocurre sin depender de ninguna enzima —útil para
        probar este módulo de forma aislada, o como "metabolismo de
        fondo" mientras no haya un gen metabólico instalado.
        """
        if self._installation is not None:
            raise MetabolismModuleError("El metabolismo ya está instalado en esta célula.")

        if enzyme_species_id is not None and not cell.species.has_species(enzyme_species_id):
            raise MetabolismModuleError(
                f"La enzima '{enzyme_species_id}' no existe en la célula "
                f"(¿se instaló primero el gen que la produce, vía GenomeModule + "
                f"TranslationModule?)."
            )
        if atp_yield_per_glucose <= 0:
            raise MetabolismModuleError(
                f"atp_yield_per_glucose debe ser positivo (recibido {atp_yield_per_glucose})."
            )

        if not cell.species.has_species(glucose_species_id):
            cell.add_species(glucose_species_id, glucose_initial_count)
        if not cell.species.has_species(atp_species_id):
            cell.add_species(atp_species_id, atp_initial_count)

        import_reaction_id = "import_glucose"
        cell.add_reaction(
            Reaction(
                import_reaction_id,
                reactants={},  # orden cero: entorno externo, no modelado como especie
                products={glucose_species_id: 1},
                rate_constant=glucose_import_rate,
            )
        )

        glycolysis_reaction_id = "glycolysis"
        glycolysis_reactants = {glucose_species_id: 1}
        glycolysis_products = {atp_species_id: atp_yield_per_glucose}
        if enzyme_species_id is not None:
            glycolysis_reactants[enzyme_species_id] = 1
            glycolysis_products[enzyme_species_id] = 1  # catalizador, no se consume

        cell.add_reaction(
            Reaction(
                glycolysis_reaction_id,
                reactants=glycolysis_reactants,
                products=glycolysis_products,
                rate_constant=glycolysis_rate,
            )
        )

        self._installation = MetabolismInstallation(
            glucose_species_id=glucose_species_id,
            atp_species_id=atp_species_id,
            import_reaction_id=import_reaction_id,
            glycolysis_reaction_id=glycolysis_reaction_id,
            enzyme_species_id=enzyme_species_id,
        )
        return self._installation

    def get_installation(self) -> MetabolismInstallation:
        if self._installation is None:
            raise MetabolismModuleError("El metabolismo no se ha instalado todavía.")
        return self._installation
