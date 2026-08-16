"""
gui/app_state.py

Estado de la aplicación GUI, completamente independiente de PySide6.

Responsabilidad única: mantener la Cell activa y los módulos
biológicos instalados sobre ella, y exponer operaciones de alto nivel
(cargar genoma, ejecutar, depurar, dividir, exportar/importar SBML)
como métodos normales de Python — sin ninguna dependencia de Qt.

Esto permite testear toda la lógica de la GUI con `assert` puro, igual
que el resto del proyecto, y mantiene los archivos de Qt (theme.py,
*_panel.py, main_window.py) como envoltorios finos que solo conectan
widgets a estos métodos, sin lógica de negocio propia.

No modifica ningún módulo existente: solo los COMPONE, en el mismo
orden y con las mismas restricciones que ya exigían sus propios tests
(p. ej. la especie de ATP debe existir antes de acoplar transcripción
a energía). Ver `GenomeSetupOptions.metabolism_mode` para el único
punto donde estas restricciones obligan a elegir entre dos
combinaciones de funcionalidades que no pueden activarse ambas a la
vez con el código actual (documentado también en el whitepaper).
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from cell import Cell
from data_io.genome_loader import load_genome, GenomeData
from data_io.sbml_io import export_cell_to_sbml_file, build_cell_from_sbml
from biology.genome import GenomeModule, InstalledGene
from biology.translation import TranslationModule, TranslatedGene
from biology.degradation import DegradationModule
from biology.metabolism import MetabolismModule
from biology.regulation import RegulationModule, RepressionLink
from biology.replication import ReplicationModule, total_count_condition
from engine.debugger import Debugger
from engine.forensics import ForensicsAnalyzer, ForensicReport

METABOLISM_MODES = ("none", "background_with_atp_coupling", "enzyme_catalyzed")


class AppStateError(ValueError):
    """Se lanza ante cualquier operación inválida sobre el estado de la aplicación."""


@dataclass
class GenomeSetupOptions:
    """
    Parámetros de configuración al cargar e instalar un genoma.

    `metabolism_mode` controla cómo se instala el metabolismo,
    porque el código actual no permite combinar libremente enzima
    catalizadora y acoplamiento a ATP en transcripción/traducción (la
    enzima exige que la traducción ya haya ocurrido; el acoplamiento a
    ATP exige que el ATP exista ANTES de instalar la transcripción):

      - "none": sin metabolismo, comportamiento de siempre.
      - "background_with_atp_coupling": metabolismo "de fondo" (sin
        enzima), con la posibilidad de acoplar el gasto de ATP a
        transcripción/traducción.
      - "enzyme_catalyzed": metabolismo catalizado por la proteína de
        un gen concreto ya traducido, SIN acoplamiento a ATP en
        transcripción/traducción (mismo orden que el test de
        integración de biology/metabolism.py).
    """

    fasta_path: Union[str, Path]
    annotation_path: Union[str, Path]
    seed: Optional[int] = None
    rnap_initial_count: int = 30
    ribosome_initial_count: int = 20
    mrna_half_life: float = 300.0
    protein_half_life: float = 3600.0
    promoter_strengths: Optional[Dict[str, float]] = None
    rbs_strengths: Optional[Dict[str, float]] = None

    metabolism_mode: str = "none"
    enzyme_gene_id: Optional[str] = None
    glucose_initial_count: int = 0
    atp_initial_count: int = 100
    glucose_import_rate: float = 1.0
    glycolysis_rate: float = 0.1
    atp_yield_per_glucose: int = 2
    atp_cost_per_transcription: int = 2
    atp_cost_per_translation: int = 4


class AppState:
    """
    Estado central de la aplicación: la Cell activa, sus módulos
    biológicos, y el número de generación tras divisiones celulares.
    """

    def __init__(self) -> None:
        self.cell: Optional[Cell] = None
        self.genome: Optional[GenomeData] = None
        self.genome_module: Optional[GenomeModule] = None
        self.translation_module: Optional[TranslationModule] = None
        self.degradation_module: Optional[DegradationModule] = None
        self.metabolism_module: Optional[MetabolismModule] = None
        self.regulation_module: RegulationModule = RegulationModule()
        self.replication_module: Optional[ReplicationModule] = None
        self.installed_genes: List[InstalledGene] = []
        self.translated_genes: List[TranslatedGene] = []
        self.generation: int = 0

    # ------------------------------------------------------------------
    # Carga e instalación del genoma
    # ------------------------------------------------------------------

    def load_and_install_genome(self, options: GenomeSetupOptions) -> None:
        """Carga un genoma desde archivo e instala transcripción, traducción, degradación y (opcional) metabolismo."""
        if options.metabolism_mode not in METABOLISM_MODES:
            raise AppStateError(
                f"metabolism_mode debe ser uno de {METABOLISM_MODES} "
                f"(recibido '{options.metabolism_mode}')."
            )

        genome = load_genome(options.fasta_path, options.annotation_path)
        cell = Cell(seed=options.seed)
        genome_module = GenomeModule(genome)
        metabolism_module: Optional[MetabolismModule] = None

        genome_atp_kwargs: Dict[str, object] = {}
        if options.metabolism_mode == "background_with_atp_coupling":
            metabolism_module = MetabolismModule()
            metabolism_module.install(
                cell,
                glucose_initial_count=options.glucose_initial_count,
                atp_initial_count=options.atp_initial_count,
                glucose_import_rate=options.glucose_import_rate,
                glycolysis_rate=options.glycolysis_rate,
                atp_yield_per_glucose=options.atp_yield_per_glucose,
            )
            genome_atp_kwargs = {
                "atp_species_id": "ATP",
                "atp_cost_per_transcription": options.atp_cost_per_transcription,
            }

        installed_genes = genome_module.install(
            cell,
            rnap_initial_count=options.rnap_initial_count,
            promoter_strengths=options.promoter_strengths,
            **genome_atp_kwargs,
        )

        translation_module = TranslationModule(installed_genes)
        translation_atp_kwargs: Dict[str, object] = {}
        if options.metabolism_mode == "background_with_atp_coupling":
            translation_atp_kwargs = {
                "atp_species_id": "ATP",
                "atp_cost_per_translation": options.atp_cost_per_translation,
            }
        translated_genes = translation_module.install(
            cell,
            ribosome_initial_count=options.ribosome_initial_count,
            rbs_strengths=options.rbs_strengths,
            **translation_atp_kwargs,
        )

        if options.metabolism_mode == "enzyme_catalyzed":
            if not options.enzyme_gene_id:
                raise AppStateError("metabolism_mode='enzyme_catalyzed' requiere indicar enzyme_gene_id.")
            enzyme_species_id = f"Prot_{options.enzyme_gene_id}"
            if not cell.species.has_species(enzyme_species_id):
                raise AppStateError(
                    f"El gen '{options.enzyme_gene_id}' no tiene proteína traducida "
                    f"('{enzyme_species_id}'); comprueba que esté en el genoma cargado."
                )
            metabolism_module = MetabolismModule()
            metabolism_module.install(
                cell,
                enzyme_species_id=enzyme_species_id,
                glucose_initial_count=options.glucose_initial_count,
                atp_initial_count=options.atp_initial_count,
                glucose_import_rate=options.glucose_import_rate,
                glycolysis_rate=options.glycolysis_rate,
                atp_yield_per_glucose=options.atp_yield_per_glucose,
            )

        degradation_module = DegradationModule(installed_genes, translated_genes)
        degradation_module.install(
            cell,
            mrna_half_life=options.mrna_half_life,
            protein_half_life=options.protein_half_life,
        )

        self.cell = cell
        self.genome = genome
        self.genome_module = genome_module
        self.translation_module = translation_module
        self.degradation_module = degradation_module
        self.metabolism_module = metabolism_module
        self.installed_genes = installed_genes
        self.translated_genes = translated_genes
        self.regulation_module = RegulationModule()
        self.replication_module = None
        self.generation = 0

    def require_cell(self) -> Cell:
        if self.cell is None:
            raise AppStateError("No hay ninguna célula cargada todavía. Carga un genoma primero.")
        return self.cell

    def gene_ids(self) -> List[str]:
        return [g.gene_id for g in self.installed_genes]

    # ------------------------------------------------------------------
    # Regulación
    # ------------------------------------------------------------------

    def add_repression(
        self,
        target_gene_id: str,
        repressor_gene_id: str,
        binding_rate: float = 1.0,
        unbinding_rate: float = 1.0,
    ) -> RepressionLink:
        cell = self.require_cell()
        target = next((g for g in self.installed_genes if g.gene_id == target_gene_id), None)
        repressor = next((g for g in self.translated_genes if g.gene_id == repressor_gene_id), None)
        if target is None:
            raise AppStateError(f"Gen diana '{target_gene_id}' no encontrado entre los genes instalados.")
        if repressor is None:
            raise AppStateError(f"Gen represor '{repressor_gene_id}' no encontrado entre los genes traducidos.")
        return self.regulation_module.add_repression(cell, target, repressor, binding_rate, unbinding_rate)

    # ------------------------------------------------------------------
    # Ejecución
    # ------------------------------------------------------------------

    def run_steps(self, n: int):
        return self.require_cell().run(max_steps=n)

    def run_time(self, seconds: float):
        cell = self.require_cell()
        return cell.run(max_time=cell.time + seconds)

    def species_snapshot(self) -> Dict[str, int]:
        return self.require_cell().species.snapshot()

    def trajectory(self):
        return self.require_cell().get_trajectory()

    # ------------------------------------------------------------------
    # Depuración
    # ------------------------------------------------------------------

    def new_debugger(self, enable_snapshots: bool = True, snapshot_interval: int = 20) -> Debugger:
        return self.require_cell().debug(enable_snapshots=enable_snapshots, snapshot_interval=snapshot_interval)

    # ------------------------------------------------------------------
    # Análisis forense
    # ------------------------------------------------------------------

    def analyze_species(self, species_id: str, window_start: float, window_end: float) -> ForensicReport:
        cell = self.require_cell()
        analyzer = ForensicsAnalyzer(cell.get_events(), cell.reactions)
        return analyzer.analyze_species(species_id, window_start, window_end)

    # ------------------------------------------------------------------
    # División celular
    # ------------------------------------------------------------------

    def divide(
        self,
        threshold_species_ids: List[str],
        threshold: int,
        gene_species_ids: List[str],
    ) -> Tuple[Cell, Cell]:
        """
        Divide la célula activa si se cumple la condición (suma de
        `threshold_species_ids` >= `threshold`), replicando antes las
        copias de gen indicadas en `gene_species_ids` (whitepaper,
        aviso de `biology/replication.py`: sin esto, una hija puede
        quedarse sin ninguna copia de un gen presente en una sola
        copia). Devuelve las dos células hijas sin adoptar ninguna
        todavía — usar adopt_daughter() para continuar la simulación
        con una de ellas.
        """
        cell = self.require_cell()
        condition = total_count_condition(threshold_species_ids, threshold)
        replication_module = ReplicationModule(condition)
        if not replication_module.should_divide(cell):
            raise AppStateError(
                f"Todavía no se cumple la condición de división "
                f"(suma de {threshold_species_ids} < {threshold})."
            )
        replication_module.replicate_gene_copies(cell, gene_species_ids)
        daughter_a, daughter_b = replication_module.divide(cell)
        self.replication_module = replication_module
        return daughter_a, daughter_b

    def adopt_daughter(self, daughter: Cell) -> None:
        """Sustituye la célula activa por una de las hijas de una división ya ocurrida."""
        self.cell = daughter
        self.generation += 1

    # ------------------------------------------------------------------
    # SBML
    # ------------------------------------------------------------------

    def export_sbml(self, path: Union[str, Path]) -> Path:
        return export_cell_to_sbml_file(self.require_cell(), path)

    def import_sbml(self, path: Union[str, Path], seed: Optional[int] = None) -> None:
        """
        Sustituye la célula activa por una reconstruida desde SBML.
        Los módulos biológicos (genome_module, etc.) se pierden — el
        SBML solo sabe representar especies y reacciones, no de qué
        módulo biológico vinieron (whitepaper, sección 7).
        """
        cell = build_cell_from_sbml(path, seed=seed)
        self.cell = cell
        self.genome = None
        self.genome_module = None
        self.translation_module = None
        self.degradation_module = None
        self.metabolism_module = None
        self.regulation_module = RegulationModule()
        self.replication_module = None
        self.installed_genes = []
        self.translated_genes = []
        self.generation = 0
