"""
biology/genome.py

GenomeModule: primer módulo biológico (Capa 2), construido sobre un
GenomeData ya cargado (data_io/genome_loader.py). Por cada gen anotado,
registra en una Cell: una especie que representa la copia del gen (no
se consume al transcribirse, actúa como catalizador), una especie de
ARNm, y una reacción de transcripción.

Simplificación deliberada de la Fase 1: la transcripción se modela
como un único paso estocástico (RNAP + gen -> RNAP + gen + ARNm), en
vez de como una cadena explícita de unión/elongación/terminación con
estados intermedios. Es el modelo de "dos etapas" estándar en
expresión génica estocástica (transcripción de un paso + traducción de
un paso), con amplio precedente en la literatura (p.ej. Thattai & van
Oudenaarden, 2001). La elongación explícita queda como refinamiento de
biology/transcription.py sin romper esta interfaz: seguiría
generándose "una reacción de transcripción por gen", solo que pasaría
a ser el paso final de una cadena más larga en vez de la única.

No conoce nada de traducción, metabolismo ni depuración — solo genera
las especies y reacciones de transcripción a partir de la anotación.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional

from cell import Cell
from data_io.genome_loader import GenomeData
from engine.reactions import Reaction

DEFAULT_RNAP_SPECIES_ID = "RNAP"


class GenomeModuleError(ValueError):
    """Se lanza si el genoma no puede instalarse sobre la célula (conflictos de nombres, etc.)."""


@dataclass(frozen=True)
class InstalledGene:
    """Referencia a las especies/reacción generadas para un gen concreto al instalarlo."""

    gene_id: str
    gene_species_id: str
    mrna_species_id: str
    transcription_reaction_id: str


class GenomeModule:
    """
    Instala un genoma cargado (GenomeData) sobre una Cell: una especie
    de ARNm y una reacción de transcripción por cada gen anotado.

    La ARN-polimerasa se trata como maquinaria preexistente de la
    célula, no codificada por gen en esta fase (whitepaper, sección 5,
    nota "opcional" sobre simplificar la Fase 1) — su cantidad inicial
    se indica al instalar, no se deriva del genoma.
    """

    def __init__(self, genome: GenomeData) -> None:
        self.genome = genome
        self._installed: Dict[str, InstalledGene] = {}

    def install(
        self,
        cell: Cell,
        rnap_species_id: str = DEFAULT_RNAP_SPECIES_ID,
        rnap_initial_count: int = 30,
        base_transcription_rate: float = 1.0,
        promoter_strengths: Optional[Dict[str, float]] = None,
        atp_species_id: Optional[str] = None,
        atp_cost_per_transcription: int = 0,
    ) -> List[InstalledGene]:
        """
        Registra en `cell` las especies y reacciones de transcripción
        de todos los genes del genoma.

        `promoter_strengths`: multiplicador de tasa por gen, indexado
        por `promoter_id` si el gen lo tiene, o por `gene_id` si no.
        Los genes sin entrada en este diccionario usan multiplicador 1.0
        (es decir, tasa = base_transcription_rate).

        `atp_species_id` / `atp_cost_per_transcription`: acoplamiento
        energético opcional (whitepaper, Fase 3). Si se indican ambos
        (con coste > 0), cada transcripción consume esa cantidad de
        ATP además de RNAP y el gen — la especie de ATP debe existir
        ya en la célula (instalada por biology/metabolism.py). Por
        defecto no hay coste de ATP, para no romper el comportamiento
        ya validado en fases anteriores.

        Devuelve la lista de genes instalados, con los IDs de especie
        y reacción generados para cada uno — útiles para poner
        breakpoints o inspeccionar resultados después.
        """
        promoter_strengths = promoter_strengths or {}
        couple_to_atp = atp_species_id is not None and atp_cost_per_transcription > 0

        if couple_to_atp and not cell.species.has_species(atp_species_id):
            raise GenomeModuleError(
                f"Se pidió acoplar la transcripción a ATP ('{atp_species_id}'), pero esa "
                f"especie no existe en la célula (¿se instaló primero MetabolismModule?)."
            )

        if not cell.species.has_species(rnap_species_id):
            cell.add_species(rnap_species_id, rnap_initial_count)

        installed: List[InstalledGene] = []
        for gene in self.genome.genes:
            gene_species_id = f"gene_{gene.gene_id}"
            mrna_species_id = f"mRNA_{gene.gene_id}"
            reaction_id = f"transcribe_{gene.gene_id}"

            if cell.species.has_species(gene_species_id) or cell.species.has_species(mrna_species_id):
                raise GenomeModuleError(
                    f"Conflicto al instalar el gen '{gene.gene_id}': "
                    f"'{gene_species_id}' o '{mrna_species_id}' ya existen en la célula "
                    f"(¿se ha instalado este genoma dos veces?)."
                )

            cell.add_species(gene_species_id, 1)  # una copia del gen, no se consume
            cell.add_species(mrna_species_id, 0)

            lookup_key = gene.promoter_id if gene.promoter_id else gene.gene_id
            strength_multiplier = promoter_strengths.get(lookup_key, 1.0)
            rate = base_transcription_rate * strength_multiplier

            reactants = {rnap_species_id: 1, gene_species_id: 1}
            products = {rnap_species_id: 1, gene_species_id: 1, mrna_species_id: 1}
            if couple_to_atp:
                reactants[atp_species_id] = atp_cost_per_transcription  # se consume, no se repone

            cell.add_reaction(Reaction(reaction_id, reactants=reactants, products=products, rate_constant=rate))

            record = InstalledGene(
                gene_id=gene.gene_id,
                gene_species_id=gene_species_id,
                mrna_species_id=mrna_species_id,
                transcription_reaction_id=reaction_id,
            )
            installed.append(record)
            self._installed[gene.gene_id] = record

        return installed

    def get_installed(self, gene_id: str) -> InstalledGene:
        try:
            return self._installed[gene_id]
        except KeyError:
            raise GenomeModuleError(
                f"El gen '{gene_id}' no se ha instalado todavía (¿llamaste a install()?)."
            ) from None
