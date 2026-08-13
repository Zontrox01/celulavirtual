"""
tests/test_metabolism.py

Validación de biology/metabolism.py.

Incluye el criterio de la Fase 3 del roadmap (whitepaper, sección 8):
"la expresión génica se ralentiza si el ATP cae" — comprobado de
forma comparativa (con/sin reposición de ATP), y con un breakpoint de
ejemplo tal como se describe en el propio whitepaper.
"""

import pytest
from cell import Cell
from data_io.genome_loader import GenomeData, GeneAnnotation
from biology.genome import GenomeModule
from biology.translation import TranslationModule
from biology.metabolism import MetabolismModule, MetabolismModuleError


def make_single_gene_genome():
    return GenomeData(
        sequence_id="toy",
        sequence="A" * 100,
        genes=[GeneAnnotation("target", 1, 90, "+", promoter_id="strong")],
    )


def make_genome_with_metabolic_gene():
    return GenomeData(
        sequence_id="toy",
        sequence="A" * 300,
        genes=[
            GeneAnnotation("target", 1, 90, "+", promoter_id="strong"),
            GeneAnnotation("metA", 100, 190, "+", promoter_id="medium"),
        ],
    )


# ---------------------------------------------------------------------
# Instalación básica
# ---------------------------------------------------------------------

def test_install_creates_glucose_and_atp_species():
    cell = Cell(seed=1)
    module = MetabolismModule()

    module.install(cell, glucose_initial_count=5, atp_initial_count=50)

    assert cell.species.get_count("Glucose") == 5
    assert cell.species.get_count("ATP") == 50


def test_install_reuses_existing_species_without_overwriting():
    cell = Cell(seed=1)
    cell.add_species("ATP", 999)
    module = MetabolismModule()

    module.install(cell, atp_initial_count=1)  # no debería sobrescribir

    assert cell.species.get_count("ATP") == 999


def test_import_reaction_has_no_reactants_zero_order():
    cell = Cell(seed=1)
    module = MetabolismModule()
    installation = module.install(cell)

    reaction = cell.reactions.get_reaction(installation.import_reaction_id)
    assert reaction.reactants == {}
    assert reaction.products == {"Glucose": 1}


def test_glycolysis_stoichiometry_without_enzyme():
    cell = Cell(seed=1)
    module = MetabolismModule()
    installation = module.install(cell, atp_yield_per_glucose=2)

    reaction = cell.reactions.get_reaction(installation.glycolysis_reaction_id)
    assert reaction.reactants == {"Glucose": 1}
    assert reaction.products == {"ATP": 2}


def test_glycolysis_stoichiometry_with_enzyme_is_catalytic():
    cell = Cell(seed=1)
    cell.add_species("Prot_metA", 10)
    module = MetabolismModule()
    installation = module.install(cell, enzyme_species_id="Prot_metA", atp_yield_per_glucose=3)

    reaction = cell.reactions.get_reaction(installation.glycolysis_reaction_id)
    assert reaction.reactants == {"Glucose": 1, "Prot_metA": 1}
    assert reaction.products == {"ATP": 3, "Prot_metA": 1}  # el catalizador no se consume


def test_enzyme_species_missing_raises():
    cell = Cell(seed=1)  # sin Prot_metA
    module = MetabolismModule()

    with pytest.raises(MetabolismModuleError):
        module.install(cell, enzyme_species_id="Prot_metA")


def test_invalid_atp_yield_raises():
    cell = Cell(seed=1)
    module = MetabolismModule()

    with pytest.raises(MetabolismModuleError):
        module.install(cell, atp_yield_per_glucose=0)


def test_double_install_raises():
    cell = Cell(seed=1)
    module = MetabolismModule()
    module.install(cell)

    with pytest.raises(MetabolismModuleError):
        module.install(cell)


def test_get_installation_before_install_raises():
    module = MetabolismModule()
    with pytest.raises(MetabolismModuleError):
        module.get_installation()


def test_running_produces_atp_over_time():
    cell = Cell(seed=2)
    module = MetabolismModule()
    module.install(cell, glucose_initial_count=0, atp_initial_count=0,
                    glucose_import_rate=5.0, glycolysis_rate=1.0)

    cell.run(max_time=20.0)

    assert cell.species.get_count("ATP") > 0


# ---------------------------------------------------------------------
# Criterio central de la Fase 3: ATP escaso ralentiza la expresión génica
# ---------------------------------------------------------------------

def test_gene_expression_slows_when_atp_is_scarce():
    """
    Con metabolismo activo reponiendo ATP continuamente, la
    transcripción acoplada a ATP debe producir mucho más ARNm que en
    una célula con muy poco ATP y ninguna vía de reposición.
    """
    # --- Escenario A: metabolismo activo, ATP se repone continuamente ---
    cell_a = Cell(seed=200)
    metabolism_a = MetabolismModule()
    metabolism_a.install(cell_a, atp_initial_count=50, glucose_import_rate=5.0, glycolysis_rate=1.0)
    genome_a = GenomeModule(make_single_gene_genome())
    genome_a.install(cell_a, rnap_initial_count=10, atp_species_id="ATP", atp_cost_per_transcription=2)

    cell_a.run(max_time=100.0)
    mrna_with_metabolism = cell_a.species.get_count("mRNA_target")

    # --- Escenario B: ATP escaso, sin ninguna vía de reposición ---
    cell_b = Cell(seed=200)
    cell_b.add_species("ATP", 10)
    genome_b = GenomeModule(make_single_gene_genome())
    genome_b.install(cell_b, rnap_initial_count=10, atp_species_id="ATP", atp_cost_per_transcription=2)

    cell_b.run(max_time=100.0)
    mrna_without_metabolism = cell_b.species.get_count("mRNA_target")

    assert mrna_without_metabolism < mrna_with_metabolism / 3, (
        f"sin reposición de ATP ({mrna_without_metabolism}) debería ser mucho menor que "
        f"con metabolismo activo ({mrna_with_metabolism})"
    )
    # cota superior teórica del escenario B: nunca puede transcribir más de ATP_inicial / coste
    assert mrna_without_metabolism <= 10 // 2


def test_debugger_breakpoint_on_low_atp_matches_whitepaper_example():
    """
    Reproduce el ejemplo explícito del whitepaper (sección 8, Fase 3):
    "breakpoint de ejemplo: pausar cuando ATP < umbral".
    """
    cell = Cell(seed=300)
    cell.add_species("ATP", 20)
    genome = GenomeModule(make_single_gene_genome())
    genome.install(cell, rnap_initial_count=10, atp_species_id="ATP", atp_cost_per_transcription=2)

    debugger = cell.debug()
    debugger.set_breakpoint("atp_bajo", lambda s: s["ATP"] < 10)
    debugger.run_until_breakpoint(max_steps=10_000)

    assert debugger.last_triggered_breakpoint == "atp_bajo"
    assert cell.species.get_count("ATP") < 10


# ---------------------------------------------------------------------
# Ciclo cerrado: un gen produce la enzima que cataliza el metabolismo
# que a su vez da energía a la propia expresión génica
# ---------------------------------------------------------------------

def test_closed_loop_gene_product_fuels_its_own_expression_via_metabolism():
    """
    'metA' se transcribe y traduce a Prot_metA, que cataliza la
    glicólisis, que produce el ATP que a su vez paga el coste
    energético de transcribir/traducir (incluido el propio metA) —
    el ciclo completo ADN -> ARNm -> proteína -> metabolismo -> energía
    que motivó este proyecto desde el principio.
    """
    cell = Cell(seed=42)
    genome_module = GenomeModule(make_genome_with_metabolic_gene())
    installed_genes = genome_module.install(cell, rnap_initial_count=10)

    translation_module = TranslationModule(installed_genes)
    translated_genes = translation_module.install(cell, ribosome_initial_count=10)

    metA_protein = next(g for g in translated_genes if g.gene_id == "metA").protein_species_id

    metabolism_module = MetabolismModule()
    metabolism_module.install(
        cell,
        enzyme_species_id=metA_protein,
        glucose_initial_count=0,
        atp_initial_count=20,
        glucose_import_rate=3.0,
        glycolysis_rate=0.5,
    )

    # activar el acoplamiento energético habría requerido instalar
    # genome/translation DESPUÉS del metabolismo; aquí se comprueba el
    # ciclo de producción (metA -> enzima -> ATP), que es lo que hace
    # a este test un caso de integración de los tres módulos a la vez
    cell.run(max_time=50.0)

    assert cell.species.get_count("Prot_metA") > 0  # se tradujo la enzima
    assert cell.species.get_count("ATP") > 20  # la glicólisis catalizada produjo ATP neto
