"""
tests/test_replication.py

Validación de biology/replication.py.
"""

import pytest
from cell import Cell
from engine.reactions import Reaction
from data_io.genome_loader import GenomeData, GeneAnnotation
from biology.genome import GenomeModule
from biology.translation import TranslationModule
from biology.replication import (
    ReplicationModule,
    ReplicationModuleError,
    total_count_condition,
)


def test_total_count_condition_triggers_at_threshold():
    condition = total_count_condition(["A", "B"], threshold=10)
    assert condition({"A": 4, "B": 5}) is False  # 9 < 10
    assert condition({"A": 5, "B": 5}) is True   # 10 >= 10
    assert condition({"A": 20, "B": 0}) is True


def test_total_count_condition_ignores_missing_species():
    condition = total_count_condition(["A", "B"], threshold=5)
    assert condition({"A": 5}) is True  # "B" ausente cuenta como 0


def test_total_count_condition_invalid_threshold_raises():
    with pytest.raises(ReplicationModuleError):
        total_count_condition(["A"], threshold=0)


def test_should_divide_reflects_condition():
    cell = Cell(seed=1)
    cell.add_species("Prot_X", 5)
    module = ReplicationModule(total_count_condition(["Prot_X"], threshold=10), seed=1)

    assert module.should_divide(cell) is False

    cell.species.change_count("Prot_X", 10)
    assert module.should_divide(cell) is True


def test_check_and_divide_returns_none_if_not_ready():
    cell = Cell(seed=1)
    cell.add_species("Prot_X", 1)
    module = ReplicationModule(total_count_condition(["Prot_X"], threshold=100), seed=1)

    result = module.check_and_divide(cell)
    assert result is None


def test_divide_conserves_total_count_per_species():
    cell = Cell(seed=1)
    cell.add_species("A", 100)
    cell.add_species("B", 37)
    module = ReplicationModule(lambda s: True, seed=42)

    daughter_a, daughter_b = module.divide(cell)

    assert daughter_a.species.get_count("A") + daughter_b.species.get_count("A") == 100
    assert daughter_a.species.get_count("B") + daughter_b.species.get_count("B") == 37


def test_divide_zero_count_species_stays_zero_in_both_daughters():
    cell = Cell(seed=1)
    cell.add_species("A", 0)
    module = ReplicationModule(lambda s: True, seed=1)

    daughter_a, daughter_b = module.divide(cell)

    assert daughter_a.species.get_count("A") == 0
    assert daughter_b.species.get_count("A") == 0


def test_divide_split_is_roughly_balanced_over_many_species():
    """
    No se comprueba un valor exacto (es estocástico), sino que el
    reparto de una cantidad grande no está sistemáticamente sesgado
    hacia una de las dos hijas.
    """
    cell = Cell(seed=1)
    cell.add_species("A", 100_000)
    module = ReplicationModule(lambda s: True, seed=7)

    daughter_a, daughter_b = module.divide(cell)

    fraction_a = daughter_a.species.get_count("A") / 100_000
    assert 0.45 < fraction_a < 0.55  # debería rondar 0.5 con una muestra tan grande


def test_daughters_inherit_the_same_reactions():
    cell = Cell(seed=1)
    cell.add_species("A", 10)
    cell.add_species("B", 0)
    cell.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))
    module = ReplicationModule(lambda s: True, seed=1)

    daughter_a, daughter_b = module.divide(cell)

    assert len(daughter_a.reactions) == 1
    assert len(daughter_b.reactions) == 1
    assert daughter_a.reactions.get_reaction("A_to_B") is daughter_b.reactions.get_reaction("A_to_B")


def test_daughters_are_independent_cells():
    cell = Cell(seed=1)
    cell.add_species("A", 100)
    module = ReplicationModule(lambda s: True, seed=1)

    daughter_a, daughter_b = module.divide(cell)
    count_a_before = daughter_a.species.get_count("A")

    daughter_b.species.change_count("A", -1)  # solo debería afectar a daughter_b

    assert daughter_a.species.get_count("A") == count_a_before


def test_replicate_gene_copies_doubles_the_count():
    cell = Cell(seed=1)
    cell.add_species("gene_target", 1)

    ReplicationModule.replicate_gene_copies(cell, ["gene_target"])

    assert cell.species.get_count("gene_target") == 2


def test_divide_records_history():
    cell = Cell(seed=1)
    cell.add_species("A", 10)
    module = ReplicationModule(lambda s: True, seed=1)

    assert len(module.history) == 0
    module.divide(cell)
    assert len(module.history) == 1
    assert module.history[0].daughter_a_state["A"] + module.history[0].daughter_b_state["A"] == 10


# ---------------------------------------------------------------------
# Ciclo celular completo: transcripción + traducción -> división ->
# ambas hijas siguen funcionando de forma independiente
# ---------------------------------------------------------------------

def test_full_cell_cycle_until_division_and_independent_daughters():
    genome = GenomeData(
        sequence_id="toy",
        sequence="A" * 100,
        genes=[GeneAnnotation("target", 1, 90, "+", promoter_id="strong")],
    )

    cell = Cell(seed=17)
    genome_module = GenomeModule(genome)
    installed_genes = genome_module.install(cell, rnap_initial_count=10)

    translation_module = TranslationModule(installed_genes)
    translation_module.install(cell, ribosome_initial_count=10)

    replication_module = ReplicationModule(
        total_count_condition(["Prot_target"], threshold=30), seed=99
    )

    # crecer hasta que se cumpla la condición de división
    while not replication_module.should_divide(cell):
        events = cell.run(max_steps=50)
        if not events:  # protección: no debería quedarse atascada
            break

    assert replication_module.should_divide(cell)

    # replicar el gen antes de dividir, para que ambas hijas puedan
    # seguir transcribiéndolo de forma independiente
    replication_module.replicate_gene_copies(cell, ["gene_target"])
    assert cell.species.get_count("gene_target") == 2

    daughter_a, daughter_b = replication_module.divide(cell)

    # ambas hijas deben poder seguir corriendo su propia simulación,
    # de forma independiente, sin heredar el histórico de la madre
    assert len(daughter_a.get_events()) == 0
    assert len(daughter_b.get_events()) == 0

    events_a = daughter_a.run(max_steps=20)
    events_b = daughter_b.run(max_steps=20)

    # si alguna hija heredó al menos una copia del gen, debería poder
    # seguir produciendo ARNm de forma independiente
    if daughter_a.species.get_count("gene_target") > 0:
        assert len(events_a) > 0
    if daughter_b.species.get_count("gene_target") > 0:
        assert len(events_b) > 0
