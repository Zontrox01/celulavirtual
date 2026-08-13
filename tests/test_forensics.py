"""
tests/test_forensics.py

Validación de engine/forensics.py.
"""

import pytest
from engine.reactions import Reaction, ReactionManager
from engine.ssa import SSAEvent
from engine.forensics import ForensicsAnalyzer, ForensicReport


def make_event(step, time, reaction_id):
    return SSAEvent(step=step, time=time, reaction_id=reaction_id, propensities={}, total_propensity=0.0)


def build_simple_reactions():
    rm = ReactionManager()
    rm.add_reaction(Reaction("produce_X", {}, {"X": 1}, rate_constant=1.0))
    rm.add_reaction(Reaction("consume_X", {"X": 1}, {}, rate_constant=1.0))
    rm.add_reaction(Reaction("unrelated", {"Y": 1}, {"Z": 1}, rate_constant=1.0))
    return rm


def test_events_in_window_filters_by_time():
    events = [make_event(1, 1.0, "a"), make_event(2, 5.0, "b"), make_event(3, 9.0, "c")]
    analyzer = ForensicsAnalyzer(events, ReactionManager())

    result = analyzer.events_in_window(2.0, 6.0)
    assert [e.reaction_id for e in result] == ["b"]


def test_events_in_window_is_inclusive_at_boundaries():
    events = [make_event(1, 1.0, "a"), make_event(2, 5.0, "b")]
    analyzer = ForensicsAnalyzer(events, ReactionManager())

    result = analyzer.events_in_window(1.0, 5.0)
    assert len(result) == 2


def test_last_event_before_returns_the_most_recent_one():
    events = [make_event(1, 1.0, "a"), make_event(2, 3.0, "b"), make_event(3, 7.0, "c")]
    analyzer = ForensicsAnalyzer(events, ReactionManager())

    result = analyzer.last_event_before(5.0)
    assert result.reaction_id == "b"


def test_last_event_before_returns_none_if_no_events_yet():
    analyzer = ForensicsAnalyzer([make_event(1, 10.0, "a")], ReactionManager())
    assert analyzer.last_event_before(5.0) is None


def test_reactions_fired_count_in_window():
    events = [
        make_event(1, 1.0, "a"), make_event(2, 2.0, "a"),
        make_event(3, 3.0, "b"), make_event(4, 20.0, "a"),  # fuera de ventana
    ]
    analyzer = ForensicsAnalyzer(events, ReactionManager())

    counts = analyzer.reactions_fired_count(0.0, 10.0)
    assert counts == {"a": 2, "b": 1}


def test_analyze_species_net_change_matches_direct_computation():
    events = [
        make_event(1, 1.0, "produce_X"),
        make_event(2, 2.0, "produce_X"),
        make_event(3, 3.0, "consume_X"),
    ]
    analyzer = ForensicsAnalyzer(events, build_simple_reactions())

    report = analyzer.analyze_species("X", 0.0, 10.0)
    assert report.net_change == 1  # +1 +1 -1


def test_analyze_species_ignores_reactions_that_do_not_touch_the_species():
    events = [
        make_event(1, 1.0, "produce_X"),
        make_event(2, 2.0, "unrelated"),  # no toca X
    ]
    analyzer = ForensicsAnalyzer(events, build_simple_reactions())

    report = analyzer.analyze_species("X", 0.0, 10.0)
    reaction_ids = {c.reaction_id for c in report.contributions}
    assert "unrelated" not in reaction_ids
    assert reaction_ids == {"produce_X"}


def test_analyze_species_top_producer_and_consumer():
    events = (
        [make_event(i, float(i), "produce_X") for i in range(1, 6)]      # +5
        + [make_event(i, float(i), "consume_X") for i in range(6, 8)]    # -2
    )
    analyzer = ForensicsAnalyzer(events, build_simple_reactions())

    report = analyzer.analyze_species("X", 0.0, 10.0)

    assert report.net_change == 3
    assert report.top_producer().reaction_id == "produce_X"
    assert report.top_producer().net_change == 5
    assert report.top_consumer().reaction_id == "consume_X"
    assert report.top_consumer().net_change == -2


def test_analyze_species_catalytic_reaction_has_zero_net_and_is_excluded():
    rm = ReactionManager()
    rm.add_reaction(Reaction("bind_and_release", {"R": 1}, {"R": 1}, rate_constant=1.0))  # neto 0 sobre R
    events = [make_event(1, 1.0, "bind_and_release")] * 3
    analyzer = ForensicsAnalyzer(events, rm)

    report = analyzer.analyze_species("R", 0.0, 10.0)
    assert report.contributions == []
    assert report.net_change == 0


def test_analyze_species_no_top_producer_or_consumer_when_none_exist():
    report = ForensicReport("X", 0.0, 1.0, net_change=0, contributions=[])
    assert report.top_producer() is None
    assert report.top_consumer() is None


def test_unknown_reaction_in_log_is_ignored_gracefully():
    events = [make_event(1, 1.0, "produce_X"), make_event(2, 2.0, "reaccion_fantasma")]
    analyzer = ForensicsAnalyzer(events, build_simple_reactions())

    report = analyzer.analyze_species("X", 0.0, 10.0)  # no debe lanzar excepción
    assert report.net_change == 1


# ---------------------------------------------------------------------
# Integración real: transcripción + degradación
# ---------------------------------------------------------------------

def test_forensics_identifies_producer_and_consumer_in_real_simulation():
    from cell import Cell
    from data_io.genome_loader import GenomeData, GeneAnnotation
    from biology.genome import GenomeModule
    from biology.degradation import DegradationModule

    genome = GenomeData(
        sequence_id="toy", sequence="A" * 100,
        genes=[GeneAnnotation("target", 1, 90, "+", promoter_id="strong")],
    )
    cell = Cell(seed=11)
    genome_module = GenomeModule(genome)
    installed_genes = genome_module.install(cell, rnap_initial_count=10)
    degradation_module = DegradationModule(installed_genes)
    degradation_module.install(cell, mrna_half_life=10.0)

    cell.run(max_time=100.0)

    analyzer = ForensicsAnalyzer(cell.get_events(), cell.reactions)
    report = analyzer.analyze_species("mRNA_target", 0.0, cell.time)

    producer = report.top_producer()
    consumer = report.top_consumer()

    assert producer is not None and producer.reaction_id == "transcribe_target"
    assert consumer is not None and consumer.reaction_id == "degrade_mRNA_target"
    assert report.net_change == cell.species.get_count("mRNA_target")
