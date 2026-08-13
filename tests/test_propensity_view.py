"""
tests/test_propensity_view.py

Validación de engine/propensity_view.py.
"""

import pytest
from engine.species import SpeciesManager
from engine.reactions import Reaction, ReactionManager
from engine.ssa import GillespieSSA
from engine.debugger import Debugger
from engine.propensity_view import (
    PropensityBar,
    prepare_propensity_bars,
    render_propensity_bars_as_text,
    show_propensities,
)


def test_prepare_bars_sorted_by_propensity_descending():
    bars = prepare_propensity_bars({"slow": 1.0, "fast": 10.0, "medium": 5.0})
    assert [b.reaction_id for b in bars] == ["fast", "medium", "slow"]


def test_prepare_bars_fractions_sum_to_one():
    bars = prepare_propensity_bars({"a": 3.0, "b": 1.0})
    total_fraction = sum(b.fraction for b in bars)
    assert total_fraction == pytest.approx(1.0)
    assert bars[0].fraction == pytest.approx(0.75)
    assert bars[1].fraction == pytest.approx(0.25)


def test_prepare_bars_handles_zero_total_without_division_error():
    bars = prepare_propensity_bars({"a": 0.0, "b": 0.0})
    assert all(b.fraction == 0.0 for b in bars)


def test_prepare_bars_empty_input():
    bars = prepare_propensity_bars({})
    assert bars == []


def test_render_text_contains_all_reaction_ids():
    bars = prepare_propensity_bars({"transcribe_x": 4.0, "degrade_x": 1.0})
    text = render_propensity_bars_as_text(bars)
    assert "transcribe_x" in text
    assert "degrade_x" in text


def test_render_text_highlights_the_fired_reaction():
    bars = prepare_propensity_bars({"a": 2.0, "b": 1.0})
    text = render_propensity_bars_as_text(bars, highlight_reaction_id="b")
    lines = text.splitlines()
    b_line = next(line for line in lines if line.startswith("b "))
    a_line = next(line for line in lines if line.startswith("a "))
    assert "disparó" in b_line
    assert "disparó" not in a_line


def test_render_text_empty_bars_returns_placeholder():
    text = render_propensity_bars_as_text([])
    assert text == "(sin reacciones)"


def test_render_text_bar_width_is_respected():
    bars = prepare_propensity_bars({"only": 5.0})
    text = render_propensity_bars_as_text(bars, bar_width=10)
    # peso 100% -> la barra debe estar completamente rellena (10 símbolos '#')
    assert "#" * 10 in text


def test_show_propensities_without_rich_falls_back_to_text_and_prints(capsys):
    propensities = {"a": 2.0, "b": 1.0}
    returned_text = show_propensities(propensities, use_rich=False)

    captured = capsys.readouterr()
    assert "a" in captured.out
    assert "b" in captured.out
    assert returned_text == render_propensity_bars_as_text(prepare_propensity_bars(propensities))


# ---------------------------------------------------------------------
# Integración real con el depurador
# ---------------------------------------------------------------------

def test_integration_with_debugger_pending_propensities():
    sm = SpeciesManager()
    sm.add_species("A", 10)
    sm.add_species("B", 10)
    sm.add_species("C", 0)
    rm = ReactionManager()
    rm.add_reaction(Reaction("A_to_C", {"A": 1}, {"C": 1}, rate_constant=3.0))
    rm.add_reaction(Reaction("B_to_C", {"B": 1}, {"C": 1}, rate_constant=1.0))
    ssa = GillespieSSA(sm, rm, seed=1)
    debugger = Debugger(ssa)

    event = debugger.step()[0]
    propensities = debugger.pending_propensities()

    bars = prepare_propensity_bars(propensities)
    assert {b.reaction_id for b in bars} == {"A_to_C", "B_to_C"}

    text = render_propensity_bars_as_text(bars, highlight_reaction_id=event.reaction_id)
    assert event.reaction_id in text
    assert "disparó" in text
