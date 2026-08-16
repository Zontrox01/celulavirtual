"""
tests/test_app_state.py

Validación de gui/app_state.py. Ningún test de este archivo necesita
PyQt6: toda la lógica de la GUI vive separada de Qt, precisamente para
poder testearla igual que el resto del proyecto.
"""

import pytest
from pathlib import Path
from gui.app_state import AppState, AppStateError, GenomeSetupOptions
from data_io.genome_loader import load_genome
from tools.generate_synthetic_genome import generate_synthetic_genome


def test_require_cell_raises_before_any_genome_loaded():
    state = AppState()
    with pytest.raises(AppStateError):
        state.require_cell()


def test_load_and_install_basic(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=1)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=1))

    assert state.cell is not None
    assert len(state.installed_genes) == 4
    assert len(state.translated_genes) == 4
    # transcripción + traducción + degradación(ARNm) + degradación(proteína), por gen
    assert len(state.cell.reactions) == 4 * 4


def test_run_steps(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=2)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=2))
    events = state.run_steps(20)
    assert len(events) == 20


def test_species_snapshot_contains_machinery(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=3)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=3))
    snapshot = state.species_snapshot()
    assert "RNAP" in snapshot
    assert "Ribosome" in snapshot


def test_trajectory_matches_events_run(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=4)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=4))
    state.run_steps(10)
    df = state.trajectory()
    assert len(df) == 11  # fila inicial + 10 eventos


def test_metabolism_background_with_atp_coupling(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=5)
    state = AppState()
    options = GenomeSetupOptions(
        fasta_path=fasta, annotation_path=csv, seed=5,
        metabolism_mode="background_with_atp_coupling",
        atp_cost_per_transcription=2, atp_cost_per_translation=3,
    )
    state.load_and_install_genome(options)

    assert state.metabolism_module is not None
    assert state.cell.species.has_species("ATP")
    reaction = state.cell.reactions.get_reaction(state.installed_genes[0].transcription_reaction_id)
    assert reaction.reactants.get("ATP") == 2


def test_metabolism_enzyme_catalyzed(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=6)
    genome = load_genome(fasta, csv)
    enzyme_gene = genome.genes[0].gene_id

    state = AppState()
    options = GenomeSetupOptions(
        fasta_path=fasta, annotation_path=csv, seed=6,
        metabolism_mode="enzyme_catalyzed", enzyme_gene_id=enzyme_gene,
    )
    state.load_and_install_genome(options)

    assert state.metabolism_module is not None
    installation = state.metabolism_module.get_installation()
    assert installation.enzyme_species_id == f"Prot_{enzyme_gene}"


def test_metabolism_enzyme_catalyzed_without_gene_id_raises(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=7)
    state = AppState()
    options = GenomeSetupOptions(
        fasta_path=fasta, annotation_path=csv, seed=7,
        metabolism_mode="enzyme_catalyzed", enzyme_gene_id=None,
    )
    with pytest.raises(AppStateError):
        state.load_and_install_genome(options)


def test_invalid_metabolism_mode_raises(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=8)
    state = AppState()
    options = GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, metabolism_mode="no_existe")
    with pytest.raises(AppStateError):
        state.load_and_install_genome(options)


def test_add_repression(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=9)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=9))
    gene_ids = state.gene_ids()

    link = state.add_repression(gene_ids[0], gene_ids[1], binding_rate=2.0, unbinding_rate=0.5)
    assert link.target_gene_id == gene_ids[0]
    assert link.repressor_gene_id == gene_ids[1]


def test_add_repression_unknown_gene_raises(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=10)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=10))
    with pytest.raises(AppStateError):
        state.add_repression("fantasma", state.gene_ids()[0])


def test_new_debugger_and_forensics_integration(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=11)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=11))

    debugger = state.new_debugger(enable_snapshots=True, snapshot_interval=5)
    debugger.step(n=30)

    target_species = f"mRNA_{state.gene_ids()[0]}"
    report = state.analyze_species(target_species, 0.0, state.cell.time)
    assert report.species_id == target_species


def test_divide_and_adopt(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=12)
    state = AppState()
    state.load_and_install_genome(
        GenomeSetupOptions(
            fasta_path=fasta, annotation_path=csv, seed=12,
            rnap_initial_count=10, ribosome_initial_count=10,
        )
    )
    protein_species = [f"Prot_{g}" for g in state.gene_ids()]
    gene_species = [f"gene_{g}" for g in state.gene_ids()]

    for _ in range(200):
        total = sum(state.cell.species.get_count(p) for p in protein_species)
        if total >= 5:
            break
        state.run_steps(20)

    daughter_a, daughter_b = state.divide(protein_species, 5, gene_species)
    assert daughter_a is not None
    assert daughter_b is not None

    state.adopt_daughter(daughter_a)
    assert state.cell is daughter_a
    assert state.generation == 1


def test_divide_before_threshold_raises(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=13)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=13))
    protein_species = [f"Prot_{g}" for g in state.gene_ids()]
    gene_species = [f"gene_{g}" for g in state.gene_ids()]

    with pytest.raises(AppStateError):
        state.divide(protein_species, 1_000_000, gene_species)


def test_export_and_import_sbml_round_trip(tmp_path):
    fasta, csv = generate_synthetic_genome(4, tmp_path, seed=14)
    state = AppState()
    state.load_and_install_genome(GenomeSetupOptions(fasta_path=fasta, annotation_path=csv, seed=14))
    state.run_steps(10)

    sbml_path = Path(tmp_path) / "export.xml"
    state.export_sbml(sbml_path)
    assert sbml_path.exists()

    reimported = AppState()
    reimported.import_sbml(sbml_path, seed=1)
    assert reimported.cell is not None
    assert reimported.genome_module is None  # documentado: se pierde la info de módulos
    events = reimported.run_steps(5)
    assert len(events) == 5
