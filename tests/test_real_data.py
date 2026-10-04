"""Tests against real RFdiffusion output.

Everything else in the suite runs on synthetic structures. These two fixtures are
genuine C5 backbones from RFdiffusion (60 aa subunits, symmetry config, no guiding
potentials), so they pin the parser and the metric layer against the actual file
format rather than against assumptions about it.

Both files were produced by the calibration run described in
configs/presets/compact_subunit.yaml.
"""
from __future__ import annotations

from pathlib import Path

from cyclicnano.config import load_config
from cyclicnano.filters import FilterSet
from cyclicnano.geometry import backbone_metrics
from cyclicnano.pdbio import read_pdb

DATA = Path(__file__).resolve().parent / "data"
PASS_PDB = DATA / "rfdiffusion_c5_60aa_pass.pdb"
FAIL_PDB = DATA / "rfdiffusion_c5_60aa_fail.pdb"


def test_real_output_parses_into_five_chains():
    """RFdiffusion symmetric mode writes one chain per symmetry copy, not one long chain."""
    struct = read_pdb(str(PASS_PDB))
    assert struct.n_chains == 5
    assert struct.chain_ids == ["A", "B", "C", "D", "E"]
    for cid in struct.chain_ids:
        assert struct.chains[cid].n_res == 60


def test_real_output_residue_numbering_is_continuous_across_chains():
    """Numbering runs 1..300 rather than restarting per chain; the parser must not assume."""
    struct = read_pdb(str(PASS_PDB))
    assert int(struct.chains["A"].resnums[0]) == 1
    assert int(struct.chains["E"].resnums[-1]) == 300


def test_real_output_has_only_backbone_atoms():
    struct = read_pdb(str(PASS_PDB))
    assert set(struct.chains["A"].coords) == {"N", "CA", "C", "O"}


def test_symmetry_is_recovered_from_real_output():
    m = backbone_metrics(read_pdb(str(PASS_PDB)), expected_sym=5)
    assert abs(m["sym_order_detected"] - 5) < 0.01
    assert m["sym_rmsd"] < 0.01
    assert m["sym_order_ok"] is True


def test_metrics_on_real_output_are_in_plausible_ranges():
    m = backbone_metrics(read_pdb(str(PASS_PDB)), expected_sym=5)
    assert m["n_res_monomer"] == 60
    assert 0.0 <= m["helix_frac"] <= 1.0
    assert abs(m["helix_frac"] + m["strand_frac"] + m["loop_frac"] - 1.0) < 1e-9
    assert 0.5 < m["rg_ratio"] < 2.5
    assert 0.0 <= m["shape_anisotropy"] <= 1.0
    assert 0.0 <= m["assembly_shape_anisotropy"] <= 1.0
    assert m["iface_contacts_per_chain"] > 0


def test_preset_separates_the_two_real_fixtures():
    """The calibrated preset must accept the good backbone and reject the poor one."""
    cfg = load_config(profile="configs/presets/compact_subunit.yaml")
    rules = FilterSet("backbone", cfg.get("backbone_filter.rules"))
    good = backbone_metrics(read_pdb(str(PASS_PDB)), expected_sym=5)
    poor = backbone_metrics(read_pdb(str(FAIL_PDB)), expected_sym=5)
    rules.check_syntax(good)
    assert rules.apply(good).passed, rules.apply(good).failed_rules
    assert not rules.apply(poor).passed


def test_default_rules_accept_real_output():
    """The shipped defaults must not reject ordinary RFdiffusion output wholesale."""
    cfg = load_config()
    rules = FilterSet("backbone", cfg.get("backbone_filter.rules"))
    assert rules.apply(backbone_metrics(read_pdb(str(PASS_PDB)), expected_sym=5)).passed


# ---------------------------------------------------------------- ProteinMPNN output
MPNN_FA = DATA / "proteinmpnn_c5_tied.fa"


def test_real_fasta_skips_the_input_sequence():
    """Entry 0 is the poly-glycine input backbone and must not be treated as a design."""
    from cyclicnano.s02_design import parse_fasta
    entries = parse_fasta(MPNN_FA)
    assert len(entries) == 4
    assert not any(set(e["sequence"]) == {"G"} for e in entries)


def test_real_fasta_confirms_tied_positions_took_effect():
    """All five chains must carry the same sequence, or the oligomer cannot assemble."""
    from cyclicnano.s02_design import parse_fasta
    for e in parse_fasta(MPNN_FA):
        assert e["n_chains_in_fasta"] == 5
        assert e["chains_identical"] is True
        assert len(e["sequence"]) == 60          # monomer, not the 300-residue concatenation


def test_real_fasta_metadata_is_parsed():
    """The header fields survive the comma split, including on the first design entry."""
    from cyclicnano.s02_design import parse_fasta
    first = parse_fasta(MPNN_FA)[0]
    assert first["mpnn_score"] == 0.6490
    assert first["mpnn_global_score"] == 0.6490
    # seq_recovery is 0 by construction: RFdiffusion writes poly-glycine backbones
    assert first["mpnn_seq_recovery"] == 0.0


def test_real_fasta_header_brackets_do_not_break_parsing():
    """The skipped input header contains designed_chains=['A', 'B', ...] with commas."""
    from cyclicnano.s02_design import parse_fasta
    assert "designed_chains=['A', 'B'" in MPNN_FA.read_text(encoding="utf-8")
    assert len(parse_fasta(MPNN_FA)) == 4        # parsed anyway


# ---------------------------------------------------------------- AlphaFold2 output
PRED_ID = "batch000_5__seq000"
PRED_PDB = DATA / f"{PRED_ID}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
SCORES = DATA / f"{PRED_ID}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json"


def test_collect_matches_the_real_colabfold_filenames():
    """The identifiers carry a double underscore; colabfold keeps it, so the glob must too."""
    from cyclicnano.s03_fold import collect
    assert "__" in PRED_ID
    row = collect(DATA, PRED_ID)
    assert row["fold_ok"] is True
    assert Path(row["pred_path"]).name == PRED_PDB.name


def test_collect_reads_the_confidence_scores():
    from cyclicnano.s03_fold import collect
    row = collect(DATA, PRED_ID)
    assert abs(row["mean_plddt"] - 94.5157) < 1e-3
    assert abs(row["min_plddt"] - 72.75) < 1e-6
    assert abs(row["ptm"] - 0.73) < 1e-6


def test_collect_reports_a_missing_prediction_instead_of_raising():
    from cyclicnano.s03_fold import collect
    row = collect(DATA, "no_such_sequence")
    assert row["fold_ok"] is False and row["pred_path"] is None


def test_real_scores_json_has_one_plddt_per_residue():
    import json
    data = json.loads(SCORES.read_text(encoding="utf-8"))
    assert set(data) >= {"plddt", "ptm"}
    assert len(data["plddt"]) == 60                  # matches the 60-residue monomer


def test_real_prediction_is_a_single_chain_with_full_backbone():
    """AlphaFold2 predicts the monomer only; the design oligomer has five chains."""
    struct = read_pdb(str(PRED_PDB))
    assert struct.n_chains == 1
    assert struct.chains["A"].n_res == 60
    assert set(struct.chains["A"].coords) >= {"N", "CA", "C", "O"}


def test_real_prediction_has_real_residue_names_unlike_the_design():
    """The design backbone is poly-glycine, which is why sequence alignment cannot pair them."""
    design = [ln[17:20] for ln in
              (DATA / "rfdiffusion_c5_60aa_pass.pdb").read_text(encoding="utf-8").splitlines()
              if ln.startswith("ATOM")]
    pred = [ln[17:20] for ln in PRED_PDB.read_text(encoding="utf-8").splitlines()
            if ln.startswith("ATOM")]
    assert set(design) == {"GLY"}
    assert len(set(pred)) > 5


def test_self_consistency_rmsd_runs_on_real_files():
    from cyclicnano.s04_validate import monomer_rmsd
    assert monomer_rmsd(PRED_PDB, PRED_PDB) < 1e-9   # identical structures
    value = monomer_rmsd(DATA / "rfdiffusion_c5_60aa_pass.pdb", PRED_PDB)
    assert value > 0                                 # unrelated design and prediction


# ---------------------------------------------------------------- axis-aligned preset
def test_axis_aligned_preset_parses_against_real_metrics():
    from cyclicnano.config import load_config
    from cyclicnano.filters import FilterSet
    cfg = load_config(profile="configs/presets/axis_aligned.yaml")
    rules = FilterSet("backbone", cfg.get("backbone_filter.rules"))
    rules.check_syntax(backbone_metrics(read_pdb(str(PASS_PDB)), expected_sym=5))


def test_axis_aligned_preset_rejects_a_strand_heavy_backbone():
    """Both shipped fixtures are alpha/beta; one carries enough strand to be rejected."""
    from cyclicnano.config import load_config
    from cyclicnano.filters import FilterSet
    cfg = load_config(profile="configs/presets/axis_aligned.yaml")
    rules = FilterSet("backbone", cfg.get("backbone_filter.rules"))
    verdicts = [rules.apply(backbone_metrics(read_pdb(str(p)), expected_sym=5)).passed
                for p in (PASS_PDB, FAIL_PDB)]
    assert not all(verdicts)


def test_axis_aligned_preset_does_not_filter_on_clashes():
    """Stage 03 settles whether a close contact matters; screening twice loses designs."""
    from cyclicnano.config import load_config
    rules = load_config(profile="configs/presets/axis_aligned.yaml").get("backbone_filter.rules")
    assert not any("clash" in r for r in rules)
