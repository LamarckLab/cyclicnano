"""Tests against real RFdiffusion output.

Everything else in the suite runs on synthetic structures. These three fixtures are
genuine C5 backbones from RFdiffusion (60 aa subunits, symmetry config), so they pin
the parser and the metric layer against the actual file format rather than against
assumptions about it.

They are named for their fold, because that is what the stage-01 rules select on and
each one is here to exercise a different verdict:

    alpha_beta    two helices and two strands. Rejected: it carries beta.
    two_helix     all helix, but a hairpin of a 26-mer and a 30-mer. Rejected: this
                  is the shape that dominates the generated C5 set and the one the
                  protocol rules out, the helices being too long to design against.
    three_helix   three helices, no beta. Accepted. One backbone in the hundred
                  generated looked like this, which is why a run generates 200.
"""
from __future__ import annotations

from pathlib import Path

from cyclicnano.config import load_config
from cyclicnano.filters import FilterSet
from cyclicnano.geometry import backbone_metrics
from cyclicnano.pdbio import read_pdb

DATA = Path(__file__).resolve().parent / "data"
ALPHA_BETA_PDB = DATA / "rfdiffusion_c5_60aa_alpha_beta.pdb"
TWO_HELIX_PDB = DATA / "rfdiffusion_c5_60aa_two_helix.pdb"
THREE_HELIX_PDB = DATA / "rfdiffusion_c5_60aa_three_helix.pdb"


def test_real_output_parses_into_five_chains():
    """RFdiffusion symmetric mode writes one chain per symmetry copy, not one long chain."""
    struct = read_pdb(str(ALPHA_BETA_PDB))
    assert struct.n_chains == 5
    assert struct.chain_ids == ["A", "B", "C", "D", "E"]
    for cid in struct.chain_ids:
        assert struct.chains[cid].n_res == 60


def test_real_output_residue_numbering_is_continuous_across_chains():
    """Numbering runs 1..300 rather than restarting per chain; the parser must not assume."""
    struct = read_pdb(str(ALPHA_BETA_PDB))
    assert int(struct.chains["A"].resnums[0]) == 1
    assert int(struct.chains["E"].resnums[-1]) == 300


def test_real_output_has_only_backbone_atoms():
    struct = read_pdb(str(ALPHA_BETA_PDB))
    assert set(struct.chains["A"].coords) == {"N", "CA", "C", "O"}


def test_symmetry_is_recovered_from_real_output():
    m = backbone_metrics(read_pdb(str(ALPHA_BETA_PDB)), expected_sym=5)
    assert abs(m["sym_order_detected"] - 5) < 0.01
    assert m["sym_rmsd"] < 0.01
    assert m["sym_order_ok"] is True


def test_metrics_on_real_output_are_in_plausible_ranges():
    m = backbone_metrics(read_pdb(str(ALPHA_BETA_PDB)), expected_sym=5)
    assert m["n_res_monomer"] == 60
    assert 0.0 <= m["helix_frac"] <= 1.0
    assert abs(m["helix_frac"] + m["strand_frac"] + m["loop_frac"] - 1.0) < 1e-9
    assert 0.5 < m["rg_ratio"] < 2.5
    assert 0.0 <= m["shape_anisotropy"] <= 1.0
    assert 0.0 <= m["assembly_shape_anisotropy"] <= 1.0
    assert m["iface_contacts_per_chain"] > 0


def _preset_rules():
    cfg = load_config(profile="configs/presets/pure_three_helix.yaml")
    return FilterSet("backbone", cfg.get("backbone_filter.rules"))


def test_each_preset_rule_rejects_on_its_own_grounds():
    """One fixture per verdict, and the rejections do not overlap by accident."""
    rules = _preset_rules()

    strandy = backbone_metrics(read_pdb(str(ALPHA_BETA_PDB)), expected_sym=5)
    assert strandy["strand_frac"] > 0
    assert "strand_frac == 0" in rules.apply(strandy).failed_rules

    hairpin = backbone_metrics(read_pdb(str(TWO_HELIX_PDB)), expected_sym=5)
    assert hairpin["strand_frac"] == 0                  # the beta rule has nothing to say
    assert rules.apply(hairpin).failed_rules == ["n_helices == 3"]

    bundle = backbone_metrics(read_pdb(str(THREE_HELIX_PDB)), expected_sym=5)
    assert rules.apply(bundle).passed, bundle["ss_string"]


def test_default_rules_accept_real_output():
    """The shipped defaults must not reject ordinary RFdiffusion output wholesale."""
    cfg = load_config()
    rules = FilterSet("backbone", cfg.get("backbone_filter.rules"))
    assert rules.apply(backbone_metrics(read_pdb(str(ALPHA_BETA_PDB)), expected_sym=5)).passed


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
              (ALPHA_BETA_PDB).read_text(encoding="utf-8").splitlines()
              if ln.startswith("ATOM")]
    pred = [ln[17:20] for ln in PRED_PDB.read_text(encoding="utf-8").splitlines()
            if ln.startswith("ATOM")]
    assert set(design) == {"GLY"}
    assert len(set(pred)) > 5


def test_self_consistency_rmsd_runs_on_real_files():
    from cyclicnano.s04_validate import monomer_rmsd
    assert monomer_rmsd(PRED_PDB, PRED_PDB) < 1e-9   # identical structures
    value = monomer_rmsd(ALPHA_BETA_PDB, PRED_PDB)
    assert value > 0                                 # unrelated design and prediction


# ------------------------------------------------------------ pure-three-helix preset
def test_preset_parses_against_real_metrics():
    _preset_rules().check_syntax(
        backbone_metrics(read_pdb(str(ALPHA_BETA_PDB)), expected_sym=5))


def test_preset_accepts_only_the_three_helix_fixture():
    rules = _preset_rules()
    verdicts = {p.stem: rules.apply(backbone_metrics(read_pdb(str(p)), expected_sym=5)).passed
                for p in (ALPHA_BETA_PDB, TWO_HELIX_PDB, THREE_HELIX_PDB)}
    assert sum(verdicts.values()) == 1, verdicts
    assert verdicts[THREE_HELIX_PDB.stem]


def test_three_helix_fixture_really_is_three_helices_and_no_beta():
    """Guards the fixture itself: if the SS assignment drifts, the preset tests above
    would start passing or failing for a reason that has nothing to do with the rules."""
    m = backbone_metrics(read_pdb(str(THREE_HELIX_PDB)), expected_sym=5)
    assert m["strand_frac"] == 0.0
    assert m["n_helices"] == 3
    assert m["helix_frac"] > 0.8
    assert m["ss_string"].count("E") == 0


def test_preset_does_not_filter_on_clashes():
    """Stage 03 settles whether a close contact matters; screening twice loses designs."""
    assert not any("clash" in r for r in _preset_rules().rules)


def test_preset_does_not_filter_on_size():
    """The radius ceiling was dropped: the only qualifying C5 backbone measured 27.9,
    above every threshold under consideration, so keeping it took the yield to zero."""
    rules = _preset_rules().rules
    assert not any("centre_distance" in r or "height" in r for r in rules)
    m = backbone_metrics(read_pdb(str(THREE_HELIX_PDB)), expected_sym=5)
    assert m["max_centre_distance"] > 27                 # and it is still measured
