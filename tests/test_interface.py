"""Stage 05: interface evaluation.

AlphaFold2-multimer and PyRosetta are not available to the suite, so the prediction
and the energy are synthesised. What is tested here is everything around them: which
sequences are chosen, how the oligomer is spelled for colabfold, how the results are
carried into the delivery, and that none of it touches the verdicts.
"""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from cyclicnano import s04_delivery, s05_interface
from cyclicnano.config import load_config, parse_overrides
from cyclicnano.runner import Runner


def _frames():
    """Four backbones, two of which passed stage 04, three sequences each."""
    backbones = pd.DataFrame([
        {"backbone_id": "design_1", "bb_pass": True, "final_pass": True},
        {"backbone_id": "design_2", "bb_pass": True, "final_pass": True},
        {"backbone_id": "design_3", "bb_pass": True, "final_pass": False},
        {"backbone_id": "design_4", "bb_pass": False, "final_pass": False},
    ])
    rows = []
    for bid in ("design_1", "design_2", "design_3", "design_4"):
        for k, passing in enumerate((True, False, True)):
            rows.append({"sequence_id": f"{bid}__seq{k:03d}", "backbone_id": bid,
                         "seq_index": k, "sequence": "AC" * 30, "sc_pass": passing})
    return backbones, pd.DataFrame(rows)


def test_only_passing_sequences_of_passing_backbones_are_evaluated():
    """Both filters matter. Every sequence of a passing backbone would spend most of
    the compute on sequences already known not to fold back; every passing sequence
    regardless of backbone would describe designs that were never delivered."""
    backbones, sequences = _frames()
    chosen = {r["sequence_id"] for r in s05_interface.targets(backbones, sequences)}
    assert chosen == {"design_1__seq000", "design_1__seq002",
                      "design_2__seq000", "design_2__seq002"}
    assert not any(s.startswith("design_3") or s.startswith("design_4") for s in chosen)


def test_evaluation_needs_a_stage_04_verdict():
    backbones, sequences = _frames()
    try:
        s05_interface.targets(backbones.drop(columns=["final_pass"]), sequences)
        raise AssertionError("expected a RuntimeError")
    except RuntimeError as exc:
        assert "stage 04" in str(exc)


def test_the_oligomer_is_spelled_with_one_chain_per_symmetry_copy():
    """colabfold is told to fold a homo-oligomer by repeating the chain, so the count
    has to come from the symmetry of the run: three for a C3, five for a C5."""
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_iface_"))
    try:
        rows = [{"sequence_id": "design_1__seq000", "sequence": "AAA"},
                {"sequence_id": "design_1__seq002", "sequence": "CCC"}]
        text = s05_interface.multimer_fasta(rows, 5, tmp / "x.fasta").read_text()
        assert text.splitlines()[1] == "AAA:AAA:AAA:AAA:AAA"
        assert text.count(">") == 2

        text = s05_interface.multimer_fasta(rows, 3, tmp / "y.fasta").read_text()
        assert text.splitlines()[1] == "AAA:AAA:AAA"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _cfg(*extra):
    return load_config(overrides=parse_overrides(
        ["target.symmetry=C5", "target.monomer_length=60", *extra]))


def test_the_prediction_uses_the_multimer_model():
    """ipTM is not defined for one chain, so stage 03's model cannot produce it."""
    argv = [str(a) for a in s05_interface.build_command(_cfg(), Path("/in.fasta"),
                                                        Path("/out"))]
    assert "--model-type" in argv
    assert argv[argv.index("--model-type") + 1] == "alphafold2_multimer_v3"
    assert "single_sequence" in argv              # still no MSA; five chains change nothing


def test_the_prediction_model_can_be_overridden():
    argv = [str(a) for a in s05_interface.build_command(
        _cfg("interface.model_type=alphafold2_multimer_v2"), Path("/in.fasta"), Path("/out"))]
    assert argv[argv.index("--model-type") + 1] == "alphafold2_multimer_v2"


def test_predictions_are_regrouped_under_short_names():
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_iface_"))
    try:
        sid = "design_1__seq000"
        pdb = tmp / f"{sid}_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"
        pdb.write_text("ATOM", encoding="utf-8")
        (tmp / f"{sid}_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
         ).write_text("{}", encoding="utf-8")

        out = s05_interface.regroup(tmp, [{"sequence_id": sid, "interface_pred_path": str(pdb)}])
        assert sorted(p.name for p in (tmp / "design_1").iterdir()) == [
            "seq000_af2multimer.pdb", "seq000_af2multimer_scores.json"]
        assert Path(out[0]["interface_pred_path"]).name == "seq000_af2multimer.pdb"
        assert list(tmp.glob("*.pdb")) == []
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_ipTM_is_read_from_the_multimer_scores():
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_iface_"))
    try:
        sid = "design_1__seq000"
        (tmp / f"{sid}_unrelaxed_rank_001_x.pdb").write_text("ATOM", encoding="utf-8")
        (tmp / f"{sid}_scores_rank_001_x.json").write_text(
            json.dumps({"iptm": 0.88, "ptm": 0.9, "plddt": [90.0, 92.0]}), encoding="utf-8")

        row = s05_interface.collect_prediction(tmp, sid)
        assert row["interface_ok"] is True
        assert row["iptm"] == 0.88
        assert row["multimer_ptm"] == 0.9
        assert abs(row["multimer_plddt"] - 91.0) < 1e-9

        absent = s05_interface.collect_prediction(tmp, "design_9__seq000")
        assert absent["interface_ok"] is False
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_the_energy_is_one_interface_between_two_chains():
    """Chains A and B, taken out of the relaxed oligomer so that they sit in the
    conformation the assembly puts them in. Chain A against the other four would be
    the cost of removing a subunit from the ring, which counts both of its contacts
    and measured 2.34 times one interface on a real structure."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("iface_worker", s05_interface.WORKER)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)            # no pyrosetta import at module level

    source = s05_interface.WORKER.read_text(encoding="utf-8")
    assert "two chains" in (worker.score_one.__doc__ or "")
    # the partner is found, not assumed to be chain 2
    assert "contact_partner" in source
    assert "pair_pose(pose, 1, partner)" in source


def test_the_contacting_partner_is_not_assumed_to_be_chain_two():
    """AlphaFold numbers chains as it pleases. On design_205__seq001, a prediction
    that failed to close into a ring (detected C2.53, ipTM 0.56), chain A sat 12.7 A
    from B and 5.9 A from C, so scoring chains 1 and 2 measured the gap across the
    ring and returned dSASA 0. That happens on exactly the predictions worth reading.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location("iface_worker2", s05_interface.WORKER)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)

    class FakePose:
        """Chain 1 near chain 3, far from 2: the design_205__seq001 arrangement."""
        coords = {1: [(0, 0, 0)], 2: [(20, 0, 0)], 3: [(6, 0, 0)], 4: [(30, 0, 0)]}

        def num_chains(self):
            return 4

        def chain_begin(self, c):
            return c

        def chain_end(self, c):
            return c

        def residue(self, r):
            xyz = self.coords[r][0]
            return type("R", (), {"xyz": lambda self, _n, v=xyz: v})()

    partner, gap = worker.contact_partner(FakePose())
    assert partner == 3, partner
    assert abs(gap - 6.0) < 1e-6


def test_an_unevaluated_sequence_gets_a_null_block_not_a_block_of_nulls():
    """The distinction is the point: not measured must not read as measured and poor.

    The row comes back from the table, not from the stage, so an unevaluated sequence
    carries NaN rather than False. bool(NaN) is True in Python, so the obvious test
    passed for exactly the rows it existed to catch, and all three unevaluated
    sequences of design_27 were handed a block of nulls.
    """
    import math

    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_iface_"))
    try:
        assert s04_delivery._interface_block({"sequence_id": "x"}, tmp) is None
        assert s04_delivery._interface_block(
            {"sequence_id": "x", "interface_scored": float("nan"),
             "iptm": float("nan")}, tmp) is None
        assert s04_delivery._interface_block(
            {"sequence_id": "x", "interface_scored": False, "iptm": None}, tmp) is None

        scored = {"sequence_id": "x", "interface_scored": True, "iptm": 0.88,
                  "dG_interface": -38.5, "dSASA_interface": 1849.0}
        block = s04_delivery._interface_block(scored, tmp)
        assert block["iptm"] == 0.88
        assert block["dG_interface"] == -38.5
        assert block["prediction"] is None            # nothing was copied into the folder
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_a_dry_run_shows_both_commands_and_writes_nothing():
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_iface_"))
    try:
        import contextlib
        import io
        from cyclicnano.manifest import RunState

        cfg = _cfg(f"run.outdir={tmp.as_posix()}", "backend.kind=local")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            s05_interface.run(cfg, RunState(cfg.outdir), Runner(dry_run=True), dry_run=True)
        out = buf.getvalue()
        assert "alphafold2_multimer_v3" in out, out
        assert "interface_energy.py" in out, out
        assert not (tmp / ".done").exists()
        assert not (tmp / "05_interface").exists()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_the_worker_script_is_shipped_and_standalone():
    """It runs inside the PyRosetta environment, where cyclicnano is not installed."""
    assert s05_interface.WORKER.exists(), s05_interface.WORKER
    source = s05_interface.WORKER.read_text(encoding="utf-8")
    assert "import cyclicnano" not in source
    assert "from cyclicnano" not in source


def test_interface_metrics_are_not_part_of_any_rule():
    """The whole point: these describe a design, they do not select one."""
    cfg = _cfg()
    rules = (cfg.get("backbone_filter.rules") or []) + (cfg.get("validate.rules") or [])
    for metric in ("iptm", "dG_interface", "dSASA_interface", "dG_dSASA_100"):
        assert not any(metric in rule for rule in rules), (metric, rules)
