"""End-to-end test of the stages that need no external tool.

Stages 00, 02 and 03 require RFdiffusion / ProteinMPNN / AlphaFold2, so their outputs
are synthesised here. Stages 01 and 04 then run unmodified, which exercises the metric
layer, the filter engine, the manifest merge and the self-consistency rollup together.

Stage 03 is faked by perturbing the design monomer by a known amount, so the expected
per-backbone success rate is known exactly and can be asserted on.
"""
from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np

from cyclicnano import s01_backbone_filter, s04_validate
from cyclicnano.config import load_config, parse_overrides
from cyclicnano.manifest import RunState
from cyclicnano.pdbio import extract_monomer, read_pdb, write_pdb
from cyclicnano.runner import Runner
from cyclicnano.s00_generate import build_command
from cyclicnano.s02_design import build_commands, parse_fasta
from cyclicnano.s03_fold import build_command as fold_command
from synthetic import clashing_oligomer, cyclic_oligomer

N_GOOD, N_BAD, N_SEQ = 4, 3, 10


class Harness:
    """A temporary run directory pre-populated with synthetic stage-00 output."""

    def __init__(self, rules: str = '["n_clash == 0"]'):
        self.tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_pipeline_"))
        # Rules are stated here rather than inherited: filtering is off in the shipped
        # defaults, so a test that wants rejection must ask for it.
        self.cfg = load_config(overrides=parse_overrides([
            "target.symmetry=C5", "target.monomer_length=60",
            f"run.outdir={self.tmp.as_posix()}",
            f"backbone_filter.rules={rules}",
            "validate.rmsd_cutoff=1.0", "validate.plddt_cutoff=80",
        ]))
        self.state = RunState(self.cfg.outdir)

    def seed_backbones(self):
        outdir = self.state.stage_dir("00_backbones")
        rows = []
        for i in range(N_GOOD):
            p = outdir / f"good_{i:03d}.pdb"
            write_pdb(str(p), cyclic_oligomer(n_sym=5, ring_radius=16.0 + 0.5 * i))
            rows.append({"backbone_id": p.stem, "backbone_path": str(p), "symmetry": "C5"})
        for i in range(N_BAD):
            p = outdir / f"clash_{i:03d}.pdb"
            write_pdb(str(p), clashing_oligomer(n_sym=5))
            rows.append({"backbone_id": p.stem, "backbone_path": str(p), "symmetry": "C5"})
        self.state.backbones.upsert(rows)
        self.state.mark_done("00_backbones")

    def seed_predictions(self, seed: int = 0):
        """Half the sequences per backbone fold back tightly, half do not."""
        rng = np.random.default_rng(seed)
        preddir = self.state.stage_dir("03_predictions")
        df = self.state.backbones.load()
        df = df[df["bb_pass"].astype(bool)]
        rows = []
        for record in df.to_dict("records"):
            bid = record["backbone_id"]
            for k in range(N_SEQ):
                good = k % 2 == 0
                pred = extract_monomer(read_pdb(record["backbone_path"]))
                pred.chains["A"].coords["CA"] = (
                    pred.chains["A"].ca + rng.normal(0.0, 0.3 if good else 3.0,
                                                     pred.chains["A"].ca.shape))
                sid = f"{bid}__seq{k:03d}"
                path = preddir / f"{sid}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
                write_pdb(str(path), pred)
                rows.append({"sequence_id": sid, "backbone_id": bid, "seq_index": k,
                             "sequence": "A" * 60, "pred_path": str(path), "fold_ok": True,
                             "mean_plddt": 92.0 if good else 60.0})
        self.state.sequences.upsert(rows)
        self.state.mark_done("02_sequences")
        self.state.mark_done("03_predictions")

    def close(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


def _run_through_stage04():
    h = Harness()
    h.seed_backbones()
    summary = s01_backbone_filter.run(h.cfg, h.state)
    h.seed_predictions()
    s04_validate.run(h.cfg, h.state)
    return h, summary


# ------------------------------------------------------------------ stage 01
def test_stage01_scores_every_backbone():
    h = Harness()
    h.seed_backbones()
    try:
        summary = s01_backbone_filter.run(h.cfg, h.state)
        assert summary["n_input"] == N_GOOD + N_BAD
        assert summary["n_pass"] == N_GOOD
        assert summary["rejected_by_rule"]["n_clash == 0"] == N_BAD
    finally:
        h.close()


def test_stage01_keeps_rejected_candidates_with_reasons():
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        df = h.state.backbones.load()
        assert len(df) == N_GOOD + N_BAD                    # nothing is ever deleted
        assert df["helix_frac"].notna().all()               # metrics stored for failures too
        rejected = df[~df["bb_pass"].astype(bool)]
        assert rejected["bb_fail_reason"].str.len().gt(0).all()
    finally:
        h.close()


def test_stage01_writes_a_filter_summary():
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        assert (Path(h.cfg.outdir) / "01_backbone_filter" / "filter_summary.json").exists()
    finally:
        h.close()


# ------------------------------------------------------------------ stage 04
def test_stage04_recovers_the_planted_success_rate():
    h, _ = _run_through_stage04()
    try:
        bb = h.state.backbones.load()
        passing = bb[bb["bb_pass"].astype(bool)]
        assert passing["success_rate"].notna().all()
        assert bool(((passing["success_rate"] - 0.5).abs() < 1e-9).all()), \
            passing["success_rate"].tolist()
        assert bool((passing["best_rmsd"] < 1.0).all())
    finally:
        h.close()


def test_stage04_writes_a_ranking():
    h, _ = _run_through_stage04()
    try:
        path = Path(h.cfg.outdir) / "04_validate" / "ranked_backbones.csv"
        assert path.exists()
        import pandas as pd
        ranked = pd.read_csv(path)
        assert list(ranked["success_rate"]) == sorted(ranked["success_rate"], reverse=True)
    finally:
        h.close()


def test_stage04_preserves_stage01_columns():
    h, _ = _run_through_stage04()
    try:
        bb = h.state.backbones.load()
        assert "helix_frac" in bb.columns and "success_rate" in bb.columns
    finally:
        h.close()


def test_stage04_records_alternative_cutoffs():
    h, _ = _run_through_stage04()
    try:
        seqs = h.state.sequences.load()
        assert "sc_pass_at_2A" in seqs.columns               # from validate.also_record
        assert seqs["sc_pass_at_2A"].sum() >= seqs["sc_pass"].sum()
    finally:
        h.close()


# ------------------------------------------------------------------ resume
def test_completed_stage_is_skipped_and_can_be_forced():
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        assert h.state.is_done("01_backbone_filter")
        assert s01_backbone_filter.run(h.cfg, h.state) == {}         # skipped
        assert s01_backbone_filter.run(h.cfg, h.state, force=True)["n_input"] == N_GOOD + N_BAD
    finally:
        h.close()


# ------------------------------------------------------------------ command construction
def _rfd_cfg(*extra):
    return load_config(overrides=parse_overrides([
        "target.symmetry=C5", "target.monomer_length=60", *extra]))


def test_rfdiffusion_command_uses_total_contig_length():
    argv = build_command(_rfd_cfg(), Path("/run/batch000"), 10)
    joined = " ".join(str(a) for a in argv)
    assert "contigmap.contigs=[300-300]" in joined            # 5 x 60, not 60
    assert "--config-name symmetry" in joined


def test_rfdiffusion_symmetry_flag_is_lowercase():
    """RFdiffusion is invoked with c5, matching the hand-run command that works."""
    argv = build_command(_rfd_cfg(), Path("/run/batch000"), 10)
    assert "inference.symmetry=c5" in argv
    assert "inference.symmetry=C5" not in argv


def test_rfdiffusion_command_sets_hydra_run_dir():
    argv = build_command(_rfd_cfg(), Path("/run/batch000"), 10, Path("/run/hydra"))
    assert "hydra.run.dir=/run/hydra" in " ".join(str(a) for a in argv).replace("\\", "/")


def test_rfdiffusion_command_omits_unverified_hydra_keys():
    """An unknown hydra key aborts the run, so nothing speculative is emitted."""
    joined = " ".join(str(a) for a in build_command(_rfd_cfg(), Path("/run/b"), 10, seed=7))
    assert "diffuser.T" not in joined
    assert "inference.seed" not in joined


def test_rfdiffusion_seed_is_opt_in():
    cfg = _rfd_cfg("generate.seed_key=inference.seed")
    assert "inference.seed=7" in build_command(cfg, Path("/run/b"), 10, seed=7)


def test_rfdiffusion_potential_weights_are_written_as_integers():
    """weight_intra:1, not 1.0, matching the known-good invocation."""
    joined = " ".join(str(a) for a in
                      build_command(_rfd_cfg("generate.potentials.enabled=true"),
                                    Path("/run/b"), 10))
    assert "weight_intra:1,weight_inter:0.1" in joined
    assert "potentials.guide_scale=2" in joined
    assert "potentials.guide_scale=2.0" not in joined


def test_rfdiffusion_potentials_can_be_disabled():
    joined = " ".join(str(a) for a in
                      build_command(_rfd_cfg("generate.potentials.enabled=false"),
                                    Path("/run/b"), 10))
    assert "guiding_potentials" not in joined
    assert "contigmap.contigs=[300-300]" in joined


def test_proteinmpnn_command_ties_chains():
    cfg = load_config(overrides=parse_overrides(["paths.proteinmpnn=/opt/ProteinMPNN"]))
    cmds = build_commands(cfg, Path("/in"), Path("/work"))
    joined = " ".join(" ".join(argv) for _, argv in cmds)
    assert "make_tied_positions_dict.py" in joined
    assert "--homooligomer 1" in joined
    assert "--tied_positions_jsonl" in joined


def test_colabfold_command_disables_msa():
    cfg = load_config()
    argv = fold_command(cfg, Path("/in.fasta"), Path("/out"))
    assert "--msa-mode" in argv and "single_sequence" in argv
    assert "--templates" not in argv


def test_runner_sets_deterministic_gpu_ordering():
    env = Runner(kind="local", gpus=[2, 3]).env_vars()
    assert env["CUDA_DEVICE_ORDER"] == "PCI_BUS_ID"
    assert env["CUDA_VISIBLE_DEVICES"] == "2,3"


def test_runner_conda_wrapping():
    r = Runner(kind="conda", envs={"proteinmpnn": "lmk_ProteinMPNN"})
    assert r.wrap("proteinmpnn", ["python", "x.py"])[:4] == ["conda", "run", "-n", "lmk_ProteinMPNN"]


def test_runner_rejects_unconfigured_tool():
    try:
        Runner(kind="conda", envs={}).wrap("proteinmpnn", ["python"])
        raise AssertionError("expected a KeyError")
    except KeyError:
        pass


# ------------------------------------------------------------------ parsing
def test_parse_fasta_skips_the_input_sequence_and_splits_tied_chains():
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "x.fa"
        path.write_text(
            ">input, score=0.0\nAAAA/AAAA\n"
            ">T=0.1, sample=1, score=1.23, global_score=1.50, seq_recovery=0.31\n"
            "MKVL/MKVL\n"
            ">T=0.1, sample=2, score=1.45, global_score=1.60, seq_recovery=0.29\n"
            "MKWL/MKWL\n",
            encoding="utf-8")
        entries = parse_fasta(path)
    assert len(entries) == 2                                  # input entry dropped
    assert entries[0]["sequence"] == "MKVL"                   # tied chains collapsed
    assert entries[0]["chains_identical"] is True
    assert entries[0]["n_chains_in_fasta"] == 2
    assert entries[0]["mpnn_score"] == 1.23
    assert entries[1]["mpnn_seq_recovery"] == 0.29


def test_rfdiffusion_contact_params_are_omitted_unless_set():
    """r_0 and d_0 exist on olig_contacts but stay at RFdiffusion defaults by default."""
    joined = " ".join(str(a) for a in build_command(_rfd_cfg(), Path("/run/b"), 10))
    assert "r_0:" not in joined and "d_0:" not in joined


def test_rfdiffusion_contact_params_are_emitted_when_set():
    cfg = _rfd_cfg("generate.potentials.enabled=true",
                   "generate.potentials.r_0=6", "generate.potentials.d_0=2")
    joined = " ".join(str(a) for a in build_command(cfg, Path("/run/b"), 10))
    assert "weight_intra:1,weight_inter:0.1,r_0:6,d_0:2" in joined


def test_rfdiffusion_never_emits_monomer_rog():
    """monomer_ROG has no chain_lengths and would shrink the whole ring, not the subunits."""
    joined = " ".join(str(a) for a in build_command(_rfd_cfg(), Path("/run/b"), 10))
    assert "monomer_ROG" not in joined


def test_default_command_matches_the_verified_baseline():
    """The shipped default reproduces the hand-run command exactly, nothing added.

    The baseline is the invocation with oligomer-contact potentials, which measurably
    removes backbone clashes and gives more compact subunits than the bare form.
    """
    argv = build_command(_rfd_cfg(), Path("/out/output"), 50, Path("/logs"))
    joined = " ".join(str(a) for a in argv).replace("\\", "/")
    for expected in (
        "--config-name symmetry",
        "inference.symmetry=c5",
        "contigmap.contigs=[300-300]",
        'potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.1"]',
        "potentials.olig_intra_all=True",
        "potentials.olig_inter_all=True",
        "potentials.guide_scale=2",
        "potentials.guide_decay=quadratic",
        "inference.output_prefix=/out/output",
        "hydra.run.dir=/logs",
        "inference.num_designs=50",
    ):
        assert expected in joined, expected
    # still nothing speculative: these keys do not exist in the shipped hydra configs
    for absent in ("diffuser.T", "inference.seed", "monomer_ROG", "guide_scale=2.0"):
        assert absent not in joined, absent


# ------------------------------------------------------------------ stage 03 wiring
def _fold_cfg(*extra):
    return load_config(profile="configs/profiles/amax.yaml",
                       overrides=parse_overrides(list(extra)))


def test_colabfold_command_passes_the_weight_directory():
    """Without --data colabfold_batch looks elsewhere and re-downloads the parameters."""
    argv = fold_command(_fold_cfg(), Path("/in.fasta"), Path("/out"))
    assert "--data" in argv
    assert argv[argv.index("--data") + 1] == "/data/lmk/colabfold_parameters"


def test_colabfold_positionals_come_last():
    """Matches the verified invocation: flags first, then input, then output."""
    argv = [str(a).replace("\\", "/") for a in
            fold_command(_fold_cfg(), Path("/in.fasta"), Path("/out"))]
    assert argv[-2:] == ["/in.fasta", "/out"]
    assert argv[0].endswith("colabfold_batch")


def test_colabfold_omits_unverified_flags_by_default():
    argv = fold_command(_fold_cfg(), Path("/in.fasta"), Path("/out"))
    assert "--model-type" not in argv
    assert "--random-seed" not in argv
    assert "--templates" not in argv


def test_colabfold_optional_flags_are_emitted_when_set():
    argv = fold_command(_fold_cfg("fold.model_type=alphafold2_ptm", "fold.random_seed=0"),
                        Path("/in.fasta"), Path("/out"))
    assert "--model-type" in argv and "--random-seed" in argv


def test_pixi_managed_tool_is_not_wrapped_in_conda_run():
    """LocalColabFold registers no conda env; it is reached through PATH instead."""
    from cyclicnano.runner import runner_from_config
    r = runner_from_config(_fold_cfg(), Path("/tmp"))
    assert r.wrap("colabfold", ["colabfold_batch", "x"]) == ["colabfold_batch", "x"]
    assert r.wrap("proteinmpnn", ["python", "x"])[:3] == ["conda", "run", "-n"]


def test_tool_env_prepends_rather_than_replaces_path():
    from cyclicnano.runner import runner_from_config
    r = runner_from_config(_fold_cfg(), Path("/tmp"))
    env = r.env_vars("colabfold")
    assert env["PATH"].startswith("/data/lmk/localcolabfold/.pixi/envs/default/bin")
    assert len(env["PATH"]) > len("/data/lmk/localcolabfold/.pixi/envs/default/bin")
    assert env["LD_LIBRARY_PATH"].startswith("/data/lmk/localcolabfold/.pixi/envs/default/lib")


def test_tool_env_is_scoped_to_its_tool():
    from cyclicnano.runner import runner_from_config
    r = runner_from_config(_fold_cfg(), Path("/tmp"))
    assert "/localcolabfold/" not in r.env_vars("rfdiffusion").get("PATH", "")


def test_tool_with_neither_conda_env_nor_tool_env_is_rejected():
    from cyclicnano.runner import Runner
    try:
        Runner(kind="conda", envs={}, tool_env={}).wrap("colabfold", ["x"])
        raise AssertionError("expected a KeyError")
    except KeyError:
        pass


def test_dry_run_leaves_no_completion_marker():
    """A dry run that marked stages done would make the following real run skip them."""
    import tempfile
    from cyclicnano.manifest import RunState
    from cyclicnano.runner import Runner
    from cyclicnano import s00_generate
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_dry_"))
    cfg = load_config(overrides=parse_overrides([
        f"run.outdir={tmp.as_posix()}", "target.n_backbones=2"]))
    state = RunState(cfg.outdir)
    s00_generate.run(cfg, state, Runner(kind="local", dry_run=True))
    assert not state.is_done("00_backbones")


def test_filtering_is_off_in_the_shipped_defaults():
    """Every backbone must reach sequence design unless rules are added explicitly."""
    cfg = load_config()
    assert cfg.get("backbone_filter.rules") == []


def test_metrics_are_still_recorded_when_nothing_is_filtered():
    """Turning filtering off must not turn measurement off."""
    h = Harness(rules="[]")
    h.seed_backbones()
    try:
        summary = s01_backbone_filter.run(h.cfg, h.state)
        assert summary["n_pass"] == N_GOOD + N_BAD          # nothing rejected
        df = h.state.backbones.load()
        assert df["bb_pass"].astype(bool).all()
        for column in ("helix_frac", "n_clash", "shape_anisotropy", "sym_order_detected"):
            assert df[column].notna().all(), column          # measured anyway
    finally:
        h.close()


# ------------------------------------------------------------------ output hygiene
def test_rfdiffusion_command_disables_trajectory_output():
    """Trajectory files were 98 percent of stage 00 output and nothing reads them."""
    argv = build_command(_rfd_cfg(), Path("/run/b"), 10)
    assert "inference.write_trajectory=False" in argv


def test_trajectory_output_can_be_re_enabled():
    cfg = _rfd_cfg("generate.write_trajectory=true")
    assert "inference.write_trajectory=False" not in build_command(cfg, Path("/run/b"), 10)


def test_proteinmpnn_intermediates_are_written_outside_the_results_dir():
    cfg = load_config(overrides=parse_overrides(["paths.proteinmpnn=/opt/ProteinMPNN"]))
    cmds = build_commands(cfg, Path("/in"), Path("/work"), Path("/logs"))
    joined = " ".join(" ".join(str(a) for a in argv) for _, argv in cmds).replace("\\", "/")
    assert "/logs/parsed_chains.jsonl" in joined
    assert "/logs/tied_positions.jsonl" in joined
    assert "/work/parsed_chains.jsonl" not in joined


# ------------------------------------------------------------------ stage 05
def _run_through_delivery():
    from cyclicnano import s05_delivery
    h, _ = _run_through_stage04()
    s05_delivery.run(h.cfg, h.state)
    return h, Path(h.cfg.outdir) / "05_delivery"


def test_delivery_writes_a_pdb_and_a_fasta_per_passing_design():
    h, out = _run_through_delivery()
    try:
        bb = h.state.backbones.load()
        n = int(bb["final_pass"].fillna(False).astype(bool).sum())
        assert n > 0
        assert len(list(out.glob("*.pdb"))) == n
        assert len(list(out.glob("*.fasta"))) == n
        assert (out / "delivery.csv").exists()
    finally:
        h.close()


def test_delivery_excludes_designs_that_failed_validation():
    h, out = _run_through_delivery()
    try:
        bb = h.state.backbones.load()
        failed = bb[~bb["final_pass"].fillna(False).astype(bool)]["backbone_id"]
        for bid in failed:
            assert not (out / f"{bid}.pdb").exists()
    finally:
        h.close()


def test_delivery_fasta_carries_the_best_sequence():
    """The delivered sequence must be the one that folded back closest."""
    h, out = _run_through_delivery()
    try:
        seqs = h.state.sequences.load()
        import pandas as pd
        table = pd.read_csv(out / "delivery.csv")
        for row in table.to_dict("records"):
            group = seqs[seqs.backbone_id == row["backbone_id"]]
            assert abs(row["best_rmsd"] - group["sc_rmsd"].min()) < 1e-9
            text = (out / f"{row['backbone_id']}.fasta").read_text()
            assert row["sequence"] in text
            assert text.startswith(">")
    finally:
        h.close()


def test_delivery_is_sorted_best_first():
    h, out = _run_through_delivery()
    try:
        import pandas as pd
        table = pd.read_csv(out / "delivery.csv")
        assert list(table["success_rate"]) == sorted(table["success_rate"], reverse=True)
    finally:
        h.close()


def test_delivery_rerun_clears_stale_hits():
    from cyclicnano import s05_delivery
    h, out = _run_through_delivery()
    try:
        (out / "zz_stale.pdb").write_text("junk", encoding="utf-8")
        s05_delivery.run(h.cfg, h.state, force=True)
        assert not (out / "zz_stale.pdb").exists()
    finally:
        h.close()


def test_delivery_requires_a_stage_04_verdict():
    h = Harness()
    h.seed_backbones()
    try:
        from cyclicnano import s05_delivery
        s01_backbone_filter.run(h.cfg, h.state)
        try:
            s05_delivery.run(h.cfg, h.state)
            raise AssertionError("expected a RuntimeError")
        except RuntimeError:
            pass
    finally:
        h.close()


def test_rfdiffusion_batches_continue_each_others_numbering():
    """Batch two starts where batch one stopped, instead of restarting at 0."""
    argv = [str(a) for a in build_command(_rfd_cfg(), Path("/out/design"), 10, start=10)]
    assert "inference.design_startnum=10" in argv
    assert "inference.output_prefix=/out/design" in " ".join(argv).replace("\\", "/")


def test_rfdiffusion_startnum_can_be_pinned_in_config():
    cfg = _rfd_cfg("generate.design_startnum=20")
    assert "inference.design_startnum=20" in build_command(cfg, Path("/out/design"), 10, start=0)


def test_runner_forces_a_headless_plotting_backend():
    """colabfold plots with matplotlib; an inherited DISPLAY makes the job die on logout."""
    import os
    from cyclicnano.runner import Runner
    os.environ["DISPLAY"] = "localhost:10.0"
    try:
        env = Runner(kind="local").env_vars("colabfold")
        assert env["MPLBACKEND"] == "Agg"
        assert "DISPLAY" not in env
    finally:
        os.environ.pop("DISPLAY", None)


def test_headless_backend_applies_to_every_tool():
    from cyclicnano.runner import Runner
    for tool in ("rfdiffusion", "proteinmpnn", "colabfold"):
        assert Runner(kind="local").env_vars(tool)["MPLBACKEND"] == "Agg"


def test_rfdiffusion_is_seeded_so_a_run_can_be_repeated():
    """Unseeded, the same command returns a different hundred backbones every time."""
    assert "inference.deterministic=True" in build_command(_rfd_cfg(), Path("/run/b"), 10)


def test_deterministic_seeding_can_be_turned_off():
    cfg = _rfd_cfg("generate.deterministic=false")
    assert "inference.deterministic=True" not in build_command(cfg, Path("/run/b"), 10)


def test_each_batch_gets_a_distinct_seed_range():
    """Seeding is per design index, so batches must not restart it or they repeat."""
    first = build_command(_rfd_cfg(), Path("/run/design"), 10, start=0)
    second = build_command(_rfd_cfg(), Path("/run/design"), 10, start=10)
    assert "inference.design_startnum=0" in first
    assert "inference.design_startnum=10" in second
    assert all("inference.deterministic=True" in c for c in (first, second))
