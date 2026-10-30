"""End-to-end test of the stages that need no external tool.

Stages 00, 02 and 03 require RFdiffusion / ProteinMPNN / AlphaFold2, so their outputs
are synthesised here. Stages 01 and 04 then run unmodified, which exercises the metric
layer, the filter engine, the manifest merge and the self-consistency rollup together.

Stage 03 is faked by perturbing the design monomer by a known amount, so the expected
per-backbone success rate is known exactly and can be asserted on.
"""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

from cyclicnano import s01_backbone_filter, s03_fold, s04_delivery, s04_validate
from cyclicnano.config import load_config, parse_overrides
from cyclicnano.manifest import RunState
from cyclicnano.pdbio import extract_monomer, read_pdb, write_pdb
from cyclicnano.runner import Runner
from cyclicnano import s00_generate
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
            p = outdir / s00_generate.backbone_filename(f"good_{i:03d}")
            write_pdb(str(p), cyclic_oligomer(n_sym=5, ring_radius=16.0 + 0.5 * i))
            rows.append({"backbone_id": s00_generate.backbone_id_from(p),
                         "backbone_path": str(p), "symmetry": "C5"})
        for i in range(N_BAD):
            p = outdir / s00_generate.backbone_filename(f"clash_{i:03d}")
            write_pdb(str(p), clashing_oligomer(n_sym=5))
            rows.append({"backbone_id": s00_generate.backbone_id_from(p),
                         "backbone_path": str(p), "symmetry": "C5"})
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
                # The layout stage 03 leaves behind: one folder per backbone, two
                # files per sequence, the scores beside the structure.
                folder = preddir / bid
                folder.mkdir(parents=True, exist_ok=True)
                path = folder / f"seq{k:03d}{s03_fold.PRED_SUFFIX}"
                write_pdb(str(path), pred)
                plddt = 92.0 if good else 60.0
                s03_fold.scores_path_for(path).write_text(
                    json.dumps({"plddt": [plddt] * 60, "ptm": 0.8}), encoding="utf-8")
                rows.append({"sequence_id": sid, "backbone_id": bid, "seq_index": k,
                             "sequence": "A" * 60, "pred_path": str(path), "fold_ok": True,
                             "mean_plddt": plddt})
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
        path = Path(h.cfg.outdir) / "01_backbone_filter" / "filter_summary.json"
        assert path.exists()
        summary = json.loads(path.read_text(encoding="utf-8"))
        # The rules are restated beside the results: a reader of the output tree
        # should not have to find the config that produced it.
        assert summary["rules"] == ["n_clash == 0"]
        assert len(summary["passed_backbones"]) == summary["n_pass"] == N_GOOD
    finally:
        h.close()


def test_a_backbone_file_says_it_is_a_backbone_at_every_stage():
    """The id stays design_N; only the file carries the tag, so that a structure sitting
    beside an AlphaFold prediction says which of the two it is. The round trip through
    ProteinMPNN is the part that can go wrong: it names its FASTA after the PDB it read,
    so the tag would otherwise ride into every sequence id."""
    from cyclicnano.s00_generate import BACKBONE_TAG, backbone_filename, backbone_id_from

    assert backbone_filename("design_27") == "design_27_bb.pdb"
    assert backbone_id_from("design_27_bb.pdb") == "design_27"
    assert backbone_id_from(Path("/x/design_27_bb.fa")) == "design_27"   # ProteinMPNN
    assert backbone_id_from("design_27.pdb") == "design_27"                    # untagged
    assert BACKBONE_TAG in backbone_filename("design_27")

    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        staged = s01_backbone_filter.passed_dir(h.state)
        names = sorted(p.name for p in staged.glob("*.pdb"))
        assert names and all(n.endswith("_bb.pdb") for n in names), names
        # the ids in the table are untagged
        bb = h.state.backbones.load()
        ids = bb["backbone_id"].tolist()
        assert not any(str(b).endswith(BACKBONE_TAG) for b in ids), ids
    finally:
        h.close()


def test_stage01_copies_the_survivors_into_their_own_folder():
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        staged = s01_backbone_filter.passed_dir(h.state)
        names = sorted(s00_generate.backbone_id_from(p) for p in staged.glob("*.pdb"))
        assert names == sorted(h.state.backbones.load()
                               .query("bb_pass")["backbone_id"].tolist())
        assert len(names) == N_GOOD
    finally:
        h.close()


def test_stage01_rebuilds_the_survivor_folder_when_the_rules_tighten():
    """A backbone the new rules reject must not be left behind for stage 02 to pick
    up: the folder is the handover, so a stale file there is a design nobody chose."""
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        staged = s01_backbone_filter.passed_dir(h.state)
        assert len(list(staged.glob("*.pdb"))) == N_GOOD

        h.cfg.data["backbone_filter"]["rules"] = ["n_clash == 0", "n_res_monomer == 61"]
        s01_backbone_filter.run(h.cfg, h.state, force=True)
        assert list(staged.glob("*.pdb")) == []
    finally:
        h.close()


def test_sequence_files_say_what_made_them_and_what_they_hold():
    """ProteinMPNN names its FASTA after the PDB it read, so the file lands as
    design_N_bb.fa. Left alone the backbone tag would ride into every sequence
    id; the rename swaps it for the one describing the file, and the id is unchanged
    on both sides."""
    from cyclicnano.s02_design import (SEQ_TAG, backbone_id_from_fasta,
                                       sequences_filename)

    assert sequences_filename("design_27") == "design_27_mpnn_seqs.fa"
    assert backbone_id_from_fasta("design_27_mpnn_seqs.fa") == "design_27"
    # ProteinMPNN's own output name, before the rename
    assert backbone_id_from_fasta("design_27_bb.fa") == "design_27"
    assert backbone_id_from_fasta(Path("/x/design_27.fa")) == "design_27"
    assert "mpnn" in SEQ_TAG


def test_stage02_refuses_to_run_without_the_survivor_folder():
    """It reads the folder, not the table; the two could otherwise disagree silently."""
    from cyclicnano import s02_design
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        shutil.rmtree(s01_backbone_filter.passed_dir(h.state))
        try:
            s02_design.run(h.cfg, h.state, Runner())
            raise AssertionError("expected a RuntimeError")
        except RuntimeError as exc:
            assert "stage 01" in str(exc)
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


def test_a_nonzero_success_rate_already_bounds_the_best_rmsd():
    """Why the shipped validate rules carry no best_rmsd threshold.

    A sequence passes only when its RMSD is at or below validate.rmsd_cutoff, so the
    minimum across a backbone is below that as soon as one sequence passes, and a rule
    restating the cutoff can reject nothing. Measured on real output it rejected none
    of 16 C5 and none of 45 C3 backbones, and on C3 it would have admitted 36 where
    the success rate admitted 27.

    The redundancy is invisible while the two numbers agree, which is the reason this
    is a test and not a comment: raising rmsd_cutoff would turn such a rule into a real
    and quite different constraint without anyone deciding that.
    """
    h, _ = _run_through_stage04()
    try:
        cutoff = float(h.cfg.get("validate.rmsd_cutoff"))
        bb = h.state.backbones.load()
        scored = bb[bb["success_rate"].notna() & (bb["success_rate"] > 0)]
        assert len(scored), "the harness produced no scored backbone"
        assert (scored["best_rmsd"] <= cutoff).all(), scored[["backbone_id", "best_rmsd"]]

        shipped = load_config().get("validate.rules", []) or []
        assert not any("best_rmsd" in rule for rule in shipped), shipped
    finally:
        h.close()


def test_best_plddt_follows_the_best_rmsd_sequence_not_the_loudest_one():
    """best_plddt is the confidence OF the best sequence, not the best confidence of
    any sequence. The two are different numbers, and the second flatters a failure:
    on real output the maximum was usually set by a sequence AlphaFold was sure about
    and had folded wrong, 95.3 at 8.8 A. Sitting beside best_rmsd and best_sequence_id
    it read as though all three described one sequence; now they do."""
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        h.seed_predictions()

        # a sequence AlphaFold is very sure about. The harness gives it a large
        # perturbation, so it is nowhere near the best by RMSD.
        seqs = h.state.sequences.load()
        bid = _filtered_ids(h)[0]
        loud = seqs[(seqs.backbone_id == bid) & (seqs.seq_index == 1)].iloc[0]
        h.state.sequences.upsert([{"sequence_id": loud.sequence_id, "mean_plddt": 99.0}])

        s04_validate.run(h.cfg, h.state, force=True)

        seqs = h.state.sequences.load()
        group = seqs[seqs.backbone_id == bid].dropna(subset=["sc_rmsd"])
        closest = group.nsmallest(1, "sc_rmsd").iloc[0]
        row = h.state.backbones.load().set_index("backbone_id").loc[bid]

        assert closest.sequence_id != loud.sequence_id, "the planted row must not be the best"
        assert abs(row["best_plddt"] - closest.mean_plddt) < 1e-9, row["best_plddt"]
        assert row["best_plddt"] != 99.0, "it took the maximum instead of the best"

        # and the delivery row keeps all three columns on that one sequence
        s04_delivery.run(h.cfg, h.state, force=True)
        import pandas as pd
        table = pd.read_csv(Path(h.cfg.outdir) / "04_delivery" / "delivery.csv")
        d = table[table.backbone_id == bid].iloc[0]
        assert d.best_sequence_id == closest.sequence_id
        assert abs(d.best_rmsd - closest.sc_rmsd) < 1e-9
        assert abs(d.best_plddt - closest.mean_plddt) < 1e-9
    finally:
        h.close()


def test_stage04_writes_no_directory_of_its_own():
    """Validation only adds columns. A second ranking file beside the delivery table
    would be the same numbers under another name, free to disagree with it."""
    h, _ = _run_through_stage04()
    try:
        assert not (Path(h.cfg.outdir) / "04_validate").exists()
        bb = h.state.backbones.load()
        assert "final_pass" in bb.columns and "success_rate" in bb.columns
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


def test_noise_scale_is_absent_unless_asked_for():
    """Unset, the command must stay byte-identical to the verified baseline."""
    joined = " ".join(str(a) for a in build_command(_rfd_cfg(), Path("/run/b"), 10))
    assert "noise_scale" not in joined


def test_noise_scale_sets_translations_and_rotations_together():
    """RFdiffusion exposes the two halves separately and advises tying them, so the
    config carries one number and the stage emits both flags from it."""
    argv = [str(a) for a in build_command(_rfd_cfg("generate.noise_scale=0.5"),
                                          Path("/run/b"), 10)]
    assert "denoiser.noise_scale_ca=0.5" in argv
    assert "denoiser.noise_scale_frame=0.5" in argv


def test_noise_scale_of_zero_is_emitted_rather_than_treated_as_unset():
    """0 is a meaningful setting, and `if noise:` would silently drop it."""
    argv = [str(a) for a in build_command(_rfd_cfg("generate.noise_scale=0"),
                                          Path("/run/b"), 10)]
    assert "denoiser.noise_scale_ca=0" in argv


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


def test_stage03_regroups_colabfold_output_by_backbone():
    """colabfold writes everything flat under names that repeat the same rank, model
    and seed on every file. Only rank 1 of one model is produced here, so those 41
    characters distinguish nothing while the directory reaches 640 entries."""
    from cyclicnano.s03_fold import regroup

    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_regroup_"))
    try:
        rows, backbone_of = [], {}
        for bid in ("design_7", "design_12"):
            for k in range(3):
                sid = f"{bid}__seq{k:03d}"
                pdb = tmp / f"{sid}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
                pdb.write_text("ATOM", encoding="utf-8")
                (tmp / f"{sid}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json"
                 ).write_text("{}", encoding="utf-8")
                rows.append({"sequence_id": sid, "pred_path": str(pdb)})
                backbone_of[sid] = bid
        # a sequence that never folded must survive the pass untouched
        rows.append({"sequence_id": "design_7__seq099", "pred_path": None})
        backbone_of["design_7__seq099"] = "design_7"

        out = regroup(tmp, rows, backbone_of)

        assert sorted(d.name for d in tmp.iterdir() if d.is_dir()) == ["design_12", "design_7"]
        assert sorted(p.name for p in (tmp / "design_7").iterdir()) == [
            "seq000_af2.pdb", "seq000_af2_scores.json",
            "seq001_af2.pdb", "seq001_af2_scores.json",
            "seq002_af2.pdb", "seq002_af2_scores.json"]
        assert list(tmp.glob("*.pdb")) == []                 # nothing left flat
        for row in out:
            if row["sequence_id"].endswith("seq099"):
                assert row["pred_path"] is None
            else:
                assert Path(row["pred_path"]).exists()       # rows point at the new place
                assert Path(row["pred_path"]).parent.name in ("design_7", "design_12")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


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


# ------------------------------------------------------------------ stage 04b
def _run_through_delivery():
    h, _ = _run_through_stage04()
    s04_delivery.run(h.cfg, h.state)
    return h, Path(h.cfg.outdir) / "04_delivery"


def _filtered_ids(h):
    bb = h.state.backbones.load()
    return bb[bb["bb_pass"].fillna(False).astype(bool)]["backbone_id"].tolist()


def _delivered(out: Path) -> dict:
    """Map backbone_id to its folder across both groups, rank prefix stripped."""
    found = {}
    for group, strip in ((s04_delivery.PASSED, True), (s04_delivery.BELOW, False)):
        for d in sorted((out / group).glob("*")):
            if d.is_dir():
                found[d.name.split("_", 1)[1] if strip else d.name] = d
    return found


def test_delivery_writes_one_folder_per_backbone_that_passed_the_shape_filter():
    h, out = _run_through_delivery()
    try:
        expected = _filtered_ids(h)
        assert expected
        assert sorted(_delivered(out)) == sorted(expected)
        assert (out / "delivery.csv").exists()
        assert (out / "delivery.json").exists()
    finally:
        h.close()


def test_delivery_splits_the_folders_by_verdict_and_numbers_the_passing_ones():
    """A directory listing should answer which cleared stage 04, and in what order,
    without opening the table."""
    h, out = _run_through_delivery()
    try:
        import pandas as pd
        table = pd.read_csv(out / "delivery.csv")
        passed = table[table["final_pass"]]["backbone_id"].tolist()
        assert passed

        names = sorted(d.name for d in (out / s04_delivery.PASSED).glob("*") if d.is_dir())
        assert len(names) == len(passed)
        # the prefix is the rank, counting from 1 in the order of the table
        for i, bid in enumerate(passed, start=1):
            assert (out / s04_delivery.PASSED / f"{i}_{bid}").is_dir(), names
        # the folder column resolves to the real directory, which is what scripts use
        for row in table.to_dict("records"):
            assert (out / row["folder"]).is_dir(), row["folder"]
    finally:
        h.close()


def test_delivery_table_carries_only_the_delivery_columns():
    """The stage-01 metrics and the looser success rates live in the backbone table.
    Copying them here made a row too wide to read, and a second copy that can age."""
    h, out = _run_through_delivery()
    try:
        import pandas as pd
        cols = set(pd.read_csv(out / "delivery.csv").columns)
        assert cols == {"backbone_id", "folder", "symmetry", "n_res_monomer", "n_seqs",
                        "n_predictions", "n_interface_evaluated", "success_rate",
                        "best_rmsd", "best_plddt", "final_pass", "final_fail_reason",
                        "best_sequence_id", "best_sequence"}, sorted(cols)
        # the looser cutoffs are recorded, but this table answers the 1.0 A question
        assert not any(c.startswith("success_rate_at_") for c in cols), sorted(cols)
        # ... and they are still recorded, one level up
        bb = set(h.state.backbones.load().columns)
        for dropped in ("n_helices", "strand_frac", "max_centre_distance"):
            assert dropped in bb, dropped
    finally:
        h.close()


def test_delivery_excludes_backbones_the_shape_filter_rejected():
    """Stage 01 rejects on geometry and that verdict is final; stage 04 does not."""
    h, out = _run_through_delivery()
    try:
        bb = h.state.backbones.load()
        delivered = _delivered(out)
        for bid in bb[~bb["bb_pass"].fillna(False).astype(bool)]["backbone_id"]:
            assert str(bid) not in delivered
    finally:
        h.close()


def test_delivery_keeps_the_near_misses_and_records_the_verdict():
    """A backbone that misses the RMSD cutoff is still laid out, flagged rather than
    dropped: the cutoffs are uncalibrated and a near miss is worth looking at."""
    h, _ = _run_through_stage04()
    try:
        demoted = _filtered_ids(h)[0]
        h.state.backbones.upsert([{"backbone_id": demoted, "final_pass": False,
                                   "final_fail_reason": "planted"}])
        s04_delivery.run(h.cfg, h.state)
        out = Path(h.cfg.outdir) / "04_delivery"
        assert (out / s04_delivery.BELOW / demoted).is_dir()

        import pandas as pd
        table = pd.read_csv(out / "delivery.csv").set_index("backbone_id")
        assert not bool(table.loc[demoted, "final_pass"])
        assert table["final_pass"].sum() == len(table) - 1
        assert pd.read_csv(out / "delivery.csv")["backbone_id"].iloc[-1] == demoted
    finally:
        h.close()


def test_delivery_folder_holds_the_backbone_and_every_prediction():
    h, out = _run_through_delivery()
    try:
        bid = _filtered_ids(h)[0]
        folder = _delivered(out)[bid]
        assert (folder / s00_generate.backbone_filename(bid)).exists()
        assert (folder / "sequences.json").exists()
        assert (folder / "summary.csv").exists()
        preds = folder / "predictions"
        assert len(list(preds.glob("*.pdb"))) == N_SEQ          # all of them, not the best
        assert len(list(preds.glob("*_af2_scores.json"))) == N_SEQ

        import pandas as pd
        per_seq = pd.read_csv(folder / "summary.csv")
        assert len(per_seq) == N_SEQ
        assert per_seq["sc_rmsd"].notna().all()
    finally:
        h.close()


def test_delivery_keeps_the_names_stage_03_gave_the_predictions():
    """The two layouts should stay recognisable as the same files, and a structure
    carried out of either one should still say what produced it."""
    h, out = _run_through_delivery()
    try:
        bid = _filtered_ids(h)[0]
        preds = _delivered(out)[bid] / "predictions"
        names = sorted(p.name for p in preds.glob("*.pdb"))
        assert names[0] == f"seq000{s03_fold.PRED_SUFFIX}", names
        assert all("af2" in n for n in names)
        assert not any("rank_001" in n for n in names)
        # and the score file is found through the stage-03 helper, not by guesswork
        assert s03_fold.scores_path_for(preds / names[0]).exists()
    finally:
        h.close()


def test_delivery_says_so_when_a_score_file_is_missing(capsys=None):
    """Predictions written under the old stage-03 naming still copy, because the
    structure is found through pred_path, but their scores are looked up by name and
    quietly do not arrive. The run must say that rather than deliver confidence
    columns backed by nothing."""
    import contextlib
    import io

    h, _ = _run_through_stage04()
    try:
        bid = _filtered_ids(h)[0]
        preds = Path(h.cfg.outdir) / "03_predictions" / bid
        for scores in list(preds.glob("*_af2_scores.json"))[:3]:
            scores.unlink()

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            s04_delivery.run(h.cfg, h.state, force=True)
        out = buf.getvalue()
        assert "WARNING" in out and "3 predictions had no score file" in out, out

        # the structures still arrive; it is only the confidence that is absent
        folder = _delivered(Path(h.cfg.outdir) / "04_delivery")[bid] / "predictions"
        assert len(list(folder.glob("*.pdb"))) == N_SEQ
        assert len(list(folder.glob("*_af2_scores.json"))) == N_SEQ - 3
    finally:
        h.close()


def test_delivery_is_quiet_when_every_score_file_is_present():
    h, _ = _run_through_stage04()
    try:
        import contextlib
        import io
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            s04_delivery.run(h.cfg, h.state)
        assert "WARNING" not in buf.getvalue(), buf.getvalue()
    finally:
        h.close()


def test_delivery_sequences_json_records_every_sequence():
    """Every sequence is listed, exactly one is flagged best, and each points at the
    prediction it produced by a path relative to the folder, so the folder can be
    copied somewhere else and still resolve."""
    h, out = _run_through_delivery()
    try:
        seqs = h.state.sequences.load()
        import pandas as pd
        table = pd.read_csv(out / "delivery.csv")
        delivered = _delivered(out)
        for row in table.to_dict("records"):
            folder = delivered[row["backbone_id"]]
            doc = json.loads((folder / "sequences.json").read_text(encoding="utf-8"))

            assert doc["backbone_id"] == row["backbone_id"]
            assert doc["n_sequences"] == N_SEQ
            assert len(doc["sequences"]) == N_SEQ
            assert doc["best_sequence_id"] == row["best_sequence_id"]
            assert doc["criteria"]["rmsd_cutoff"] == 1.0
            assert doc["criteria"]["plddt_cutoff"] == 80.0

            flagged = [e for e in doc["sequences"] if e["best"]]
            assert len(flagged) == 1
            assert flagged[0]["sequence_id"] == row["best_sequence_id"]

            group = seqs[seqs.backbone_id == row["backbone_id"]]
            assert abs(row["best_rmsd"] - group["sc_rmsd"].min()) < 1e-9
            assert abs(flagged[0]["sc_rmsd"] - group["sc_rmsd"].min()) < 1e-9
            assert doc["n_pass"] == int(group["sc_pass"].fillna(False).astype(bool).sum())

            for entry in doc["sequences"]:
                assert (folder / entry["prediction"]).exists(), entry["prediction"]
                assert (folder / entry["scores"]).exists(), entry["scores"]
                assert entry["sequence"]
    finally:
        h.close()


def test_delivery_json_is_valid_json_when_values_are_missing():
    """json.dumps writes a bare NaN for a missing float, which is a Python extension
    that strict parsers reject. A backbone whose sequences never folded has several."""
    h, _ = _run_through_stage04()
    try:
        bid = _filtered_ids(h)[0]
        seqs = h.state.sequences.load()
        blanked = [{"sequence_id": s, "mean_plddt": None, "min_plddt": None,
                    "ptm": None, "mpnn_score": None}
                   for s in seqs[seqs.backbone_id == bid]["sequence_id"]]
        h.state.sequences.upsert(blanked)
        s04_delivery.run(h.cfg, h.state, force=True)

        out = Path(h.cfg.outdir) / "04_delivery"

        def strict(path):
            def boom(const):
                raise AssertionError(f"{path.name} carries a bare {const}")
            return json.loads(path.read_text(encoding="utf-8"), parse_constant=boom)

        doc = strict(out / "delivery.json")
        assert doc["n_backbones"]
        folder = _delivered(out)[bid]
        entry = strict(folder / "sequences.json")["sequences"][0]
        assert entry["mean_plddt"] is None
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
    h, out = _run_through_delivery()
    try:
        (out / "zz_stale.pdb").write_text("junk", encoding="utf-8")
        (out / s04_delivery.PASSED / "99_gone").mkdir()
        s04_delivery.run(h.cfg, h.state, force=True)
        assert not (out / "zz_stale.pdb").exists()
        assert not (out / s04_delivery.PASSED / "99_gone").exists()
    finally:
        h.close()


def test_a_dry_run_works_on_an_empty_run_directory(capsys=None):
    """The one moment --dry-run is actually wanted is before the first real run, when
    nothing upstream exists. Stages 02 and 03 used to raise there, so the flag failed
    precisely when it was needed, and the GPU commands were never shown at all."""
    import io
    import contextlib
    from cyclicnano import cli

    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_dryrun_"))
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = cli.main(["run", "--dry-run",
                             "--set", "target.symmetry=C5",
                             "--set", "target.monomer_length=60",
                             "--set", "generate.noise_scale=0.5",
                             "--set", "backend.kind=local",
                             "--set", f"run.outdir={tmp.as_posix()}"])
        out = buf.getvalue()
        assert code == 0, out

        # every tool command is shown, not just filed away in commands.jsonl
        assert "run_inference.py" in out, out
        assert "protein_mpnn_run.py" in out, out
        assert "colabfold" in out, out
        assert "denoiser.noise_scale_ca=0.5" in out, out

        # and nothing was actually produced
        assert not (tmp / ".done").exists()
        assert not (tmp / "logs" / "designs.fasta").exists()
        assert list((tmp / "00_backbones").glob("*.pdb")) == []
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_delivery_requires_a_stage_04_verdict():
    h = Harness()
    h.seed_backbones()
    try:
        s01_backbone_filter.run(h.cfg, h.state)
        try:
            s04_delivery.run(h.cfg, h.state)
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
