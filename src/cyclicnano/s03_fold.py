"""Stage 03 - AlphaFold2 monomer prediction via localcolabfold.

Per protocol only the monomer is predicted. MSA mode is forced to single_sequence:
a de novo design has no homologues, so a search either returns nothing or returns
unrelated hits, and in both cases the resulting confidence is not meaningful.
"""
from __future__ import annotations

import json
from pathlib import Path

from .config import Config
from .manifest import RunState
from .runner import Runner

STAGE = "03_predictions"


def write_fasta(rows: list[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(f">{r['sequence_id']}\n{r['sequence']}\n")
    return path


def build_command(cfg: Config, fasta: Path, outdir: Path) -> list[str]:
    """Mirror the verified colabfold_batch invocation: --data, flags, input, output.

    `--data` points at the weight directory and is not optional: without it
    colabfold_batch looks in its own default location and will try to download the
    parameters again.
    """
    exe = cfg.get("paths.colabfold") or "colabfold_batch"
    argv = [str(exe)]

    params = cfg.get("paths.colabfold_params")
    if params:
        argv += ["--data", str(params)]

    argv += [
        "--msa-mode", str(cfg.get("fold.msa_mode", "single_sequence")),
        "--num-recycle", str(cfg.get("fold.num_recycle", 3)),
        "--num-models", str(cfg.get("fold.num_models", 1)),
    ]
    # Opt-in only: the verified invocation sets neither, and colabfold_batch already
    # defaults to alphafold2_ptm for monomers.
    if cfg.get("fold.model_type"):
        argv += ["--model-type", str(cfg.get("fold.model_type"))]
    if cfg.get("fold.random_seed") is not None:
        argv += ["--random-seed", str(int(cfg.get("fold.random_seed")))]
    if cfg.get("fold.use_templates", False):
        argv.append("--templates")

    return argv + [str(fasta), str(outdir)]                 # positionals last, as verified


def collect(outdir: Path, sequence_id: str) -> dict:
    """Pick the rank-1 prediction and its score file for one sequence."""
    pdbs = sorted(outdir.glob(f"{sequence_id}_unrelaxed_rank_001_*.pdb"))
    if not pdbs:
        pdbs = sorted(outdir.glob(f"{sequence_id}_*rank_001*.pdb"))
    if not pdbs:
        return {"pred_path": None, "fold_ok": False}

    row = {"pred_path": str(pdbs[0]), "fold_ok": True}
    scores = sorted(outdir.glob(f"{sequence_id}_scores_rank_001_*.json"))
    if scores:
        data = json.loads(scores[0].read_text(encoding="utf-8"))
        plddt = data.get("plddt") or []
        if plddt:
            row["mean_plddt"] = float(sum(plddt) / len(plddt))
            row["min_plddt"] = float(min(plddt))
        if data.get("ptm") is not None:
            row["ptm"] = float(data["ptm"])
    return row


def run(cfg: Config, state: RunState, runner: Runner, force: bool = False) -> list[dict]:
    if state.is_done(STAGE) and not force:
        print(f"[{STAGE}] already done, skipping")
        return state.sequences.load().to_dict("records")

    if cfg.get("fold.msa_mode") != "single_sequence":
        print(f"[{STAGE}] WARNING: msa_mode={cfg.get('fold.msa_mode')!r} on de novo designs")

    df = state.sequences.load()
    work = state.stage_dir(STAGE)
    logs = state.outdir / "logs"

    # As in stage 02: the command needs the FASTA path, not the FASTA, so it can be
    # shown before anything upstream has run. Writing the file would leave a trace.
    if runner.dry_run:
        runner.run("colabfold", build_command(cfg, logs / "designs.fasta", work))
        return []

    if df.empty:
        raise RuntimeError("no sequences in the manifest; run stage 02 first")

    records = df.to_dict("records")
    fasta = write_fasta(records, logs / "designs.fasta")     # an input, not a result
    runner.run("colabfold", build_command(cfg, fasta, work))

    rows = [{"sequence_id": r["sequence_id"], **collect(work, r["sequence_id"])} for r in records]
    state.sequences.upsert(rows)

    # colabfold_batch writes eight files per sequence; two are read downstream. The
    # rest are plots, the single-sequence a3m and PAE dumps, and they go. The two
    # run-level files it leaves in the output directory, its resolved config and the
    # citation list, are provenance rather than results and move to logs/, so the
    # stage directory ends up holding one structure and one score file per sequence
    # and nothing else.
    if not cfg.get("fold.keep_all_outputs", False):
        keep = {Path(r["pred_path"]).name for r in rows if r.get("pred_path")}
        keep |= {p.name for p in work.glob("*_scores_rank_001_*.json")}
        archive = logs / "colabfold"
        removed = 0
        for path in work.iterdir():
            if not path.is_file() or path.name in keep:
                continue
            if path.name in ("config.json", "cite.bibtex"):
                archive.mkdir(parents=True, exist_ok=True)
                path.replace(archive / path.name)
            else:
                path.unlink()
                removed += 1
        print(f"[{STAGE}] removed {removed} unused colabfold products")
    if not runner.dry_run:                                   # a dry run must leave no trace
        state.mark_done(STAGE)
    n_ok = sum(1 for r in rows if r.get("fold_ok"))
    print(f"[{STAGE}] {n_ok}/{len(rows)} predictions recovered")
    return rows
