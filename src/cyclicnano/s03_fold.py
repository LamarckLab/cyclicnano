"""Stage 03 - AlphaFold2 monomer prediction via localcolabfold.

Per protocol only the monomer is predicted. MSA mode is forced to single_sequence:
a de novo design has no homologues, so a search either returns nothing or returns
unrelated hits, and in both cases the resulting confidence is not meaningful.

The stage leaves one folder per backbone holding two files per sequence,
seqNNN_af2.pdb and seqNNN_af2_scores.json, rather than the flat directory of long
names colabfold writes. The tool is named in the file rather than only in the stage
directory, so a structure stays identifiable once it is moved or sent on.
"""
from __future__ import annotations

import json
from pathlib import Path

from .config import Config
from .manifest import RunState
from .runner import Runner

STAGE = "03_predictions"

# The names this stage gives its two products. Defined here and imported by the
# delivery stage rather than spelled out twice: the two have drifted apart before,
# and a mismatch is silent, costing only the score file nobody notices is missing.
PRED_SUFFIX = "_af2.pdb"
SCORES_SUFFIX = "_af2_scores.json"


def scores_path_for(prediction: Path) -> Path:
    """The score file this stage writes beside a prediction."""
    name = prediction.name
    stem = name[:-len(PRED_SUFFIX)] if name.endswith(PRED_SUFFIX) else prediction.stem
    return prediction.parent / f"{stem}{SCORES_SUFFIX}"


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


def regroup(work: Path, rows: list[dict], backbone_of: dict) -> list[dict]:
    """Move each rank-1 prediction and its scores into a folder named for its backbone.

    colabfold names every file after the rank, the model and the seed as well as the
    sequence. Only rank 1 of a single model is ever produced here, so those 41
    characters are identical on all of them and distinguish nothing, while the flat
    directory reaches 640 entries on a 16-backbone run. One folder per backbone with
    the sequence index as the name says everything the long form did.

    Returns the rows with pred_path pointing at the new location.
    """
    out = []
    for row in rows:
        sid, pred = row["sequence_id"], row.get("pred_path")
        bid = backbone_of.get(sid)
        if not bid or not pred or not Path(str(pred)).exists():
            out.append(row)
            continue
        src = Path(str(pred))
        short = sid.split("__", 1)[1] if "__" in sid else sid
        dest = work / str(bid)
        dest.mkdir(parents=True, exist_ok=True)
        new_pdb = dest / f"{short}{PRED_SUFFIX}"
        src.replace(new_pdb)
        for scores in sorted(src.parent.glob(f"{sid}_scores_rank_001_*.json")):
            scores.replace(scores_path_for(new_pdb))
        out.append({**row, "pred_path": str(new_pdb)})
    return out


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

    # colabfold_batch writes eight files per sequence; two are read downstream. The
    # rest are plots, the single-sequence a3m and PAE dumps, and they go. The two
    # run-level files it leaves in the output directory, its resolved config and the
    # citation list, are provenance rather than results and move to logs/.
    if not cfg.get("fold.keep_all_outputs", False):
        rows = regroup(work, rows, {r["sequence_id"]: r.get("backbone_id") for r in records})
        archive = logs / "colabfold"
        removed = 0
        for path in work.iterdir():
            if not path.is_file():                          # the regrouped folders stay
                continue
            if path.name in ("config.json", "cite.bibtex"):
                archive.mkdir(parents=True, exist_ok=True)
                path.replace(archive / path.name)
            else:
                path.unlink()
                removed += 1
        print(f"[{STAGE}] removed {removed} unused colabfold products")
    state.sequences.upsert(rows)
    if not runner.dry_run:                                   # a dry run must leave no trace
        state.mark_done(STAGE)
    n_ok = sum(1 for r in rows if r.get("fold_ok"))
    print(f"[{STAGE}] {n_ok}/{len(rows)} predictions recovered")
    return rows
