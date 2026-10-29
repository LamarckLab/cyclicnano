"""Stage 04b - lay out the results one folder per backbone.

Every backbone that survived stage 01 gets a folder, whether or not it cleared the
stage-04 cutoffs. That is deliberate. The cutoffs are uncalibrated, the survivors of
a strict shape filter are few, and a backbone whose twenty sequences all land at
1.4 A is a different kind of near-miss from one whose best is 6 A -- discarding both
as "failed" throws away the distinction. The verdict is recorded as a column instead,
and the ranking puts the ones that passed on top.

Each folder holds the design backbone, every sequence written for it, and the
prediction and score file for each of those sequences. Nothing here is new: it is the
output of stages 02 and 03 regrouped by candidate, so that one can be looked at
without cross-referencing three directories.

The folders are split into passed/ and below_cutoff/, and the passing ones carry
their rank as a prefix, so a directory listing answers which cleared stage 04 and in
what order without opening the table. The cost is that a folder name encodes a
verdict: rerunning with a different cutoff renames things, and a path written down
elsewhere can go stale.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

from .config import Config
from .manifest import RunState
from .s00_generate import backbone_filename
from .s03_fold import scores_path_for

STAGE = "04_delivery"
PASSED = "passed"          # cleared the stage-04 cutoffs
BELOW = "below_cutoff"     # kept and flagged, not discarded


def _f(value, spec: str = ".2f", missing: str = "na") -> str:
    """Format a number that may be absent, without letting a NaN reach the header."""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return missing
    return missing if x != x else format(x, spec)


def best_sequence(sequences: pd.DataFrame, backbone_id: str) -> dict | None:
    """The sequence of this backbone that folded back closest to the design."""
    group = sequences[sequences["backbone_id"] == backbone_id]
    group = group[group["sc_rmsd"].notna()]
    return group.nsmallest(1, "sc_rmsd").iloc[0].to_dict() if len(group) else None


def _copy_predictions(group: list[dict], dest: Path) -> tuple[int, int]:
    """Copy each sequence prediction and its scores, keeping the stage-03 names.

    Stage 03 has already reduced these to seqNNN_af2.pdb and seqNNN_af2_scores.json
    inside a folder named for the backbone. The delivery folder is named for the backbone too, so the
    same short names carry here unchanged and the two layouts stay recognisable as
    the same files.
    """
    n = missing = 0
    for row in group:
        pred = row.get("pred_path")
        if not pred or not Path(str(pred)).exists():
            continue
        src = Path(str(pred))
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest / src.name)
        n += 1
        scores = scores_path_for(src)                       # named by stage 03, not guessed
        if scores.exists():
            shutil.copy2(scores, dest / scores.name)
        else:
            missing += 1
    return n, missing


def _write_fasta(group: list[dict], best_id: str | None, path: Path) -> None:
    lines = []
    for row in sorted(group, key=lambda r: r.get("seq_index", 0)):
        tag = " BEST" if row["sequence_id"] == best_id else ""
        lines.append(
            f">{row['sequence_id']} rmsd={_f(row.get('sc_rmsd'))} "
            f"plddt={_f(row.get('mean_plddt'), '.1f')} "
            f"mpnn_score={_f(row.get('mpnn_score'), '.3f')}{tag}\n{row['sequence']}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


_SEQ_COLUMNS = ["sequence_id", "seq_index", "sequence", "mpnn_score", "mpnn_global_score",
                "mean_plddt", "min_plddt", "ptm", "sc_rmsd", "sc_pass"]


def run(cfg: Config, state: RunState, force: bool = False, dry_run: bool = False) -> dict:
    if dry_run:                                             # no subprocess to simulate
        print(f"[{STAGE}] dry run, skipping")
        return {}
    if state.is_done(STAGE) and not force:
        print(f"[{STAGE}] already done, skipping")
        return {}

    backbones = state.backbones.load()
    sequences = state.sequences.load()
    if backbones.empty:
        raise RuntimeError("no backbones in the manifest; run the earlier stages first")
    if "final_pass" not in backbones.columns:
        raise RuntimeError("no stage 04 verdict in the manifest; run stage 04 first")

    outdir = state.stage_dir(STAGE)
    for stale in outdir.iterdir():                          # a rerun must not leave old hits behind
        stale.unlink() if stale.is_file() else shutil.rmtree(stale, ignore_errors=True)

    # Everything the shape filter kept, not only what the RMSD cutoff kept.
    delivered = backbones[backbones["bb_pass"].fillna(False).astype(bool)]
    delivered = delivered.sort_values(
        ["final_pass", "success_rate", "best_rmsd"], ascending=[False, False, True])

    extra = [c for c in backbones.columns if c.startswith("success_rate_at_")]
    rows, rank, missing_scores = [], 0, 0
    for record in delivered.to_dict("records"):
        bid = record["backbone_id"]
        # Grouped and numbered so that the verdict and the ranking are visible from a
        # directory listing. The near misses keep their own folder rather than being
        # dropped; what they lose is only the top of the list.
        if bool(record.get("final_pass")):
            rank += 1
            rel = Path(PASSED) / f"{rank}_{bid}"
        else:
            rel = Path(BELOW) / str(bid)
        folder = outdir / rel
        folder.mkdir(parents=True, exist_ok=True)

        source = Path(str(record["backbone_path"]))
        if source.exists():
            shutil.copy2(source, folder / backbone_filename(bid))

        group = sequences[sequences["backbone_id"] == bid].to_dict("records")
        best = best_sequence(sequences, bid)
        best_id = best["sequence_id"] if best else None

        n_pred, n_missing = _copy_predictions(group, folder / "predictions")
        missing_scores += n_missing
        if group:
            _write_fasta(group, best_id, folder / "sequences.fasta")
            cols = [c for c in _SEQ_COLUMNS if c in sequences.columns]
            per_seq = pd.DataFrame(group)[cols].sort_values("seq_index")
            per_seq.to_csv(folder / "summary.csv", index=False)

        rows.append({
            "backbone_id": bid,
            "folder": rel.as_posix(),
            "symmetry": record.get("symmetry"),
            "n_res_monomer": record.get("n_res_monomer"),
            "n_helices": record.get("n_helices"),
            "strand_frac": record.get("strand_frac"),
            "max_centre_distance": record.get("max_centre_distance"),
            "n_seqs": record.get("n_seqs"),
            "n_predictions": n_pred,
            "success_rate": record.get("success_rate"),
            **{c: record.get(c) for c in extra},
            "best_rmsd": best["sc_rmsd"] if best else None,
            # From the same row as best_rmsd and best_sequence_id, so the three cannot
            # describe different sequences.
            "best_plddt": best.get("mean_plddt") if best else None,
            "final_pass": bool(record.get("final_pass")),
            "final_fail_reason": record.get("final_fail_reason"),
            "best_sequence_id": best_id,
            "best_sequence": best["sequence"] if best else None,
        })

    summary = outdir / "delivery.csv"
    pd.DataFrame(rows).to_csv(summary, index=False)
    # The same table as JSON, for anything reading the run without pandas.
    (outdir / "delivery.json").write_text(
        json.dumps({"n_backbones": len(rows),
                    "n_final_pass": sum(1 for r in rows if r["final_pass"]),
                    "backbones": rows}, indent=2, default=str),
        encoding="utf-8")

    # A structure whose scores file is absent is delivered without its confidence, and
    # nothing downstream notices. It means the predictions were written by a different
    # version of stage 03 than the one naming them now, which is a migration left half
    # done rather than a missing file, so it is said out loud rather than raised.
    if missing_scores:
        print(f"[{STAGE}] WARNING: {missing_scores} predictions had no score file beside "
              f"them. The stage-03 output predates the current naming; re-run stage 03 "
              f"or migrate it before trusting the confidence columns.")

    n_pass = sum(1 for r in rows if r["final_pass"])
    print(f"[{STAGE}] {len(rows)} backbones delivered to {outdir} "
          f"({n_pass} of them pass stage 04)")
    if rows:
        top = rows[0]
        print(f"[{STAGE}] best: {top['backbone_id']} "
              f"(success_rate {_f(top['success_rate'])}, rmsd {_f(top['best_rmsd'])} A)")
    state.mark_done(STAGE)
    return {"n_delivered": len(rows), "n_final_pass": n_pass, "summary": str(summary)}
