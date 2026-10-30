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

import numpy as np
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


def _jsonable(value):
    """NaN and pandas NA are not JSON. json.dumps writes a bare NaN token for them,
    which is a Python extension that strict parsers reject, so they become null."""
    if value is None:
        return None
    if isinstance(value, float) and value != value:                 # NaN
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if float(value) != float(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


_SEQ_FIELDS = ["seq_index", "sc_pass", "sc_rmsd", "mean_plddt", "min_plddt", "ptm",
               "mpnn_score", "mpnn_global_score"]

# Written by stage 05 when it runs. A sequence that was never evaluated carries null
# rather than a missing key, so that "not measured" is visible and cannot be mistaken
# for "measured and poor".
_INTERFACE_FIELDS = ["iptm", "multimer_ptm", "multimer_plddt", "dG_interface",
                     "dSASA_interface", "dG_dSASA_100"]
INTERFACE_DIR = "interface"


def _interface_block(row: dict, folder: Path) -> dict | None:
    """The stage-05 measurements for one sequence, or None if it was not evaluated.

    These describe the contact between subunits, which nothing in stages 03 and 04
    looks at: RMSD and pLDDT are both properties of a single chain. They are reported
    and never filtered on, so a poor interface marks a design rather than removing it.
    """
    # Through _jsonable first. A sequence that was never evaluated reads back from the
    # table as NaN, and bool(NaN) is True in Python, so the obvious test passes for
    # exactly the rows it was written to catch and every unevaluated sequence was
    # handed a block of nulls instead of a null block.
    if not _jsonable(row.get("interface_scored")) and _jsonable(row.get("iptm")) is None:
        return None
    block = {k: _jsonable(row.get(k)) for k in _INTERFACE_FIELDS}
    pred = row.get("interface_pred_path")
    if pred:
        name = Path(str(pred)).name
        block["prediction"] = (f"{INTERFACE_DIR}/{name}"
                               if (folder / INTERFACE_DIR / name).exists() else None)
    else:
        block["prediction"] = None
    return block


def _copy_interface(group: list[dict], dest: Path) -> int:
    """Carry the predicted oligomer into the folder beside the monomer predictions."""
    n = 0
    for row in group:
        pred = row.get("interface_pred_path")
        if not pred or not Path(str(pred)).exists():
            continue
        src = Path(str(pred))
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest / src.name)
        n += 1
    return n


def _write_sequences(group: list[dict], best_id: str | None, backbone_id: str,
                     folder: Path, criteria: dict) -> None:
    """Every sequence written for this backbone, with where its prediction landed.

    This replaces the FASTA the folder used to carry. A FASTA can hold the sequence
    and whatever is crammed into its header, and that header had grown into an ad-hoc
    record nothing could parse. The sequences are still a FASTA one directory up, in
    02_sequences, for any tool that wants one.
    """
    entries = []
    for row in sorted(group, key=lambda r: r.get("seq_index", 0)):
        raw = row.get("pred_path")
        pred = Path(str(raw)) if raw else None
        entry = {"sequence_id": row["sequence_id"],
                 "best": row["sequence_id"] == best_id}
        entry.update({k: _jsonable(row.get(k)) for k in _SEQ_FIELDS})
        # Relative to this folder, so the folder survives being copied elsewhere.
        has_pred = bool(pred and pred.exists())
        entry["prediction"] = f"predictions/{pred.name}" if has_pred else None
        scores = scores_path_for(pred) if has_pred else None
        entry["scores"] = f"predictions/{scores.name}" if scores and scores.exists() else None
        entry["interface"] = _interface_block(row, folder)
        entry["sequence"] = row["sequence"]                         # last: it is the long one
        entries.append(entry)

    passed = sum(1 for e in entries if e["sc_pass"])
    (folder / "sequences.json").write_text(
        json.dumps({
            "backbone_id": backbone_id,
            "n_sequences": len(entries),
            "n_pass": passed,
            "success_rate": round(passed / len(entries), 4) if entries else None,
            "best_sequence_id": best_id,
            "criteria": criteria,
            "sequences": entries,
        }, indent=2),
        encoding="utf-8")


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

    # Repeated into every sequences.json so a folder states the rule it was judged by.
    criteria = {"rmsd_cutoff": float(cfg.get("validate.rmsd_cutoff", 1.0)),
                "plddt_cutoff": float(cfg.get("validate.plddt_cutoff", 0.0))}
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
        n_iface = _copy_interface(group, folder / INTERFACE_DIR)
        if group:
            _write_sequences(group, best_id, bid, folder, criteria)
            cols = [c for c in _SEQ_COLUMNS if c in sequences.columns]
            per_seq = pd.DataFrame(group)[cols].sort_values("seq_index")
            per_seq.to_csv(folder / "summary.csv", index=False)

        # Deliberately narrow. Everything the shape filter measured, and the success
        # rate at the looser cutoffs, stays in tables/backbones.csv; repeating it here
        # made a table nobody could read across. The protocol settles on 1.0 A, so the
        # looser readings are not what this table is answering. What is left is which
        # design, where it is, how well it folded back, and the verdict.
        rows.append({
            "backbone_id": bid,
            "folder": rel.as_posix(),                       # scripts resolve the path from this
            "symmetry": record.get("symmetry"),
            "n_res_monomer": record.get("n_res_monomer"),
            "n_seqs": record.get("n_seqs"),
            "n_predictions": n_pred,
            "n_interface_evaluated": n_iface,
            "success_rate": record.get("success_rate"),
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
                    "backbones": [{k: _jsonable(v) for k, v in r.items()}
                                  for r in rows]}, indent=2, default=str),
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
