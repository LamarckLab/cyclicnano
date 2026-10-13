"""Stage 04 - self-consistency scoring and per-backbone rollup.

Per sequence: CA RMSD between the predicted monomer and the chain-A monomer of the
design backbone. Per backbone: the fraction of its sequences that fold back within
the cutoff. That success rate is the ranking quantity the protocol ends on.

Extra cutoffs from `validate.also_record` are stored alongside the primary one, so a
stricter or looser definition can be applied later without re-predicting anything.

The stage writes no directory of its own. Everything it produces is a column on the
backbone or sequence table, and the ranked view of those columns is what the delivery
stage lays out; a second ranking file here would only be the same numbers under a
different name, free to disagree with the first.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import Config
from .filters import FilterSet
from .geometry import rmsd
from .manifest import RunState
from .pdbio import read_pdb

STAGE = "04_validate"


def monomer_rmsd(design_pdb: str | Path, predicted_pdb: str | Path) -> float:
    """CA RMSD of the predicted monomer against chain A of the design oligomer."""
    design = read_pdb(str(design_pdb))
    pred = read_pdb(str(predicted_pdb))
    a = design.chains[design.chain_ids[0]].ca
    b = pred.chains[pred.chain_ids[0]].ca
    if len(a) != len(b):
        n = min(len(a), len(b))                             # tolerate a trimmed terminus
        a, b = a[:n], b[:n]
    return rmsd(a, b)


def run(cfg: Config, state: RunState, force: bool = False, dry_run: bool = False) -> dict:
    if dry_run:                                             # no subprocess to simulate
        print(f"[{STAGE}] dry run, skipping")
        return {}
    if state.is_done(STAGE) and not force:
        print(f"[{STAGE}] already done, skipping")
        return {}

    seqs = state.sequences.load()
    bbs = state.backbones.load()
    if seqs.empty:
        raise RuntimeError("no sequences in the manifest; run stages 02-03 first")

    bb_path = dict(zip(bbs["backbone_id"], bbs["backbone_path"]))
    cutoff = float(cfg.get("validate.rmsd_cutoff", 1.0))
    extra = [float(x) for x in (cfg.get("validate.also_record", []) or [])]
    plddt_cut = float(cfg.get("validate.plddt_cutoff", 0.0))

    seq_rows = []
    for r in seqs.to_dict("records"):
        sid, bid, pred = r["sequence_id"], r["backbone_id"], r.get("pred_path")
        if not pred or not Path(str(pred)).exists():
            seq_rows.append({"sequence_id": sid, "sc_rmsd": None, "sc_pass": False})
            continue
        try:
            value = monomer_rmsd(bb_path[bid], pred)
        except Exception as exc:
            seq_rows.append({"sequence_id": sid, "sc_rmsd": None, "sc_pass": False,
                             "sc_error": str(exc)})
            continue
        plddt = r.get("mean_plddt")
        row = {
            "sequence_id": sid,
            "sc_rmsd": value,
            "sc_pass": bool(value <= cutoff and (plddt is None or float(plddt) >= plddt_cut)),
        }
        for c in extra:
            row[f"sc_pass_at_{c:g}A"] = bool(value <= c)
        seq_rows.append(row)
    state.sequences.upsert(seq_rows)

    seqs = state.sequences.load()
    rules = FilterSet("validate", cfg.get("validate.rules", []) or [])
    bb_rows, metrics = [], []
    for bid, group in seqs.groupby("backbone_id"):
        rm = group["sc_rmsd"].astype(float)
        valid = rm.dropna()
        m = {
            "backbone_id": bid,
            "n_seqs": int(len(group)),
            "n_seqs_folded": int(len(valid)),
            "success_rate": float(group["sc_pass"].fillna(False).astype(bool).mean()),
            "best_rmsd": float(valid.min()) if len(valid) else float("inf"),
            "median_rmsd": float(valid.median()) if len(valid) else float("inf"),
            "best_plddt": float(group["mean_plddt"].max()) if "mean_plddt" in group else float("nan"),
        }
        for c in extra:
            col = f"sc_pass_at_{c:g}A"
            if col in group:
                m[f"success_rate_at_{c:g}A"] = float(group[col].fillna(False).astype(bool).mean())
        verdict = rules.apply(m)
        metrics.append(m)
        bb_rows.append({**m, "final_pass": verdict.passed, "final_fail_reason": verdict.reason})
    state.backbones.upsert(bb_rows)

    summary = {
        "n_backbones": len(bb_rows),
        "n_final_pass": sum(1 for r in bb_rows if r["final_pass"]),
        "mean_success_rate": float(np.mean([r["success_rate"] for r in bb_rows])) if bb_rows else 0.0,
        "rule_summary": rules.summary(metrics) if metrics else {},
    }
    print(f"[{STAGE}] {summary['n_final_pass']}/{summary['n_backbones']} backbones pass; "
          f"mean success rate {summary['mean_success_rate']:.2f}")
    state.mark_done(STAGE)
    return summary
