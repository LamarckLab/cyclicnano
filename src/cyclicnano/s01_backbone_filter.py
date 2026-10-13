"""Stage 01 - objective backbone filtering.

This is the stage that replaces "eyeball the backbones and throw out the bad ones".
It computes every metric for every backbone unconditionally, then applies the rules
from the config to set a pass flag. Nothing is deleted, so raising or lowering a
threshold costs one CPU-minute instead of a GPU-day.

The stage leaves two things behind: `filter_summary.json`, which states the rules and
how many backbones each one rejected, and `passed/`, holding a copy of the structures
that survived. The copy exists so that the selection is inspectable on its own -- the
pass flag also lives in the backbone table, but a folder of PDBs is what gets opened
in PyMOL -- and it is what stage 02 hands to ProteinMPNN.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .config import Config
from .filters import FilterSet
from .geometry import backbone_metrics
from .manifest import RunState
from .pdbio import read_pdb

STAGE = "01_backbone_filter"
PASSED = "passed"


def passed_dir(state: RunState) -> Path:
    """Where the surviving backbones are copied. Stage 02 reads this directory."""
    return state.outdir / STAGE / PASSED


def metrics_for(path: str | Path, expected_sym: int) -> dict:
    return backbone_metrics(read_pdb(str(path)), expected_sym=expected_sym)


def run(cfg: Config, state: RunState, force: bool = False, dry_run: bool = False) -> dict:
    # The GPU stages take a Runner and honour its dry_run; this one has no subprocess
    # to simulate, so it has to be told. Without it a dry run rewrites the folder that
    # stage 02 reads, which is a real change made by a command that promises none.
    if dry_run:
        print(f"[{STAGE}] dry run, skipping")
        return {}
    if state.is_done(STAGE) and not force:
        print(f"[{STAGE}] already done, skipping")
        return {}

    df = state.backbones.load()
    if df.empty:
        raise RuntimeError("no backbones in the manifest; run stage 00 first")

    rules = FilterSet("backbone", cfg.get("backbone_filter.rules", []) or [])
    expected = cfg.sym_order

    rows, all_metrics, survivors = [], [], []
    for record in df.to_dict("records"):
        bid, path = record["backbone_id"], record["backbone_path"]
        try:
            m = metrics_for(path, expected)
        except Exception as exc:                            # a malformed PDB is a fail, not a crash
            rows.append({"backbone_id": bid, "bb_pass": False,
                         "bb_fail_reason": f"metric error: {exc}"})
            continue

        verdict = rules.apply(m)
        all_metrics.append(m)
        row = {"backbone_id": bid, "bb_pass": verdict.passed, "bb_fail_reason": verdict.reason}
        row.update({k: v for k, v in m.items() if k != "ss_string"})
        row["ss_string"] = m["ss_string"]
        rows.append(row)
        if verdict.passed:
            survivors.append((bid, path))

    state.backbones.upsert(rows)

    # Rebuilt from scratch every time, because a tightened rule must not leave the
    # backbones it now rejects sitting in the folder that stage 02 reads.
    staged = passed_dir(state)
    if staged.exists():
        shutil.rmtree(staged)
    staged.mkdir(parents=True)
    for bid, path in survivors:
        shutil.copy2(path, staged / f"{bid}.pdb")

    summary = rules.summary(all_metrics) if all_metrics else {"rules": list(rules.rules),
                                                              "n_input": 0, "n_pass": 0,
                                                              "rejected_by_rule": {}}
    summary["n_metric_errors"] = sum(1 for r in rows if not r.get("bb_pass")
                                     and str(r.get("bb_fail_reason", "")).startswith("metric error"))
    summary["passed_backbones"] = [bid for bid, _ in survivors]
    out = state.stage_dir(STAGE) / "filter_summary.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[{STAGE}] {summary['n_pass']}/{summary['n_input']} backbones pass")
    for rule, n in sorted(summary["rejected_by_rule"].items(), key=lambda kv: -kv[1]):
        if n:
            print(f"          rejected by {rule!r}: {n}")
    print(f"[{STAGE}] survivors copied to {staged}")
    state.mark_done(STAGE)
    return summary
