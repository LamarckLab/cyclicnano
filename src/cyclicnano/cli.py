"""Command-line entry point.

    cyclicnano run       --config my.yaml --profile config/profiles/amax.yaml
    cyclicnano metrics   some_backbone.pdb --symmetry C5
    cyclicnano check     --config my.yaml           # validate config and rule syntax only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import (s00_generate, s01_backbone_filter, s02_design, s03_fold, s04_validate,
               s05_delivery)
from .config import load_config, parse_overrides
from .filters import FilterSet, RuleError
from .geometry import backbone_metrics
from .manifest import RunState
from .pdbio import read_pdb
from .runner import runner_from_config

STAGES = ["00", "01", "02", "03", "04", "05"]
_SAMPLE_METRICS = {
    "n_chains": 5, "n_res_monomer": 60, "helix_frac": 0.5, "strand_frac": 0.1,
    "loop_frac": 0.4, "n_sse": 3, "loop_max_len": 8, "nterm_helix_len": 5,
    "cterm_helix_len": 5, "rg": 12.0, "rg_ratio": 1.0, "contact_order": 0.2,
    "n_helices": 3, "asphericity": 0.4, "acylindricity": 0.1,
    "shape_anisotropy": 0.2, "axis_ratio": 1.8,
    "assembly_asphericity": 0.1, "assembly_acylindricity": 0.05,
    "assembly_shape_anisotropy": 0.08, "assembly_axis_ratio": 1.4,
    "iface_contacts_per_chain": 40.0, "pore_radius": 8.0, "max_radius": 20.0,
    "height": 30.0, "sym_angle_deg": 72.0, "sym_order_detected": 5.0, "sym_rmsd": 0.1,
    "n_clash_intra": 0, "n_clash_inter": 0, "n_clash": 0, "sym_order_ok": True,
    "min_interchain_dist": 3.0,
    "max_axis_angle": 25.0, "mean_axis_angle": 15.0,
    "max_inter_helix_angle": 20.0, "max_helix_len": 18, "max_strand_len": 6,
    "success_rate": 0.5, "best_rmsd": 0.8, "median_rmsd": 1.2, "n_seqs": 10,
    "n_seqs_folded": 10, "best_plddt": 90.0,
}


def _check(cfg) -> int:
    problems = cfg.validate()
    for group, key in (("backbone", "backbone_filter.rules"), ("validate", "validate.rules")):
        try:
            FilterSet(group, cfg.get(key, []) or []).check_syntax(_SAMPLE_METRICS)
        except RuleError as exc:
            problems.append(f"{key}: {exc}")
    if problems:
        print("config problems:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    print(f"config OK: {cfg.symmetry} x {cfg.monomer_length} aa "
          f"(RFdiffusion contig total = {cfg.total_length})")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="cyclicnano", description="de novo cyclic oligomer design pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)

    for name in ("run", "check"):
        p = sub.add_parser(name)
        p.add_argument("--config", action="append", dest="config",
                       help="may be given several times; later files win")
        p.add_argument("--profile")
        p.add_argument("--set", dest="overrides", action="append", default=[],
                       help="dotted override, e.g. --set target.symmetry=C3")
        if name == "run":
            p.add_argument("--stages", default=",".join(STAGES))
            p.add_argument("--dry-run", action="store_true")
            p.add_argument("--force", action="store_true")

    pm = sub.add_parser("metrics")
    pm.add_argument("pdb")
    pm.add_argument("--symmetry", default=None, help="expected Cn, e.g. C5")

    args = ap.parse_args(argv)

    if args.cmd == "metrics":
        expected = int(args.symmetry[1:]) if args.symmetry else None
        m = backbone_metrics(read_pdb(args.pdb), expected_sym=expected)
        print(json.dumps(m, indent=2, default=str))
        return 0

    cfg = load_config(args.config, args.profile, parse_overrides(args.overrides))
    if args.cmd == "check":
        return _check(cfg)

    if _check(cfg) != 0:
        return 1

    state = RunState(cfg.outdir)
    Path(cfg.outdir).mkdir(parents=True, exist_ok=True)

    # Record the recipe beside the results: a run made with --set overrides is
    # otherwise only reconstructable from shell history.
    if not args.dry_run:
        logs = Path(cfg.outdir) / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        (logs / "resolved_config.yaml").write_text(cfg.resolved_yaml(), encoding="utf-8")
    runner = runner_from_config(cfg, cfg.outdir, dry_run=args.dry_run)
    wanted = [s.strip() for s in args.stages.split(",") if s.strip()]

    if "00" in wanted:
        s00_generate.run(cfg, state, runner, force=args.force)
    if "01" in wanted:
        s01_backbone_filter.run(cfg, state, force=args.force)
    if "02" in wanted:
        s02_design.run(cfg, state, runner, force=args.force)
    if "03" in wanted:
        s03_fold.run(cfg, state, runner, force=args.force)
    if "04" in wanted:
        s04_validate.run(cfg, state, force=args.force)
    if "05" in wanted:
        s05_delivery.run(cfg, state, force=args.force)

    print(f"\ntables: {state.backbones.path}\n        {state.sequences.path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
