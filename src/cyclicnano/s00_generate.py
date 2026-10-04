"""Stage 00 - symmetric backbone generation with RFdiffusion.

The emitted command mirrors a hand-run invocation that is known to work, and adds
nothing to it by default. Extra hydra keys are opt-in, because hydra rejects keys the
installed config does not define and an unsupported one fails the whole run.

The one conversion this stage performs is the contig length: the symmetric config
takes the TOTAL residue count across all chains, so a C5 of 60-residue subunits is
[300-300], not [60-60]. Passing the per-chain number silently produces a much smaller
oligomer, so `Config.total_length` computes it instead.
"""
from __future__ import annotations

import shutil
from pathlib import Path

from .config import Config
from .manifest import RunState
from .runner import Runner

STAGE = "00_backbones"


def _num(value) -> str:
    """Format 1.0 as 1 and 0.1 as 0.1, matching how the flags are written by hand."""
    f = float(value)
    return str(int(f)) if f == int(f) else str(f)


def build_command(cfg: Config, out_prefix: Path, n_designs: int,
                  hydra_dir: Path | None = None, seed: int | None = None,
                  start: int | None = None) -> list[str]:
    root = cfg.get("paths.rfdiffusion")
    script = f"{root}/scripts/run_inference.py" if root else "run_inference.py"

    argv = [
        "python", script,
        "--config-name", "symmetry",
        f"inference.symmetry={cfg.symmetry.lower()}",       # RFdiffusion expects c5, not C5
        f"contigmap.contigs=[{cfg.total_length}-{cfg.total_length}]",
    ]

    if cfg.get("generate.potentials.enabled", True):
        pot = cfg.get("generate.potentials.olig_contacts", {}) or {}
        terms = [
            "type:olig_contacts",
            f"weight_intra:{_num(pot.get('weight_intra', 1))}",
            f"weight_inter:{_num(pot.get('weight_inter', 0.1))}",
        ]
        # r_0 / d_0 exist on the potential but are left at the RFdiffusion defaults
        # unless set, so the emitted command stays identical to the verified one.
        for name in ("r_0", "d_0"):
            value = cfg.get(f"generate.potentials.{name}")
            if value is not None:
                terms.append(f"{name}:{_num(value)}")
        argv += [
            'potentials.guiding_potentials=["' + ",".join(terms) + '"]',
            "potentials.olig_intra_all=True",
            "potentials.olig_inter_all=True",
            f"potentials.guide_scale={_num(cfg.get('generate.potentials.guide_scale', 2))}",
            f"potentials.guide_decay={cfg.get('generate.potentials.guide_decay', 'quadratic')}",
        ]

    # Trajectory files are the bulk of the output and nothing downstream reads them:
    # 50 C5 backbones produced 337 MB of them against 4 MB of actual structures.
    if not cfg.get("generate.write_trajectory", False):
        argv += ["inference.write_trajectory=False"]

    # Without this RFdiffusion seeds itself from entropy and the same command gives a
    # different hundred backbones every time. It seeds per design from the design
    # index, which the batch loop already advances, so every design in a run gets a
    # distinct seed and two runs of the same config produce the same structures.
    if cfg.get("generate.deterministic", True):
        argv += ["inference.deterministic=True"]

    # Batches share one prefix and continue each other's numbering, so the stage ends
    # with design_0 .. design_N-1 rather than a restart inside every batch. The start
    # index is counted by the stage rather than left to RFdiffusion's -1 autodetect,
    # so the caller stays in control of it. It also advances the per-design seed
    # across batches, which matters if inference.deterministic is ever turned on.
    override = cfg.get("generate.design_startnum")
    start = override if override is not None else start
    if start is not None:
        argv += [f"inference.design_startnum={int(start)}"]

    argv += [f"inference.output_prefix={out_prefix}"]
    if hydra_dir is not None:
        argv += [f"hydra.run.dir={hydra_dir}"]              # keep hydra logs out of the cwd
    argv += [f"inference.num_designs={n_designs}"]

    # Opt-in only: the seed key differs between RFdiffusion builds and an unknown
    # hydra key aborts the run, so it is emitted solely when configured explicitly.
    seed_key = cfg.get("generate.seed_key")
    if seed_key and seed is not None:
        argv += [f"{seed_key}={seed}"]

    return argv + list(cfg.get("generate.extra_args", []) or [])


def run(cfg: Config, state: RunState, runner: Runner, force: bool = False) -> list[dict]:
    if state.is_done(STAGE) and not force:
        print(f"[{STAGE}] already done, skipping")
        return state.backbones.load().to_dict("records")

    outdir = state.stage_dir(STAGE)
    configured = cfg.get("generate.hydra_run_dir")
    # Hydra logs belong with the other logs, not among the structures.
    hydra_dir = Path(configured) if configured else state.outdir / "logs" / "hydra"

    total = int(cfg.require("target.n_backbones"))
    batch = int(cfg.get("generate.batch_size", 10))
    base_seed = int(cfg.get("run.seed", 0))

    prefix = outdir / str(cfg.get("generate.output_prefix", "design"))
    if force:
        for stale in outdir.glob(f"{prefix.name}_*.pdb"):   # otherwise auto-numbering appends
            stale.unlink()

    # Batching exists for checkpointing, not for naming: each call continues the
    # numbering of the last, so an interrupted run resumes at the next free index.
    b = 0
    while True:
        done = len(list(outdir.glob(f"{prefix.name}_*.pdb")))
        n = min(batch, total - done)
        if n <= 0:
            break
        runner.run("rfdiffusion",
                   build_command(cfg, prefix, n, hydra_dir / f"batch{b:03d}",
                                 base_seed + b, start=done))
        if runner.dry_run:
            break                                           # nothing lands, so the count never moves
        if len(list(outdir.glob(f"{prefix.name}_*.pdb"))) <= done:
            raise RuntimeError(f"RFdiffusion produced no new backbones in {outdir}")
        b += 1

    # Leave only structures behind: .trb holds RFdiffusion metadata that nothing
    # downstream reads, and a stray traj/ survives from runs made before the flag.
    if not cfg.get("generate.keep_trb", False):
        for junk in outdir.glob("*.trb"):
            junk.unlink()
    traj = outdir / "traj"
    if traj.is_dir() and not cfg.get("generate.write_trajectory", False):
        shutil.rmtree(traj, ignore_errors=True)

    rows = []
    for path in sorted(outdir.glob("*.pdb")):
        rows.append({
            "backbone_id": path.stem,
            "backbone_path": str(path),
            "symmetry": cfg.symmetry,
            "monomer_length_requested": cfg.monomer_length,
            "stage": STAGE,
        })
    state.backbones.upsert(rows)
    if not runner.dry_run:                                   # a dry run must leave no trace
        state.mark_done(STAGE)
    print(f"[{STAGE}] {len(rows)} backbones in {outdir}")
    return rows
