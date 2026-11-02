"""Stage 05 - interface quality for the sequences that passed.

Stages 03 and 04 answer one question: does the designed sequence fold into the
designed subunit. Both RMSD and pLDDT are properties of a single chain, so a design
can clear them while saying nothing about whether the subunits assemble. This stage
adds the part they leave out.

One measurement:

  dG_separated, from Rosetta's InterfaceAnalyzerMover, for the contact between two
  subunits: the score of the complex minus the two separated and repacked chains, in
  REU. It is not measured on a prediction of the oligomer. The monomer that stage 03
  predicted for this sequence is superposed onto two adjacent chains of the design
  backbone, which gives the designed geometry carrying real side chains, and that
  dimer is relaxed under coordinate restraints and scored.

  A sequence counts as having a usable interface below `interface.ddg_cutoff`, -25 REU
  by default. That number came from the group; nothing here calibrated it.

No GPU is involved. The superposition is linear algebra and the scoring is Rosetta, so
this stage runs on CPU alone.

**What this stage no longer does.** It used to fold the whole oligomer with
AlphaFold2-multimer for an ipTM. That was dropped because the multimer predictions
showed clashes often enough not to be trusted. The removal is worth stating plainly,
because it removes something: ipTM was the only number in this pipeline that asked
whether the sequence would assemble into the ring at all. dG asks a narrower question
- whether the interface that was designed scores well once the side chains the
sequence implies are placed on it - and it asks it on a geometry that is assumed
rather than predicted. Nothing here now tests that assumption. A clashing multimer
prediction is also not unambiguous evidence against the predictor; it can be evidence
against the design. The two cannot be told apart from the prediction alone, which is
why this is a judgement call and not a measurement.

Scope is deliberately narrow. Only backbones that passed stage 04, and within them
only sequences that passed, because this is an evaluation of the designs that got
through and not a second filter applied to everything.
"""
from __future__ import annotations

import json
from pathlib import Path

from .config import Config
from .manifest import RunState
from .runner import Runner

STAGE = "05_interface"

DIMER_SUFFIX = "_dimer.pdb"
DIMER_RELAXED_SUFFIX = "_dimer_relaxed.pdb"

WORKER = Path(__file__).resolve().parent / "scripts" / "interface_energy.py"


def targets(backbones, sequences) -> list[dict]:
    """The sequences to evaluate: the passing ones of the backbones that passed.

    Both filters are needed. Taking every sequence of a passing backbone would spend
    most of the compute on sequences already known not to fold back, and taking every
    passing sequence regardless of its backbone would describe interfaces of designs
    that were not delivered.
    """
    if backbones.empty or sequences.empty:
        return []
    if "final_pass" not in backbones.columns:
        raise RuntimeError("no stage 04 verdict in the manifest; run stage 04 first")
    passed = set(backbones[backbones["final_pass"].fillna(False).astype(bool)]["backbone_id"])
    rows = sequences[sequences["backbone_id"].isin(passed)]
    if "sc_pass" in rows.columns:
        rows = rows[rows["sc_pass"].fillna(False).astype(bool)]
    return rows.to_dict("records")


def _transformed_atoms(pdb_path: Path, rot, shift, chain_id: str) -> list[str]:
    """Every ATOM line of a file, moved by (rot, shift) and relabelled to one chain.

    The transform is derived from CA positions but applied to the raw text, because
    the PDB reader keeps only the backbone and the side chains are the entire point:
    Rosetta cannot score an interface made of glycine, and an AlphaFold monomer is the
    only place in this pipeline where real side chains exist.
    """
    import numpy as np

    out = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines(keepends=True):
        if not line.startswith("ATOM"):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        x, y, z = xyz @ rot.T + shift
        out.append(f"{line[:21]}{chain_id}{line[22:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
    return out


def build_dimer(design_pdb: Path, monomer_pdb: Path, out_path: Path,
                pair: tuple[str, str] = ("A", "B")) -> dict:
    """Place two copies of the predicted monomer at two adjacent subunit positions.

    RFdiffusion writes a cyclic oligomer in ring order, so chains A and B are
    neighbours; measured on real backbones they sit 5.4 A apart while A-C and A-D are
    at 11.3 A. That is unlike a prediction, where the chain labels are arbitrary.
    """
    from .geometry import kabsch, rmsd
    from .pdbio import read_pdb

    design = read_pdb(str(design_pdb))
    monomer = read_pdb(str(monomer_pdb))
    mono = monomer.chains[monomer.chain_ids[0]]

    fits, lines = [], []
    for target, label in zip(pair, ("A", "B")):
        if target not in design.chains:
            raise ValueError(f"{design_pdb} has no chain {target}")
        anchor = design.chains[target].ca
        n = min(len(mono.ca), len(anchor))
        rot, shift = kabsch(mono.ca[:n], anchor[:n])
        fits.append(rmsd(mono.ca[:n] @ rot.T + shift, anchor[:n], superpose=False))
        lines += _transformed_atoms(monomer_pdb, rot, shift, label)
        lines.append("TER\n")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return {"dimer_path": str(out_path), "dimer_fit_rmsd": round(max(fits), 3)}


def backbone_summary(rows: list[dict], cutoff: float) -> list[dict]:
    """The per-backbone interface verdict, from the per-sequence energies.

    A second success rate beside the one stage 04 computes, and deliberately beside
    rather than folded into it: `final_pass` still rests on the monomer alone, so a
    backbone can be read as having failed on the fold or on the interface instead of
    merely as having failed.

    The denominator is the sequences that were scored, not all of them. Those that
    never reached this stage were already rejected on the monomer, and counting them
    here would report the monomer filter a second time under another name.
    """
    per_backbone: dict[str, list[dict]] = {}
    for row in rows:
        bid = str(row.get("backbone_id") or str(row["sequence_id"]).split("__", 1)[0])
        per_backbone.setdefault(bid, []).append(row)

    out = []
    for bid, group in sorted(per_backbone.items()):
        scored = [r for r in group if r.get("dG_interface") is not None]
        passing = [r for r in scored if float(r["dG_interface"]) < cutoff]
        out.append({
            "backbone_id": bid,
            "n_interface_scored": len(scored),
            "n_interface_pass": len(passing),
            "interface_success_rate": (round(len(passing) / len(scored), 4)
                                       if scored else None),
            "best_dG_interface": (min(float(r["dG_interface"]) for r in scored)
                                  if scored else None),
        })
    return out


def run(cfg: Config, state: RunState, runner: Runner, force: bool = False,
        dry_run: bool = False) -> dict:
    cutoff = float(cfg.get("interface.ddg_cutoff", -25.0))
    relax_cycles = str(int(cfg.get("interface.relax_cycles", 1)))

    # As in the other stages: the command is built from paths, so it can be shown on
    # an empty run directory, which is when someone wants to read it.
    if dry_run:
        runner.run("rosetta", ["python", str(WORKER),
                               str(state.outdir / "logs" / "interface_jobs.jsonl"),
                               str(state.outdir / "logs" / "interface_results.jsonl"),
                               "--relax-cycles", relax_cycles])
        return {}
    if state.is_done(STAGE) and not force:
        print(f"[{STAGE}] already done, skipping")
        return {}

    backbones, sequences = state.backbones.load(), state.sequences.load()
    rows = targets(backbones, sequences)
    if not rows:
        print(f"[{STAGE}] no sequence passed stage 04; nothing to evaluate")
        state.mark_done(STAGE)
        return {"n_evaluated": 0}

    work = state.stage_dir(STAGE)
    # Created here rather than assumed. It used to appear as a side effect of
    # writing the multimer FASTA, and removing that prediction left this stage
    # opening a manifest inside a directory nothing had made.
    logs = state.outdir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    print(f"[{STAGE}] {len(rows)} sequences from "
          f"{len({r['backbone_id'] for r in rows})} backbones, CPU only")

    # A sequence that already has an energy is not scored again. The relax is the
    # expensive half of this stage, about a minute and a half a structure, and
    # re-running it to apply a changed cutoff to numbers already on disk would
    # cost hours to produce the same energies. So `--stages 05 --force` re-reads
    # the cutoff cheaply and sends only the unscored sequences to Rosetta.
    existing = {}
    if "dG_interface" in sequences.columns:
        for row in sequences.to_dict("records"):
            energy = row.get("dG_interface")
            if energy is None or energy != energy:          # NaN
                continue
            existing[row["sequence_id"]] = float(energy)

    # The energy is measured on the designed geometry: two copies of the monomer
    # stage 03 predicted, placed at two adjacent subunit positions of the backbone
    # that was designed.
    bb_path = dict(zip(backbones["backbone_id"], backbones["backbone_path"]))
    mono_path = dict(zip(sequences["sequence_id"], sequences.get("pred_path", [])))
    scored, n_built, n_kept = [], 0, 0
    for row in rows:
        sid = row["sequence_id"]
        bid = sid.split("__", 1)[0]
        entry = {"sequence_id": sid, "backbone_id": row.get("backbone_id", bid)}
        if sid in existing:
            entry["dG_interface"] = existing[sid]
            n_kept += 1
            scored.append(entry)
            continue
        design, monomer = bb_path.get(bid), mono_path.get(sid)
        if not design or not monomer or not Path(str(monomer)).exists():
            entry["dimer_error"] = "no stage 03 prediction for this sequence"
            scored.append(entry)
            continue
        short = sid.split("__", 1)[1] if "__" in sid else sid
        try:
            entry.update(build_dimer(Path(str(design)), Path(str(monomer)),
                                     work / bid / f"{short}{DIMER_SUFFIX}"))
            n_built += 1
        except Exception as exc:                            # one bad pair is not a run
            entry["dimer_error"] = f"{type(exc).__name__}: {exc}"
        scored.append(entry)
    if n_kept:
        print(f"[{STAGE}] {n_kept} sequences already had an energy, kept as they are")
    print(f"[{STAGE}] {n_built}/{len(rows)} design-geometry dimers built")

    # Rosetta runs in its own environment, so the structures go over as a manifest and
    # come back as one result per line.
    jobs = logs / "interface_jobs.jsonl"
    results = logs / "interface_results.jsonl"
    with open(jobs, "w", encoding="utf-8") as fh:
        for row in scored:
            dimer = row.get("dimer_path")
            if not dimer:
                continue
            relaxed = str(Path(dimer).with_name(
                Path(dimer).name.replace(DIMER_SUFFIX, DIMER_RELAXED_SUFFIX)))
            row["dimer_relaxed_path"] = relaxed      # the delivery copies this one
            fh.write(json.dumps({"id": row["sequence_id"], "pdb": dimer,
                                 "relaxed_pdb": relaxed}) + "\n")

    if n_built:
        runner.run("rosetta", ["python", str(WORKER), str(jobs), str(results),
                               "--relax-cycles", relax_cycles],
                   check=False)         # a structure Rosetta chokes on must not end the run

    # Only when Rosetta was actually asked. The results file survives a previous
    # run, so reading it unconditionally counted energies carried over from the
    # table as having been computed here, and the log said "43 this run" about a
    # run that invoked Rosetta zero times. Nothing is lost by not reading it: a
    # sequence whose energy is in the table was reused above, and one whose is
    # not had its dimer built, which is what makes n_built non-zero.
    energies = {}
    if n_built and results.exists():
        for line in results.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            if rec.get("ok"):
                energies[rec["id"]] = {k: v for k, v in rec.items()
                                       if k not in ("id", "ok", "seconds")}

    merged = []
    for row in scored:
        entry = dict(row)
        entry.update(energies.get(row["sequence_id"], {}))
        energy = entry.get("dG_interface")
        # Scored means "has an energy", from this run or an earlier one. Testing
        # only this run's results would mark every sequence carried over from a
        # previous stage 05 as unevaluated and throw away its number.
        entry["interface_scored"] = energy is not None
        # None rather than False when there is no energy: "not measured" and
        # "measured and above the cutoff" are different facts, and the delivery
        # renders them apart.
        entry["interface_pass"] = (None if energy is None
                                   else bool(float(energy) < cutoff))
        merged.append(entry)
    state.sequences.upsert(merged)
    state.backbones.upsert(backbone_summary(merged, cutoff))

    # Two counts, because they answer different questions: how many Rosetta
    # returned this time, and how many sequences now hold an energy at all.
    n_new = len(energies)
    n_scored = sum(1 for r in merged if r.get("interface_scored"))
    n_pass = sum(1 for r in merged if r.get("interface_pass"))
    if n_new < n_built:
        print(f"[{STAGE}] WARNING: {n_built - n_new} dimers have no interface "
              f"energy; see {results}")
    print(f"[{STAGE}] {n_scored}/{len(rows)} sequences scored ({n_new} this run), "
          f"{n_pass} below {cutoff} REU")
    state.mark_done(STAGE)
    return {"n_targets": len(rows), "n_scored": n_scored, "n_new": n_new,
            "n_interface_pass": n_pass, "ddg_cutoff": cutoff}
