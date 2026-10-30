"""Stage 05 - interface quality for the sequences that passed.

Stages 03 and 04 answer one question: does the designed sequence fold into the
designed subunit. Both RMSD and pLDDT are properties of a single chain, so a design
can clear them while saying nothing at all about whether the subunits assemble. This
stage adds the part they leave out, and adds it as description rather than as another
gate: nothing here changes success_rate or final_pass.

Two measurements, from one structure:

  ipTM, from an AlphaFold2-multimer prediction of the whole oligomer. The monomer run
  in stage 03 cannot produce it, because an interface predicted TM-score is not
  defined for a single chain. This is a separate, heavier prediction: for C5 it folds
  300 residues rather than 60.

  dG_separated, from Rosetta, for taking one subunit out of the assembled ring. It is
  computed on that same predicted oligomer, after a restrained relax, so the two
  numbers describe the same object and can be read against each other.

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

PRED_SUFFIX = "_af2multimer.pdb"
SCORES_SUFFIX = "_af2multimer_scores.json"
RELAXED_SUFFIX = "_relaxed.pdb"

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


def multimer_fasta(rows: list[dict], n_chains: int, path: Path) -> Path:
    """One record per sequence, the chain repeated and joined by colons.

    That is how colabfold is told to fold a homo-oligomer. The count comes from the
    symmetry of the run, so a C3 gets three copies and a C5 five.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(f">{row['sequence_id']}\n" + ":".join([row["sequence"]] * n_chains) + "\n")
    return path


def build_command(cfg: Config, fasta: Path, outdir: Path) -> list[str]:
    exe = cfg.get("paths.colabfold") or "colabfold_batch"
    argv = [str(exe)]
    params = cfg.get("paths.colabfold_params")
    if params:
        argv += ["--data", str(params)]
    argv += [
        # Still single-sequence. A de novo design has no homologues whether it is
        # folded as one chain or five, and the multimer model does not change that.
        "--msa-mode", str(cfg.get("interface.msa_mode", "single_sequence")),
        "--model-type", str(cfg.get("interface.model_type", "alphafold2_multimer_v3")),
        "--num-recycle", str(cfg.get("interface.num_recycle", 3)),
        "--num-models", str(cfg.get("interface.num_models", 1)),
    ]
    return argv + [str(fasta), str(outdir)]


def collect_prediction(outdir: Path, sequence_id: str) -> dict:
    """ipTM and its companions from the rank-1 multimer prediction."""
    pdbs = sorted(outdir.glob(f"{sequence_id}_unrelaxed_rank_001_*.pdb"))
    if not pdbs:
        return {"interface_pred_path": None, "interface_ok": False}

    row = {"interface_pred_path": str(pdbs[0]), "interface_ok": True}
    scores = sorted(outdir.glob(f"{sequence_id}_scores_rank_001_*.json"))
    if scores:
        data = json.loads(scores[0].read_text(encoding="utf-8"))
        if data.get("iptm") is not None:
            row["iptm"] = float(data["iptm"])
        if data.get("ptm") is not None:
            row["multimer_ptm"] = float(data["ptm"])
        plddt = data.get("plddt") or []
        if plddt:
            row["multimer_plddt"] = float(sum(plddt) / len(plddt))
    return row


def regroup(work: Path, rows: list[dict]) -> list[dict]:
    """Same treatment stage 03 gives its output: one folder per backbone, short names."""
    out = []
    for row in rows:
        sid, pred = row["sequence_id"], row.get("interface_pred_path")
        if not pred or not Path(str(pred)).exists():
            out.append(row)
            continue
        src = Path(str(pred))
        bid = sid.split("__", 1)[0]
        short = sid.split("__", 1)[1] if "__" in sid else sid
        dest = work / bid
        dest.mkdir(parents=True, exist_ok=True)
        new_pdb = dest / f"{short}{PRED_SUFFIX}"
        src.replace(new_pdb)
        for scores in sorted(src.parent.glob(f"{sid}_scores_rank_001_*.json")):
            scores.replace(dest / f"{short}{SCORES_SUFFIX}")
        out.append({**row, "interface_pred_path": str(new_pdb)})
    return out


def run(cfg: Config, state: RunState, runner: Runner, force: bool = False,
        dry_run: bool = False) -> dict:
    # As in stages 02 and 03: the command is built from paths, so it can be shown on
    # an empty run directory, which is when someone wants to read it.
    if dry_run:
        runner.run("colabfold", build_command(cfg, state.outdir / "logs" / "interface.fasta",
                                              state.outdir / STAGE))
        runner.run("rosetta", ["python", str(WORKER),
                               str(state.outdir / "logs" / "interface_jobs.jsonl"),
                               str(state.outdir / "logs" / "interface_results.jsonl"),
                               "--relax-cycles", str(int(cfg.get("interface.relax_cycles", 1)))])
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
    logs = state.outdir / "logs"
    n_chains = cfg.sym_order
    print(f"[{STAGE}] {len(rows)} sequences from "
          f"{len({r['backbone_id'] for r in rows})} backbones, as {n_chains}-mers")

    fasta = multimer_fasta(rows, n_chains, logs / "interface.fasta")
    runner.run("colabfold", build_command(cfg, fasta, work))

    scored = [{"sequence_id": r["sequence_id"], **collect_prediction(work, r["sequence_id"])}
              for r in rows]
    scored = regroup(work, scored)
    n_pred = sum(1 for r in scored if r.get("interface_ok"))
    print(f"[{STAGE}] {n_pred}/{len(rows)} multimer predictions recovered")

    # Rosetta runs in its own environment, so the structures go over as a manifest and
    # come back as one result per line.
    jobs = logs / "interface_jobs.jsonl"
    results = logs / "interface_results.jsonl"
    with open(jobs, "w", encoding="utf-8") as fh:
        for row in scored:
            pred = row.get("interface_pred_path")
            if not pred:
                continue
            relaxed = str(Path(pred).with_name(Path(pred).name.replace(PRED_SUFFIX,
                                                                      RELAXED_SUFFIX)))
            fh.write(json.dumps({"id": row["sequence_id"], "pdb": pred,
                                 "relaxed_pdb": relaxed}) + "\n")

    if n_pred:
        runner.run("rosetta", ["python", str(WORKER), str(jobs), str(results),
                               "--relax-cycles", str(int(cfg.get("interface.relax_cycles", 1)))],
                   check=False)         # a structure Rosetta chokes on must not end the run

    energies = {}
    if results.exists():
        for line in results.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            if rec.get("ok"):
                energies[rec["id"]] = {k: v for k, v in rec.items()
                                       if k not in ("id", "ok", "seconds")}

    merged = []
    for row in scored:
        sid = row["sequence_id"]
        entry = dict(row)
        entry.update(energies.get(sid, {}))
        entry["interface_scored"] = sid in energies
        merged.append(entry)
    state.sequences.upsert(merged)

    n_energy = len(energies)
    if n_energy < n_pred:
        print(f"[{STAGE}] WARNING: {n_pred - n_energy} predictions have no interface "
              f"energy; see {results}")
    print(f"[{STAGE}] {n_energy}/{len(rows)} sequences fully evaluated")
    state.mark_done(STAGE)
    return {"n_targets": len(rows), "n_predicted": n_pred, "n_scored": n_energy}
