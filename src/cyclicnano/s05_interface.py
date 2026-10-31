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

  dG_interface, from Rosetta, for the contact between two subunits. It is NOT measured
  on the prediction. The monomer that stage 03 predicted for this sequence is
  superposed onto two adjacent chains of the design backbone, giving the designed
  geometry carrying real side chains, and that dimer is relaxed and scored.

  The two therefore describe different structures, deliberately. ipTM asks whether the
  sequence would assemble this way at all; dG asks whether the interface we designed is
  any good once the side chains the sequence implies are put on it. A design can do well
  on one and badly on the other, and that difference is information rather than noise.

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


DIMER_SUFFIX = "_dimer.pdb"
DIMER_RELAXED_SUFFIX = "_dimer_relaxed.pdb"


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
    at 11.3 A. That is unlike the predictions, where the chain labels are arbitrary.
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

    # The same tidying stage 03 does, which this stage was missing: colabfold leaves
    # six files per sequence beside the two that are kept, and on 328 sequences that
    # was 1971 loose files and 260 of the 319 MB. The coverage plot and the a3m are
    # particularly empty here, since single-sequence mode means the alignment is the
    # input sequence and nothing else.
    if not cfg.get("interface.keep_all_outputs", False):
        archive = logs / "colabfold_multimer"
        removed = 0
        for path in work.iterdir():
            if not path.is_file():                          # the regrouped folders stay
                continue
            if path.name in ("config.json", "cite.bibtex", "log.txt"):
                archive.mkdir(parents=True, exist_ok=True)
                path.replace(archive / path.name)
            else:
                path.unlink()
                removed += 1
        if removed:
            print(f"[{STAGE}] removed {removed} unused colabfold products")

    # The energy is measured on the designed geometry, not on the prediction: two
    # copies of the monomer stage 03 predicted, placed at two adjacent subunit
    # positions of the backbone that was designed. Nothing here needs a GPU.
    bb_path = dict(zip(backbones["backbone_id"], backbones["backbone_path"]))
    mono_path = dict(zip(sequences["sequence_id"], sequences.get("pred_path", [])))
    n_built = 0
    for row in scored:
        sid = row["sequence_id"]
        bid = sid.split("__", 1)[0]
        design, monomer = bb_path.get(bid), mono_path.get(sid)
        if not design or not monomer or not Path(str(monomer)).exists():
            continue
        short = sid.split("__", 1)[1] if "__" in sid else sid
        try:
            row.update(build_dimer(Path(str(design)), Path(str(monomer)),
                                   work / bid / f"{short}{DIMER_SUFFIX}"))
            n_built += 1
        except Exception as exc:                            # one bad pair is not a run
            row["dimer_error"] = f"{type(exc).__name__}: {exc}"
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
    if n_energy < n_built:
        print(f"[{STAGE}] WARNING: {n_built - n_energy} dimers have no interface "
              f"energy; see {results}")
    print(f"[{STAGE}] {n_energy}/{len(rows)} sequences fully evaluated")
    state.mark_done(STAGE)
    return {"n_targets": len(rows), "n_predicted": n_pred, "n_scored": n_energy}
