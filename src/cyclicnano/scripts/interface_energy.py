"""Rosetta interface scoring for one predicted oligomer, run out of process.

PyRosetta lives in its own conda environment, so this cannot be imported by the
pipeline and is invoked as a script instead. It takes a JSON Lines manifest of jobs
and writes one result line per job, so the PyRosetta import and the scorefunction
setup are paid once for a whole run rather than once per structure.

Nothing here imports cyclicnano. The only contract with the pipeline is the shape of
the manifest and of the results, and a job that fails is reported rather than allowed
to take the rest of the batch down with it.

    python interface_energy.py jobs.jsonl results.jsonl [--relax-cycles 1]

Each job: {"id": ..., "pdb": ..., "relaxed_pdb": ...}
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path


def adjacent_pairs(n_chains: int) -> list[tuple[int, int]]:
    """The contacting pairs of a cyclic oligomer: 1-2, 2-3, ... n-1.

    In a Cn ring a subunit touches only its two neighbours; the chains across the
    ring are 13 A away and contribute nothing. Two chains make one pair rather than
    two, since there is only one interface between them.
    """
    if n_chains < 2:
        return []
    if n_chains == 2:
        return [(1, 2)]
    return [(i, i % n_chains + 1) for i in range(1, n_chains + 1)]


def contact_partner(pose) -> tuple[int, float]:
    """The chain that chain 1 actually touches, and how close it gets.

    Not chain 2. AlphaFold numbers chains as it pleases, and on a prediction that
    failed to close into a ring it put chain A next to C and E rather than B; scoring
    1 against 2 then measured the 12.7 A gap across the ring and reported no
    interface at all. That happens precisely on the poor predictions, which are the
    ones worth looking at.

    Ties are broken towards the lower chain index, so a clean Cn ring, where chain 1
    has two equally close neighbours, still returns the same pair every time.
    """
    import numpy as np

    def ca(chain_index):
        begin, end = pose.chain_begin(chain_index), pose.chain_end(chain_index)
        return np.array([[pose.residue(r).xyz("CA")[k] for k in range(3)]
                         for r in range(begin, end + 1)])

    first = ca(1)
    best, best_d = None, float("inf")
    for other in range(2, pose.num_chains() + 1):
        d = float(np.linalg.norm(first[:, None, :] - ca(other)[None, :, :], axis=-1).min())
        if d < best_d - 1e-9:
            best, best_d = other, d
    return best, best_d


def pair_pose(pose, i: int, j: int):
    """A two-chain pose holding chains i and j in their oligomer conformation."""
    from pyrosetta.rosetta.core.pose import append_pose_to_pose
    from pyrosetta.rosetta.protocols.grafting import return_region

    out = return_region(pose.clone(), pose.chain_begin(i), pose.chain_end(i))
    append_pose_to_pose(out, return_region(pose.clone(), pose.chain_begin(j),
                                           pose.chain_end(j)), True)
    return out


def score_one(pdb: str, relaxed_out: str, relax_cycles: int) -> dict:
    """Relax under coordinate restraints, then score the interface between two chains.

    The restraints matter. An unrelaxed AlphaFold model scores around +7000 Rosetta
    energy units on a 300-residue pentamer and around -1000 after relaxing, a gap
    made almost entirely of small clashes the force field hates and the predictor
    never considered. Scoring without relaxing therefore ranks designs by how much
    the two conventions disagree rather than by anything about the interface. The
    restraints keep the relax from walking away from the structure that was predicted.
    """
    from pyrosetta import get_fa_scorefxn, pose_from_pdb
    from pyrosetta.rosetta.core.scoring import ScoreType
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.std import set_int_t

    pose = pose_from_pdb(pdb)
    n_chains = pose.num_chains()
    if n_chains < 2:
        raise ValueError(f"{pdb} has {n_chains} chain(s); an interface needs at least two")

    plain = get_fa_scorefxn()
    before = float(plain(pose))

    # Zero means skip, which FastRelax does not do on its own: given 0 repeats it
    # falls back to a default script that cost 450 seconds a structure, five times
    # the one-cycle price, and still moved an already-relaxed pose by 40 units.
    if relax_cycles > 0:
        restrained = get_fa_scorefxn()
        restrained.set_weight(ScoreType.coordinate_constraint, 1.0)
        relax = FastRelax(restrained, relax_cycles)
        relax.constrain_relax_to_start_coords(True)
        relax.ramp_down_constraints(False)
        relax.apply(pose)
    after = float(plain(pose))
    if relaxed_out:
        Path(relaxed_out).parent.mkdir(parents=True, exist_ok=True)
        pose.dump_pdb(relaxed_out)

    # One interface, between the first two chains.
    #
    # They are taken out of the relaxed oligomer rather than folded on their own, so
    # the two subunits sit in the conformation the assembly puts them in. Scoring
    # chain A against the other four would answer a different question, the cost of
    # taking one subunit out of the ring, which counts both of its contacts: measured
    # on a real structure that came to 2.34 times one interface.
    #
    # The equivalent interfaces of a Cn ring are not quite identical, because the
    # relax runs on the whole oligomer without symmetry restraints and the repacking
    # is stochastic; A-B and A-E differed by 1 to 4 energy units on real structures.
    # That spread is the precision of this number, and it is not averaged away here.
    partner, gap = contact_partner(pose)
    fixed = set_int_t()
    fixed.add(1)
    analyzer = InterfaceAnalyzerMover(fixed, False, plain, True, False, True, False)
    analyzer.apply(pair_pose(pose, 1, partner))
    dg = float(analyzer.get_separated_interface_energy())
    dsasa = float(analyzer.get_interface_delta_sasa())
    return {
        "n_chains": int(n_chains),
        "interface_chains": f"{chr(64 + 1)}_{chr(64 + partner)}",
        "interface_gap": round(gap, 2),
        "dG_interface": dg,
        "dSASA_interface": dsasa,
        "dG_dSASA_100": (dg / dsasa * 100.0) if dsasa else None,
        "total_score_before_relax": before,
        "total_score_after_relax": after,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("jobs")
    ap.add_argument("results")
    ap.add_argument("--relax-cycles", type=int, default=1)
    args = ap.parse_args(argv)

    from pyrosetta import init
    init("-ex1 -ex2aro -use_input_sc -ignore_unrecognized_res -mute all")

    jobs = [json.loads(line) for line in Path(args.jobs).read_text(encoding="utf-8").splitlines()
            if line.strip()]
    out = Path(args.results)
    out.parent.mkdir(parents=True, exist_ok=True)

    done = failed = 0
    with open(out, "w", encoding="utf-8") as fh:
        for i, job in enumerate(jobs, start=1):
            start = time.time()
            try:
                row = score_one(job["pdb"], job.get("relaxed_pdb"), args.relax_cycles)
                row.update({"id": job["id"], "ok": True,
                            "seconds": round(time.time() - start, 1)})
                done += 1
            except Exception as exc:                        # one bad structure is not a run
                row = {"id": job["id"], "ok": False, "error": f"{type(exc).__name__}: {exc}"}
                failed += 1
                print(f"[{i}/{len(jobs)}] {job['id']} FAILED: {exc}", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)
            fh.write(json.dumps(row) + "\n")
            fh.flush()                                      # so a kill leaves what finished
            if row.get("ok"):
                print(f"[{i}/{len(jobs)}] {job['id']} {row['interface_chains']} "
                      f"dG={row['dG_interface']:.1f} "
                      f"dSASA={row['dSASA_interface']:.0f} ({row['seconds']}s)",
                      flush=True)

    print(f"scored {done}, failed {failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
