"""Copies of a predicted monomer placed at the designed subunit positions.

The pipeline never predicts the oligomer. What it has for each sequence is an
AlphaFold prediction of one subunit, and a design backbone that says where the
subunits sit. Superposing the first onto the second gives the assembly the design
asks for, carrying the side chains the sequence implies.

Two callers want that, for different reasons and at different sizes. Stage 05
places two copies to get an interface Rosetta can relax and score; the web page
places one per subunit so the whole ring can be looked at. The placement is the
same operation either way, so it lives here rather than twice.

What this is not: a prediction. The subunit positions are assumed, the fit is
rigid-body, and nothing relaxes the result, so neighbouring subunits may clash.
`fit_rmsd` is the only thing that reports on the placement - how far the
prediction had to be from a design chain to be laid onto it.
"""
from __future__ import annotations

from pathlib import Path


def place_monomer(design_pdb: Path | str, monomer_pdb: Path | str,
                  chains: list[tuple[str, str]] | None = None) -> dict:
    """Place the predicted monomer at the named subunit positions of the design.

    `chains` pairs a chain of the design backbone with the chain id to write the
    copy out as; `None` means every chain of the design, keeping its own label.
    The fit is on CA atoms over the shorter of the two chains, and every heavy
    and light atom of the prediction is carried through the same transform.

    Returns the assembled PDB as text rather than writing it: one caller saves it
    for Rosetta, the other serves it to a browser and has nowhere to put a file.
    """
    import numpy as np

    from .geometry import kabsch, rmsd
    from .pdbio import read_pdb

    design = read_pdb(str(design_pdb))
    monomer = read_pdb(str(monomer_pdb))
    mono = monomer.chains[monomer.chain_ids[0]]

    if chains is None:
        chains = [(cid, cid) for cid in design.chain_ids]
    if not chains:
        raise ValueError(f"{design_pdb} has no chains to place onto")

    # The prediction's atom records, read once: the same lines are rewritten per
    # subunit with new coordinates and a new chain id, and nothing else changes.
    atoms = [line for line in Path(monomer_pdb).read_text(encoding="utf-8").splitlines()
             if line.startswith("ATOM")]
    coords = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                       for line in atoms])

    fits, out = [], []
    for target, label in chains:
        if target not in design.chains:
            raise ValueError(f"{design_pdb} has no chain {target}")
        anchor = design.chains[target].ca
        n = min(len(mono.ca), len(anchor))
        rot, shift = kabsch(mono.ca[:n], anchor[:n])
        fits.append(rmsd(mono.ca[:n] @ rot.T + shift, anchor[:n], superpose=False))
        for line, (x, y, z) in zip(atoms, coords @ rot.T + shift):
            out.append(f"{line[:21]}{label}{line[22:30]}"
                       f"{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}\n")
        out.append("TER\n")

    return {
        "pdb": "".join(out) + "END\n",
        "chains": [label for _, label in chains],
        # The worst of the placements, not the mean: one badly seated subunit is
        # the thing worth knowing about and an average would hide it.
        "fit_rmsd": round(max(fits), 3),
    }
