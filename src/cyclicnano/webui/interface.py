"""Which residues of a dimer are at the interface.

Stage 05 records how good the interface is and not which residues make it: the
Rosetta worker returns dG, dSASA and a residue count, but no identities. They
are computed here instead, from the relaxed dimer the energy was measured on, so
the answer describes the same structure as the number beside it.

Reported as a list of residue numbers, not drawn. An earlier version also served
the side chains as a structure for the viewer to lay over the dimer; that was
removed as unwanted detail, and the numbers are what is left.

Geometry, so this module uses numpy. `store` stays standard-library-only on
purpose and this is why the work is not in it.

Heavy atoms only, and all of them rather than the backbone. An interface is made
of side chains, and the relaxed dimer is the one structure in this pipeline that
has them; a backbone-only distance would miss most of the contacts it is looking
for.
"""
from __future__ import annotations

from pathlib import Path

# 5 angstrom between heavy atoms. Wide enough to take in a salt bridge or a
# hydrogen bond with its geometry imperfect after a restrained relax, narrow
# enough not to collect the second shell.
DEFAULT_CUTOFF = 5.0
MAX_CUTOFF = 12.0


def parse_atoms(text: str) -> dict:
    """Heavy-atom coordinates per chain, with the residue each one belongs to.

    Read from the raw text rather than through `pdbio`, which keeps only the
    four backbone atoms by name. Hydrogens are dropped by element where the
    column is filled and by name otherwise, since a relaxed Rosetta structure
    may carry them and they would widen every contact by a bond length.
    """
    import numpy as np

    chains: dict[str, dict] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip() if len(line) >= 78 else ""
        name = line[12:16].strip()
        if element == "H" or (not element and name[:1] == "H"):
            continue
        chain = line[21]
        try:
            resnum = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        bucket = chains.setdefault(chain, {"resnum": [], "xyz": []})
        bucket["resnum"].append(resnum)
        bucket["xyz"].append(xyz)

    return {chain: {"resnum": np.array(b["resnum"]), "xyz": np.array(b["xyz"])}
            for chain, b in chains.items() if b["xyz"]}


def contacts(text: str, cutoff: float = DEFAULT_CUTOFF) -> dict:
    """The residues of each chain within `cutoff` of any heavy atom of another.

    Pairwise over every heavy atom, which for a dimer of this size is a few
    hundred by a few hundred and not worth a spatial index. Reported per chain
    and sorted, because the answer is read as a list of residue numbers.
    """
    import numpy as np

    cutoff = max(0.5, min(float(cutoff), MAX_CUTOFF))
    chains = parse_atoms(text)
    labels = sorted(chains)

    found = {label: set() for label in labels}
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            left, right = chains[a], chains[b]
            # (n_a, n_b) distances, then back to residues through the atom index.
            delta = left["xyz"][:, None, :] - right["xyz"][None, :, :]
            close = (delta * delta).sum(-1) <= cutoff * cutoff
            if not close.any():
                continue
            found[a].update(int(r) for r in np.unique(left["resnum"][close.any(axis=1)]))
            found[b].update(int(r) for r in np.unique(right["resnum"][close.any(axis=0)]))

    per_chain = {label: sorted(found[label]) for label in labels if found[label]}
    return {
        "cutoff": cutoff,
        "chains": per_chain,
        "n_residues": sum(len(v) for v in per_chain.values()),
        # Reported so a reader can tell "no contacts" from "only one chain here":
        # a monomer has no interface, which is not the same as a dimer that does
        # not touch.
        "n_chains": len(labels),
    }


def of_file(path: Path | str, cutoff: float = DEFAULT_CUTOFF) -> dict:
    return contacts(Path(path).read_text(encoding="utf-8"), cutoff)
