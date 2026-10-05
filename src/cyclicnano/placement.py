"""Placing a cyclic component onto an icosahedral axis.

A C5 component goes on a fivefold axis and a C3 on a threefold axis. The axis
directions are fixed by the geometry, which leaves each component two freedoms:
how far along its axis it sits, and how far it is turned about that axis. Four
numbers in total, and they determine the entire 60-chain assembly.

The work here is to put a component into a known frame first. An oligomer as
RFdiffusion writes it sits wherever the diffusion left it, so its own symmetry
axis is recovered, moved to +z and its centroid to the origin. From that canonical
form a placement is one turn and one shift.
"""
from __future__ import annotations

import numpy as np

from .geometry import symmetry_frame
from .icosahedral import rotation_matrix
from .pdbio import Structure


def _align(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotation carrying unit vector `source` onto unit vector `target`."""
    v = np.cross(source, target)
    c = float(np.dot(source, target))
    if np.linalg.norm(v) < 1e-12:                       # already parallel or antiparallel
        if c > 0:
            return np.eye(3)
        perp = np.array([1.0, 0.0, 0.0])
        if abs(source[0]) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        axis = np.cross(source, perp)
        return rotation_matrix(axis / np.linalg.norm(axis), np.pi)
    k = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3) + k + k @ k / (1.0 + c)


def canonical_frame(structure: Structure) -> tuple[np.ndarray, list[int]]:
    """Per-chain CA coordinates with the symmetry axis on +z and the centroid at 0.

    Returns the stacked coordinates and the chain boundaries, so a single subunit
    can be sliced out later without re-deriving which atoms belong to it.
    """
    frame = symmetry_frame(structure)
    rot = _align(frame["axis"], np.array([0.0, 0.0, 1.0]))

    chains, sizes = [], []
    for cid in structure.chain_ids:
        ca = structure.chains[cid].ca
        chains.append((ca - frame["centre"]) @ rot.T)
        sizes.append(len(ca))
    return np.concatenate(chains), sizes


def place(canonical: np.ndarray, axis: np.ndarray, distance: float, turn_deg: float) -> np.ndarray:
    """Turn about the component axis, then carry it onto `axis` at `distance`.

    The turn is applied in the canonical frame, where the component axis is +z, so
    it is a plain rotation about z and does not interact with the alignment.
    """
    spun = canonical @ rotation_matrix(np.array([0.0, 0.0, 1.0]), np.radians(turn_deg)).T
    onto = _align(np.array([0.0, 0.0, 1.0]), axis)
    return spun @ onto.T + distance * axis


def subunit(coords: np.ndarray, sizes: list[int], index: int) -> np.ndarray:
    """One chain out of a stacked component."""
    start = sum(sizes[:index])
    return coords[start:start + sizes[index]]


def termini(coords: np.ndarray, sizes: list[int]) -> np.ndarray:
    """First and last CA of every chain, as (n_chains, 2, 3).

    These are the only atoms the search needs: a connector is built between one
    terminus of a C5 subunit and one of a C3 subunit.
    """
    out = []
    for i in range(len(sizes)):
        ca = subunit(coords, sizes, i)
        out.append([ca[0], ca[-1]])
    return np.array(out)


def component_radius(canonical: np.ndarray) -> float:
    """Farthest CA from the origin, used to bound the search and to prune clashes."""
    return float(np.linalg.norm(canonical, axis=1).max())
