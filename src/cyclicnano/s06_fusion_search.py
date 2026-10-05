"""Stage 06 - find placements where a C5 and a C3 subunit can be joined.

The icosahedral frame fixes the two axis directions, which leaves four numbers: how
far each component sits along its axis and how far each is turned about it. Those
four determine the whole 60-chain assembly, so the search is over them.

A placement is kept when the two components do not interpenetrate and some pair of
subunit termini, one from each, lies close enough to be bridged by a short
connector. Both directions are tried, since a connector can run from the C5 chain
into the C3 chain or the other way.

The output is a motif PDB per hit: the two subunits alone, in their placed
positions, which is what stage 07 hands to RFdiffusion to build the connector
between.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

from .icosahedral import FIVE_THREE_ANGLE_DEG, axis_pair
from .pdbio import Chain, Structure, read_pdb, write_pdb
from .placement import canonical_frame, component_radius, place, subunit, termini

STAGE = "06_fusion_search"
TERMINUS_NAMES = ("N", "C")


def centre_distance(d5: float, d3: float) -> float:
    """Distance between the two component centres, by the law of cosines."""
    return float(np.sqrt(d5 ** 2 + d3 ** 2
                         - 2 * d5 * d3 * np.cos(np.radians(FIVE_THREE_ANGLE_DEG))))


def distance_bounds(r5: float, r3: float, reach: float) -> tuple[float, float]:
    """Centre separations worth examining.

    Closer than the radii sum and the bodies interpenetrate; further than that plus
    the connector reach and no pair of termini can meet. Bounding the 2D grid this
    way removes most of the 4D search before any coordinate is transformed.
    """
    return r5 + r3 - 0.5 * (r5 + r3), r5 + r3 + reach


def search(c5: Structure, c3: Structure, *,
           link_min: float = 8.0, link_max: float = 30.0,
           turn_step: float = 5.0, shift_step: float = 2.0,
           clash_cutoff: float = 4.0, contact_cutoff: float = 9.0,
           min_contacts: int = 10, max_hits: int = 20000) -> list[dict]:
    """Scan the four degrees of freedom and return the placements worth building on."""
    five, three = axis_pair()
    k5, s5 = canonical_frame(c5)
    k3, s3 = canonical_frame(c3)
    rad5, rad3 = component_radius(k5), component_radius(k3)

    # One turn covers all subunits: turning the C5 by 72 degrees is the same as
    # choosing the next subunit, so the sweep stops there and every subunit is
    # considered at each step.
    turns5 = np.arange(0.0, 360.0 / len(s5), turn_step)
    turns3 = np.arange(0.0, 360.0 / len(s3), turn_step)

    lo, hi = distance_bounds(rad5, rad3, link_max)
    shifts = np.arange(rad5 + rad3 - 10.0, 4.0 * (rad5 + rad3), shift_step)

    # Termini in the canonical frame, before any shift along the axis.
    term5 = np.array([termini(place(k5, five, 0.0, t), s5) for t in turns5])   # (T5, n5, 2, 3)
    term3 = np.array([termini(place(k3, three, 0.0, t), s3) for t in turns3])  # (T3, n3, 2, 3)
    body5 = np.array([place(k5, five, 0.0, t) for t in turns5])
    body3 = np.array([place(k3, three, 0.0, t) for t in turns3])

    # Requiring the components to touch bounds the search by itself: beyond the two
    # radii plus the contact cutoff no contact is possible, so those placements are
    # dropped without transforming a single body coordinate.
    contact_possible_within = rad5 + rad3 + contact_cutoff

    hits: list[dict] = []
    for d5, d3 in itertools.product(shifts, shifts):
        separation = centre_distance(d5, d3)
        if not (lo <= separation <= hi):
            continue
        a = term5 + d5 * five                                 # (T5, n5, 2, 3)
        b = term3 + d3 * three                                # (T3, n3, 2, 3)
        d = np.linalg.norm(a.reshape(-1, 1, 3) - b.reshape(1, -1, 3), axis=-1)
        close = np.argwhere((d >= link_min) & (d <= link_max))
        if not len(close):
            continue

        # Group by placement before touching the bodies: every terminus pair at the
        # same (turn5, turn3) shares one clash check, and there can be thousands of
        # pairs per placement.
        by_placement: dict[tuple[int, int], list] = {}
        for flat5, flat3 in close:
            i5, rest5 = divmod(flat5, len(s5) * 2)
            sub5, end5 = divmod(rest5, 2)
            i3, rest3 = divmod(flat3, len(s3) * 2)
            sub3, end3 = divmod(rest3, 2)
            if end5 == end3:
                continue                                      # N to N or C to C cannot be joined
            by_placement.setdefault((i5, i3), []).append(
                (sub5, end5, sub3, end3, float(d[flat5, flat3])))

        if separation > contact_possible_within:
            continue                                           # the two bodies cannot meet

        for (i5, i3), pairs in by_placement.items():
            p5 = body5[i5] + d5 * five
            p3 = body3[i3] + d3 * three
            gap = np.linalg.norm(p5[:, None, :] - p3[None, :, :], axis=-1)
            min_gap = float(gap.min())
            if min_gap < clash_cutoff:
                continue
            # The covalent connector alone does not make a rigid particle. Without an
            # interface between the two components the assembly is a string of beads
            # held by a linker, so contact is required rather than merely permitted.
            n_contacts = int((gap < contact_cutoff).sum())
            if n_contacts < min_contacts:
                continue

            sub5, end5, sub3, end3, link = min(pairs, key=lambda x: x[4])
            hits.append({
                "d5": float(d5), "turn5": float(turns5[i5]),
                "d3": float(d3), "turn3": float(turns3[i3]),
                "subunit5": int(sub5), "subunit3": int(sub3),
                "terminus5": TERMINUS_NAMES[end5], "terminus3": TERMINUS_NAMES[end3],
                "link_distance": link,
                "centre_distance": separation,
                "min_body_gap": min_gap,
                "n_contacts": n_contacts,
                "n_terminus_pairs": len(pairs),
            })
            if len(hits) >= max_hits:
                return hits
    return hits


def motif_structure(c5: Structure, c3: Structure, hit: dict) -> Structure:
    """The two subunits alone, placed, as the motif RFdiffusion builds between."""
    five, three = axis_pair()
    k5, s5 = canonical_frame(c5)
    k3, s3 = canonical_frame(c3)
    a = subunit(place(k5, five, hit["d5"], hit["turn5"]), s5, hit["subunit5"])
    b = subunit(place(k3, three, hit["d3"], hit["turn3"]), s3, hit["subunit3"])
    return Structure(chains={
        "A": Chain("A", np.arange(1, len(a) + 1), {"CA": a}),
        "B": Chain("B", np.arange(1, len(b) + 1), {"CA": b}),
    })


def write_hits(outdir: Path, c5: Structure, c3: Structure, hits: list[dict]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for i, hit in enumerate(hits):
        write_pdb(str(outdir / f"motif_{i:04d}.pdb"), motif_structure(c5, c3, hit))
        hit["motif"] = f"motif_{i:04d}.pdb"
    (outdir / "hits.json").write_text(json.dumps(hits, indent=2), encoding="utf-8")
