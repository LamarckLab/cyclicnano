"""Objective backbone metrics for cyclic homo-oligomers.

Every quantity here is deterministic and dependency-free (numpy only). The point of
this module is that a human judgement like "this backbone looks bad" becomes a number
that can be thresholded in a config file and re-tuned without recomputing anything.
"""
from __future__ import annotations

import numpy as np

from .pdbio import Structure

# --- P-SEA secondary-structure criteria (Labesse et al. 1997), CA-trace only ---
PSEA = {
    "helix": {"d2": (5.5, 0.5), "d3": (5.3, 0.5), "d4": (6.4, 0.6),
              "theta": (89.0, 12.0), "tau": (50.0, 20.0)},
    "strand": {"d2": (6.7, 0.6), "d3": (9.9, 0.9), "d4": (12.4, 1.1),
               "theta": (124.0, 14.0), "tau": (-170.0, 45.0)},
}
MIN_RUN = {"H": 4, "E": 3}          # shorter runs get demoted to loop


# ---------------------------------------------------------------- superposition
def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (R, t) such that target ~= mobile @ R.T + t."""
    mc, tc = mobile.mean(0), target.mean(0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, tc - r @ mc


def rmsd(a: np.ndarray, b: np.ndarray, superpose: bool = True) -> float:
    """Cartesian RMSD between two equal-length coordinate sets."""
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    if superpose:
        r, t = kabsch(a, b)
        a = a @ r.T + t
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


# ---------------------------------------------------------------- local geometry
def _angles(ca: np.ndarray) -> np.ndarray:
    """Virtual bond angle at each CA(i) using CA(i-1), CA(i), CA(i+1), in degrees."""
    out = np.full(len(ca), np.nan)
    v1, v2 = ca[:-2] - ca[1:-1], ca[2:] - ca[1:-1]
    cos = (v1 * v2).sum(1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
    out[1:-1] = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    return out


def _torsions(ca: np.ndarray) -> np.ndarray:
    """Virtual CA torsion at each CA(i) over CA(i-1..i+2), in degrees."""
    out = np.full(len(ca), np.nan)
    p0, p1, p2, p3 = ca[:-3], ca[1:-2], ca[2:-1], ca[3:]
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - (b0 * b1).sum(1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(1, keepdims=True) * b1
    out[1:-2] = np.degrees(np.arctan2((np.cross(b1, v) * w).sum(1), (v * w).sum(1)))
    return out


def _within(value: float, spec: tuple[float, float]) -> bool:
    centre, tol = spec
    return bool(abs(value - centre) <= tol)


def _enforce_min_run(labels: list[str]) -> list[str]:
    """Demote secondary-structure runs shorter than the P-SEA minimum to loop."""
    out, i = list(labels), 0
    while i < len(out):
        j = i
        while j < len(out) and out[j] == out[i]:
            j += 1
        if out[i] in MIN_RUN and (j - i) < MIN_RUN[out[i]]:
            for k in range(i, j):
                out[k] = "L"
        i = j
    return out


def secondary_structure(ca: np.ndarray) -> str:
    """P-SEA-style H/E/L string from a CA trace.

    Approximate by construction: it needs no external binary, which keeps the metric
    layer reproducible across machines. A DSSP backend can be slotted in later.
    """
    n = len(ca)
    if n < 6:
        return "L" * n

    def gap(k: int) -> np.ndarray:
        d = np.full(n, np.nan)
        d[: n - k] = np.linalg.norm(ca[k:] - ca[: n - k], axis=1)
        return d

    d2, d3, d4 = gap(2), gap(3), gap(4)
    theta, tau = _angles(ca), _torsions(ca)

    labels = []
    for i in range(n):
        lab = "L"
        for code, key in (("H", "helix"), ("E", "strand")):
            c = PSEA[key]
            dist_ok = all(
                not np.isnan(v) and _within(v, c[k])
                for v, k in ((d2[i], "d2"), (d3[i], "d3"), (d4[i], "d4"))
            )
            ang_ok = (
                not np.isnan(theta[i]) and not np.isnan(tau[i])
                and _within(theta[i], c["theta"]) and _within(tau[i], c["tau"])
            )
            if dist_ok or ang_ok:
                lab = code
                break
        labels.append(lab)
    return "".join(_enforce_min_run(labels))


def _runs(ss: str, code: str) -> list[int]:
    lengths, count = [], 0
    for ch in list(ss) + [""]:
        if ch == code:
            count += 1
        else:
            if count:
                lengths.append(count)
            count = 0
    return lengths


def sse_segments(ss: str, min_len: int = 4) -> list[tuple[str, int, int]]:
    """Contiguous helix and strand runs as (code, start, stop), loops excluded."""
    out, i = [], 0
    while i < len(ss):
        j = i
        while j < len(ss) and ss[j] == ss[i]:
            j += 1
        if ss[i] in "HE" and j - i >= min_len:
            out.append((ss[i], i, j))
        i = j
    return out


def _principal_axis(points: np.ndarray) -> np.ndarray:
    """Direction of greatest extent, i.e. the long axis of a helix or strand."""
    centred = points - points.mean(0)
    _, _, vt = np.linalg.svd(centred)
    return vt[0]


def axis_alignment(ca: np.ndarray, ss: str, sym_axis: np.ndarray) -> dict:
    """How closely each secondary structure element runs along the symmetry axis.

    This is the criterion that decides whether a subunit can form a stable
    interface. A helix lying parallel to the axis extends along it and can pack
    against the matching helix of the neighbouring subunit; one that is splayed
    meets its neighbour at a point instead of along a line, and the assembly has
    nothing to hold it together.

    Angles are unsigned and folded into 0-90 degrees, since a helix pointing up
    the axis and one pointing down are equally well aligned.
    """
    segments = sse_segments(ss)
    if not segments:
        return {"max_axis_angle": float("nan"), "mean_axis_angle": float("nan"),
                "max_inter_helix_angle": float("nan"),
                "max_helix_len": 0, "max_strand_len": 0}

    angles, helix_dirs = [], []
    for code, start, stop in segments:
        direction = _principal_axis(ca[start:stop])
        angle = np.degrees(np.arccos(min(1.0, abs(float(np.dot(direction, sym_axis))))))
        angles.append(angle)
        if code == "H":
            helix_dirs.append(direction)

    inter = [
        np.degrees(np.arccos(min(1.0, abs(float(np.dot(helix_dirs[i], helix_dirs[j]))))))
        for i in range(len(helix_dirs)) for j in range(i + 1, len(helix_dirs))
    ]
    return {
        "max_axis_angle": float(max(angles)),
        "mean_axis_angle": float(np.mean(angles)),
        "max_inter_helix_angle": float(max(inter)) if inter else 0.0,
        "max_helix_len": max((b - a for c, a, b in segments if c == "H"), default=0),
        "max_strand_len": max((b - a for c, a, b in segments if c == "E"), default=0),
    }


def terminal_run(ss: str, code: str, end: str) -> int:
    """Length of the leading (N) or trailing (C) run of the given SS code, else 0."""
    seq = ss if end == "N" else ss[::-1]
    n = 0
    for ch in seq:
        if ch != code:
            break
        n += 1
    return n


# ---------------------------------------------------------------- symmetry
def symmetry_frame(structure: Structure) -> dict:
    """Recover the rotation relating chain 0 to chain 1, and the resulting Cn axis.

    Also reports how cleanly the chains superpose, which catches generation runs where
    the enforced symmetry drifted.
    """
    ids = structure.chain_ids
    if len(ids) < 2:
        return {"sym_angle_deg": 0.0, "sym_order_detected": 1.0, "sym_rmsd": 0.0,
                "axis": np.array([0.0, 0.0, 1.0]), "centre": structure.all_ca().mean(0)}

    a, b = structure.chains[ids[0]].ca, structure.chains[ids[1]].ca
    n = min(len(a), len(b))
    r, _ = kabsch(a[:n], b[:n])

    cos_angle = float(np.clip((np.trace(r) - 1.0) / 2.0, -1.0, 1.0))
    angle = float(np.degrees(np.arccos(cos_angle)))
    order = 360.0 / angle if angle > 1e-6 else float("inf")

    eigval, eigvec = np.linalg.eig(r)                       # rotation axis = eigenvalue-1 vector
    axis = np.real(eigvec[:, int(np.argmin(np.abs(eigval - 1.0)))])
    axis = axis / np.linalg.norm(axis)

    ref = structure.chains[ids[0]].ca
    devs = [rmsd(structure.chains[c].ca[: len(ref)], ref) for c in ids[1:]]

    return {
        "sym_angle_deg": angle,
        "sym_order_detected": order,
        "sym_rmsd": float(np.max(devs)) if devs else 0.0,
        "axis": axis,
        "centre": structure.all_ca().mean(0),
    }


# ---------------------------------------------------------------- clashes & contacts
def _pair_dists(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)


def _backbone_stack(chain) -> tuple[np.ndarray, np.ndarray]:
    """All available backbone atoms of one chain, plus the residue index of each."""
    names = [n for n in ("N", "CA", "C", "O") if n in chain.coords]
    return (np.concatenate([chain.coords[n] for n in names]),
            np.tile(np.arange(chain.n_res), len(names)))


def count_clashes(structure: Structure, cutoff: float = 2.5) -> dict:
    """Backbone heavy-atom clashes, split into intra- and inter-chain.

    One criterion is used throughout: two backbone atoms closer than `cutoff` overlap.
    A CA-CA proxy is deliberately avoided for the inter-chain case, because ordinary
    interface contacts sit at 4-6 A and would otherwise all register as clashes.
    """
    intra = 0
    for chain in structure.chains.values():
        stacked, res_idx = _backbone_stack(chain)
        d = _pair_dists(stacked, stacked)
        sep = np.abs(res_idx[:, None] - res_idx[None, :])
        intra += int(np.triu((d < cutoff) & (sep >= 2), k=1).sum())

    inter, min_inter, ids = 0, float("inf"), structure.chain_ids
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, _ = _backbone_stack(structure.chains[ids[i]])
            b, _ = _backbone_stack(structure.chains[ids[j]])
            d = _pair_dists(a, b)
            inter += int((d < cutoff).sum())
            min_inter = min(min_inter, float(d.min()))
    return {
        "n_clash_intra": intra,
        "n_clash_inter": inter,
        "n_clash": intra + inter,
        "min_interchain_dist": min_inter if np.isfinite(min_inter) else float("nan"),
    }


def interface_contacts(structure: Structure, cutoff: float = 8.0) -> float:
    """Mean number of inter-chain CA-CA contacts per chain."""
    ids = structure.chain_ids
    if len(ids) < 2:
        return 0.0
    total = 0
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            total += int((_pair_dists(structure.chains[ids[i]].ca,
                                      structure.chains[ids[j]].ca) < cutoff).sum())
    return float(2 * total / len(ids))


def gyration_shape(ca: np.ndarray) -> dict:
    """Shape descriptors from the gyration tensor eigenvalues.

    `shape_anisotropy` (kappa squared) is the direct answer to "is this subunit
    globular": 0 for a perfect sphere, 1 for a straight rod. Compact globular folds
    sit near 0.05-0.15, a two-helix hairpin well above that. `rg_ratio` only measures
    overall size against a reference and cannot distinguish a compact ball from an
    equally sized disc, which is why this is computed separately.
    """
    centred = ca - ca.mean(0)
    tensor = (centred[:, :, None] * centred[:, None, :]).mean(axis=0)
    lam = np.sort(np.linalg.eigvalsh(tensor))            # lam[0] <= lam[1] <= lam[2]
    rg_sq = float(lam.sum())
    if rg_sq <= 0:
        return {"asphericity": 0.0, "acylindricity": 0.0,
                "shape_anisotropy": 0.0, "axis_ratio": 1.0}

    b = float(lam[2] - 0.5 * (lam[0] + lam[1]))          # asphericity
    c = float(lam[1] - lam[0])                           # acylindricity
    return {
        "asphericity": b / rg_sq,
        "acylindricity": c / rg_sq,
        "shape_anisotropy": (b ** 2 + 0.75 * c ** 2) / (rg_sq ** 2),
        "axis_ratio": float(np.sqrt(lam[2] / lam[0])) if lam[0] > 1e-9 else float("inf"),
    }


def relative_contact_order(ca: np.ndarray, cutoff: float = 8.0, min_sep: int = 3) -> float:
    """Relative contact order; high values flag topologies that fold poorly."""
    d = _pair_dists(ca, ca)
    idx = np.arange(len(ca))
    sep = np.abs(idx[:, None] - idx[None, :])
    mask = np.triu((d < cutoff) & (sep >= min_sep), k=1)
    n = int(mask.sum())
    return float(sep[mask].sum() / (n * len(ca))) if n else 0.0


# ---------------------------------------------------------------- top-level
def backbone_metrics(structure: Structure, expected_sym: int | None = None) -> dict:
    """Full objective metric set for one generated oligomer backbone."""
    ids = structure.chain_ids
    ca = structure.chains[ids[0]].ca
    n_res = len(ca)

    ss = secondary_structure(ca)
    loops = _runs(ss, "L")
    frame = symmetry_frame(structure)

    axis, centre = frame["axis"], frame["centre"]
    rel = structure.all_ca() - centre
    along = rel @ axis
    radial = np.linalg.norm(rel - np.outer(along, axis), axis=1)

    rg = float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(axis=1).mean()))
    rg_ref = 2.2 * n_res ** 0.38                            # compact-globular reference

    m = {
        "n_chains": structure.n_chains,
        "n_res_monomer": n_res,
        "ss_string": ss,
        "helix_frac": ss.count("H") / n_res,
        "strand_frac": ss.count("E") / n_res,
        "loop_frac": ss.count("L") / n_res,
        "n_sse": len(_runs(ss, "H")) + len(_runs(ss, "E")),
        "loop_max_len": max(loops) if loops else 0,
        "nterm_helix_len": terminal_run(ss, "H", "N"),      # recorded now, used by the fusion stage later
        "cterm_helix_len": terminal_run(ss, "H", "C"),
        "rg": rg,
        "rg_ratio": rg / rg_ref,
        "n_helices": len(_runs(ss, "H")),                   # subunit topology, e.g. 3 for a 3-helix fold
        "contact_order": relative_contact_order(ca),
        "iface_contacts_per_chain": interface_contacts(structure),
        "pore_radius": float(radial.min()),
        "max_radius": float(radial.max()),
        "height": float(along.max() - along.min()),
        "sym_angle_deg": frame["sym_angle_deg"],
        "sym_order_detected": frame["sym_order_detected"],
        "sym_rmsd": frame["sym_rmsd"],
    }
    m.update(axis_alignment(ca, ss, axis))
    m.update(gyration_shape(ca))
    # The stated design target is that the whole ring is globular, not just one
    # subunit, so the same descriptors are computed over the full assembly.
    m.update({f"assembly_{k}": v for k, v in gyration_shape(structure.all_ca()).items()})
    m.update(count_clashes(structure))
    if expected_sym is not None:
        m["sym_order_ok"] = bool(abs(m["sym_order_detected"] - expected_sym) < 0.1)
    return m
