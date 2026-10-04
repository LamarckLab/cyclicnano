"""Metric-layer tests against synthetic oligomers with known ground truth."""
from __future__ import annotations

import numpy as np

from cyclicnano.geometry import (backbone_metrics, count_clashes, gyration_shape, rmsd,
                            secondary_structure, symmetry_frame, terminal_run)
from synthetic import (clashing_oligomer, cyclic_oligomer, extended_monomer,
                       helical_hairpin, spherical_cloud)


# ------------------------------------------------------------------ symmetry
def test_symmetry_order_recovered():
    for n in (3, 5, 6, 8):
        detected = symmetry_frame(cyclic_oligomer(n_sym=n))["sym_order_detected"]
        assert abs(detected - n) < 0.01, f"C{n}: detected {detected}"


def test_symmetry_chains_superpose_exactly():
    for n in (3, 5, 6):
        assert symmetry_frame(cyclic_oligomer(n_sym=n))["sym_rmsd"] < 1e-6


def test_symmetry_axis_is_z():
    axis = symmetry_frame(cyclic_oligomer(n_sym=5))["axis"]
    assert abs(abs(float(np.dot(axis, [0, 0, 1]))) - 1.0) < 1e-6


def test_expected_symmetry_mismatch_is_flagged():
    m = backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=3)
    assert m["sym_order_ok"] is False
    assert backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=5)["sym_order_ok"] is True


# ------------------------------------------------------------------ rmsd
def test_rmsd_is_zero_for_identical_coordinates():
    ca = cyclic_oligomer(n_sym=5).chains["A"].ca
    assert rmsd(ca, ca) < 1e-9


def test_rmsd_invariant_under_rigid_motion():
    struct = cyclic_oligomer(n_sym=5)
    a = struct.chains["A"].ca
    assert rmsd(struct.chains["B"].ca, a) < 1e-6           # symmetry mate is a pure rotation
    assert rmsd(a + np.array([3.0, -2.0, 7.0]), a) < 1e-6


def test_rmsd_rejects_length_mismatch():
    ca = cyclic_oligomer(n_sym=5).chains["A"].ca
    try:
        rmsd(ca, ca[:-5])
        raise AssertionError("expected a shape mismatch error")
    except ValueError:
        pass


# ------------------------------------------------------------------ secondary structure
def test_helix_bundle_is_mostly_helix():
    ss = secondary_structure(cyclic_oligomer(n_sym=5).chains["A"].ca)
    assert ss.count("H") / len(ss) > 0.70
    assert ss.count("E") == 0                              # an all-alpha bundle has no strand
    assert "L" in ss                                       # loops must not be called helix


def test_terminal_run():
    assert terminal_run("HHHLLLHH", "H", "N") == 3
    assert terminal_run("HHHLLLHH", "H", "C") == 2
    assert terminal_run("LHHH", "H", "N") == 0


# ------------------------------------------------------------------ clashes
def test_separated_ring_is_clash_free():
    m = backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=5)
    assert m["n_clash"] == 0
    assert m["iface_contacts_per_chain"] > 0               # but the subunits do touch


def test_interpenetrating_ring_is_flagged():
    m = backbone_metrics(clashing_oligomer(n_sym=5), expected_sym=5)
    assert m["n_clash_inter"] > 0
    assert m["min_interchain_dist"] < 2.5


def test_interface_contacts_are_not_counted_as_clashes():
    """Regression: a 4 A CA-CA proxy flagged ordinary interface contacts as overlap."""
    struct = cyclic_oligomer(n_sym=5, ring_radius=16.0)
    clashes = count_clashes(struct)
    assert clashes["n_clash_inter"] == 0
    assert clashes["min_interchain_dist"] > 2.5


# ------------------------------------------------------------------ shape
def test_extended_backbone_is_flagged_by_compactness():
    assert backbone_metrics(extended_monomer(60))["rg_ratio"] > 1.5


def test_compact_bundle_passes_compactness():
    assert backbone_metrics(cyclic_oligomer(n_sym=5))["rg_ratio"] < 1.3


def test_ring_geometry_is_self_consistent():
    m = backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=5)
    assert 0 < m["pore_radius"] < m["max_radius"]
    assert m["height"] > 0


def test_metric_set_is_complete():
    m = backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=5)
    required = {
        "n_chains", "n_res_monomer", "helix_frac", "strand_frac", "loop_frac", "n_sse",
        "loop_max_len", "nterm_helix_len", "cterm_helix_len", "rg", "rg_ratio",
        "contact_order", "iface_contacts_per_chain", "pore_radius", "max_radius",
        "height", "sym_angle_deg", "sym_order_detected", "sym_rmsd", "n_clash",
        "n_clash_intra", "n_clash_inter", "min_interchain_dist", "sym_order_ok",
        "n_helices", "asphericity", "acylindricity", "shape_anisotropy", "axis_ratio",
        "assembly_asphericity", "assembly_shape_anisotropy", "assembly_axis_ratio",
        "max_axis_angle", "mean_axis_angle", "max_inter_helix_angle",
        "max_helix_len", "max_strand_len",
    }
    assert required <= set(m), f"missing: {sorted(required - set(m))}"


# ------------------------------------------------------------------ shape
def test_shape_anisotropy_ranks_sphere_bundle_hairpin_rod():
    """kappa squared must order the four reference shapes monotonically."""
    k = lambda ca: gyration_shape(ca)["shape_anisotropy"]        # noqa: E731
    sphere = k(spherical_cloud().chains["A"].ca)
    bundle = k(cyclic_oligomer(n_sym=5).chains["A"].ca)
    hairpin = k(helical_hairpin())
    rod = k(extended_monomer(60).chains["A"].ca)
    assert sphere < bundle < hairpin < rod, (sphere, bundle, hairpin, rod)


def test_ideal_sphere_is_near_zero_anisotropy():
    assert gyration_shape(spherical_cloud().chains["A"].ca)["shape_anisotropy"] < 0.05


def test_rod_is_near_one_anisotropy():
    assert gyration_shape(extended_monomer(60).chains["A"].ca)["shape_anisotropy"] > 0.9


def test_three_helix_bundle_is_compact_and_all_alpha():
    m = backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=5)
    assert m["n_helices"] == 3
    assert m["shape_anisotropy"] <= 0.30
    assert m["strand_frac"] <= 0.05


def test_two_helix_hairpin_is_more_elongated_than_a_bundle():
    from cyclicnano.pdbio import Chain, Structure
    mono = helical_hairpin()
    struct = Structure(chains={"A": Chain("A", np.arange(1, len(mono) + 1), {"CA": mono})})
    m = backbone_metrics(struct)
    assert m["shape_anisotropy"] > 0.30 or m["axis_ratio"] > 3.0


def test_assembly_shape_is_measured_separately_from_the_subunit():
    """The design target is a globular ring, which is not the same as a globular subunit."""
    m = backbone_metrics(cyclic_oligomer(n_sym=5), expected_sym=5)
    assert m["assembly_shape_anisotropy"] != m["shape_anisotropy"]
    assert 0.0 <= m["assembly_shape_anisotropy"] <= 1.0


# ------------------------------------------------------------------ axis alignment
def test_axis_angle_is_zero_for_helices_built_along_the_axis():
    """The synthetic bundle runs its helices parallel to z, which is the Cn axis."""
    from cyclicnano.geometry import axis_alignment, secondary_structure, symmetry_frame
    struct = cyclic_oligomer(n_sym=5)
    ca = struct.chains["A"].ca
    a = axis_alignment(ca, secondary_structure(ca), symmetry_frame(struct)["axis"])
    assert a["max_axis_angle"] < 10            # built parallel, so near zero
    assert a["max_inter_helix_angle"] < 10     # and parallel to each other


def test_axis_angle_is_ninety_for_a_helix_across_the_axis():
    """A helix perpendicular to the symmetry axis is the failure case the rule targets."""
    import numpy as np
    from cyclicnano.geometry import axis_alignment
    from synthetic import ideal_helix
    ca = ideal_helix(20)                        # runs along z
    ss = "H" * 20
    across = axis_alignment(ca, ss, np.array([1.0, 0.0, 0.0]))
    assert across["max_axis_angle"] > 80
    along = axis_alignment(ca, ss, np.array([0.0, 0.0, 1.0]))
    assert along["max_axis_angle"] < 10


def test_axis_angle_ignores_helix_direction():
    """Up the axis and down it are equally well aligned, so angles fold into 0-90."""
    import numpy as np
    from cyclicnano.geometry import axis_alignment
    from synthetic import ideal_helix
    ca, ss = ideal_helix(20), "H" * 20
    up = axis_alignment(ca, ss, np.array([0.0, 0.0, 1.0]))["max_axis_angle"]
    down = axis_alignment(ca, ss, np.array([0.0, 0.0, -1.0]))["max_axis_angle"]
    assert abs(up - down) < 1e-6


def test_segment_lengths_are_reported():
    from cyclicnano.geometry import sse_segments
    segs = sse_segments("LLHHHHHHHHLLEEEEEELLHHH")   # 8-helix, 6-strand, 3-helix below min_len
    assert [(c, b - a) for c, a, b in segs] == [("H", 8), ("E", 6)]
