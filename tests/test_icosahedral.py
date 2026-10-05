"""The icosahedral rotation group.

Everything in the assembly stages rests on these 60 matrices and on the two axes
drawn from them, and a mistake here is quiet: a slightly wrong axis still yields a
plausible-looking shell that simply fails to close. The group is therefore checked
against its own algebraic structure rather than against stored numbers.
"""
from __future__ import annotations

from collections import Counter

import numpy as np

from cyclicnano.icosahedral import (FIVE_THREE_ANGLE_DEG, axis_pair,
                                    icosahedral_rotations, icosahedron_vertices,
                                    rotation_angle_deg, rotation_matrix)

G = icosahedral_rotations()


def test_group_has_sixty_elements():
    assert G.shape == (60, 3, 3)


def test_every_element_is_a_proper_rotation():
    for m in G:
        assert np.allclose(m @ m.T, np.eye(3), atol=1e-9)
        assert abs(np.linalg.det(m) - 1.0) < 1e-9        # +1, not a reflection


def test_elements_are_distinct():
    flat = G.reshape(60, 9)
    for i in range(60):
        for j in range(i + 1, 60):
            assert not np.allclose(flat[i], flat[j], atol=1e-8)


def test_rotation_angles_match_the_conjugacy_classes():
    """1 identity, 12 at 72, 20 at 120, 12 at 144, 15 at 180.

    This is the strongest single check available: a wrong axis or a failed closure
    changes the histogram immediately.
    """
    assert Counter(round(rotation_angle_deg(m)) for m in G) == {0: 1, 72: 12, 120: 20, 144: 12, 180: 15}


def test_group_is_closed_under_multiplication():
    flat = G.reshape(60, 9)
    for a in G:
        for b in G:
            product = (a @ b).ravel()
            assert np.isclose(np.abs(flat - product).max(axis=1), 0, atol=1e-8).any()


def test_every_element_has_an_inverse_in_the_group():
    for a in G:
        assert any(np.allclose(a @ b, np.eye(3), atol=1e-8) for b in G)


def test_vertices_are_twelve_unit_vectors_in_opposite_pairs():
    v = icosahedron_vertices()
    assert v.shape == (12, 3)
    assert np.allclose(np.linalg.norm(v, axis=1), 1.0)
    for x in v:                                           # every vertex has an antipode
        assert any(np.allclose(x, -y, atol=1e-9) for y in v)


def test_the_two_axes_meet_at_the_icosahedral_angle():
    five, three = axis_pair()
    angle = np.degrees(np.arccos(abs(float(np.dot(five, three)))))
    assert abs(angle - FIVE_THREE_ANGLE_DEG) < 1e-9
    assert abs(FIVE_THREE_ANGLE_DEG - 37.3774) < 1e-3     # the value quoted in the literature


def test_the_axes_have_the_orders_they_claim():
    five, three = axis_pair()
    for axis, order in ((five, 5), (three, 3)):
        r = rotation_matrix(axis, 2 * np.pi / order)
        assert any(np.allclose(r, m, atol=1e-8) for m in G)          # the rotation is in the group
        assert np.allclose(np.linalg.matrix_power(r, order), np.eye(3), atol=1e-9)
        for k in range(1, order):                                     # and no smaller power is identity
            assert not np.allclose(np.linalg.matrix_power(r, k), np.eye(3), atol=1e-6)


def test_sixty_copies_of_one_subunit_make_twelve_pentamers():
    """Why 60 rotations times one chain gives 12 x 5 and not something else.

    The stabiliser of the fivefold axis is a C5 subgroup, so the 60 elements fall
    into 12 cosets of 5. Each coset carries the pentamer to one vertex, and within
    a coset the five elements fill that pentamer.
    """
    five, _ = axis_pair()
    stabiliser = [m for m in G if np.allclose(m @ five, five, atol=1e-8)]
    assert len(stabiliser) == 5
    assert len(G) // len(stabiliser) == 12


def test_sixty_copies_of_one_subunit_make_twenty_trimers():
    _, three = axis_pair()
    stabiliser = [m for m in G if np.allclose(m @ three, three, atol=1e-8)]
    assert len(stabiliser) == 3
    assert len(G) // len(stabiliser) == 20
