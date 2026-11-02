"""Placing a predicted monomer at the designed subunit positions.

One operation with two callers: stage 05 places two copies to get an interface
Rosetta can score, and the web page places one per subunit so the ring can be
looked at. What is checked here is that the placement is a rigid motion onto the
design chains, that every atom of the prediction comes along and not only the
four the metric layer reads, and that the labels come out in ring order.

No Rosetta and no AlphaFold: a synthetic ring, and a subunit taken out of it, so
the right answer is known exactly.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

import synthetic
from cyclicnano.assembly import place_monomer
from cyclicnano.pdbio import extract_monomer, read_pdb, write_pdb
from cyclicnano.s05_interface import build_dimer


def _ring(tmp: Path, n_sym: int = 5, twist: float = 0.0):
    """A design backbone and a prediction of its subunit, as two files.

    `twist` rotates the prediction out of the design's frame, which is what an
    AlphaFold prediction always is: it comes back in its own coordinates and has
    to be superposed before it means anything here.
    """
    design = synthetic.cyclic_oligomer(n_sym=n_sym, per_helix=12, n_helices=2)
    mono = extract_monomer(design)
    if twist:
        rot = synthetic._rotz(twist)
        for chain in mono.chains.values():
            chain.coords["CA"] = chain.coords["CA"] @ rot.T + np.array([30.0, -5.0, 2.0])

    design_pdb, mono_pdb = tmp / "design_bb.pdb", tmp / "pred.pdb"
    write_pdb(str(design_pdb), design)
    write_pdb(str(mono_pdb), mono)
    return design, design_pdb, mono_pdb


def _with_side_chains(path: Path) -> int:
    """Add a CB beside every CA, and say how many atoms the file now has.

    The prediction is the only structure in this pipeline carrying side chains,
    and carrying them is the entire reason the placement reads the file as text
    instead of going through `read_pdb`, which keeps four atoms per residue.
    """
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    out = []
    for line in lines:
        out.append(line)
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            x = float(line[30:38]) + 1.5
            out.append(f"{line[:12]} CB {line[16:30]}{x:8.3f}{line[38:]}")
    path.write_text("".join(out), encoding="utf-8")
    return sum(1 for line in out if line.startswith("ATOM"))


def test_every_subunit_of_the_design_gets_a_copy():
    """Five chains out for a C5 design, labelled as the design labels them.

    Ring order matters downstream: chains A and B are neighbours in an
    RFdiffusion oligomer and the dimer relies on it, so the labels are carried
    over rather than reassigned.
    """
    with tempfile.TemporaryDirectory() as tmp:
        design, design_pdb, mono_pdb = _ring(Path(tmp), n_sym=5, twist=37.0)
        built = place_monomer(design_pdb, mono_pdb)

        assert built["chains"] == ["A", "B", "C", "D", "E"]
        out = read_pdb_text(built["pdb"], Path(tmp))
        assert out.chain_ids == design.chain_ids
        for cid in design.chain_ids:
            assert out.chains[cid].n_res == design.chains[cid].n_res


def test_the_copies_land_on_the_design_chains():
    """The fit is the whole point, so it is measured against the design itself.

    The prediction here IS a design subunit, moved somewhere else, so a correct
    placement puts it back exactly and `fit_rmsd` is zero to rounding. A real
    prediction differs from the design and the number says by how much.
    """
    with tempfile.TemporaryDirectory() as tmp:
        design, design_pdb, mono_pdb = _ring(Path(tmp), n_sym=5, twist=119.0)
        built = place_monomer(design_pdb, mono_pdb)

        assert built["fit_rmsd"] < 0.01
        out = read_pdb_text(built["pdb"], Path(tmp))
        for cid in design.chain_ids:
            moved = out.chains[cid].ca
            target = design.chains[cid].ca
            assert np.abs(moved - target).max() < 0.01, cid


def test_the_side_chains_come_along():
    """Not only the atoms `read_pdb` keeps.

    Rosetta cannot score an interface made of glycine main chain, and a ring
    drawn without side chains is not what anyone wants to look at either. The
    placement reads the prediction as text for exactly this reason.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _, design_pdb, mono_pdb = _ring(Path(tmp), n_sym=3)
        n_atoms = _with_side_chains(mono_pdb)

        built = place_monomer(design_pdb, mono_pdb)
        placed = [line for line in built["pdb"].splitlines() if line.startswith("ATOM")]
        assert len(placed) == 3 * n_atoms
        # Half the prediction's atoms are CB, and there are three copies of it.
        assert sum(1 for line in placed if line[12:16].strip() == "CB") == 3 * (n_atoms // 2)
        # Every chain id in the output is one the design asked for, side chains
        # included: the relabelling is applied to the line, not to an atom list.
        assert {line[21] for line in placed} == {"A", "B", "C"}


def test_a_chain_the_design_does_not_have_is_refused():
    """A pairing that cannot be honoured is an error and not an empty ring.

    The web route turns this into a 403 rather than serving half an assembly.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _, design_pdb, mono_pdb = _ring(Path(tmp), n_sym=3)
        for bad in ([("A", "A"), ("Z", "B")], []):
            try:
                place_monomer(design_pdb, mono_pdb, bad)
            except ValueError:
                continue
            raise AssertionError(f"expected ValueError for {bad!r}")


def test_the_dimer_is_the_same_placement_on_two_chains():
    """Stage 05 and the page must not disagree about where a subunit sits.

    They share the placement precisely so that the energy and the picture
    describe one geometry; this is what says so.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _, design_pdb, mono_pdb = _ring(tmp, n_sym=5, twist=64.0)

        out = tmp / "dimer.pdb"
        info = build_dimer(design_pdb, mono_pdb, out)
        pair = place_monomer(design_pdb, mono_pdb, [("A", "A"), ("B", "B")])

        assert info["dimer_fit_rmsd"] == pair["fit_rmsd"]
        assert out.read_text(encoding="utf-8") == pair["pdb"]
        dimer = read_pdb_text(pair["pdb"], tmp)
        assert dimer.chain_ids == ["A", "B"]


def read_pdb_text(text: str, tmp: Path):
    """`read_pdb` takes a path, and what is under test returns text."""
    path = tmp / "_read_back.pdb"
    path.write_text(text, encoding="utf-8")
    return read_pdb(str(path))
