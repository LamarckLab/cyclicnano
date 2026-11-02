"""The web UI's reading layer.

The HTTP layer is a thin wrapper and FastAPI is not a dependency of the suite,
so what is tested here is `store`: which directories count as runs, where the
funnel numbers come from, how a backbone is found when its folder carries a rank
prefix, and that a path arriving from a browser cannot leave the run.

No pytest fixtures, so that `tests/run_all.py` can run this module too: it skips
any test that takes an argument.
"""
from __future__ import annotations

import contextlib
import json
import shutil
import tempfile
from pathlib import Path

from cyclicnano.pdbio import extract_monomer, write_pdb
from cyclicnano.webui import store

import synthetic


def _tree(root: Path, run: str = "c5", *, tables: bool = True,
          interface: bool = True, done: bool = True) -> Path:
    """One run laid out the way stage 04 lays it out.

    Two backbones: `design_27` passed and so sits under a rank prefix, while
    `design_14` did not and has no interface folder, because stage 05 only
    evaluates the sequences of passing backbones.
    """
    base = root / run / store.DELIVERY
    (root / run).mkdir(parents=True, exist_ok=True)
    if done:
        (root / run / ".done").touch()

    # A real C5 backbone and real predictions: the assembled ring is built from
    # both at request time, and a one-line stub cannot be assembled.
    ring = synthetic.cyclic_oligomer(n_sym=5, per_helix=10, n_helices=2)
    monomer = extract_monomer(ring)

    folders = {"design_27": "passed/1_design_27", "design_14": "below_cutoff/design_14"}
    for bid, folder in folders.items():
        passing = folder.startswith("passed")
        d = base / folder
        (d / "predictions").mkdir(parents=True, exist_ok=True)
        write_pdb(str(d / f"{bid}_bb.pdb"), ring)

        sequences = []
        for k in range(2):
            write_pdb(str(d / "predictions" / f"seq{k:03d}_af2.pdb"), monomer)
            block = None
            if passing and interface and k == 0:
                (d / "interface").mkdir(exist_ok=True)
                (d / "interface" / f"seq{k:03d}_dimer_relaxed.pdb").write_text(
                    "ATOM\n", encoding="utf-8")
                block = {"dG_interface": -31.4, "interface_pass": True,
                         "dimer": f"interface/seq{k:03d}_dimer_relaxed.pdb"}
            sequences.append({
                "sequence_id": f"{bid}__seq{k:03d}", "best": k == 0, "seq_index": k,
                "sc_pass": k == 0, "sc_rmsd": 0.46 if k == 0 else 3.3,
                "mean_plddt": 95.1 if k == 0 else 74.0,
                "prediction": f"predictions/seq{k:03d}_af2.pdb",
                "scores": f"predictions/seq{k:03d}_af2_scores.json",
                "interface": block, "sequence": "AC" * 30,
            })
        (d / store.SEQUENCES).write_text(json.dumps({
            "backbone_id": bid, "n_sequences": 2, "n_pass": 1, "success_rate": 0.5,
            "best_sequence_id": f"{bid}__seq000",
            "criteria": {"rmsd_cutoff": 1.0, "plddt_cutoff": 80.0,
                         "ddg_cutoff": -25.0},
            "sequences": sequences,
        }), encoding="utf-8")

    (base / store.MANIFEST).write_text(json.dumps({
        "n_backbones": 2, "n_final_pass": 1,
        "backbones": [
            {"backbone_id": "design_14", "folder": folders["design_14"], "symmetry": "C5",
             "n_res_monomer": 60, "final_pass": False, "success_rate": 0.1,
             "n_interface_scored": 0, "n_interface_pass": 0,
             "interface_success_rate": None},
            {"backbone_id": "design_27", "folder": folders["design_27"], "symmetry": "C5",
             "n_res_monomer": 60, "final_pass": True, "success_rate": 0.5,
             "n_interface_scored": 1, "n_interface_pass": 1,
             "interface_success_rate": 1.0, "best_dG_interface": -31.4},
        ],
    }), encoding="utf-8")

    if tables:
        t = root / run / "tables"
        t.mkdir(parents=True, exist_ok=True)
        rows = ["backbone_id,symmetry,bb_pass"]
        rows += [f"design_{i},C5,{i < 2}" for i in range(5)]
        (t / "backbones.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        rows = ["sequence_id,fold_ok,sc_pass,interface_scored,interface_pass"]
        rows += ["design_27__seq000,True,True,True,True",
                 "design_27__seq001,True,False,,", "design_14__seq000,False,False,,"]
        (t / "sequences.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    return root


@contextlib.contextmanager
def outputs(**kwargs):
    """A throwaway output root holding one run called `c5`."""
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_webui_"))
    try:
        yield _tree(tmp, **kwargs)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def raises(exc, fn, *args):
    want = " or ".join(e.__name__ for e in (exc if isinstance(exc, tuple) else (exc,)))
    try:
        fn(*args)
    except exc:
        return
    except Exception as other:
        raise AssertionError(f"expected {want}, got {other!r}") from other
    raise AssertionError(f"expected {want} from {fn.__name__}{args}")


# ----------------------------------------------------------------------- runs

def test_only_real_runs_are_listed():
    """A run is a directory with a delivery manifest, and an underscore directory
    is never one.

    The output root also holds `_archive`, `_experiments` and pilot directories,
    some of which do contain a delivery manifest from an earlier attempt.
    Listing those beside the real runs would be the most misleading thing this
    page could do, so the name rule lives here and is not left to the caller.
    """
    with outputs() as root:
        _tree(root, "c3")
        (root / "_iptm_pilot").mkdir()
        _tree(root, "_archive")                         # a complete old run
        (root / "notes").mkdir()                        # a directory, but not a run
        (root / "c5.log").write_text("", encoding="utf-8")

        assert sorted(r["run"] for r in store.list_runs(root)) == ["c3", "c5"]
        raises(store.NotFound, store.run_dir, root, "_archive")
        raises(store.NotFound, store.run_dir, root, "notes")


def test_a_browser_cannot_name_a_run_outside_the_root():
    with outputs() as root:
        for bad in ("..", "../c5", "/etc", "_archive", "", "c5/../.."):
            raises((store.Forbidden, store.NotFound), store.run_dir, root, bad)


def test_the_funnel_comes_from_the_tables_not_the_delivery():
    """The dropped backbones are the whole point of a funnel.

    The delivery folder holds only what survived stage 01, so counting it would
    report 2 generated backbones instead of 5 and make the filter look free.
    """
    with outputs() as root:
        f = store.funnel(root, "c5")
        assert f["n_generated"] == 5                    # tables/backbones.csv
        assert f["n_passed_filter"] == 2                # bb_pass
        assert f["n_sequences"] == 3
        assert f["n_predictions"] == 2                  # fold_ok
        assert f["n_sequences_pass"] == 1               # sc_pass
        assert f["n_delivered"] == 2                    # delivery.json
        assert f["n_final_pass"] == 1
        assert f["n_interface"] == 1
        assert f["n_interface_pass"] == 1
        assert f["symmetry"] == "C5"


def test_a_stage_that_has_not_run_is_unknown_rather_than_zero():
    """None and 0 are different claims, and the page renders them differently.

    A run whose tables are missing has not been counted; saying "0 generated"
    about a run that delivered two backbones would be a visible contradiction.
    """
    with outputs(tables=False) as root:
        f = store.funnel(root, "c5")
        assert f["n_generated"] is None
        assert f["n_passed_filter"] is None
        assert f["n_interface"] is None
        assert f["n_final_pass"] == 1                   # the manifest is still there


# ------------------------------------------------------------------ backbones

def test_passing_backbones_come_first_and_rank_sorts_as_a_number():
    """`10_design_169` must not sort before `1_design_222`.

    Stage 04 prefixes passing folders with their delivery rank, so the obvious
    string sort puts ranks 10 to 19 ahead of rank 2, silently reordering the
    table that gets read top-down.
    """
    with outputs() as root:
        base = root / "c5" / store.DELIVERY
        for rank, bid in ((2, "design_a"), (10, "design_b")):
            folder = f"passed/{rank}_{bid}"
            (base / folder).mkdir(parents=True)
            (base / folder / store.SEQUENCES).write_text(
                json.dumps({"backbone_id": bid, "sequences": []}), encoding="utf-8")
        manifest = json.loads((base / store.MANIFEST).read_text(encoding="utf-8"))
        manifest["backbones"] += [
            {"backbone_id": "design_a", "folder": "passed/2_design_a", "final_pass": True},
            {"backbone_id": "design_b", "folder": "passed/10_design_b", "final_pass": True},
        ]
        (base / store.MANIFEST).write_text(json.dumps(manifest), encoding="utf-8")

        rows = store.backbones(root, "c5")
        assert [r["backbone_id"] for r in rows] == [
            "design_27", "design_a", "design_b", "design_14"]
        assert [r["group"] for r in rows] == ["passed", "passed", "passed", "below_cutoff"]
        assert rows[0]["has_interface"] is True
        assert rows[-1]["has_interface"] is False


def test_a_backbone_is_found_by_id_whatever_its_rank_prefix():
    """The caller knows `design_27`, not `passed/1_design_27`.

    The prefix moves whenever a run is re-ranked, so the id is the only stable
    handle, and a URL bookmarked last week has to keep working.
    """
    with outputs() as root:
        assert store.backbone(root, "c5", "design_27")["folder"] == "passed/1_design_27"
        assert store.backbone(root, "c5", "design_14")["folder"] == "below_cutoff/design_14"
        raises(store.NotFound, store.backbone, root, "c5", "design_999")
        raises(store.Forbidden, store.backbone, root, "c5", "../../etc")


def test_interface_metrics_are_lifted_out_of_their_block():
    """The energy and its verdict sit inside `interface`, null for a failing sequence.

    The table needs them as plain columns, and a sequence that was never evaluated
    has to come out as None rather than raise or read as zero. `interface_pass` in
    particular: False means "measured and above the cutoff", None means "not
    measured", and those are different facts.
    """
    with outputs() as root:
        detail = store.backbone(root, "c5", "design_27")
        good, bad = detail["sequences"]
        assert good["dG_interface"] == -31.4
        assert good["interface_pass"] is True
        assert bad["interface"] is None
        assert bad["dG_interface"] is None and bad["interface_pass"] is None
        assert detail["criteria"]["rmsd_cutoff"] == 1.0
        assert detail["criteria"]["ddg_cutoff"] == -25.0
        assert detail["row"]["n_interface_pass"] == 1
        assert detail["row"]["interface_success_rate"] == 1.0


def test_the_viewable_structures_are_the_ones_the_pipeline_named():
    """Four kinds: the design backbone, and per sequence a ring, a monomer and a dimer.

    The file-backed ones are read out of `sequences.json` rather than globbed, because
    that file is the pipeline's statement of which structure belongs to which sequence.
    A sequence with no interface block contributes no interface structure, and the
    oligomer *prediction* is gone along with ipTM - the ring listed here is assembled
    from the monomer on request and is not a prediction of anything.
    """
    with outputs() as root:
        detail = store.backbone(root, "c5", "design_27")
        assert [(s["key"], s["sequence_id"]) for s in detail["structures"]] == [
            ("backbone", None),
            ("oligomer", "design_27__seq000"),
            ("monomer", "design_27__seq000"),
            ("dimer", "design_27__seq000"),
            ("oligomer", "design_27__seq001"),
            ("monomer", "design_27__seq001"),
        ]
        assert not any(s["key"] == "multimer" for s in detail["structures"])
        # Which of them are in the designed coordinate frame, which is what the note
        # under the viewer has to get right.
        assert {s["key"]: s["design_frame"] for s in detail["structures"]} == {
            "backbone": True, "oligomer": True, "monomer": False, "dimer": True}
        assert {s["key"]: s["n_chains"] for s in detail["structures"]} == {
            "backbone": 5, "oligomer": 5, "monomer": 1, "dimer": 2}
        assert all(not s["missing"] for s in detail["structures"])


def test_the_assembled_ring_names_the_two_files_it_is_built_from():
    """It has no file of its own, so `source` is what the route is given.

    `path` is a name and nothing else: it identifies the selection and titles a
    download, and resolving it as a file would fail. Everything with a file behind
    it says `source: None`, so one page-side rule tells the two apart.
    """
    with outputs() as root:
        detail = store.backbone(root, "c5", "design_27")
        ring = next(s for s in detail["structures"] if s["key"] == "oligomer")
        assert ring["path"] == "passed/1_design_27/seq000_af2_C5.pdb"
        assert not (root / "c5" / store.DELIVERY / ring["path"]).exists()
        assert ring["source"] == "passed/1_design_27/predictions/seq000_af2.pdb"
        assert ring["design"] == "passed/1_design_27/design_27_bb.pdb"
        for named in (ring["source"], ring["design"]):
            assert store.resolve(root, "c5", named).is_file()
        assert all(s["source"] is None
                   for s in detail["structures"] if s["key"] != "oligomer")


def test_a_single_chain_design_has_no_ring_to_assemble():
    """A one-chain backbone is one subunit, and one subunit is the monomer already.

    The guard is the chain count of the file, not the declared symmetry: the count
    is what an assembly would actually be built from.
    """
    with outputs() as root:
        folder = root / "c5" / store.DELIVERY / "passed" / "1_design_27"
        write_pdb(str(folder / "design_27_bb.pdb"),
                  extract_monomer(synthetic.cyclic_oligomer(n_sym=5, per_helix=10,
                                                            n_helices=2)))
        detail = store.backbone(root, "c5", "design_27")
        assert not any(s["key"] == "oligomer" for s in detail["structures"])
        assert [s["key"] for s in detail["structures"]][0] == "backbone"


def test_a_path_the_manifest_names_but_disk_lacks_is_reported_missing():
    """Quietly dropping it would hide a truncated run.

    A prediction can be named by `sequences.json` and absent from disk if a copy
    was interrupted; the page says so instead of showing a shorter list.
    """
    with outputs() as root:
        (root / "c5" / store.DELIVERY / "passed" / "1_design_27"
         / "predictions" / "seq001_af2.pdb").unlink()
        detail = store.backbone(root, "c5", "design_27")
        missing = [s for s in detail["structures"] if s["missing"]]
        assert [s["sequence_id"] for s in missing] == ["design_27__seq001"]


# ----------------------------------------------------------------------- files

def test_a_file_request_cannot_leave_the_delivery_folder():
    """The path comes from a query string, so this is the one real attack surface.

    Checked after resolution rather than by inspecting the string, so that a
    symlink pointing out of the tree is caught as well.
    """
    with outputs() as root:
        assert store.resolve(root, "c5", "passed/1_design_27/design_27_bb.pdb").is_file()

        (root / "secret.pdb").write_text("ATOM\n", encoding="utf-8")
        for bad in ("../../secret.pdb", "../../../etc/passwd", "/etc/passwd",
                    "passed/../../../secret.pdb"):
            raises(store.Forbidden, store.resolve, root, "c5", bad)


def test_only_structures_and_scores_are_served():
    """A read-only viewer, not a file server.

    The delivery tree sits next to logs and configs on the same machine, and a
    suffix allowlist is what keeps one careless route from exposing them.
    """
    with outputs() as root:
        base = root / "c5" / store.DELIVERY
        (base / "notes.txt").write_text("hi", encoding="utf-8")
        (base / "script.sh").write_text("hi", encoding="utf-8")
        for bad in ("notes.txt", "script.sh"):
            raises(store.Forbidden, store.resolve, root, "c5", bad)

        assert store.resolve(root, "c5", store.MANIFEST).name == store.MANIFEST
        raises(store.NotFound, store.resolve, root, "c5", "passed/1_design_27/nope.pdb")


# ---------------------------------------------------------------------- export

def test_the_export_scopes_answer_three_different_questions():
    """Narrow by default, because the default is the one people act on.

    The fixture has one delivered backbone and one below the cutoff, each with
    one passing and one failing sequence. `delivered` is what you would order:
    a passing sequence on a passing backbone. `passed` adds the sequence that
    folded back on a backbone whose overall rate did not, which is where a lucky
    single design hides. `all` is the record, failures included.
    """
    with outputs() as root:
        counts = {scope: store.export_rows(root, "c5", scope) for scope in store.SCOPES}
        assert [r["sequence_id"] for r in counts["delivered"]] == ["design_27__seq000"]
        assert sorted(r["sequence_id"] for r in counts["passed"]) == [
            "design_14__seq000", "design_27__seq000"]
        assert len(counts["all"]) == 4
        # Ranked backbones keep the delivery order, so the sheet reads top-down
        # the same way the table does.
        assert counts["delivered"][0]["rank"] == 1
        assert counts["passed"][1]["rank"] is None        # below_cutoff has no rank
        raises(store.Forbidden, store.export_rows, root, "c5", "everything")


def test_a_fasta_record_starts_with_a_marker_and_carries_its_metrics():
    """A FASTA that gets mailed to a vendor arrives without this page.

    So the numbers that tell two records apart ride in the description line, and
    the record marker has to be there or the file is not a FASTA at all.
    """
    with outputs() as root:
        rows = store.export_rows(root, "c5", "delivered")
        text = store.export_fasta(rows, "C5")
        lines = text.splitlines()

        assert len(lines) == 2 * len(rows)
        assert all(line.startswith(">") for line in lines[0::2])
        assert all(not line.startswith(">") for line in lines[1::2])
        assert text.endswith("\n")

        header = lines[0]
        assert header.startswith(">design_27__seq000 C5 ")
        for key in ("sc_rmsd=", "mean_plddt=", "dG_interface=", "len=60"):
            assert key in header, header
        assert "iptm" not in header
        assert lines[1] == "AC" * 30

        assert store.export_fasta([], "C5") == ""        # no stray newline


def test_the_csv_opens_in_excel_and_names_the_stored_field():
    """The byte-order mark is the whole reason a spreadsheet was asked for.

    Without it Excel reads a UTF-8 file as the system code page on a
    double-click. The energy column also keeps the name the pipeline stores it
    under rather than the label the table shows, so a file that outlives the
    page is still traceable to the field it came from.
    """
    with outputs() as root:
        rows = store.export_rows(root, "c5", "all")
        text = store.export_csv(rows)

        assert text.startswith("\ufeff")
        header = text.lstrip("\ufeff").splitlines()[0]
        assert header.split(",") == store.EXPORT_FIELDS
        assert "dG_interface" in header and "ddG" not in header
        assert len(text.splitlines()) == len(rows) + 1

        # A sequence that was never evaluated leaves the cell empty, not zero.
        body = [line.split(",") for line in text.lstrip("\ufeff").splitlines()[1:]]
        unevaluated = [r for r in body if r[2].endswith("seq001")]
        assert unevaluated and all(r[8] == "" and r[9] == "" for r in unevaluated)


def test_a_run_stopped_after_the_filter_is_still_a_run():
    """Stage 01 is the cheap look, and it has no delivery manifest.

    Running 00 and 01 alone is how you find out whether the filter leaves
    anything before spending GPU hours on ProteinMPNN and AlphaFold. Requiring
    the manifest to recognise a run made exactly that run invisible on the page
    that would have shown its funnel, which is the only thing it has to show.
    """
    with outputs() as root:
        early = root / "c9_early"
        (early / "tables").mkdir(parents=True)
        rows = ["backbone_id,symmetry,bb_pass"]
        rows += [f"design_{i},C9,{i < 1}" for i in range(4)]
        (early / "tables" / "backbones.csv").write_text(
            "\n".join(rows) + "\n", encoding="utf-8")

        listed = {r["run"]: r for r in store.list_runs(root)}
        assert set(listed) == {"c5", "c9_early"}
        assert listed["c9_early"]["delivered"] is False
        assert listed["c5"]["delivered"] is True

        # The funnel is there; the delivery numbers are unknown, not zero.
        assert listed["c9_early"]["n_generated"] == 4
        assert listed["c9_early"]["n_passed_filter"] == 1
        assert listed["c9_early"]["n_final_pass"] is None

        # And the delivery table is empty rather than an exception.
        assert store.backbones(root, "c9_early") == []


def test_the_table_carries_the_whole_distribution_not_a_best_value():
    """A single best RMSD says almost nothing about twenty sequences.

    The table draws a histogram instead, so it needs every value. Raw numbers are
    sent rather than counts: the browser holds the cutoffs and bins against them, and
    sending bins would put the bin edges in one place and the thresholds they are read
    against in another. They come from the same `sequences.json` the detail page
    shows, so the two cannot disagree.
    """
    with outputs() as root:
        rows = {r["backbone_id"]: r for r in store.backbones(root, "c5")}
        dist = rows["design_27"]["distribution"]

        assert set(dist) == set(store.DISTRIBUTION_FIELDS) | {"sc_pass"}
        assert dist["sc_rmsd"] == [0.46, 3.3]
        assert dist["mean_plddt"] == [95.1, 74.0]
        assert dist["sc_pass"] == [True, False]
        # The unevaluated sequence holds the place rather than being dropped, so the
        # three arrays stay aligned with each other and with the sequence order.
        assert dist["dG_interface"] == [-31.4, None]
        assert all(len(v) == 2 for v in dist.values())

        # A backbone whose sequences were never evaluated still has the monomer arrays.
        below = rows["design_14"]["distribution"]
        assert below["sc_rmsd"] == [0.46, 3.3]
        assert below["dG_interface"] == [None, None]
