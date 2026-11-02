"""Reading a finished run for the web UI.

Everything here is standard library over `04_delivery/`, which stage 04 already
made self-contained: each backbone folder holds its design backbone, every
prediction, the interface structures, and a `sequences.json` that names them by
relative path. The browser therefore never needs to reach into 00-03 or 05, and
this module needs neither a web framework nor pandas.

Nothing is cached. A run directory is a few hundred kilobytes of CSV and JSON,
and a stage that has just finished should show up on the next page load rather
than whenever a cache decides to expire.
"""
from __future__ import annotations

import csv
import io
import json
import re
from pathlib import Path

DELIVERY = "04_delivery"
MANIFEST = "delivery.json"
SEQUENCES = "sequences.json"
PASSED = "passed"
BELOW = "below_cutoff"

# A run directory is named by the user, so it is also a path component arriving
# from the browser. Only names matching this reach the filesystem, and they are
# checked against the discovered runs on top of that.
#
# A leading underscore is allowed here on purpose, even though those directories
# are not served: they are hidden rather than illegal, so `run_dir` can answer
# "no such run" for them instead of "forbidden", and not confirm to a browser
# that `_archive` exists.
RUN_NAME = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")

# What a browser may ask for out of a run: structures and scores, nothing else.
# The point is to keep a read-only viewer from becoming a file server.
SERVABLE = {".pdb", ".json", ".csv"}


class NotFound(Exception):
    """A run, backbone or file that is not there."""


class Forbidden(Exception):
    """A path that resolves outside the run, or a file type we do not serve."""


# ----------------------------------------------------------------- small helpers

def _truthy(value: object) -> bool:
    """Whether a CSV cell written by pandas means True."""
    return str(value).strip().lower() in {"true", "1", "1.0", "yes"}


def _rows(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _count(rows: list[dict], column: str) -> int | None:
    """How many rows hold a true value in a column, or None if there is no column."""
    if not rows or column not in rows[0]:
        return None
    return sum(1 for r in rows if _truthy(r.get(column)))


def _rank_of(folder: str) -> tuple[int, str]:
    """Sort key for a delivery folder name.

    Stage 04 prefixes passing backbones with their rank, so as strings
    `10_design_169` sorts before `1_design_222`. The rank is the number.
    """
    name = folder.rsplit("/", 1)[-1]
    head = name.split("_", 1)[0]
    return (int(head), name) if head.isdigit() else (10**6, name)


# ----------------------------------------------------------------------- runs

TABLES = "tables"
BACKBONE_TABLE = "backbones.csv"


PROJECT = "project.json"


def is_run(path: Path) -> bool:
    """Whether a directory is a project.

    Any one of three markers is enough, because a project is visible from the
    moment it is named and gains the others as stages run.

    `project.json` alone is a project with nothing computed yet, which is a
    useful state: it says what is intended before anything has been spent on it.
    `tables/backbones.csv` alone is a run stopped after stage 01, the cheap way
    to see how many backbones survive the filter — requiring the delivery
    manifest made exactly that run invisible on the page that would show its
    funnel. The manifest alone covers runs that predate both files.
    """
    return any((path / rel).is_file() for rel in
               (PROJECT, f"{TABLES}/{BACKBONE_TABLE}", f"{DELIVERY}/{MANIFEST}"))


def list_runs(root: Path | str) -> list[dict]:
    """Every project under the output root, most recently touched first.

    Directories whose name begins with an underscore are skipped. That is where
    pilots, archived attempts and the job records live, and showing those beside
    the real projects would be the most misleading thing this page could do.
    """
    from . import projects                 # projects reads store; imported here

    root = Path(root)
    if not root.is_dir():
        raise NotFound(f"output root does not exist: {root}")

    runs = []
    for child in sorted(root.iterdir()):
        if not child.is_dir() or child.name.startswith("_") or not is_run(child):
            continue
        # Whichever record this project actually has: a run stopped after stage
        # 01 has only the table, and one just created has only its manifest.
        marker = next((p for p in (child / DELIVERY / MANIFEST,
                                   child / TABLES / BACKBONE_TABLE,
                                   child / PROJECT) if p.is_file()),
                      child)
        runs.append(_describe(root, child, marker))
    return sorted(runs, key=lambda r: r["updated"], reverse=True)


def _describe(root: Path, folder: Path, marker: Path) -> dict:
    """One project, as both the list and its own page need it.

    Assembled in one place so the card and the page cannot disagree about a
    project's symmetry or whether it has started.
    """
    from . import projects

    measured = funnel(root, folder.name)
    meta = projects.read(root, folder.name)
    return {
        "run": folder.name,
        "updated": int(marker.stat().st_mtime),
        "delivered": (folder / DELIVERY / MANIFEST).is_file(),
        "complete": (folder / ".done").exists(),
        # False until a stage has run, which is how the page knows to offer the
        # submit form in place of a funnel.
        "started": bool(measured.get("n_generated")),
        **meta,
        **measured,
        # The measured symmetry wins where there is one; the declaration is all
        # an empty project has. Both are kept, so a disagreement between what
        # was declared and what the tables say stays visible.
        "symmetry": measured.get("symmetry") or meta.get("symmetry_declared") or "",
    }


def describe(root: Path | str, run: str) -> dict:
    """One project by name, having checked the browser may ask for it."""
    folder = run_dir(root, run)
    marker = next((p for p in (folder / DELIVERY / MANIFEST,
                               folder / TABLES / BACKBONE_TABLE,
                               folder / PROJECT) if p.is_file()), folder)
    return _describe(Path(root), folder, marker)


def run_dir(root: Path | str, run: str) -> Path:
    """The directory of one run, having checked that the browser may name it."""
    if not RUN_NAME.match(run or ""):
        raise Forbidden(f"not a run name: {run!r}")
    path = Path(root) / run
    if run.startswith("_") or not path.is_dir() or not is_run(path):
        raise NotFound(f"no such run: {run}")
    return path


def funnel(root: Path | str, run: str) -> dict:
    """The numbers at each stage of one run.

    Read from the run tables rather than recounted from the delivery folder: the
    tables are the pipeline's own record and they still hold the backbones and
    sequences that were dropped, which is the entire interest of a funnel. A
    stage that has not run yields None, not zero, so the page can tell "nothing
    passed" apart from "not attempted".
    """
    path = Path(root) / run
    backbone_rows = _rows(path / "tables" / "backbones.csv")
    sequence_rows = _rows(path / "tables" / "sequences.csv")

    manifest_path = path / DELIVERY / MANIFEST
    manifest = {}
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    symmetry = ""
    for row in backbone_rows:
        if row.get("symmetry"):
            symmetry = row["symmetry"]
            break

    return {
        "symmetry": symmetry,
        "n_generated": len(backbone_rows) or None,
        "n_passed_filter": _count(backbone_rows, "bb_pass"),
        "n_sequences": len(sequence_rows) or None,
        "n_predictions": _count(sequence_rows, "fold_ok"),
        "n_sequences_pass": _count(sequence_rows, "sc_pass"),
        "n_delivered": manifest.get("n_backbones"),
        "n_final_pass": manifest.get("n_final_pass"),
        "n_interface": _count(sequence_rows, "interface_scored"),
        "n_interface_pass": _count(sequence_rows, "interface_pass"),
    }


# ------------------------------------------------------------------ backbones

def backbones(root: Path | str, run: str) -> list[dict]:
    """The delivery table of one run, passing backbones first.

    Each row is `delivery.json` as stage 04 wrote it, plus the group it landed
    in and whether its interface structures are on disk, so that the table can
    be rendered without one request per row.
    """
    base = run_dir(root, run) / DELIVERY
    if not (base / MANIFEST).is_file():
        return []                                   # a run that has not reached 04
    manifest = json.loads((base / MANIFEST).read_text(encoding="utf-8"))

    rows = []
    for entry in manifest.get("backbones", []):
        folder = str(entry.get("folder", ""))
        row = dict(entry)
        row["group"] = PASSED if folder.startswith(PASSED) else BELOW
        row["has_interface"] = (base / folder / "interface").is_dir()
        row["distribution"] = _distribution(base / folder / SEQUENCES)
        rows.append(row)

    rows.sort(key=lambda r: (r["group"] != PASSED, _rank_of(str(r.get("folder", "")))))
    return rows


# The columns a histogram is drawn from. Raw values rather than counts, because the
# browser holds the cutoffs and bins against them; sending bins would put the bin
# edges in one place and the thresholds they are read against in another.
DISTRIBUTION_FIELDS = ("sc_rmsd", "mean_plddt", "dG_interface")


def _distribution(path: Path) -> dict:
    """Every sequence's metrics for one backbone, as parallel arrays.

    Read from the same `sequences.json` the detail page shows, so the table and the
    detail cannot disagree. A backbone has tens of sequences, not thousands, so this
    is a few kilobytes and no summarising is needed before it reaches the browser.
    """
    if not path.is_file():
        return {}
    detail = json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, list] = {field: [] for field in DISTRIBUTION_FIELDS}
    out["sc_pass"] = []
    for seq in detail.get("sequences", []):
        iface = seq.get("interface") or {}
        out["sc_rmsd"].append(seq.get("sc_rmsd"))
        out["mean_plddt"].append(seq.get("mean_plddt"))
        out["dG_interface"].append(iface.get("dG_interface"))
        out["sc_pass"].append(bool(seq.get("sc_pass")))
    return out


def _find_folder(base: Path, backbone_id: str) -> str:
    """The delivery folder of one backbone, found by id rather than by rank prefix.

    The caller knows `design_27`; on disk that is `passed/1_design_27` or
    `below_cutoff/design_27`, and the rank moves whenever the run is re-ranked.
    """
    for group in (PASSED, BELOW):
        parent = base / group
        if not parent.is_dir():
            continue
        for child in sorted(parent.iterdir()):
            if child.is_dir() and (child.name == backbone_id
                                   or child.name.split("_", 1)[-1] == backbone_id):
                return f"{group}/{child.name}"
    raise NotFound(f"no such backbone: {backbone_id}")


def _chain_ids(path: Path) -> list[str]:
    """The chain ids of a PDB, in the order they first appear.

    Column 22, read straight off the atom records. RFdiffusion writes a cyclic
    oligomer in ring order, so the count is the number of subunits the design
    asks for and the order is the ring - both of which the page needs to say
    "pentamer" rather than "oligomer".
    """
    seen = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                cid = line[21]
                if cid not in seen:
                    seen.append(cid)
    return seen


def backbone(root: Path | str, run: str, backbone_id: str) -> dict:
    """One backbone: its delivery row, its sequences, and what can be viewed.

    The structures come from `sequences.json` rather than from globbing the
    folder, because that file is the pipeline's statement of which structure
    belongs to which sequence. A path it names that is missing on disk is
    reported as missing rather than quietly dropped.
    """
    if not RUN_NAME.match(backbone_id or ""):
        raise Forbidden(f"not a backbone id: {backbone_id!r}")

    base = run_dir(root, run) / DELIVERY
    folder = _find_folder(base, backbone_id)
    detail = json.loads((base / folder / SEQUENCES).read_text(encoding="utf-8"))

    row = next((r for r in backbones(root, run)
                if r.get("backbone_id") == detail.get("backbone_id")), {})

    # No display text here. A structure is described by its kind, the sequence it
    # belongs to, and whether it is in the designed coordinate frame; the label and
    # the explanation are the page's business, and the page has two languages.
    design = f"{detail['backbone_id']}_bb.pdb"
    design_rel = f"{folder}/{design}"
    subunits = (_chain_ids(base / folder / design)
                if (base / folder / design).is_file() else [])
    structures = []
    # Listed on the strength of the file being there, not of this reader having
    # understood it: a design backbone is the pipeline's own output and the
    # viewer can be left to parse it. Only the assembled ring needs the chain
    # count, so only the ring is withheld when it comes back empty.
    if (base / folder / design).is_file():
        structures.append({
            "key": "backbone",
            "path": design_rel,
            "source": None,
            "n_chains": len(subunits) or None,
            "sequence_id": None,
            "seq_short": None,
            "design_frame": True,
            "missing": False,
        })

    sequences = []
    for seq in detail.get("sequences", []):
        entry = dict(seq)
        iface = seq.get("interface") or {}
        entry["dG_interface"] = iface.get("dG_interface")
        entry["interface_pass"] = iface.get("interface_pass")
        sequences.append(entry)

        sid = str(seq.get("sequence_id", ""))
        short = sid.split("__", 1)[-1] if "__" in sid else sid
        prediction = seq.get("prediction")
        # The whole ring for this sequence, built on request from the prediction
        # and the design backbone rather than stored: it is the one structure
        # here that no stage produces, and a file per sequence per subunit would
        # be tens of megabytes of something derivable in milliseconds. The path
        # is a name for the selection and for a download, and no file is behind
        # it; `source` is what the route is actually given.
        if prediction and len(subunits) > 1 and (base / folder / prediction).is_file():
            structures.append({
                "key": "oligomer",
                "path": f"{folder}/{Path(prediction).stem}_C{len(subunits)}.pdb",
                "source": f"{folder}/{prediction}",
                "design": design_rel,
                "n_chains": len(subunits),
                "sequence_id": sid,
                "seq_short": short,
                "design_frame": True,
                "missing": False,
            })

        for rel, key, design_frame, n_chains in (
            (prediction, "monomer", False, 1),
            (iface.get("dimer"), "dimer", True, 2),
        ):
            if not rel:
                continue
            structures.append({
                "key": key,
                "path": f"{folder}/{rel}",
                "source": None,
                "n_chains": n_chains,
                "sequence_id": sid,
                "seq_short": short,
                "design_frame": design_frame,
                "missing": not (base / folder / rel).is_file(),
            })

    return {
        "run": run,
        "folder": folder,
        "group": PASSED if folder.startswith(PASSED) else BELOW,
        "row": row,
        "criteria": detail.get("criteria", {}),
        "n_sequences": detail.get("n_sequences"),
        "n_pass": detail.get("n_pass"),
        "success_rate": detail.get("success_rate"),
        "best_sequence_id": detail.get("best_sequence_id"),
        "sequences": sequences,
        "structures": structures,
    }


# ---------------------------------------------------------------------- export

# What a sequence is worth ordering on. Deliberately not every column of
# sequences.csv: this is the sheet someone reads before spending money, so it
# carries the verdict, the two monomer numbers it rests on, the two interface
# numbers that describe but do not gate, and the sequence itself.
EXPORT_FIELDS = [
    "backbone_id", "rank", "sequence_id", "seq_index", "best",
    "sc_pass", "sc_rmsd", "mean_plddt", "dG_interface", "interface_pass",
    "length", "sequence",
]

SCOPES = ("delivered", "passed", "all")


def export_rows(root: Path | str, run: str, scope: str = "delivered") -> list[dict]:
    """Every sequence worth exporting, best backbone first.

    Three scopes, because three different questions get asked of this file:

    `delivered` is the default and the narrowest - sequences that passed, from
    backbones whose success rate passed. These are the ones you would actually
    order, and ranking them by their backbone keeps the sheet in the same order
    as the delivery table.

    `passed` adds the sequences that folded back correctly on a backbone whose
    overall success rate did not, which is where a lucky single design hides.

    `all` is everything in the delivery, failures included, for the record.
    """
    if scope not in SCOPES:
        raise Forbidden(f"not an export scope: {scope!r}")

    base = run_dir(root, run) / DELIVERY
    rows = []
    for entry in backbones(root, run):
        if scope == "delivered" and not entry.get("final_pass"):
            continue
        folder = str(entry.get("folder", ""))
        detail_path = base / folder / SEQUENCES
        if not detail_path.is_file():
            continue
        detail = json.loads(detail_path.read_text(encoding="utf-8"))

        rank = _rank_of(folder)[0]
        for seq in detail.get("sequences", []):
            if scope in ("delivered", "passed") and not seq.get("sc_pass"):
                continue
            iface = seq.get("interface") or {}
            sequence = str(seq.get("sequence", ""))
            rows.append({
                "backbone_id": detail.get("backbone_id"),
                "rank": rank if rank < 10**6 else None,
                "sequence_id": seq.get("sequence_id"),
                "seq_index": seq.get("seq_index"),
                "best": bool(seq.get("best")),
                "sc_pass": bool(seq.get("sc_pass")),
                "sc_rmsd": seq.get("sc_rmsd"),
                "mean_plddt": seq.get("mean_plddt"),
                "dG_interface": iface.get("dG_interface"),
                "interface_pass": iface.get("interface_pass"),
                "length": len(sequence),
                "sequence": sequence,
            })
    return rows


def _round(value: object, digits: int) -> object:
    return round(float(value), digits) if isinstance(value, (int, float)) else value


def export_fasta(rows: list[dict], symmetry: str = "") -> str:
    """The rows as FASTA, one record per sequence.

    The metrics ride in the description line rather than in a separate file,
    because a FASTA that gets mailed to a synthesis vendor arrives on its own
    and the numbers are what let the reader tell two records apart.
    """
    out = []
    for row in rows:
        parts = [">" + str(row["sequence_id"])]
        if symmetry:
            parts.append(symmetry)
        for key, digits in (("sc_rmsd", 2), ("mean_plddt", 1),
                            ("dG_interface", 1)):
            value = row.get(key)
            if value is not None:
                parts.append(f"{key}={_round(value, digits)}")
        parts.append(f"len={row['length']}")
        if row.get("best"):
            parts.append("best")
        out.append(" ".join(parts))
        out.append(str(row["sequence"]))
    return "\n".join(out) + ("\n" if out else "")


def export_csv(rows: list[dict]) -> str:
    """The rows as CSV, with a byte-order mark.

    The BOM is there so Excel opens it as UTF-8 by double-click instead of
    mangling it, which is the whole reason a spreadsheet was asked for. Values
    are rounded here rather than in the reader: this file is read by people.
    """
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=EXPORT_FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        clean = dict(row)
        for key, digits in (("sc_rmsd", 3), ("mean_plddt", 2),
                            ("dG_interface", 2)):
            clean[key] = _round(clean.get(key), digits)
        writer.writerow({k: ("" if clean.get(k) is None else clean.get(k))
                         for k in EXPORT_FIELDS})
    return "\ufeff" + buffer.getvalue()


# ----------------------------------------------------------------------- files

def resolve(root: Path | str, run: str, rel: str) -> Path:
    """A file inside one run's delivery folder, for a path that came from a browser.

    Three things have to hold, and all three are checked here rather than at the
    route: the run is a real run, the path stays inside its delivery folder once
    symlinks have been followed, and the suffix is one we serve.
    """
    base = (run_dir(root, run) / DELIVERY).resolve()
    target = (base / (rel or "")).resolve()
    if not target.is_relative_to(base):
        raise Forbidden(f"path escapes the run: {rel!r}")
    if target.suffix.lower() not in SERVABLE:
        raise Forbidden(f"not a servable file type: {target.suffix!r}")
    if not target.is_file():
        raise NotFound(f"no such file: {rel}")
    return target
