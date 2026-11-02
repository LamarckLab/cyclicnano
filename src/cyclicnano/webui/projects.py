"""Projects: the directory a design target lives in.

A project is one run directory. The page calls it a project and the pipeline
calls it a run, and they are the same folder: `<root>/<name>/`, holding the
stage directories the pipeline writes plus a `project.json` this module adds.
One name for the person choosing a target, one for the tool executing stages.

Keeping them one-to-one is what lets the pipeline's own resume semantics do the
work. Submitting a job inside a project runs stages into that same directory, so
the stages already done are skipped and the project accumulates rather than
being rebuilt.

Creating one is a write and is gated by the same `--allow-submit` flag as
starting a job. Everything else here reads.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

from . import store

MANIFEST = "project.json"

# Cn for n from 2 to 12, the range the pipeline's contig arithmetic covers.
SYMMETRY = tuple(f"C{n}" for n in range(2, 13))

TITLE = re.compile(r"^[^\x00-\x1f]{1,80}$")


class Invalid(Exception):
    """A project that cannot be created as asked."""


class Exists(Exception):
    """A directory of that name is already there."""


def read(root: Path | str, run: str) -> dict:
    """The declared metadata of one project, or what can be inferred without it.

    A run that predates this file is still a project: its name is the directory
    and its symmetry is whatever its own tables measured. So the page does not
    need a migration to show the runs that already exist.
    """
    path = Path(root) / run / MANIFEST
    declared = {}
    if path.is_file():
        try:
            declared = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            declared = {}

    return {
        "name": run,
        "title": str(declared.get("title") or run),
        "note": str(declared.get("note") or ""),
        # Declared when it was created, measured when the run has output. Both
        # are reported: a project whose tables disagree with its own declaration
        # is worth seeing rather than quietly resolving.
        "symmetry_declared": declared.get("symmetry") or None,
        "monomer_length_declared": declared.get("monomer_length") or None,
        "created": declared.get("created") or None,
        "declared": bool(declared),
    }


def validate(form: dict) -> dict:
    """The fields of a new project, or a refusal.

    The name becomes a directory and a path component in a URL, so it is held to
    the same pattern a discovered run is; the title is only ever displayed, so it
    only has to be one line of text.
    """
    form = form or {}

    name = str(form.get("name") or "").strip()
    if not store.RUN_NAME.match(name) or name.startswith("_"):
        raise Invalid("A project name may use letters, digits, underscore, dot and "
                      "hyphen, and may not start with a symbol.")
    if len(name) > 64:
        raise Invalid("A project name may be at most 64 characters.")

    symmetry = str(form.get("symmetry") or "")
    if symmetry not in SYMMETRY:
        raise Invalid(f"Symmetry must be one of {', '.join(SYMMETRY)}.")

    title = str(form.get("title") or name).strip() or name
    if not TITLE.match(title):
        raise Invalid("A title must be a single line of at most 80 characters.")

    try:
        monomer_length = int(form.get("monomer_length") or 60)
    except (TypeError, ValueError):
        raise Invalid("Monomer length must be a number.") from None
    if not 20 <= monomer_length <= 400:
        raise Invalid("Monomer length must be between 20 and 400 residues.")

    note = str(form.get("note") or "").strip()[:2000]

    return {"name": name, "title": title, "symmetry": symmetry,
            "monomer_length": monomer_length, "note": note}


def create(root: Path | str, form: dict) -> dict:
    """Make the directory and write its manifest. No pipeline stage runs.

    Deliberately separate from submitting a job: naming a target and spending
    GPU hours on it are different decisions, and a project that exists with
    nothing in it is a useful state - it says what is intended before anything
    has been computed.
    """
    fields = validate(form)
    folder = Path(root) / fields["name"]
    if folder.exists():
        raise Exists(f"{fields['name']} already exists under the output root.")

    folder.mkdir(parents=True)
    record = {**fields, "created": int(time.time())}
    (folder / MANIFEST).write_text(
        json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8")
    return record
