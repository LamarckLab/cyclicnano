"""Launching pipeline runs from the browser.

This is the one part of the web UI that is not read-only, and it is off unless
the service was started with `--allow-submit`. There is no authentication, so
whoever can reach the port can spend GPU hours on a shared machine; that is a
deployment decision and the flag is there to make it one.

Three rules hold everything together.

Every field is declared with a type and a range, and nothing reaches the command
line that has not been through `validate`. The command is built as a list and
spawned without a shell, so a value cannot become an argument.

One job at a time. The design server is shared with other people's work, and two
concurrent pipelines would contend for the same cards.

A job gets its own session, so cancelling it kills the whole tree. The pipeline
spawns RFdiffusion, ProteinMPNN and colabfold as children; signalling only the
parent would leave a GPU occupied by an orphan.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

from . import projects

JOBS_DIR = "_jobs"          # under the output root; the underscore keeps it out of list_runs
RECORD = "job.json"
LOG = "log.txt"

# A run name becomes a directory under the output root and a path component in a
# URL, so it is held to the same pattern a discovered run is.
RUN_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")

JOB_ID = re.compile(r"^[0-9]{8}-[0-9]{6}-[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")


class Invalid(Exception):
    """A submitted form that cannot be turned into a run."""


class Busy(Exception):
    """A job is already running, and there is only one slot."""


class Disabled(Exception):
    """Submission was not enabled on this service."""


# --------------------------------------------------------------------- fields

@dataclass
class Field:
    """One form field, and the config key it overrides.

    `key` is the dotted path the pipeline knows it by, which is also what gets
    written into the job record: a reader of that record can reproduce the run
    from the CLI without consulting this file.
    """
    name: str
    key: str
    label: str
    kind: str                                   # int | float | name | choice
    default: object
    minimum: float | None = None
    maximum: float | None = None
    choices: tuple = ()
    help: str = ""
    warn_above: float | None = None


# The `label` of each field is an i18n key, not a label. Returning display text
# from here would let the API speak only one language, which is the mistake the
# structure labels made before they were moved out of `store`. The page looks
# each key up in static/i18n.js; `help` is the same key with `.help` appended
# when a field has one.
FIELDS: tuple[Field, ...] = (
    Field("name", "run.name", "jobs.field.name", "name", ""),
    Field("symmetry", "target.symmetry", "jobs.field.symmetry", "choice", "C5",
          # One list, in projects.py, so a project cannot declare a symmetry a
          # job would refuse or the other way round.
          choices=projects.SYMMETRY),
    Field("monomer_length", "target.monomer_length", "jobs.field.monomer_length",
          "int", 60, minimum=20, maximum=400),
    Field("n_backbones", "target.n_backbones", "jobs.field.n_backbones", "int", 500,
          minimum=1, maximum=2000, warn_above=500),
    Field("noise_scale", "generate.noise_scale", "jobs.field.noise_scale", "float",
          0.5, minimum=0.0, maximum=2.0),
    Field("batch_size", "generate.batch_size", "jobs.field.batch_size", "int", 10,
          minimum=1, maximum=200),
    Field("n_seq_per_backbone", "design.n_seq_per_backbone",
          "jobs.field.n_seq_per_backbone", "int", 20, minimum=1, maximum=200),
    Field("rmsd_cutoff", "validate.rmsd_cutoff", "jobs.field.rmsd_cutoff", "float",
          1.0, minimum=0.1, maximum=10.0),
    Field("plddt_cutoff", "validate.plddt_cutoff", "jobs.field.plddt_cutoff", "float",
          80.0, minimum=0.0, maximum=100.0),
    Field("seed", "run.seed", "jobs.field.seed", "int", 1,
          minimum=0, maximum=2**31 - 1),
)

PRESETS = ("pure_three_helix", "none")

# What to stop after. Seeing how many backbones survive stage 01 costs minutes;
# finding out after ProteinMPNN and AlphaFold have run costs hours.
STAGE_SETS = {
    "all": ("00", "01", "02", "03", "04", "05"),
    "through_delivery": ("00", "01", "02", "03", "04"),
    "backbones_only": ("00", "01"),
}


def form_spec() -> dict:
    """The fields, for the browser to render. Nothing secret, and no paths."""
    return {
        "fields": [{
            "name": f.name, "key": f.key, "label": f.label, "kind": f.kind,
            "default": f.default, "minimum": f.minimum, "maximum": f.maximum,
            "choices": list(f.choices), "help": f.help, "warn_above": f.warn_above,
        } for f in FIELDS],
        "presets": list(PRESETS),
        "stage_sets": {k: list(v) for k, v in STAGE_SETS.items()},
    }


def validate(form: dict) -> dict:
    """Turn a submitted form into overrides, or refuse it.

    Refusing is the point. Everything here ends up on a command line that starts
    a job on a shared GPU server, so a value that is not a number in range, or a
    name that is not a plain name, does not get a best effort.

    The messages name the field by its key, not by its label, and are English on
    a page that is also Chinese. A label here is an i18n key, so interpolating
    one would print `jobs.field.seed` at the reader; and a refusal from this
    layer is read while debugging a request, where the config key the value maps
    to is the useful name. The form prevents the ordinary cases client-side.
    """
    form = form or {}
    overrides: dict[str, object] = {}

    for f in FIELDS:
        raw = form.get(f.name, f.default)
        if raw is None or raw == "":
            if f.kind == "name":
                raise Invalid(f"{f.key} must not be empty")
            raw = f.default

        if f.kind == "name":
            value = str(raw).strip()
            if not RUN_NAME.match(value):
                raise Invalid(f"{f.key} may use letters, digits, underscore, dot and "
                              f"hyphen, and may not start with a symbol: {value!r}")
        elif f.kind == "choice":
            value = str(raw)
            if value not in f.choices:
                raise Invalid(f"{f.key} must be one of "
                              f"{', '.join(map(str, f.choices))}: {value!r}")
        else:
            try:
                value = int(raw) if f.kind == "int" else float(raw)
            except (TypeError, ValueError):
                raise Invalid(f"{f.key} must be a number: {raw!r}") from None
            if f.minimum is not None and value < f.minimum:
                raise Invalid(f"{f.key} must be at least {f.minimum}")
            if f.maximum is not None and value > f.maximum:
                raise Invalid(f"{f.key} must be at most {f.maximum}")

        overrides[f.key] = value

    preset = str(form.get("preset", PRESETS[0]))
    if preset not in PRESETS:
        raise Invalid(f"preset must be one of {', '.join(PRESETS)}: {preset!r}")

    stages = str(form.get("stages", "all"))
    if stages not in STAGE_SETS:
        raise Invalid(f"stages must be one of {', '.join(STAGE_SETS)}: {stages!r}")

    return {"overrides": overrides, "preset": preset, "stages": stages}


# -------------------------------------------------------------------- command

def build_command(repo: Path, parsed: dict, profile: str | None) -> list[str]:
    """The CLI invocation for a validated form.

    A list, never a string: the values came from a browser, and a shell would
    give them a second chance to mean something. `-u` because a buffered child
    writes a log that looks empty for the first half hour.
    """
    repo = Path(repo)
    symmetry = str(parsed["overrides"]["target.symmetry"]).lower()

    configs = [repo / "configs" / "runs" / "base.yaml"]
    # The per-symmetry file exists for C5 and C3 and sets the run name as well,
    # so it is only used when it matches, and every value it holds is overridden
    # below anyway.
    per_symmetry = repo / "configs" / "runs" / f"{symmetry}.yaml"
    if per_symmetry.is_file():
        configs.append(per_symmetry)
    if parsed["preset"] != "none":
        configs.append(repo / "configs" / "presets" / f"{parsed['preset']}.yaml")

    for path in configs:
        if not path.is_file():
            raise Invalid(f"config file is missing: {path}")

    command = [sys.executable, "-u", "-m", "cyclicnano.cli", "run"]
    for path in configs:
        command += ["--config", str(path)]
    if profile:
        if not Path(profile).is_file():
            raise Invalid(f"machine profile is missing: {profile}")
        command += ["--profile", str(profile)]
    command += ["--stages", ",".join(STAGE_SETS[parsed["stages"]])]
    for key, value in parsed["overrides"].items():
        command += ["--set", f"{key}={value}"]
    return command


# ---------------------------------------------------------------------- state

def jobs_dir(root: Path | str) -> Path:
    return Path(root) / JOBS_DIR


def _alive(pid: int | None) -> bool:
    """Whether a process id is still a process.

    Signal 0 is a probe and not a signal: it checks that the process exists
    without touching it. PermissionError is caught first and means yes, because
    a process owned by someone else is still a process holding a card.

    The remaining OSError is what a pid that cannot exist raises, and the type
    differs by platform: ProcessLookupError on POSIX, a plain OSError on
    Windows. For a liveness probe, any other failure to find it means it is not
    there.
    """
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except PermissionError:
        return True
    except (OSError, ValueError, OverflowError):
        return False
    return True


def _read(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _write(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)                                # never a half-written record


def _settle(root: Path, record: dict) -> dict:
    """Bring a stored record up to date with the operating system.

    A job marked running whose process is gone was not reaped by this service:
    the service restarted while it was running. The exit code is unrecoverable
    in that case, so the state says `unknown` rather than guessing, and points
    at the pipeline's own completion marker, which is the only evidence left.
    """
    if record.get("state") != "running" or _alive(record.get("pid")):
        return record

    run = record.get("run")
    done = bool(run) and (Path(root) / str(run) / ".done").exists()
    record["state"] = "unknown"
    record["note"] = (
        "The service restarted while this job was running, so its exit code was "
        "never collected. "
        + ("The run directory has a .done marker, so it looks finished."
           if done else
           "The run directory has no .done marker, so it may not have finished; "
           "read the log."))
    record["finished"] = record.get("finished") or int(time.time())
    _write(jobs_dir(root) / str(record["id"]) / RECORD, record)
    return record


def listing(root: Path | str, limit: int = 25) -> list[dict]:
    """Every job, newest first, with its state checked against the system."""
    base = jobs_dir(root)
    if not base.is_dir():
        return []
    out = []
    for child in sorted(base.iterdir(), reverse=True):
        record = _read(child / RECORD) if child.is_dir() else None
        if record:
            out.append(_settle(Path(root), record))
        if len(out) >= limit:
            break
    return out


def running(root: Path | str) -> dict | None:
    return next((j for j in listing(root) if j.get("state") == "running"), None)


def get(root: Path | str, job_id: str) -> dict:
    if not JOB_ID.match(job_id or ""):
        raise Invalid(f"not a job id: {job_id!r}")
    record = _read(jobs_dir(root) / job_id / RECORD)
    if not record:
        raise Invalid(f"no such job: {job_id}")
    return _settle(Path(root), record)


def log_tail(root: Path | str, job_id: str, offset: int = 0,
             limit: int = 64 * 1024) -> dict:
    """The log from a byte offset, so the page can follow it without refetching.

    Reading bytes rather than lines because the offset has to survive a partial
    line: the pipeline is writing to this file while it is read.
    """
    get(root, job_id)                                # validates the id
    path = jobs_dir(root) / job_id / LOG
    if not path.is_file():
        return {"offset": 0, "text": "", "size": 0}

    size = path.stat().st_size
    offset = max(0, min(int(offset), size))
    with open(path, "rb") as fh:
        fh.seek(offset)
        blob = fh.read(limit)
    return {"offset": offset + len(blob), "size": size,
            "text": blob.decode("utf-8", errors="replace")}


# --------------------------------------------------------------------- submit

@dataclass
class Launcher:
    """What it takes to start a run, held apart from the routes.

    The repository and the machine profile are decided when the service starts,
    not per request: a browser chooses parameters, never which code or which
    machine configuration runs.
    """
    root: Path
    repo: Path
    profile: str | None = None
    enabled: bool = False
    _threads: list = field(default_factory=list)

    def submit(self, form: dict) -> dict:
        if not self.enabled:
            raise Disabled("Job submission is not enabled on this service; "
                        "start it with --allow-submit")

        current = running(self.root)
        if current:
            raise Busy(f"a job is already running: {current['id']} "
                       f"({current.get('run')}). This machine is shared, so one "
                       "job runs at a time.")

        parsed = validate(form)
        command = build_command(self.repo, parsed, self.profile)
        run = str(parsed["overrides"]["run.name"])
        job_id = f"{time.strftime('%Y%m%d-%H%M%S')}-{run}"

        folder = jobs_dir(self.root) / job_id
        folder.mkdir(parents=True, exist_ok=True)
        log_path = folder / LOG

        record = {
            "id": job_id,
            "run": run,
            "state": "running",
            "started": int(time.time()),
            "finished": None,
            "returncode": None,
            "pid": None,
            "stages": parsed["stages"],
            "preset": parsed["preset"],
            "overrides": parsed["overrides"],
            # The exact invocation, so the run can be repeated from a terminal
            # without reconstructing it from the form.
            "command": " ".join(shlex.quote(c) for c in command),
            "resumed": (Path(self.root) / run).is_dir(),
        }

        with open(log_path, "ab") as log:
            log.write(f"$ {record['command']}\n\n".encode())
            log.flush()
            process = subprocess.Popen(
                command, cwd=str(self.repo), stdout=log, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                # Its own session, so cancelling signals the whole tree: the
                # pipeline's own children hold the GPU, not this process.
                start_new_session=True,
            )

        record["pid"] = process.pid
        _write(folder / RECORD, record)

        thread = threading.Thread(target=self._reap, args=(process, folder, job_id),
                                  daemon=True)
        thread.start()
        self._threads.append(thread)
        return record

    def _reap(self, process: subprocess.Popen, folder: Path, job_id: str) -> None:
        code = process.wait()
        record = _read(folder / RECORD) or {"id": job_id}
        # Cancel marks the record before signalling, so a cancelled job is not
        # relabelled as a failure by the exit code the signal produced.
        if record.get("state") != "cancelled":
            record["state"] = "finished" if code == 0 else "failed"
        record["returncode"] = code
        record["finished"] = int(time.time())
        _write(folder / RECORD, record)

    def cancel(self, job_id: str) -> dict:
        if not self.enabled:
            raise Disabled("Job submission is not enabled on this service; "
                        "start it with --allow-submit")

        record = get(self.root, job_id)
        if record.get("state") != "running":
            raise Invalid(f"this job is not running ({record.get('state')})")

        pid = record.get("pid")
        record["state"] = "cancelled"
        record["note"] = "Cancelled from the page."
        _write(jobs_dir(self.root) / job_id / RECORD, record)

        if _alive(pid):
            # The group, not the process: RFdiffusion and colabfold are children
            # and they are the ones holding a card.
            try:
                os.killpg(os.getpgid(int(pid)), signal.SIGTERM)
            except AttributeError:                      # no process groups here
                os.kill(int(pid), signal.SIGTERM)
            except (ProcessLookupError, PermissionError, OSError):
                pass
        return record
