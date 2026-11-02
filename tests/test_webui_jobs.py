"""Submitting pipeline runs from the browser.

This module is the only part of the web UI that spends anything, so what is
tested here is the refusals: which forms are rejected, what the command comes
out as, and that the state machine never reports a job as having succeeded when
it does not know.

Nothing is spawned. The three lines that start a process are exercised against
the real pipeline on the server; everything around them is pure and lives here.
"""
from __future__ import annotations

import contextlib
import json
import shutil
import tempfile
import time
from pathlib import Path

from cyclicnano.webui import jobs


@contextlib.contextmanager
def workspace():
    """An output root and a stand-in checkout holding the configs we name."""
    tmp = Path(tempfile.mkdtemp(prefix="cyclicnano_jobs_"))
    try:
        root = tmp / "outputs"
        root.mkdir()
        repo = tmp / "repo"
        for rel in ("configs/runs/base.yaml", "configs/runs/c5.yaml",
                    "configs/runs/c3.yaml", "configs/presets/pure_three_helix.yaml"):
            path = repo / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("{}\n", encoding="utf-8")
        yield root, repo
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def raises(exc, fn, *args, **kwargs):
    want = " or ".join(e.__name__ for e in (exc if isinstance(exc, tuple) else (exc,)))
    try:
        fn(*args, **kwargs)
    except exc:
        return
    except Exception as other:
        raise AssertionError(f"expected {want}, got {other!r}") from other
    raise AssertionError(f"expected {want} from {fn.__name__}")


def record(root: Path, job_id: str, **fields) -> dict:
    row = {"id": job_id, "run": fields.pop("run", "c5"), "state": "running",
           "started": int(time.time()), "finished": None, "returncode": None,
           "pid": None, "stages": "all", "preset": "pure_three_helix",
           "overrides": {}, "command": "python -m cyclicnano.cli run"}
    row.update(fields)
    folder = jobs.jobs_dir(root) / job_id
    folder.mkdir(parents=True, exist_ok=True)
    (folder / jobs.RECORD).write_text(json.dumps(row), encoding="utf-8")
    return row


# ------------------------------------------------------------------ validation

def test_nothing_reaches_the_command_line_without_passing_a_range():
    """Every field is declared with a type and bounds, and this is the gate.

    The values come from a browser on an unauthenticated port and end up
    starting a job on a shared GPU server, so a bad one gets a refusal rather
    than a best effort.
    """
    good = jobs.validate({"name": "c5_test", "symmetry": "C5", "n_backbones": 50})
    assert good["overrides"]["run.name"] == "c5_test"
    assert good["overrides"]["target.n_backbones"] == 50
    assert good["overrides"]["target.symmetry"] == "C5"
    # An omitted field takes the declared default rather than being left out,
    # so the job record always states the whole configuration.
    assert set(good["overrides"]) == {f.key for f in jobs.FIELDS}

    for bad in (
        {"name": ""},                                  # a run must be named
        {"name": "../evil"},                           # a path, not a name
        {"name": "/abs"},
        {"name": "has space"},
        {"name": "-leading"},
        {"name": "x", "n_backbones": 99999},           # over the cap
        {"name": "x", "n_backbones": 0},
        {"name": "x", "symmetry": "C99"},
        {"name": "x", "symmetry": "C5; rm -rf /"},
        {"name": "x", "noise_scale": "abc"},
        {"name": "x", "noise_scale": -1},
        {"name": "x", "plddt_cutoff": 101},
        {"name": "x", "preset": "/etc/passwd"},        # not an allowed preset
        {"name": "x", "preset": "../../secret"},
        {"name": "x", "stages": "rm -rf"},
    ):
        raises(jobs.Invalid, jobs.validate, bad)


def test_the_form_spec_and_the_gate_cannot_disagree():
    """The browser renders from the spec, so the spec has to be the same rules.

    A field present in the form but absent from the spec would be enforced and
    never shown; one in the spec but not enforced would be shown and ignored.
    """
    spec = jobs.form_spec()
    assert [f["name"] for f in spec["fields"]] == [f.name for f in jobs.FIELDS]
    assert spec["presets"] == list(jobs.PRESETS)
    assert set(spec["stage_sets"]) == set(jobs.STAGE_SETS)
    for shown, declared in zip(spec["fields"], jobs.FIELDS):
        assert shown["minimum"] == declared.minimum
        assert shown["maximum"] == declared.maximum
        assert shown["choices"] == list(declared.choices)


# --------------------------------------------------------------------- command

def test_the_command_is_a_list_and_the_overrides_come_last():
    """Config layering decides which value wins, so the order is the behaviour.

    `configs/runs/c5.yaml` sets `run.name: c5`. A run the user named something
    else only gets that name because `--set` is applied after every `--config`,
    so the arguments have to be emitted in that order.
    """
    with workspace() as (root, repo):
        parsed = jobs.validate({"name": "c5_quick", "symmetry": "C5",
                                "stages": "backbones_only"})
        command = jobs.build_command(repo, parsed, None)

        assert isinstance(command, list)
        assert all(isinstance(c, str) for c in command)
        assert command[1:5] == ["-u", "-m", "cyclicnano.cli", "run"]

        configs = [command[i + 1] for i, c in enumerate(command) if c == "--config"]
        assert [Path(c).name for c in configs] == [
            "base.yaml", "c5.yaml", "pure_three_helix.yaml"]

        assert command[command.index("--stages") + 1] == "00,01"
        assert "run.name=c5_quick" in command
        assert command.index("--set") > max(
            i for i, c in enumerate(command) if c == "--config")


def test_a_symmetry_with_no_run_config_still_works_and_no_preset_adds_no_file():
    """Only C5 and C3 have a run config; the rest are the same pipeline.

    Naming a file that is not there would fail the run minutes later inside the
    CLI, so the per-symmetry config is used only when it exists and every value
    it would have set is overridden anyway.
    """
    with workspace() as (root, repo):
        parsed = jobs.validate({"name": "c7_try", "symmetry": "C7", "preset": "none"})
        command = jobs.build_command(repo, parsed, None)
        configs = [command[i + 1] for i, c in enumerate(command) if c == "--config"]
        assert [Path(c).name for c in configs] == ["base.yaml"]
        assert "target.symmetry=C7" in command
        assert not any("pure_three_helix" in c for c in command)

        # A profile that is not there is refused here rather than by the CLI.
        raises(jobs.Invalid, jobs.build_command, repo, parsed, str(repo / "nope.yaml"))


# ----------------------------------------------------------------- the refusals

def test_submission_is_off_until_it_is_switched_on():
    """The page has no authentication, so this is the whole access control.

    Anyone who can reach the port could otherwise spend GPU hours on a machine
    other people are using, so the capability is a deployment decision.
    """
    with workspace() as (root, repo):
        off = jobs.Launcher(root=root, repo=repo, enabled=False)
        raises(jobs.Disabled, off.submit, {"name": "c5"})
        raises(jobs.Disabled, off.cancel, "20260101-000000-c5")


def test_one_slot_only():
    """Two pipelines on this machine would contend for the same cards.

    The refusal names the job that holds the slot, because the next question is
    always which one.
    """
    with workspace() as (root, repo):
        launcher = jobs.Launcher(root=root, repo=repo, enabled=True)
        # A record whose pid is this process: alive, so the slot is taken.
        import os
        record(root, "20260101-000000-busy", run="busy", pid=os.getpid())
        assert jobs.running(root)["run"] == "busy"
        raises(jobs.Busy, launcher.submit, {"name": "c5_other"})


# -------------------------------------------------------------- the state machine

def test_an_orphaned_job_is_unknown_rather_than_guessed():
    """A restart loses the exit code, and inventing one would be a lie.

    The pipeline's own `.done` marker is the only evidence left, so the note
    says what it shows and the state says the code is not known. Calling it
    "finished" because the process is gone would report a crash as a success.
    """
    with workspace() as (root, repo):
        # pid 1 is not ours to signal but is alive, so use a pid that cannot be.
        record(root, "20260101-000001-gone", run="gone", pid=2 ** 22)
        settled = jobs.get(root, "20260101-000001-gone")
        assert settled["state"] == "unknown"
        assert settled["returncode"] is None
        assert "may not have finished" in settled["note"]
        assert settled["finished"]

        (root / "done_run").mkdir()
        (root / "done_run" / ".done").touch()
        record(root, "20260101-000002-done", run="done_run", pid=2 ** 22)
        assert "looks finished" in jobs.get(root, "20260101-000002-done")["note"]

        # The settled state is written back, so it is not recomputed each read.
        stored = json.loads((jobs.jobs_dir(root) / "20260101-000001-gone"
                             / jobs.RECORD).read_text(encoding="utf-8"))
        assert stored["state"] == "unknown"


def test_cancelling_kills_the_process_and_the_state_survives_the_exit_code():
    """SIGTERM gives a non-zero exit code, and that is not a failure.

    Cancel marks the record before signalling precisely so the reaper, which
    only sees the code, leaves the state alone. Run against a real throwaway
    child rather than a fabricated pid, because the signalling is the part worth
    testing: in production the signal goes to the process group, since
    RFdiffusion and colabfold are children and they are the ones holding a card.
    """
    import subprocess
    import sys
    import time as clock

    with workspace() as (root, repo):
        launcher = jobs.Launcher(root=root, repo=repo, enabled=True)
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"],
                                 start_new_session=True)
        try:
            record(root, "20260101-000003-x", run="x", pid=child.pid)
            assert jobs.get(root, "20260101-000003-x")["state"] == "running"

            cancelled = launcher.cancel("20260101-000003-x")
            assert cancelled["state"] == "cancelled"

            for _ in range(50):                        # it does not die instantly
                if child.poll() is not None:
                    break
                clock.sleep(0.1)
            assert child.poll() is not None, "the child outlived its cancellation"

            launcher._reap(child, jobs.jobs_dir(root) / "20260101-000003-x",
                           "20260101-000003-x")
            after = jobs.get(root, "20260101-000003-x")
            assert after["state"] == "cancelled"       # not "failed"
            assert after["returncode"] == child.returncode != 0

            raises(jobs.Invalid, launcher.cancel, "20260101-000003-x")  # already done
            raises(jobs.Invalid, launcher.cancel, "../../etc/passwd")
        finally:
            if child.poll() is None:
                child.kill()


def test_the_log_is_followed_by_byte_offset():
    """The pipeline is writing to this file while the page reads it.

    Line-based reads would lose or duplicate a partially written line at the
    boundary, so the offset counts bytes and a decode error degrades to a
    replacement character instead of failing the request.
    """
    with workspace() as (root, repo):
        record(root, "20260101-000004-log", run="log", pid=2 ** 22)
        path = jobs.jobs_dir(root) / "20260101-000004-log" / jobs.LOG
        path.write_bytes(b"first\nsecond\n")

        head = jobs.log_tail(root, "20260101-000004-log", 0)
        assert head["text"] == "first\nsecond\n" and head["offset"] == 13

        with open(path, "ab") as fh:
            fh.write("第三行\n".encode())
        more = jobs.log_tail(root, "20260101-000004-log", head["offset"])
        assert more["text"] == "第三行\n"
        assert jobs.log_tail(root, "20260101-000004-log", more["offset"])["text"] == ""

        # An offset past the end is clamped rather than raising.
        assert jobs.log_tail(root, "20260101-000004-log", 10 ** 9)["text"] == ""
        raises(jobs.Invalid, jobs.log_tail, root, "no such job", 0)
