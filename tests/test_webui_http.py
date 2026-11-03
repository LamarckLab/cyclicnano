"""The web UI's HTTP layer.

Separate from `test_webui.py` on purpose. That module tests the reading layer
with the standard library alone, so the suite runs on a machine with no web
stack; this one needs FastAPI and its test client.

The skip is deliberately narrow. If FastAPI is missing the machine cannot serve
the page at all and there is nothing to test, so the module stands down. But if
FastAPI is installed and only the test client is missing, that is a missing dev
dependency rather than a machine without a web stack, and the module says so
instead of passing quietly: a test that skips itself looks like coverage.
"""
from __future__ import annotations

import importlib.util

from test_webui import outputs

HAVE_FASTAPI = importlib.util.find_spec("fastapi") is not None
HAVE_CLIENT = importlib.util.find_spec("httpx") is not None


def client_for(root):
    from fastapi.testclient import TestClient

    from cyclicnano.webui.app import create_app
    return TestClient(create_app(root))


_announced = False


def _unavailable() -> bool:
    """Whether to stand down, saying so out loud.

    `tests/run_all.py` counts any test that does not raise as a pass, so a
    silent early return here would be indistinguishable from coverage. It
    prints instead, once, so the run says what it did not test.
    """
    global _announced
    if not HAVE_FASTAPI:
        if not _announced:
            print("  [skip] 本机没有 fastapi，HTTP 层未测试"
                  "（装上 pip install -e \".[webui,dev]\" 再跑）")
            _announced = True
        return True
    assert HAVE_CLIENT, ("fastapi is installed but httpx is not, so the HTTP layer "
                         "is untested; install with pip install -e \".[dev]\"")
    return False


# ----------------------------------------------------------------- the routes

def test_the_three_levels_answer():
    if _unavailable():
        return
    with outputs() as root:
        c = client_for(root)

        health = c.get("/api/health").json()
        assert health["ok"] is True and health["root"] == str(root)
        # Off unless the service was started with --allow-submit. The page reads
        # this to leave the submit form out entirely, so the default matters.
        assert health["can_submit"] is False

        runs = c.get("/api/runs").json()
        assert [r["run"] for r in runs] == ["c5"]
        assert runs[0]["n_generated"] == 5 and runs[0]["n_final_pass"] == 1

        rows = c.get("/api/runs/c5/backbones").json()
        assert [r["backbone_id"] for r in rows] == ["design_27", "design_14"]

        detail = c.get("/api/runs/c5/backbones/design_27").json()
        assert detail["folder"] == "passed/1_design_27"
        assert len(detail["sequences"]) == 2

        assert c.get("/").status_code == 200
        assert "text/html" in c.get("/").headers["content-type"]


def test_a_refused_path_keeps_its_status():
    """403 for a path that leaves the run, 404 for something that is not there.

    The distinction is the point: a hidden run must not be confirmed to exist,
    and a traversal attempt must not read as a typo.
    """
    if _unavailable():
        return
    with outputs() as root:
        c = client_for(root)
        assert c.get("/api/runs/c5/file", params={"path": "../../etc/passwd"}).status_code == 403
        assert c.get("/api/runs/c5/file", params={"path": "notes.txt"}).status_code == 403
        assert c.get("/api/runs/_archive/backbones").status_code == 404
        assert c.get("/api/runs/c5/backbones/design_999").status_code == 404


# -------------------------------------------------------------- the download

def test_a_structure_is_viewed_as_text_and_saved_as_a_pdb():
    """The viewer and a save want different things from the same bytes.

    Serving a PDB as text/plain to a download is what made browsers write
    `.txt` files: the URL carries no filename, so the browser invents one from
    the media type. The download therefore has to name the file and use a type
    the browser will not rename, while the viewer keeps text/plain because it
    parses the bytes itself.
    """
    if _unavailable():
        return
    with outputs() as root:
        c = client_for(root)
        path = "passed/1_design_27/design_27_bb.pdb"

        shown = c.get("/api/runs/c5/file", params={"path": path})
        assert shown.status_code == 200
        assert shown.headers["content-type"].startswith("text/plain")
        assert "content-disposition" not in shown.headers
        assert shown.text.startswith("ATOM")

        saved = c.get("/api/runs/c5/file", params={"path": path, "download": 1})
        assert saved.status_code == 200
        assert saved.headers["content-type"] == "application/octet-stream"
        disposition = saved.headers["content-disposition"]
        assert "attachment" in disposition
        assert 'filename="design_27_bb.pdb"' in disposition
        assert saved.content == shown.content

        # Both forms are uncached: a re-run overwrites these paths in place.
        assert shown.headers["cache-control"] == "no-store"
        assert saved.headers["cache-control"] == "no-store"


def test_the_ring_is_assembled_on_request_and_not_read_off_disk():
    """It has no file, so the route is given the two files it is built from.

    Both go through the same resolution as any other path, which is what keeps a
    browser from pairing something outside the delivery folder. The result is a
    structure the viewer can load directly and a save can name.
    """
    if _unavailable():
        return
    with outputs() as root:
        c = client_for(root)
        source = "passed/1_design_27/predictions/seq000_af2.pdb"
        design = "passed/1_design_27/design_27_bb.pdb"

        shown = c.get("/api/runs/c5/assembly", params={"path": source, "design": design})
        assert shown.status_code == 200
        assert shown.headers["content-type"].startswith("text/plain")
        assert shown.headers["cache-control"] == "no-store"
        chains = {line[21] for line in shown.text.splitlines() if line.startswith("ATOM")}
        assert chains == {"A", "B", "C", "D", "E"}

        # Five copies of the one-chain prediction, so five times its atoms.
        monomer = c.get("/api/runs/c5/file", params={"path": source}).text
        n_atoms = sum(1 for line in monomer.splitlines() if line.startswith("ATOM"))
        assert sum(1 for line in shown.text.splitlines()
                   if line.startswith("ATOM")) == 5 * n_atoms

        saved = c.get("/api/runs/c5/assembly",
                      params={"path": source, "design": design, "download": 1})
        assert saved.headers["content-type"] == "application/octet-stream"
        assert 'filename="seq000_af2_C5.pdb"' in saved.headers["content-disposition"]
        assert saved.text == shown.text

        # Neither path may leave the run, and neither may name something that is
        # not a structure.
        for bad in ({"path": "../../../etc/passwd", "design": design},
                    {"path": source, "design": "../../secret.pdb"},
                    {"path": "manifest.json", "design": design}):
            assert c.get("/api/runs/c5/assembly", params=bad).status_code in (403, 404)


def test_the_page_is_never_served_from_a_stale_cache():
    """A deploy that the browser does not see is a deploy that did not happen.

    The API routes said no-store from the start and these said nothing, so the
    browser applied its own heuristic: a CSS fix and a JavaScript feature were
    copied to the server, verified byte for byte over HTTP, and still did not
    reach the page. The vendored bundle is the one exception, pinned by version
    and hash, so it cannot change under a URL that stays the same.
    """
    if _unavailable():
        return
    with outputs() as root:
        c = client_for(root)
        for name in ("app.js", "style.css", "i18n.js", "index.html"):
            r = c.get(f"/static/{name}")
            if r.status_code != 200:
                continue                       # not every file is in the checkout
            assert r.headers["cache-control"] == "no-store", name


def test_the_export_route_hands_over_a_named_file():
    """It exists to hand a sheet to whoever orders the synthesis.

    So it is always an attachment with a filename that says which run and which
    scope it holds; a FASTA rendered in a browser tab is one step away from
    being pasted by hand.
    """
    if _unavailable():
        return
    with outputs() as root:
        c = client_for(root)

        fasta = c.get("/api/runs/c5/export", params={"format": "fasta", "scope": "delivered"})
        assert fasta.status_code == 200
        assert fasta.headers["content-type"] == "application/octet-stream"
        assert 'filename="c5_sequences_delivered.fasta"' in fasta.headers["content-disposition"]
        assert fasta.text.startswith(">design_27__seq000 C5 ")

        csv = c.get("/api/runs/c5/export", params={"format": "csv", "scope": "all"})
        assert csv.status_code == 200
        assert csv.headers["content-type"].startswith("text/csv")
        assert 'filename="c5_sequences_all.csv"' in csv.headers["content-disposition"]
        assert csv.content.startswith(b"\xef\xbb\xbf")      # the BOM Excel needs

        assert c.get("/api/runs/c5/export", params={"format": "xlsx"}).status_code == 403
        assert c.get("/api/runs/c5/export", params={"scope": "whatever"}).status_code == 403


def test_the_job_routes_refuse_everything_until_submission_is_enabled():
    """Read-only by default, including the routes that would start a run.

    The page has no authentication, so this is the whole of its access control:
    a service started without --allow-submit answers 403 to a submit and to a
    cancel, and says so in /api/jobs so the form is never drawn.
    """
    if _unavailable():
        return
    with outputs() as root:
        c = client_for(root)

        info = c.get("/api/jobs").json()
        assert info["enabled"] is False
        assert info["running"] is None and info["jobs"] == []
        # The spec is served either way, so the page can explain what it would ask.
        assert [f["name"] for f in info["spec"]["fields"]][0] == "name"

        assert c.post("/api/jobs", json={"name": "c9"}).status_code == 403
        assert c.post("/api/jobs/20260101-000000-c9/cancel").status_code == 403

        # A bad job id is a 400 whether or not submission is on.
        assert c.get("/api/jobs/not-an-id").status_code == 400
        assert c.get("/api/jobs/not-an-id/log").status_code == 400
