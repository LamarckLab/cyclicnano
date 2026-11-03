"""The HTTP layer of the web UI.

Read-only on purpose. Every route is a thin wrapper over `store`, which holds
all the knowledge of the output layout, so the interesting behaviour can be
tested without starting a server. Nothing here launches a pipeline run, touches
a GPU or writes to the output tree.

There is no authentication, because the service is meant to be reachable only
from inside the lab network. That is a deployment assumption, not a property of
this code: do not put it on a public address.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

from fastapi import Body, FastAPI, Query, Response
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from ..assembly import place_monomer
from . import interface, jobs, projects, store

STATIC = Path(__file__).resolve().parent / "static"


class _PageFiles(StaticFiles):
    """The page's own files, served with a cache policy of their own.

    The API routes have said `no-store` from the start; these had said nothing,
    so a browser applied its heuristic and a deploy could leave someone looking
    at the previous page. That happened: a CSS fix and a JavaScript feature were
    copied to the server, verified byte for byte over HTTP, and still did not
    appear in the browser.

    The vendored bundle is the exception. It is pinned by version and checked
    against a hash, so it cannot change under a URL that stays the same, and
    re-sending five megabytes on every page load to protect against a change
    that cannot happen is not a trade worth making.
    """

    def file_response(self, full_path, stat_result, scope, status_code=200):
        response = super().file_response(full_path, stat_result, scope, status_code)
        cached = "vendor" in Path(full_path).parts
        response.headers["Cache-Control"] = ("public, max-age=31536000, immutable"
                                             if cached else "no-store")
        return response


DEFAULT_ROOT = "/data/lmk/cyclicnano_outputs"
ENV_ROOT = "CYCLICNANO_OUTPUTS"


def _installed_version() -> str:
    """The version this service is running, for the footer.

    Read from the installed distribution rather than kept in a constant here, so
    pyproject.toml stays the one place it is written. A checkout that was never
    installed - which is how the test suite runs - has no metadata to read, and
    the answer is an empty string: the footer then shows a name without a
    version instead of a name with a wrong one.
    """
    try:
        from importlib.metadata import PackageNotFoundError, version
        return version("cyclicnano")
    except (PackageNotFoundError, ImportError, ValueError):
        return ""



def output_root() -> Path:
    return Path(os.environ.get(ENV_ROOT) or DEFAULT_ROOT)


def repo_root() -> Path:
    """The checkout this package was installed from.

    src/cyclicnano/webui/app.py -> three parents up is the repository, which is
    where `configs/` lives. Needed only to launch a run; the reading side never
    touches the code directory.
    """
    return Path(__file__).resolve().parents[3]


def create_app(root: Path | str | None = None, *, allow_submit: bool = False,
               repo: Path | str | None = None, profile: str | None = None) -> FastAPI:
    """The application, reading one output root.

    The root is captured here rather than read per request so that a running
    service cannot be pointed somewhere else by changing the environment under
    it, and so tests can hold several roots at once. The same goes for the
    repository and the machine profile: a browser picks parameters, never which
    code or which machine configuration a run uses.
    """
    root = Path(root) if root is not None else output_root()
    launcher = jobs.Launcher(root=root, repo=Path(repo) if repo else repo_root(),
                             profile=profile, enabled=bool(allow_submit))

    app = FastAPI(
        title="cyclicnano",
        description="Read-only browser for cyclicnano design runs.",
        docs_url="/api/docs",
        openapi_url="/api/openapi.json",
    )
    app.state.root = root

    @app.exception_handler(store.NotFound)
    def _not_found(request, exc):
        return JSONResponse({"detail": str(exc)}, status_code=404)

    @app.exception_handler(store.Forbidden)
    def _forbidden(request, exc):
        return JSONResponse({"detail": str(exc)}, status_code=403)

    @app.exception_handler(jobs.Invalid)
    def _invalid(request, exc):
        return JSONResponse({"detail": str(exc)}, status_code=400)

    @app.exception_handler(jobs.Busy)
    def _busy(request, exc):
        # 409, not 400: the form was fine, the machine was not free.
        return JSONResponse({"detail": str(exc)}, status_code=409)

    @app.exception_handler(jobs.Disabled)
    def _disabled(request, exc):
        return JSONResponse({"detail": str(exc)}, status_code=403)

    @app.exception_handler(projects.Invalid)
    def _bad_project(request, exc):
        return JSONResponse({"detail": str(exc)}, status_code=400)

    @app.exception_handler(projects.Exists)
    def _project_exists(request, exc):
        # 409, not 400: the form was fine, the name was taken.
        return JSONResponse({"detail": str(exc)}, status_code=409)

    @app.get("/api/health")
    def health() -> dict:
        """Whether the output root is readable, and what this service allows.

        Worth having: the single most likely failure of this service is a root
        that is not mounted. `can_submit` is here so the page can leave the
        submit form out entirely rather than offering a button that 403s.
        """
        return {"ok": root.is_dir(), "root": str(root),
                "can_submit": launcher.enabled, "version": _installed_version()}

    @app.get("/api/runs")
    def runs() -> list[dict]:
        return store.list_runs(root)

    @app.post("/api/runs")
    def create_run(form: dict = Body(default_factory=dict)) -> dict:
        """Create a project: the directory and its manifest, and nothing else.

        Naming a target and spending GPU hours on it are different decisions, so
        no stage runs here. A project that exists with nothing in it says what is
        intended before anything has been computed.

        Gated by the same flag as submitting a job, because it writes.
        """
        if not launcher.enabled:
            raise jobs.Disabled("Project creation is not enabled on this service; "
                                "start it with --allow-submit")
        return projects.create(root, form)

    @app.get("/api/runs/{run}")
    def describe_run(run: str) -> dict:
        """One project: its declared metadata and its measured funnel."""
        return store.describe(root, run)

    @app.get("/api/runs/{run}/backbones")
    def backbones(run: str) -> list[dict]:
        return store.backbones(root, run)

    @app.get("/api/runs/{run}/backbones/{backbone_id}")
    def backbone(run: str, backbone_id: str) -> dict:
        return store.backbone(root, run, backbone_id)

    @app.get("/api/runs/{run}/file")
    def file(run: str,
             path: str = Query(..., description="Path inside 04_delivery/"),
             download: bool = Query(False, description="Send as a named attachment")):
        """One structure or score file.

        Two callers with different needs, so the flag rather than two routes.
        The viewer fetches the bytes and parses them itself, which wants
        text/plain; a save wants a filename and a type the browser will not
        rename. Serving a PDB as text/plain to a download is what produced
        `.txt` files: the URL carries no filename, so the browser invents one
        from the media type.

        Caching is off either way, because a re-run overwrites these paths in
        place and a stale structure beside fresh numbers is worse than a refetch.
        """
        target = store.resolve(root, run, path)
        headers = {"Cache-Control": "no-store"}
        if download:
            # `filename` is what sets Content-Disposition, and it encodes a
            # non-ASCII name properly, which a hand-written header would not.
            return FileResponse(target, media_type="application/octet-stream",
                                headers=headers, filename=target.name)
        return FileResponse(target, media_type="text/plain", headers=headers)

    @app.get("/api/runs/{run}/export")
    def export(run: str,
               fmt: str = Query("fasta", alias="format", description="fasta or csv"),
               scope: str = Query("delivered", description="|".join(store.SCOPES))):
        """The sequences of one run, as a file to send on.

        Always an attachment: this route exists to hand a sheet to whoever
        orders the synthesis, and a FASTA rendered in a browser tab is a step
        away from being pasted by hand.

        The exported `dG_interface` column keeps the name the pipeline stores it
        under, not the label the table shows. A file that outlives this page has
        to stay traceable to the field it came from.
        """
        rows = store.export_rows(root, run, scope)
        symmetry = store.funnel(root, run).get("symmetry", "")

        if fmt == "csv":
            body, media, suffix = store.export_csv(rows), "text/csv", "csv"
        elif fmt == "fasta":
            # No registered media type for FASTA, and text/* invites the browser
            # to rename the file after its type. Octet-stream keeps the name.
            body, media, suffix = (store.export_fasta(rows, symmetry),
                                   "application/octet-stream", "fasta")
        else:
            raise store.Forbidden(f"not an export format: {fmt!r}")

        name = f"{run}_sequences_{scope}.{suffix}"
        return Response(content=body.encode("utf-8"), media_type=media,
                        headers={"Cache-Control": "no-store",
                                 "Content-Disposition": f'attachment; filename="{name}"'})

    @app.get("/api/runs/{run}/assembly")
    def assembly_pdb(run: str,
                     path: str = Query(..., description="A prediction inside 04_delivery/"),
                     design: str = Query(..., description="The design backbone of its folder"),
                     download: bool = Query(False, description="Send as a named attachment")):
        """The whole designed ring for one sequence, built here and not stored.

        The pipeline predicts one subunit per sequence and never the oligomer, so
        this is the prediction laid onto every chain of the design backbone - the
        designed geometry carrying the side chains the sequence implies, which is
        the closest thing to a pentamer this pipeline has. It is assembled per
        request because it is derivable in milliseconds and a file per sequence
        would be tens of megabytes of the same thing.

        Both paths come from the browser and both go through the same resolution
        as any other file, so neither can leave the delivery folder. They are
        pairings the caller could have made anyway by downloading the two files.
        """
        target = store.resolve(root, run, path)
        backbone_pdb = store.resolve(root, run, design)
        for chosen in (target, backbone_pdb):
            if chosen.suffix.lower() != ".pdb":
                raise store.Forbidden(f"not a structure: {chosen.name}")

        try:
            built = place_monomer(backbone_pdb, target)
        except ValueError as err:
            raise store.Forbidden(f"cannot assemble: {err}") from err

        name = f"{target.stem}_C{len(built['chains'])}.pdb"
        headers = {"Cache-Control": "no-store"}
        if download:
            headers["Content-Disposition"] = f'attachment; filename="{name}"'
        return Response(content=built["pdb"],
                        media_type=("application/octet-stream" if download
                                    else "text/plain"),
                        headers=headers)

    @app.get("/api/runs/{run}/interface")
    def interface_residues(run: str,
                           path: str = Query(..., description="Path inside 04_delivery/"),
                           cutoff: float = Query(interface.DEFAULT_CUTOFF, gt=0,
                                                 le=interface.MAX_CUTOFF)):
        """Which residues of a dimer are at the interface.

        Stage 05 records how good the interface is, not which residues make it,
        so this is computed from the structure the energy was measured on. The
        path goes through the same resolution as any other file, so it cannot
        leave the delivery folder.
        """
        target = store.resolve(root, run, path)
        if target.suffix.lower() != ".pdb":
            raise store.Forbidden(f"not a structure: {target.name}")
        return interface.of_file(target, cutoff)

    @app.post("/api/clientlog")
    def client_log(payload: dict = Body(default_factory=dict)) -> dict:
        """What the browser could not tell anyone, written to the server log.

        This service is developed on a machine that cannot open the page, so when
        the viewer fails to draw there is no console to read. The page gathers
        what would have been in it - whether Mol* loaded, whether WebGL exists,
        the host element's measured size, the caught exception - and posts it
        here, where it can be read without asking anyone to copy text out of a
        browser.

        It writes to the log and nowhere else, and the payload is truncated: an
        unauthenticated endpoint that appends to a file is a way to fill a disk.
        """
        body = json.dumps(payload, ensure_ascii=False, default=str)[:8000]
        print(f"[webui] client report: {body}", flush=True)
        return {"logged": True}

    @app.get("/api/jobs")
    def job_list() -> dict:
        """The queue, such as it is: one slot, plus the recent history.

        The form spec rides along so the browser does not need a second request
        and cannot disagree with the server about the ranges.
        """
        return {"enabled": launcher.enabled, "spec": jobs.form_spec(),
                "running": jobs.running(root), "jobs": jobs.listing(root)}

    @app.post("/api/jobs")
    def job_submit(form: dict = Body(default_factory=dict)) -> dict:
        return launcher.submit(form)

    @app.get("/api/jobs/{job_id}")
    def job_get(job_id: str) -> dict:
        return jobs.get(root, job_id)

    @app.get("/api/jobs/{job_id}/log")
    def job_log(job_id: str, offset: int = Query(0, ge=0)) -> dict:
        """The log from a byte offset, so the page can follow a running job."""
        return jobs.log_tail(root, job_id, offset)

    @app.post("/api/jobs/{job_id}/cancel")
    def job_cancel(job_id: str) -> dict:
        return launcher.cancel(job_id)

    if not STATIC.is_dir():                            # pragma: no cover - packaging
        raise RuntimeError(f"static assets are missing: {STATIC}")

    @app.get("/")
    def index():
        return FileResponse(STATIC / "index.html", headers={"Cache-Control": "no-store"})

    app.mount("/static", _PageFiles(directory=STATIC), name="static")
    return app


def main(argv: list[str] | None = None) -> int:
    import argparse

    import uvicorn

    parser = argparse.ArgumentParser(
        prog="python -m cyclicnano.webui",
        description="Serve the read-only browser for cyclicnano design runs.")
    parser.add_argument("--root", default=None,
                        help=f"output root to browse (default: ${ENV_ROOT} or {DEFAULT_ROOT})")
    parser.add_argument("--host", default="0.0.0.0",
                        help="address to bind; 127.0.0.1 to keep it on this machine")
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument("--reload", action="store_true", help="reload on code changes")
    parser.add_argument("--allow-submit", action="store_true",
                        help="let the page start pipeline runs on this machine. "
                             "There is no authentication, so anyone who can reach "
                             "the port can spend GPU hours; off by default.")
    parser.add_argument("--repo", default=None,
                        help="the checkout holding configs/ (default: this package's)")
    parser.add_argument("--profile", default=None,
                        help="machine profile passed to every submitted run, "
                             "e.g. configs/profiles/amax.yaml")
    args = parser.parse_args(argv)

    root = Path(args.root) if args.root else output_root()
    if not root.is_dir():
        parser.error(f"output root does not exist: {root}")

    repo = Path(args.repo) if args.repo else repo_root()
    if args.allow_submit and not (repo / "configs" / "runs" / "base.yaml").is_file():
        parser.error(f"--allow-submit needs configs/ under {repo}; pass --repo")

    print(f"[webui] serving {root}")
    print(f"[webui] submit: {'on' if args.allow_submit else 'off'}"
          + (f", repo {repo}" if args.allow_submit else ""))
    if args.allow_submit:
        print("[webui] WARNING: no authentication. Anyone who can reach this port "
              "can start a run on this machine.")
    print(f"[webui] http://{args.host}:{args.port}/")

    if args.reload:
        os.environ[ENV_ROOT] = str(root)
        uvicorn.run("cyclicnano.webui.app:create_app", factory=True,
                    host=args.host, port=args.port, reload=True)
    else:
        uvicorn.run(create_app(root, allow_submit=args.allow_submit, repo=repo,
                               profile=args.profile),
                    host=args.host, port=args.port)
    return 0


if __name__ == "__main__":                             # pragma: no cover
    raise SystemExit(main())
