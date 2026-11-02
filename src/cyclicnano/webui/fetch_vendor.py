"""Fetch the browser libraries the viewer needs.

Mol* is not committed: it is five megabytes of minified bundle, which does not
belong in a research repository. It is fetched once, pinned by version and
checked against a hash, into `static/vendor/`, which git ignores.

The page falls back to the CDN when a file is missing, so this is an
optimisation for a lab machine whose browsers may not reach the open internet,
not a build step. Run it once after cloning:

    python -m cyclicnano.webui.fetch_vendor
"""
from __future__ import annotations

import hashlib
import sys
import urllib.request
from pathlib import Path

VERSION = "5.13.1"

# Two mirrors, because on this network cdnjs is unreachable while these two are
# not. Tried in order; the hashes make it safe not to care which one answers.
MIRRORS = (
    "https://cdn.jsdelivr.net/npm/molstar@{version}/build/viewer/{name}",
    "https://unpkg.com/molstar@{version}/build/viewer/{name}",
)

ASSETS = {
    "molstar.js": "4fc6493c2471db8834dcc391151b344d7017e4eb8a4eb76c152be5f9dbd23874",
    "molstar.css": "22fbf73b75c11ef57687d08437c95dbfe8430ea895264f6da2b9b01e28bacc49",
}

VENDOR = Path(__file__).resolve().parent / "static" / "vendor"


def fetch(name: str, digest: str, force: bool = False) -> bool:
    target = VENDOR / name

    if target.is_file() and not force:
        if hashlib.sha256(target.read_bytes()).hexdigest() == digest:
            print(f"[vendor] {name} is already correct")
            return True
        print(f"[vendor] {name} has the wrong hash, refetching")

    errors = []
    for mirror in MIRRORS:
        url = mirror.format(version=VERSION, name=name)
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                blob = response.read()
        except Exception as exc:
            errors.append(f"    {url}: {type(exc).__name__}: {exc}")
            continue

        got = hashlib.sha256(blob).hexdigest()
        if got != digest:
            errors.append(f"    {url}: hash mismatch, got {got}")
            continue

        VENDOR.mkdir(parents=True, exist_ok=True)
        target.write_bytes(blob)
        print(f"[vendor] wrote {name} ({len(blob):,} bytes) from {url}")
        return True

    print(f"[vendor] could not fetch {name}; the page will fall back to the CDN")
    print("\n".join(errors))
    return False


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    force = "--force" in argv
    ok = all(fetch(name, digest, force) for name, digest in ASSETS.items())

    # A bundle left over from an earlier version of this viewer. Removing it
    # keeps the vendor directory from quietly growing a second copy of a 3D
    # engine that nothing loads any more.
    stale = VENDOR / "3Dmol-min.js"
    if stale.is_file():
        stale.unlink()
        print(f"[vendor] removed {stale.name}, no longer used")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
