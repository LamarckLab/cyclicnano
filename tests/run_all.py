"""Dependency-free test runner.

The test modules are ordinary pytest files; this runner exists so the suite can also
be executed on a machine without pytest installed.

    python tests/run_all.py
"""
from __future__ import annotations

import importlib
import inspect
import sys
import traceback
from pathlib import Path

import conftest  # noqa: F401  (sets up sys.path as a side effect)

MODULES = ["test_pdbio", "test_geometry", "test_filters", "test_config",
           "test_manifest", "test_pipeline", "test_real_data", "test_icosahedral"]


def main() -> int:
    passed, failed = 0, []
    for name in MODULES:
        module = importlib.import_module(name)
        print(f"\n{name}")
        for fn_name, fn in sorted(vars(module).items()):
            if not fn_name.startswith("test_") or not inspect.isfunction(fn):
                continue
            if inspect.signature(fn).parameters:
                continue                                    # fixture-dependent, pytest only
            try:
                fn()
                passed += 1
                print(f"  [ok  ] {fn_name}")
            except Exception:
                failed.append(f"{name}.{fn_name}")
                print(f"  [FAIL] {fn_name}")
                traceback.print_exc(limit=3)
    print("\n" + "=" * 66)
    print(f"{passed} passed, {len(failed)} failed")
    for name in failed:
        print(f"  FAILED: {name}")
    print("=" * 66)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
