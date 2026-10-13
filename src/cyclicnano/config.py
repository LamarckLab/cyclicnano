"""Config loading: default.yaml, optionally overlaid with a machine profile and CLI overrides."""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "configs" / "default.yaml"

SYMMETRY_ORDERS = {f"C{i}": i for i in range(2, 13)}


def _is_absolute(path: Path) -> bool:
    """Absolute on the host, or a POSIX path written for the machine that will run it.

    Path.is_absolute() is False for "/data/lmk/runs" on Windows, since it carries no
    drive, so a config authored for a Linux host would otherwise be resolved against
    the repository root when inspected from a workstation.
    """
    # as_posix(), not str(): on Windows str() renders the separators as backslashes,
    # so the leading slash of a POSIX path would not be visible.
    return path.is_absolute() or path.as_posix().startswith("/")


def _deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


class Config:
    """Dotted read access over the merged config dict."""

    def __init__(self, data: dict, sources: list | None = None, overrides: dict | None = None):
        self.data = data
        self.sources = sources or []                        # files merged, in precedence order
        self.overrides = overrides or {}                    # command-line --set values

    def resolved_yaml(self) -> str:
        """The merged config, annotated with where it came from.

        Written into every run so the directory records its own recipe: a run made
        with --set overrides is otherwise only reconstructable from shell history.
        """
        lines = ["# Fully merged configuration, written when the run started.",
                 "# Later sources override earlier ones.", "#", "# Sources:"]
        lines += [f"#   {src}" for src in self.sources]
        if self.overrides:
            lines.append(f"#   command line: {self.overrides}")
        body = yaml.safe_dump(self.data, sort_keys=False, allow_unicode=True,
                              default_flow_style=False)
        return "\n".join(lines) + "\n\n" + body

    def get(self, dotted: str, default: Any = None) -> Any:
        node: Any = self.data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def require(self, dotted: str) -> Any:
        value = self.get(dotted, None)
        if value is None:
            raise KeyError(f"config is missing required key: {dotted}")
        return value

    # -- derived values -----------------------------------------------------
    @property
    def symmetry(self) -> str:
        return str(self.require("target.symmetry")).upper()

    @property
    def sym_order(self) -> int:
        sym = self.symmetry
        if sym not in SYMMETRY_ORDERS:
            raise ValueError(f"unsupported symmetry {sym!r}; expected one of {sorted(SYMMETRY_ORDERS)}")
        return SYMMETRY_ORDERS[sym]

    @property
    def monomer_length(self) -> int:
        return int(self.require("target.monomer_length"))

    @property
    def total_length(self) -> int:
        """RFdiffusion symmetric mode takes the TOTAL contig length, not per-chain.

        Getting this wrong silently produces oligomers of the wrong size, so the
        conversion lives here rather than in the caller.
        """
        return self.sym_order * self.monomer_length

    @property
    def outdir(self) -> Path:
        """Where this run writes, as `run.outdir` or else `run.outroot` / `run.name`.

        Splitting the two lets the machine profile say where results live on this
        host while the run config says only what the run is called, so switching
        between the C5 and C3 variants needs no path on the command line.
        """
        explicit = self.get("run.outdir")
        if explicit:
            path = Path(str(explicit))
        else:
            path = Path(str(self.get("run.outroot", "runs"))) / str(self.require("run.name"))
        return path if _is_absolute(path) else (REPO_ROOT / path).resolve()

    def unknown_keys(self) -> list[str]:
        """Keys the merged config carries that `default.yaml` does not define.

        default.yaml is documented as the complete parameter surface, so anything not
        in it is a typo or a misplacement, and both are silent: a key under the wrong
        section is simply never read. That cost a run. A dropped `generate:` header put
        `noise_scale` under `target:`, nothing complained, and five hundred backbones
        were generated at the default noise while the config said otherwise.

        Where default.yaml holds a non-dict or an empty dict, the subtree below is
        free-form (machine-specific env vars, image names) and is not descended into.
        """
        defaults = yaml.safe_load(DEFAULT_CONFIG.read_text(encoding="utf-8")) or {}

        def walk(got: dict, want: dict, prefix: str) -> list[str]:
            out = []
            for key, value in got.items():
                path = f"{prefix}{key}"
                if key not in want:
                    out.append(path)
                elif isinstance(value, dict) and isinstance(want[key], dict) and want[key]:
                    out += walk(value, want[key], f"{path}.")
            return out

        return sorted(walk(self.data, defaults, ""))

    def validate(self) -> list[str]:
        """Cheap startup checks; returns a list of human-readable problems."""
        problems = []
        for key in self.unknown_keys():
            problems.append(f"{key} is not a key defined in configs/default.yaml; it is "
                            f"being ignored. Check the spelling and the section it sits under")
        try:
            self.sym_order
        except ValueError as exc:
            problems.append(str(exc))
        if self.monomer_length < 20:
            problems.append(f"target.monomer_length={self.monomer_length} is implausibly short")
        if int(self.get("design.n_seq_per_backbone", 0)) < 1:
            problems.append("design.n_seq_per_backbone must be >= 1")
        if int(self.get("run.seed", 1)) == 0:
            problems.append("run.seed=0 makes ProteinMPNN pick a random seed, so the "
                            "same backbones give different sequences on every run; "
                            "use any non-zero value")
        if not self.get("design.tie_chains", True):
            problems.append("design.tie_chains=false will give each chain a different "
                            "sequence, which cannot form a homo-oligomer")
        if self.get("fold.msa_mode") != "single_sequence":
            problems.append("fold.msa_mode is not single_sequence; de novo designs have no "
                            "homologues and an MSA yields misleading confidence")
        gpus = self.get("compute.gpus")
        if gpus is not None and not isinstance(gpus, list):
            problems.append("compute.gpus must be null or a list of device indices")
        return problems


def load_config(path: str | Path | list | None = None,
                profile: str | Path | None = None,
                overrides: dict | None = None) -> Config:
    """default.yaml <- user config(s) <- machine profile <- CLI overrides.

    `path` takes several files so shared settings can live in one of them and the
    per-run differences in another, instead of being copied between near-identical
    files that then drift apart.
    """
    data = yaml.safe_load(DEFAULT_CONFIG.read_text(encoding="utf-8")) or {}
    sources = [str(DEFAULT_CONFIG)]

    if path is None:
        given: list = []
    elif isinstance(path, (str, Path)):
        given = [path]
    else:
        given = list(path)

    for extra in [*given, profile]:
        if not extra:
            continue
        p = Path(extra)
        if not p.is_absolute():
            p = REPO_ROOT / p
        data = _deep_merge(data, yaml.safe_load(p.read_text(encoding="utf-8")) or {})
        sources.append(str(p))

    if overrides:
        data = _deep_merge(data, overrides)
    return Config(data, sources=sources, overrides=overrides)


def parse_overrides(pairs: list[str]) -> dict:
    """Turn ['target.symmetry=C3', 'target.monomer_length=80'] into a nested dict."""
    out: dict = {}
    for pair in pairs or []:
        if "=" not in pair:
            raise ValueError(f"override must be key=value, got {pair!r}")
        key, raw = pair.split("=", 1)
        node = out
        parts = key.split(".")
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = yaml.safe_load(raw)
    return out
