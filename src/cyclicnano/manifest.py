"""Run state: two tables that every stage reads and appends to.

The pipeline fans out one backbone into many sequences, so the state is kept at two
grains rather than one wide table:

    backbones.csv   one row per generated backbone      (stages 00, 01, 04-rollup)
    sequences.csv   one row per (backbone, seq_index)   (stages 02, 03, 04)

Both are written as CSV for inspection and Parquet for typed round-trips. Stages are
idempotent: re-running one overwrites its own columns and leaves the rest alone.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

BACKBONE_KEY = "backbone_id"
SEQUENCE_KEY = "sequence_id"


class Table:
    """A CSV/Parquet-backed frame keyed by a single id column."""

    def __init__(self, path: Path, key: str):
        self.path = Path(path)
        self.key = key
        self.path.parent.mkdir(parents=True, exist_ok=True)

    @property
    def parquet_path(self) -> Path:
        return self.path.with_suffix(".parquet")

    def load(self) -> pd.DataFrame:
        if self.parquet_path.exists():
            return pd.read_parquet(self.parquet_path)
        if self.path.exists():
            return pd.read_csv(self.path)
        return pd.DataFrame(columns=[self.key])

    def save(self, df: pd.DataFrame) -> None:
        df = df.sort_values(self.key).reset_index(drop=True)
        df.to_csv(self.path, index=False)
        try:
            df.to_parquet(self.parquet_path, index=False)
        except Exception:                                   # pyarrow is optional; CSV is the source of truth
            pass

    def upsert(self, rows: list[dict]) -> pd.DataFrame:
        """Merge rows in by key, overwriting only the columns each row carries.

        Columns a row does not mention are left untouched, which is what makes stages
        idempotent: re-running stage 01 rewrites its metrics without disturbing the
        sequence-level results written later.
        """
        if not rows:
            return self.load()
        new = pd.DataFrame(rows)
        if self.key not in new.columns:
            raise ValueError(f"rows are missing the key column {self.key!r}")
        new = new.drop_duplicates(subset=self.key, keep="last").set_index(self.key)

        old = self.load()
        if old.empty:
            merged = new
        else:
            old = old.set_index(self.key)
            merged = old.reindex(old.index.union(new.index))
            for col in new.columns:
                if col not in merged.columns:
                    merged[col] = pd.NA
                values = new[col]
                # None is not assignable into a float column on pandas 2; it raises
                # rather than treating it as missing. A stage writes None for a value
                # it could not compute, so a re-run after a prediction went missing
                # crashed here. NaN means the same thing and keeps the column numeric.
                if values.isna().any() and pd.api.types.is_numeric_dtype(merged[col].dtype):
                    values = pd.to_numeric(values, errors="coerce")
                merged.loc[new.index, col] = values.values
        self.save(merged.reset_index())
        return self.load()


class RunState:
    """Bundles the two tables plus the run directory layout."""

    def __init__(self, outdir: str | Path):
        self.outdir = Path(outdir)
        self.backbones = Table(self.outdir / "tables" / "backbones.csv", BACKBONE_KEY)
        self.sequences = Table(self.outdir / "tables" / "sequences.csv", SEQUENCE_KEY)

    def stage_dir(self, name: str) -> Path:
        d = self.outdir / name
        d.mkdir(parents=True, exist_ok=True)
        return d

    def done_marker(self, stage: str) -> Path:
        return self.outdir / ".done" / f"{stage}.done"

    def is_done(self, stage: str) -> bool:
        return self.done_marker(stage).exists()

    def mark_done(self, stage: str) -> None:
        marker = self.done_marker(stage)
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text("ok\n", encoding="utf-8")

    def clear_done(self, stage: str) -> None:
        self.done_marker(stage).unlink(missing_ok=True)
