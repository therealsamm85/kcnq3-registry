#!/usr/bin/env python3
"""Build cohort aggregates with k-anonymity guard.

Reads `data/registry.jsonl`, partitions submissions into cohort cells,
suppresses cells with n < K_MIN, and writes
`releases/v1/aggregates.json` for the desktop app to download.

Cell hierarchy (finest → coarsest)
----------------------------------
1. (gene, variant_protein, age_bucket, sex)
2. (gene, variant_protein, age_bucket)
3. (gene, variant_protein)
4. (gene, age_bucket)
5. (gene)

We publish EVERY cell at every level that meets n ≥ K_MIN. The desktop
app's lookup picks the finest cell available for a given child.

Stats per cell
--------------
For each quantitative finding field present in ≥ K_MIN submissions:
- n_with_value
- mean, sd, median, p10, p25, p75, p90
- min, max (only published if n ≥ 10 — extremes are more identifying)

For categorical fields (csws_criterion_met, spindle_interpretation,
activation_label, quality_grade) we publish the count distribution.

K_MIN = 5
---------
A standard k-anonymity threshold for rare-disease registries.
Configurable via env var REGISTRY_K_MIN for testing.
"""

from __future__ import annotations

import json
import math
import os
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts._lib_validate import validate_submission
else:
    from ._lib_validate import validate_submission


K_MIN = int(os.environ.get("REGISTRY_K_MIN", "5"))
K_MIN_FOR_EXTREMES = 10

REPO_ROOT = Path(__file__).resolve().parent.parent
REGISTRY_PATH = REPO_ROOT / "data" / "registry.jsonl"
OUTPUT_PATH = REPO_ROOT / "releases" / "v1" / "aggregates.json"


# ── Field catalog ──────────────────────────────────────────────────────────

# Quantitative findings — we publish numeric distributions for these.
QUANT_FIELDS: tuple[str, ...] = (
    "background_pdr_hz",
    "csws_threshold_pct",
    "spindle_density_per_min_central",
    "activation_factor",
    "morphology_events_per_min",
    "morphology_spike_wave_pct",
    "n_sleep_cycles",
)

# Per-stage percentage dicts — flatten to "{base}__{STAGE}".
STAGE_FIELDS: tuple[str, ...] = (
    "swi_pct_by_stage",
    "sleep_stages_pct",
)
STAGE_KEYS = ("WAKE", "N1", "N2", "N3", "REM")

# Categorical (count distributions).
CATEGORICAL_FIELDS: tuple[str, ...] = (
    "csws_criterion_met",        # bool → "true"/"false"
    "spindle_interpretation",
    "activation_label",
    "quality_grade",
)


# ── Cell key extraction ────────────────────────────────────────────────────

def _cell_keys(sub: dict) -> list[tuple]:
    """Return the list of cell keys this submission belongs to, finest first."""
    s = sub.get("subject", {})
    gene = s.get("variant_gene", "?")
    protein = s.get("variant_protein", "?")
    age = s.get("age_years_bucket", "?")
    sex = s.get("sex", "?")
    return [
        ("gene_protein_age_sex", gene, protein, age, sex),
        ("gene_protein_age",      gene, protein, age),
        ("gene_protein",          gene, protein),
        ("gene_age",              gene, age),
        ("gene",                  gene),
    ]


def _cell_key_to_dict(key: tuple) -> dict:
    """Render a cell key as a JSON-friendly dict."""
    level = key[0]
    rest = key[1:]
    fields_by_level = {
        "gene_protein_age_sex":
            ["variant_gene", "variant_protein", "age_years_bucket", "sex"],
        "gene_protein_age":
            ["variant_gene", "variant_protein", "age_years_bucket"],
        "gene_protein": ["variant_gene", "variant_protein"],
        "gene_age":     ["variant_gene", "age_years_bucket"],
        "gene":         ["variant_gene"],
    }
    fields = fields_by_level[level]
    out = {"level": level}
    for f, v in zip(fields, rest):
        out[f] = v
    return out


# ── Stat helpers ───────────────────────────────────────────────────────────

def _percentile(sorted_vals: list[float], p: float) -> float:
    """Linear-interpolated percentile (p in [0, 100])."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    k = (len(sorted_vals) - 1) * (p / 100.0)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return sorted_vals[int(k)]
    return sorted_vals[f] + (sorted_vals[c] - sorted_vals[f]) * (k - f)


def _summarize_numeric(values: list[float]) -> dict[str, float | int] | None:
    """Return a stat dict, or None if n < K_MIN."""
    vals = [v for v in values
            if isinstance(v, (int, float)) and not math.isnan(v)
            and not math.isinf(v)]
    n = len(vals)
    if n < K_MIN:
        return None
    vals_sorted = sorted(vals)
    out: dict[str, float | int] = {
        "n": n,
        "mean": round(statistics.fmean(vals), 4),
        "median": round(_percentile(vals_sorted, 50), 4),
        "sd": round(statistics.pstdev(vals), 4) if n >= 2 else 0.0,
        "p10": round(_percentile(vals_sorted, 10), 4),
        "p25": round(_percentile(vals_sorted, 25), 4),
        "p75": round(_percentile(vals_sorted, 75), 4),
        "p90": round(_percentile(vals_sorted, 90), 4),
    }
    if n >= K_MIN_FOR_EXTREMES:
        out["min"] = round(min(vals), 4)
        out["max"] = round(max(vals), 4)
    return out


def _summarize_categorical(values: Iterable[Any]) -> dict[str, int] | None:
    counts: dict[str, int] = defaultdict(int)
    for v in values:
        if isinstance(v, bool):
            counts["true" if v else "false"] += 1
        elif v is None:
            continue
        else:
            counts[str(v)] += 1
    total = sum(counts.values())
    if total < K_MIN:
        return None
    return dict(sorted(counts.items()))


# ── Main aggregation ───────────────────────────────────────────────────────

def aggregate(submissions: list[dict], k_min: int | None = None) -> dict:
    """Compute the full aggregates payload.

    `k_min` overrides the module-level K_MIN if given (used in tests).
    """
    global K_MIN  # used by helpers in this module
    if k_min is not None:
        K_MIN = k_min

    # Partition submissions into all the cell hierarchies they fit.
    cells: dict[tuple, list[dict]] = defaultdict(list)
    for sub in submissions:
        for key in _cell_keys(sub):
            cells[key].append(sub)

    out_cells: list[dict] = []
    for key, members in cells.items():
        n = len(members)
        if n < K_MIN:
            continue  # k-anonymity guard

        cell_summary: dict = {
            "cell": _cell_key_to_dict(key),
            "n": n,
            "stats": {},
            "categorical": {},
        }

        # Quantitative
        for fld in QUANT_FIELDS:
            vals = [m["findings"][fld]
                    for m in members
                    if fld in (m.get("findings") or {})]
            summ = _summarize_numeric(vals)
            if summ is not None:
                cell_summary["stats"][fld] = summ

        # Per-stage dicts
        for base in STAGE_FIELDS:
            for stage in STAGE_KEYS:
                vals = []
                for m in members:
                    d = (m.get("findings") or {}).get(base) or {}
                    if stage in d:
                        vals.append(d[stage])
                summ = _summarize_numeric(vals)
                if summ is not None:
                    cell_summary["stats"][f"{base}__{stage}"] = summ

        # Categorical
        for fld in CATEGORICAL_FIELDS:
            vals = [m["findings"][fld]
                    for m in members
                    if fld in (m.get("findings") or {})]
            summ = _summarize_categorical(vals)
            if summ is not None:
                cell_summary["categorical"][fld] = summ

        # Suppress empty cells (nothing met threshold inside them)
        if not cell_summary["stats"] and not cell_summary["categorical"]:
            continue

        out_cells.append(cell_summary)

    # Sort: finest level first, then by cell key for stability.
    level_order = {
        "gene_protein_age_sex": 0, "gene_protein_age": 1,
        "gene_protein": 2, "gene_age": 3, "gene": 4,
    }
    out_cells.sort(key=lambda c: (
        level_order.get(c["cell"]["level"], 99),
        json.dumps(c["cell"], sort_keys=True),
    ))

    return {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"
        ),
        "k_min": K_MIN,
        "n_submissions": len(submissions),
        "n_cells_published": len(out_cells),
        "cells": out_cells,
    }


def _load_jsonl(path: Path) -> list[dict]:
    if not path.exists():
        return []
    out = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            ok, _ = validate_submission(obj)
            if ok:
                out.append(obj)
    return out


def main() -> int:
    subs = _load_jsonl(REGISTRY_PATH)
    agg = aggregate(subs)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(agg, indent=2, sort_keys=False))
    print(f"  registry: {len(subs)} valid submissions")
    print(f"  cells published (n ≥ {K_MIN}): {agg['n_cells_published']}")
    print(f"  output: {OUTPUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
