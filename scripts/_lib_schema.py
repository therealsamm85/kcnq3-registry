"""Registry submission schema v1 — defines what families can contribute.

The schema is the privacy boundary. The submission builder constructs
JSON from scratch using ONLY keys defined here. Anything not listed is
structurally absent from the output — there is no way for a free-text
field, a filename, or an unforeseen findings key to leak.

Versioning
----------
- `SCHEMA_VERSION` is bumped any time the contract changes.
- Old submissions remain valid against their version; the aggregator
  always knows which schema_version a given record was built against.
- DO NOT add free-text fields without an explicit allowlist on the
  value-side. If a clinical concept can't be expressed via a controlled
  vocabulary, add a new enum here.

Why no JSON-Schema dependency
-----------------------------
Pure-Python definitions stay in lockstep with the validator + builder
without runtime deps. The trade-off is we hand-code the validation
logic, but we test it exhaustively (see tests/test_registry.py).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

SCHEMA_VERSION = 1


# ─── Controlled vocabularies ────────────────────────────────────────────────

# Gene symbol: 2-16 chars, uppercase alphanumeric. Matches HGNC convention.
GENE_SYMBOL_RE = re.compile(r"^[A-Z][A-Z0-9]{1,15}$")

# Variant protein notation. Matches "p.Arg230His", "p.Arg230*", "p.Val1Met",
# "p.Arg230fs", etc. Keeps the field too short and too structured to fit PHI.
VARIANT_PROTEIN_RE = re.compile(
    r"^p\.[A-Z][a-z]{2}\d{1,5}(?:[A-Z][a-z]{2}|\*|fs|del|dup)$"
)

# Variant type — kept small and clinically meaningful.
VARIANT_TYPES: frozenset[str] = frozenset({
    "missense_GoF",
    "missense_LoF",
    "missense_unknown",
    "truncating",
    "splice",
    "deletion",
    "duplication",
    "regulatory",
    "unknown",
})

# Age buckets — coarse enough that {variant, age_bucket, country}
# can't usually re-identify in cohorts of < 10. The k-anonymity guard
# in the aggregator enforces n>=5 per cell anyway.
AGE_BUCKETS: tuple[str, ...] = (
    "0-1", "1-2", "2-3", "3-5", "5-7", "7-10",
    "10-13", "13-18", "18-30", "30+",
)

DURATION_BUCKETS: tuple[str, ...] = (
    "<1", "1-4", "4-12", "12-24", "24-48", "48+",
)

SEX_VALUES: frozenset[str] = frozenset({"F", "M", "X", "unknown"})

# ISO 3166-1 alpha-2. We allow any two uppercase letters — validating
# against a full country list would not catch hostile input and adds no
# real privacy value (countries are by design coarse).
COUNTRY_RE = re.compile(r"^[A-Z]{2}$")

MONTAGE_VALUES: frozenset[str] = frozenset({
    "10-20_monopolar",
    "10-10_monopolar",
    "10-20_bipolar",
    "double_banana",
    "other_monopolar",
    "other_bipolar",
    "unknown",
})

SPINDLE_INTERPRETATIONS: frozenset[str] = frozenset({
    "below", "in", "above", "no_age_provided", "no_data",
})

ACTIVATION_LABELS: frozenset[str] = frozenset({
    "none", "mild", "moderate", "strong",
})

QUALITY_GRADES: frozenset[str] = frozenset({"A", "B", "C", "D", "unknown"})

INTERVENTION_TYPES: frozenset[str] = frozenset({
    "medication", "diet", "stimulation", "behavioral", "other",
})

INTERVENTION_RECORD_KINDS: frozenset[str] = frozenset({
    "baseline", "pre", "post", "followup",
})

# Free-text fields we accept ONLY from controlled string sets above plus:
# variant_gene (regex), variant_protein (regex), intervention.name.
# intervention.name is loose free-text because medication names vary
# globally — we cap length and pass through the PHI scan, but families
# CAN type whatever they want here. Keep that in mind.
INTERVENTION_NAME_MAX_LEN = 64

# UUID4 regex (8-4-4-4-12 lowercase hex)
UUID4_RE = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$"
)

# YYYY-MM (no day — submission month is reported coarsely)
SUBMITTED_AT_MONTH_RE = re.compile(r"^\d{4}-(?:0[1-9]|1[0-2])$")


# ─── Allowlist tree — the canonical shape of a submission ───────────────────
# Each leaf is (python_type, validator_callable_or_None).
# Validators return True/False; None means "any value of the right type".

def _is_uuid(v: object) -> bool:
    return isinstance(v, str) and bool(UUID4_RE.match(v))

def _is_month(v: object) -> bool:
    return isinstance(v, str) and bool(SUBMITTED_AT_MONTH_RE.match(v))

def _is_pct(v: object) -> bool:
    return isinstance(v, (int, float)) and 0.0 <= float(v) <= 100.0

def _is_nonneg_finite(v: object) -> bool:
    if not isinstance(v, (int, float)):
        return False
    f = float(v)
    return f >= 0.0 and f == f and f not in (float("inf"), float("-inf"))

def _is_finite(v: object) -> bool:
    if not isinstance(v, (int, float)):
        return False
    f = float(v)
    return f == f and f not in (float("inf"), float("-inf"))


# Per-stage SWI / sleep_stages dict — keys are constrained.
SLEEP_STAGE_KEYS = frozenset({"WAKE", "N1", "N2", "N3", "REM"})


@dataclass(frozen=True)
class FieldSpec:
    """One leaf in the allowlist tree."""
    required: bool
    py_type: type | tuple[type, ...]
    validator: object  # callable(v) -> bool, or None
    description: str = ""
