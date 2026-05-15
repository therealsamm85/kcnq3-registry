"""PHI pattern scanner — belt-and-suspenders defense.

The submission builder's allowlist is the primary defense: fields not
in `schema.SCHEMA_VERSION 1` simply don't exist in the output. This
module is the second layer: it scans every string VALUE in a built
submission for patterns that look like PHI (dates, paths, emails,
names-with-numbers, free-form sentences) and flags them.

A flagged submission MUST NOT be uploaded. The UI should refuse, the
CI in the registry repo should reject the PR, and `validate_submission`
calls `scan_for_phi` before returning success.

Patterns checked
- Exact ISO date 'YYYY-MM-DD' (we only allow 'YYYY-MM' anywhere)
- DMY/MDY date strings with dots or slashes
- Email addresses
- Phone-like number sequences (7+ consecutive digits, or known formats)
- File path separators ('/', '\\\\') with multiple segments
- Long free-text strings (> 200 chars) — there is no schema field that
  requires that much text; any such string is suspicious
- Words that look like person names (Capitalized + space + Capitalized
  sequences) — heuristic, can false-positive on legitimate place names,
  but every flag is reviewed manually before merge

This list is conservative on purpose: we'd rather block a legitimate
submission and have the family adjust than leak something.
"""

from __future__ import annotations

import re
from typing import Any


# ── Patterns ────────────────────────────────────────────────────────────────

# YYYY-MM-DD (full date). The schema only allows YYYY-MM at the top level;
# anything else is suspicious.
_PAT_ISO_DATE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")

# DD.MM.YYYY (EU) or DD/MM/YYYY or MM/DD/YYYY
_PAT_NUMERIC_DATE = re.compile(
    r"\b\d{1,2}[./-]\d{1,2}[./-]\d{2,4}\b"
)

# Email
_PAT_EMAIL = re.compile(
    r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}"
)

# 7+ consecutive digits (MRNs, long IDs, phone-without-separator)
_PAT_LONG_NUMBER = re.compile(r"\b\d{7,}\b")

# Phone-like with explicit phone markers — must include either '+' or
# whitespace or parens. Pure dashed-digit sequences (uuid hex blocks,
# date-like all-digit ids) are intentionally NOT matched here; they're
# caught by the more specific date / long-number patterns instead.
_PAT_PHONE = re.compile(
    r"(?:\+\d{1,3}[\s.-]?\d{2,4}[\s.-]?\d{3,4}[\s.-]?\d{3,4}"
    r"|\(\d{2,4}\)[\s.-]?\d{3,4}[\s.-]?\d{3,4}"
    r"|\b\d{2,4}\s\d{3,4}\s\d{3,4}\b)"
)

# Multi-segment path
_PAT_PATH = re.compile(r"(?:[A-Za-z]:[\\/]|/)[\w./\\-]{4,}")

# Two consecutive Capitalized words (heuristic person name).
# Excludes ALL-CAPS gene symbols which would not match (only first letter cap).
_PAT_NAME_LIKE = re.compile(
    r"\b[A-Z][a-z]{1,}\s+[A-Z][a-z]{1,}\b"
)

# Common name-y connectors in narrative free text
_PAT_NARRATIVE = re.compile(
    r"\b(?:my (?:son|daughter|child)|patient|named|called)\b",
    re.IGNORECASE,
)

# Hard cap for any string value in the submission. The longest legitimate
# string in v1 schema is intervention.name (capped at 64). Use 80 as the
# scanner threshold to leave margin.
_MAX_STRING_LEN = 80

# Paths that are already validated as opaque IDs (uuid4 regex in the
# schema). The PHI scanner skips them — a uuid will contain digit runs
# that look phone-like or MRN-like otherwise, and the dedicated regex
# is the real authority on that field.
_SKIP_PATHS = frozenset({
    "$.submission_id",
    "$.intervention.linked_pre_submission_id",
})


def scan_for_phi(obj: Any, *, path: str = "$") -> list[str]:
    """Walk a submission JSON and return a list of human-readable findings.

    Each finding is a one-line description of where a suspicious value
    sits. Empty list = clean.

    The walk is structural: dicts, lists, tuples, and primitives only.
    We don't follow arbitrary attributes — the input is always a plain
    JSON-shaped object.
    """
    findings: list[str] = []

    if isinstance(obj, dict):
        for k, v in obj.items():
            findings.extend(scan_for_phi(v, path=f"{path}.{k}"))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            findings.extend(scan_for_phi(v, path=f"{path}[{i}]"))
    elif isinstance(obj, str):
        if path not in _SKIP_PATHS:
            findings.extend(_scan_string(obj, path))
    # Numbers, bools, None: no PHI risk.

    return findings


def _scan_string(s: str, path: str) -> list[str]:
    found: list[str] = []

    if len(s) > _MAX_STRING_LEN:
        found.append(
            f"{path}: string longer than {_MAX_STRING_LEN} chars "
            f"({len(s)} chars) — free-text not allowed"
        )

    if _PAT_ISO_DATE.search(s):
        found.append(f"{path}: contains ISO date 'YYYY-MM-DD'")
    if _PAT_NUMERIC_DATE.search(s):
        found.append(f"{path}: contains numeric date (DMY/MDY)")
    if _PAT_EMAIL.search(s):
        found.append(f"{path}: contains email-like pattern")
    if _PAT_LONG_NUMBER.search(s):
        found.append(f"{path}: contains 7+ consecutive digits")
    if _PAT_PHONE.search(s):
        found.append(f"{path}: contains phone-like pattern")
    if _PAT_PATH.search(s):
        found.append(f"{path}: contains file path pattern")
    if _PAT_NAME_LIKE.search(s):
        found.append(f"{path}: contains 'Capitalized Capitalized' "
                     f"sequence (name-like)")
    if _PAT_NARRATIVE.search(s):
        found.append(f"{path}: contains narrative free text "
                     f"('my child' / 'patient' / ...)")

    return found


def is_clean(obj: Any) -> bool:
    """True iff scan_for_phi(obj) returns no findings."""
    return not scan_for_phi(obj)
