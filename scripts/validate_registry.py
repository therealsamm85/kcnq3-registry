#!/usr/bin/env python3
"""Validate every line in data/registry.jsonl.

Run in CI on every PR. Exits 1 if any line fails — block the merge.

Additionally checks:
- Every submission_id is unique across the file
- No duplicate (variant_gene, variant_protein, submitted_at_month,
  intervention.name) tuples — a heuristic against accidental
  double-submission of the same record. Real edge cases (true two
  pre-treatment recordings same month) can be approved manually.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts._lib_validate import validate_submission
else:
    from ._lib_validate import validate_submission


REGISTRY_PATH = Path(__file__).resolve().parent.parent / "data" / "registry.jsonl"


def main() -> int:
    if not REGISTRY_PATH.exists():
        print(f"  registry file does not exist: {REGISTRY_PATH}")
        print("  (empty registry is valid — nothing to check)")
        return 0

    n_ok = 0
    n_fail = 0
    ids_seen: Counter[str] = Counter()
    near_dups: Counter[tuple] = Counter()

    with open(REGISTRY_PATH, "r", encoding="utf-8") as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as e:
                print(f"\x1b[31m✗\x1b[0m line {line_no}: JSON parse error: {e}")
                n_fail += 1
                continue
            ok, errors = validate_submission(obj)
            if not ok:
                print(f"\x1b[31m✗\x1b[0m line {line_no}: INVALID")
                for err in errors:
                    print(f"    {err}")
                n_fail += 1
                continue

            sid = obj.get("submission_id", "")
            ids_seen[sid] += 1

            subj = obj.get("subject", {})
            it = obj.get("intervention") or {}
            near_key = (
                subj.get("variant_gene"),
                subj.get("variant_protein"),
                obj.get("submitted_at_month"),
                (it.get("name") if isinstance(it, dict) else None),
            )
            near_dups[near_key] += 1
            n_ok += 1

    # Duplicate submission_id is a hard fail.
    dups = [sid for sid, c in ids_seen.items() if c > 1]
    for sid in dups:
        print(f"\x1b[31m✗\x1b[0m duplicate submission_id: {sid}")
        n_fail += 1

    # Near-duplicates are a soft warning.
    suspicious = [k for k, c in near_dups.items() if c > 1]
    for k in suspicious:
        print(f"\x1b[33m⚠\x1b[0m near-duplicate "
              f"(gene/protein/month/intervention): {k}")

    print(f"\n{n_ok} valid lines, {n_fail} invalid, "
          f"{len(suspicious)} near-duplicate warnings")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
