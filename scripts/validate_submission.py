#!/usr/bin/env python3
"""Validate one or more submission JSON files.

Usage:
    python3 scripts/validate_submission.py path/to/sub1.json [sub2.json ...]

Exits 0 if all valid, 1 otherwise. Prints per-file pass/fail and the
human-readable error list for any failure.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Allow running both as `python3 scripts/validate_submission.py` and
# `python3 -m scripts.validate_submission`.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts._lib_validate import validate_submission
else:
    from ._lib_validate import validate_submission


def main(argv: list[str]) -> int:
    if not argv:
        print("usage: validate_submission.py <file.json> [file2.json ...]",
              file=sys.stderr)
        return 2

    n_ok = 0
    n_fail = 0
    for arg in argv:
        path = Path(arg)
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as e:
            print(f"\x1b[31m✗\x1b[0m {path}: cannot read/parse: {e}")
            n_fail += 1
            continue
        ok, errors = validate_submission(obj)
        if ok:
            print(f"\x1b[32m✓\x1b[0m {path}: valid")
            n_ok += 1
        else:
            print(f"\x1b[31m✗\x1b[0m {path}: INVALID")
            for err in errors:
                print(f"    {err}")
            n_fail += 1

    print(f"\n{n_ok} valid, {n_fail} invalid")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
