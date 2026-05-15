# kcnq3-registry

A federated, privacy-preserving registry of quantitative EEG findings for
families affected by rare-variant KCNQ-spectrum epilepsies (KCNQ2, KCNQ3,
KCNQ5) and related channelopathies.

The registry exists so that n=1 case reports can become n=many cohort
data — accelerating recognition of variant-specific EEG signatures,
treatment-response patterns, and developmental trajectories that
otherwise take decades to surface in the published literature.

## How it works

1. A family runs [KCNQ3-Lens](https://github.com/your-org/kcnq3-lens)
   on a recording.
2. The app produces a quantitative findings JSON (locally, no upload).
3. The family clicks **Contribute** → the de-identification builder
   constructs a registry-shaped submission (no PHI, by construction).
4. The app opens a pre-filled GitHub PR adding **one JSON line** to
   `data/registry.jsonl`.
5. CI validates the line against schema v1 and runs a PHI-pattern scan.
6. A maintainer reviews the PR and merges. CI rebuilds
   `releases/v1/aggregates.json` (k≥5 per cohort cell).
7. KCNQ3-Lens downloads aggregates and shows peer-comparison
   percentiles in subsequent recordings.

There is no backend. There is no database server. The "database" is a
public, version-controlled JSONL file. Every submission is auditable.

## Privacy model

- **Allowlist by construction.** Submissions are built from a fixed
  schema. Fields not in the schema cannot appear in the output —
  there is no scrub pass that might miss something.
- **No exact dates.** Submission month only, never per-day timestamps.
- **No exact ages.** Bucketed (e.g., "5-7" not "5y 4m").
- **No exact durations.** Bucketed.
- **No filenames, paths, or free-text labels.**
- **PHI-pattern scan.** Belt-and-suspenders regex sweep on every
  submission — ISO dates, emails, MRN-like numbers, paths,
  name-like sequences, narrative text. CI rejects matches.
- **k-anonymity guard.** Aggregates only publish cells with n ≥ 5.
  Cells below threshold are merged with the next-coarser bucket or
  suppressed. Prevents identifying "the only Italian 4-year-old".
- **Right to erasure.** Families keep their `submission_id` locally.
  To withdraw: open an issue with the id; the maintainer removes the
  line in the next merge.

See `data/consent_v1.md` for the consent text families must affirm.

## Schema

Defined in code in [`scripts/_lib_schema.py`](scripts/_lib_schema.py).
A human-readable summary lives in [`schema/v1.md`](schema/v1.md).

## Running the tools locally

```bash
# Validate every line in data/registry.jsonl
python3 scripts/validate_registry.py

# Validate a single submission JSON
python3 scripts/validate_submission.py path/to/submission.json

# Rebuild aggregates (produces releases/v1/aggregates.json)
python3 scripts/build_aggregates.py

# Run aggregator + validator tests
python3 scripts/test_aggregator.py
```

No dependencies beyond the Python stdlib.

## Governance

This registry is a community resource. The maintainer's only role is
to review PRs for PHI leakage. Maintainers do NOT:
- accept submissions privately (only through public PRs)
- modify submitted records (errors are corrected by appending a
  superseding record)
- delete records except on documented right-to-erasure requests

## License

- **Code** (scripts/, workflows): MIT. See `LICENSE`.
- **Data** (`data/registry.jsonl`, `releases/`): CC0 1.0 (public
  domain dedication) — contributors waive copyright so the data can
  be combined with other open datasets without friction.

## Not medical advice

This registry is a research aggregation tool. It is not a medical
device. Findings here do not constitute diagnosis or treatment
recommendations.
