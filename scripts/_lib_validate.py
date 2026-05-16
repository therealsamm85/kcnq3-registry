"""Validate a registry submission JSON.

Used in two places with identical logic:
1. The desktop app, right before showing the preview to the family
   (a built submission must validate or the build itself is buggy)
2. The registry repo's CI on every PR — each new JSONL line is
   parsed and run through this validator before merge

By colocating the validator with the builder we keep them in lockstep.
If a future schema version adds fields, validator + builder + tests
all move together.

Returns (ok: bool, errors: list[str]). Errors are human-readable.
"""

from __future__ import annotations

from typing import Any

from . import _lib_schema as _schema
from . import _lib_phi_check as phi_check


def validate_submission(obj: Any) -> tuple[bool, list[str]]:
    """Top-level entry point."""
    errors: list[str] = []
    if not isinstance(obj, dict):
        return False, ["top-level must be a JSON object"]

    # ── Required top-level keys ────────────────────────────────────────
    required = {
        "submission_id", "schema_version", "submitted_at_month",
        "consent", "subject", "recording", "findings", "tool_version",
    }
    missing = required - set(obj.keys())
    if missing:
        errors.append(f"missing required keys: {sorted(missing)}")

    # ── Reject unknown top-level keys (forward-compat strictness) ──────
    allowed_top = required | {"intervention"}
    extra = set(obj.keys()) - allowed_top
    if extra:
        errors.append(f"unknown top-level keys: {sorted(extra)}")

    # If we already have structural errors, don't try to dig deeper —
    # the messages would be noisy and unhelpful.
    if errors:
        return False, errors

    # ── Field-by-field ─────────────────────────────────────────────────
    sv = obj["schema_version"]
    # Accept schema_version 1 (legacy) and 2 (current) — additive bump.
    _VALID_SCHEMA_VERSIONS = {1, 2}
    if sv not in _VALID_SCHEMA_VERSIONS:
        errors.append(
            f"schema_version is {sv}, expected one of "
            f"{sorted(_VALID_SCHEMA_VERSIONS)}"
        )

    if not _schema._is_uuid(obj["submission_id"]):
        errors.append("submission_id is not a uuid4 string")

    if not _schema._is_month(obj["submitted_at_month"]):
        errors.append("submitted_at_month must match YYYY-MM")

    errors.extend(_validate_consent(obj["consent"]))
    errors.extend(_validate_subject(obj["subject"]))
    errors.extend(_validate_recording(obj["recording"]))
    errors.extend(_validate_findings(obj["findings"]))

    if obj.get("intervention") is not None:
        errors.extend(_validate_intervention(obj["intervention"]))

    if not isinstance(obj["tool_version"], str) or not obj["tool_version"]:
        errors.append("tool_version must be a non-empty string")

    # ── PHI scan ───────────────────────────────────────────────────────
    phi_findings = phi_check.scan_for_phi(obj)
    if phi_findings:
        errors.append("PHI scan: " + "; ".join(phi_findings))

    return (len(errors) == 0), errors


# ─── Sub-validators ────────────────────────────────────────────────────

def _validate_consent(c: Any) -> list[str]:
    errs: list[str] = []
    if not isinstance(c, dict):
        return ["consent must be an object"]
    if not isinstance(c.get("version"), int):
        errs.append("consent.version must be an int")
    if not isinstance(c.get("given"), bool):
        errs.append("consent.given must be a bool")
    elif c.get("given") is False:
        errs.append("consent.given must be True for valid submissions")
    if not _schema._is_month(c.get("given_at_month")):
        errs.append("consent.given_at_month must match YYYY-MM")
    extra = set(c.keys()) - {"version", "given", "given_at_month"}
    if extra:
        errs.append(f"consent has unknown keys: {sorted(extra)}")
    return errs


def _validate_subject(s: Any) -> list[str]:
    errs: list[str] = []
    if not isinstance(s, dict):
        return ["subject must be an object"]
    req = {"variant_gene", "variant_protein", "variant_type",
           "age_years_bucket", "sex"}
    for k in req:
        if k not in s:
            errs.append(f"subject missing required key {k!r}")
    if errs:
        return errs
    if not _schema.GENE_SYMBOL_RE.match(s["variant_gene"] or ""):
        errs.append("subject.variant_gene fails gene-symbol regex")
    if not _schema.VARIANT_PROTEIN_RE.match(s["variant_protein"] or ""):
        errs.append("subject.variant_protein fails protein-notation regex")
    if s["variant_type"] not in _schema.VARIANT_TYPES:
        errs.append("subject.variant_type not in allowed set")
    if s["age_years_bucket"] not in _schema.AGE_BUCKETS:
        errs.append("subject.age_years_bucket not in allowed set")
    if s["sex"] not in _schema.SEX_VALUES:
        errs.append("subject.sex not in allowed set")
    if "country_region" in s and s["country_region"] is not None:
        if not _schema.COUNTRY_RE.match(s["country_region"]):
            errs.append("subject.country_region must be ISO 3166-1 alpha-2")
    extra = set(s.keys()) - (
        req | {"country_region"}
    )
    if extra:
        errs.append(f"subject has unknown keys: {sorted(extra)}")
    return errs


def _validate_recording(r: Any) -> list[str]:
    errs: list[str] = []
    if not isinstance(r, dict):
        return ["recording must be an object"]
    req = {"duration_hours_bucket", "had_sleep", "montage", "n_channels"}
    for k in req:
        if k not in r:
            errs.append(f"recording missing required key {k!r}")
    if errs:
        return errs
    if r["duration_hours_bucket"] not in _schema.DURATION_BUCKETS:
        errs.append("recording.duration_hours_bucket not in allowed set")
    if not isinstance(r["had_sleep"], bool):
        errs.append("recording.had_sleep must be a bool")
    if r["montage"] not in _schema.MONTAGE_VALUES:
        errs.append("recording.montage not in allowed set")
    if not isinstance(r["n_channels"], int) or not (0 <= r["n_channels"] <= 256):
        errs.append("recording.n_channels must be int in [0, 256]")
    extra = set(r.keys()) - req
    if extra:
        errs.append(f"recording has unknown keys: {sorted(extra)}")
    return errs


def _validate_findings(f: Any) -> list[str]:
    """Every key must be in the allowlist; every value must match its
    type/range. Unknown keys are rejected outright (forward strictness).

    v1 keys and v2 keys are both included in the allowlist — any submission
    (v1 or v2) may contain either subset. New v2 keys are all optional.
    """
    errs: list[str] = []
    if not isinstance(f, dict):
        return ["findings must be an object"]

    # v1 scalar fields
    allowed_keys: dict[str, Any] = {
        "background_pdr_hz": _schema._is_finite,
        "csws_criterion_met": lambda v: isinstance(v, bool),
        "csws_threshold_pct": _schema._is_pct,
        "spindle_density_per_min_central": _schema._is_nonneg_finite,
        "spindle_interpretation": (
            lambda v: v in _schema.SPINDLE_INTERPRETATIONS
        ),
        "activation_factor": _schema._is_nonneg_finite,
        "activation_label": lambda v: v in _schema.ACTIVATION_LABELS,
        "morphology_events_per_min": _schema._is_nonneg_finite,
        "morphology_spike_wave_pct": _schema._is_pct,
        "n_sleep_cycles": (
            lambda v: isinstance(v, int) and 0 <= v <= 100
        ),
        "quality_grade": lambda v: v in _schema.QUALITY_GRADES,
        # v2 scalar fields (all optional)
        "coupling_plv_bucket": (
            lambda v: v in _schema.PLV_BUCKETS
        ),
        "coupling_preferred_phase_octant": (
            lambda v: v in _schema.PHASE_OCTANTS
        ),
        "coupling_n_events_bucket": (
            lambda v: v in _schema.COUPLED_EVENTS_BUCKETS
        ),
        "coupling_rayleigh_significant": lambda v: isinstance(v, bool),
        "sw_density_bucket": (
            lambda v: v in _schema.SW_DENSITY_BUCKETS
        ),
        "sw_mean_ptp_bucket": (
            lambda v: v in _schema.SW_PTP_BUCKETS
        ),
        "sw_method": lambda v: v in _schema.SW_METHODS,
        "hfo_rate_bucket": (
            lambda v: v in _schema.HFO_RATE_BUCKETS
        ),
        "hfo_available": lambda v: isinstance(v, bool),
        "hfo_pct_on_spike_bucket": (
            lambda v: v in _schema.HFO_PCT_ON_SPIKE_BUCKETS
        ),
        # v0.13.3 — IED detection (all optional)
        "ied_method": lambda v: v in _schema.IED_METHODS,
        "ied_rate_bucket": lambda v: v in _schema.IED_RATE_BUCKETS,
        "ied_age_flag": lambda v: v in _schema.IED_AGE_FLAGS,
        "ied_agreement_bucket": lambda v: v in _schema.IED_AGREEMENT_BUCKETS,
        "ied_n_rolandic_benign_bucket": (
            lambda v: v in _schema.IED_ROLANDIC_BUCKETS
        ),
        "ied_nrem_rate_bucket": lambda v: v in _schema.IED_NREM_RATE_BUCKETS,
    }

    # Keys handled by dedicated sub-validators (excluded from unknown-key check)
    _complex_keys = {"swi_pct_by_stage", "sleep_stages_pct", "spindle_age_norm_range"}

    extra = set(f.keys()) - (set(allowed_keys.keys()) | _complex_keys)
    if extra:
        errs.append(f"findings has unknown keys: {sorted(extra)}")

    for k, validator in allowed_keys.items():
        if k not in f:
            continue
        if not validator(f[k]):
            errs.append(f"findings.{k} failed validation: value={f[k]!r}")

    if "swi_pct_by_stage" in f:
        errs.extend(_validate_stage_pct_dict(
            f["swi_pct_by_stage"], "findings.swi_pct_by_stage"
        ))
    if "sleep_stages_pct" in f:
        errs.extend(_validate_stage_pct_dict(
            f["sleep_stages_pct"], "findings.sleep_stages_pct"
        ))
    if "spindle_age_norm_range" in f:
        r = f["spindle_age_norm_range"]
        if not (
            isinstance(r, list) and len(r) == 2
            and all(_schema._is_nonneg_finite(x) for x in r)
        ):
            errs.append(
                "findings.spindle_age_norm_range must be [low, high] floats"
            )

    return errs


def _validate_stage_pct_dict(d: Any, label: str) -> list[str]:
    errs: list[str] = []
    if not isinstance(d, dict):
        return [f"{label} must be an object"]
    extra = set(d.keys()) - _schema.SLEEP_STAGE_KEYS
    if extra:
        errs.append(f"{label} has unknown stage keys: {sorted(extra)}")
    for k, v in d.items():
        if k in _schema.SLEEP_STAGE_KEYS and not _schema._is_pct(v):
            errs.append(f"{label}.{k} must be a percentage [0, 100]")
    return errs


def _validate_intervention(it: Any) -> list[str]:
    errs: list[str] = []
    if not isinstance(it, dict):
        return ["intervention must be an object"]
    req = {"type", "name", "record_kind"}
    for k in req:
        if k not in it:
            errs.append(f"intervention missing required key {k!r}")
    if errs:
        return errs
    if it["type"] not in _schema.INTERVENTION_TYPES:
        errs.append("intervention.type not in allowed set")
    if not isinstance(it["name"], str) or not it["name"]:
        errs.append("intervention.name must be a non-empty string")
    elif len(it["name"]) > _schema.INTERVENTION_NAME_MAX_LEN:
        errs.append(
            f"intervention.name longer than "
            f"{_schema.INTERVENTION_NAME_MAX_LEN}"
        )
    if it["record_kind"] not in _schema.INTERVENTION_RECORD_KINDS:
        errs.append("intervention.record_kind not in allowed set")
    if it.get("linked_pre_submission_id") is not None:
        if not _schema._is_uuid(it["linked_pre_submission_id"]):
            errs.append("intervention.linked_pre_submission_id not a uuid4")
    extra = set(it.keys()) - (req | {"linked_pre_submission_id"})
    if extra:
        errs.append(f"intervention has unknown keys: {sorted(extra)}")
    return errs
