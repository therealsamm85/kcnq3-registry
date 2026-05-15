"""Registry validation + aggregation toolchain.

Vendored standalone — no dependencies beyond the Python stdlib. Mirrors
the equivalent modules in the KCNQ3-Lens app (src/registry/) so the
two ends of the pipeline cannot drift apart silently. When schema
changes, both copies must move in lockstep — CI runs the round-trip
test that catches drift.
"""
