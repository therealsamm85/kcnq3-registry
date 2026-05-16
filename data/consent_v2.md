# Registry consent — version 2

By contributing a submission to this registry, you affirm that:

1. **You have the right to contribute.** You are the patient, the
   patient's legal guardian, or have explicit consent from the
   patient/guardian to contribute these data.

2. **You understand what is uploaded.** The submission contains:
   - The variant (gene, protein change, type)
   - Age bucket (e.g., "5-7"), not exact age
   - Sex (F / M / X / unknown)
   - Optional ISO country code
   - Bucketed recording duration
   - Quantitative EEG findings (numbers, percentages, sleep stage
     fractions, etc.) — no raw EEG signal
   - Optional intervention metadata (medication name, "pre"/"post")

   **Added in v2 (schema v2, shipped with kcnq3-lens v0.13.2+):**
   - HFO (high-frequency oscillation) rate bucket — a research metric
     derived from recordings sampled at ≥500 Hz
   - SO-spindle coupling strength and angle buckets — descriptive
     measure of memory-consolidation circuitry maturation
   - Slow-wave density bucket — NREM3 slow-oscillation marker
   - IED detection method flag — indicates whether automated IED
     count used the ensemble heuristic or the SpikeNet path

   All v2 fields are bucketed (coarse bins) and contain no raw signal,
   no waveforms, no timestamps, and no free-text.

3. **You understand what is NOT uploaded.**
   - No name, date of birth, exact age, address, contact information
   - No filename, file path, recording date
   - No raw EEG signal, no waveforms, no images
   - No free-text narrative

4. **Your data may be re-shared.** Aggregated, k-anonymized cohort
   statistics may be downloaded and used by other families,
   clinicians, and researchers. Individual submission records remain
   in the public `data/registry.jsonl` file, but contain no
   identifying information by construction.

5. **You may withdraw your submission.** You keep your
   `submission_id` locally. To withdraw: open a GitHub issue
   referencing the id. The maintainer will remove the line in the
   next merge.

6. **This is not medical advice.** The registry is a research
   aggregation. Findings here do not constitute diagnosis or
   treatment recommendations. Do not change medications based on
   peer-comparison percentiles without consulting your clinician.

7. **No warranty.** The maintainer cannot guarantee against errors,
   downtime, or data loss. Treat the registry as research infrastructure,
   not as a clinical record.

If you understand and agree, check the consent box in the app. The app
will record your affirmation as `{"version": 2, "given": true,
"given_at_month": "YYYY-MM"}`. The exact day is not recorded.

If you change your mind later, you may withdraw any submission at any
time without explanation (see point 5).

---

*Consent v1 was valid through kcnq3-lens v0.13.1. Families who
contributed under v1 do not need to re-affirm for existing records;
new contributions from v0.13.2 onwards require v2 affirmation.*
