# attic

Superseded code kept for the record. Nothing here is called by the pipeline.

- `background_model_v1/`: the original background builder (Pisces-VCF-based
  Beta fit and helpers). Replaced in 0.3.0 by `build_background.py`, which
  works from pileups.
- `superseded/`: the long-format error-model scorer, the artefact mask
  builder (mask retired in 0.3.0), the v1 indel blocklist builder, and the
  nohup batch launcher from the single-server days.

`patches/` (top level) holds the exact-anchor patch scripts that produced
the 0.3.0 engine, estimator and HGVS changes; they are the audit trail for
those edits and are already applied.
