# Background error model

Site- and substitution-specific error model built from biological negative
controls, after Waalkes et al. Haematologica 2017;102(9):1549-1557.

See `error_model.md` for the full procedure. Summary:

1. Pisces on the 8 BNC consensus BAMs, `--minbq 30 --minvf 0.0001 -c 1`
2. `remove_variants_gtr_20.pl`  -- drop VAF > 0.2 (germline)
3. `print_multiple_variants_at_same_location.pl` -- split multi-allelic
4. `fill_empty_mips.pl`  -- write zero counts at uncalled positions
5. `beta_distribution.py` -- fit per-position Beta, emit the matrix

Run once per consensus track. `default_error_rate` in `beta_distribution.py`
must be set per track: 1/15000 for SSCS, 1/200000 for DCS, both from Wang et al.
NAR 2019 supplementary table 5 (SmallDeep panel).

## Known issues

- `default_error_rate` is a module-level constant, so the two tracks currently
  require two copies of the script. It should become a command-line argument.

- `fill_empty_mips.pl` iterates `$start <= $i <= $end` against a 0-based BED
  start, so each probe beginning a contiguous block emits one extra position at
  its 5' end. 97 of 182 probes are block starts, giving 21,937 rows where the
  panel spans 21,840 positions. The surplus rows carry zero counts and take the
  default, so the model is unaffected, but the matrix cannot be joined to other
  tables on position without dropping them.

- Depth at filled positions is carried forward from the last matched variant
  line rather than being the true depth at that position. Error rates are
  unaffected (0/anything = 0) but the depth column is not usable, which rules
  out any depth-aware test built from this file.

- Pisces applies `MinimumVariantQScore` (default 20) beneath the requested
  `--minvf`, censoring roughly everything below 7-12 alt reads. Positions below
  that are recorded as zeros and take the default. The default is therefore
  doing the work at most positions, and it is a literature constant rather than
  a measurement from this panel.
