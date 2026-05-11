# RuleMonkey Benchmark Report

**Date:** 2026-05-11
**Commit:** `6d7f240 docs(readme,quickstart): document the cancellation hook for embedders`
**Reps per model:** 10
**NFsim reference:** 100-rep ensemble

## Correctness

- **screen**: max |RM_mean - NFsim_mean| / (NFsim_std / sqrt(n_reps)) over all (time, obs) pairs. Fast early-warning metric. Flags values ≥ 5.0 as "suspicious" but does not determine verdict. Historical benchmarks used this as the pass/fail criterion; for models with many rare-event Size_N observables it is statistically unreliable.
- **tz_max**: max over observables of the z-score of the per-rep trapezoidal time integral, computed against precomputed NFsim stats in `ensemble/{model}.tint.tsv`. Collapses each 1001-point trajectory to one number per observable per rep, eliminating single-time-point coincidences.
- **T**: per-model verdict threshold = max(5.0, 1.2 × tz_p99), where tz_p99 is the 99th percentile of `tz_max` from self-splitting the NFsim replicates at n=10 (see `tests/reference/nfsim/noise_floor.tsv`; provenance and regen recipe in the same directory's `PROVENANCE.md`). Adaptive to each model's intrinsic rare-event noise floor.
- **std_ratio**: max(RM_std / NFsim_std) across observables with nontrivial variance. Diagnostic for variance consistency; not part of verdict.
- **verdict**: PASS if `tz_max < T`, FAIL otherwise. Degenerate-observable mismatches (both stds zero, values differ) fail unconditionally.

## Efficiency

- **nfsim_s**: NFsim mean wall time (100-rep, from reference campaign).
- **rm_s**: RM mean wall time (10-rep average).
- **ev/s**: SSA events per wall-second (throughput).
- **sel/fire/obs/upd**: Phase breakdown as % of RM engine time.

## Results

| model | nfsim_s | rm_s | events | ev/s | sel% | fire% | obs% | upd% | screen | tz_max | T | std_ratio | worst_obs | verdict |
|-------|--------:|-----:|-------:|-----:|-----:|------:|-----:|-----:|-------:|-------:|----:|----------:|-----------|---------|
| r01 | — | 1.9 | — | — | 0.0 | 0.0 | 0.0 | 0.0 | 10.22 | 12.66 | 5.00 | 2.00 | Xp_free | FAIL |
| r05 | — | 0.3 | — | — | 0.0 | 0.0 | 0.0 | 0.0 | 7.83 | 8.51 | 5.00 | 1.86 | Complex | FAIL |
| r20 | — | 2.1 | — | — | 0.0 | 0.0 | 0.0 | 0.0 | 10.22 | 12.66 | 5.00 | 2.00 | Xp_free | FAIL |
| r22 | — | 0.3 | — | — | 0.0 | 0.0 | 0.0 | 0.0 | 6.87 | 10.28 | 5.00 | 1.00 | rib_elong | FAIL |
| r32 | — | 0.2 | — | — | 0.0 | 0.0 | 0.0 | 0.0 | 6.22 | 9.56 | 5.00 | 2.00 | Mintra | FAIL |

## Summary

| Metric | Count |
|--------|------:|
| PASS | 0 |
| FAIL | 5 |
| TIMEOUT | 0 |
| SKIP | 0 |
| **Total** | **5** |
