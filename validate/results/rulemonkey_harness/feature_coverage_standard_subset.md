# Feature Coverage Benchmark Report

Generated: 2026-05-11 15:29:02

**Summary: 1 PASS / 2 FAIL / 0 SKIP**

## Feature Coverage

| Model | Tier | Features | RM vs NFsim | RM vs ODE | Verdict |
|-------|------|----------|-------------|-----------|--------|
| combo_synth_degrade_equilibrium | combinations | Feature combination: synthesis + degradation + bin | FAIL | - | **FAIL** |
| ft_delete_molecules | base | Feature: DeleteMolecules keyword on degradation ru | PASS | - | PASS |
| ft_multi_product | base | Feature: rule producing multiple new molecules (on | FAIL | - | **FAIL** |

## Detailed Results

### combo_synth_degrade_equilibrium
- Tier: combinations
- RM reps: 5, wall time: 0.318s
- vs NFsim: max_z=9.79 (A_free), tz_max=5.94 — **FAIL**
- **Overall: FAIL**

### ft_delete_molecules
- Tier: base
- RM reps: 5, wall time: 0.284s
- vs NFsim: max_z=3.62 (A_total), tz_max=4.82 — **PASS**
- **Overall: PASS**

### ft_multi_product
- Tier: base
- RM reps: 5, wall time: 0.285s
- vs NFsim: max_z=7.42 (BC), tz_max=516.40 — **FAIL**
- **Overall: FAIL**

