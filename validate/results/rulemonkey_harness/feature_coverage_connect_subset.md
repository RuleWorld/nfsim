# Feature Coverage Benchmark Report

Generated: 2026-05-11 15:28:57

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
- RM reps: 5, wall time: 0.293s
- vs NFsim: max_z=17.33 (AB), tz_max=17.22 — **FAIL**
- **Overall: FAIL**

### ft_delete_molecules
- Tier: base
- RM reps: 5, wall time: 0.288s
- vs NFsim: max_z=3.31 (B_free), tz_max=4.83 — **PASS**
- **Overall: PASS**

### ft_multi_product
- Tier: base
- RM reps: 5, wall time: 0.265s
- vs NFsim: max_z=25.97 (BC), tz_max=2329.70 — **FAIL**
- **Overall: FAIL**

