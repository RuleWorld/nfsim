# RuleMonkey Validation Framework (TLBR)

This validator runs paired NFsim simulations in standard mode and RuleMonkey mode using the TLBR model, then writes reproducible metrics and plots.

## What it measures

- Observable trajectory agreement between modes across replicate seeds.
- Null-event counts parsed from NFsim console output.
- Wall-clock runtime comparison and speedup.

## Inputs

Default model:

- `test/tlbr/tlbr.bngl`

Reduced-count benchmark preset:

- `validate/tlbr_tiny.bngl` via `--tiny`

Default executable:

- `build/NFsim.exe` (Windows)
- `build/NFsim` (Linux/macOS)

## Run

From repository root:

```bash
python validate/rulemonkey_validate.py --replicates 30 --sim 3000 --osteps 300 --verbose
```

For a quick speed check, use the tiny preset:

```bash
python validate/rulemonkey_validate.py --tiny --replicates 5 --sim 10 --osteps 10 --skip-plots
```

Optional flags:

- `--bngl <path>`: use a different BNGL model.
- `--tiny`: use the reduced TLBR benchmark model.
- `--nfsim <path>`: explicit NFsim executable path.
- `--seed-start <int>`: first seed in paired sequence.
- `--extra-args "..."`: extra NFsim args for both modes (example: `-cb`).
- `--outdir <path>`: output directory.
- `--skip-plots`: skip PNG generation.

## Outputs

In the output directory (`validate/results/rulemonkey_<timestamp>/` by default):

- `runs_standard.csv`: per-seed standard mode metrics.
- `runs_rulemonkey.csv`: per-seed RuleMonkey mode metrics.
- `observable_summary.csv`: aggregate observable comparison stats.
- `summary.json`: machine-readable report with runtime/null-event summary.
- `observable_<name>.png`: mean ± std overlay plots (unless `--skip-plots`).
- `tlbr_*_seed*.log`: raw NFsim output logs used for parsing metrics.

## Notes

- The validator enforces matched headers, shapes, and time grids across runs.
- If a run fails, the error includes the exact mode/seed and log file path.
