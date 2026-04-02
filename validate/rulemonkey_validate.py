import argparse
import csv
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import bionetgen
except ImportError as exc:
    raise RuntimeError(
        "bionetgen is required for this validator. Install validate/requirements.txt first."
    ) from exc


NULL_EVENTS_RE = re.compile(r"Null events:\s*(\d+)")
REACTIONS_RE = re.compile(r"simulated\s+(\d+)\s+reactions\s+in\s+([0-9eE+\-.]+)s")


@dataclass
class RunMetrics:
    mode: str
    seed: int
    wall_time_sec: float
    null_events: Optional[int]
    reactions: Optional[int]
    reported_sim_time_sec: Optional[float]
    gdat_path: str
    log_path: str


def _workspace_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def _default_nfsim_path() -> str:
    root = _workspace_root()
    exe = "NFsim.exe" if os.name == "nt" else "NFsim"
    return os.path.join(root, "build", exe)


def _default_tlbr_path() -> str:
    root = _workspace_root()
    return os.path.join(root, "test", "tlbr", "tlbr.bngl")


def _tiny_tlbr_path() -> str:
    root = _workspace_root()
    return os.path.join(root, "validate", "tlbr_tiny.bngl")


def _aa_dimerization_path() -> str:
    root = _workspace_root()
    return os.path.join(root, "validate", "aa_dimerization.bngl")


def _load_gdat(path: str) -> Tuple[List[str], np.ndarray]:
    with open(path, "r", encoding="utf-8") as f:
        raw_header = re.sub(r"\s+", " ", f.readline().strip())
    header = [h for h in raw_header.split(" ") if h and h != "#"]
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if len(header) > data.shape[1]:
        header = header[-data.shape[1] :]
    return header, data


def _safe_mean(values: List[float]) -> float:
    return float(statistics.mean(values)) if values else float("nan")


def _safe_stdev(values: List[float]) -> float:
    if len(values) < 2:
        return 0.0
    return float(statistics.stdev(values))


def _extract_metrics(text: str) -> Tuple[Optional[int], Optional[int], Optional[float]]:
    null_events = None
    reactions = None
    sim_time = None

    null_match = NULL_EVENTS_RE.search(text)
    if null_match:
        null_events = int(null_match.group(1))

    rxn_match = REACTIONS_RE.search(text)
    if rxn_match:
        reactions = int(rxn_match.group(1))
        sim_time = float(rxn_match.group(2))

    return null_events, reactions, sim_time


def _run_bngl_to_xml(bngl_path: str, out_dir: str) -> str:
    bng2 = os.path.join(bionetgen.defaults.bng_path, "BNG2.pl")
    if not os.path.exists(bng2):
        raise FileNotFoundError(f"Could not find BNG2.pl at: {bng2}")

    model_name = os.path.basename(bngl_path)
    local_bngl = os.path.join(out_dir, model_name)
    shutil.copy2(bngl_path, local_bngl)

    cmd = ["perl", bng2, "-outdir", out_dir, "-log", local_bngl]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "BNG XML generation failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    xml_path = os.path.join(out_dir, os.path.splitext(model_name)[0] + ".xml")
    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"Expected XML file not found: {xml_path}")
    return xml_path


def _run_nfsim(
    nfsim_path: str,
    xml_path: str,
    out_dir: str,
    mode: str,
    seed: int,
    sim_time: float,
    osteps: int,
    extra_args: List[str],
    verbose: bool,
) -> RunMetrics:
    model_tag = os.path.splitext(os.path.basename(xml_path))[0]
    mode_suffix = "rm" if mode == "rulemonkey" else "std"
    gdat_path = os.path.join(out_dir, f"{model_tag}_{mode_suffix}_seed{seed}.gdat")
    log_path = os.path.join(out_dir, f"{model_tag}_{mode_suffix}_seed{seed}.log")

    cmd = [
        nfsim_path,
        "-xml",
        xml_path,
        "-o",
        gdat_path,
        "-sim",
        str(sim_time),
        "-oSteps",
        str(osteps),
        "-seed",
        str(seed),
    ]
    if verbose:
        cmd.append("-v")
    if mode == "rulemonkey":
        cmd.append("-rulemonkey")
    cmd.extend(extra_args)

    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - start

    full_log = (proc.stdout or "") + "\n" + (proc.stderr or "")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(full_log)

    if proc.returncode != 0:
        raise RuntimeError(
            "NFsim run failed.\n"
            f"Mode={mode} Seed={seed}\n"
            f"Command: {' '.join(cmd)}\n"
            f"Log: {log_path}"
        )

    if not os.path.exists(gdat_path):
        raise FileNotFoundError(f"Expected output not found: {gdat_path}")

    null_events, reactions, reported_sim = _extract_metrics(full_log)
    if verbose:
        print(
            f"[{mode}] seed={seed} wall={wall:.3f}s null={null_events} "
            f"reactions={reactions}"
        )

    return RunMetrics(
        mode=mode,
        seed=seed,
        wall_time_sec=wall,
        null_events=null_events,
        reactions=reactions,
        reported_sim_time_sec=reported_sim,
        gdat_path=gdat_path,
        log_path=log_path,
    )


def _compute_summary(
    std_runs: List[RunMetrics],
    rm_runs: List[RunMetrics],
    std_data: Dict[int, np.ndarray],
    rm_data: Dict[int, np.ndarray],
    headers: List[str],
) -> Dict[str, object]:
    obs_names = headers[1:]
    seeds = sorted(std_data.keys())

    std_stack = np.stack([std_data[s] for s in seeds], axis=0)
    rm_stack = np.stack([rm_data[s] for s in seeds], axis=0)

    std_mean = np.mean(std_stack, axis=0)
    rm_mean = np.mean(rm_stack, axis=0)
    std_std = np.std(std_stack, axis=0, ddof=1 if len(seeds) > 1 else 0)
    rm_std = np.std(rm_stack, axis=0, ddof=1 if len(seeds) > 1 else 0)

    obs_summary = []
    for idx, name in enumerate(obs_names, start=1):
        mean_abs_diff = float(np.mean(np.abs(std_mean[:, idx] - rm_mean[:, idx])))
        final_std = float(std_mean[-1, idx])
        final_rm = float(rm_mean[-1, idx])
        final_delta = final_rm - final_std

        pooled = std_std[:, idx] + rm_std[:, idx]
        pooled = np.where(pooled < 1e-12, 1e-12, pooled)
        max_z = float(np.max(np.abs((rm_mean[:, idx] - std_mean[:, idx]) / pooled)))

        obs_summary.append(
            {
                "observable": name,
                "mean_abs_diff": mean_abs_diff,
                "final_mean_standard": final_std,
                "final_mean_rulemonkey": final_rm,
                "final_delta": final_delta,
                "max_normalized_mean_diff": max_z,
            }
        )

    std_wall = [r.wall_time_sec for r in std_runs]
    rm_wall = [r.wall_time_sec for r in rm_runs]

    std_null = [r.null_events for r in std_runs if r.null_events is not None]
    rm_null = [r.null_events for r in rm_runs if r.null_events is not None]

    summary = {
        "n_replicates": len(seeds),
        "seeds": seeds,
        "runtime": {
            "standard_mean_wall_sec": _safe_mean(std_wall),
            "rulemonkey_mean_wall_sec": _safe_mean(rm_wall),
            "speedup_rm_vs_standard": _safe_mean(std_wall) / max(_safe_mean(rm_wall), 1e-12),
        },
        "null_events": {
            "standard_mean": _safe_mean(std_null),
            "rulemonkey_mean": _safe_mean(rm_null),
            "standard_stdev": _safe_stdev(std_null),
            "rulemonkey_stdev": _safe_stdev(rm_null),
        },
        "observables": obs_summary,
    }
    return summary


def _write_runs_csv(path: str, runs: List[RunMetrics], headers: List[str], data_by_seed: Dict[int, np.ndarray]) -> None:
    obs_names = headers[1:]
    fieldnames = [
        "mode",
        "seed",
        "wall_time_sec",
        "null_events",
        "reactions",
        "reported_sim_time_sec",
    ] + [f"final_{o}" for o in obs_names]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for run in runs:
            row = {
                "mode": run.mode,
                "seed": run.seed,
                "wall_time_sec": run.wall_time_sec,
                "null_events": run.null_events,
                "reactions": run.reactions,
                "reported_sim_time_sec": run.reported_sim_time_sec,
            }
            arr = data_by_seed[run.seed]
            for idx, obs in enumerate(obs_names, start=1):
                row[f"final_{obs}"] = float(arr[-1, idx])
            writer.writerow(row)


def _write_observable_summary_csv(path: str, summary: Dict[str, object]) -> None:
    rows = summary["observables"]
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _maybe_plot(
    out_dir: str,
    headers: List[str],
    seeds: List[int],
    std_data: Dict[int, np.ndarray],
    rm_data: Dict[int, np.ndarray],
) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    std_stack = np.stack([std_data[s] for s in seeds], axis=0)
    rm_stack = np.stack([rm_data[s] for s in seeds], axis=0)

    std_mean = np.mean(std_stack, axis=0)
    rm_mean = np.mean(rm_stack, axis=0)
    std_std = np.std(std_stack, axis=0, ddof=1 if len(seeds) > 1 else 0)
    rm_std = np.std(rm_stack, axis=0, ddof=1 if len(seeds) > 1 else 0)

    t = std_mean[:, 0]
    for idx, obs in enumerate(headers[1:], start=1):
        plt.figure(figsize=(8, 4.5))
        plt.plot(t, std_mean[:, idx], label="NFsim standard", linewidth=2)
        plt.fill_between(t, std_mean[:, idx] - std_std[:, idx], std_mean[:, idx] + std_std[:, idx], alpha=0.2)

        plt.plot(t, rm_mean[:, idx], label="RuleMonkey", linewidth=2)
        plt.fill_between(t, rm_mean[:, idx] - rm_std[:, idx], rm_mean[:, idx] + rm_std[:, idx], alpha=0.2)

        plt.xlabel("time")
        plt.ylabel(obs)
        plt.title(f"RuleMonkey validation: {obs}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"observable_{obs}.png"), dpi=150)
        plt.close()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate NFsim standard mode vs RuleMonkey mode with paired seeds."
    )
    parser.add_argument("--bngl", default=_default_tlbr_path(), help="Path to BNGL model.")
    parser.add_argument("--tiny", action="store_true", help="Use the reduced-count TLBR benchmark preset.")
    parser.add_argument(
        "--aa-dimerization",
        action="store_true",
        help="Use the A+A dimerization benchmark preset (catches identical-reactant molecularity issues).",
    )
    parser.add_argument("--nfsim", default=_default_nfsim_path(), help="Path to NFsim executable.")
    parser.add_argument("--replicates", type=int, default=30, help="Number of paired seeds to run.")
    parser.add_argument("--seed-start", type=int, default=1, help="Starting seed value.")
    parser.add_argument("--sim", type=float, default=3000.0, help="Simulation end time (-sim).")
    parser.add_argument("--osteps", type=int, default=300, help="Output steps (-oSteps).")
    parser.add_argument(
        "--extra-args",
        default="",
        help='Extra NFsim args passed to both modes, e.g. "-cb -utl 0".',
    )
    parser.add_argument(
        "--outdir",
        default=os.path.join(_workspace_root(), "validate", "results", f"rulemonkey_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
        help="Directory for outputs (logs, CSV, plots, summary).",
    )
    parser.add_argument("--skip-plots", action="store_true", help="Skip PNG plot generation.")
    parser.add_argument("--verbose", action="store_true", help="Print per-run progress.")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()

    if args.tiny:
        args.bngl = _tiny_tlbr_path()
    if args.aa_dimerization:
        args.bngl = _aa_dimerization_path()

    if not os.path.exists(args.bngl):
        raise FileNotFoundError(f"BNGL model not found: {args.bngl}")
    if not os.path.exists(args.nfsim):
        raise FileNotFoundError(f"NFsim executable not found: {args.nfsim}")

    os.makedirs(args.outdir, exist_ok=True)
    run_dir = os.path.join(args.outdir, "generated")
    os.makedirs(run_dir, exist_ok=True)

    print(f"[setup] outdir={args.outdir}")
    print(f"[setup] bngl={args.bngl}")
    print(f"[setup] nfsim={args.nfsim}")

    xml_path = _run_bngl_to_xml(args.bngl, run_dir)
    print(f"[setup] generated xml={xml_path}")

    extra = [x for x in args.extra_args.split(" ") if x.strip()]
    seeds = [args.seed_start + i for i in range(args.replicates)]

    std_runs: List[RunMetrics] = []
    rm_runs: List[RunMetrics] = []
    std_data: Dict[int, np.ndarray] = {}
    rm_data: Dict[int, np.ndarray] = {}

    headers_ref: Optional[List[str]] = None
    time_ref: Optional[np.ndarray] = None

    for seed in seeds:
        std = _run_nfsim(
            nfsim_path=args.nfsim,
            xml_path=xml_path,
            out_dir=args.outdir,
            mode="standard",
            seed=seed,
            sim_time=args.sim,
            osteps=args.osteps,
            extra_args=extra,
            verbose=args.verbose,
        )
        rm = _run_nfsim(
            nfsim_path=args.nfsim,
            xml_path=xml_path,
            out_dir=args.outdir,
            mode="rulemonkey",
            seed=seed,
            sim_time=args.sim,
            osteps=args.osteps,
            extra_args=extra,
            verbose=args.verbose,
        )

        h_std, d_std = _load_gdat(std.gdat_path)
        h_rm, d_rm = _load_gdat(rm.gdat_path)

        if h_std != h_rm:
            raise RuntimeError(f"Header mismatch at seed {seed}: standard={h_std}, rulemonkey={h_rm}")
        if d_std.shape != d_rm.shape:
            raise RuntimeError(f"Shape mismatch at seed {seed}: standard={d_std.shape}, rulemonkey={d_rm.shape}")

        if headers_ref is None:
            headers_ref = h_std
            time_ref = d_std[:, 0]
        else:
            if h_std != headers_ref:
                raise RuntimeError(f"Header drift detected at seed {seed}")
            if not np.allclose(d_std[:, 0], time_ref):
                raise RuntimeError(f"Standard time grid mismatch at seed {seed}")
            if not np.allclose(d_rm[:, 0], time_ref):
                raise RuntimeError(f"RuleMonkey time grid mismatch at seed {seed}")

        std_runs.append(std)
        rm_runs.append(rm)
        std_data[seed] = d_std
        rm_data[seed] = d_rm

    assert headers_ref is not None

    summary = _compute_summary(std_runs, rm_runs, std_data, rm_data, headers_ref)

    _write_runs_csv(os.path.join(args.outdir, "runs_standard.csv"), std_runs, headers_ref, std_data)
    _write_runs_csv(os.path.join(args.outdir, "runs_rulemonkey.csv"), rm_runs, headers_ref, rm_data)
    _write_observable_summary_csv(os.path.join(args.outdir, "observable_summary.csv"), summary)

    with open(os.path.join(args.outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    if not args.skip_plots:
        _maybe_plot(args.outdir, headers_ref, seeds, std_data, rm_data)

    print("\nValidation complete")
    print(f"  Replicates: {summary['n_replicates']}")
    print(f"  Mean wall time (standard): {summary['runtime']['standard_mean_wall_sec']:.4f}s")
    print(f"  Mean wall time (rulemonkey): {summary['runtime']['rulemonkey_mean_wall_sec']:.4f}s")
    print(f"  Speedup (standard/rulemonkey): {summary['runtime']['speedup_rm_vs_standard']:.4f}x")
    print(f"  Mean null events (standard): {summary['null_events']['standard_mean']}")
    print(f"  Mean null events (rulemonkey): {summary['null_events']['rulemonkey_mean']}")
    print(f"  Output directory: {args.outdir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
