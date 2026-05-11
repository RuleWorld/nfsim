#!/usr/bin/env python3
"""Adapt NFsim CLI to RuleMonkey's rm_driver interface.

This lets RuleMonkey's Python harnesses validate an NFsim binary directly
against the vendored NFsim reference ensembles.

Expected rm_driver CLI:
    rulemonkey_nfsim_driver.py <model.xml> <t_end> <n_steps> <seed> [rm_flags...]

Useful environment variables:
    NFSIM_BIN          Path to NFsim executable (default: <repo>/build/NFsim)
    NFSIM_SIM_PARAMS   Path to RuleMonkey sim_params.tsv for model-specific flags
    NFSIM_EXTRA_FLAGS  Extra NFsim CLI flags, e.g. "-connect"
"""

from __future__ import annotations

import csv
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NFSIM_BIN = REPO_ROOT / "build" / "NFsim"


def _usage() -> str:
    return (
        "Usage: rulemonkey_nfsim_driver.py <model.xml> <t_end> <n_steps> <seed> "
        "[rm_flags...]"
    )


def _strip_header_hash(fieldnames: list[str] | None) -> list[str] | None:
    if fieldnames and fieldnames[0].startswith("#"):
        fieldnames = list(fieldnames)
        fieldnames[0] = fieldnames[0].lstrip("#")
    return fieldnames


def _load_model_flags(sim_params_path: Path | None, model_name: str) -> list[str]:
    if sim_params_path is None or not sim_params_path.exists():
        return []

    with sim_params_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        reader.fieldnames = _strip_header_hash(reader.fieldnames)
        for row in reader:
            model = (row.get("model") or "").strip()
            if not model or model.startswith("#"):
                continue
            if model != model_name:
                continue
            raw_flags = (row.get("nfsim_flags") or "").strip()
            return _normalize_nfsim_flags(shlex.split(raw_flags))
    return []


def _normalize_nfsim_flags(tokens: list[str]) -> list[str]:
    """Drop flags that the harness already supplies explicitly.

    Keep model-specific behavior flags such as -cb, -bscb, -gml, and -utl.
    """

    drop_with_value = {
        "-xml",
        "-o",
        "-sim",
        "-eq",
        "-oSteps",
        "-oTimes",
        "-seed",
        "-ss",
        "-rxnlog",
        "-logbuffer",
        "-maxcputime",
    }

    keep_with_value = {
        "-gml",
        "-utl",
    }

    out: list[str] = []
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in keep_with_value:
            if i + 1 < len(tokens):
                out.extend([tok, tokens[i + 1]])
            i += 2
            continue
        if tok in drop_with_value:
            i += 2
            continue
        out.append(tok)
        i += 1
    return out


def _dedupe_flags(tokens: list[str]) -> list[str]:
    """Keep first occurrence of standalone flags; last occurrence of value flags."""

    value_flags = {"-gml", "-utl"}
    standalone_seen: set[str] = set()
    value_map: dict[str, str] = {}
    ordered_values: list[str] = []
    out: list[str] = []

    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in value_flags:
            if i + 1 < len(tokens):
                if tok not in value_map:
                    ordered_values.append(tok)
                value_map[tok] = tokens[i + 1]
            i += 2
            continue
        if tok not in standalone_seen:
            standalone_seen.add(tok)
            out.append(tok)
        i += 1

    for tok in ordered_values:
        out.extend([tok, value_map[tok]])
    return out


def _normalize_gdat_text(raw_text: str) -> str:
    """Rewrite NFsim gdat text as clean tab-separated output."""

    out_lines: list[str] = []
    for raw_line in raw_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            parts = line.lstrip("#").strip().split()
            out_lines.append("#" + "\t".join(parts))
            continue
        parts = line.split()
        out_lines.append("\t".join(parts))
    return "\n".join(out_lines) + ("\n" if out_lines else "")


def main() -> int:
    if len(sys.argv) < 5:
        print(_usage(), file=sys.stderr)
        return 2

    xml_path = Path(sys.argv[1]).resolve()
    t_end = sys.argv[2]
    n_steps = sys.argv[3]
    seed = sys.argv[4]
    passthrough_flags = sys.argv[5:]

    nfsim_bin = Path(os.environ.get("NFSIM_BIN", str(DEFAULT_NFSIM_BIN))).resolve()
    sim_params_env = os.environ.get("NFSIM_SIM_PARAMS")
    sim_params_path = Path(sim_params_env).resolve() if sim_params_env else None
    extra_flags = shlex.split(os.environ.get("NFSIM_EXTRA_FLAGS", ""))

    if not nfsim_bin.exists():
        print(f"NFsim binary not found: {nfsim_bin}", file=sys.stderr)
        return 2
    if not xml_path.exists():
        print(f"XML not found: {xml_path}", file=sys.stderr)
        return 2

    model_name = xml_path.stem
    model_flags = _load_model_flags(sim_params_path, model_name)
    merged_flags = _dedupe_flags(model_flags + passthrough_flags + extra_flags)

    with tempfile.TemporaryDirectory(prefix=f"rm_nfsim_{model_name}_") as td_raw:
        out_gdat = Path(td_raw) / f"{model_name}.gdat"
        cmd = [
            str(nfsim_bin),
            "-xml",
            str(xml_path),
            "-sim",
            str(t_end),
            "-oSteps",
            str(n_steps),
            "-seed",
            str(seed),
            *merged_flags,
            "-o",
            str(out_gdat),
        ]
        result = subprocess.run(
            cmd,
            cwd=str(xml_path.parent),
            capture_output=True,
            text=True,
        )

        if result.returncode != 0 or not out_gdat.exists():
            print("NFsim driver wrapper failed.", file=sys.stderr)
            print(f"Command: {' '.join(shlex.quote(x) for x in cmd)}", file=sys.stderr)
            if result.stdout.strip():
                print("--- stdout ---", file=sys.stderr)
                print(result.stdout.strip(), file=sys.stderr)
            if result.stderr.strip():
                print("--- stderr ---", file=sys.stderr)
                print(result.stderr.strip(), file=sys.stderr)
            return 1

        sys.stdout.write(_normalize_gdat_text(out_gdat.read_text()))
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
