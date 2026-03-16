"""Small Python API for launching NFsim and reading GDAT output.

This module is intentionally lightweight so users can script NFsim runs
without writing subprocess boilerplate in every project.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence
import csv
import os
import subprocess


@dataclass
class NFsimResult:
    """Result metadata for a completed NFsim run."""

    command: List[str]
    returncode: int
    xml_path: Path
    output_path: Path
    stdout: str
    stderr: str


def _repo_root_from_here() -> Path:
    return Path(__file__).resolve().parents[2]


def default_nfsim_binary(repo_root: Optional[Path] = None) -> Path:
    """Return the default local build path to the NFsim executable."""

    root = Path(repo_root) if repo_root is not None else _repo_root_from_here()
    if os.name == "nt":
        return root / "build" / "NFsim.exe"
    return root / "build" / "NFsim"


def run_nfsim(
    xml_path: Path,
    output_path: Optional[Path] = None,
    *,
    nfsim_binary: Optional[Path] = None,
    extra_args: Optional[Sequence[str]] = None,
    working_dir: Optional[Path] = None,
    check: bool = True,
    capture_output: bool = True,
) -> NFsimResult:
    """Run NFsim on an XML model and return execution metadata.

    Args:
        xml_path: Path to the input model XML.
        output_path: Destination GDAT path. Defaults to <xml_stem>_nf.gdat.
        nfsim_binary: Path to NFsim executable. Defaults to local build output.
        extra_args: Additional command line arguments (e.g. ['-sim', '100']).
        working_dir: Optional process working directory.
        check: Raise CalledProcessError on non-zero exit code when True.
        capture_output: Capture stdout/stderr when True.
    """

    xml_path = Path(xml_path)
    if output_path is None:
        output_path = xml_path.with_name(f"{xml_path.stem}_nf.gdat")
    output_path = Path(output_path)

    binary = Path(nfsim_binary) if nfsim_binary is not None else default_nfsim_binary()
    cmd = [str(binary), "-xml", str(xml_path), "-o", str(output_path)]
    if extra_args:
        cmd.extend([str(a) for a in extra_args])

    completed = subprocess.run(
        cmd,
        cwd=str(working_dir) if working_dir is not None else None,
        check=check,
        text=True,
        capture_output=capture_output,
    )

    return NFsimResult(
        command=cmd,
        returncode=completed.returncode,
        xml_path=xml_path,
        output_path=output_path,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
    )


def read_gdat(path: Path) -> List[dict]:
    """Read a space-delimited GDAT file into a list of row dictionaries."""

    rows: List[dict] = []
    path = Path(path)
    with path.open("r", newline="") as handle:
        reader = csv.reader(handle, delimiter=" ")
        headers: List[str] = []
        for raw in reader:
            row = [col for col in raw if col != ""]
            if not row:
                continue
            if not headers:
                headers = row
                continue
            rows.append({headers[i]: row[i] for i in range(min(len(headers), len(row)))})
    return rows


__all__ = [
    "NFsimResult",
    "default_nfsim_binary",
    "run_nfsim",
    "read_gdat",
]
