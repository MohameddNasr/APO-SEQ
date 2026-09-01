"""Slurm submission helpers for APO-SEQ."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

from aposeq.config import LoadedConfig
from aposeq.exceptions import ApoSeqError


@dataclass(frozen=True)
class SlurmScript:
    """A generated APO-SEQ Slurm script."""

    path: Path
    command: tuple[str, ...]


def write_sbatch_script(
    loaded_config: LoadedConfig,
    script_path: Path | None = None,
    steps: tuple[str, ...] = ("process_bams", "coverage", "call_variants"),
    reference_fasta: Path | None = None,
    submit_directory: Path | None = None,
) -> SlurmScript:
    """Write an sbatch script for heavy APO-SEQ processing steps."""

    output_directory = Path(str(loaded_config.run["output"]["directory"])).expanduser()
    execution = loaded_config.execution or {}
    resources = execution.get("resources", {}) if isinstance(execution.get("resources"), dict) else {}
    logs_dir = output_directory / "logs"
    scripts_dir = output_directory / "slurm"
    logs_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    script_path = script_path or scripts_dir / f"{loaded_config.run['name']}.sbatch"
    threads = resources.get("threads", 16)
    memory = resources.get("memory", "64G")
    time = resources.get("time", "24:00:00")
    job_name = sanitize_job_name(f"aposeq_{loaded_config.run['name']}")
    modules = execution.get("modules", [])
    module_lines = []
    if isinstance(modules, list):
        module_lines = [f"module load {module}" for module in modules]

    submit_directory = submit_directory or Path.cwd()
    run_config = _relative_to_or_absolute(loaded_config.run_config_path, submit_directory)
    reference = _relative_to_or_absolute(reference_fasta, submit_directory) if reference_fasta else None

    command = [
        "python3",
        "-m",
        "aposeq.cli",
        "run",
        "--run-config",
        str(run_config),
        "--steps",
        *steps,
    ]
    if reference is not None:
        command.extend(["--reference-fasta", str(reference)])

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --cpus-per-task={threads}",
        f"#SBATCH --mem={memory}",
        f"#SBATCH --time={time}",
        f"#SBATCH --output={logs_dir}/%x-%j.out",
        f"#SBATCH --error={logs_dir}/%x-%j.err",
        "",
        "set -euo pipefail",
        "",
        *module_lines,
        "",
        'cd "${SLURM_SUBMIT_DIR:-' + str(submit_directory) + '}"',
        " ".join(command),
        "",
    ]
    script_path.write_text("\n".join(lines), encoding="utf-8")
    return SlurmScript(path=script_path, command=("sbatch", str(script_path)))


def submit_sbatch(script: SlurmScript) -> str:
    """Submit a generated Slurm script with sbatch."""

    try:
        completed = subprocess.run(
            script.command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise ApoSeqError("sbatch executable was not found") from exc
    except subprocess.CalledProcessError as exc:
        raise ApoSeqError(f"sbatch submission failed: {exc.stderr.strip()}") from exc
    return completed.stdout.strip()


def sanitize_job_name(name: str) -> str:
    """Return a Slurm-safe job name."""

    safe = "".join(char if char.isalnum() or char in {"_", "-"} else "_" for char in name)
    return safe[:128]


def _relative_to_or_absolute(path: Path | None, root: Path) -> Path | None:
    if path is None:
        return None
    try:
        return path.resolve().relative_to(root.resolve())
    except ValueError:
        return path
