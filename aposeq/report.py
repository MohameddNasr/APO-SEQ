"""Run manifest and text report generation."""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

from aposeq import __version__
from aposeq.config import LoadedConfig
from aposeq.state import utc_now


def write_run_manifest(
    loaded_config: LoadedConfig,
    output_directory: Path,
    dry_run: bool,
    steps: list[str],
) -> Path:
    """Write a machine-readable run manifest."""

    manifest = {
        "pipeline": "APO-SEQ",
        "version": __version__,
        "run": loaded_config.run["name"],
        "assay": loaded_config.assay["name"],
        "dry_run": dry_run,
        "steps": steps,
        "created_at": utc_now(),
        "run_config_path": str(loaded_config.run_config_path),
        "assay_config_path": str(loaded_config.assay_config_path),
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
    path = output_directory / "run_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


def write_text_report(
    loaded_config: LoadedConfig,
    output_directory: Path,
    state: dict[str, Any],
    dry_run: bool,
) -> Path:
    """Write a compact human-readable run report."""

    path = output_directory / "Report.txt"
    lines = [
        "APO-SEQ Run Report",
        "==================",
        "",
        f"Run: {loaded_config.run['name']}",
        f"Assay: {loaded_config.assay['name']}",
        f"Dry run: {str(dry_run).lower()}",
        f"Generated: {utc_now()}",
        "",
        "Steps:",
    ]
    for name, step in state.get("steps", {}).items():
        status = step.get("status", "unknown")
        lines.append(f"- {name}: {status}")
        if step.get("error"):
            lines.append(f"  error: {step['error']}")
    lines.extend(
        [
            "",
            "Key output roots:",
            f"- BAM: {output_directory / 'BAM'}",
            f"- Coverage: {output_directory / 'Coverage'}",
            f"- Variants: {output_directory / 'Variants'}",
            f"- Analysis: {output_directory / 'Analysis'}",
            f"- IGV: {output_directory / 'IGV'}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
