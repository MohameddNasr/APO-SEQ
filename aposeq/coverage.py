"""Coverage and depth workflow for APO-SEQ."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from aposeq.command import PlannedCommand, execute_command, format_command
from aposeq.config import LoadedConfig, resolve_run_samples
from aposeq.exceptions import ApoSeqError


@dataclass(frozen=True)
class CoveragePlan:
    """Planned coverage outputs and commands for one sample."""

    barcode: str
    genotype: str
    replicate: str | int | None
    sorted_bam: Path
    depth_tsv: Path
    summary_tsv: Path
    command: PlannedCommand


@dataclass(frozen=True)
class CoverageResult:
    """Result of planning or executing coverage analysis."""

    plans: tuple[CoveragePlan, ...]
    command_manifest_path: Path
    summary_manifest_path: Path
    dry_run: bool


def build_coverage_plans(
    loaded_config: LoadedConfig,
    samtools: str = "samtools",
) -> tuple[CoveragePlan, ...]:
    """Create `samtools depth` plans for all configured samples."""

    run = loaded_config.run
    defaults = loaded_config.defaults
    assay = loaded_config.assay
    output_directory = Path(str(run["output"]["directory"])).expanduser()
    depth_root = output_directory / "Coverage" / "depth"
    summary_root = output_directory / "Coverage" / "summaries"
    sorted_root = output_directory / "BAM" / "sorted"

    alignment_config = defaults.get("alignment", {})
    variant_config = defaults.get("variant_calling", {})
    mapq = str(alignment_config.get("mapq", 20))
    baseq = str(variant_config.get("mpileup_base_quality", 10))
    chromosome = str(assay["reference"]["chromosome"])

    plans: list[CoveragePlan] = []
    samples = resolve_run_samples(run, output_directory)
    if not samples:
        raise ApoSeqError("No samples are available for coverage planning")

    for sample in samples:
        barcode = str(sample["barcode"])
        sorted_bam = sorted_root / f"{barcode}.sorted.bam"
        depth_tsv = depth_root / f"{barcode}.depth.tsv"
        summary_tsv = summary_root / f"{barcode}.coverage_summary.tsv"
        command = PlannedCommand(
            step="depth",
            argv=(
                samtools,
                "depth",
                "-aa",
                "-q",
                mapq,
                "-Q",
                baseq,
                "-r",
                chromosome,
                str(sorted_bam),
            ),
            stdout_path=depth_tsv,
        )
        replicate = sample.get("replicate")
        plans.append(
            CoveragePlan(
                barcode=barcode,
                genotype=str(sample["genotype"]),
                replicate=replicate if isinstance(replicate, (str, int)) else None,
                sorted_bam=sorted_bam,
                depth_tsv=depth_tsv,
                summary_tsv=summary_tsv,
                command=command,
            )
        )
    return tuple(plans)


def run_coverage(
    loaded_config: LoadedConfig,
    dry_run: bool = False,
    samtools: str = "samtools",
) -> CoverageResult:
    """Plan and optionally run coverage analysis."""

    plans = build_coverage_plans(loaded_config, samtools=samtools)
    output_directory = Path(str(loaded_config.run["output"]["directory"])).expanduser()
    create_output_directories(output_directory)

    command_manifest_path = output_directory / "Coverage" / "coverage_commands.tsv"
    summary_manifest_path = output_directory / "Coverage" / "coverage_manifest.tsv"
    write_coverage_command_manifest(plans, command_manifest_path, dry_run=dry_run)
    write_coverage_manifest(plans, summary_manifest_path)

    if not dry_run:
        for plan in plans:
            execute_command(plan.command, context=plan.barcode)
            summarize_depth_file(plan.depth_tsv, plan.summary_tsv, plan)

    return CoverageResult(
        plans=plans,
        command_manifest_path=command_manifest_path,
        summary_manifest_path=summary_manifest_path,
        dry_run=dry_run,
    )


def create_output_directories(output_directory: Path) -> None:
    for path in (
        output_directory / "Coverage",
        output_directory / "Coverage" / "depth",
        output_directory / "Coverage" / "summaries",
    ):
        path.mkdir(parents=True, exist_ok=True)


def summarize_depth_file(depth_tsv: Path, summary_tsv: Path, plan: CoveragePlan) -> None:
    """Summarize a samtools depth file."""

    total_positions = 0
    covered_positions = 0
    total_depth = 0
    minimum_depth: int | None = None
    maximum_depth = 0

    with depth_tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            depth = int(fields[2])
            total_positions += 1
            total_depth += depth
            maximum_depth = max(maximum_depth, depth)
            minimum_depth = depth if minimum_depth is None else min(minimum_depth, depth)
            if depth > 0:
                covered_positions += 1

    mean_depth = total_depth / total_positions if total_positions else 0
    fraction_covered = covered_positions / total_positions if total_positions else 0
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "barcode",
                "genotype",
                "replicate",
                "total_positions",
                "covered_positions",
                "fraction_covered",
                "mean_depth",
                "minimum_depth",
                "maximum_depth",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "barcode": plan.barcode,
                "genotype": plan.genotype,
                "replicate": "" if plan.replicate is None else plan.replicate,
                "total_positions": total_positions,
                "covered_positions": covered_positions,
                "fraction_covered": f"{fraction_covered:.6f}",
                "mean_depth": f"{mean_depth:.6f}",
                "minimum_depth": 0 if minimum_depth is None else minimum_depth,
                "maximum_depth": maximum_depth,
            }
        )


def write_coverage_command_manifest(
    plans: Iterable[CoveragePlan],
    path: Path,
    dry_run: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["barcode", "step", "dry_run", "command"],
            delimiter="\t",
        )
        writer.writeheader()
        for plan in plans:
            writer.writerow(
                {
                    "barcode": plan.barcode,
                    "step": plan.command.step,
                    "dry_run": str(dry_run).lower(),
                    "command": format_command(plan.command),
                }
            )


def write_coverage_manifest(plans: Iterable[CoveragePlan], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "barcode",
                "genotype",
                "replicate",
                "sorted_bam",
                "depth_tsv",
                "summary_tsv",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for plan in plans:
            writer.writerow(
                {
                    "barcode": plan.barcode,
                    "genotype": plan.genotype,
                    "replicate": "" if plan.replicate is None else plan.replicate,
                    "sorted_bam": plan.sorted_bam,
                    "depth_tsv": plan.depth_tsv,
                    "summary_tsv": plan.summary_tsv,
                }
            )
