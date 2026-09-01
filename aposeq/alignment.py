"""BAM discovery and preprocessing for APO-SEQ."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from aposeq.command import PlannedCommand, execute_command, format_command
from aposeq.config import LoadedConfig
from aposeq.exceptions import ApoSeqError


@dataclass(frozen=True)
class BamRecord:
    """A discovered BAM file assigned to a barcode."""

    path: Path
    barcode: str


@dataclass(frozen=True)
class SampleBams:
    """All discovered BAMs for one configured sample."""

    barcode: str
    genotype: str
    replicate: str | int | None
    bams: tuple[Path, ...]


@dataclass(frozen=True)
class BamProcessingPlan:
    """Planned BAM-processing outputs and commands for one sample."""

    barcode: str
    genotype: str
    replicate: str | int | None
    input_bams: tuple[Path, ...]
    merged_bam: Path
    sorted_bam: Path
    flagstat_path: Path
    commands: tuple[PlannedCommand, ...]


@dataclass(frozen=True)
class BamProcessingResult:
    """Result of planning or executing BAM preprocessing."""

    plans: tuple[BamProcessingPlan, ...]
    manifest_path: Path
    command_manifest_path: Path
    dry_run: bool


def discover_bams(
    bam_directory: Path,
    barcode_regex: str = r"(barcode[0-9]+)",
    recursive: bool = True,
) -> tuple[BamRecord, ...]:
    """Find BAM files and assign each to a barcode using a regular expression."""

    if not bam_directory.exists():
        raise ApoSeqError(f"BAM directory does not exist: {bam_directory}")
    if not bam_directory.is_dir():
        raise ApoSeqError(f"BAM input path is not a directory: {bam_directory}")

    try:
        pattern = re.compile(barcode_regex)
    except re.error as exc:
        raise ApoSeqError(f"Invalid barcode regex {barcode_regex!r}: {exc}") from exc

    paths = bam_directory.rglob("*.bam") if recursive else bam_directory.glob("*.bam")
    records: list[BamRecord] = []
    for path in sorted(paths):
        match = pattern.search(str(path))
        if match is None:
            continue
        barcode = match.group(1) if match.groups() else match.group(0)
        records.append(BamRecord(path=path.resolve(), barcode=barcode))
    return tuple(records)


def group_bams_by_sample(
    records: Iterable[BamRecord],
    samples: list[dict[str, object]],
) -> tuple[SampleBams, ...]:
    """Group discovered BAM records according to run-config samples."""

    by_barcode: dict[str, list[Path]] = {}
    for record in records:
        by_barcode.setdefault(record.barcode, []).append(record.path)

    grouped: list[SampleBams] = []
    for sample in samples:
        barcode = str(sample["barcode"])
        genotype = str(sample["genotype"])
        replicate = sample.get("replicate")
        grouped.append(
            SampleBams(
                barcode=barcode,
                genotype=genotype,
                replicate=replicate if isinstance(replicate, (str, int)) else None,
                bams=tuple(sorted(by_barcode.get(barcode, []))),
            )
        )
    return tuple(grouped)


def build_bam_processing_plans(
    grouped: Iterable[SampleBams],
    output_directory: Path,
    samtools: str = "samtools",
) -> tuple[BamProcessingPlan, ...]:
    """Build merge, sort, and index commands for each sample."""

    bam_root = output_directory / "BAM"
    merged_root = bam_root / "merged"
    sorted_root = bam_root / "sorted"
    metrics_root = bam_root / "metrics"

    plans: list[BamProcessingPlan] = []
    for sample in grouped:
        if not sample.bams:
            raise ApoSeqError(f"No BAM files found for configured barcode: {sample.barcode}")
        merged_bam = merged_root / f"{sample.barcode}.merged.bam"
        sorted_bam = sorted_root / f"{sample.barcode}.sorted.bam"
        flagstat_path = metrics_root / f"{sample.barcode}.flagstat.txt"
        commands = (
            PlannedCommand(
                step="merge",
                argv=(samtools, "merge", "-f", str(merged_bam), *map(str, sample.bams)),
            ),
            PlannedCommand(
                step="sort",
                argv=(samtools, "sort", "-o", str(sorted_bam), str(merged_bam)),
            ),
            PlannedCommand(
                step="index",
                argv=(samtools, "index", str(sorted_bam)),
            ),
            PlannedCommand(
                step="flagstat",
                argv=(samtools, "flagstat", str(sorted_bam)),
                stdout_path=flagstat_path,
            ),
        )
        plans.append(
            BamProcessingPlan(
                barcode=sample.barcode,
                genotype=sample.genotype,
                replicate=sample.replicate,
                input_bams=sample.bams,
                merged_bam=merged_bam,
                sorted_bam=sorted_bam,
                flagstat_path=flagstat_path,
                commands=commands,
            )
        )
    return tuple(plans)


def process_bams(
    loaded_config: LoadedConfig,
    dry_run: bool = False,
    samtools: str = "samtools",
) -> BamProcessingResult:
    """Discover, group, and optionally process BAMs for a run configuration."""

    run = loaded_config.run
    defaults = loaded_config.defaults
    input_config = run["input"]
    output_config = run["output"]
    alignment_config = defaults.get("alignment", {})

    bam_directory = Path(str(input_config["bam_directory"])).expanduser()
    barcode_regex = str(input_config.get("barcode_regex", r"(barcode[0-9]+)"))
    recursive = bool(alignment_config.get("recursive_bam_discovery", True))
    output_directory = Path(str(output_config["directory"])).expanduser()

    records = discover_bams(
        bam_directory=bam_directory,
        barcode_regex=barcode_regex,
        recursive=recursive,
    )
    samples = run.get("samples", [])
    if input_config.get("auto_discover_samples"):
        samples = infer_samples_from_bams(
            records,
            genotype=str(input_config["default_genotype"]),
        )
        run["samples"] = samples
    grouped = group_bams_by_sample(records, samples)
    plans = build_bam_processing_plans(grouped, output_directory, samtools=samtools)

    create_output_directories(output_directory)
    manifest_path = output_directory / "BAM" / "bam_manifest.tsv"
    command_manifest_path = output_directory / "BAM" / "bam_commands.tsv"
    write_bam_manifest(plans, manifest_path)
    write_command_manifest(plans, command_manifest_path, dry_run=dry_run)

    if not dry_run:
        for plan in plans:
            execute_plan(plan)

    return BamProcessingResult(
        plans=plans,
        manifest_path=manifest_path,
        command_manifest_path=command_manifest_path,
        dry_run=dry_run,
    )


def infer_samples_from_bams(records: Iterable[BamRecord], genotype: str) -> list[dict[str, object]]:
    """Infer sample metadata for every discovered barcode."""

    barcodes = sorted({record.barcode for record in records})
    return [
        {
            "barcode": barcode,
            "genotype": genotype,
            "replicate": index,
        }
        for index, barcode in enumerate(barcodes, start=1)
    ]


def create_output_directories(output_directory: Path) -> None:
    """Create BAM output directories."""

    for path in (
        output_directory / "BAM",
        output_directory / "BAM" / "merged",
        output_directory / "BAM" / "sorted",
        output_directory / "BAM" / "metrics",
        output_directory / "logs",
    ):
        path.mkdir(parents=True, exist_ok=True)


def execute_plan(plan: BamProcessingPlan) -> None:
    """Run all shell commands for one BAM-processing plan."""

    for command in plan.commands:
        execute_command(command, context=plan.barcode)


def write_bam_manifest(plans: Iterable[BamProcessingPlan], path: Path) -> None:
    """Write one row per input BAM with its planned outputs."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "barcode",
                "genotype",
                "replicate",
                "input_bam",
                "merged_bam",
                "sorted_bam",
                "flagstat",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for plan in plans:
            for bam in plan.input_bams:
                writer.writerow(
                    {
                        "barcode": plan.barcode,
                        "genotype": plan.genotype,
                        "replicate": "" if plan.replicate is None else plan.replicate,
                        "input_bam": bam,
                        "merged_bam": plan.merged_bam,
                        "sorted_bam": plan.sorted_bam,
                        "flagstat": plan.flagstat_path,
                    }
                )


def write_command_manifest(
    plans: Iterable[BamProcessingPlan],
    path: Path,
    dry_run: bool,
) -> None:
    """Write planned command lines for reproducibility."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["barcode", "step", "dry_run", "command"],
            delimiter="\t",
        )
        writer.writeheader()
        for plan in plans:
            for command in plan.commands:
                command_text = format_command(command)
                writer.writerow(
                    {
                        "barcode": plan.barcode,
                        "step": command.step,
                        "dry_run": str(dry_run).lower(),
                        "command": command_text,
                    }
                )
