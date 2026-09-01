"""Variant calling and VCF parsing for APO-SEQ."""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from aposeq.command import PlannedCommand, execute_command, format_command
from aposeq.config import LoadedConfig, resolve_run_samples
from aposeq.exceptions import ApoSeqError


@dataclass(frozen=True)
class VariantPlan:
    """Planned variant-calling outputs and commands for one sample."""

    barcode: str
    genotype: str
    replicate: str | int | None
    sorted_bam: Path
    raw_bcf: Path
    vcf_gz: Path
    mutation_tsv: Path
    commands: tuple[PlannedCommand, ...]


@dataclass(frozen=True)
class VariantCallingResult:
    """Result of planning or executing variant calling."""

    plans: tuple[VariantPlan, ...]
    command_manifest_path: Path
    master_mutation_table: Path
    dry_run: bool


def build_variant_plans(
    loaded_config: LoadedConfig,
    bcftools: str = "bcftools",
) -> tuple[VariantPlan, ...]:
    """Create bcftools mpileup/call plans for all configured samples."""

    run = loaded_config.run
    defaults = loaded_config.defaults
    assay = loaded_config.assay
    output_directory = Path(str(run["output"]["directory"])).expanduser()
    sorted_root = output_directory / "BAM" / "sorted"
    variant_root = output_directory / "Variants"
    raw_root = variant_root / "raw"
    vcf_root = variant_root / "vcf"
    table_root = variant_root / "tables"

    alignment_config = defaults.get("alignment", {})
    variant_config = defaults.get("variant_calling", {})
    mapq = str(alignment_config.get("mapq", 20))
    baseq = str(variant_config.get("mpileup_base_quality", 10))
    reference_fasta = str(assay["reference"]["fasta"])

    plans: list[VariantPlan] = []
    samples = resolve_run_samples(run, output_directory)
    if not samples:
        raise ApoSeqError("No samples are available for variant-calling planning")

    for sample in samples:
        barcode = str(sample["barcode"])
        sorted_bam = sorted_root / f"{barcode}.sorted.bam"
        raw_bcf = raw_root / f"{barcode}.raw.bcf"
        vcf_gz = vcf_root / f"{barcode}.vcf.gz"
        mutation_tsv = table_root / f"{barcode}.mutations.tsv"
        commands = (
            PlannedCommand(
                step="mpileup",
                argv=(
                    bcftools,
                    "mpileup",
                    "-f",
                    reference_fasta,
                    "-q",
                    mapq,
                    "-Q",
                    baseq,
                    "-a",
                    "FORMAT/AD,FORMAT/DP",
                    "-Ob",
                    "-o",
                    str(raw_bcf),
                    str(sorted_bam),
                ),
            ),
            PlannedCommand(
                step="call",
                argv=(bcftools, "call", "-mv", "-Oz", "-o", str(vcf_gz), str(raw_bcf)),
            ),
            PlannedCommand(
                step="index_vcf",
                argv=(bcftools, "index", "-t", str(vcf_gz)),
            ),
        )
        replicate = sample.get("replicate")
        plans.append(
            VariantPlan(
                barcode=barcode,
                genotype=str(sample["genotype"]),
                replicate=replicate if isinstance(replicate, (str, int)) else None,
                sorted_bam=sorted_bam,
                raw_bcf=raw_bcf,
                vcf_gz=vcf_gz,
                mutation_tsv=mutation_tsv,
                commands=commands,
            )
        )
    return tuple(plans)


def run_variant_calling(
    loaded_config: LoadedConfig,
    dry_run: bool = False,
    bcftools: str = "bcftools",
) -> VariantCallingResult:
    """Plan and optionally run variant calling."""

    plans = build_variant_plans(loaded_config, bcftools=bcftools)
    output_directory = Path(str(loaded_config.run["output"]["directory"])).expanduser()
    create_output_directories(output_directory)

    command_manifest_path = output_directory / "Variants" / "variant_commands.tsv"
    master_mutation_table = output_directory / "Variants" / "master_mutation_table.tsv"
    write_variant_command_manifest(plans, command_manifest_path, dry_run=dry_run)

    if not dry_run:
        mutation_rows: list[dict[str, object]] = []
        chromosome = str(loaded_config.assay["reference"]["chromosome"])
        for plan in plans:
            for command in plan.commands:
                execute_command(command, context=plan.barcode)
            rows = parse_vcf_to_rows(
                plan.vcf_gz,
                barcode=plan.barcode,
                genotype=plan.genotype,
                replicate=plan.replicate,
                expected_chromosome=chromosome,
            )
            write_mutation_table(rows, plan.mutation_tsv)
            mutation_rows.extend(rows)
        write_mutation_table(mutation_rows, master_mutation_table)

    return VariantCallingResult(
        plans=plans,
        command_manifest_path=command_manifest_path,
        master_mutation_table=master_mutation_table,
        dry_run=dry_run,
    )


def create_output_directories(output_directory: Path) -> None:
    for path in (
        output_directory / "Variants",
        output_directory / "Variants" / "raw",
        output_directory / "Variants" / "vcf",
        output_directory / "Variants" / "tables",
    ):
        path.mkdir(parents=True, exist_ok=True)


def parse_vcf_to_rows(
    vcf_path: Path,
    barcode: str,
    genotype: str,
    replicate: str | int | None = None,
    expected_chromosome: str | None = None,
) -> list[dict[str, object]]:
    """Parse SNV records from a VCF into APO-SEQ mutation-table rows."""

    rows: list[dict[str, object]] = []
    sample_name: str | None = None
    for line in _open_text(vcf_path):
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            sample_name = header[9] if len(header) > 9 else None
            continue
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue
        chrom, pos, _record_id, ref, alt_text, qual, filt, info = fields[:8]
        if expected_chromosome and chrom != expected_chromosome:
            continue
        alts = alt_text.split(",")
        format_keys = fields[8].split(":") if len(fields) > 8 else []
        sample_values = fields[9].split(":") if len(fields) > 9 else []
        sample_data = dict(zip(format_keys, sample_values))
        depth = _extract_depth(info, sample_data)
        alt_depths = _extract_alt_depths(sample_data, len(alts))

        for alt_index, alt in enumerate(alts):
            if len(ref) != 1 or len(alt) != 1:
                continue
            alt_count = alt_depths[alt_index] if alt_index < len(alt_depths) else 0
            allele_frequency = alt_count / depth if depth else 0
            rows.append(
                {
                    "BARCODE": barcode,
                    "GENOTYPE": genotype,
                    "REPLICATE": "" if replicate is None else replicate,
                    "CHROM": chrom,
                    "POS": int(pos),
                    "REF": ref,
                    "ALT": alt,
                    "MUTATION": f"{ref}>{alt}",
                    "QUAL": qual,
                    "FILTER": filt,
                    "DEPTH": depth,
                    "ALT_COUNT": alt_count,
                    "AF": f"{allele_frequency:.6f}",
                    "VCF_SAMPLE": "" if sample_name is None else sample_name,
                }
            )
    return rows


def write_mutation_table(rows: Iterable[dict[str, object]], path: Path) -> None:
    """Write APO-SEQ mutation rows as TSV."""

    fieldnames = [
        "BARCODE",
        "GENOTYPE",
        "REPLICATE",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "MUTATION",
        "QUAL",
        "FILTER",
        "DEPTH",
        "ALT_COUNT",
        "AF",
        "VCF_SAMPLE",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_variant_command_manifest(
    plans: Iterable[VariantPlan],
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
            for command in plan.commands:
                writer.writerow(
                    {
                        "barcode": plan.barcode,
                        "step": command.step,
                        "dry_run": str(dry_run).lower(),
                        "command": format_command(command),
                    }
                )


def _open_text(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def _extract_depth(info: str, sample_data: dict[str, str]) -> int:
    if sample_data.get("DP", "").isdigit():
        return int(sample_data["DP"])
    for entry in info.split(";"):
        if entry.startswith("DP="):
            value = entry.split("=", 1)[1]
            if value.isdigit():
                return int(value)
    ad = sample_data.get("AD")
    if ad:
        try:
            return sum(int(value) for value in ad.split(",") if value != ".")
        except ValueError as exc:
            raise ApoSeqError(f"Could not parse FORMAT/AD value: {ad}") from exc
    return 0


def _extract_alt_depths(sample_data: dict[str, str], alt_count: int) -> list[int]:
    ad = sample_data.get("AD")
    if not ad:
        return [0] * alt_count
    try:
        values = [int(value) for value in ad.split(",") if value != "."]
    except ValueError as exc:
        raise ApoSeqError(f"Could not parse FORMAT/AD value: {ad}") from exc
    return values[1:]
