"""Automated IGV snapshot selection and batch generation."""

from __future__ import annotations

import csv
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from aposeq.config import LoadedConfig
from aposeq.exceptions import ApoSeqError
from aposeq.filters import filter_mutations
from aposeq.table import load_mutation_table, save_table


@dataclass(frozen=True)
class IgvSnapshot:
    """One selected mutation and its planned IGV snapshot output."""

    assay: str
    run: str
    genotype: str
    barcode: str
    chromosome: str
    position: int
    ref: str
    alt: str
    mutation: str
    af: float
    alt_count: int
    depth: int
    motif: str
    bam: Path
    snapshot: Path
    locus: str
    selection_reason: str
    status: str = "planned"


@dataclass(frozen=True)
class IgvResult:
    """Outputs written by the IGV workflow."""

    selected_mutations: Path
    snapshot_manifest: Path
    failed_snapshots: Path
    batch_files: tuple[Path, ...]
    dry_run: bool


def run_igv_snapshots(
    loaded_config: LoadedConfig,
    mutation_table: Path | None = None,
    dry_run: bool = True,
    igv_executable: str = "igv.sh",
) -> IgvResult:
    """Select mutations, write IGV batch files, and optionally execute IGV."""

    output_directory = Path(str(loaded_config.run["output"]["directory"])).expanduser()
    mutation_table = mutation_table or _default_mutation_table(output_directory)
    table = load_mutation_table(mutation_table)

    igv_root = output_directory / "IGV"
    snapshot_root = igv_root / "snapshots"
    batch_root = igv_root / "batches"
    for path in (igv_root, snapshot_root, batch_root):
        path.mkdir(parents=True, exist_ok=True)

    selected = select_mutations_for_igv(table, loaded_config)
    selected_path = igv_root / "selected_mutations.tsv"
    save_table(selected, selected_path)

    snapshots = build_snapshot_manifest(selected, loaded_config, output_directory)
    manifest_path = igv_root / "snapshot_manifest.tsv"
    failed_path = igv_root / "failed_snapshots.tsv"
    write_snapshot_manifest(snapshots, manifest_path)
    write_failed_snapshots([], failed_path)

    batch_files = write_igv_batch_files(snapshots, loaded_config, batch_root, snapshot_root)
    if not dry_run:
        failures = execute_igv_batches(batch_files, igv_executable=igv_executable)
        write_failed_snapshots(failures, failed_path)

    return IgvResult(
        selected_mutations=selected_path,
        snapshot_manifest=manifest_path,
        failed_snapshots=failed_path,
        batch_files=batch_files,
        dry_run=dry_run,
    )


def select_mutations_for_igv(table: pd.DataFrame, loaded_config: LoadedConfig) -> pd.DataFrame:
    """Select mutation rows for IGV according to assay/run IGV settings."""

    igv_config = _effective_igv_config(loaded_config)
    selection = igv_config.get("selection", {})
    if not isinstance(selection, dict):
        selection = {}

    minimum_af = float(selection.get("minimum_af", table["AF"].min() if not table.empty else 0))
    minimum_alt_count = int(selection.get("minimum_alt_count", 0))
    apobec_only = bool(selection.get("apobec_only", False))
    maximum_snapshots = int(selection.get("maximum_snapshots", igv_config.get("maximum_snapshots", 500)))
    modes = selection.get("mode", ["passing_mutations"])
    if isinstance(modes, str):
        modes = [modes]

    filtered = filter_mutations(
        table,
        minimum_af=minimum_af,
        minimum_alt_count=minimum_alt_count,
        apobec_only=apobec_only,
        deduplicate_positions=False,
    )

    if filtered.empty:
        selected = filtered.copy()
        selected["SELECTION_REASON"] = []
        return selected

    reason_by_index: dict[int, set[str]] = {}
    if "passing_mutations" in modes:
        for index in filtered.index:
            reason_by_index.setdefault(index, set()).add("passing_mutations")

    if "recurrent" in modes:
        recurrent_positions = (
            filtered.groupby(["CHROM", "POS"], dropna=False)["BARCODE"].nunique()
            >= 2
        )
        recurrent_keys = set(recurrent_positions[recurrent_positions].index)
        for index, row in filtered.iterrows():
            if (row["CHROM"], row["POS"]) in recurrent_keys:
                reason_by_index.setdefault(index, set()).add("recurrent")

    if "outside_expected_region" in modes:
        inside_regions = _expected_region_ranges(loaded_config)
        for index, row in filtered.iterrows():
            if not _position_in_any_region(int(row["POS"]), inside_regions):
                reason_by_index.setdefault(index, set()).add("outside_expected_region")

    selected = filtered.loc[sorted(reason_by_index)].copy()
    selected["SELECTION_REASON"] = [
        ",".join(sorted(reason_by_index[index])) for index in selected.index
    ]
    selected = selected.sort_values(["SELECTION_REASON", "GENOTYPE", "BARCODE", "POS"])
    if maximum_snapshots > 0:
        selected = selected.head(maximum_snapshots)
    return selected.reset_index(drop=True)


def build_snapshot_manifest(
    selected: pd.DataFrame,
    loaded_config: LoadedConfig,
    output_directory: Path,
) -> tuple[IgvSnapshot, ...]:
    """Convert selected mutations into snapshot manifest records."""

    igv_config = _effective_igv_config(loaded_config)
    flank_size = int(igv_config.get("flank_size", 100))
    image_format = str(igv_config.get("image_format", "png"))
    assay_name = str(loaded_config.assay["name"])
    run_name = str(loaded_config.run["name"])
    chromosome = str(loaded_config.assay["reference"]["chromosome"])
    reference_length = int(loaded_config.assay["reference"]["length"])
    sorted_root = output_directory / "BAM" / "sorted"
    snapshot_root = output_directory / "IGV" / "snapshots"

    snapshots: list[IgvSnapshot] = []
    for row in selected.itertuples(index=False):
        barcode = str(row.BARCODE)
        position = int(row.POS)
        start = max(1, position - flank_size)
        end = min(reference_length, position + flank_size)
        ref = str(row.REF)
        alt = str(row.ALT)
        filename = sanitize_snapshot_name(
            f"{chromosome}_{position}_{ref}_{alt}.{image_format}"
        )
        snapshot = snapshot_root / barcode / filename
        motif = str(getattr(row, "MOTIF", ""))
        reason = str(getattr(row, "SELECTION_REASON", "passing_mutations"))
        snapshots.append(
            IgvSnapshot(
                assay=assay_name,
                run=run_name,
                genotype=str(row.GENOTYPE),
                barcode=barcode,
                chromosome=chromosome,
                position=position,
                ref=ref,
                alt=alt,
                mutation=str(row.MUTATION),
                af=float(row.AF),
                alt_count=int(row.ALT_COUNT),
                depth=int(row.DEPTH),
                motif="" if motif == "nan" else motif,
                bam=sorted_root / f"{barcode}.sorted.bam",
                snapshot=snapshot,
                locus=f"{chromosome}:{start}-{end}",
                selection_reason=reason,
            )
        )
    return tuple(snapshots)


def write_igv_batch_files(
    snapshots: tuple[IgvSnapshot, ...],
    loaded_config: LoadedConfig,
    batch_root: Path,
    snapshot_root: Path,
) -> tuple[Path, ...]:
    """Write one IGV batch file per barcode."""

    if not snapshots:
        return tuple()

    igv_config = _effective_igv_config(loaded_config)
    display = igv_config.get("display", {})
    if not isinstance(display, dict):
        display = {}
    reference_fasta = str(loaded_config.assay["reference"]["fasta"])

    by_barcode: dict[str, list[IgvSnapshot]] = {}
    for snapshot in snapshots:
        by_barcode.setdefault(snapshot.barcode, []).append(snapshot)

    batch_files: list[Path] = []
    for barcode, barcode_snapshots in sorted(by_barcode.items()):
        batch_path = batch_root / f"{barcode}.batch"
        first = barcode_snapshots[0]
        lines = [
            "new",
            f"genome {reference_fasta}",
            f"load {first.bam}",
            f"snapshotDirectory {snapshot_root / barcode}",
        ]
        sort_by = display.get("sort_alignments_by")
        color_by = display.get("color_alignments_by")
        if sort_by:
            lines.append(f"sort {sort_by}")
        if color_by:
            lines.append(f"colorBy {color_by}")
        for snapshot in barcode_snapshots:
            lines.extend(
                [
                    f"goto {snapshot.locus}",
                    "collapse",
                    f"snapshot {snapshot.snapshot.name}",
                ]
            )
        lines.append("exit")
        batch_path.parent.mkdir(parents=True, exist_ok=True)
        batch_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        batch_files.append(batch_path)
    return tuple(batch_files)


def write_snapshot_manifest(snapshots: tuple[IgvSnapshot, ...], path: Path) -> None:
    """Write the IGV snapshot manifest."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "assay",
                "run",
                "genotype",
                "barcode",
                "chromosome",
                "position",
                "ref",
                "alt",
                "mutation",
                "AF",
                "ALT_COUNT",
                "DEPTH",
                "motif",
                "bam",
                "snapshot",
                "locus",
                "selection_reason",
                "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for snapshot in snapshots:
            writer.writerow(
                {
                    "assay": snapshot.assay,
                    "run": snapshot.run,
                    "genotype": snapshot.genotype,
                    "barcode": snapshot.barcode,
                    "chromosome": snapshot.chromosome,
                    "position": snapshot.position,
                    "ref": snapshot.ref,
                    "alt": snapshot.alt,
                    "mutation": snapshot.mutation,
                    "AF": f"{snapshot.af:.6f}",
                    "ALT_COUNT": snapshot.alt_count,
                    "DEPTH": snapshot.depth,
                    "motif": snapshot.motif,
                    "bam": snapshot.bam,
                    "snapshot": snapshot.snapshot,
                    "locus": snapshot.locus,
                    "selection_reason": snapshot.selection_reason,
                    "status": snapshot.status,
                }
            )


def write_failed_snapshots(failures: list[dict[str, object]], path: Path) -> None:
    """Write IGV execution failures."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["batch_file", "return_code", "error"],
            delimiter="\t",
        )
        writer.writeheader()
        for failure in failures:
            writer.writerow(failure)


def execute_igv_batches(batch_files: tuple[Path, ...], igv_executable: str) -> list[dict[str, object]]:
    """Execute IGV batch files and return failure records."""

    failures: list[dict[str, object]] = []
    for batch_file in batch_files:
        completed = subprocess.run(
            (igv_executable, "-b", str(batch_file)),
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            failures.append(
                {
                    "batch_file": batch_file,
                    "return_code": completed.returncode,
                    "error": completed.stderr.strip(),
                }
            )
    return failures


def sanitize_snapshot_name(name: str) -> str:
    """Make a snapshot filename safe for common filesystems."""

    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def _default_mutation_table(output_directory: Path) -> Path:
    motif_table = output_directory / "Analysis" / "motifs" / "annotated_motifs.tsv"
    filtered_tables = sorted((output_directory / "Analysis" / "filtered").glob("*.tsv"))
    if motif_table.exists():
        return motif_table
    if filtered_tables:
        return filtered_tables[0]
    return output_directory / "Variants" / "master_mutation_table.tsv"


def _effective_igv_config(loaded_config: LoadedConfig) -> dict[str, object]:
    config = dict(loaded_config.defaults.get("igv", {}))
    assay_igv = loaded_config.assay.get("igv", {})
    if isinstance(assay_igv, dict):
        config.update(assay_igv)
    overrides = loaded_config.run.get("overrides", {})
    if isinstance(overrides, dict) and isinstance(overrides.get("igv"), dict):
        config.update(overrides["igv"])
    return config


def _expected_region_ranges(loaded_config: LoadedConfig) -> list[tuple[int, int]]:
    igv_config = _effective_igv_config(loaded_config)
    expected_ids = igv_config.get("expected_region_ids")
    if not expected_ids:
        expected_ids = ["repair_interval", "gap"]
    if isinstance(expected_ids, str):
        expected_ids = [expected_ids]

    ranges: list[tuple[int, int]] = []
    for region in loaded_config.assay.get("regions", []):
        if isinstance(region, dict) and region.get("id") in expected_ids:
            ranges.append((int(region["start"]), int(region["end"])))
    return ranges


def _position_in_any_region(position: int, ranges: list[tuple[int, int]]) -> bool:
    if not ranges:
        return True
    return any(start <= position <= end for start, end in ranges)
