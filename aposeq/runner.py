"""One-command APO-SEQ pipeline runner."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from aposeq.alignment import process_bams
from aposeq.analysis import run_analysis
from aposeq.config import LoadedConfig
from aposeq.coverage import run_coverage
from aposeq.exceptions import ApoSeqError
from aposeq.igv import run_igv_snapshots
from aposeq.report import write_run_manifest, write_text_report
from aposeq.state import (
    is_completed,
    load_state,
    mark_completed,
    mark_failed,
    mark_skipped,
    mark_started,
    save_state,
    state_path,
)
from aposeq.variants import run_variant_calling


DEFAULT_STEPS = ("process_bams", "coverage", "call_variants", "analyze", "igv_snapshots")


class PipelineStepSkipped(Exception):
    """Raised internally when a dry-run step cannot run without upstream data."""


@dataclass(frozen=True)
class RunnerResult:
    """Outputs from a full APO-SEQ run."""

    output_directory: Path
    state_file: Path
    run_manifest: Path
    report: Path
    completed_steps: tuple[str, ...]
    skipped_steps: tuple[str, ...]
    dry_run: bool


def run_pipeline(
    loaded_config: LoadedConfig,
    dry_run: bool = False,
    resume: bool = True,
    steps: tuple[str, ...] = DEFAULT_STEPS,
    reference_fasta: Path | None = None,
    samtools: str = "samtools",
    bcftools: str = "bcftools",
    igv_executable: str = "igv.sh",
) -> RunnerResult:
    """Run APO-SEQ steps in order with state tracking."""

    output_directory = Path(str(loaded_config.run["output"]["directory"])).expanduser()
    output_directory.mkdir(parents=True, exist_ok=True)
    filename = str(loaded_config.defaults.get("state", {}).get("filename", "pipeline_state.json"))
    path = state_path(output_directory, filename)
    state = load_state(path) if resume else {"steps": {}}
    _hydrate_auto_discovered_samples(loaded_config, output_directory)
    state["run"] = loaded_config.run["name"]
    state["assay"] = loaded_config.assay["name"]
    state["dry_run"] = dry_run
    save_state(state, path)

    completed_steps: list[str] = []
    skipped_steps: list[str] = []
    manifest = write_run_manifest(loaded_config, output_directory, dry_run=dry_run, steps=list(steps))

    for step in steps:
        if resume and is_completed(state, step):
            skipped_steps.append(step)
            continue
        try:
            mark_started(state, path, step)
            _run_step(
                step,
                loaded_config=loaded_config,
                dry_run=dry_run,
                reference_fasta=reference_fasta,
                samtools=samtools,
                bcftools=bcftools,
                igv_executable=igv_executable,
            )
            mark_completed(state, path, step)
            completed_steps.append(step)
        except PipelineStepSkipped:
            mark_skipped(state, path, step)
            skipped_steps.append(step)
        except Exception as exc:
            mark_failed(state, path, step, exc)
            raise

    state = load_state(path)
    report = write_text_report(loaded_config, output_directory, state, dry_run=dry_run)
    return RunnerResult(
        output_directory=output_directory,
        state_file=path,
        run_manifest=manifest,
        report=report,
        completed_steps=tuple(completed_steps),
        skipped_steps=tuple(skipped_steps),
        dry_run=dry_run,
    )


def _run_step(
    step: str,
    loaded_config: LoadedConfig,
    dry_run: bool,
    reference_fasta: Path | None,
    samtools: str,
    bcftools: str,
    igv_executable: str,
) -> None:
    if step == "process_bams":
        process_bams(loaded_config, dry_run=dry_run, samtools=samtools)
    elif step == "coverage":
        run_coverage(loaded_config, dry_run=dry_run, samtools=samtools)
    elif step == "call_variants":
        run_variant_calling(loaded_config, dry_run=dry_run, bcftools=bcftools)
    elif step == "analyze":
        if dry_run and not _default_master_table(loaded_config).exists():
            raise PipelineStepSkipped()
        run_analysis(loaded_config, reference_fasta=reference_fasta)
    elif step == "igv_snapshots":
        if dry_run and not _default_igv_input_table(loaded_config).exists():
            raise PipelineStepSkipped()
        run_igv_snapshots(
            loaded_config,
            dry_run=True if dry_run else False,
            igv_executable=igv_executable,
        )
    else:
        raise ValueError(f"Unknown APO-SEQ run step: {step}")


def _output_directory(loaded_config: LoadedConfig) -> Path:
    return Path(str(loaded_config.run["output"]["directory"])).expanduser()


def _default_master_table(loaded_config: LoadedConfig) -> Path:
    return _output_directory(loaded_config) / "Variants" / "master_mutation_table.tsv"


def _default_igv_input_table(loaded_config: LoadedConfig) -> Path:
    output_directory = _output_directory(loaded_config)
    motif_table = output_directory / "Analysis" / "motifs" / "annotated_motifs.tsv"
    if motif_table.exists():
        return motif_table
    filtered_dir = output_directory / "Analysis" / "filtered"
    filtered_tables = sorted(filtered_dir.glob("*.tsv")) if filtered_dir.exists() else []
    if filtered_tables:
        return filtered_tables[0]
    return _default_master_table(loaded_config)


def _hydrate_auto_discovered_samples(loaded_config: LoadedConfig, output_directory: Path) -> None:
    from aposeq.config import resolve_run_samples

    if loaded_config.run.get("samples"):
        return
    samples = resolve_run_samples(loaded_config.run, output_directory)
    if samples:
        loaded_config.run["samples"] = samples
    elif not loaded_config.run.get("input", {}).get("auto_discover_samples"):
        raise ApoSeqError("Run config has no samples")
