"""End-to-end analysis workflow for APO-SEQ mutation tables."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from aposeq.config import LoadedConfig
from aposeq.figures import write_analysis_figures
from aposeq.filters import filter_mutations, run_threshold_sweep
from aposeq.motif import annotate_motifs, load_fasta, summarize_motifs
from aposeq.qc import summarize_mutation_qc
from aposeq.statistics import summarize_by_genotype, summarize_positions
from aposeq.table import load_mutation_table, save_table


@dataclass(frozen=True)
class AnalysisResult:
    """Paths written by the APO-SEQ analysis workflow."""

    filtered_table: Path
    threshold_sweep: Path
    qc_summary: Path
    genotype_summary: Path
    position_summary: Path
    motif_table: Path | None
    motif_summary: Path | None
    figures: tuple[Path, ...]


def run_analysis(
    loaded_config: LoadedConfig,
    mutation_table: Path | None = None,
    reference_fasta: Path | None = None,
    motif_flank: int = 5,
) -> AnalysisResult:
    """Run filtering, summaries, and optional motif annotation."""

    output_directory = Path(str(loaded_config.run["output"]["directory"])).expanduser()
    mutation_table = mutation_table or output_directory / "Variants" / "master_mutation_table.tsv"
    analysis_root = output_directory / "Analysis"
    filtered_root = analysis_root / "filtered"
    summary_root = analysis_root / "summaries"
    motif_root = analysis_root / "motifs"
    plot_root = analysis_root / "plot_inputs"
    figure_root = analysis_root / "figures"

    for path in (filtered_root, summary_root, motif_root, plot_root, figure_root):
        path.mkdir(parents=True, exist_ok=True)

    table = load_mutation_table(mutation_table)
    analysis_config = _effective_analysis_config(loaded_config)
    af_thresholds = [float(value) for value in analysis_config["af_thresholds"]]
    alt_count_thresholds = [int(value) for value in analysis_config["alt_count_thresholds"]]
    primary_af = af_thresholds[0]
    primary_alt_count = alt_count_thresholds[0]
    apobec_only = bool(analysis_config.get("apobec_only", True))
    deduplicate_positions = bool(analysis_config.get("deduplicate_positions", True))

    filtered = filter_mutations(
        table,
        minimum_af=primary_af,
        minimum_alt_count=primary_alt_count,
        apobec_only=apobec_only,
        deduplicate_positions=deduplicate_positions,
    )
    sweep = run_threshold_sweep(
        table,
        af_thresholds=af_thresholds,
        alt_count_thresholds=alt_count_thresholds,
        apobec_only=apobec_only,
        deduplicate_positions=deduplicate_positions,
    )

    filtered_table = filtered_root / f"mutations_AF{primary_af:g}_ALT{primary_alt_count}.tsv"
    threshold_sweep = summary_root / "threshold_sweep.tsv"
    qc_summary = summary_root / "qc_summary.tsv"
    genotype_summary = summary_root / "genotype_summary.tsv"
    position_summary = plot_root / "position_summary.tsv"

    save_table(filtered, filtered_table)
    save_table(sweep, threshold_sweep)
    qc = summarize_mutation_qc(filtered)
    genotype = summarize_by_genotype(filtered)
    positions = summarize_positions(filtered)
    save_table(qc, qc_summary)
    save_table(genotype, genotype_summary)
    save_table(positions, position_summary)

    motif_table_path = None
    motif_summary_path = None
    if reference_fasta is not None:
        reference = load_fasta(reference_fasta)
        annotated = annotate_motifs(
            filtered,
            fasta_records=reference,
            reference_chromosome=str(loaded_config.assay["reference"]["chromosome"]),
            flank=motif_flank,
        )
        motif_table_path = motif_root / "annotated_motifs.tsv"
        motif_summary_path = motif_root / "motif_summary.tsv"
        save_table(annotated, motif_table_path)
        save_table(summarize_motifs(annotated), motif_summary_path)

    figures = write_analysis_figures(
        loaded_config=loaded_config,
        filtered=filtered,
        threshold_sweep=sweep,
        genotype_summary=genotype,
        position_summary=positions,
        output_directory=figure_root,
    )

    return AnalysisResult(
        filtered_table=filtered_table,
        threshold_sweep=threshold_sweep,
        qc_summary=qc_summary,
        genotype_summary=genotype_summary,
        position_summary=position_summary,
        motif_table=motif_table_path,
        motif_summary=motif_summary_path,
        figures=figures,
    )


def _effective_analysis_config(loaded_config: LoadedConfig) -> dict[str, object]:
    config = dict(loaded_config.defaults.get("analysis", {}))
    overrides = loaded_config.run.get("overrides", {})
    if isinstance(overrides, dict) and isinstance(overrides.get("analysis"), dict):
        config.update(overrides["analysis"])
    assay_analysis = loaded_config.assay.get("analysis", {})
    if "af_thresholds" not in config and isinstance(assay_analysis, dict):
        config["af_thresholds"] = assay_analysis.get("default_af_thresholds", [0.95])
    if "alt_count_thresholds" not in config and isinstance(assay_analysis, dict):
        config["alt_count_thresholds"] = [assay_analysis.get("default_alt_count", 10)]
    return config
