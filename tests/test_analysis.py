from pathlib import Path

import pandas as pd

from aposeq.analysis import run_analysis
from aposeq.cli import main
from aposeq.config import load_run_config
from aposeq.filters import filter_mutations, run_threshold_sweep
from aposeq.motif import annotate_motifs, reverse_complement, summarize_motifs
from aposeq.qc import summarize_mutation_qc
from aposeq.statistics import summarize_by_genotype, summarize_positions
from aposeq.table import save_table


def make_run_config(tmp_path: Path) -> Path:
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "\n".join(
            [
                "name: test_run",
                "assay: ebony",
                "execution: longleaf",
                "input:",
                f"  bam_directory: {tmp_path / 'input'}",
                '  barcode_regex: "(barcode[0-9]+)"',
                "output:",
                f"  directory: {tmp_path / 'Results'}",
                "samples:",
                "  - barcode: barcode01",
                "    genotype: Cas9-negative",
                "    replicate: 1",
                "  - barcode: barcode02",
                "    genotype: edited",
                "    replicate: 1",
                "overrides:",
                "  analysis:",
                "    af_thresholds:",
                "      - 0.5",
                "      - 0.9",
                "    alt_count_thresholds:",
                "      - 10",
            ]
        ),
        encoding="utf-8",
    )
    return run_config


def make_mutation_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "BARCODE": "barcode01",
                "GENOTYPE": "Cas9-negative",
                "REPLICATE": 1,
                "CHROM": "ebony4.2kb_ref",
                "POS": 3,
                "REF": "C",
                "ALT": "T",
                "MUTATION": "C>T",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 20,
                "ALT_COUNT": 12,
                "AF": 0.6,
            },
            {
                "BARCODE": "barcode01",
                "GENOTYPE": "Cas9-negative",
                "REPLICATE": 1,
                "CHROM": "ebony4.2kb_ref",
                "POS": 5,
                "REF": "A",
                "ALT": "G",
                "MUTATION": "A>G",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 20,
                "ALT_COUNT": 15,
                "AF": 0.75,
            },
            {
                "BARCODE": "barcode02",
                "GENOTYPE": "edited",
                "REPLICATE": 1,
                "CHROM": "ebony4.2kb_ref",
                "POS": 7,
                "REF": "G",
                "ALT": "A",
                "MUTATION": "G>A",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 25,
                "ALT_COUNT": 20,
                "AF": 0.8,
            },
            {
                "BARCODE": "barcode02",
                "GENOTYPE": "edited",
                "REPLICATE": 1,
                "CHROM": "ebony4.2kb_ref",
                "POS": 9,
                "REF": "C",
                "ALT": "T",
                "MUTATION": "C>T",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 25,
                "ALT_COUNT": 2,
                "AF": 0.08,
            },
        ]
    )


def test_filter_mutations_applies_af_alt_and_apobec() -> None:
    table = make_mutation_table()

    filtered = filter_mutations(table, 0.5, 10, apobec_only=True)

    assert filtered["MUTATION"].tolist() == ["C>T", "G>A"]


def test_threshold_sweep_counts_mutations() -> None:
    sweep = run_threshold_sweep(
        make_mutation_table(),
        af_thresholds=[0.5, 0.9],
        alt_count_thresholds=[10],
        apobec_only=True,
        deduplicate_positions=True,
    )

    assert sweep["mutation_count"].tolist() == [2, 0]


def test_summaries_are_grouped_for_qc_stats_and_positions() -> None:
    filtered = filter_mutations(make_mutation_table(), 0.5, 10, apobec_only=True)

    qc = summarize_mutation_qc(filtered)
    genotype = summarize_by_genotype(filtered)
    positions = summarize_positions(filtered)

    assert set(qc["BARCODE"]) == {"barcode01", "barcode02"}
    assert set(genotype["GENOTYPE"]) == {"Cas9-negative", "edited"}
    assert positions["mutation_count"].sum() == 2


def test_motif_annotation_handles_tc_and_ga_contexts() -> None:
    table = filter_mutations(make_mutation_table(), 0.5, 10, apobec_only=True)
    fasta_records = {"ebony4.2kb_ref": "ATCGAGAATC"}

    annotated = annotate_motifs(table, fasta_records, "ebony4.2kb_ref", flank=1)
    summary = summarize_motifs(annotated)

    assert annotated["IS_APOBEC_CONTEXT"].tolist() == [True, True]
    assert annotated["MOTIF"].tolist() == ["TC", "GA"]
    assert reverse_complement("TGA") == "TCA"
    assert summary["mutation_count"].sum() == 2


def test_run_analysis_writes_expected_outputs(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))
    mutation_table = tmp_path / "mutations.tsv"
    fasta = tmp_path / "reference.fa"
    save_table(make_mutation_table(), mutation_table)
    fasta.write_text(">ebony4.2kb_ref\nATCGAGAATC\n", encoding="utf-8")

    result = run_analysis(
        loaded,
        mutation_table=mutation_table,
        reference_fasta=fasta,
        motif_flank=1,
    )

    assert result.filtered_table.exists()
    assert result.threshold_sweep.exists()
    assert result.qc_summary.exists()
    assert result.genotype_summary.exists()
    assert result.position_summary.exists()
    assert result.motif_table is not None and result.motif_table.exists()
    assert result.motif_summary is not None and result.motif_summary.exists()
    assert len(result.figures) == 6
    assert all(path.exists() for path in result.figures)
    assert result.figures[0].name == "apobec3a_locus_map.png"
    assert result.figures[1].name == "apobec3a_locus_map.pdf"
    assert result.figures[2].name == "apobec3a_locus_map.svg"


def test_analyze_cli_writes_outputs(tmp_path: Path) -> None:
    run_config = make_run_config(tmp_path)
    mutation_table = tmp_path / "mutations.tsv"
    fasta = tmp_path / "reference.fa"
    save_table(make_mutation_table(), mutation_table)
    fasta.write_text(">ebony4.2kb_ref\nATCGAGAATC\n", encoding="utf-8")

    exit_code = main(
        [
            "analyze",
            "--run-config",
            str(run_config),
            "--mutation-table",
            str(mutation_table),
            "--reference-fasta",
            str(fasta),
            "--motif-flank",
            "1",
        ]
    )

    assert exit_code == 0
    assert (tmp_path / "Results" / "Analysis" / "summaries" / "qc_summary.tsv").exists()
    assert (tmp_path / "Results" / "Analysis" / "figures" / "apobec3a_locus_map.png").exists()
    assert (tmp_path / "Results" / "Analysis" / "figures" / "apobec3a_locus_map.pdf").exists()
    assert (tmp_path / "Results" / "Analysis" / "figures" / "apobec3a_locus_map.svg").exists()
    assert (tmp_path / "Results" / "Analysis" / "figures" / "mutation_positions.svg").exists()


def test_run_analysis_writes_empty_state_figures(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))
    mutation_table = tmp_path / "mutations.tsv"
    save_table(make_mutation_table().iloc[0:0], mutation_table)

    result = run_analysis(loaded, mutation_table=mutation_table)

    assert len(result.figures) == 6
    locus_text = result.figures[2].read_text(encoding="utf-8")
    position_text = result.figures[3].read_text(encoding="utf-8")
    assert "No high-confidence APOBEC3A mutations passed the active filters." in locus_text
    assert "No mutations passed the active filters." in position_text
