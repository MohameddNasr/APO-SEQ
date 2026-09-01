from pathlib import Path

import pandas as pd

from aposeq.cli import main
from aposeq.config import load_run_config
from aposeq.igv import (
    build_snapshot_manifest,
    run_igv_snapshots,
    select_mutations_for_igv,
    sanitize_snapshot_name,
    write_igv_batch_files,
)
from aposeq.table import save_table


def make_run_config(tmp_path: Path, assay: str = "pwa") -> Path:
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "\n".join(
            [
                "name: test_run",
                f"assay: {assay}",
                "execution: longleaf",
                "input:",
                f"  bam_directory: {tmp_path / 'input'}",
                '  barcode_regex: "(barcode[0-9]+)"',
                "output:",
                f"  directory: {tmp_path / 'Results'}",
                "samples:",
                "  - barcode: barcode01",
                "    genotype: control",
                "    replicate: 1",
                "  - barcode: barcode02",
                "    genotype: edited",
                "    replicate: 1",
                "overrides:",
                "  igv:",
                "    selection:",
                "      mode:",
                "        - passing_mutations",
                "        - recurrent",
                "        - outside_expected_region",
                "      minimum_af: 0.5",
                "      minimum_alt_count: 10",
                "      apobec_only: true",
                "      maximum_snapshots: 10",
                "    expected_region_ids:",
                "      - repair_interval",
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
                "GENOTYPE": "control",
                "REPLICATE": 1,
                "CHROM": "pwa_reference",
                "POS": 14000,
                "REF": "C",
                "ALT": "T",
                "MUTATION": "C>T",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 20,
                "ALT_COUNT": 12,
                "AF": 0.6,
                "MOTIF": "TC",
            },
            {
                "BARCODE": "barcode02",
                "GENOTYPE": "edited",
                "REPLICATE": 1,
                "CHROM": "pwa_reference",
                "POS": 14000,
                "REF": "C",
                "ALT": "T",
                "MUTATION": "C>T",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 30,
                "ALT_COUNT": 20,
                "AF": 0.7,
                "MOTIF": "TC",
            },
            {
                "BARCODE": "barcode02",
                "GENOTYPE": "edited",
                "REPLICATE": 1,
                "CHROM": "pwa_reference",
                "POS": 1000,
                "REF": "G",
                "ALT": "A",
                "MUTATION": "G>A",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 25,
                "ALT_COUNT": 22,
                "AF": 0.88,
                "MOTIF": "GA",
            },
            {
                "BARCODE": "barcode02",
                "GENOTYPE": "edited",
                "REPLICATE": 1,
                "CHROM": "pwa_reference",
                "POS": 15000,
                "REF": "A",
                "ALT": "G",
                "MUTATION": "A>G",
                "QUAL": 60,
                "FILTER": "PASS",
                "DEPTH": 25,
                "ALT_COUNT": 22,
                "AF": 0.88,
                "MOTIF": "",
            },
        ]
    )


def test_select_mutations_for_igv_marks_recurrent_and_outside_region(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))

    selected = select_mutations_for_igv(make_mutation_table(), loaded)

    assert len(selected) == 3
    assert "recurrent" in selected["SELECTION_REASON"].str.cat(sep=",")
    assert "outside_expected_region" in selected["SELECTION_REASON"].str.cat(sep=",")
    assert "A>G" not in selected["MUTATION"].tolist()


def test_build_snapshot_manifest_uses_sorted_bam_and_flanked_locus(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))
    selected = select_mutations_for_igv(make_mutation_table(), loaded)

    snapshots = build_snapshot_manifest(selected, loaded, tmp_path / "Results")

    assert snapshots[0].bam.name.endswith(".sorted.bam")
    assert snapshots[0].snapshot.suffix == ".png"
    assert snapshots[0].locus.startswith("GridION_pwa_scalloped_actual_copia_actual_flanks:")


def test_write_igv_batch_files_contains_genome_load_goto_and_snapshot(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))
    selected = select_mutations_for_igv(make_mutation_table(), loaded)
    snapshots = build_snapshot_manifest(selected, loaded, tmp_path / "Results")

    batch_files = write_igv_batch_files(
        snapshots,
        loaded,
        tmp_path / "Results" / "IGV" / "batches",
        tmp_path / "Results" / "IGV" / "snapshots",
    )

    assert len(batch_files) == 2
    batch_text = batch_files[0].read_text(encoding="utf-8")
    assert "genome reference_files/pwa_reference.fa" in batch_text
    assert "load " in batch_text
    assert "goto GridION_pwa_scalloped_actual_copia_actual_flanks:" in batch_text
    assert "snapshot " in batch_text


def test_run_igv_snapshots_dry_run_writes_outputs(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))
    table_path = tmp_path / "mutations.tsv"
    save_table(make_mutation_table(), table_path)

    result = run_igv_snapshots(loaded, mutation_table=table_path, dry_run=True)

    assert result.selected_mutations.exists()
    assert result.snapshot_manifest.exists()
    assert result.failed_snapshots.exists()
    assert len(result.batch_files) == 2


def test_igv_cli_dry_run_writes_outputs(tmp_path: Path) -> None:
    run_config = make_run_config(tmp_path)
    table_path = tmp_path / "mutations.tsv"
    save_table(make_mutation_table(), table_path)

    exit_code = main(
        [
            "igv-snapshots",
            "--run-config",
            str(run_config),
            "--mutation-table",
            str(table_path),
            "--dry-run",
        ]
    )

    assert exit_code == 0
    assert (tmp_path / "Results" / "IGV" / "snapshot_manifest.tsv").exists()


def test_sanitize_snapshot_name() -> None:
    assert sanitize_snapshot_name("chr:1 C>T.png") == "chr_1_C_T.png"
