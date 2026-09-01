from pathlib import Path

from aposeq.cli import main
from aposeq.config import load_run_config
from aposeq.runner import run_pipeline
from aposeq.table import save_table

from tests.test_analysis import make_mutation_table


def make_run_config(tmp_path: Path) -> Path:
    input_dir = tmp_path / "input" / "barcode01"
    input_dir.mkdir(parents=True)
    (input_dir / "read_1.bam").write_text("", encoding="utf-8")
    (input_dir / "read_2.bam").write_text("", encoding="utf-8")

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
                "overrides:",
                "  analysis:",
                "    af_thresholds:",
                "      - 0.5",
                "    alt_count_thresholds:",
                "      - 10",
                "  igv:",
                "    selection:",
                "      minimum_af: 0.5",
                "      minimum_alt_count: 10",
                "      apobec_only: true",
            ]
        ),
        encoding="utf-8",
    )
    return run_config


def prepare_master_table(tmp_path: Path) -> None:
    path = tmp_path / "Results" / "Variants" / "master_mutation_table.tsv"
    save_table(make_mutation_table(), path)


def test_run_pipeline_dry_run_writes_state_manifest_and_report(tmp_path: Path) -> None:
    run_config = make_run_config(tmp_path)
    prepare_master_table(tmp_path)
    loaded = load_run_config(run_config)

    result = run_pipeline(loaded, dry_run=True)

    assert result.state_file.exists()
    assert result.run_manifest.exists()
    assert result.report.exists()
    assert result.completed_steps == (
        "process_bams",
        "coverage",
        "call_variants",
        "analyze",
        "igv_snapshots",
    )
    assert (tmp_path / "Results" / "Analysis" / "summaries" / "qc_summary.tsv").exists()


def test_run_pipeline_full_dry_run_skips_data_dependent_steps_without_master_table(tmp_path: Path) -> None:
    run_config = make_run_config(tmp_path)
    loaded = load_run_config(run_config)

    result = run_pipeline(loaded, dry_run=True)

    assert result.completed_steps == ("process_bams", "coverage", "call_variants")
    assert result.skipped_steps == ("analyze", "igv_snapshots")
    assert result.report.exists()


def test_run_pipeline_resume_skips_completed_steps_without_overwriting_state(tmp_path: Path) -> None:
    run_config = make_run_config(tmp_path)
    prepare_master_table(tmp_path)
    loaded = load_run_config(run_config)

    first = run_pipeline(loaded, dry_run=True, steps=("coverage",))
    second = run_pipeline(loaded, dry_run=True, steps=("coverage",))

    assert first.completed_steps == ("coverage",)
    assert second.completed_steps == tuple()
    assert second.skipped_steps == ("coverage",)
    assert '"status": "completed"' in second.state_file.read_text(encoding="utf-8")


def test_run_pipeline_resume_recovers_auto_discovered_samples_from_manifest(tmp_path: Path) -> None:
    input_dir = tmp_path / "input" / "barcode01"
    input_dir.mkdir(parents=True)
    (input_dir / "read_1.bam").write_text("", encoding="utf-8")
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "\n".join(
            [
                "name: auto_resume",
                "assay: ebony",
                "execution: longleaf",
                "input:",
                f"  bam_directory: {tmp_path / 'input'}",
                '  barcode_regex: "(barcode[0-9]+)"',
                "  auto_discover_samples: true",
                "  default_genotype: Cas9-negative",
                "output:",
                f"  directory: {tmp_path / 'Results'}",
            ]
        ),
        encoding="utf-8",
    )
    loaded = load_run_config(run_config)
    first = run_pipeline(loaded, dry_run=True, steps=("process_bams",))

    reloaded = load_run_config(run_config)
    second = run_pipeline(reloaded, dry_run=True, steps=("process_bams", "coverage"))

    assert first.completed_steps == ("process_bams",)
    assert second.skipped_steps == ("process_bams",)
    assert second.completed_steps == ("coverage",)


def test_run_cli_dry_run_writes_report(tmp_path: Path) -> None:
    run_config = make_run_config(tmp_path)
    prepare_master_table(tmp_path)

    exit_code = main(["run", "--run-config", str(run_config), "--dry-run"])

    assert exit_code == 0
    assert (tmp_path / "Results" / "Report.txt").exists()
