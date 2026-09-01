import json
from pathlib import Path

from aposeq.cli import main
from aposeq.validation import create_smoke_test_fixture, run_smoke_test, validate_environment


def test_validate_environment_reports_missing_tools(tmp_path: Path) -> None:
    report = tmp_path / "env.json"

    result = validate_environment(
        output_path=report,
        samtools="definitely_missing_samtools",
        bcftools="definitely_missing_bcftools",
        igv="definitely_missing_igv",
    )

    assert report.exists()
    assert all(not check.available for check in result.checks)
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["checks"][0]["name"] == "samtools"


def test_create_smoke_test_fixture_writes_core_files(tmp_path: Path) -> None:
    fixture = create_smoke_test_fixture(tmp_path / "smoke")

    assert fixture.run_config.exists()
    assert fixture.reference_fasta.exists()
    assert fixture.mutation_table.exists()
    assert fixture.report.exists()
    assert (fixture.root / "input" / "barcode01" / "read_1.bam").exists()


def test_run_smoke_test_writes_runner_outputs(tmp_path: Path) -> None:
    fixture = run_smoke_test(tmp_path / "smoke")

    assert (fixture.root / "Results" / "pipeline_state.json").exists()
    assert (fixture.root / "Results" / "run_manifest.json").exists()
    assert (fixture.root / "Results" / "Report.txt").exists()
    assert (fixture.root / "Results" / "Analysis" / "summaries" / "qc_summary.tsv").exists()
    assert (fixture.root / "Results" / "IGV" / "snapshot_manifest.tsv").exists()


def test_validate_env_cli_writes_report(tmp_path: Path) -> None:
    report = tmp_path / "env.json"

    exit_code = main(
        [
            "validate-env",
            "--output",
            str(report),
            "--samtools",
            "definitely_missing_samtools",
            "--bcftools",
            "definitely_missing_bcftools",
            "--igv",
            "definitely_missing_igv",
        ]
    )

    assert exit_code == 0
    assert report.exists()


def test_smoke_test_cli_runs(tmp_path: Path) -> None:
    output_dir = tmp_path / "smoke"

    exit_code = main(["smoke-test", "--output-dir", str(output_dir)])

    assert exit_code == 0
    assert (output_dir / "Results" / "Report.txt").exists()
