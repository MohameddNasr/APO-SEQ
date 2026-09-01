from pathlib import Path

import pytest

from aposeq.config import load_run_config
from aposeq.coverage import build_coverage_plans, run_coverage, summarize_depth_file
from aposeq.exceptions import ApoSeqError


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
            ]
        ),
        encoding="utf-8",
    )
    return run_config


def test_build_coverage_plans_uses_sorted_bams_and_reference_contig(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))

    plans = build_coverage_plans(loaded)

    assert len(plans) == 1
    assert plans[0].sorted_bam.name == "barcode01.sorted.bam"
    assert plans[0].command.argv[:3] == ("samtools", "depth", "-aa")
    assert "ebony4.2kb_ref" in plans[0].command.argv


def test_run_coverage_dry_run_writes_manifests(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))

    result = run_coverage(loaded, dry_run=True)

    assert result.dry_run is True
    assert result.command_manifest_path.exists()
    assert result.summary_manifest_path.exists()
    assert "samtools depth" in result.command_manifest_path.read_text(encoding="utf-8")


def test_summarize_depth_file_writes_depth_metrics(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))
    plan = build_coverage_plans(loaded)[0]
    depth_tsv = tmp_path / "barcode01.depth.tsv"
    summary_tsv = tmp_path / "barcode01.coverage_summary.tsv"
    depth_tsv.write_text(
        "\n".join(
            [
                "ebony4.2kb_ref\t1\t0",
                "ebony4.2kb_ref\t2\t5",
                "ebony4.2kb_ref\t3\t10",
            ]
        ),
        encoding="utf-8",
    )

    summarize_depth_file(depth_tsv, summary_tsv, plan)

    summary = summary_tsv.read_text(encoding="utf-8")
    assert "total_positions" in summary
    assert "0.666667" in summary
    assert "5.000000" in summary


def test_coverage_requires_resolved_samples(tmp_path: Path) -> None:
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "\n".join(
            [
                "name: auto_without_manifest",
                "assay: pwa",
                "execution: longleaf",
                "input:",
                f"  bam_directory: {tmp_path / 'input'}",
                "  auto_discover_samples: true",
                "  default_genotype: test",
                "output:",
                f"  directory: {tmp_path / 'Results'}",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ApoSeqError, match="No samples"):
        build_coverage_plans(load_run_config(run_config))
