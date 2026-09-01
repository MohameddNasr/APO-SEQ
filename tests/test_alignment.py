from pathlib import Path

from aposeq.alignment import (
    build_bam_processing_plans,
    discover_bams,
    group_bams_by_sample,
    infer_samples_from_bams,
    process_bams,
)
from aposeq.config import load_run_config


def test_discover_bams_groups_nested_barcode_paths(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    barcode01 = input_dir / "barcode01"
    barcode02 = input_dir / "nested" / "barcode02"
    barcode01.mkdir(parents=True)
    barcode02.mkdir(parents=True)
    (barcode01 / "pass_1.bam").write_text("", encoding="utf-8")
    (barcode01 / "pass_2.bam").write_text("", encoding="utf-8")
    (barcode02 / "pass_1.bam").write_text("", encoding="utf-8")
    (input_dir / "unmatched.bam").write_text("", encoding="utf-8")

    records = discover_bams(input_dir)

    assert [record.barcode for record in records] == ["barcode01", "barcode01", "barcode02"]


def test_build_bam_processing_plans_contains_samtools_commands(tmp_path: Path) -> None:
    bams = (tmp_path / "a.bam", tmp_path / "b.bam")
    samples = [
        {"barcode": "barcode01", "genotype": "Cas9-negative", "replicate": 1},
    ]
    grouped = group_bams_by_sample(
        records=[
            type("Record", (), {"barcode": "barcode01", "path": bams[0]}),
            type("Record", (), {"barcode": "barcode01", "path": bams[1]}),
        ],
        samples=samples,
    )

    plans = build_bam_processing_plans(grouped, tmp_path / "Results")

    assert len(plans) == 1
    assert plans[0].commands[0].argv[:3] == ("samtools", "merge", "-f")
    assert plans[0].commands[1].step == "sort"
    assert plans[0].commands[2].step == "index"
    assert plans[0].commands[3].step == "flagstat"


def test_process_bams_dry_run_writes_manifests(tmp_path: Path) -> None:
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
            ]
        ),
        encoding="utf-8",
    )

    loaded = load_run_config(run_config)
    result = process_bams(loaded, dry_run=True)

    assert result.dry_run is True
    assert result.manifest_path.exists()
    assert result.command_manifest_path.exists()
    assert "barcode01" in result.manifest_path.read_text(encoding="utf-8")
    assert "samtools merge" in result.command_manifest_path.read_text(encoding="utf-8")
    assert "samtools flagstat" in result.command_manifest_path.read_text(encoding="utf-8")


def test_process_bams_can_auto_discover_samples(tmp_path: Path) -> None:
    for barcode in ("barcode01", "barcode02"):
        input_dir = tmp_path / "input" / barcode
        input_dir.mkdir(parents=True)
        (input_dir / "read_1.bam").write_text("", encoding="utf-8")

    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "\n".join(
            [
                "name: auto_run",
                "assay: pwa",
                "execution: longleaf",
                "input:",
                f"  bam_directory: {tmp_path / 'input'}",
                '  barcode_regex: "(barcode[0-9]+)"',
                "  auto_discover_samples: true",
                "  default_genotype: delta2-3 transposase negative",
                "output:",
                f"  directory: {tmp_path / 'Results'}",
            ]
        ),
        encoding="utf-8",
    )

    loaded = load_run_config(run_config)
    result = process_bams(loaded, dry_run=True)

    assert [plan.barcode for plan in result.plans] == ["barcode01", "barcode02"]
    assert all(plan.genotype == "delta2-3 transposase negative" for plan in result.plans)


def test_infer_samples_from_bams_assigns_replicates() -> None:
    records = [
        type("Record", (), {"barcode": "barcode02"}),
        type("Record", (), {"barcode": "barcode01"}),
    ]

    samples = infer_samples_from_bams(records, genotype="test")

    assert [sample["barcode"] for sample in samples] == ["barcode01", "barcode02"]
    assert [sample["replicate"] for sample in samples] == [1, 2]
