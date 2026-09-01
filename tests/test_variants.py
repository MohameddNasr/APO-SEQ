import gzip
from pathlib import Path

import pytest

from aposeq.config import load_run_config
from aposeq.exceptions import ApoSeqError
from aposeq.variants import build_variant_plans, parse_vcf_to_rows, run_variant_calling, write_mutation_table


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


def test_build_variant_plans_contains_bcftools_commands(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))

    plans = build_variant_plans(loaded)

    assert len(plans) == 1
    assert plans[0].commands[0].argv[:2] == ("bcftools", "mpileup")
    assert plans[0].commands[1].argv[:2] == ("bcftools", "call")
    assert plans[0].commands[2].step == "index_vcf"


def test_run_variant_calling_dry_run_writes_command_manifest(tmp_path: Path) -> None:
    loaded = load_run_config(make_run_config(tmp_path))

    result = run_variant_calling(loaded, dry_run=True)

    assert result.dry_run is True
    assert result.command_manifest_path.exists()
    manifest = result.command_manifest_path.read_text(encoding="utf-8")
    assert "bcftools mpileup" in manifest
    assert "bcftools call" in manifest


def test_parse_vcf_to_rows_extracts_snv_depth_and_af(tmp_path: Path) -> None:
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbarcode01",
                "ebony4.2kb_ref\t10\t.\tC\tT\t60\tPASS\tDP=20\tGT:AD:DP\t0/1:12,8:20",
                "ebony4.2kb_ref\t11\t.\tG\tGA\t60\tPASS\tDP=10\tGT:AD:DP\t0/1:5,5:10",
                "other\t12\t.\tC\tT\t60\tPASS\tDP=10\tGT:AD:DP\t0/1:5,5:10",
            ]
        ),
        encoding="utf-8",
    )

    rows = parse_vcf_to_rows(
        vcf,
        barcode="barcode01",
        genotype="Cas9-negative",
        replicate=1,
        expected_chromosome="ebony4.2kb_ref",
    )

    assert len(rows) == 1
    assert rows[0]["MUTATION"] == "C>T"
    assert rows[0]["DEPTH"] == 20
    assert rows[0]["ALT_COUNT"] == 8
    assert rows[0]["AF"] == "0.400000"


def test_parse_gzipped_vcf_and_write_mutation_table(tmp_path: Path) -> None:
    vcf = tmp_path / "sample.vcf.gz"
    with gzip.open(vcf, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbarcode01",
                    "ebony4.2kb_ref\t10\t.\tG\tA\t60\tPASS\tDP=10\tGT:AD\t0/1:3,7",
                ]
            )
        )

    rows = parse_vcf_to_rows(vcf, "barcode01", "test")
    output = tmp_path / "mutations.tsv"
    write_mutation_table(rows, output)

    text = output.read_text(encoding="utf-8")
    assert "G>A" in text
    assert "0.700000" in text


def test_variant_planning_requires_resolved_samples(tmp_path: Path) -> None:
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
        build_variant_plans(load_run_config(run_config))
