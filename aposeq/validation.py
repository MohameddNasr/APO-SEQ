"""Validation helpers for APO-SEQ installations and smoke tests."""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from aposeq.config import LoadedConfig
from aposeq.exceptions import ApoSeqError
from aposeq.runner import run_pipeline
from aposeq.table import save_table


@dataclass(frozen=True)
class ToolCheck:
    """Availability and version status for one external tool."""

    name: str
    executable: str
    available: bool
    version: str


@dataclass(frozen=True)
class EnvironmentValidationResult:
    """Result of validating local external tool availability."""

    checks: tuple[ToolCheck, ...]
    report_path: Path | None


@dataclass(frozen=True)
class SmokeTestResult:
    """Files generated for a local APO-SEQ smoke test."""

    root: Path
    run_config: Path
    reference_fasta: Path
    mutation_table: Path
    report: Path


def validate_environment(
    output_path: Path | None = None,
    samtools: str = "samtools",
    bcftools: str = "bcftools",
    igv: str = "igv.sh",
) -> EnvironmentValidationResult:
    """Check whether external tools needed by APO-SEQ are available."""

    checks = (
        check_tool("samtools", samtools, ("--version",)),
        check_tool("bcftools", bcftools, ("--version",)),
        check_tool("igv", igv, ("--version",)),
    )
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(
                {
                    "checks": [
                        {
                            "name": check.name,
                            "executable": check.executable,
                            "available": check.available,
                            "version": check.version,
                        }
                        for check in checks
                    ]
                },
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")
    return EnvironmentValidationResult(checks=checks, report_path=output_path)


def check_tool(name: str, executable: str, version_args: tuple[str, ...]) -> ToolCheck:
    """Check one executable and capture a compact version string."""

    if shutil.which(executable) is None and not Path(executable).exists():
        return ToolCheck(name=name, executable=executable, available=False, version="")
    try:
        completed = subprocess.run(
            (executable, *version_args),
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except OSError:
        return ToolCheck(name=name, executable=executable, available=False, version="")
    output = (completed.stdout or completed.stderr).strip().splitlines()
    version = output[0] if output else ""
    return ToolCheck(name=name, executable=executable, available=completed.returncode == 0, version=version)


def create_smoke_test_fixture(root: Path) -> SmokeTestResult:
    """Create a tiny local APO-SEQ fixture for dry-run and analysis validation."""

    root.mkdir(parents=True, exist_ok=True)
    input_dir = root / "input" / "barcode01"
    input_dir.mkdir(parents=True, exist_ok=True)
    (input_dir / "read_1.bam").write_text("", encoding="utf-8")
    (input_dir / "read_2.bam").write_text("", encoding="utf-8")

    reference_fasta = root / "reference.fa"
    reference_fasta.write_text(">ebony4.2kb_ref\nATCGAGAATC\n", encoding="utf-8")

    mutation_table = root / "Results" / "Variants" / "master_mutation_table.tsv"
    save_table(_smoke_mutation_table(), mutation_table)

    run_config = root / "smoke_run.yaml"
    run_config.write_text(
        "\n".join(
            [
                "name: smoke_test",
                "assay: ebony",
                "execution: longleaf",
                "input:",
                f"  bam_directory: {root / 'input'}",
                '  barcode_regex: "(barcode[0-9]+)"',
                "output:",
                f"  directory: {root / 'Results'}",
                "samples:",
                "  - barcode: barcode01",
                "    genotype: smoke",
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

    report = root / "SMOKE_TEST.md"
    report.write_text(
        "\n".join(
            [
                "# APO-SEQ Smoke Test Fixture",
                "",
                "This fixture is intended for dry-run planning and analysis-layer validation.",
                "The BAM files are placeholders and are not valid sequencing files.",
                "",
                "Run:",
                "",
                "```bash",
                f"python3 -m aposeq.cli run --run-config {run_config} --dry-run --reference-fasta {reference_fasta}",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return SmokeTestResult(
        root=root,
        run_config=run_config,
        reference_fasta=reference_fasta,
        mutation_table=mutation_table,
        report=report,
    )


def run_smoke_test(root: Path, loaded_config: LoadedConfig | None = None) -> SmokeTestResult:
    """Create and run the smoke-test fixture through the dry-run pipeline."""

    fixture = create_smoke_test_fixture(root)
    if loaded_config is not None:
        raise ApoSeqError("run_smoke_test does not accept preloaded external configs")
    from aposeq.config import load_run_config

    run_pipeline(
        load_run_config(fixture.run_config),
        dry_run=True,
        resume=False,
        reference_fasta=fixture.reference_fasta,
    )
    return fixture


def _smoke_mutation_table():
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "BARCODE": "barcode01",
                "GENOTYPE": "smoke",
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
                "GENOTYPE": "smoke",
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
        ]
    )
