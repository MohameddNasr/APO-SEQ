"""Command-line interface for APO-SEQ."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from aposeq.alignment import process_bams
from aposeq.analysis import run_analysis
from aposeq.config import load_run_config
from aposeq.coverage import run_coverage
from aposeq.exceptions import ApoSeqError
from aposeq.igv import run_igv_snapshots
from aposeq.runner import DEFAULT_STEPS, run_pipeline
from aposeq.slurm import submit_sbatch, write_sbatch_script
from aposeq.validation import create_smoke_test_fixture, run_smoke_test, validate_environment
from aposeq.variants import run_variant_calling


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="aposeq")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser(
        "validate-config",
        help="Validate an APO-SEQ run configuration and its linked assay settings.",
    )
    validate.add_argument(
        "--run-config",
        required=True,
        help="Path or config/runs name for the run YAML file.",
    )
    validate.add_argument(
        "--json",
        action="store_true",
        help="Print loaded config paths as JSON.",
    )

    process = subparsers.add_parser(
        "process-bams",
        help="Discover, merge, sort, and index run BAM files.",
    )
    process.add_argument(
        "--run-config",
        required=True,
        help="Path or config/runs name for the run YAML file.",
    )
    process.add_argument(
        "--dry-run",
        action="store_true",
        help="Write manifests and planned commands without running samtools.",
    )
    process.add_argument(
        "--samtools",
        default="samtools",
        help="samtools executable path or command name.",
    )

    coverage = subparsers.add_parser(
        "coverage",
        help="Run or plan per-base depth and coverage summaries.",
    )
    coverage.add_argument("--run-config", required=True)
    coverage.add_argument("--dry-run", action="store_true")
    coverage.add_argument("--samtools", default="samtools")

    variants = subparsers.add_parser(
        "call-variants",
        help="Run or plan bcftools variant calling and mutation-table creation.",
    )
    variants.add_argument("--run-config", required=True)
    variants.add_argument("--dry-run", action="store_true")
    variants.add_argument("--bcftools", default="bcftools")

    analyze = subparsers.add_parser(
        "analyze",
        help="Filter mutation tables and write APO-SEQ analysis summaries.",
    )
    analyze.add_argument("--run-config", required=True)
    analyze.add_argument(
        "--mutation-table",
        help="Optional mutation table path. Defaults to Variants/master_mutation_table.tsv.",
    )
    analyze.add_argument(
        "--reference-fasta",
        help="Optional reference FASTA for motif-context annotation.",
    )
    analyze.add_argument(
        "--motif-flank",
        type=int,
        default=5,
        help="Bases on each side of the mutation for motif context.",
    )

    igv = subparsers.add_parser(
        "igv-snapshots",
        help="Select mutations and generate IGV snapshot batch files.",
    )
    igv.add_argument("--run-config", required=True)
    igv.add_argument(
        "--mutation-table",
        help="Optional selected/filtered mutation table. Defaults to analysis outputs.",
    )
    igv.add_argument(
        "--dry-run",
        action="store_true",
        help="Write selected mutations, manifests, and batch files without launching IGV.",
    )
    igv.add_argument(
        "--igv",
        default="igv.sh",
        help="IGV executable path or command name.",
    )

    run = subparsers.add_parser(
        "run",
        help="Run the APO-SEQ workflow with state tracking and reporting.",
    )
    run.add_argument("--run-config", required=True)
    run.add_argument("--dry-run", action="store_true")
    run.add_argument(
        "--no-resume",
        action="store_true",
        help="Ignore existing pipeline_state.json and rerun requested steps.",
    )
    run.add_argument(
        "--steps",
        nargs="+",
        choices=DEFAULT_STEPS,
        default=list(DEFAULT_STEPS),
        help="Pipeline steps to run in order.",
    )
    run.add_argument("--reference-fasta")
    run.add_argument("--samtools", default="samtools")
    run.add_argument("--bcftools", default="bcftools")
    run.add_argument("--igv", default="igv.sh")

    env = subparsers.add_parser(
        "validate-env",
        help="Check external APO-SEQ tools such as samtools, bcftools, and IGV.",
    )
    env.add_argument("--output", help="Optional JSON report path.")
    env.add_argument("--samtools", default="samtools")
    env.add_argument("--bcftools", default="bcftools")
    env.add_argument("--igv", default="igv.sh")

    smoke = subparsers.add_parser(
        "smoke-test",
        help="Create and optionally run a tiny APO-SEQ smoke-test fixture.",
    )
    smoke.add_argument("--output-dir", required=True)
    smoke.add_argument(
        "--create-only",
        action="store_true",
        help="Create the fixture without running the dry-run pipeline.",
    )

    sbatch = subparsers.add_parser(
        "write-sbatch",
        help="Write an sbatch script for heavy APO-SEQ steps.",
    )
    sbatch.add_argument("--run-config", required=True)
    sbatch.add_argument(
        "--steps",
        nargs="+",
        choices=DEFAULT_STEPS,
        default=["process_bams", "coverage", "call_variants"],
    )
    sbatch.add_argument("--reference-fasta")
    sbatch.add_argument("--script-path")
    sbatch.add_argument("--submit", action="store_true")

    return parser


def validate_config_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    if args.json:
        payload = asdict(loaded)
        payload["run_config_path"] = str(loaded.run_config_path)
        payload["assay_config_path"] = str(loaded.assay_config_path)
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("APO-SEQ configuration is valid.")
        print(f"Run config: {Path(loaded.run_config_path)}")
        print(f"Assay config: {Path(loaded.assay_config_path)}")
        if loaded.execution:
            print(f"Execution backend: {loaded.execution['backend']}")
    return 0


def process_bams_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    result = process_bams(
        loaded_config=loaded,
        dry_run=args.dry_run,
        samtools=args.samtools,
    )
    mode = "planned" if args.dry_run else "processed"
    print(f"APO-SEQ BAM processing {mode} for {len(result.plans)} sample(s).")
    print(f"BAM manifest: {result.manifest_path}")
    print(f"Command manifest: {result.command_manifest_path}")
    return 0


def coverage_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    result = run_coverage(loaded, dry_run=args.dry_run, samtools=args.samtools)
    mode = "planned" if args.dry_run else "completed"
    print(f"APO-SEQ coverage {mode} for {len(result.plans)} sample(s).")
    print(f"Coverage manifest: {result.summary_manifest_path}")
    print(f"Command manifest: {result.command_manifest_path}")
    return 0


def call_variants_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    result = run_variant_calling(loaded, dry_run=args.dry_run, bcftools=args.bcftools)
    mode = "planned" if args.dry_run else "completed"
    print(f"APO-SEQ variant calling {mode} for {len(result.plans)} sample(s).")
    print(f"Command manifest: {result.command_manifest_path}")
    print(f"Master mutation table: {result.master_mutation_table}")
    return 0


def analyze_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    result = run_analysis(
        loaded,
        mutation_table=Path(args.mutation_table) if args.mutation_table else None,
        reference_fasta=Path(args.reference_fasta) if args.reference_fasta else None,
        motif_flank=args.motif_flank,
    )
    print("APO-SEQ analysis completed.")
    print(f"Filtered mutations: {result.filtered_table}")
    print(f"Threshold sweep: {result.threshold_sweep}")
    print(f"QC summary: {result.qc_summary}")
    print(f"Genotype summary: {result.genotype_summary}")
    print(f"Position summary: {result.position_summary}")
    if result.motif_table:
        print(f"Motif table: {result.motif_table}")
        print(f"Motif summary: {result.motif_summary}")
    for figure in result.figures:
        print(f"Figure: {figure}")
    return 0


def igv_snapshots_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    result = run_igv_snapshots(
        loaded,
        mutation_table=Path(args.mutation_table) if args.mutation_table else None,
        dry_run=args.dry_run,
        igv_executable=args.igv,
    )
    mode = "planned" if args.dry_run else "completed"
    print(f"APO-SEQ IGV snapshots {mode}.")
    print(f"Selected mutations: {result.selected_mutations}")
    print(f"Snapshot manifest: {result.snapshot_manifest}")
    print(f"Failed snapshots: {result.failed_snapshots}")
    print(f"Batch files: {len(result.batch_files)}")
    return 0


def run_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    result = run_pipeline(
        loaded,
        dry_run=args.dry_run,
        resume=not args.no_resume,
        steps=tuple(args.steps),
        reference_fasta=Path(args.reference_fasta) if args.reference_fasta else None,
        samtools=args.samtools,
        bcftools=args.bcftools,
        igv_executable=args.igv,
    )
    mode = "dry run" if args.dry_run else "run"
    print(f"APO-SEQ {mode} completed.")
    print(f"Output directory: {result.output_directory}")
    print(f"State file: {result.state_file}")
    print(f"Run manifest: {result.run_manifest}")
    print(f"Report: {result.report}")
    print(f"Completed steps: {', '.join(result.completed_steps) or 'none'}")
    print(f"Skipped steps: {', '.join(result.skipped_steps) or 'none'}")
    return 0


def validate_env_command(args: argparse.Namespace) -> int:
    result = validate_environment(
        output_path=Path(args.output) if args.output else None,
        samtools=args.samtools,
        bcftools=args.bcftools,
        igv=args.igv,
    )
    for check in result.checks:
        status = "available" if check.available else "missing"
        detail = f" ({check.version})" if check.version else ""
        print(f"{check.name}: {status}{detail}")
    if result.report_path:
        print(f"Environment report: {result.report_path}")
    return 0


def smoke_test_command(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    fixture = (
        create_smoke_test_fixture(output_dir)
        if args.create_only
        else run_smoke_test(output_dir)
    )
    print("APO-SEQ smoke-test fixture is ready.")
    print(f"Root: {fixture.root}")
    print(f"Run config: {fixture.run_config}")
    print(f"Reference FASTA: {fixture.reference_fasta}")
    print(f"Mutation table: {fixture.mutation_table}")
    print(f"Smoke-test notes: {fixture.report}")
    return 0


def write_sbatch_command(args: argparse.Namespace) -> int:
    loaded = load_run_config(args.run_config)
    script = write_sbatch_script(
        loaded,
        script_path=Path(args.script_path) if args.script_path else None,
        steps=tuple(args.steps),
        reference_fasta=Path(args.reference_fasta) if args.reference_fasta else None,
    )
    print(f"SBATCH script: {script.path}")
    print(f"Submit command: {' '.join(script.command)}")
    if args.submit:
        print(submit_sbatch(script))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "validate-config":
            return validate_config_command(args)
        if args.command == "process-bams":
            return process_bams_command(args)
        if args.command == "coverage":
            return coverage_command(args)
        if args.command == "call-variants":
            return call_variants_command(args)
        if args.command == "analyze":
            return analyze_command(args)
        if args.command == "igv-snapshots":
            return igv_snapshots_command(args)
        if args.command == "run":
            return run_command(args)
        if args.command == "validate-env":
            return validate_env_command(args)
        if args.command == "smoke-test":
            return smoke_test_command(args)
        if args.command == "write-sbatch":
            return write_sbatch_command(args)
    except ApoSeqError as exc:
        print(f"aposeq: error: {exc}", file=sys.stderr)
        return 2

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
