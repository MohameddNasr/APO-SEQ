# APO-SEQ

APO-SEQ is an automated sequencing analysis pipeline for APOBEC-assisted mutation detection from Oxford Nanopore GridION amplicon data.

The v1.0 target workflow starts with unmerged BAM files and produces merged/sorted/indexed BAMs, coverage tables, variant calls, mutation tables, APOBEC motif analysis, QC/statistics, publication figures, IGV snapshots, and a final run report.

## Current Status

This repository currently contains:

- Goal 1: the configuration and metadata foundation.
- Goal 2: BAM discovery and preprocessing command planning/execution.

```bash
python -m aposeq.cli validate-config --run-config config/runs/ebony_cas9_negative.yaml
python -m aposeq.cli process-bams --run-config config/runs/ebony_cas9_negative.yaml --dry-run
python -m aposeq.cli coverage --run-config config/runs/ebony_cas9_negative.yaml --dry-run
python -m aposeq.cli call-variants --run-config config/runs/ebony_cas9_negative.yaml --dry-run
python -m aposeq.cli analyze --run-config config/runs/ebony_cas9_negative.yaml
python -m aposeq.cli igv-snapshots --run-config config/runs/ebony_cas9_negative.yaml --dry-run
python -m aposeq.cli run --run-config config/runs/ebony_cas9_negative.yaml --dry-run
python -m aposeq.cli validate-env
python -m aposeq.cli smoke-test --output-dir /tmp/aposeq_smoke
```

Before running `process-bams`, edit the selected run config so `input.bam_directory` points to the real GridION unmerged BAM directory.

Dry-run commands create manifests and exact command lines without running `samtools` or `bcftools`. This is the safest way to inspect a run before submitting work to Longleaf.

The full runner records progress in `pipeline_state.json`. If a run stops, rerunning the same command resumes by skipping completed steps. Use `--no-resume` to ignore the state file and rerun requested steps.

During a from-scratch `--dry-run`, APO-SEQ plans BAM, coverage, and variant commands. Data-dependent downstream steps such as analysis and IGV snapshots are skipped unless the expected mutation table already exists.

Main planned outputs:

```text
Results/<run>/
  BAM/
    bam_manifest.tsv
    bam_commands.tsv
    merged/
    sorted/
    metrics/
  Coverage/
    coverage_manifest.tsv
    coverage_commands.tsv
    depth/
    summaries/
  Variants/
    variant_commands.tsv
    master_mutation_table.tsv
    raw/
    vcf/
    tables/
  Analysis/
    filtered/
    summaries/
    motifs/
    plot_inputs/
    figures/
  IGV/
    selected_mutations.tsv
    snapshot_manifest.tsv
    failed_snapshots.tsv
    batches/
    snapshots/
  pipeline_state.json
  run_manifest.json
  Report.txt
```

IGV snapshot generation defaults to dry-run/batch generation. Remove `--dry-run` only on a machine or Longleaf job where the IGV executable is available.

The primary manuscript-style figure is written as a 1200 dpi PNG at
`Analysis/figures/apobec3a_locus_map.png` and as a vector PDF at
`Analysis/figures/apobec3a_locus_map.pdf`, with an SVG companion file.

See [docs/GUIDELINE.md](docs/GUIDELINE.md) for the full user-facing workflow, [docs/VALIDATION.md](docs/VALIDATION.md) before running a full experiment, and [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the exact reproducibility record for the delta2-3 transposase-negative P{wa} run.

## Goal Checkpoints

1. Scope and terminology: APO-SEQ v1.0 supports the `ebony` and `pwa` assays. Cas9-negative is a genotype/run under the ebony assay.
2. Configuration layer: default settings, assay metadata, run metadata, execution settings, and validation.
3. BAM automation: discover unmerged BAMs, merge, sort, index, and collect alignment metrics.
4. Variant workflow: run coverage and `bcftools`, then build a master mutation table.
5. Analysis workflow: filtering, QC, statistics, APOBEC motif analysis, and figures.
6. IGV snapshots: generate selected mutation snapshots plus manifests and failure logs.
7. Runner/reporting: one-command CLI, resumability, state files, software versions, and reports.
8. Validation: dry run, small data, then full Longleaf run.
9. Guideline: user-facing automated sequencing analysis guide.

## Configuration Layout

```text
config/
  default.yaml
  execution/
    longleaf.yaml
  assays/
    ebony.yaml
    pwa.yaml
  runs/
    ebony_cas9_negative.yaml
    pwa_example.yaml
```

Assay files describe the biology and coordinates. Run files describe an actual sequencing run, including BAM location, genotype/barcode mapping, output directory, and optional overrides.
