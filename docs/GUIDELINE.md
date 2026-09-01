# APO-SEQ Automated Sequencing Analysis Guideline

This guideline describes how to run APO-SEQ from unmerged Oxford Nanopore GridION BAM files through mutation tables, QC/statistics, APOBEC motif outputs, IGV snapshot batches, and final reports.

## 1. Inputs

APO-SEQ expects one run configuration and one assay configuration.

The run configuration tells APO-SEQ where the sequencing data are:

```yaml
name: ebony_cas9_negative
assay: ebony
execution: longleaf

input:
  bam_directory: /path/to/GridION/unmerged_bams
  barcode_regex: "(barcode[0-9]+)"

output:
  directory: Results/ebony_cas9_negative

samples:
  - barcode: barcode01
    genotype: Cas9-negative
    replicate: 1
```

The assay configuration tells APO-SEQ the biological coordinate system:

- reference FASTA
- reference contig name
- reference length
- assay regions
- DSB positions or repair boundaries
- plotting/shading defaults
- IGV selection defaults

Current assays:

- `ebony`
- `pwa`

Cas9-negative is treated as a genotype/sample group within an `ebony` run, not as a separate assay.

## 2. Prepare The Run Config

Start from one of these files:

```text
config/runs/ebony_cas9_negative.yaml
config/runs/pwa_example.yaml
```

Edit:

- `name`
- `assay`
- `input.bam_directory`
- `output.directory`
- `samples`
- optional AF/ALT_COUNT overrides

The barcode names in `samples` must match the barcode names found in the BAM paths, usually `barcode01`, `barcode02`, and so on.

## 3. Validate The Config

Run:

```bash
python3 -m aposeq.cli validate-config \
  --run-config config/runs/ebony_cas9_negative.yaml
```

This catches invalid YAML structure, missing samples, duplicate barcodes, bad assay names, and invalid coordinate definitions before any sequencing data are processed.

## 4. Validate The Environment

Run:

```bash
python3 -m aposeq.cli validate-env \
  --output Results/environment_report.json
```

For real data processing, `samtools` and `bcftools` must be available. IGV is required only when executing snapshot generation, not when writing dry-run batch files.

## 5. Dry-Run The Full Pipeline

Run:

```bash
python3 -m aposeq.cli run \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --dry-run
```

Dry-run mode writes manifests and planned commands without running `samtools`, `bcftools`, or IGV.

Inspect:

```text
Results/<run>/
  BAM/bam_manifest.tsv
  BAM/bam_commands.tsv
  Coverage/coverage_commands.tsv
  Variants/variant_commands.tsv
  pipeline_state.json
  run_manifest.json
  Report.txt
```

If this is a fresh run with no mutation table yet, APO-SEQ will plan BAM, coverage, and variant-calling steps, then skip analysis and IGV snapshots until real variant outputs exist.

## 6. Run The Core Pipeline

For a real first pass, run:

```bash
python3 -m aposeq.cli run \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --steps process_bams coverage call_variants
```

This performs:

```text
unmerged BAMs
  -> merged BAMs
  -> sorted BAMs
  -> BAM indexes
  -> flagstat metrics
  -> per-base depth
  -> bcftools variant calls
  -> master mutation table
```

Primary output:

```text
Results/<run>/Variants/master_mutation_table.tsv
```

## 7. Run Analysis

After the master mutation table exists, run:

```bash
python3 -m aposeq.cli analyze \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --reference-fasta reference_files/ebony4.2kb_ref.flyPatched.fa
```

This writes:

```text
Results/<run>/Analysis/
  filtered/
  summaries/
  motifs/
  plot_inputs/
  figures/
```

Key files:

- `filtered/mutations_AF<value>_ALT<count>.tsv`
- `summaries/threshold_sweep.tsv`
- `summaries/qc_summary.tsv`
- `summaries/genotype_summary.tsv`
- `motifs/annotated_motifs.tsv`
- `motifs/motif_summary.tsv`
- `plot_inputs/position_summary.tsv`
- `figures/apobec3a_locus_map.png`
- `figures/apobec3a_locus_map.pdf`
- `figures/apobec3a_locus_map.svg`
- `figures/mutation_positions.svg`
- `figures/genotype_mutation_counts.svg`
- `figures/threshold_sweep.svg`

## 8. Generate IGV Snapshot Batches

First run IGV snapshots in dry-run mode:

```bash
python3 -m aposeq.cli igv-snapshots \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --dry-run
```

This writes:

```text
Results/<run>/IGV/
  selected_mutations.tsv
  snapshot_manifest.tsv
  failed_snapshots.tsv
  batches/
  snapshots/
```

Selection can include:

- passing mutations
- recurrent mutations
- mutations outside the expected repair region
- APOBEC-only mutations
- AF and ALT_COUNT thresholds
- maximum snapshot limits

To execute IGV batch files, remove `--dry-run` on a machine where IGV is installed:

```bash
python3 -m aposeq.cli igv-snapshots \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --igv /path/to/igv.sh
```

## 9. Resume Interrupted Runs

APO-SEQ writes:

```text
Results/<run>/pipeline_state.json
```

If a job stops, rerun the same command. Completed steps are skipped automatically.

To force rerun:

```bash
python3 -m aposeq.cli run \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --no-resume
```

## 10. Longleaf Workflow

Recommended Longleaf order:

1. Activate the APO-SEQ environment.
2. Validate the config.
3. Validate external tools.
4. Run a full dry-run.
5. Run one barcode or a small subset.
6. Inspect command manifests, `Report.txt`, and `master_mutation_table.tsv`.
7. Run the full experiment.
8. Run analysis with the correct reference FASTA.
9. Generate IGV batches in dry-run mode.
10. Execute IGV snapshots only where IGV is available.

## 11. Output Interpretation

Use `Report.txt` for a quick run overview.

Use `run_manifest.json` for reproducibility metadata:

- run name
- assay name
- config paths
- dry-run status
- Python/platform metadata
- requested steps

Use `pipeline_state.json` to inspect whether each step completed, failed, or was skipped.

Use `BAM/bam_manifest.tsv` to confirm barcode-to-BAM grouping.

Use `Variants/master_mutation_table.tsv` as the core variant evidence table.

Use `Analysis/summaries/threshold_sweep.tsv` to compare AF and ALT_COUNT thresholds.

Use `Analysis/motifs/annotated_motifs.tsv` for APOBEC TC/GA context interpretation.

Use `IGV/snapshot_manifest.tsv` to connect each planned snapshot back to its mutation evidence.

## 12. Troubleshooting

If config validation fails, check:

- assay name
- duplicate barcode names
- missing `samples`
- invalid reference length
- assay regions extending beyond reference length

If no BAMs are found, check:

- `input.bam_directory`
- `input.barcode_regex`
- whether BAM files end in `.bam`
- whether barcode names in file paths match the sample table

If variant calling fails, check:

- `samtools` and `bcftools` availability
- sorted BAM paths
- reference FASTA path
- reference contig names
- permissions on the output directory

If analysis fails, check:

- `Variants/master_mutation_table.tsv` exists
- required columns are present
- `AF`, `ALT_COUNT`, `DEPTH`, and `POS` are numeric

If IGV snapshot generation fails, check:

- IGV executable path
- sorted BAM files exist
- reference FASTA is accessible
- selected mutations are not empty
- `IGV/failed_snapshots.tsv`

## 13. Minimal Command Sequence

```bash
python3 -m aposeq.cli validate-config --run-config config/runs/ebony_cas9_negative.yaml
python3 -m aposeq.cli validate-env
python3 -m aposeq.cli run --run-config config/runs/ebony_cas9_negative.yaml --dry-run
python3 -m aposeq.cli run --run-config config/runs/ebony_cas9_negative.yaml --steps process_bams coverage call_variants
python3 -m aposeq.cli analyze --run-config config/runs/ebony_cas9_negative.yaml --reference-fasta reference_files/ebony4.2kb_ref.flyPatched.fa
python3 -m aposeq.cli igv-snapshots --run-config config/runs/ebony_cas9_negative.yaml --dry-run
```
