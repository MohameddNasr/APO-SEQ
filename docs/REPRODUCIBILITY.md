# APO-SEQ Reproducibility Documentation

This document records the reproducible APO-SEQ workflow used for the delta2-3 transposase-negative P{wa} sequencing dataset.

## Dataset

Run label:

```text
pwa_delta2_3_transposase_negative
```

Biological description:

```text
delta2-3 transposase-negative P{wa} assay
```

Input data root:

```text
/path/to/GridION/run/bam_pass
```

Only BAM files underneath this directory are included in this run.

## APO-SEQ Deployment Location

APO-SEQ was deployed on Longleaf under:

```text
$HOME/APO-SEQ
```

All commands below assume:

```bash
cd $HOME/APO-SEQ
```

## Configuration Files

Run configuration:

```text
config/runs/pwa_delta2_3_transposase_negative.yaml
```

Assay configuration:

```text
config/assays/pwa_delta2_3_transposase_negative.yaml
```

The standard P{wa} assay remains in `config/assays/pwa.yaml`. This run uses a
separate control-specific assay config because the delta2-3 transposase-negative
control has a different feature map and no DSB positions.

Execution configuration:

```text
config/execution/longleaf.yaml
```

The run configuration uses barcode auto-discovery:

```yaml
input:
  bam_directory: /path/to/GridION/run/bam_pass
  barcode_regex: "(barcode[0-9]+)"
  auto_discover_samples: true
  default_genotype: delta2-3 transposase negative
```

This means every discovered barcode under the selected `bam_pass` directory is assigned:

```text
delta2-3 transposase negative
```

## Assay Parameters

Assay:

```text
pwa_delta2_3_transposase_negative
```

Reference FASTA configured in the assay file:

```text
reference_files/pwa_reference.fa
```

This FASTA was supplied from the GridION device/local reference memory as:

```text
delta2-3 neg flypatched pwa complete.fasta
```

For the APO-SEQ run it is stored under the project-standard path:

```text
$HOME/APO-SEQ/reference_files/pwa_reference.fa
```

Reference contig:

```text
GridION_pwa_scalloped_actual_copia_actual_flanks
```

Reference length:

```text
37389
```

No DSB positions are configured for this control-specific assay because this is the delta2-3 transposase-negative control.

Regions:

```text
left_outside_break    1-13533
repair_interval      13534-27632
right_outside_break   27633-37389
```

Features:

```text
sd_left       1-13533
p5            13534-14120
w_plus        14121-18133
left_ltr      18134-18409
copia         18410-23003
right_ltr     23004-23279
w_plus_right  23280-27376
p3            27377-27632
sd_right      27633-37389
```

Default analysis thresholds:

```text
AF >= 0.95
ALT_COUNT >= 10
APOBEC-only filtering enabled
```

APOBEC mutation classes:

```text
C>T
G>A
```

APOBEC motif contexts:

```text
TC
GA
```

## External Tools

The heavy sequencing steps are submitted through Slurm `sbatch`.

The core workflow requires:

```text
samtools
bcftools
python3
```

IGV snapshot execution additionally requires:

```text
igv.sh
```

Record tool availability and versions with:

```bash
python3 -m aposeq.cli validate-env \
  --output Results/pwa_delta2_3_transposase_negative/environment_report.json
```

The bootstrap run initially reported `samtools`, `bcftools`, and `igv` as missing in the login-shell environment, but the submitted Slurm script loads:

```bash
module load samtools
module load bcftools
module load minimap2
```

For final publication-grade reproducibility, archive the final `environment_report.json` and Slurm log files.

## Slurm Submission

Generate the Slurm script:

```bash
python3 -m aposeq.cli write-sbatch \
  --run-config config/runs/pwa_delta2_3_transposase_negative.yaml \
  --steps process_bams coverage call_variants
```

Submit the heavy core job:

```bash
sbatch Results/pwa_delta2_3_transposase_negative/slurm/pwa_delta2_3_transposase_negative.sbatch
```

Submitted job observed during this run:

```text
66181888
```

Monitor:

```bash
squeue -u $USER
```

After completion:

```bash
sacct -j 66181888 --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

## Core Workflow

The submitted Slurm job runs these APO-SEQ steps:

```text
process_bams
coverage
call_variants
```

Step details:

1. Discover input BAM files under the configured `bam_pass` directory.
2. Extract barcode names using `(barcode[0-9]+)`.
3. Auto-create sample metadata for discovered barcodes.
4. Merge unmerged BAMs by barcode with `samtools merge`.
5. Sort merged BAMs with `samtools sort`.
6. Index sorted BAMs with `samtools index`.
7. Collect alignment metrics with `samtools flagstat`.
8. Generate per-base depth with `samtools depth`.
9. Run `bcftools mpileup`.
10. Run `bcftools call`.
11. Parse SNVs into the APO-SEQ master mutation table.

## Expected Output Structure

Primary output directory:

```text
Results/pwa_delta2_3_transposase_negative
```

Expected outputs:

```text
Results/pwa_delta2_3_transposase_negative/
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
  logs/
  slurm/
  pipeline_state.json
  run_manifest.json
  Report.txt
```

During the submitted job, early evidence showed:

```text
BAM/bam_manifest.tsv
BAM/bam_commands.tsv
BAM/merged/barcode01.merged.bam through barcode24.merged.bam
BAM/sorted/barcode01.sorted.bam through barcode23.sorted.bam
BAM/sorted/*.bai indexes through barcode23
BAM/metrics/barcode01.flagstat.txt through barcode23.flagstat.txt
```

At that snapshot, the job was still running.

The first submitted job later failed during coverage with Slurm job ID `66181888`:

```text
State: FAILED
ExitCode: 2:0
Elapsed: 00:18:06
```

The APO-SEQ state showed `process_bams` completed and `coverage` failed at barcode01. The failed command requested:

```text
-r pwa_reference
```

The BAM header reported the actual contig:

```text
@SQ SN:GridION_pwa_scalloped_actual_copia_actual_flanks LN:37389
```

`samtools depth` therefore failed with:

```text
samtools depth: cannot parse region "pwa_reference"
```

The reproducible fix is to use the control-specific assay file,
`config/assays/pwa_delta2_3_transposase_negative.yaml`, with:

```yaml
reference:
  name: GridION_pwa_scalloped_actual_copia_actual_flanks
  fasta: reference_files/pwa_reference.fa
  length: 37389
  chromosome: GridION_pwa_scalloped_actual_copia_actual_flanks
```

The right-side control-region and feature annotations should end at `37389`.

Before resubmitting the job, verify the deployed reference FASTA:

```bash
cd $HOME/APO-SEQ
head -1 reference_files/pwa_reference.fa
samtools faidx reference_files/pwa_reference.fa
cut -f1,2 reference_files/pwa_reference.fa.fai
```

Expected FASTA index entry:

```text
GridION_pwa_scalloped_actual_copia_actual_flanks    37389
```

## Completion Checks

After the Slurm job is no longer in `squeue`, run:

```bash
cd $HOME/APO-SEQ
cat Results/pwa_delta2_3_transposase_negative/pipeline_state.json
cat Results/pwa_delta2_3_transposase_negative/Report.txt
tail -100 Results/pwa_delta2_3_transposase_negative/logs/aposeq_pwa_delta2_3_transposase_negative-66181888.err
tail -100 Results/pwa_delta2_3_transposase_negative/logs/aposeq_pwa_delta2_3_transposase_negative-66181888.out
ls -lh Results/pwa_delta2_3_transposase_negative/Variants/master_mutation_table.tsv
```

Successful core completion requires:

```text
pipeline_state.json marks process_bams, coverage, and call_variants as completed
Variants/master_mutation_table.tsv exists
Slurm exit code is 0:0
No fatal errors appear in the .err log
```

## Analysis After Core Job

After `master_mutation_table.tsv` exists, run:

```bash
python3 -m aposeq.cli analyze \
  --run-config config/runs/pwa_delta2_3_transposase_negative.yaml \
  --reference-fasta reference_files/pwa_reference.fa
```

This writes:

```text
Analysis/
  filtered/
  summaries/
  motifs/
  plot_inputs/
  figures/
```

Important analysis outputs:

```text
Analysis/filtered/mutations_AF0.95_ALT10.tsv
Analysis/summaries/threshold_sweep.tsv
Analysis/summaries/qc_summary.tsv
Analysis/summaries/genotype_summary.tsv
Analysis/motifs/annotated_motifs.tsv
Analysis/motifs/motif_summary.tsv
Analysis/plot_inputs/position_summary.tsv
Analysis/figures/apobec3a_locus_map.png  # 1200 dpi PNG
Analysis/figures/apobec3a_locus_map.pdf  # vector PDF
Analysis/figures/apobec3a_locus_map.svg
Analysis/figures/mutation_positions.svg
Analysis/figures/genotype_mutation_counts.svg
Analysis/figures/threshold_sweep.svg
```

## IGV Snapshot Reproducibility

Plan IGV snapshots without launching IGV:

```bash
python3 -m aposeq.cli igv-snapshots \
  --run-config config/runs/pwa_delta2_3_transposase_negative.yaml \
  --dry-run
```

This writes:

```text
IGV/selected_mutations.tsv
IGV/snapshot_manifest.tsv
IGV/failed_snapshots.tsv
IGV/batches/
IGV/snapshots/
```

Snapshot selection settings:

```yaml
minimum_af: 0.95
minimum_alt_count: 10
apobec_only: true
maximum_snapshots: 500
mode:
  - passing_mutations
  - recurrent
  - outside_expected_region
```

The snapshot manifest links each selected mutation to:

```text
assay
run
genotype
barcode
chromosome
position
ref
alt
mutation
AF
ALT_COUNT
DEPTH
motif
bam
snapshot
locus
selection_reason
status
```

## Resume Behavior

APO-SEQ writes:

```text
Results/pwa_delta2_3_transposase_negative/pipeline_state.json
```

If a Slurm job stops, rerun the same `sbatch` command from the APO-SEQ directory. Completed steps are skipped automatically.

To force a rerun of requested steps:

```bash
python3 -m aposeq.cli run \
  --run-config config/runs/pwa_delta2_3_transposase_negative.yaml \
  --steps process_bams coverage call_variants \
  --no-resume
```

Use `--no-resume` cautiously because it can overwrite outputs.

## Files To Archive With The Results

For reproducibility, archive these with the final results:

```text
config/runs/pwa_delta2_3_transposase_negative.yaml
config/assays/pwa_delta2_3_transposase_negative.yaml
config/default.yaml
config/execution/longleaf.yaml
Results/pwa_delta2_3_transposase_negative/run_manifest.json
Results/pwa_delta2_3_transposase_negative/pipeline_state.json
Results/pwa_delta2_3_transposase_negative/environment_report.json
Results/pwa_delta2_3_transposase_negative/BAM/bam_manifest.tsv
Results/pwa_delta2_3_transposase_negative/BAM/bam_commands.tsv
Results/pwa_delta2_3_transposase_negative/Coverage/coverage_commands.tsv
Results/pwa_delta2_3_transposase_negative/Variants/variant_commands.tsv
Results/pwa_delta2_3_transposase_negative/slurm/
Results/pwa_delta2_3_transposase_negative/logs/
```

Also record:

```text
APO-SEQ package version
Git commit or archive checksum, if available
Slurm job ID
Slurm state and exit code
Tool versions for samtools and bcftools
Reference FASTA checksum
```

## Reproducibility Checklist

Before considering the run final:

- Confirm only the intended `bam_pass` directory was used.
- Confirm all expected barcodes are present in `BAM/bam_manifest.tsv`.
- Confirm each barcode has merged, sorted, indexed, and flagstat outputs.
- Confirm coverage files exist for each barcode.
- Confirm VCF files exist for each barcode.
- Confirm `Variants/master_mutation_table.tsv` exists and is non-empty.
- Confirm `pipeline_state.json` marks core steps completed.
- Confirm Slurm exit code is `0:0`.
- Confirm analysis outputs were generated from the final master table.
- Confirm IGV snapshot batches were generated from the filtered or motif-annotated table.
