# APO-SEQ Validation Guide

This guide describes the checks to run before a full APO-SEQ analysis.

## 1. Validate Python Tests

Run this after any code change:

```bash
python3 -m pytest
```

Expected result:

```text
33 passed
```

## 2. Validate External Tools

Run:

```bash
python3 -m aposeq.cli validate-env --output Results/environment_report.json
```

This checks:

- `samtools`
- `bcftools`
- `igv.sh`

IGV may be missing on a compute node if snapshots are generated elsewhere. `samtools` and `bcftools` must be available before real BAM processing and variant calling.

## 3. Create A Local Smoke-Test Fixture

Run:

```bash
python3 -m aposeq.cli smoke-test --output-dir /tmp/aposeq_smoke
```

This creates:

```text
/tmp/aposeq_smoke/
  smoke_run.yaml
  reference.fa
  input/barcode01/
  Results/
```

The BAM files in this fixture are placeholders. The smoke test validates config parsing, dry-run planning, analysis tables, IGV batch generation, state tracking, and reporting. It does not validate real BAM parsing.

## 4. Dry-Run A Real Run

Edit a run config so:

```yaml
input:
  bam_directory: /path/to/real/GridION/unmerged_bams
```

Then run:

```bash
python3 -m aposeq.cli run \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --dry-run
```

Check:

- `BAM/bam_manifest.tsv`
- `BAM/bam_commands.tsv`
- `Coverage/coverage_commands.tsv`
- `Variants/variant_commands.tsv`
- `pipeline_state.json`
- `run_manifest.json`
- `Report.txt`

## 5. Small Real-Data Test

Before the full experiment, test one barcode or a small subset of BAMs.

Recommended command:

```bash
python3 -m aposeq.cli run \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --steps process_bams coverage call_variants
```

Then run analysis after confirming the master mutation table:

```bash
python3 -m aposeq.cli analyze \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --reference-fasta reference_files/ebony4.2kb_ref.flyPatched.fa
```

Then plan IGV snapshots:

```bash
python3 -m aposeq.cli igv-snapshots \
  --run-config config/runs/ebony_cas9_negative.yaml \
  --dry-run
```

## 6. Longleaf Validation

On Longleaf:

1. Load or activate the environment containing APO-SEQ, `samtools`, and `bcftools`.
2. Run `validate-env`.
3. Run a dry-run using the real run config.
4. Submit a small one-barcode test.
5. Inspect `Report.txt`, `pipeline_state.json`, command manifests, and the master mutation table.
6. Submit the full run.

If a job stops, rerun the same command. APO-SEQ resumes from `pipeline_state.json` and skips completed steps.
