# APO-SEQ Roadmap

This file is the pause/resume checklist for APO-SEQ v1.0.

## Goal 1: Configuration And Metadata

Status: complete

Deliverables:

- Package scaffold
- Default configuration
- `ebony` assay configuration
- `pwa` assay configuration
- Example run configurations
- Longleaf execution configuration
- Config validation CLI

Exit criteria:

- `python -m aposeq.cli validate-config --run-config config/runs/ebony_cas9_negative.yaml` succeeds
- `python -m aposeq.cli validate-config --run-config config/runs/pwa_example.yaml` succeeds

## Goal 2: Automated BAM Processing

Status: complete

Deliverables:

- Discover unmerged BAMs under an input directory
- Group BAMs by barcode/sample
- Merge per-barcode BAMs
- Sort merged BAMs
- Index sorted BAMs
- Write BAM manifest
- Write command manifest
- Collect `samtools flagstat` alignment metrics

Exit criteria:

- `python3 -m pytest` succeeds
- `python3 -m aposeq.cli process-bams --run-config <run.yaml> --dry-run` writes manifests when the run config points to an existing BAM directory

## Goal 3: Coverage And Variant Calling

Status: complete

Deliverables:

- Per-base depth tables
- Coverage summaries
- `bcftools mpileup` command generation
- `bcftools call` command generation
- VCF normalization/parsing
- Master mutation table

Exit criteria:

- `python3 -m pytest` succeeds
- `python3 -m aposeq.cli coverage --run-config <run.yaml> --dry-run` writes `Coverage/coverage_commands.tsv`
- `python3 -m aposeq.cli call-variants --run-config <run.yaml> --dry-run` writes `Variants/variant_commands.tsv`
- Non-dry-run variant calling writes `Variants/master_mutation_table.tsv`

## Goal 4: Analysis Modules

Status: complete

Deliverables:

- AF and ALT_COUNT sweeps
- APOBEC-only filtering
- Position deduplication
- QC summaries
- Statistical summaries
- Motif tables
- Publication figure inputs

Exit criteria:

- `python3 -m pytest` succeeds
- `python3 -m aposeq.cli analyze --run-config <run.yaml>` reads `Variants/master_mutation_table.tsv`
- Optional `--reference-fasta` writes APOBEC context and motif summaries

## Goal 5: IGV Snapshots

Status: complete

Deliverables:

- Select mutations for snapshotting
- Generate IGV batch files
- Run IGV in batch mode
- Write snapshot manifest
- Write failed snapshot log

Exit criteria:

- `python3 -m pytest` succeeds
- `python3 -m aposeq.cli igv-snapshots --run-config <run.yaml> --dry-run` writes `IGV/selected_mutations.tsv`
- Dry-run writes one batch file per selected barcode under `IGV/batches/`
- Snapshot manifest links assay, run, genotype, barcode, chromosome, position, ref, alt, AF, ALT_COUNT, DEPTH, motif, BAM, image path, selection reason, and status

## Goal 6: Runner, CLI, And Reporting

Status: complete

Deliverables:

- `aposeq run`
- Dry-run mode
- Resume mode
- `pipeline_state.json`
- Software version manifest
- Final HTML/text report

Exit criteria:

- `python3 -m pytest` succeeds
- `python3 -m aposeq.cli run --run-config <run.yaml> --dry-run` plans upstream command steps
- `pipeline_state.json` records completed, skipped, running, and failed states
- `run_manifest.json` records run/assay/config/software metadata
- `Report.txt` summarizes step status and output roots
- Resume skips already completed steps without erasing completed state

## Goal 7: Validation

Status: complete

Deliverables:

- Config validation tests
- Small local smoke test
- Small Longleaf test
- Full Longleaf validation run
- Bug-fix pass

Exit criteria:

- `python3 -m pytest` succeeds
- `python3 -m aposeq.cli validate-env` reports external tool availability
- `python3 -m aposeq.cli smoke-test --output-dir <dir>` creates and runs a self-contained dry-run fixture
- `docs/VALIDATION.md` documents local, small real-data, and Longleaf validation workflow

## Goal 8: User-Facing Guideline

Status: complete

Deliverables:

- Input preparation guide
- Config editing guide
- Longleaf execution guide
- Output interpretation guide
- IGV snapshot guide
- Troubleshooting section

Exit criteria:

- `docs/GUIDELINE.md` documents the complete APO-SEQ user workflow
- `README.md` links to the guideline
- `python3 -m pytest` succeeds
