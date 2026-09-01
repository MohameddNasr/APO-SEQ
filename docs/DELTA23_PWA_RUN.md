# Delta2-3 Transposase-Negative P{wa} Run

Dataset:

```text
/path/to/GridION/run/bam_pass
```

Run config:

```text
config/runs/pwa_delta2_3_transposase_negative.yaml
```

The run config uses barcode auto-discovery:

```yaml
auto_discover_samples: true
default_genotype: delta2-3 transposase negative
```

This means every barcode discovered under `bam_pass` will be processed and assigned to the same genotype.

## Heavy Core Job

Submit from the APO-SEQ project directory on Longleaf:

```bash
sbatch Results/pwa_delta2_3_transposase_negative/slurm/pwa_delta2_3_transposase_negative.sbatch
```

This runs:

```text
process_bams
coverage
call_variants
```

## After The Core Job

Check:

```text
Results/pwa_delta2_3_transposase_negative/Report.txt
Results/pwa_delta2_3_transposase_negative/pipeline_state.json
Results/pwa_delta2_3_transposase_negative/Variants/master_mutation_table.tsv
```

Then run analysis and IGV batch generation.

For the full reproducibility record, see:

```text
docs/REPRODUCIBILITY.md
```
