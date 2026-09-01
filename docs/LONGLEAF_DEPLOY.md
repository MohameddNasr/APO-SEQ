# Deploy APO-SEQ To Longleaf

This project was created in the local Codex workspace. To run it on Longleaf, deploy it under:

```text
$HOME/APO-SEQ
```

The current dataset run config analyzes only:

```text
/path/to/GridION/run/bam_pass
```

The bootstrap script checks that this folder exists and contains `.bam` files before submitting the Slurm job.

## One-Paste Deployment

Use the generated local file:

```text
longleaf_bootstrap_aposeq.sh
```

Copy or paste it into Longleaf, then run:

```bash
bash $HOME/longleaf_bootstrap_aposeq.sh
```

The script will:

- create `$HOME/APO-SEQ`
- unpack APO-SEQ there
- install the package with `python3 -m pip install --user -e .`
- validate the P{wa} delta2-3 config
- validate external tools
- generate the Slurm script
- submit the heavy core job with `sbatch`

## Monitor

```bash
squeue -u $USER
```

## Inspect After Completion

```bash
cd $HOME/APO-SEQ
cat Results/pwa_delta2_3_transposase_negative/Report.txt
cat Results/pwa_delta2_3_transposase_negative/pipeline_state.json
ls -lh Results/pwa_delta2_3_transposase_negative/Variants/master_mutation_table.tsv
```
