# Reference Files

Place assay-specific FASTA references in this directory before running APO-SEQ.

Reference FASTA files are intentionally not tracked by Git because they are often
lab-specific, large, or generated from local instrument/reference databases.

Each assay configuration points to its expected FASTA with:

```yaml
reference:
  fasta: reference_files/your_reference.fa
  chromosome: your_reference_contig_name
  length: 12345
```

The FASTA header, `reference.chromosome`, and BAM contig name must match. You can
check a BAM contig with:

```bash
samtools view -H sample.sorted.bam | grep '^@SQ'
```

Index the reference before variant calling:

```bash
samtools faidx reference_files/your_reference.fa
```
