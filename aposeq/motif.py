"""APOBEC motif annotation and motif summaries."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from aposeq.exceptions import ApoSeqError


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""

    return sequence.translate(COMPLEMENT)[::-1].upper()


def load_fasta(path: Path) -> dict[str, str]:
    """Load a small FASTA file into a contig-to-sequence dictionary."""

    if not path.exists():
        raise ApoSeqError(f"Reference FASTA does not exist: {path}")

    records: dict[str, list[str]] = {}
    current_name: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                records[current_name] = []
            elif current_name is None:
                raise ApoSeqError(f"Invalid FASTA: sequence found before header in {path}")
            else:
                records[current_name].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def extract_context(sequence: str, position: int, flank: int) -> str:
    """Extract a 1-based position-centered sequence context."""

    if position < 1 or position > len(sequence):
        raise ApoSeqError(f"Position {position} is outside reference length {len(sequence)}")
    start = max(0, position - flank - 1)
    end = min(len(sequence), position + flank)
    return sequence[start:end].upper()


def classify_apobec_context(ref: str, alt: str, context: str, flank: int) -> tuple[bool, str]:
    """Classify whether a mutation falls in canonical TC/GA APOBEC context."""

    mutation = f"{ref.upper()}>{alt.upper()}"
    center = flank if len(context) > flank else len(context) // 2
    previous_base = context[center - 1] if center - 1 >= 0 else ""
    next_base = context[center + 1] if center + 1 < len(context) else ""

    if mutation == "C>T":
        motif = f"{previous_base}C" if previous_base else "C"
        return motif == "TC", motif
    if mutation == "G>A":
        motif = f"G{next_base}" if next_base else "G"
        return motif == "GA", motif
    return False, ""


def annotate_motifs(
    table: pd.DataFrame,
    fasta_records: dict[str, str],
    reference_chromosome: str,
    flank: int = 5,
) -> pd.DataFrame:
    """Add sequence context and APOBEC motif columns to a mutation table."""

    if reference_chromosome not in fasta_records:
        raise ApoSeqError(f"Reference contig {reference_chromosome!r} not found in FASTA")

    sequence = fasta_records[reference_chromosome]
    annotated = table.copy()
    contexts: list[str] = []
    motifs: list[str] = []
    is_apobec_context: list[bool] = []
    oriented_contexts: list[str] = []

    for row in annotated.itertuples(index=False):
        context = extract_context(sequence, int(row.POS), flank)
        is_context, motif = classify_apobec_context(str(row.REF), str(row.ALT), context, flank)
        mutation = f"{str(row.REF).upper()}>{str(row.ALT).upper()}"
        contexts.append(context)
        motifs.append(motif)
        is_apobec_context.append(is_context)
        oriented_contexts.append(reverse_complement(context) if mutation == "G>A" else context)

    annotated["CONTEXT"] = contexts
    annotated["MOTIF"] = motifs
    annotated["IS_APOBEC_CONTEXT"] = is_apobec_context
    annotated["ORIENTED_CONTEXT"] = oriented_contexts
    return annotated


def summarize_motifs(table: pd.DataFrame) -> pd.DataFrame:
    """Summarize APOBEC motif annotations by genotype and mutation."""

    if table.empty:
        return pd.DataFrame(
            columns=["GENOTYPE", "MUTATION", "MOTIF", "IS_APOBEC_CONTEXT", "mutation_count"]
        )
    return (
        table.groupby(["GENOTYPE", "MUTATION", "MOTIF", "IS_APOBEC_CONTEXT"], dropna=False)
        .size()
        .reset_index(name="mutation_count")
        .sort_values(["GENOTYPE", "MUTATION", "MOTIF"])
        .reset_index(drop=True)
    )
