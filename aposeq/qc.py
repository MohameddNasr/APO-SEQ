"""QC summaries for APO-SEQ mutation tables."""

from __future__ import annotations

import pandas as pd


def summarize_mutation_qc(table: pd.DataFrame) -> pd.DataFrame:
    """Summarize depth, allele frequency, and mutation counts by barcode."""

    if table.empty:
        return pd.DataFrame(
            columns=[
                "BARCODE",
                "GENOTYPE",
                "mutation_count",
                "mean_depth",
                "median_depth",
                "mean_af",
                "max_af",
                "apobec_mutation_count",
            ]
        )

    grouped = table.groupby(["BARCODE", "GENOTYPE"], dropna=False)
    summary = grouped.agg(
        mutation_count=("POS", "count"),
        mean_depth=("DEPTH", "mean"),
        median_depth=("DEPTH", "median"),
        mean_af=("AF", "mean"),
        max_af=("AF", "max"),
    ).reset_index()
    apobec_counts = (
        table[table["MUTATION"].isin({"C>T", "G>A"})]
        .groupby(["BARCODE", "GENOTYPE"], dropna=False)
        .size()
        .reset_index(name="apobec_mutation_count")
    )
    summary = summary.merge(apobec_counts, how="left", on=["BARCODE", "GENOTYPE"])
    summary["apobec_mutation_count"] = summary["apobec_mutation_count"].fillna(0).astype(int)
    return summary
