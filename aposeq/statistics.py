"""Statistical summaries for APO-SEQ mutation tables."""

from __future__ import annotations

import pandas as pd


def summarize_by_genotype(table: pd.DataFrame) -> pd.DataFrame:
    """Summarize mutation burden by genotype."""

    if table.empty:
        return pd.DataFrame(
            columns=[
                "GENOTYPE",
                "barcode_count",
                "mutation_count",
                "unique_position_count",
                "mean_mutations_per_barcode",
                "apobec_mutation_count",
            ]
        )

    rows: list[dict[str, object]] = []
    for genotype, group in table.groupby("GENOTYPE", dropna=False):
        barcode_count = group["BARCODE"].nunique()
        mutation_count = len(group)
        rows.append(
            {
                "GENOTYPE": genotype,
                "barcode_count": barcode_count,
                "mutation_count": mutation_count,
                "unique_position_count": group["POS"].nunique(),
                "mean_mutations_per_barcode": mutation_count / barcode_count if barcode_count else 0,
                "apobec_mutation_count": group["MUTATION"].isin({"C>T", "G>A"}).sum(),
            }
        )
    return pd.DataFrame(rows).sort_values("GENOTYPE").reset_index(drop=True)


def summarize_positions(table: pd.DataFrame) -> pd.DataFrame:
    """Create plot-ready per-position mutation counts."""

    if table.empty:
        return pd.DataFrame(
            columns=["CHROM", "POS", "REF", "ALT", "MUTATION", "mutation_count", "barcode_count"]
        )

    return (
        table.groupby(["CHROM", "POS", "REF", "ALT", "MUTATION"], dropna=False)
        .agg(
            mutation_count=("POS", "count"),
            barcode_count=("BARCODE", "nunique"),
            mean_af=("AF", "mean"),
            max_af=("AF", "max"),
            mean_depth=("DEPTH", "mean"),
        )
        .reset_index()
        .sort_values(["CHROM", "POS", "MUTATION"])
        .reset_index(drop=True)
    )
