"""Mutation filtering helpers."""

from __future__ import annotations

import pandas as pd


APOBEC_MUTATIONS = {"C>T", "G>A"}


def filter_mutations(
    table: pd.DataFrame,
    minimum_af: float,
    minimum_alt_count: int,
    apobec_only: bool = False,
    deduplicate_positions: bool = False,
) -> pd.DataFrame:
    """Filter mutations by AF, ALT_COUNT, mutation class, and position."""

    filtered = table[
        (table["AF"] >= minimum_af)
        & (table["ALT_COUNT"] >= minimum_alt_count)
    ].copy()

    if apobec_only:
        filtered = filtered[filtered["MUTATION"].isin(APOBEC_MUTATIONS)].copy()

    if deduplicate_positions and not filtered.empty:
        sort_columns = ["BARCODE", "POS", "AF", "ALT_COUNT"]
        filtered = filtered.sort_values(sort_columns, ascending=[True, True, False, False])
        filtered = filtered.drop_duplicates(subset=["BARCODE", "POS", "REF", "ALT"], keep="first")

    return filtered.reset_index(drop=True)


def run_threshold_sweep(
    table: pd.DataFrame,
    af_thresholds: list[float],
    alt_count_thresholds: list[int],
    apobec_only: bool,
    deduplicate_positions: bool,
) -> pd.DataFrame:
    """Summarize mutation counts across AF and ALT_COUNT thresholds."""

    rows: list[dict[str, object]] = []
    for af in af_thresholds:
        for alt_count in alt_count_thresholds:
            filtered = filter_mutations(
                table,
                minimum_af=float(af),
                minimum_alt_count=int(alt_count),
                apobec_only=apobec_only,
                deduplicate_positions=deduplicate_positions,
            )
            rows.append(
                {
                    "minimum_af": af,
                    "minimum_alt_count": alt_count,
                    "apobec_only": apobec_only,
                    "deduplicate_positions": deduplicate_positions,
                    "mutation_count": len(filtered),
                    "barcode_count": filtered["BARCODE"].nunique() if not filtered.empty else 0,
                    "genotype_count": filtered["GENOTYPE"].nunique() if not filtered.empty else 0,
                }
            )
    return pd.DataFrame(rows)
