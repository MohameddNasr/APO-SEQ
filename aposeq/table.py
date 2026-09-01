"""Mutation-table loading and validation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from aposeq.exceptions import ApoSeqError


REQUIRED_MUTATION_COLUMNS = (
    "BARCODE",
    "GENOTYPE",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "MUTATION",
    "DEPTH",
    "ALT_COUNT",
    "AF",
)


def load_mutation_table(path: Path) -> pd.DataFrame:
    """Load and validate an APO-SEQ mutation table."""

    if not path.exists():
        raise ApoSeqError(f"Mutation table does not exist: {path}")
    table = pd.read_csv(path, sep="\t")
    validate_mutation_table(table)
    return normalize_mutation_table(table)


def validate_mutation_table(table: pd.DataFrame) -> None:
    """Validate that required mutation-table columns exist."""

    missing = [column for column in REQUIRED_MUTATION_COLUMNS if column not in table.columns]
    if missing:
        raise ApoSeqError(f"Mutation table is missing required columns: {', '.join(missing)}")


def normalize_mutation_table(table: pd.DataFrame) -> pd.DataFrame:
    """Normalize mutation-table numeric and sequence columns."""

    normalized = table.copy()
    normalized["POS"] = pd.to_numeric(normalized["POS"], errors="raise").astype(int)
    normalized["DEPTH"] = pd.to_numeric(normalized["DEPTH"], errors="coerce").fillna(0).astype(int)
    normalized["ALT_COUNT"] = pd.to_numeric(normalized["ALT_COUNT"], errors="coerce").fillna(0).astype(int)
    normalized["AF"] = pd.to_numeric(normalized["AF"], errors="coerce").fillna(0.0).astype(float)
    normalized["REF"] = normalized["REF"].astype(str).str.upper()
    normalized["ALT"] = normalized["ALT"].astype(str).str.upper()
    normalized["MUTATION"] = normalized["REF"] + ">" + normalized["ALT"]
    return normalized


def save_table(table: pd.DataFrame, path: Path) -> None:
    """Write a table as TSV."""

    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)
