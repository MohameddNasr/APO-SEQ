"""Configuration loading and validation for APO-SEQ."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import csv

import yaml

from aposeq.exceptions import ConfigError


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_ROOT = PROJECT_ROOT / "config"


@dataclass(frozen=True)
class LoadedConfig:
    """A validated set of APO-SEQ configuration layers."""

    run: dict[str, Any]
    assay: dict[str, Any]
    defaults: dict[str, Any]
    execution: dict[str, Any] | None
    run_config_path: Path
    assay_config_path: Path


def load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML file as a dictionary."""

    if not path.exists():
        raise ConfigError(f"Configuration file does not exist: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"Configuration file must contain a mapping: {path}")
    return data


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge two dictionaries without mutating inputs."""

    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def resolve_config_path(path_or_name: str | Path, section: str | None = None) -> Path:
    """Resolve either an explicit path or a config name under config/."""

    candidate = Path(path_or_name)
    if candidate.exists():
        return candidate.resolve()

    if section is None:
        section_path = CONFIG_ROOT
    else:
        section_path = CONFIG_ROOT / section

    if candidate.suffix in {".yaml", ".yml"}:
        named = section_path / candidate.name
    else:
        named = section_path / f"{candidate}.yaml"
    return named.resolve()


def load_run_config(run_config: str | Path) -> LoadedConfig:
    """Load and validate a run configuration with its assay and defaults."""

    run_path = resolve_config_path(run_config, "runs")
    defaults = load_yaml(CONFIG_ROOT / "default.yaml")
    run = load_yaml(run_path)

    assay_name = require_string(run, "assay")
    assay_path = resolve_config_path(assay_name, "assays")
    assay = load_yaml(assay_path)

    execution = None
    execution_name = run.get("execution")
    if execution_name:
        execution = load_yaml(resolve_config_path(str(execution_name), "execution"))

    validate_assay_config(assay, assay_path)
    validate_run_config(run, assay, run_path)
    validate_defaults(defaults)
    if execution is not None:
        validate_execution_config(execution)

    return LoadedConfig(
        run=run,
        assay=assay,
        defaults=defaults,
        execution=execution,
        run_config_path=run_path,
        assay_config_path=assay_path,
    )


def require_string(data: dict[str, Any], key: str) -> str:
    value = data.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Missing or invalid required string: {key}")
    return value


def require_version(data: dict[str, Any], key: str) -> str:
    value = data.get(key)
    if isinstance(value, str) and value.strip():
        return value
    if isinstance(value, (int, float)):
        return str(value)
    raise ConfigError(f"Missing or invalid required version: {key}")


def require_int(data: dict[str, Any], key: str, minimum: int | None = None) -> int:
    value = data.get(key)
    if not isinstance(value, int):
        raise ConfigError(f"Missing or invalid required integer: {key}")
    if minimum is not None and value < minimum:
        raise ConfigError(f"{key} must be >= {minimum}; got {value}")
    return value


def validate_assay_config(assay: dict[str, Any], path: Path | None = None) -> None:
    """Validate assay-level biological metadata."""

    name = require_string(assay, "name")
    require_version(assay, "version")

    reference = require_mapping(assay, "reference")
    require_string(reference, "name")
    require_string(reference, "fasta")
    length = require_int(reference, "length", minimum=1)
    require_string(reference, "chromosome")

    regions = assay.get("regions")
    if not isinstance(regions, list) or not regions:
        raise ConfigError(f"Assay {name} must define at least one region")

    seen_region_ids: set[str] = set()
    for index, region in enumerate(regions, start=1):
        if not isinstance(region, dict):
            raise ConfigError(f"Assay region #{index} must be a mapping")
        region_id = require_string(region, "id")
        if region_id in seen_region_ids:
            raise ConfigError(f"Duplicate region id in assay {name}: {region_id}")
        seen_region_ids.add(region_id)
        start = require_int(region, "start", minimum=1)
        end = require_int(region, "end", minimum=1)
        if start > end:
            raise ConfigError(f"Region {region_id} has start > end")
        if end > length:
            raise ConfigError(f"Region {region_id} ends beyond reference length {length}")

    dsb = require_mapping(assay, "dsb")
    positions = dsb.get("positions", [])
    if not isinstance(positions, list):
        raise ConfigError("dsb.positions must be a list")
    for position in positions:
        if not isinstance(position, int) or position < 1 or position > length:
            raise ConfigError(f"Invalid DSB position for assay {name}: {position}")

    _validate_optional_plotting_regions(assay, seen_region_ids)
    _validate_optional_igv(assay)

    if path is not None and path.name.split(".")[0] != name:
        raise ConfigError(f"Assay name {name!r} does not match file name {path.name!r}")


def validate_run_config(run: dict[str, Any], assay: dict[str, Any], path: Path | None = None) -> None:
    """Validate run-level metadata."""

    require_string(run, "name")
    if require_string(run, "assay") != require_string(assay, "name"):
        raise ConfigError("Run assay does not match loaded assay configuration")

    input_data = require_mapping(run, "input")
    require_string(input_data, "bam_directory")
    auto_discover_samples = bool(input_data.get("auto_discover_samples", False))
    barcode_regex = input_data.get("barcode_regex")
    if barcode_regex is not None and not isinstance(barcode_regex, str):
        raise ConfigError("input.barcode_regex must be a string when provided")
    if auto_discover_samples:
        require_string(input_data, "default_genotype")

    output = require_mapping(run, "output")
    require_string(output, "directory")

    samples = run.get("samples")
    if auto_discover_samples and samples is None:
        samples = []
    if not auto_discover_samples and (not isinstance(samples, list) or not samples):
        raise ConfigError("Run config must define at least one sample")
    if samples is not None and not isinstance(samples, list):
        raise ConfigError("Run samples must be a list")
    seen_barcodes: set[str] = set()
    for index, sample in enumerate(samples or [], start=1):
        if not isinstance(sample, dict):
            raise ConfigError(f"Sample #{index} must be a mapping")
        barcode = require_string(sample, "barcode")
        if barcode in seen_barcodes:
            raise ConfigError(f"Duplicate barcode in run config: {barcode}")
        seen_barcodes.add(barcode)
        require_string(sample, "genotype")

    overrides = run.get("overrides", {})
    if overrides is not None and not isinstance(overrides, dict):
        raise ConfigError("run.overrides must be a mapping")


def validate_defaults(defaults: dict[str, Any]) -> None:
    """Validate default pipeline settings."""

    analysis = require_mapping(defaults, "analysis")
    af_thresholds = analysis.get("af_thresholds")
    if not isinstance(af_thresholds, list) or not af_thresholds:
        raise ConfigError("analysis.af_thresholds must be a non-empty list")
    for threshold in af_thresholds:
        if not isinstance(threshold, (int, float)) or threshold < 0 or threshold > 1:
            raise ConfigError(f"Invalid AF threshold: {threshold}")

    alt_counts = analysis.get("alt_count_thresholds")
    if not isinstance(alt_counts, list) or not alt_counts:
        raise ConfigError("analysis.alt_count_thresholds must be a non-empty list")
    for threshold in alt_counts:
        if not isinstance(threshold, int) or threshold < 0:
            raise ConfigError(f"Invalid ALT_COUNT threshold: {threshold}")


def validate_execution_config(execution: dict[str, Any]) -> None:
    """Validate execution backend settings."""

    require_string(execution, "name")
    backend = require_string(execution, "backend")
    if backend not in {"local", "slurm"}:
        raise ConfigError(f"Unsupported execution backend: {backend}")


def resolve_run_samples(run: dict[str, Any], output_directory: Path | None = None) -> list[dict[str, object]]:
    """Return configured samples or recover auto-discovered samples from outputs."""

    samples = run.get("samples")
    if isinstance(samples, list) and samples:
        return samples

    input_data = run.get("input", {})
    if not isinstance(input_data, dict) or not input_data.get("auto_discover_samples"):
        return []

    genotype = str(input_data.get("default_genotype", "unknown"))
    if output_directory is not None:
        manifest = output_directory / "BAM" / "bam_manifest.tsv"
        if manifest.exists():
            return _samples_from_bam_manifest(manifest, fallback_genotype=genotype)
        sorted_dir = output_directory / "BAM" / "sorted"
        if sorted_dir.exists():
            barcodes = sorted(
                path.name.removesuffix(".sorted.bam")
                for path in sorted_dir.glob("*.sorted.bam")
            )
            return [
                {"barcode": barcode, "genotype": genotype, "replicate": index}
                for index, barcode in enumerate(barcodes, start=1)
            ]
    return []


def _samples_from_bam_manifest(path: Path, fallback_genotype: str) -> list[dict[str, object]]:
    samples: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            barcode = row.get("barcode", "")
            if not barcode or barcode in samples:
                continue
            replicate = row.get("replicate") or len(samples) + 1
            samples[barcode] = {
                "barcode": barcode,
                "genotype": row.get("genotype") or fallback_genotype,
                "replicate": replicate,
            }
    return [samples[barcode] for barcode in sorted(samples)]


def require_mapping(data: dict[str, Any], key: str) -> dict[str, Any]:
    value = data.get(key)
    if not isinstance(value, dict):
        raise ConfigError(f"Missing or invalid required mapping: {key}")
    return value


def _validate_optional_plotting_regions(assay: dict[str, Any], region_ids: set[str]) -> None:
    plotting = assay.get("plotting", {})
    if not isinstance(plotting, dict):
        raise ConfigError("plotting must be a mapping")
    shaded = plotting.get("default_shaded_regions", [])
    if not isinstance(shaded, list):
        raise ConfigError("plotting.default_shaded_regions must be a list")
    unknown = set(shaded) - region_ids
    if unknown:
        raise ConfigError(f"Unknown shaded region id(s): {sorted(unknown)}")


def _validate_optional_igv(assay: dict[str, Any]) -> None:
    igv = assay.get("igv", {})
    if not isinstance(igv, dict):
        raise ConfigError("igv must be a mapping")
    flank_size = igv.get("flank_size")
    if flank_size is not None and (not isinstance(flank_size, int) or flank_size < 0):
        raise ConfigError("igv.flank_size must be a non-negative integer")
