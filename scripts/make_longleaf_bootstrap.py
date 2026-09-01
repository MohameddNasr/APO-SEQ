"""Build the one-file Longleaf APO-SEQ bootstrap installer."""

from __future__ import annotations

import base64
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = Path("/tmp/aposeq_longleaf_deploy.tar.gz")
OUTPUT = PROJECT_ROOT / "longleaf_bootstrap_aposeq.sh"
BAM_DIR = "/path/to/GridION/run/bam_pass"


def main() -> None:
    encoded = base64.b64encode(ARCHIVE.read_bytes()).decode("ascii")
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        'APOSEQ_ROOT="${APOSEQ_ROOT:-$HOME/APO-SEQ}"',
        'ARCHIVE="${ARCHIVE:-$HOME/aposeq_longleaf_deploy.tar.gz}"',
        'RUN_CONFIG="config/runs/pwa_delta2_3_transposase_negative.yaml"',
        f'BAM_DIR="${{BAM_DIR:-{BAM_DIR}}}"',
        "",
        'if [ ! -d "$BAM_DIR" ]; then',
        '  echo "ERROR: BAM directory does not exist: $BAM_DIR" >&2',
        "  exit 1",
        "fi",
        "",
        'BAM_COUNT=$(find "$BAM_DIR" -type f -name "*.bam" | wc -l)',
        'if [ "$BAM_COUNT" -eq 0 ]; then',
        '  echo "ERROR: No .bam files were found under: $BAM_DIR" >&2',
        "  exit 1",
        "fi",
        'echo "Found $BAM_COUNT BAM file(s) under $BAM_DIR"',
        "",
        'mkdir -p "$APOSEQ_ROOT"',
        'cat > "$ARCHIVE.b64" <<\'EOF_APOSEQ_ARCHIVE\'',
    ]
    for index in range(0, len(encoded), 76):
        lines.append(encoded[index : index + 76])
    lines.extend(
        [
            "EOF_APOSEQ_ARCHIVE",
            "",
            'base64 -d "$ARCHIVE.b64" > "$ARCHIVE"',
            'tar -xzf "$ARCHIVE" -C "$APOSEQ_ROOT"',
            'cd "$APOSEQ_ROOT"',
            "",
            "python3 -m pip install --user -e .",
            'python3 -m aposeq.cli validate-config --run-config "$RUN_CONFIG"',
            'python3 -m aposeq.cli validate-env --output Results/pwa_delta2_3_transposase_negative/environment_report.json',
            'python3 -m aposeq.cli write-sbatch --run-config "$RUN_CONFIG" --steps process_bams coverage call_variants',
            "sbatch Results/pwa_delta2_3_transposase_negative/slurm/pwa_delta2_3_transposase_negative.sbatch",
            "",
            'echo "APO-SEQ submitted. Monitor with: squeue -u $USER"',
        ]
    )
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    OUTPUT.chmod(0o755)


if __name__ == "__main__":
    main()
