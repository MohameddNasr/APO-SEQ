from pathlib import Path

from aposeq.config import load_run_config
from aposeq.slurm import sanitize_job_name, write_sbatch_script


def test_write_sbatch_script_contains_heavy_steps(tmp_path: Path) -> None:
    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "\n".join(
            [
                "name: pwa_delta2_3",
                "assay: pwa",
                "execution: longleaf",
                "input:",
                f"  bam_directory: {tmp_path / 'bam_pass'}",
                "output:",
                f"  directory: {tmp_path / 'Results'}",
                "samples:",
                "  - barcode: barcode01",
                "    genotype: delta2-3 transposase negative",
            ]
        ),
        encoding="utf-8",
    )
    loaded = load_run_config(run_config)

    script = write_sbatch_script(loaded, submit_directory=Path.cwd())

    text = script.path.read_text(encoding="utf-8")
    assert "#SBATCH --job-name=aposeq_pwa_delta2_3" in text
    assert "--steps process_bams coverage call_variants" in text
    assert 'cd "${SLURM_SUBMIT_DIR:-' in text
    assert script.command[0] == "sbatch"


def test_sanitize_job_name() -> None:
    assert sanitize_job_name("delta2-3 transposase negative") == "delta2-3_transposase_negative"
