from pathlib import Path

from aposeq.config import load_run_config


def test_ebony_run_config_loads() -> None:
    loaded = load_run_config(Path("config/runs/ebony_cas9_negative.yaml"))
    assert loaded.run["assay"] == "ebony"
    assert loaded.assay["reference"]["length"] == 27383


def test_pwa_run_config_loads() -> None:
    loaded = load_run_config(Path("config/runs/pwa_example.yaml"))
    assert loaded.run["assay"] == "pwa"
    assert loaded.assay["reference"]["length"] == 37389
    assert loaded.assay["dsb"]["positions"] == [13531, 22758]


def test_pwa_delta2_3_negative_run_config_loads() -> None:
    loaded = load_run_config(Path("config/runs/pwa_delta2_3_transposase_negative.yaml"))
    assert loaded.run["assay"] == "pwa_delta2_3_transposase_negative"
    assert loaded.assay["reference"]["length"] == 37389
    assert loaded.assay["dsb"]["positions"] == []
