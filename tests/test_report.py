"""Machine-readable example output (for the BATS pipeline): numpy-safe JSON results."""

import json

import numpy as np

from ti_sim import report


def test_save_results_writes_plain_json(tmp_path):
    out = tmp_path / "r.json"
    res = {"peak": np.float64(0.26), "n": np.int64(3), "v": np.arange(3) * 0.5, "name": "x"}
    assert report.save_results(res, out) == out
    assert json.loads(out.read_text()) == {"peak": 0.26, "n": 3, "v": [0.0, 0.5, 1.0], "name": "x"}


def test_save_results_uses_env_and_is_a_noop_without_it(tmp_path, monkeypatch):
    monkeypatch.delenv("TI_SIM_OUT", raising=False)
    assert report.save_results({"a": 1}) is None
    monkeypatch.setenv("TI_SIM_OUT", str(tmp_path / "sub" / "r.json"))
    p = report.save_results({"a": 1})
    assert p is not None and json.loads(p.read_text()) == {"a": 1}


def test_figure_path_follows_results_dir(tmp_path, monkeypatch):
    monkeypatch.delenv("TI_SIM_OUT", raising=False)
    assert report.figure_path("a.png") == "a.png"
    monkeypatch.setenv("TI_SIM_OUT", str(tmp_path / "res" / "x.json"))
    assert report.figure_path("a.png") == str(tmp_path / "res" / "a.png") and (tmp_path / "res").is_dir()
