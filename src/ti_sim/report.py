"""Machine-readable results from the examples, so a pipeline (tests/pipeline.bats) can
assert on numbers instead of scraping stdout. save_results writes plain JSON to `path`,
or to $TI_SIM_OUT when no path is given, and does nothing when neither is set."""

import json
import os
from pathlib import Path

import numpy as np


def _plain(x):
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, np.ndarray):
        return _plain(x.tolist())
    if isinstance(x, np.generic):
        return x.item()
    return x


def figure_path(name):
    """Where an example saves a figure: beside $TI_SIM_OUT when set (the cwd may be a
    read-only install, e.g. /opt/ti-sim in the container), else the bare name (cwd)."""
    out = os.environ.get("TI_SIM_OUT")
    if not out:
        return name
    d = Path(out).parent
    d.mkdir(parents=True, exist_ok=True)
    return str(d / name)


def save_results(results, path=None):
    """Write `results` as JSON; returns the Path written or None when there is nowhere to."""
    path = path or os.environ.get("TI_SIM_OUT")
    if not path:
        return None
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_plain(results), indent=2) + "\n")
    return path
