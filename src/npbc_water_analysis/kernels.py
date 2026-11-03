"""Explicit loader for the archived production scripts."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

KERNEL_NAMES = {
    "whole": "analyze_radius_equivalence.py",
    "rdf": "compute_rdf_extensions.py",
    "radial": "compute_radial_layer_analysis.py",
    "figures": "make_stochkin_figure_suite.py",
}


def load_kernel(name: str) -> ModuleType:
    if name not in KERNEL_NAMES:
        raise ValueError(f"unknown kernel {name!r}; choose from {sorted(KERNEL_NAMES)}")
    path = Path(__file__).resolve().parents[2] / "analysis_scripts" / KERNEL_NAMES[name]
    spec = importlib.util.spec_from_file_location(f"npbc_kernel_{name}", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load analysis kernel: {path}")
    module = importlib.util.module_from_spec(spec)
    # Dataclasses with postponed annotations consult ``sys.modules`` while the
    # class decorator executes, so register the module before loading it.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module
