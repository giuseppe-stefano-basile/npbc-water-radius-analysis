#!/usr/bin/env python3
"""Print a few scientifically meaningful checks from the reference tables."""

from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[1] / "reference_results/data"

geometry = pd.read_csv(root / "system_geometry.csv")
comparisons = pd.read_csv(root / "npbc_pbc_rdf_comparisons.csv")
radial = pd.read_csv(root / "radial_layer_gOO_metrics.csv")

print("System geometry")
print(geometry.to_string(index=False))
print("\nNPBC--PBC intermolecular RDF errors")
print(comparisons.to_string(index=False))
print("\nOuter-layer radial O--O metrics")
print(radial.query("layer_scheme == 'coarse' and layer_index == 0").to_string(index=False))
