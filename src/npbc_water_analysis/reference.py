"""Validation of the packaged benchmark tables and figures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REQUIRED_TABLES = (
    "temporal_shell_statistics.csv",
    "potential_amplitudes_full_radius.csv",
    "potential_fields_full_radius.csv",
    "regional_distributions.csv",
    "boundary_safe_rdf.csv",
    "intermolecular_rdf_with_block_std.csv",
    "radial_layer_gOO.csv",
    "radial_layer_gOO_metrics.csv",
    "radial_layer_property_distributions_vs_pbc.csv",
    "radial_layer_distribution_metrics_vs_pbc.csv",
)


def validate_reference_results(root: Path) -> dict[str, int]:
    root = Path(root)
    missing = [name for name in REQUIRED_TABLES if not (root / "data" / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing reference tables: {missing}")

    figures = sorted((root / "figures").glob("[0-9][0-9]_*.pdf"))
    if len(figures) != 28:
        raise ValueError(f"expected 28 vector figures, found {len(figures)}")
    numbers = [int(path.name[:2]) for path in figures]
    if numbers != list(range(1, 29)):
        raise ValueError(f"figure numbering is incomplete: {numbers}")

    manifest = pd.read_csv(root / "figure_manifest.csv")
    if len(manifest) != 28 or manifest["figure"].duplicated().any():
        raise ValueError("figure manifest must contain 28 unique rows")

    distributions = pd.read_csv(root / "data/radial_layer_property_distributions_vs_pbc.csv")
    grouped = distributions.groupby(["system", "layer_label", "property"])["probability"].sum()
    if not ((grouped - 1.0).abs() < 1.0e-8).all():
        raise ValueError("one or more radial histograms do not normalize to unity")

    metrics = pd.read_csv(root / "data/radial_layer_gOO_metrics.csv")
    numeric = metrics.select_dtypes(include="number").drop(columns=["wall_hi_A"], errors="ignore")
    if numeric.isna().any().any():
        raise ValueError("radial RDF metric table contains missing numerical values")

    return {
        "tables": len(REQUIRED_TABLES),
        "figures": len(figures),
        "normalized_histograms": len(grouped),
        "rdf_metric_rows": len(metrics),
    }
