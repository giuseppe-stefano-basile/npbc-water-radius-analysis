from pathlib import Path

from npbc_water_analysis.reference import validate_reference_results


def test_packaged_reference_results_are_complete():
    root = Path(__file__).resolve().parents[1] / "reference_results"
    summary = validate_reference_results(root)
    assert summary == {
        "tables": 10,
        "figures": 28,
        "normalized_histograms": 192,
        "rdf_metric_rows": 54,
    }
