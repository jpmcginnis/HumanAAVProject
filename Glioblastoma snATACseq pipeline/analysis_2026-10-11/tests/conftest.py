"""Pytest configuration. Skips tests gracefully when the Zenodo-deposited data
hasn't been downloaded yet (bare git clone). On the author's machine (or any
clone with data/matrix/ and reports/detail/ populated from Zenodo DOI
10.5281/zenodo.23252324), all 17 tests run normally."""
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent.parent
MATRIX_PRESENT = (HERE / "data" / "matrix" / "enhancer_candidate_matrix.parquet").exists()
REPORTS_PRESENT = (HERE / "reports" / "detail" / "pan_malignant.parquet").exists()

requires_matrix = pytest.mark.skipif(
    not MATRIX_PRESENT,
    reason=("Skipping: data/matrix/enhancer_candidate_matrix.parquet not present. "
            "Download the atlas data from Zenodo DOI 10.5281/zenodo.23252324 and "
            "place it under data/matrix/ to run this test.")
)

requires_reports = pytest.mark.skipif(
    not REPORTS_PRESENT,
    reason=("Skipping: reports/detail/pan_malignant.parquet not present. "
            "Download the atlas reports from Zenodo DOI 10.5281/zenodo.23252324 and "
            "place them under reports/ to run this test.")
)


def pytest_collection_modifyitems(config, items):
    """Auto-apply the skip markers to all tests. Tests that need both markers
    (schema tests that read matrix AND reports) get both — pytest treats the
    strictest active skip as the reason."""
    for item in items:
        # Any test that touches reports/detail/ gets the reports marker.
        if "test_filters" in item.nodeid or "test_schema" in item.nodeid:
            item.add_marker(requires_reports)
        # Any test that reads from data/matrix/ gets the matrix marker.
        if "test_schema" in item.nodeid:
            item.add_marker(requires_matrix)
