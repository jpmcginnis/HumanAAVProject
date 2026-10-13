# Regeneration log — v2026-10-11 (built 2026-10-11, rc2 revision)

This file is the receipt that the deliverables and tests in this repository
were regenerated end-to-end before being pushed. A reviewer who clones does
not have to rerun the pipeline to trust the numbers — this is what the author
saw when they ran the regeneration themselves.

## Build sequence

```bash
# Download the atlas data from Zenodo first:
python3 scripts/zenodo_upload.py --download 10.5281/zenodo.23252324 --out .

# Then regenerate the deliverables:
python3 scripts/build_reports.py            # 16 report parquets → annotated CSVs + per-report HTMLs + top30 BED
python3 scripts/build_polished_reports.py   # 12 top-level polished top-100 CSVs + combined xlsx/csv/parquet + top10_per_arm_combined
python3 scripts/build_combined_html.py      # sortable + filterable + hover-definition HTML
python3 -m pytest tests/ -v                 # 17/17 regression tests
```

## pytest output on the author's full tree (2026-10-11 21:14)

```
============================= test session starts ==============================
platform darwin -- Python 3.13.5, pytest-9.1.1, pluggy-1.5.0 -- /usr/local/Caskroom/miniconda/base/bin/python3
cachedir: .pytest_cache
rootdir: /Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results
plugins: zarr-3.4.0, anyio-4.13.0
collecting ... collected 17 items

tests/test_computations.py::test_distance_category_boundaries PASSED     [  5%]
tests/test_computations.py::test_nearest_tss_picks_closest PASSED        [ 11%]
tests/test_computations.py::test_nearest_tss_handles_missing_chrom PASSED [ 17%]
tests/test_computations.py::test_daigle_z_formula PASSED                 [ 23%]
tests/test_filters.py::test_patient_level_replication PASSED             [ 29%]
tests/test_filters.py::test_distal_top30_bed_excludes_promoters PASSED   [ 35%]
tests/test_filters.py::test_top50_combined_has_800_rows PASSED           [ 41%]
tests/test_filters.py::test_annotated_distance_categories_consistent PASSED [ 47%]
tests/test_schema.py::test_peak_provenance_columns PASSED                [ 52%]
tests/test_schema.py::test_peak_provenance_sources_match_canonical_slugs PASSED [ 58%]
tests/test_schema.py::test_enhancer_matrix_has_peak_id PASSED            [ 64%]
tests/test_schema.py::test_enhancer_matrix_has_pseudobulk_columns PASSED [ 70%]
tests/test_schema.py::test_pan_malignant_report_schema PASSED            [ 76%]
tests/test_schema.py::test_annotated_csvs_have_gene_cols PASSED          [ 82%]
tests/test_schema.py::test_n_target_patients_within_atlas_bounds PASSED  [ 88%]
tests/test_schema.py::test_greenwald_fraction_in_range PASSED            [ 94%]
tests/test_schema.py::test_union_bed_well_formed PASSED                  [100%]

============================= 17 passed in 10.37s ==============================
```

## Expected pytest output on a bare git clone (Zenodo data not yet downloaded)

A bare clone of this repo does not contain the data/matrix/ or reports/detail/
parquets — they live at the Zenodo deposit. Running pytest on a bare clone
should show **4 passed, 13 skipped**, with each skip reason naming the Zenodo
DOI to download. Once the Zenodo data is in place, all 17 tests run.
