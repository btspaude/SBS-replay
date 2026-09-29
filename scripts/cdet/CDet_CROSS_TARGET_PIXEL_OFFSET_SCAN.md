# Cross-target pixel-offset scan

`CDet_CrossTargetPixelOffsetScan.C` is a first-pass diagnostic for checking
whether the Run 5710 pixel timing offsets transfer to the other admitted
cross-target runs. It does not activate, edit, or overwrite
`CDet_calibration_dt.dat`.

For each run, the macro runs the existing cross-target master at calibration
stage 0, then calls the hierarchical pixel-fit routine in diagnostic mode.
The resulting candidate is therefore a run-by-run measurement, not a proposed
production calibration. The extractor determines a relative pixel correction
using that run's valid-pixel reference median. A run-wide timing-origin shift
largely cancels in this relative quantity, but differences in illumination,
cuts, geometry, and ECal timing must still be checked before combining runs.

The output table has one row per logical pixel. Its columns are:

```text
pixel_id  run5710_offset_ns  run_<run>_offset_ns  run_<run>_minus_5710_ns ...
```

The difference is defined as

```text
run offset - Run 5710 offset
```

If a pixel has no usable Run 5710 offset, its Run 5710 value is `n/a` and all
corresponding differences are also `n/a`. If a run does not produce a valid
pixel fit, both fields for that run are `n/a`. A nonzero `retained_existing`
value in the Run 5710 results is retained because it can represent the later
Run 5710 half-bar alignment; zero retained rows are treated as unavailable.

Run the default inventory from the CDet script directory:

```cpp
root -l
.L CDet_CrossTargetPixelOffsetScan.C+
CDet_CrossTargetPixelOffsetScan();
```

This uses the cross-target inventory documented in
`CDet_CROSS_TARGET_RUN_INVENTORY.md`, excluding Run 5710. Diagnostic candidate
files, fit-result files, ROOT diagnostics, and plots are placed under
`CDet_cross_target_pixel_offset_scan_work/`. The summary table is
`CDet_cross_target_pixel_offset_scan.tsv`.

The macro also writes compact diagnostic outputs in the work directory so the
scan can be reviewed without reading the wide pixel-by-run table:

```text
CDet_cross_target_pixel_offset_run_summary.tsv
CDet_cross_target_pixel_offset_pixel_summary.tsv
CDet_cross_target_pixel_offset_diagnostics.pdf
CDet_cross_target_pixel_offset_diagnostics.png
CDet_cross_target_pixel_offset_coverage.png
```

The run summary reports the number of valid offsets, the number of valid
comparisons to Run 5710, and the median, RMS, and median absolute deviation
(MAD) of the offset differences for each run. The pixel summary reports the
same robust statistics for each pixel, along with its layer, bar, coverage,
minimum and maximum observed difference, and a simple status field. The PDF
and PNG diagnostics contain an all-comparison difference histogram, a
run-by-run median plot, a per-pixel median plot, and a run-versus-pixel
heatmap. The coverage PNG shows how many runs contributed a valid comparison
for each pixel.

To use a custom list, provide one run per line. A second column can give an
explicit configuration path:

```text
# run  configuration
5722 CDet_run5722_projection.conf
5723 CDet_run5723_projection.conf
```

```cpp
CDet_CrossTargetPixelOffsetScan(
  "my_cross_target_runs.txt",
  "CDet_pixel_timing_fit_results_run5710_reviewed_102_cut_pass.dat",
  "my_pixel_offset_scan.tsv");
```

The Run 5710 reference file is the reviewed result table, not a new fit made
by this macro. Before using any differences to update calibration constants,
inspect the per-run fit counts, fit sources, residual widths, and the
cross-target residual means. The next step should be a common-reference or
hierarchical fit that separates per-pixel offsets from per-run timing shifts.
