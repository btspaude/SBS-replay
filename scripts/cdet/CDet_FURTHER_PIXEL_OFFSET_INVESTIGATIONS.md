# Further CDet pixel-offset investigations

Updated: 2026-10-20.

## Run 5710 offset coverage

Counting literal `n/a` entries in the `run5710_offset_ns` column of
[`CDet_cross_target_pixel_offset_scan.tsv`](CDet_cross_target_pixel_offset_scan.tsv)
gives:

| Coverage | Layer 1: IDs 0–1343 | Layer 2: IDs 1344–2687 | Whole detector |
| --- | ---: | ---: | ---: |
| `n/a` for Run 5710 | 800 | 717 | **1,517** |
| Listed Run 5710 offset | 544 | 627 | **1,171** |
| Total pixels | 1,344 | 1,344 | 2,688 |

Thus **56.44%** of the pixels have `n/a` for Run 5710 in this table.
All 2,688 pixel IDs were checked for uniqueness. Bar IDs are global and
zero-based: `bar = pixel_id / 16`, with integer division.

This is a count of unavailable reference offsets under the scan's rules,
not a count of dead pixels or pixels with no recorded hits. As documented in
the [scan description](CDet_CROSS_TARGET_PIXEL_OFFSET_SCAN.md), the scan retains
nonzero `retained_existing` reference values but treats zero retained rows as
unavailable. These table counts alone do not establish what constants are
active in a particular replay.

Reproduce the count from the CDet script directory:

```python
import csv

with open("CDet_cross_target_pixel_offset_scan.tsv") as source:
    rows = list(csv.DictReader(
        (line for line in source if not line.startswith("#")), delimiter="\t"))
missing = [row for row in rows
           if row["run5710_offset_ns"].strip().lower() == "n/a"]
print(len(missing), "of", len(rows))
```

## Investigation strategy: start with LH2

The current priority is to inspect individual bars on LH2 before switching
back to cross-target runs. The user's observations are that LH2 backgrounds
make uncalibrated timing peaks difficult to identify, while cross-target
illumination leaves many bars poorly covered. These motivate the investigation;
they are not yet a quantified signal/background or coverage result here.

Start with study bar 19 and control bar 30. Inspect whole-bar and 4×4 pixel
distributions of manual dt versus rx without an ellipse centered at −26 ns.
The aim is to find localized structure before choosing a timing center.
Different pixel offsets could smear a whole-bar distribution, making the
per-pixel views especially useful. Localization alone does not establish signal.
Compare candidate structures across pixels/runs and against a control bar.

## Current PulseStudy diagnostics

[`Plot_CDet_PulseStudy.C`](Plot_CDet_PulseStudy.C) reads the run into aligned
vectors once; `plotPulseStudyBars` draws from those vectors.

| Quantity | Source and convention |
| --- | --- |
| Stored rx | `vPulseXResiduals`, from `earm.cdet.pulse.ecal_x_residual`, in metres |
| Manual dt | `vECalAdcTime[event] - vPulseLe[event][hit]`, in ns |
| Stored dt | `vPulseDt`, from `earm.cdet.pulse.ecal_residual`, in ns |

`vPulseLe` receives the 0.01 conversion to ns during loading. Manual dt omits
the additional pulse timing corrections; it is not a subtraction of untouched
hardware counts. Stored rx uses aligned channel geometry and projected ECal x,
not pulse timing. With the current database's spatial settings:

```text
x_corr = 1.07 * channel_x - 0.03 m
rx = x_corr - projected_ECal_x - 0.0 m
```

The settings actually used to replay a ROOT file determine its stored values.

The study/control routine currently produces eight canvases per bar:
summed and 4×4 pixel views of stored dt, raw LE, raw LE versus ToT, and manual
dt versus rx. It additionally draws one **stored-dt-versus-rx** canvas for the
configured control bar, for **17 canvases total**. That extra histogram is
`hPulseBarDt_resVsRx`; the user's `selected == 1` fill follows the control bar
argument/configuration, rather than being permanently fixed to bar 30.

The control histogram has been initialized with the existing 100 rx bins
over `[-xRes, +xRes]` and shared dt binning. Its canvas uses `COLZ` and saves as
`bar030_Control_dtResVsRx.pdf` with the default control bar. Repeated calls
replace the histogram/canvas. The missing `SaveAs` semicolon was also fixed.

Example, from the CDet script directory with `OUT_DIR` set:

```cpp
.L Plot_CDet_PulseStudy.C+
Plot_CDet_PulseStudy("CDet_run6077_pulse_pixel_timing.conf");
plotPulseStudyBars("CDet_run6077_pulse_pixel_timing.conf", "offset_investigation", 19, 30);
```

Outputs go to `PusleStudy/<tag>/`. The current dt display range is −90 to 0 ns;
values outside it remain in histogram flow bins. Check those bins before
concluding that a pixel has no timing structure.

### Selection context

Set `good_pulse.single_cut_radius: -1.0` to disable the ellipse. The current
configuration already does so. This does **not** make PulseStudy uncut:
loading applies ECal energy/time, replay ECal eligibility and spatial flags,
calibration validity, and strict ToT bounds. Plotting additionally applies its
LE/ToT selection and any enabled secondary ECal cuts.

In the two-bar routine, disabling the ellipse activates the standalone rx
cut. Its threshold still reads `analysis.x_difference_max`; the separately
present `analysis.x_res_max` key is unused. This naming mismatch remains to
be resolved. The current manual and stored dt-versus-rx histograms are filled
after the shared plotting selection. Enabling the ellipse therefore selects
them using **stored dt**; they are not independent before/after-cut panels.

## Separate check for pixels with no recorded activity

[`Plot_CDet_PixelActivity.C`](Plot_CDet_PixelActivity.C) and
[`CDet_run6077_pixel_activity.conf`](CDet_run6077_pixel_activity.conf) provide
an independent scan without ECal, timing, ToT, calibration, or spatial cuts.
It counts `earm.cdet.hits.TDCelemID`, including decoder slots with incomplete
pulses, and counts each pixel at most once per event for event occupancy.

The `earm.cdet.hit.tdc_le` branch used by PlotElastic's `vGoodLe` already has
replay hit selection and LE/ToT windows, so it is unsuitable as an entirely
uncut activity check.

```cpp
.L Plot_CDet_PixelActivity.C+
Plot_CDet_PixelActivity("CDet_run6077_pixel_activity.conf");
```

The default scans all available events/segment groups. Results are two layer
maps, `pixel_activity.csv`, and `zero_activity_pixels.txt`, saved under
`CDetPixelActivity/run6077/` by default. Zero recorded activity is sample-specific
and does not prove hardware failure; noise also counts as activity.

## Validation and next checks

ROOT compilation and synthetic-data checks passed for the new control
stored-dt-versus-rx histogram, its separation from manual dt, pixel sums,
redraw cleanup, and changed study/control bars. These are implementation checks,
not evidence that a physics peak has been isolated. Scratch test log:
`/private/tmp/cdet_pulse_dt_rx_residual/check.log`.

Next, inspect the LH2 bar-19/control-bar plots and pixel grids for localization,
and compare the full activity scan with the Run 5710 unavailable-offset list.
Use cross-target runs as a complementary check where they illuminate the
relevant pixels. No new offsets have been extracted or applied in this step.
Explicit before/after-ellipse comparisons and the detector-wide dt-versus-rx
diagnostic remain unimplemented.
