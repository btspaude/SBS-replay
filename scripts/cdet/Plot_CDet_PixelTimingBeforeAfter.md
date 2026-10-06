# Cross-target pixel timing before and after calibration

`Plot_CDet_PixelTimingBeforeAfter.C` makes eight canvases for DNP slides:

- Six bar canvases: LE, LE versus ToT, and ECal-minus-CDet Δt, each before and
  after calibration. Each shows **four columns and two rows**, retaining the
  middle two rows of the original 4×4 layout.
- Two full-detector canvases: CDet time versus ECal time, one with pixel timing
  offsets applied and one with the full timing calibration applied.

Any logical pixel ID selects its containing 16-pixel bar: `bar = pixel / 16`.
The bar canvases display pixel IDs `16*bar + 4` through `16*bar + 11`. The
default selected pixel **469** selects bar 29 and displays **468–475**; its
aligned bar base is 464. Layer 2 pixel IDs
are also supported. Bar numbering is global, 0–167. The selected bar does not
restrict the full-detector correlation or the detector-wide DT reference fits.

Saved per-pixel LE–ToT polygons now select **both** timing states by default.
The default file is the same reviewed 102-cut Run 5710 artifact used by
`Run_CDet_Calibration_Run5710_Accepted.C`; see **Saved polygon selection** below
for its coordinates, provenance limits, and how to choose another file.

## Run it

### Replay pairs and the Bar 30 pulse spectrum

The same file also provides an independent entry point for two replay-timing
plots and two comparisons with the optional y-propagation correction. The
comparisons are enabled by default with `n = 1.59`. This entry point does not
run the eight legacy before/after canvases:

```cpp
.L Plot_CDet_PixelTimingBeforeAfter.C+
Plot_CDet_PairDTvsDXAndBarTiming(
    "CDet_run6077_projection.conf", 30,
    "cdet_pair_timing", nullptr, -1);
```

`nullptr` uses `OUT_DIR`; replace it with the replay-file directory if needed.
The `-1` event override processes all entries in the configured segment range,
including rollover files. The default `-2` uses `analysis.events` instead.
Progress prints and flushes every 1,000 entries, showing the run, current
entry, requested entry limit, and total chain entries. This counts every
entry read, including events rejected by the selections.
Bar IDs are zero-based: **Bar 30 comprises logical pixels 480–495**, including
all instrumented pixels in that block, not just the middle eight shown by the
legacy before/after canvases.

The first canvas uses **one entry per stored `earm.cdet.pair.*` pair**, across
the full detector. It reproduces the original `hDtvsDxCDetECal` coordinates:

```text
dt = pair.ecal_residual = t_ECal - (t_L1,corr + t_L2,corr)/2
dx = [(x_L1,corr - x_ECal,projected,L1)
    + (x_L2,corr - x_ECal,projected,L2)]/2
```

Member positions come from `pulse.x_corr` and `pulse.ecal_x_proj`, indexed by
`pair.pulse_index_l1/l2`. This dx is the **pair-mean position residual**, not
`pair.dx` (the inter-layer displacement) or `pair.trajectory_residual` (the
inter-layer displacement minus the ECal prediction). The routine uses the
stored assignment as-is; any upstream pair cuts remain part of that sample.
It applies the configured inclusive ECal energy interval and no new pair cut.

The second canvas exactly reproduces the population of
`hCDetBar30ECalMinusCDet_ProjectedQuality` in
`Plot_CDet_GoodPulseCandidates_AllTDC.C`: the inclusive
`analysis.ecal_energy_min/max` interval, finite pulse ID/residual/ToT, and all
four stored flags `calib_valid`, `ecal_eligible`, `spatial_pass`, and
`broad_quality_pass`. Each qualifying pulse in the selected bar contributes
its `pulse.ecal_residual = t_ECal - t_CDet,corr`. **No stored pair or second-layer
hit is required.** Multiple pulses in one event contribute separately; the
canvas reports both pulse and event counts.

Here “ECal projection” means the replay's `spatial_pass` flag, matching that
existing pulse diagnostic. Its spatial/time/quality limits are already encoded
in the replay flags. No additional configuration time window or saved pixel
polygon is applied, and existing replay timing corrections are not reapplied. In
particular, replay broad quality is not evidence of a saved polygon cut.

The two additional canvases compare these exact same pairs/pulses before and
after the y correction from `plotCDetLayersTimeComp()` in
`PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C`. Pairing
and pulse selection are completed using the stored flags and membership;
the correction changes only the plotted timing values:

```text
b_side = -n/c for the left side, +n/c for the right side
delta_i = b_side * (pulse.ecal_y_proj[i] - pulse.y[i])
dt_pulse,y = pulse.ecal_residual + delta_i
dt_pair,y = pair.ecal_residual + (delta_L1 + delta_L2)/2
```

`c = 0.299792458 m/ns`, so the default magnitude is `5.30367 ns/m`.
The positive sign in the residual correction follows from
`t_CDet,y = t_CDet,corr - delta_i` and `dt = t_ECal - t_CDet`.
Each pair member uses its own layer projection and readout-side sign,
including pairs crossing the central seam. Channel IDs modulo 1344 below
672 use the left-side sign; the remaining IDs use the right-side sign.
For Bar 30, a projected y displacement of +0.10 m changes dt by −0.530 ns.

`pulse.y` stores the fixed half-bar center; the current database gives the
same center to every instrumented paddle within each half-bar. The correction
uses `pulse.ecal_y_proj` directly, not `pulse.ecal_y_residual`, which also
contains the spatial-selection offset. This reproduces the master's center
and projection convention without loading its geometry or calibration tables.
The two geometry branches must be present and finite for selected entries;
otherwise the enabled correction stops with a diagnostic instead of producing
a partially corrected comparison.

The pair comparison shows the two 2D distributions side by side with matching
color scales. The bar comparison overlays the spectra and separately fits
both with the same Gaussian plus linear-background model and fit interval.
This is a model-based timing diagnostic; it does not refit a propagation
speed, rerun selection, modify the ROOT input, or install new calibration.

Following the first five arguments, optional arguments are:

```cpp
// dx bin width/min/max (m), dt bin width/min/max (ns),
// Bar-spectrum fit min/max (ns), savePlots, y-correction refractive index
Plot_CDet_PairDTvsDXAndBarTiming(
    "CDet_run6077_projection.conf", 30, "cdet_pair_timing", nullptr, -1,
    0.002, -0.16, 0.16, 1.0, -60.0, 30.0, -55.0, -10.0, true, 1.59);
```

Set the final argument to `0.0` to disable the y comparison and produce only
the original two replay-timing canvases; the y geometry branches are then
not required. This routine uses that argument, not the master's separate
`display.y_correction_refractive_index` configuration key.

Ranges control the displayed histograms, not selection. The bar spectrum uses
the existing diagnostic's Gaussian plus linear-background fit and default
fit/peak-seed windows. The fit is for display and never updates calibration.
All canvases remain open and save PDF/PNG by default:

```text
cdet_pair_timing/CDet_run6077_pair_ecal_cdet_dt_vs_dx.{pdf,png}
cdet_pair_timing/CDet_run6077_bar030_ecal_cdet_dt_projected_quality.{pdf,png}
cdet_pair_timing/CDet_run6077_pair_ecal_cdet_dt_vs_dx_y_comparison.{pdf,png}
cdet_pair_timing/CDet_run6077_bar030_ecal_cdet_dt_y_comparison.{pdf,png}
```

Validated on 2026-10-06 with ROOT ACLiC and all 354,656 entries in the local
Run 6077 replay set. The configured ECal energy interval admitted 201,698
events; the plots contain 113,689 stored pairs in 58,477 events and 9,470
Bar-30 pulses in 6,656 events. Of those bar pulses, 2,965 come from events
with no stored pair. Both histograms matched independent `TTree::Draw`
selections bin-for-bin, including underflow/overflow. Both PDF/PNG pairs
exported and the PNG layouts were inspected. The existing Podd dictionary
autoload warnings appeared while opening the files; compilation, reading,
and histogram comparisons completed successfully. These results reproduce
the calibration and flags already stored in that local replay dataset.

The y-comparison revision passed ACLiC compilation and synthetic checks of
left/right signs, zero displacement, opposite-side pair members, unpaired
bar pulses, unchanged selections, and disabling the correction. On the same
full Run 6077 dataset, all four histograms matched independent `TTree::Draw`
formulas bin-for-bin, including flow bins. The two original histograms also
matched the previously verified release. Each comparison retained all 113,689
pairs or 9,470 bar pulses, respectively. With `n = 1.59`, the Bar-30 fitted
peak width changed from `2.91 +/- 0.10 ns` to `2.82 +/- 0.10 ns`; the fitted
centroid changed from `-27.78 +/- 0.08 ns` to `-27.54 +/- 0.08 ns`. These are
fits to the same events under an assumed propagation model, not a measurement
of the propagation speed or a claim of a statistically established improvement.
All four PDF/PNG pairs exported, and the comparison PNG layouts were inspected.

### Legacy hit before/after canvases

From `gocdetscripts` / `SBS-replay/scripts/cdet`, start a fresh ROOT session:

```cpp
.L Plot_CDet_PixelTimingBeforeAfter.C+
Plot_CDet_PixelTimingBeforeAfter("CDet_run5710_projection.conf", 469);
```

The existing run configuration supplies event limits, segment limits and
analysis cuts. Input defaults to `OUT_DIR`. To override the destination,
input directory and event limit:

```cpp
Plot_CDet_PixelTimingBeforeAfter(
    "CDet_run5710_projection.conf", 469,
    "dnp_run5710_bar29", "/path/to/Rootfiles", 100000);
```

The event override defaults to `-2` (use configuration); `-1` reads all selected
entries. As in the master, `analysis.first_event` is not used. This small macro
handles one run per call; it fits peaks for display but does not combine runs
or update calibration constants.
`analysis.calibration_stage` is intentionally ignored: every call produces both
timing states. The bar plots' before state is uncorrected; the detector
correlation's before state has pixel offsets applied, as labeled on its canvas.

The next arguments control display binning, not acceptance. Four optional
arguments set the before/after DT fit windows; the final optional argument
selects the polygon file:

```cpp
// LE bin width/min/max, ToT bin width/min/max, savePlots,
// dt bin width/min/max, ECal-time bin width/min/max,
// before DT fit min/max, after DT fit min/max
Plot_CDet_PixelTimingBeforeAfter(
    "CDet_run5710_projection.conf", 469, "dnp_run5710_bar29", nullptr, -2,
    1.0, 0, 60, 1.0, 0, 40, true,
    1.0, -40, 10, 1.0, 5, 40,
    -30, 10, -30, 10);
```

The dt display defaults to −40 to +10 ns. The correlation plots use ECal time
on the horizontal axis (5 to 40 ns by default) and CDet LE on the vertical
axis (the same `leMin/leMax` as the LE plots). All timing axes default to 1 ns
bins, following the 2026-10-02 ifarm update. Existing shorter calls remain valid.
The DT fits use the chosen DT binning; the 1 ns default now matches the
production extractor's default bin width. Both fit windows default to −30 to
+10 ns, matching the current cross-target extractor. Set them to the appropriate
signal region for other timing conditions; for the local Run 6077 check, both
were −45 to −10 ns.
Fit-window limits must lie within the DT histogram; either endpoint may equal
the corresponding histogram edge. ROOT maps a coordinate equal to the upper
edge to its overflow bin, so the macro clamps the fit's counting and seed-bin
indices to the visible bins. Out-of-range fit windows report the offending
before/after interval and histogram limits explicitly. Changing a fit window
changes the fit and reference estimate, not the selected hits.

These are display ranges, not event cuts: for example, the ECal-time display
does not replace `analysis.ecal_time_min/max`. For another run, choose a display
range that includes its selected ECal times.

The eight canvases remain open in interactive ROOT. Each saves as both PDF and
PNG by default; `savePlots=false` displays them without saving. Default outputs:

```text
cdet_timing_before_after/CDet_run5710_bar029_le_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar029_le_after.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar029_le_vs_tot_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar029_le_vs_tot_after.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar029_ecal_cdet_dt_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar029_ecal_cdet_dt_after.{pdf,png}
cdet_timing_before_after/CDet_run5710_detector_cdet_t_vs_ecal_t_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_detector_cdet_t_vs_ecal_t_after.{pdf,png}
```

Reusing a destination for the same run/bar replaces those figures. The two
detector correlation files are shared across selected bars. Earlier
`barNNN_cdet_t_vs_ecal_t_*` files in an existing destination are left untouched;
use the new `detector_cdet_t_vs_ecal_t_*` files for the aggregate views.
The PDF is the vector export for slides. Before/after pairs use the same
binning and x ranges. **LE and DT panels have independent count-axis maxima**,
each 5% above its own tallest bin (or its displayed DT fit, if taller).
This keeps the broader before spectra from inheriting a much taller after
peak's scale. Their axes remain counts, not normalized distributions, and
their minima remain zero. The 2D plots retain matched color scales.
The compact `N` and `SD` labels use the entire
selected sample, including values outside the display range. No underflow or
overflow boxes are drawn. LE retains those two labels. The DT panels also
show the fitted signal centroid `mu_i`, Gaussian width `sigma_fit`, their
fit uncertainties in ns, and dimensionless chi-square/NDF. SD describes the
whole distribution, including background; `sigma_fit` describes the fitted
signal peak. Neither is automatically an intrinsic CDet timing resolution.
LE statistics and DT fit labels sit in compact headers above the plot frames;
the count axes no longer reserve extra height for those labels.

The dt convention follows the cross-target pixel-calibration plots:
`dt = earm.ecal.adctime - CDet LE`, evaluated with the before/after CDet time.
It is the ECal-minus-single-hit time difference, not the inter-layer pair dt.
No additional dt cut is applied.

## DT fits and the run reference

The red DT curves use `gaus(0)+pol1(3)`: a Gaussian signal plus a linear
background, following the **automatic hierarchical** path in
`extractHierarchicalCDetPixelTimingOffsetsDiagnostic()`. This is the current
method described in [the cross-target calibration guide](CDet_CROSSTARGET_TIMING_CALIBRATION.md),
rather than the older inverse-variance-weighted reference method.

For each before/after sample separately:

1. Sum counts from the instrumented pixels in each contiguous eight-pixel
   group. Fit that sum within the configured broad window if it has at least
   100 entries in the window.
2. For a pixel with at least 35 entries in the broad window, use its valid
   group mean to seed a fit whose half-width is `clamp(2*sigma_group, 4, 8)` ns,
   clipped to the broad window. If the group fit is invalid, try the broad
   individual fit instead.
3. Apply the master's individual-fit checks: successful fit, positive finite
   centroid uncertainty no larger than 2 ns, Gaussian width 0.5–8 ns,
   centroid more than half a bin from the fit boundaries, positive amplitude,
   amplitude/error at least 1.5, positive NDF and finite chi-square/NDF no
   greater than 15. The group fit has the master's looser acceptance; its
   chi-square does not veto an individual fit.
4. Define `mu_0` as the median of **all valid individual pixel centroids across
   the selected detector sample**, not just the displayed bar. If no individual
   fits qualify, use the median of valid group means, explicitly labeled as
   group fits. If neither qualifies, show `mu_0 = n/a`.

The before and after DT canvas headers each show their own `mu_0`, the number
and type of fits defining it, and the broad fit window. Low-statistics or
rejected individual fits display `mu_i, sigma_fit = n/a`; when available,
their group mean is labeled separately. A group fallback is not presented
as an independently measured pixel centroid. Known unused pixels are not fit.

These are fresh **plot-sample** references for this run, timing stage, selected
segments, event limit, configured layers and cuts. They are not a saved
production-calibration `mu_0`, and they are not `shift_ns`. The run shift centers
the corrected CDet time; `mu_0` here summarizes `t_ECal - t_CDet`. The two need
not agree, and the before/after references need not agree either.

The peak model, group seeding, acceptance checks and reference rule match the
automatic production fitting method. The plot sample remains the requested
master good-hit sample with the common polygon mask described below; no extra
extractor-only 1–12 GeV energy cut is added. Consequently, these plots are not
a reproduction of the commissioned Run 5710 fit constants, whose calibration passes used
reviewed polygons and stage-specific inputs. Only the loaded constants are
applied; the fitted centroids never change a correction or create a `.dat` file.

All plots use the shared `CDetPlotStyle.h` fonts and axes. The full-detector
CDet-time versus ECal-time canvases follow the master's `cCDetTvsECalT`
convention: `COLZ` with black, style-20 `ProfileX` mean markers (size 0.6).
Each accepted CDet hit contributes one point at its event's ECal time. All
accepted pixels in the configured layer selection contribute, independently
of the bar chosen for the other plots. The profiles use only the displayed
CDet-LE bins; this pooled trend is a visual summary, not a new ECal calibration
fit. Empty displayed pixels are left empty and unused pixels are labeled.

## Event and hit selection

The macro has its own short event loop. It includes the cross-target master to
reuse configuration/calibration readers, geometry constants, and bad/unused
pixel lists, but does not invoke the master's main analysis or plotting suite.
Use a fresh ROOT process rather than loading another copy of the master first.

To reproduce the master's `vGoodLe`, this macro reads its legacy
`earm.cdet.hit.*` collection and `earm.cdet.tdc_mult`, along with ECal quantities.
Raw `earm.cdet.hits.*` values provide the optional reference-channel time. It
does not substitute replay `pulse.*` or `pair.*` selections for this sample.

The same selections are evaluated before calibration:

- ECal reconstruction: `-1.5 < x < 1.5 m`, `-1.2 < y < 1.2 m`, and both
  coordinates nonzero; strict ECal time and energy limits from `analysis.*`.
- Inclusive uncorrected LE and ToT limits from `analysis.*`.
- `tdc_mult[hit] < 100`, using the same indexing as the master.
- Corrected geometry `x = 1.08*xhit - 0.03 m`, with `x < 998 m`.
- Projection using `zECal = 6.144 m`: inclusive configured x-residual limit,
  and `|yhit - yECal*zhit/zECal - analysis.y_offset| <= 0.36 m`.
- Count qualifying hits after optional bad-channel suppression, require the
  configured occupancy in each layer, then apply `analysis.layer_choice`.
  Choice 3 requires at least one qualifying hit in each layer.

The master currently applies a valid `[ECalSelection]` window in the run timing
file over the configured window; this macro follows that actual behavior.
No additional elastic selection, layer pairing, or ellipse selection is
imposed. Saved LE–ToT polygons filter hits after these base event/occupancy
selections; the event occupancies are not recalculated after that hit mask.
The bar and middle-eight-pixel display restrictions
apply only to the LE, LE–ToT and DT canvases, after evaluating the full
detector's event-occupancy cuts. The CDet/ECal correlations and the temporary
DT histograms defining the detector-wide references use all accepted pixels.

## Saved polygon selection

The final `pixelCutFile` argument defaults to:

```text
CDet_run5710_halfbar_aligned_final_archive/CDet_pixel_quality_cuts_run5710_halfbar_aligned_final.root
```

This file contains 102 reviewed Run 5710 cuts, including cuts for displayed
pixels 468 and 471. The macro reads and clones each
`pixel_NNNN/cut_le_vs_tot` without modifying the file. A missing/unreadable file
or a file containing no pixel cuts stops before plotting. When a cut contains
`source_run` metadata, a different input run also stops: Run 5710 cuts are not
silently transferred to LH2. Supply a matching file for another run.

Each saved polygon is evaluated once as:

```text
cut.IsInside(ToT_ns, raw_LE_ns - reference_ns + pixel_offset_ns)
```

That LE coordinate follows the actual `editCDetPixelLeTotCut()` drawing code:
its `vPaddleGoodLe` values include reference subtraction and pixel offsets,
but are not updated by ECal, time-walk, or run-shift corrections. The archived
file records `calibration_stage = 7`; that metadata describes the analysis
invocation and does not make the editor's per-pixel vector fully corrected.
The same editor behavior is present in the archive's source revision
`3dedf214`.

The existing hierarchical extractor instead tests its polygons against
`vGoodLe`, which is fully corrected at stage 7. This is a pre-existing
editor/extractor coordinate discrepancy. This plot macro follows the drawing
coordinates; it does not change the extractor or refit any calibration.
Accordingly, the selected population is not asserted to reproduce the
historical polygon-gated extraction. The active master also differs from the
archive-era offsets, so physical acceptance must be reviewed on Run 5710
before interpreting the revised plots as a validated scientific result.

For a pixel with a saved polygon, an outside hit is removed from **all**
before/after LE, LE–ToT, DT, and detector-correlation histograms. Both displayed
times are then filled from the same surviving hit. Pixels without a polygon
keep the base selections. The polygons supplement the configured LE/ToT cuts;
they do not replace those cuts. Detector-wide reference fits use the same
filtered sample. Plot headers identify polygon selection, and the terminal
reports the file, number of loaded polygons, tested hits, and rejected hits
for the detector and selected bar.

To use freshly reviewed cuts, append their filename after the fit windows:

```cpp
Plot_CDet_PixelTimingBeforeAfter(
    "CDet_run5710_projection.conf", 469, "dnp_run5710_bar29", nullptr, -2,
    1, 0, 60, 1, 0, 40, true,
    1, -40, 10, 1, 5, 40, -30, 10, -30, 10,
    "CDet_pixel_quality_cuts.root");
```

An explicit final `""` or `nullptr` disables polygons for an uncut comparison.
Such plots and their terminal output are labeled **Pixel polygons disabled**.
Changing the cut file changes the sample; use a separate output directory
when retaining both versions.

## Timing correction and plotting-only scope

Within each before/after comparison, both views contain identical events and
hits. For the **bar LE, LE–ToT and DT plots**, before means the stored good-hit
LE multiplied by 0.01 to obtain ns, with no macro timing correction applied.
For the **full-detector CDet/ECal correlation only**, before includes pixel
offsets and the configured reference treatment, but no ECal correction,
time-walk correction or run `shift_ns`:

```text
t_before_bar         = 0.01 * hit.tdc_le
ToT                  = 0.01 * hit.tdc_tot
t_before_correlation = t_before_bar - reference + pixel_offset
t_after              = t_before_correlation - (p0 + p1*tECal) + delta
                       - walk_p1_layer*(1/sqrt(ToT) - 1/sqrt(ToTref_layer))
                       + shift_ns
```

Reference subtraction defaults to disabled. When enabled in the configuration,
the last valid raw channel-2696 reference is used, matching the master; no valid
reference leaves a zero subtraction. It affects the correlation's before
state and all after states, while the bar before states remain raw. ToT is
unchanged. The time-walk fit-domain
limits are metadata and do not cut the application of that correction. Hits
are not reselected after correction.

The macro reads `CDet_calibration_dt.dat` and `CDet_run<run>.dat` from the current
working directory, including any run-specific ECal `p0/p1` overrides. It requires
complete pixel, ECal (including fixed delta), time-walk and run-shift constants
before labeling a view fully corrected. It never updates those files.

The macro is solely a plot producer. It uses the master's base hit selection
and timing formula, adds the common polygon mask, fills the histograms
immediately, and retains no run-wide event vectors. Per-event candidate indices are local and
discarded after each event. Only the canvases and optional PDF/PNG exports
remain; no event tree, analysis-vector interface, or calibration output is
created. Shared calibration readers still load the existing constants into
the master's calibration variables.

Two differences from the old plotting shortcut are deliberate: the old
`plotPaddles()` reads `vPaddleGoodLe`, which is filled with pixel-offset-corrected
times but is not updated by the later ECal/time-walk/global corrections. Also,
the master's stage 0 still adds the run's global shift. Neither is suitable as
the respective fully corrected or completely uncorrected view requested here.

## Validation history (before polygon selection)

The checks in this section describe revisions without polygon filtering.
They remain the timing-formula and base-selection parity record.

The initial four-canvas version passed ROOT ACLiC compilation. On the first
20,000 entries of the locally
available Run 6077 dataset, the new routine and the original stage-7 master
selected identical input entries, hit IDs and ToTs: 4,685 accepted events and
41,627 hits. The largest LE difference was `1.42e-14 ns` (floating-point
roundoff). Before times also matched stage 0 after removing its run shift.
Those four canvases had 16 pads, identical before/after entry counts and matched
per-pixel scales; PNG exports were visually checked. These are software parity
checks on the available LH2 sample, not cross-target physics results. Run 5710
ROOT input was not available locally, so its slide plots still need to be run
on the user's cross-target dataset. ROOT emitted pre-existing analyzer-library
autoload warnings during the check, but tree reading and all comparisons
completed successfully.

The eight-canvas plotting-only revision also compiled successfully and was
checked on the same 20,000 entries. All 128 pixel histograms matched reference
histograms filled from the original master's times **bin-for-bin**, including
underflow/overflow. This checked the new dt sign, the ECal/CDet axis order,
identical before/after populations, and matched display scales. All eight PDFs
and eight PNGs were produced; the added plot layouts were visually reviewed.

The DT-fit revision passed ACLiC compilation and repeated the 128-histogram
parity check. With the same Run 6077 sample, 0.5 ns bins, and a −45 to −10 ns
fit window, it matched the original automatic hierarchical extractor on all
5,376 before/after pixel-fit acceptance classifications. All 74 accepted
individual fits agreed in centroid, centroid uncertainty and width to within
the reference text file's precision (0.0001 ns tolerance). Both detector
medians agreed: 31 valid individual fits before and 43 after. An empty sample
produced no fits and no finite reference. These checks used temporary
diagnostic outputs and did not activate or overwrite calibration constants.
Saved DT layouts were checked for Bar 30's low-statistics labels and Bar 27's
accepted fits; the compact statistics and header reference remain clear of
the histogram peaks.

On 2026-10-02, the ifarm update changed the timing bins to 1 ns, the DT display
to −40 to +10 ns, and the ECal-time display to 5–40 ns. The old validation
incorrectly rejected the default fit endpoint of +10 ns when it equaled the
histogram maximum. The fix allows matching endpoints and excludes overflow
from the fit-entry count, peak search and background seed. A regression check
failed before this fix and passed afterward: six fitted pixel signals and
their detector reference stayed unchanged after adding one million counts
to each underflow and overflow bin. A 20,000-event Run 6077 check produced all
128 histograms with the updated binning/ranges and identical before/after
populations. Invalid windows outside the histogram still fail with explicit
range messages. A forced ACLiC rebuild also passed without the eight compiler
warnings previously inherited from the master; index types and intentionally
unused legacy arguments were cleaned up without changing selections.

The subsequent middle-two-row and full-detector-correlation revision passed
ACLiC compilation and a check on the same 20,000 Run 6077 entries. All 50
displayed histograms (48 pixel spectra and two detector correlations) matched
the original master's reference times bin-for-bin, including flow bins.
The references used stage 0 for raw bar times, stage 2 with run constants
disabled for pixel-aligned correlation times, and stage 7 for fully corrected
times. All three selected the same 4,685 events and 41,627 detector hits;
317 hits fell in the displayed pixels 484–491. The check also verified the
eight canvases, pixel labels/order, four-column/two-row layout, independent
LE maxima, and matched DT/color scales. This check used 1 ns bins, DT display
−60 to +40 ns, DT fit windows −45 to −10 ns, and ECal display −10 to +50 ns
to include the available LH2 sample. It does not establish the cross-target
physics result. All eight PDFs and eight PNGs were exported; compact LE
statistics now sit on one line above the peak, and the correlation margins
allow space for the color-axis label.

The following display revision replaced the LE 15% headroom and shared DT
65% headroom with independent 5% headroom for both plot types, moved statistics
above the frames, and changed the default selected pixel to 469 (bar 29).
ACLiC compilation and 20,000-entry Run 6077 checks completed for bars 30 and 29,
each retaining 4,685 accepted events and 41,627 detector hits. Both DT reference
medians and their fit counts stayed unchanged. The bar-29 run selected 847
hits across its full bar; its canvases display pixels 468–475. PDF/PNG exports
were produced and the LE/DT layouts were visually checked, including accepted
fits and low-statistics labels. Event selection, histogram contents and the
fitting procedure were not changed by this display revision.

## Polygon-selection validation (2026-10-02)

ROOT ACLiC compilation passed. A synthetic tree with known raw, pixel-aligned,
and final times checked that the polygon uses the drawing coordinate, including
optional reference subtraction. All 50 displayed histograms matched expected
contents bin-for-bin, including flow bins, with cuts enabled and disabled.
The test also checked pixels without cuts, unchanged base event occupancies,
and stopping on missing, empty, or wrong-run cut files.

On the first 20,000 local Run 6077 entries, reference hit times were captured
from the unchanged `c916d26c` macro. With polygons explicitly disabled, the
revised macro retained its 4,685 base-admitted events, 41,627 detector hits,
and 847 bar-29 hits. Two **test-only** Run 6077 rectangles then rejected 33
detector hits, including 19 bar-29 hits, leaving 41,594 and 828 respectively.
All 50 displayed histograms matched independently filtered reference times
bin-for-bin in both cases, and before/after populations remained identical.
All eight PDF/PNG pairs exported, and the LE, LE–ToT, DT-fit, and detector
correlation layouts were visually checked. Existing Podd dictionary/autoload
warnings occurred while opening real data; tree reading and comparisons
completed successfully.

The default Run 5710 archive loaded all 102 cuts, and its run metadata
correctly prevented implicit use on Run 6077. Local Run 5710 input is absent:
these checks validate software enforcement, not the physical acceptance of
the archived polygons with the current Run 5710 master. No calibration or
cut file was changed, and no new production physics plots were generated.
