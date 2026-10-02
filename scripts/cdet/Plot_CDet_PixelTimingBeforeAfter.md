# Cross-target pixel timing before and after calibration

`Plot_CDet_PixelTimingBeforeAfter.C` makes eight canvases for DNP slides:

- Six bar canvases: LE, LE versus ToT, and ECal-minus-CDet Δt, each before and
  after calibration. Each shows **four columns and two rows**, retaining the
  middle two rows of the original 4×4 layout.
- Two full-detector canvases: CDet time versus ECal time, one with pixel timing
  offsets applied and one with the full timing calibration applied.

Any logical pixel ID selects its containing 16-pixel bar: `bar = pixel / 16`.
The bar canvases display pixel IDs `16*bar + 4` through `16*bar + 11`. For
example, pixel 485 selects bar 30 and displays **484–491**. Layer 2 pixel IDs
are also supported. Bar numbering is global, 0–167. The selected bar does not
restrict the full-detector correlation or the detector-wide DT reference fits.

## Run it

From `gocdetscripts` / `SBS-replay/scripts/cdet`, start a fresh ROOT session:

```cpp
.L Plot_CDet_PixelTimingBeforeAfter.C+
Plot_CDet_PixelTimingBeforeAfter("CDet_run5710_projection.conf", 485);
```

The existing run configuration supplies event limits, segment limits and
analysis cuts. Input defaults to `OUT_DIR`. To override the destination,
input directory and event limit:

```cpp
Plot_CDet_PixelTimingBeforeAfter(
    "CDet_run5710_projection.conf", 485,
    "dnp_run5710_bar30", "/path/to/Rootfiles", 100000);
```

The event override defaults to `-2` (use configuration); `-1` reads all selected
entries. As in the master, `analysis.first_event` is not used. This small macro
handles one run per call; it fits peaks for display but does not combine runs
or update calibration constants.
`analysis.calibration_stage` is intentionally ignored: every call produces both
timing states. The bar plots' before state is uncorrected; the detector
correlation's before state has pixel offsets applied, as labeled on its canvas.

The next arguments control display binning, not acceptance. Four optional
arguments at the end set the before/after DT fit windows:

```cpp
// LE bin width/min/max, ToT bin width/min/max, savePlots,
// dt bin width/min/max, ECal-time bin width/min/max,
// before DT fit min/max, after DT fit min/max
Plot_CDet_PixelTimingBeforeAfter(
    "CDet_run5710_projection.conf", 485, "dnp_run5710_bar30", nullptr, -2,
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
cdet_timing_before_after/CDet_run5710_bar030_le_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar030_le_after.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar030_le_vs_tot_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar030_le_vs_tot_after.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar030_ecal_cdet_dt_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_bar030_ecal_cdet_dt_after.{pdf,png}
cdet_timing_before_after/CDet_run5710_detector_cdet_t_vs_ecal_t_before.{pdf,png}
cdet_timing_before_after/CDet_run5710_detector_cdet_t_vs_ecal_t_after.{pdf,png}
```

Reusing a destination for the same run/bar replaces those figures. The two
detector correlation files are shared across selected bars. Earlier
`barNNN_cdet_t_vs_ecal_t_*` files in an existing destination are left untouched;
use the new `detector_cdet_t_vs_ecal_t_*` files for the aggregate views.
The PDF is the vector export for slides. Before/after pairs use the same
binning and x ranges. **LE panels have independent count-axis maxima**, each
15% above its own tallest bin, to reduce empty space when calibration sharpens
the peak. Their axes remain counts, not normalized distributions. DT panels
retain matched count scales, and the 2D plots retain matched color scales.
The compact `N` and `SD` labels use the entire
selected sample, including values outside the display range. No underflow or
overflow boxes are drawn. LE retains those two labels. The DT panels also
show the fitted signal centroid `mu_i`, Gaussian width `sigma_fit`, their
fit uncertainties in ns, and dimensionless chi-square/NDF. SD describes the
whole distribution, including background; `sigma_fit` describes the fitted
signal peak. Neither is automatically an intrinsic CDet timing resolution.
The DT panels reserve space above the peaks for the labels.

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
master good-hit sample: no extra extractor-only 1–12 GeV energy cut or manual
LE–ToT polygons are added. Consequently, these plots are not a reproduction
of the commissioned Run 5710 fit constants, whose calibration passes used
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
No additional elastic selection, layer pairing, ellipse selection, or saved
LE–ToT polygon is imposed. The bar and middle-eight-pixel display restrictions
apply only to the LE, LE–ToT and DT canvases, after evaluating the full
detector's event-occupancy cuts. The CDet/ECal correlations and the temporary
DT histograms defining the detector-wide references use all accepted pixels.

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

The macro is solely a plot producer. It uses the same hit selection and timing
formula as the master's `vGoodLe`, but fills the histograms immediately and
retains no run-wide event vectors. Per-event candidate indices are local and
discarded after each event. Only the canvases and optional PDF/PNG exports
remain; no event tree, analysis-vector interface, or calibration output is
created. Shared calibration readers still load the existing constants into
the master's calibration variables.

Two differences from the old plotting shortcut are deliberate: the old
`plotPaddles()` reads `vPaddleGoodLe`, which is filled with pixel-offset-corrected
times but is not updated by the later ECal/time-walk/global corrections. Also,
the master's stage 0 still adds the run's global shift. Neither is suitable as
the respective fully corrected or completely uncorrected view requested here.

## Validation (2026-10-01)

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
