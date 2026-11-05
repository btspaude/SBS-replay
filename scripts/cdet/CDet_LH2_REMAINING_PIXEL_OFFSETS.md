# LH2 methodology for remaining CDet pixel offsets

Updated 2026-11-04 to describe the current `extractBarOffsets()` implementation
in [Plot_CDet_PulseStudy.C](Plot_CDet_PulseStudy.C). Results are diagnostic
proposals, not approved calibration constants. The routine never updates the
master `.dat` or replay database.

## Purpose and analysis framework

Use LH2 statistics to investigate pixels inadequately constrained by cross-target
data. Inspect spectra first and select local pixel IDs with coherent peaks as
candidate donors. Work bar by bar; the selected donors and grouping can differ
between bars.

**Successful donors retain their individual fitted centers and offsets. Only
pixels without usable individual fits inherit a weighted donor result.**
Inheritance can use the whole bar or each eight-pixel half independently.

The main `Plot_CDet_PulseStudy()` reads the ROOT files once and saves event-aligned
vectors. The secondary `extractBarOffsets()` uses those vectors without rereading
the ROOT files. Its timing variable is `vPulseDt`, the stored branch
`earm.cdet.pulse.ecal_residual`:

$$
dt_i = t_{\mathrm{ECal}} - t_{\mathrm{CDet},i}^{\mathrm{corrected}}.
$$

This retains the replay's time-walk, run timing, and existing pixel-offset
corrections. LH2 uses ECal timing slope `p1 = 0`. A pixel lacking an independent
measured offset can still have a nonzero existing constant.

## Calling the routine

The current argument order is:

```cpp
extractBarOffsets(
    bool viewOnly = true,
    int myBar = 0,
    const std::vector<int>& goodPixels = {},
    bool fullNINO = true,
    const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
    const char* calibrationFile = "CDet_calibration_dt.dat",
    bool writeCSV = false,
    const char* csvFile = nullptr);
```

For example, after loading the run with the main macro:

```cpp
// Fit local pixels 3, 5, and 10; inherit separately within each half; save CSV.
extractBarOffsets(false, 19, {3, 5, 10}, false,
    "CDet_run6077_pulse_pixel_timing.conf",
    "CDet_calibration_dt.dat", true);
```

| Argument | Current behavior |
| --- | --- |
| `viewOnly` | `true`: histogram-only 4x4 canvases/PDFs, with no calibration-file read, fits, offset calculation, or CSV. `false`: fit and calculate proposals, overlay accepted fits, and optionally write CSV. |
| `myBar` | Global, zero-based bar ID, 0–167. |
| `goodPixels` | Local pixel IDs, 0–15, to attempt fitting. Empty means attempt all available pixel histograms. |
| `fullNINO = true` | `fitMode = 2`: one donor group covering pixels 0–15. |
| `fullNINO = false` | `fitMode = 1`: separate donor groups for pixels 0–7 and 8–15. |
| `configFile` | Shared configuration; requires `config.version: 1`. |
| `calibrationFile` | Master calibration file containing `[PixelOffsets]`, from the replay that produced the ROOT files. |
| `writeCSV` | Enables CSV output; defaults to false. |
| `csvFile` | Optional explicit `.csv` output path. Otherwise uses the default described below. |

`goodPixels` selects donors in the study bar only. It does not select the control
pixels: the routine attempts all available histograms in the configured control
bar. A single selected pixel can supply the full bar. In half-bar mode, each
half needs its own usable donor; no result is borrowed across the boundary.

## Execution and current cuts

1. Read the configuration and validate bar/pixel IDs.
2. In fit mode, read existing constants once from `[PixelOffsets]`. Each row is global pixel
   ID, offset in ns, and hit count; the hit count is not a validity flag. The
   mapping is `globalID = 16 * myBar + localPixel`. Missing constants for the
   study bar cause a return instead of silently substituting zero. The routine
   cannot verify that the supplied constants match the ROOT-file replay.
3. Fill separate dt histograms for the study and control pixels in the saved
   event/hit loops. All study pixels are filled, even when only some are selected
   for fitting. Both study and control histograms use `display.ecal_cdet_dt_min/max`
   and `pulse_study.dt_bin_width_ns` from the configuration. The 6077 config
   currently gives 90 bins over -90 to 0 ns; 5710 gives 300 bins over -150 to
   +150 ns. Both use 1 ns bins.
4. Draw all 16 study spectra and all 16 control spectra on separate 4x4
   canvases and save PDFs. View mode returns here.
5. In fit mode, fit control histograms and determine the reference `mu0`, then
   fit the selected study pixels, retaining each successful individual result.
6. Within each group, use successful donors to assign results to the other
   pixels. Print proposals, overlay accepted fits, update the PDFs, and
   optionally write the CSV.

The saved sample already passed the main macro's ECal energy/time cuts, replay
ECal eligibility, spatial and calibration flags, and strict ToT bounds. The 6077
configuration uses `-10 < ECalAdcTime < 10` ns; 5710 uses `10 < ECalAdcTime < 35`
ns, unless overridden in the main call.

The offset routine itself requires finite dt and rx, and
`abs(rx) <= analysis.x_difference_max` (default 0.08 m). It currently reads
secondary ECal and LE/ToT settings into variables but **does not apply those
additional cuts**. It does not apply an ellipse or a multiplicity selection.
Changing secondary cut settings alone therefore does not change this extractor's
selection. Rejected loading events cannot be recovered by a secondary routine.

The number of dt bins is `ceil((DtMax - DtMin) / BinWidth)`, matching the other
PulseStudy routines. Bounds stay exactly as configured; if their span is not an
integer multiple of the requested width, the actual bin width is slightly
smaller. Use an integer-nanosecond span for exact 1 ns bins. Display bounds do not
reject pulses: values outside them enter underflow/overflow bins. Fits only see
the histogram's available range, so keep the desired fit windows within it.

## Peak fitting and acceptance

Each pixel histogram is fitted with `gaus(0)+pol1(3)`:

$$
f(t) = A\exp\!\left[-\frac{(t-\mu)^2}{2\sigma^2}\right] + b_0 + b_1 t.
$$

`GetParameter(1)` supplies the peak center, and `GetParError(1)` supplies its fit
uncertainty. The Gaussian width `GetParameter(2)` is **not** the uncertainty used
for weighting.

| Fit | Configuration keys | Default window (ns) |
| --- | --- | --- |
| Control pixels | `pulse_study.control_fit_min/max` | -40 to 0 |
| Study donors | `pulse_study.fit_min/max` | -50 to 0 |

Both shared 6077/5710 configs now explicitly set these keys to the default
windows above. They are starting settings, not validated peak windows for every
run/bar; inspect the spectra and adjust before fitting. Initial estimates come from
the histogram maximum/minimum inside the fit window, with an initial Gaussian
width equal to one tenth of the window width. The ROOT fit option is `QRS0`:
quiet, use the function range, retain the fit result, and do not draw the fit.
Accepted attached fit curves are subsequently made visible on the canvases;
failed fits and inherited assignments are not drawn as fitted curves.

A fit contributes only if its returned status is zero, its result exists,
its covariance status is at least 2, its center and positive center uncertainty
are finite, its center lies strictly inside the fit window, its amplitude is
positive, and its Gaussian width has finite nonzero magnitude. There is currently
no minimum-entry, significance, chi-square, inter-pixel agreement, or systematic
stability requirement. These numerical checks do not establish a physical peak;
manual donor selection and spectrum inspection remain essential.

A selected pixel whose fit fails does not become a donor. It can inherit from
other successful donors in its group. A group with no usable donor remains
unresolved. If no usable control reference exists, the routine returns before
proposing offsets or writing a CSV. The histogram PDFs have already been saved.

## Reference and calculation formulas

All centers, offsets, corrections, and reported uncertainties are in ns.
Let `s_i` denote the fitted uncertainty on a center, not its Gaussian width.
The formulas below are the independent-fit approximation implemented in the code.

### Control reference

The configured `pulse_pixel_timing.control_bar` defaults to bar 30. For accepted
control fits in the set $C$, define inverse-variance weights:

$$
w_k = \frac{1}{s_k^2},\qquad
\mu_0 = \frac{\sum_{k\in C} w_k\mu_k}{\sum_{k\in C} w_k},\qquad
s_0 = \sqrt{\frac{1}{\sum_{k\in C}w_k}}.
$$

Here `mu0` is the empirical LH2 reference from already-calibrated control pixels.
It is distinct from the historical cross-target reference. The study donors do
not redefine it. Inspect consistency among control peaks, and where practical
check a few known calibrated pixels in other bars/layers.

### Successful individual donor

For donor $i$, retain its own fitted $\mu_i$ and $s_i$:

$$
c_i = \mu_i - \mu_0,\qquad
s_{c_i} = \sqrt{s_i^2+s_0^2},\qquad
O_i^{\mathrm{new}} = O_i^{\mathrm{existing}} + c_i.
$$

The correction is added to CDet time, so a positive correction moves dt downward.
A peak at -23 ns against a -26 ns reference needs `c = +3 ns`. If its existing
offset is +5 ns, the proposed final constant is +8 ns, not +3 ns.

### Pixels inheriting a group result

Let $D_g$ be the accepted donors in the full-bar or half-bar group $g$:

$$
W_g = \sum_{j\in D_g}\frac{1}{s_j^2},\qquad
\bar\mu_g = \frac{\sum_{j\in D_g}\mu_j/s_j^2}{W_g},\qquad
s_g = \sqrt{\frac{1}{W_g}}.
$$

For a pixel $i$ without a usable individual fit in that group:

$$
\mu_i^{\mathrm{assigned}}=\bar\mu_g,\qquad
s_i^{\mathrm{assigned}}=s_g,\qquad
c_i=\bar\mu_g-\mu_0,
$$

$$
s_{c_i}=\sqrt{s_g^2+s_0^2},\qquad
O_i^{\mathrm{new}}=O_i^{\mathrm{existing}}+c_i.
$$

The reference uncertainty is added once after combining donors. It does not
shrink with the donor count. Nor does an inherited error shrink because the
result is copied to more recipient pixels. A single donor gives its own center
and center uncertainty as the inherited group estimate.

The current method averages the **corrected dt centers**, not the donors' final
database constants. Recipients share a residual correction, but their final
constants can differ if their existing constants differ. This supersedes the
earlier plan to force one identical final constant across a group. It assumes
that the remaining timing residuals are comparable within the chosen group;
existing unequal corrections can affect that assumption and must be considered
when selecting donors.

### Worked example

Suppose `mu0 = -26.0 ± 0.2 ns` and two accepted donors in one group have:

| Pixel | Fitted mu | Center error | Own correction c | Correction error |
| --- | --- | --- | --- | --- |
| 3 | -23.0 | 0.4 | +3.0 | 0.447 |
| 5 | -21.0 | 0.8 | +5.0 | 0.825 |

Their weights are 6.25 and 1.5625, giving `groupMu = -22.6 ns` and
`groupMuError = 0.358 ns`. Other pixels in the group inherit `c = +3.4 ns`
with error `0.410 ns`. Pixels 3 and 5 retain their separate +3.0 and +5.0 ns
corrections; they are not replaced by +3.4 ns. Each pixel then adds its own
existing constant to obtain its proposed final offset.

## Meaning and limits of the uncertainties

The reported errors include fitted-center statistics and the reference error.
They omit uncertainty from background/fit-model choice, fit-window and cut
variation, donor selection, genuine pixel-to-pixel variation, and the assumption
that unmeasured pixels share the donors' residual. There is no scatter-based
error inflation in the code.

Existing calibration constants are treated as fixed inputs; their uncertainties
and correlations are not read or propagated. Consequently `barOffsetError` is
the same as the correction error under this fixed-input convention, not a full
uncertainty budget for the final calibration.

More generally, correlations would enter as:

$$
\operatorname{Var}(\mu_i-\mu_0)
= s_i^2+s_0^2-2\operatorname{Cov}(\mu_i,\mu_0).
$$

If the study bar equals the control bar, the code leaves correction errors
undefined (`NAN`, blank in CSV), because the fits overlap the reference and that
covariance is not implemented. Distinct bars are treated as independent even
though they may share events or systematic effects; that is an approximation.

Corrections share the reference uncertainty: for independent donor centers,
$\operatorname{Cov}(c_i,c_j)=s_0^2$ for distinct donors. Recipients inheriting
the same group estimate also share the donor-mean uncertainty. Their CSV rows
are not independent measurements and must not be averaged as independent data
to claim a smaller uncertainty.

## Histogram canvases and PDFs

Every valid call with saved events creates two 4x4 canvases: the study bar and
its configured control bar. Pads 1-16 show local pixels 0-15 in row-major order.
All pixels are shown, including empty histograms and pixels outside `goodPixels`.
Histograms are blue; in fit mode, accepted individual fits are overlaid in red.
An inherited correction does not generate a fit curve for its recipient.

Unused channels use the same global ID list as `PlotElastic.C`. A solid black
rectangle covers their statistics box, including when the histogram is empty.
Empty instrumented pixels are not marked unused just because they have no hits.
The statistics boxes display entries, mean, and standard deviation. This visual
unused-channel marker does not itself exclude a channel from fitting or CSV
inheritance.

Only PDF plot files are saved:

```text
PusleStudy/run6077_bar019_dt_pixels_view.pdf
PusleStudy/run6077_bar019_control030_dt_pixels_view.pdf
PusleStudy/run6077_bar019_dt_pixels_fit.pdf
PusleStudy/run6077_bar019_control030_dt_pixels_fit.pdf
```

Each call writes the two files for its mode; run and bar IDs follow the call and
configuration. View and fit filenames are separate. Repeating a call in the same
mode replaces its PDFs and the two previous offset canvases. In an interactive
ROOT session, canvases remain open; the routine does not force batch mode.

```cpp
extractBarOffsets(true, 19);                    // Histograms only, all 16 pixels.
extractBarOffsets(false, 19, {0, 4, 8}, false); // Fit chosen donors, separate halves.
```

`viewOnly = true` takes precedence over `writeCSV = true`: viewing never writes
an offset CSV. Fit mode saves initial histogram PDFs before fitting, so spectra
remain available even if a usable control reference cannot be established.
Calibration/database inputs remain read-only.

## CSV output and proposal handling

When `viewOnly = false` and `writeCSV = true`, a successful reference calculation produces all 16 local
pixel rows. The default filename is:

```text
PusleStudy/run<run>_bar<zero-padded-bar>_<full-or-halves>_offsets.csv
```

For the example call above: `PusleStudy/run6077_bar019_halves_offsets.csv`.
An explicit final `csvFile` argument overrides the path and must end in `.csv`;
its parent directory must exist. Reusing a filename replaces that CSV.

| CSV column | Meaning |
| --- | --- |
| `ID` | For example, `Bar 19 pixel 3`; local pixel numbering. |
| `mu` | Individual fitted center for successful donors; weighted group center for recipients. |
| `mu error` | Individual center error or inherited donor-mean error, excluding the reference error. |
| `c` | Additional correction `mu - mu0`, not the final database constant. |
| `c error` | Quadrature center/reference error as above; blank when covariance treatment is missing. |

Unresolved pixels have an ID and four blank numeric fields. A numeric zero is
a valid result, not the unresolved marker. All 16 IDs are emitted, including
unused channels; receiving a proposed value does not establish that a channel
is instrumented or working. The routine has no mask protecting established
calibrations from appearing in the proposal.

The five-column CSV does not record donor/inherited provenance, `mu0`, run cuts,
existing constants, or final constants. Preserve the call, console output,
configuration, and matching calibration file alongside it so the result is
interpretable. The console distinguishes donors from inherited assignments and
prints proposed final constants. CSV writing does not alter the saved event
vectors, calibration master, or replay database. Review which proposals should
actually be applied before a separate calibration update.

## Current database offset coverage

Audit date: 2026-11-04. The replay loads `earm.cdet.timing.pixel_offset`
from [the replay database](../../DB/db_earm.cdet.dat). All 2688 entries are
present and finite, and match the `[PixelOffsets]` section of
[CDet_calibration_dt.dat](CDet_calibration_dt.dat) exactly. The third column
in that master section is a hit count, not a validity flag.

The following **47 bars contain zero-valued offsets for instrumented pixels**.
Bar numbers use global, zero-based numbering (L1: 0–83; L2: 84–167).

| Layer | Bars containing zero-valued offsets |
| --- | --- |
| L1 (22 bars) | 0–3, 7, 22, 24, 39–44, 53, 70–74, 81–83 |
| L2 (25 bars) | 84–87, 89, 91, 106, 123–128, 130, 136–137, 154–159, 165–167 |

There are **724 zero-valued entries**: **630 instrumented pixels** and
**94 unused channels**. Unused channels are identified by the `unused` source
label in the [archived Run 5710 fit results](CDet_run5710_halfbar_aligned_final_archive/CDet_pixel_timing_fit_results_run5710_fixed_effects_final.dat).
Of these 47 bars, 42 have all 16 entries zero (14 instrumented and two unused).
The five partially zero bars are:

| Bar | Instrumented pixels with zero offsets | Local pixel IDs with zero offsets (unused channels excluded) |
| --- | --- | --- |
| 22 | 13 | 0, 1, 3–11, 13, 14 |
| 72 | 8 | 0–7 |
| 74 | 6 | 1, 3–7 |
| 89 | 7 | 8–14 |
| 91 | 8 | 0–7 |

**There are no absent, nonfinite, or zero-valued offset entries outside these
47 bars.** All 1936 entries in the other 121 bars are finite and nonzero.

This is a numerical coverage audit, not a certification of independent pixel
fits. Nonzero constants can include inherited group offsets or retained
existing values, so pixels needing further investigation can exist outside
this list. Conversely, zero alone is not a fit-validity flag. The archived
fit results describe an earlier calibration pass and do not reproduce every
current constant; use them for provenance with that limitation, rather than
equating their fit labels with the current database values.

## Diagnostic idea: dt spectra in different rx ranges

Compare apparently background-dominated dt spectra in several ranges of the
stored pulse x residual, `rx = earm.cdet.pulse.ecal_x_residual`, for individual
pixels and bars. Use `vPulseDt` for dt. A peak localized in rx could be obscured
when the full accepted rx range is combined; this is a hypothesis to test,
not an assumption that the background contains recoverable signal.

Keep the other selections and dt binning fixed, record the rx boundaries and
counts, and compare both yields and normalized shapes. Initially leave the
dt–rx ellipse disabled so it does not impose a timing peak. Any candidate peak
should persist under modest changes in rx boundaries and be checked on another
sample. The saved vectors only cover the spatial acceptance retained at loading.

## Selection, validation, and records

- Use background studies to select credible signal: ECal cuts, geometry,
  ToT, and single-/multiple-pulse comparisons. Count multiplicity across the
  bar after pulse cuts, then select the donor pixel.
- Keep the candidate timing search sufficiently broad. Check that loading,
  corrected-LE cuts, or pairing requirements have not already rejected the
  displaced peak. Secondary routines cannot recover rejected pulses.
- Check peak/reference stability against reasonable selection and fit-window
  changes. A successful numerical fit alone does not establish signal.
- Inspect corrected spectra and validate on an independent sample or run when
  possible. A donor fit uncertainty does not quantify the additional uncertainty
  from assuming unmeasured recipients share the donors' timing residual.
- Record run/segments, cuts, reference pixels and estimator, donor fits, old
  constants, proposed final constants, and whether each value is individually
  fitted, inherited from a half-bar, inherited from a bar, or unresolved.

Keep diagnostic proposals separate from approved constants; database updates
remain a separate reviewed step.

Related notes: [cross-target calibration](CDet_CROSSTARGET_TIMING_CALIBRATION.md)
and [earlier offset investigations](CDet_FURTHER_PIXEL_OFFSET_INVESTIGATIONS.md).
