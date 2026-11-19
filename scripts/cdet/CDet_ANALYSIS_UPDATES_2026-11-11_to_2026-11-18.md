# CDet analysis updates: November 11–18, 2026

Seven-day lookback ending November 18, 2026 (America/New_York), based on the local SBS-replay `tdccalib` history and the changes committed with this document. This records implemented analysis behavior, not conclusions about detector performance. Existing offset extraction and the original vector-loading architecture predate this period unless a change is explicitly identified below.

## Files and workflow

The work is concentrated in:

- [Plot_CDet_PulseStudy.C](Plot_CDet_PulseStudy.C): load the run once, keep event-indexed vectors, then call secondary plotting routines.
- [Plot_CDet_PulseSideband.C](Plot_CDet_PulseSideband.C): independent loading and background fits for spatial sidebands.
- [CDet_run6077_pulse_pixel_timing.conf](CDet_run6077_pulse_pixel_timing.conf): LH2 configuration.
- [CDet_run5710_pulse_pixel_timing.conf](CDet_run5710_pulse_pixel_timing.conf): cross-target configuration.

Secondary routines reread their configuration on each call. They can narrow the saved sample, but cannot recover events or pulses rejected during loading. Reload the main macro when changing a loading cut and needing a different saved sample. Repeated plotting calls replace each routine's previous histograms and canvases.

## 1. Signal fits and fit-quality information

`getBarSigma` fits each bar and the detector sum with a Gaussian signal plus a linear background, `gaus(0) + pol1(3)`, within `pulse_study.fit_signal_min/max`. The current signal fit is not a double Gaussian.

Implemented outputs include:

- Control-bar dt spectrum with its fit, signal sigma, uncertainty, and reduced chi-square.
- Detector-summed dt spectrum with its own fit and sigma.
- Sigma versus global bar number: valid layer-1 points red, layer-2 points blue, and unusable/insufficient-statistics fits marked black at zero. Those black markers are not measured zero widths.
- An inverse-variance weighted mean of the usable per-bar sigmas in the selected graph range, also drawn on the sigma plot.

For fitted widths `sigma_i` with fit errors `s_i`:

\[
w_i=1/s_i^2,\qquad
\bar\sigma=\frac{\sum_i w_i\sigma_i}{\sum_i w_i},\qquad
\delta\bar\sigma=(\sum_i w_i)^{-1/2}.
\]

The mean excludes failed fits and their zero placeholders. Its quoted uncertainty assumes independent fit errors and excludes systematics. It is distinct from the width fitted to the detector-summed spectrum.

The signal and sideband fit displays both gained reduced chi-square, `chi2 / ndf`, with unavailable cases identified when appropriate.

## 2. Independent spatial-sideband study

Added `Plot_CDet_PulseSideband.C` and its secondary `getSigmaBars` routine, using the same run configuration filenames as PulseStudy without requiring a separate sideband config.

- The default inner sideband boundary is `analysis.x_difference_max`, normally 0.08 m.
- The default outer boundary is twice that value: `0.08 < abs(rx) <= 0.16 m` for the usual settings.
- The loading call supports an outer-bound override; the default y-residual bound is 0.36 m.
- Sidebands use explicit spatial selections rather than requiring the replay signal-region spatial flag. They do not apply the dt ellipse.
- Background fits compare quadratic (`pol2`), fourth-order polynomial (`pol4`), Gaussian, and Lorentzian shapes.
- Fits use `pulse_study.fit_min/max`, excluding histogram-bin centers in the inclusive `fit_signal_min/max` interval. The full histogram and full-range model overlays remain visible; the excluded data are not replaced by zeros.
- Fit viability requires usable statistics and bins on both sides of the excluded interval, plus fit-quality checks.
- The sideband sigma-versus-bar plot reports the Gaussian **background** width, not signal timing resolution. Polynomial fits have no sigma; Lorentzian gamma is a half-width at half-maximum.

Output is under `PusleStudy/sidebands/<tag>/`, with interactive canvases. Lagrange interpolation and additional background-shape models were discussed but were not implemented in this period.

The later optional-upper-cut switches and text-report cleanup described below apply to `Plot_CDet_PulseStudy.C`; they were not propagated to the separate sideband macro. The sideband loader still applies its configured energy and ToT upper bounds.

## 3. Peak pulse yields for the control bar and all bars

`getBarSigma` now calculates the number of selected pulses inside each usable bar fit's peak window, initially for the control bar and subsequently for every bar.

\[
W_b=[\mu_b-k\sigma_b,\ \mu_b+k\sigma_b],\qquad
Y_b=\frac{N_{\text{selected pulses in bar }b\text{ within }W_b}}
{N_{\text{ECal-selected events}}}.
\]

`k` is `pulse_study.peak_nsigma` (default 2), with a call override available. The numerator counts saved pulse values in the fitted window. The dedicated ECal denominator includes events with zero accepted CDet pulses, and receives the secondary ECal cuts. Secondary CDet cuts do not reduce this denominator.

Added `peak_yield_vs_bar.pdf`, alongside the control-bar yield annotation. These are raw pulse yields per ECal-selected event, with background included. They are not rates in Hz or background-subtracted detector efficiencies; multiple pulses in an event can contribute separately.

## 4. LE plots, comparison canvases, and run labels

`plotPulseStudyBars` gained corrected-LE versions of the existing raw-LE plots:

- Summed bar LE spectra and 4×4 per-pixel LE spectra.
- Summed and per-pixel LE-versus-ToT heat maps.
- Fixed bar-30 raw-LE versus ECal ADC-time and corrected-LE versus ECal ADC-time heat maps.

The LE selection continues to use **raw LE**, even when plotting corrected LE. Nonfinite corrected LE prevents filling its own corrected-LE plots, without discarding the pulse from otherwise valid raw-LE/dt plots.

The LE/ToT comparison canvas now has a wide top panel overlaying dt before and after the secondary raw-LE window, with two lower panels showing dt versus ToT before and after the secondary ToT window. The comparison omits only the indicated secondary window; loading cuts remain in all samples.

Added `analysis.target` for titles such as `Run 6077 (LH2)` or `Run 5710 (Cross target)`. The target label does not change selections or add a target string to output filenames.

## 5. Y-position timing correction: current implementation

Implemented an optional correction based on the distance from the projected ECal hit to the readout end of the 50 cm bar. It is referenced to the bar center, preserving a zero correction there.

For each pixel:

\[
y_{\rm proj}=y_{\rm ECal}\frac{z_{\rm pixel}}{z_{\rm ECal}},\qquad
 y_{\rm readout}=y_{\rm DB}+s(0.25\,\mathrm{m}),\qquad
 d=|y_{\rm readout}-y_{\rm proj}|,
\]

where `s = +1` on the LHS and `s = -1` on the RHS. The helper identifies the side using `pixelID % 1344 < 672` for the LHS, repeating the mapping in each layer.

The initial `n/c` model was replaced on November 18 by an effective propagation delay that accounts for indirect light paths:

\[
\alpha=7\,\mathrm{ns/m},\qquad
\Delta t_y=\alpha(d-0.25\,\mathrm{m}),\qquad
 dt_y=dt+\Delta t_y,\qquad
 LE_{\rm corr,y}=LE_{\rm corr}-\Delta t_y.
\]

- Current key: `pulse_study.y_propagation_ns_per_m: 7.0` in both configs.
- This replaces `pulse_study.y_correction_refractive_index`; the old key is not used.
- An explicit numeric propagation argument now has units **ns/m**, not refractive index.
- The correction is **off by default** and enabled through the routine's boolean argument.
- `getBarSigma`: applies to bar 30 and its contribution to the detector sum.
- `plotPulseStudyBars`: applies to both selected bars and the fixed bar-30 LE/ECal-time diagnostic.
- Selections use original stored dt and raw LE. Saved vectors are unchanged.

The 7 ns/m value is the chosen effective model, not a newly fitted propagation constant.

Also added two diagnostic heat maps to `plotPulseStudyBars`, one per selected bar: **dt before y correction versus projected ECal y**. They always fill from saved `vPulseDt` (`pulse.ecal_residual`), even when y correction is enabled elsewhere. Thus “before correction” here means before this optional y correction, not before the replay's existing timing corrections. Projected y is horizontal; dt is vertical. The projected-y display uses the position binning and zooms to occupied bins plus two bins of padding.

## 6. Cut controls and current defaults

Upper ECal-energy and ToT limits can now be disabled in PulseStudy while retaining their numerical settings for future scans.

| Setting | Default/current value | Effect |
| --- | --- | --- |
| `analysis.apply_ecal_energy_max` | 0 | Loading energy upper cut off |
| `analysis.apply_tot_max` | 0 | Loading ToT upper cut off |
| `pulse_study.apply_ecal_energy_cut` | 1 | Secondary energy lower cut on |
| `pulse_study.apply_tot_cut` | 1 | Secondary ToT lower cut on |
| `pulse_study.apply_ecal_energy_max` | 0 | Secondary energy upper cut off |
| `pulse_study.apply_tot_max` | 0 | Secondary ToT upper cut off |
| `analysis.le_max` | **30 ns in both configs** | Raw-LE maximum, explicitly adopted in this update |

PulseStudy loading now replaces replay `spatial_pass` with explicit finite-residual checks: `abs(rx) <= analysis.x_difference_max` and `abs(ry) <= 0.36 m`. The x bound defaults to 0.08 m and can be widened in the config; rerun the main analysis to refill vectors. Replay ECal eligibility and calibration validity remain required. The captured loading-cut report records the actual rx/ry limits, even if the config is edited before plotting.

Loading ECal-energy and ToT lower bounds remain active. Secondary lower windows have their existing toggles. Secondary upper cuts require the corresponding secondary window to be enabled. Display maxima do not enable selection cuts.

For a ToT-upper-cut scan without reloading, load with `analysis.apply_tot_max: 0`, enable `pulse_study.apply_tot_cut` and `pulse_study.apply_tot_max`, then vary `analysis.tot_max` between plotting calls.

## 7. Position correlations and coordinate-difference spectra

Expanded `plotPulseStudyPositionCorrelations` to cover:

- The selected control bar (default 30), the full detector, layer 1, and layer 2.
- Channel CDet x, aligned/corrected CDet x, and DB channel y versus both native ECal coordinates and ECal coordinates projected to each pulse's z.
- Separate per-layer 1D spectra of `x - x_ECal_projected` and `x_corr - x_ECal_projected`.
- The same two 1D subtraction spectra specifically for **fixed bar 30**, independent of the control-bar argument.
- Full-detector 1D channel-CDet-x and native/unprojected ECal-x spectra. CDet counts accepted pulses; native ECal counts each event once if it has an accepted pulse.

The redundant study-bar group was removed from this routine. Its third argument is the control bar. Control-bar correlation displays auto-zoom with padding; detector/layer plots retain the configured common bounds. Raw/corrected control-x plots share occupied bounds; L1/L2 comparisons share axes and color scales.

Current position binning is 0.005 m (5 mm), with default full-coordinate bounds -1.6 to +1.6 m. The ECal-subtraction spectra have separate display limits, currently -0.20 to +0.20 m. These direct subtractions do not apply an extra DB residual-offset subtraction and should not be treated as aliases of stored rx. DB y is a channel center, not a reconstructed interaction position.

### Optional signal dt window for the single-pulse position plots

```text
pulse_study.position_apply_dt_cut: 0
pulse_study.position_dt_min_ns: -36.0
pulse_study.position_dt_max_ns: -20.0
```

Those are the current 6077 values. Run 5710 starts at -20 to +10 ns. Set the toggle to 1 to apply the inclusive stored-dt window to the single-pulse position plots, including fixed-bar-30 subtraction spectra and the native ECal event spectrum. These are adjustable starting windows, not guaranteed optimal signal bounds.

### New same-event L1/L2 plots, triggered by bar-30 signal

Added channel and aligned versions of:

- `x2 vs x1`, with bar-30 L1 x horizontal and L2 x vertical.
- `x2 - x1`, with that exact subtraction sign and no ECal projection subtraction.

An event qualifies when it has at least one bar-30 pulse passing the secondary ECal, raw-LE, ToT, and rx/ellipse cuts **and** the configured signal dt window. This bar-30 window is always required for these layer plots, even when `position_apply_dt_cut` is 0 for the other plots.

Every qualifying bar-30 pulse is combined with every saved L2 pulse in that same event. L2 gets loading cuts only: no additional secondary dt, rx, LE, or ToT selection. Finite coordinates are required separately for each coordinate version. Multiple qualifying bar-30 pulses produce all combinations; these are not replay-matched pairs or confirmed tracks, and counts are combinations rather than events.

The new layer-difference display keys are `pulse_study.layer_difference_min_m/max_m`, currently -1.6 to +1.6 m (narrowed from the initial ±3.2 m), using the position bin width. Display boundaries do not discard combinations. Older configs missing the position dt bounds fall back to `fit_signal_min/max`.

### Bar-30 multiplicity

Added `bar030_hits_per_event.pdf`: one entry per event passing secondary ECal cuts with at least one saved bar-30 pulse. It counts **all saved bar-30 pulses before secondary pulse/dt cuts**, including out-of-peak pulses. The axis is total N; additional pulses equal N-1. Unit-width integer bins cover the observed maximum. Zero-bar-30 events are excluded, so this is a conditional multiplicity distribution.

The position-correlation routine now produces 37 canvases: 24 coordinate correlations, four layer/ECal differences, two fixed-bar-30/ECal differences, two detector x spectra, four layer-combination plots, and one multiplicity plot.

## 8. Plot cleanup and cut reports

Removed crowded cut boxes from PulseStudy canvases, reclaimed plotting space, and retained run/bar identifiers, useful legends, fit results, and the unused-pixel masks in offset displays.

Each PulseStudy plotting routine now writes a routine-specific `*_cuts.txt` beside its PDFs. The reports distinguish:

- Actual loading-cut metadata captured when the vectors were filled, including loading-call overrides.
- The secondary cuts used by the plotting routine.
- Display and fit settings, optional y correction, and before/after comparison exceptions.
- Special selections for layer combinations, multiplicity, offset extraction, and the legacy background routine.

This avoids incorrectly labeling the saved sample with loading cuts edited after analysis. Missing loading metadata is explicitly reported. A failed report write produces an error rather than silently proceeding to save unaccompanied new plots.

Existing selection differences are documented, not silently changed: `extractBarOffsets` applies its ECal-time and rx selections but no additional energy/LE/ToT/ellipse cuts; the legacy background routine applies its energy/time/ToT windows directly and uses corrected LE. The separate sideband macro has not received this text-report/plot-box cleanup.

## 9. Example calls and output locations

From `scripts/cdet`, with the normal ROOT/analyzer environment and input directory configured:

```cpp
.L Plot_CDet_PulseStudy.C+
Plot_CDet_PulseStudy("CDet_run6077_pulse_pixel_timing.conf");
plotPulseStudyBars("CDet_run6077_pulse_pixel_timing.conf", "bars_baseline");
plotPulseStudyBars("CDet_run6077_pulse_pixel_timing.conf", "bars_y", true);
getBarSigma("CDet_run6077_pulse_pixel_timing.conf", "sigma_y", true);
plotPulseStudyPositionCorrelations("CDet_run6077_pulse_pixel_timing.conf", "positions");
```

Tagged PulseStudy outputs go to `PusleStudy/<tag>/` (existing directory spelling preserved). Examples of new filenames:

- `plotPulseStudyBars_cuts.txt`
- `plotPulseStudyPositionCorrelations_cuts.txt`
- `bar030_Control_dt_before_y_vs_ecal_y_projected.pdf`
- `bar030_raw_x_minus_ecal_projected.pdf`
- `bar030_corrected_x_minus_ecal_projected.pdf`
- `bar030_signal_raw_x2_vs_x1.pdf`
- `bar030_signal_corrected_x2_vs_x1.pdf`
- `bar030_signal_raw_x2_minus_x1.pdf`
- `bar030_signal_corrected_x2_minus_x1.pdf`
- `bar030_hits_per_event.pdf`

Use a new tag to preserve a previous study. Reusing a routine and tag replaces its corresponding plots and report. The loading ECal plot/report and offset-extraction outputs retain their existing untagged `PusleStudy/` locations with run/bar names.

For the separate sideband study:

```cpp
.L Plot_CDet_PulseSideband.C+
Plot_CDet_PulseSideband("CDet_run6077_pulse_pixel_timing.conf");
getSigmaBars("CDet_run6077_pulse_pixel_timing.conf", "background_shapes");
```

## 10. Validation and limits

Synthetic ROOT checks during this work covered repeated-call cleanup, preserved event/pulse weighting, loading-cut provenance after config edits, enabled/disabled windows and upper cuts, nonfinite coordinates, empty samples, flow-bin retention, and unchanged saved vectors. Y tests checked both readout sides, correction signs, and pre-y-correction diagnostics. New layer tests checked all combinations, absence of cross-event associations, L2 acceptance outside secondary windows, fixed bar-30 gating, and multiplicity counted independently of signal selection. Representative exported PDFs were rendered and visually inspected.

These are implementation checks. Real-run fits, the best background model, the preferred signal windows, and the effectiveness of 7 ns/m still require evaluation on the farm data. No replay DB pixel offsets were changed by these plotting additions.

## Commit trail for the lookback

| Date | SBS-replay commit | Implemented change |
| --- | --- | --- |
| Nov 11 | `00e7a45c` | Signal-window Gaussian plus linear-background fits; run configs updated |
| Nov 11 | `5c61e702` | Independent spatial-sideband study |
| Nov 12 | `cee20e52` | Reduced chi-square on signal/sideband fit plots |
| Nov 12 | `5b9c885e` | Control-bar peak yield and ECal opportunity denominator |
| Nov 12 | `83e706aa` | Exclude signal interval from sideband fits; add Lorentzian |
| Nov 13 | `3758f0b5` | Corrected-LE bar and pixel plots |
| Nov 16 | `329f4d4c` | LE-cut overlays and target labels |
| Nov 16 | `57ce976a` | All-bar peak yields and initial optional bar-30 y correction |
| Nov 17 | `1b766c45` | Bar-30 LE/ECal-time plots and weighted sigma summary |
| Nov 17 | `f3d6d27d` | Endpoint-based y correction, position correlations, optional upper cuts |
| Nov 17 | `53aa3844` | Detector/layer position groups and control-bar auto-zoom |
| Nov 18 | `6a0e642e` | Per-layer CDet-minus-ECal-projection spectra |
| Nov 18 | `9449de66` | Detector CDet-x and native ECal-x spectra |
| Nov 18, committed with this document | — | Cut reports/cleaner plots, 7 ns/m model, pre-y dt/y plots, position dt window, fixed-bar-30 differences, 30 ns LE maxima, signal-gated layer combinations, bar-30 multiplicity |

The parent CDet_replay repository tracks these changes through its SBS-replay submodule reference. No commits in this lookback were found in the locally available SBS-offline or jlab-HPC histories; this recap does not attribute unrelated pre-existing working-tree edits to this week's implementation work.
