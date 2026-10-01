# CDet database timing eras and LH2 shift survey

## Purpose

This document records how the timestamp boundaries and fallback `shift_ns`
values in `DB/db_earm.cdet.dat` were established for the five-pass GEp data,
Run 3573 through Run 6088. It also reports the first direct test of the
additional CDet timing shift needed by representative LH2 runs relative to the
already-established cross-target shift in the same era.

The distinction is important:

- the **era fallback** is the mean `shift_ns` measured from the cross-target
  calibration runs in that stable acquisition period;
- an **LH2 increment** is an additional offset relative to that fallback and
  is permitted to depend on the target/trigger condition;
- the Ben-derived cross-target ECal slope is `p1 = 0.810602`, whereas the
  physically motivated LH2 value is `p1 = 0`.

Pixel offsets, time-walk constants, the TDC conversion, ECal `p0` and `delta`,
geometry, and pulse-quality definitions remain common across the experiment.

## Inputs and reproducibility

The era construction uses:

- `CDet_cross_target_shift_calibrations.tsv`: the accepted cross-target
  `shift_ns` fits;
- `plot_cross_target_shift_vs_run.py`: the checked-in four-group plot;
- `CDet_5pass_CODA_timestamps_updated.csv`: CODA timestamps extracted from
  `gep5_RUN.evio.0.0`;
- `CDet_runs4042_4043_CODA_timestamps.csv`: the explicitly recovered adjacent
  timestamps that fix the otherwise ambiguous Era 2/Era 3 boundary;
- `CDet_run6088_CODA_timestamps.csv`: the timestamp of the final five-pass run;
- `DB/db_earm.cdet.dat`: the timestamp-valid implementation.

The LH2 comparison uses the analyzer-backed histogram
`hCDetSelectedPairMeanCorrectedLE` in each
`CDet_runRUN_good_pulse_tdc/CDetGoodPulse_AllTDC.root`. Regenerate the summary
table and figures with:

```console
python3 plot_cdet_lh2_shift_relative_to_era.py
```

The script writes `CDet_lh2_shift_relative_to_era.tsv`, PNG, and PDF.

## Cross-target measurements and era means

Thirty-eight accepted cross-target runs divide naturally into four timing
groups. The fallback for each group is the unweighted mean of its fitted
run-specific shifts.

| Era | Cross-target calibration runs | ECal ADC-time window | Mean `shift_ns` | Run-to-run sample SD |
|---:|---|---|---:|---:|
| 1 | 3573, 3575 | -15 to 30 ns | -7.937335 ns | 0.020 ns |
| 2 | 3602 through 4006 (17 accepted runs) | 0 to 35 ns | -9.359456 ns | 0.093 ns |
| 3 | 4345 through 5423 (13 accepted runs) | 10 to 35 ns | -11.060945 ns | 0.333 ns |
| 4 | 5722 through 5795 (6 accepted runs) | 10 to 35 ns | 1.095437 ns | 0.087 ns |

![Cross-target shift versus run](CDet_cross_target_shift_vs_run.png)

The dashed horizontal segments are the four means used as database fallbacks.
The large changes between groups are much larger than the within-group scatter,
which is the empirical basis for a piecewise-constant era model.

## Timestamp boundaries used by the database

The analyzer database selects constants by event timestamp, not by a run-number
comparison. The run ranges below are therefore explanatory labels for the
timestamp-valid intervals.

| Era | Descriptive run interval | Database start timestamp | Boundary evidence |
|---:|---|---|---|
| 1 | 3573--3601 | 2025-05-09 10:09:25 (Run 3573) | First five-pass CDet run and first cross-target group |
| 2 | 3602--4042 | 2025-05-10 02:02:13 (Run 3602) | First run of the second stable cross-target group |
| 3 | 4043--5709 | 2025-05-24 07:09:25 (Run 4043) | Explicit adjacent CODA timestamps for Runs 4042 and 4043; transfers the third-group fallback into the uncalibrated gap before Run 4345 |
| 4 | 5710--6088 | 2025-08-07 01:26:21 (Run 5710) | First run of the fourth timing configuration; Run 6088 at 2025-08-25 03:15:23 closes the requested five-pass scope |

Run 4042 began at 2025-05-24 04:09:43 and Run 4043 began at
2025-05-24 07:09:25. Recovering both timestamps avoided placing the Era 3
validity boundary at a guessed time. The assignment of the common Era 3 LH2
shift to the non-cross-target interval beginning with Run 4043 is an explicit
operational interpolation; it is not a direct `shift_ns` measurement of Run
4043.

Every individually calibrated cross-target run overrides its era mean with its
own fitted `shift_ns` and uses `p1 = 0.810602`. Every non-cross-target interval
resets to `p1 = 0`, the common LH2 shift for that era, the LH2 ECal timing
interval, and the LH2 pairing switches. Runs 5711 and 6077 now use the common
Era 4 LH2 shift rather than retaining separate overrides.

## Representative LH2 survey and fallback-start test

Eight representative LH2 runs are now compared:

```text
3649  4346  4724  5295  5711  5727  5886  6077
```

Runs 5711 and 6077 originally had separately determined shifts of `3.422000
ns` and `2.913449 ns`. To test whether those starting values affected pair
selection, both runs were subsequently replayed with the same Era 4 fallback
used for Runs 5727 and 5886:

```text
p1       = 0
shift_ns = 1.095437 ns
```

The standard `CDet_run5711_good_pulse_tdc` and
`CDet_run6077_good_pulse_tdc` directories now contain products from these
fallback-start replays. The older independent shift values remain useful
historical results, but they are not the shifts applied to the ROOT files used
in the table below.

The comparison applies the same deliberately narrow diagnostic selection to
every run:

```text
3.0 < E_ECal < 4.5 GeV
-5 < t_ECal < 5 ns
(trajectory residual / 0.020 m)^2
  + ((t_ECal - pair mean + 26 ns) / 5 ns)^2 <= 1
```

The last condition is the established trajectory-time ellipse with its radius
reduced from 2 to 1. The selected corrected-pair distributions have RMS values
of about `3.2--3.6 ns`, rather than the `5.3--6.8 ns` widths in the superseded
broad comparison. Medians are retained as robust location estimators so the
result does not depend on the detailed Gaussian fit model.

During this study, the plotting macro was corrected so that the configured
ECal-time interval is applied when filling the selected-pair residual and
absolute-time histograms. Earlier broad figures did not apply that event-time
gate to those particular histograms and must not be compared numerically with
the results below.

![LH2 timing relative to cross-target fallbacks](CDet_lh2_shift_relative_to_era.png)

The top panel shows the selected corrected-pair spectra after replay with the
cross-target-era fallback. The lower-left panel compares that fallback with the
candidate LH2 shift required to move the measured pair-time median to 30 ns.
The lower-right panel shows the same result as the additional LH2 shift above
the cross-target fallback.

## Fallback-start results

All four Era 4 rows below were replayed with exactly the same cross-target
`shift_ns`. Runs 5711 and 6077 therefore provide the controlled test that was
missing from the earlier comparison.

The established calibration convention is

```text
candidate LH2 shift = cross-target fallback + (30 ns - pair-time location)

LH2 increment above era = candidate LH2 shift - cross-target fallback
```

The production calibration driver uses the Gaussian-core centroid of the
projected-half-bar timing distribution as the pair-time location. The table
below uses the median of the narrow detector-wide selected-pair distribution
as a robust diagnostic proxy. These values should therefore be confirmed with
the production projected-half-bar estimator before changing the database.

| Run | Era | Selected pairs | Cross-target fallback | Pair-time median | Correction to 30 ns | Candidate LH2 shift | Increment above era |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3649 | 2 | 5,790 | -9.359456 | 26.072027 | +3.927973 | -5.431483 | +3.927973 |
| 4346 | 3 | 7,089 | -11.060945 | 25.785221 | +4.214779 | -6.846166 | +4.214779 |
| 4724 | 3 | 12,487 | -11.060945 | 25.515637 | +4.484363 | -6.576582 | +4.484363 |
| 5295 | 3 | 8,741 | -11.060945 | 25.515312 | +4.484688 | -6.576257 | +4.484688 |
| 5711 | 4 | 5,278 | +1.095437 | 26.636672 | +3.363328 | +4.458765 | +3.363328 |
| 5727 | 4 | 7,413 | +1.095437 | 26.901685 | +3.098315 | +4.193752 | +3.098315 |
| 5886 | 4 | 11,590 | +1.095437 | 26.716835 | +3.283165 | +4.378602 | +3.283165 |
| 6077 | 4 | 10,313 | +1.095437 | 26.576543 | +3.423457 | +4.518894 | +3.423457 |

The Era 4 increments are `3.363`, `3.098`, `3.283`, and `3.423 ns`. Their mean
is `3.292 ns` and their sample standard deviation is `0.141 ns`. Thus all four
LH2 runs require a common additional shift of approximately `+3.3 ns` relative
to the Era 4 cross-target fallback according to this diagnostic. Runs 5711 and
6077 differ by only `0.060 ns`, despite having started from the same fallback;
their agreement is no longer a consequence of assigning them different input
shifts.

The three Era 3 increments are `4.215`, `4.484`, and `4.485 ns`, with mean
`4.395 ns` and sample standard deviation `0.156 ns`. Run 3649 gives an Era 2
increment of `3.928 ns`. The data therefore continue to favor an LH2
correction defined relative to each cross-target era rather than one absolute
`shift_ns` for the entire experiment.

## Provisional LH2 shift policy by era

For operational use, the run measurements within each era are averaged and
one common LH2 `shift_ns` is assigned to that era. This is intended to place
the hydrogen timing population consistently within the same trajectory-time
ellipse throughout the experiment.

| Era | Cross-target average | Mean measured LH2 increment | Provisional LH2 `shift_ns` | Basis |
|---:|---:|---:|---:|---|
| 1 | -7.937335 ns | not yet measured | **-5.431483 ns** | Temporarily use the Era 2 absolute LH2 shift |
| 2 | -9.359456 ns | +3.927973 ns | **-5.431483 ns** | Run 3649 |
| 3 | -11.060945 ns | +4.394610 ns | **-6.666335 ns** | Mean of Runs 4346, 4724, and 5295 |
| 4 | +1.095437 ns | +3.292066 ns | **+4.387503 ns** | Mean of Runs 5711, 5727, 5886, and 6077 |

The Era 1 value is deliberately an absolute-shift inheritance, not an inferred
Era 1 increment. Relative to the Era 1 cross-target average it corresponds to
`+2.505852 ns`, but that difference has not been measured with Era 1 hydrogen
data. It must be replaced if a future Era 1 LH2 analysis supports a different
value.

## Interpretation and next step

The controlled fallback-start comparison supports a simple operational model:
within a measured acquisition era, the LH2 runs require an approximately
common increment above the cross-target timing origin. For Era 4 that increment
is about `+3.3 ns`; the corresponding diagnostic candidate shift is about
`4.39 ns`. Era 1 provisionally inherits the Era 2 absolute LH2 shift pending a
direct measurement.

The authoritative database now implements the common LH2 shifts in the table
above for every non-cross-target interval. Individually fitted cross-target
runs retain their own shifts. The older Run 5711 and Run 6077 values of
`3.422000 ns` and `2.913449 ns` have been superseded in the database by the
common Era 4 value `4.387503 ns`.

The present calculation uses the narrow selected-pair median, whereas the
original calibration procedure used a local Gaussian fit to the
projected-half-bar core. A future validation should replay representative runs
with the installed era constants and verify that the common ellipse accepts
the hydrogen timing population consistently. The selected-pair population
should be recorded before and after any subsequent refinement so selection
changes remain visible.
