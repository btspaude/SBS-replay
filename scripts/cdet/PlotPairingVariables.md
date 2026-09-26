# `PlotPairingVariables.C`

`PlotPairingVariables.C` is a pulse-and-pair analysis macro for the refined
CDet ROOT trees. It reads the newer `earm.cdet.pulse.*` and
`earm.cdet.pair.*` branches in one event loop. It does not use the legacy
`hit.*` or `hits.*` branches, and it does not rerun the pulse matching or pair
assignment.

## Running the macro

From ROOT, load the macro and use the TEnv-style configuration file:

```cpp
.L PlotPairingVariables.C+
PlotPairingVariables("PlotPairingVariables.conf");
```

The input directory can be overridden in the second argument:

```cpp
PlotPairingVariables("PlotPairingVariables.conf", "/path/to/Rootfiles");
```

The configuration selects the run, event limit, segment range, Layer-1 bar,
cuts, histogram ranges, and output behavior. `analysis.min_segment` and
`analysis.max_segment` set an inclusive segment range. Set both to `-1` to use
all available segments. The dataset helper groups rollover parts and chains
all matching segment files.

## Event and pair selection

The event must pass both ECal event-level requirements before any pair is
filled:

```text
-10 ns <= earm.ecal.adctime <= 4 ns
3.0 GeV <= earm.ecal.e <= 4.5 GeV
```

For each stored pair, the macro validates the pulse indices, Layer-1/Layer-2
pixel identities, finite timing values, and the pulse-to-pair array sizes. The
configured member ToT interval is inclusive for both pulses. The current
Run-6077 configuration uses `8 <= ToT <= 35 ns`.

The pair ellipse is enabled by default:

```ini
cuts.pair_radius_max: 2
```

The replay stores its normalized squared radius in `earm.cdet.pair.ecal_score`.
The macro therefore accepts a pair when `pair.ecal_score <= 4`. This is the
Run 6077 ECal trajectory-time ellipse from
`CDet_RUN6077_PAIR_SELECTION_STUDY.md`:

\[
R^2 = \left(\frac{r_x}{0.020\,\mathrm{m}}\right)^2
    + \left(\frac{\Delta t_{pair}+26\,\mathrm{ns}}{5\,\mathrm{ns}}\right)^2
\]

where

\[
r_x = (x_{L2}-x_{L1}) - \frac{x_{ECal}}{z_{ECal}}(z_{L2}-z_{L1}).
\]

The optional `cuts.layer_dt_max_ns` and `cuts.pair_radius_max` values can be
set to `-1` to disable those additional macro-level cuts. They do not undo
selections already applied upstream when the ROOT file was produced.

## Timing canvases

The first canvas contains three spectra for the selected Layer-1 bar:

1. Stored pair-mean corrected LE, `(tL1 + tL2)/2`.
2. Inter-layer timing, `tL2 - tL1`.
3. ECal-minus-pair timing residual, `tECal - (tL1 + tL2)/2`.

The second canvas overlays the corrected LE spectra for the Layer-1 member and
its Layer-2 partner. The macro reports the full-sample standard deviation and
its ROOT moment-based error; it does not perform a Gaussian fit.

The current display range for the ECal pair residual is `-40` to `0 ns`. This
is a display range, while the ellipse uses the residual centered near `-26 ns`.

## Pair geometry canvases

The focused geometry canvas has three panels, all using accepted pairs from
the selected Layer-1 bar:

1. A two-dimensional plot with `x1 - x2` on the horizontal axis and `x1` on
   the vertical axis.
2. The CDet pair position residual,
   \[
   \langle x\rangle_{pair} - x_{ECal}\frac{\langle z_{pair}\rangle}{6.144\,\mathrm{m}}.
   \]
3. A rough out-of-plane angle in transport coordinates.

For the angle diagnostic, the macro uses transport `y` and `z` and assumes the
trajectory begins at `(y,z)=(0,0)`. It fits a line constrained through the
origin and the three measured points: CDet layer 1, CDet layer 2, and ECal.
For the points `(zi, yi)`, the fitted slope is

\[
m = \frac{\sum_i z_i y_i}{\sum_i z_i^2},
\qquad \theta_y = \arctan(m).
\]

The displayed angle is `1000*theta_y` in mrad. The current display range is
`-60` to `80 mrad` with `5 mrad` bins. This is a rough geometric diagnostic,
not a replacement for a full track or optics reconstruction.

The standard deviation of the second panel is printed as a rough CDet x
position-resolution estimate in millimeters.

A second, otherwise identical geometry canvas is filled before the selected
bar restriction and therefore uses accepted pairs from all detector bars. Its
output names end in `_geometry_all`.

## Bar-width canvas

The fourth canvas shows the standard deviation of
`tECal - tCDet,member` versus bar number for both detector layers. Bars with
fewer than `plots.min_entries_per_bar` accepted members are omitted from the
graph. The companion CSV retains every bar, including low-statistics bars,
with its entries, mean, standard deviation, error, underflow, overflow, and
status.

## Progress and output

The terminal prints progress every 1,000 processed events, using the total
number of entries in the chained run files. It reports the cumulative number
of events containing at least one pair that passes the macro’s checks and
configured cuts.

With `output.save_plots: 1`, the macro saves five thesis-ready canvases as
both PDF and PNG, plus the bar-width CSV. Files are written under
`output.directory` with run and selected-bar names. Plot saving is disabled
when `output.save_plots: 0`.
