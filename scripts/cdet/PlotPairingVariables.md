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

An empty `analysis.input_directory` uses `OUT_DIR`. The input directory can be
overridden in the second argument:

```cpp
PlotPairingVariables("PlotPairingVariables.conf", "/path/to/Rootfiles");
```

The configuration selects the run, event limit, segment range, Layer-1 pixel,
cuts, histogram ranges, and output behavior. `analysis.layer1_pixel` is a
zero-based Layer-1 logical pixel ID; its containing 16-channel bar is derived
automatically. The supplied configuration uses `-1`, which selects all 16
pixels in `analysis.layer1_bar` (bar 30 by default). In that mode, the macro
accepts Layer-1 pixel IDs from `16*bar` through `16*bar+15`. Set
`analysis.layer1_pixel` to a pixel ID for a pixel-specific study. `analysis.min_segment` and
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

The existing upper row of the focused plots uses every stored pair whose
Layer-1 member matches the selected pixel and that passes the configured study
cuts. The lower row is an event-level best-pair view. For
each ECal-admitted event, the macro chooses the accepted pair with the smallest
`earm.cdet.pair.ecal_score`. If that branch is nonfinite, it uses the explicitly
recomputed ellipse radius squared as the ranking value. The focused plots choose the best accepted pair among candidates whose Layer-1
pixel is the requested pixel; the all-detector geometry plots choose the best accepted pair
anywhere in the detector. This produces one pair per event for the best-pair
row, while preserving the all-pair population above it.

This is a diagnostic ranking of pairs already stored in the ROOT file. It does
not rerun the replay assignment or create a new pair collection. The best-pair
selection is made after the configured member-ToT, optional layer-dt, and
pair-radius cuts.

The pair ellipse is enabled by default:

```ini
cuts.pair_radius_max: 2
```

The macro explicitly recomputes the normalized radius from
`earm.cdet.pair.trajectory_residual` and
`earm.cdet.pair.ecal_residual`, then accepts a pair when `R^2 <= 4`. The
replay also stores a diagnostic `earm.cdet.pair.ecal_score`; the explicit
calculation makes the timing and trajectory terms visible and checks the same
selection independently. This is the Run 6077 ECal trajectory-time ellipse from
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

The first canvas contains two rows of three spectra for the selected Layer-1
pixel. The upper row contains every accepted pair and the lower row contains the
best accepted pair per event:

1. Stored pair-mean corrected LE, `(tL1 + tL2)/2`.
2. Inter-layer timing, `tL2 - tL1`.
3. ECal-minus-pair timing residual, `tECal - (tL1 + tL2)/2`.

The lower row has the same three quantities after the best-pair ranking.

The second canvas has two rows. Each row overlays the corrected LE spectra for
the Layer-1 member and its Layer-2 partner. The upper row uses all accepted
pairs; the lower row uses the best pair per event for the selected Layer-1 pixel.
The macro reports the full-sample standard deviation and its ROOT moment-based
error; it does not perform a Gaussian fit.

Timing annotations and terminal summaries are labeled in ns. Position-residual
annotations are labeled in m, and out-of-plane angle annotations are labeled in
degrees. This avoids applying the timing unit to the geometry panels.
The `N` and standard-deviation boxes are intentionally compact. Each canvas
also includes a small two-line cut box showing the stored-pair ToT, optional
layer-dt, pair-radius, ECal adctime, and ECal energy selections.

The current display range for the ECal pair residual is `-40` to `0 ns`. This
is a display range, while the ellipse uses the residual centered near `-26 ns`.
The focused-pixel LE spectra use `plots.bar_le_min_ns` and
`plots.bar_le_max_ns`. The focused x1-versus-x1-x2 panels use
`plots.bar_x1_min_m` and `plots.bar_x1_max_m`. These ranges are separate from
the detector-wide LE and x1 ranges because their peak locations can move with
pixel or bar position. Layer timing, ECal-pair timing, pair-position residuals,
CDet delta-x residuals, and angle plots retain the detector-wide ranges because their peak
locations should be comparable across the detector.

The older `plots.le_min_ns` and `plots.le_max_ns` keys remain accepted for
compatibility with the positional interface. When using the configuration
interface, use the `plots.bar_le_*` keys for the focused LE panels.

## Pair geometry canvases

The focused geometry canvas has two rows of three panels. The upper row uses
all accepted pairs from the selected Layer-1 pixel; the lower row uses the best
accepted pair per event for that pixel. Each row contains:

1. A two-dimensional plot with `x1 - x2` on the horizontal axis and `x1` on
   the vertical axis.
2. The CDet pair position residual,
   \[
   \langle x\rangle_{pair} - x_{ECal}\frac{\langle z_{pair}\rangle}{6.144\,\mathrm{m}}.
   \]
3. A rough out-of-plane angle in transport coordinates.

### Out-of-plane angle derivation

The out-of-plane direction is the vertical transport-coordinate direction
measured by the CDet's 5 mm paddles. In this convention that is transport
`x`, with `z` pointing into the detector plane. The macro therefore uses the
pulse `x` and `z` coordinates. For each accepted pair, it constructs four
points in the transport `x-z` plane:

```text
point 0: (z0, x0) = (0, 0)                         assumed target origin
point 1: (z1, x1) = (pulse.z[i1], pulse.x[i1])     CDet Layer 1
point 2: (z2, x2) = (pulse.z[i2], pulse.x[i2])     CDet Layer 2
point 3: (z3, x3) = (6.144 m, earm.ecal.x)         ECal
```

The value `6.144 m` is the ECal distance used by the existing CDet projection
macros. The ECal x coordinate is treated as a measured point at that z. The
origin point is an analysis assumption for this rough diagnostic; it is not a
new reconstructed vertex.

The macro assumes that these four points are described approximately by a
straight line through the origin,

```text
x(z) = m z,
```

where `m` is the vertical slope. For a single point the least-squares residual
would be `xi - m zi`. The macro chooses `m` to minimize the sum of squared
residuals for all three measured points (the origin has zero residual for any
finite m):

\[
\chi^2(m) = \sum_{i=1}^{3}(x_i-mz_i)^2.
\]

Differentiating with respect to `m` and setting the result to zero gives

\[
\frac{d\chi^2}{dm} = -2\sum_{i=1}^{3}z_i(x_i-mz_i)=0,
\]

so the fitted slope is

\[
m = \frac{\sum_{i=1}^{3}z_i x_i}{\sum_{i=1}^{3}z_i^2}
  = \frac{z_1x_1+z_2x_2+6.144\,x_{ECal}}
         {z_1^2+z_2^2+6.144^2}.
\]

The physical angle of a line with slope `m` relative to the transport z axis
is

\[
\theta_x=\arctan(m).
\]

The macro converts this angle from radians to degrees:

```cpp
angle_degrees = std::atan(m) * 180.0 / \pi;
```

This calculation is repeated independently for every accepted pair in the
selected-bar geometry canvas and for every accepted pair in the all-detector
geometry canvas. A pair contributes only when both pulse x/z coordinates and
the ECal x coordinate are finite and the denominator
`z1*z1 + z2*z2 + 6.144*6.144` is positive. The pair’s two CDet points are not
averaged for this angle; both measurements enter the fit separately, along
with the ECal point.

The displayed angle is `theta_x` in degrees. The current display range is
`-40` to `40 degrees` with `0.1 degree` bins. This is a rough geometric diagnostic,
not a replacement for a full track or optics reconstruction. It does not use
uncertainties or perform a weighted fit, does not account for magnetic
transport or detector alignment, and does not fit a free vertex. A nonzero
beam-spot or vertex offset, curved transport, or significant x calibration
error can therefore bias this angle. The result should be interpreted as the
best straight-line slope through the assumed-origin four-point construction,
not as a precision scattering-angle measurement.

The standard deviation of the second panel is printed as a rough CDet x
position-resolution estimate in millimeters.

A second, otherwise identical two-row geometry canvas is filled before the
selected-bar restriction and therefore uses all accepted pairs in the detector
in its upper row and the best accepted pair per event in its lower row. Its
output names end in `_geometry_all`.

The macro also writes two-row `_geometry_residual` and `_geometry_residual_all` copies.
Their middle panels use the CDet separation residual, defined directly as

`x_residual = Δx_CDet − (xECal / zECal) × Δz_CDet`

with `Δx_CDet = x1 − x2` and `Δz_CDet = z1 − z2`. This is the ECal-guided
trajectory residual used by the ellipse, with the same sign convention as the
macro. Its standard deviation is printed as a separate Layer-1 x-position
estimate.

## Pre-ellipse diagnostic

The `_ellipse` canvas has two rows. The upper row shows the all-detector
candidate pair population before the ECal trajectory-time radius cut. The
lower row shows the best pre-ellipse candidate per event, ranked by the same
ECal score. Its horizontal axis is the trajectory residual `rₓ`, and its
vertical axis is `Δt_pair = tECal − <tCDet>pair`. The event-level ECal cuts,
pair validity, member-ToT cut, and optional inter-layer timing cut have already
been applied, but the pair-radius cut has not.

The configured ellipse is drawn on top of the candidate population. With the
Run 6077 settings it is centered at `(rₓ, Δt_pair) = (0 m, −26 ns)` and has
radius-2 semiaxes of `0.040 m` and `10 ns`. This makes it possible to see
which candidate pairs are removed by the spatial-plus-timing ellipse.

## Bar-width canvas

The fourth canvas has two rows. The upper row shows the standard deviation of
`tECal - tCDet,member` versus bar number for all accepted pairs; the lower row
shows the same quantity for the best pair per event. Both rows show the two
detector layers. Bars with fewer than `plots.min_entries_per_bar` accepted
members are omitted from the graph. The companion CSV retains every bar,
including low-statistics bars,
with its entries, mean, standard deviation, error, underflow, overflow, and
status.

## Progress and output

The terminal prints progress every 1,000 processed events, using the total
number of entries in the chained run files. It reports the cumulative number
of events containing at least one pair that passes the macro’s checks and
configured cuts.

With `output.save_plots: 1`, the macro saves eight thesis-ready canvases as
both PDF and PNG, plus the bar-width CSV. Files are written under
`output.directory` with run, containing-bar, and selected-pixel names when a
pixel is selected. Plot saving is disabled when `output.save_plots: 0`.
