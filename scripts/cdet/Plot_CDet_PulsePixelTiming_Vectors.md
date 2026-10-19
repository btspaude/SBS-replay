# Load-once pulse-vector analysis

Implemented 2026-10-19 from the current `Plot_CDet_PulsePixelTiming.C` in SBS-replay
`88fea7e6a136d287ed135c774b676850f1dbf7df` (working-file contents, not an older
historical version). The original macro and configuration are unchanged. This
copy retains the existing fits, plots and physics definitions while separating
file reading from analysis. No replay, database, calibration or pairing changes
are made.

## Start in ROOT

From the directory selected by `gocdetscripts`:

```cpp
.L Plot_CDet_PulsePixelTiming_Vectors.C+
const char* cfg = "CDet_run6077_pulse_pixel_timing_vectors.conf";

// Read all fields, initially using a small sample. OUT_DIR supplies the directory.
Load_CDet_PulsePixelTiming_Vectors(cfg, nullptr, 1000);

// These calls only read the cached vectors and reread the configuration.
Plot_CDet_PulseTiming_Vectors(cfg);
Plot_CDet_PulseLEToT_Vectors(cfg);
Plot_CDet_PulseOccupancy_Vectors(cfg);
Plot_CDet_PulseProjectionResponse_Vectors(cfg);

// Or produce the full suite. The existing cache is reused.
Plot_CDet_PulsePixelTiming_Vectors(cfg);
```

All six functions return `bool`; check for `true`. The full-suite entry point
loads on its first call if there is no complete cache. Individual routines
require a successful load first. Compilation and the first call take extra
time; they are separate from the printed loader/analysis times.

```cpp
// Edit analysis.le_max, fit settings, display ranges, etc. in cfg, then:
Plot_CDet_PulseTiming_Vectors(cfg); // no file reread or implicit recalibration

// Explicit display overrides, retaining the current cuts:
Plot_CDet_PulseTiming_Vectors(cfg, 19, 0.5, -90, 0,
                             0.005, -0.10, 0.10, -90, 0, "fine_timing");
Plot_CDet_PulseLEToT_Vectors(cfg, 19, 1, -20, 100, 1, 0, 100, "le_tot");

// Explicit reload replaces, rather than appends to, the cache.
Load_CDet_PulsePixelTiming_Vectors(cfg, nullptr, 5000);
Plot_CDet_PulseOccupancy_Vectors(cfg);

// Respect analysis.events from the configuration (-1 means all entries).
Load_CDet_PulsePixelTiming_Vectors(cfg);

// Explicit directory override:
Load_CDet_PulsePixelTiming_Vectors(cfg, "/Users/spaubt/CDet_replay/sbs/data", 1000);
```

Use an intentional small event limit until the measured storage is acceptable.
The supplied configuration inherits the original 10,000,000-event limit; it is
not automatically a small test. Loading all fields is intentional. An allocation
failure clears the cache and reports failure; no fields are dropped or narrowed.

`Plot_CDet_PulsePixelTiming_Vectors(6077)` selects
`CDet_run6077_pulse_pixel_timing_vectors.conf`. The full configuration-path
signature keeps the original argument order through the final output tag.
The `(cfg, "tag", bar, directory, events)` and `(run, "tag", ...)` conveniences
are available too. The numeric-run overload also retains the full explicit display signature.

## Function layout and ownership

- `Load_CDet_PulsePixelTiming_Vectors` discovers files through `CDetRunDataset`,
  inspects the pulse/ECal schema in every selected file, then fills all vectors
  in one `TTreeReader` pass. It applies no physics selection.
- `Plot_CDet_PulseTiming_Vectors` creates timing histograms, loops over
  `[event][pulse]`, then fits and draws seven timing canvases.
- `Plot_CDet_PulseLEToT_Vectors` fills the same final pulse selection into six
  LE and LE-versus-ToT canvases.
- `Plot_CDet_PulseOccupancy_Vectors` uses separate event/pixel masks and draws
  one occupancy canvas.
- `Plot_CDet_PulseProjectionResponse_Vectors` constructs the independent ECal
  denominator and draws three canvases plus a per-pixel CSV.
- `Plot_CDet_PulsePixelTiming_Vectors` invokes the four routines in sequence.

The functions and histogram/cache variables live in `CDetPulseVectors`.
Global `using` declarations expose the six entry-point names at the ROOT prompt.
The original macro can be loaded alongside the copy. Inspect, for example,
`CDetPulseVectors::hPixelDt[0]` versus the original `::hPixelDt[0]`.
ROOT object names have the `pv_` prefix. Output goes to
`PulsePixelTimingVectors/<tag>/` under the current directory, with the same tag
rules as the original. Only the rerun routine's canvases/histograms are replaced.
Repeated calls overwrite that routine's PDFs rather than accumulating counts.
The original `PulsePixelTiming/` products remain separate.

The cache is session-local RAM, not a persistent disk format. Treat all cache
containers as read-only from the prompt. Named `vPulse*` references are const;
the underlying typed maps exist to support discovered fields. `ClearCache()`
releases cached vector storage; plots remain until rerun or closed. A failed
reload leaves no usable old or partial cache, although old displayed plots can
remain visible. `cacheReady` and `loadGeneration` identify successful loads.

## Dataset settings versus analysis settings

Every routine rereads the supplied `.conf`. Dataset identity includes
`analysis.run_number`, `analysis.events`, `analysis.min_segment`,
`analysis.max_segment` and the resolved configuration input directory. The new
`analysis.input_directory` key is empty by default, meaning `OUT_DIR`.
An explicit loader directory overrides that key; an explicit loader event limit
other than `-2` overrides `analysis.events` (`-1` means all). Those effective
overrides remain attached to the cache. A full-suite call with differing
explicit dataset overrides is rejected.

Changing any dataset key requires an explicit load, even when that key was
previously overridden. An `OUT_DIR` change is also detected when it supplies
the configuration directory. Paths are expanded and made absolute; symlink
aliases are not canonicalized. This deliberately errs toward asking for a
reload rather than silently accepting a different dataset. File replacement
on disk under unchanged filenames is **not** automatically detected by plot
calls: reload explicitly after a replay. `loadedFiles`, `loadedSchema`,
`loadedDirectory`, `loadedEventLimit`, and `vEntry` record the loaded selection.
File grouping/rollover selection follows the existing `CDetRunDataset` policy.

Cuts, fit settings, bar choices, output tags, geometry settings and display
axes do not invalidate the cache. Projection rereads its date-selected geometry
DB. This does not change stored pulse projections or replay flags. Its independent
reference calculation uses that supplied geometry, exactly as in the original.
Relative geometry paths resolve relative to the configuration file.

New display keys in the copied configuration are:

| Keys (prefix `vectors.display.`) | Defaults |
| --- | --- |
| `dt_bin_width_ns` | 1 ns |
| `dx_bin_width_m`, `dx_min_m`, `dx_max_m` | 0.005, -0.10, 0.10 m |
| `dt_vs_dx_min_ns`, `dt_vs_dx_max_ns` (optional) | Same as 1D dt bounds |
| `le_bin_width_ns`, `le_min_ns`, `le_max_ns` | 1, -20, 100 ns |
| `tot_bin_width_ns`, `tot_min_ns`, `tot_max_ns` | 1, 0, 100 ns |

`display.ecal_cdet_dt_min/max` retain their existing role. `NaN` argument
defaults select configured values; supplied display arguments override them.
Histogram bin counts use `ceil((max-min)/width)`, as in the original. Display
ranges do not impose selection cuts; entries can appear in flow bins.

## Complete branch mapping

The actual Run 6077 input `cdet_6077_stream0_2_seg0_0.root` contains the documented
**33 pulse fields: 29 arrays and four scalars**. All are `Double_t`, including
IDs, flags and counters. This matches `SBSCDet::DefineVariables()` and the
[field reference](../../../../../docs/btspaude/CDET_PULSE_BRANCHES.md).
`Ndata.earm.cdet.pulse.*` count metadata are read through `TTreeReaderArray`;
they are not additional physics fields beginning with `earm.cdet.pulse`.

Every name below is in `CDetPulseVectors`. Pulse branch suffixes omit
`earm.cdet.pulse.`. All array views are
`const std::vector<std::vector<Double_t>>&`, indexed `[event][pulse]`.
Values, NaNs, flags and pulse order are copied without cuts or calibration.

| Branch suffix | Named vector | Units and correction state |
| --- | --- | --- |
| `pmtnum` | `vPulseID` | Logical pixel ID, stored as double; 0–2687 |
| `index` | `vPulseIndex` | Accepted-pulse index within its channel, not event-wide index |
| `le_index` | `vPulseLEIndex` | Original LE decoder slot within channel |
| `te_index` | `vPulseTEIndex` | Original TE decoder slot within channel |
| `row` | `vPulseRow` | Geometry row label |
| `col` | `vPulseCol` | Geometry column label |
| `layer` | `vPulseLayer` | Stored geometry layer label; may be zero for both physical layers in older builder |
| `x` | `vPulseX` | Geometry x, m |
| `y` | `vPulseY` | Geometry half-bar y, m; not measured interaction y |
| `z` | `vPulseZ` | Geometry z, m |
| `tdc_le` | `vPulseLe` | Legacy TDC units, low-level offset/gain already applied |
| `tdc_te` | `vPulseTe` | Legacy TDC units, low-level offset/gain already applied |
| `tdc_tot` | `vPulseToT` | Legacy TE minus LE, same units |
| `tdc_le_raw` | `vPulseLeRaw` | Decoder raw LE before low-level calibration |
| `tdc_te_raw` | `vPulseTeRaw` | Decoder raw TE before low-level calibration |
| `tdc_tot_raw` | `vPulseToTRaw` | Decoder raw TE minus LE |
| `tdc_le_corr` | `vPulseLeCorr` | Fully corrected LE, ns |
| `tdc_te_corr` | `vPulseTeCorr` | Same additive timing correction as LE, ns |
| `tdc_tot_ns` | `vPulseToTNs` | Converted ToT, ns; common timing offsets cancel |
| `ecal_residual` | `vPulseDt` | ECal ADC time minus corrected CDet LE, ns |
| `calib_valid` | `vPulseCalibValid` | Stored timing-calibration validity flag |
| `x_corr` | `vPulseXCorr` | Replay scaled/aligned CDet x, m |
| `ecal_x_proj` | `vPulseECalXProj` | Target-origin ECal ray projected to pulse z, m |
| `ecal_y_proj` | `vPulseECalYProj` | Same projection in y, m |
| `ecal_x_residual` | `vPulseDx` | Corrected CDet x minus projected ECal x minus replay residual offset, m |
| `ecal_y_residual` | `vPulseDy` | CDet y minus projected ECal y minus replay residual offset, m |
| `broad_quality_pass` | `vPulseBroadQualityPass` | Stored LE/ToT-window flag; retained but not required by these plots |
| `ecal_eligible` | `vPulseECalEligible` | Stored cluster/time/acceptance/projection decision |
| `spatial_pass` | `vPulseSpatialPass` | Stored x **and** y residual decision |

The four scalar views are `const std::vector<Double_t>&`, indexed `[event]`:

| Branch suffix | Named vector | Meaning |
| --- | --- | --- |
| `n` | `vPulseN` | Number of complete stored pulses |
| `n_le_only` | `vPulseNLEOnly` | Ignored extra/unmatched leading edges |
| `n_te_only` | `vPulseNTEOnly` | Trailing edges without an open LE |
| `n_invalid_pair` | `vPulseNInvalidPair` | Rejected nonpositive raw/calibrated ToT combinations |

| Other branch | Vector | Type / units |
| --- | --- | --- |
| `earm.ecal.x` | `vECalX` | Double_t, m |
| `earm.ecal.y` | `vECalY` | Double_t, m |
| `earm.ecal.adctime` | `vECalTime` | Double_t, selected cluster ADC time in ns |
| `earm.ecal.nclus` | `vECalNclus` | Double_t, cluster count |
| `earm.ecal.e` | `vECalEnergy` | Double_t, selected cluster energy in GeV |
| TChain entry index | `vEntry` | Long64_t, original chain index |

Additional pulse-prefixed branches are dynamically discovered. Arrays go into
`pulseArrays<T>[fullBranchName][event][pulse]`, scalars into
`pulseScalars<T>[fullBranchName][event]`, retaining the inspected primitive type.
Supported types are ROOT's Double_t, Float_t, Int_t, UInt_t, Long64_t, ULong64_t,
Short_t, UShort_t, Char_t, UChar_t, Bool_t, Long_t and ULong_t. The documented
baseline currently requires Double_t, so a changed baseline type fails clearly
rather than converting it. Unsupported types, multi-leaf branches, fixed-length
arrays, missing baseline/ECal fields, per-file schema differences, and array
lengths differing from `pulse.n` also fail explicitly. Extend the loader with a
corresponding typed reader if the replay schema evolves; nothing is silently
skipped. A successful cache has exactly one entry in every field per event,
including empty arrays for events with zero pulses.

## Selections and an event walkthrough

Timing and LE/ToT call the same ordinary `GoodBeforeEllipse` and `GoodEllipse`
functions inside direct event/pulse loops. The former requires finite ID, dx,
dt, converted LE and positive finite ToT, stored `calib_valid`, `ecal_eligible`
and `spatial_pass` greater than 0.5, strict configured ECal energy bounds,
inclusive configured LE bounds, and optional inclusive ToT bounds. The ellipse
uses stored dx/dt and an inclusive normalized squared-radius boundary. Radius
`-1` disables it; otherwise it must be positive. Neither helper rewrites data.

For event `e`, the loader saves `vEntry[e]`, ECal scalars and every pulse field.
Pulse `p` retains its identity across all arrays. Suppose its ID is 304, its
legacy LE is 1000, stored ToT is 12 ns, dx is 0.01 m and stored dt is -26 ns.
The LE test uses `0.01 * vPulseLe[e][p] = 10 ns`; the timing histogram fills
`vPulseDt[e][p] = -26 ns` directly. It does not multiply corrected times by 0.01
or apply offsets again. If the replay flags, ECal energy and configured cuts
pass, bar 19 pixel 0 receives one entry. A second passing pulse on pixel 304
adds another timing entry. Occupancy records that pixel only once for event e.

Raw occupancy accepts all finite physical pixel IDs. Its ECal subset adds the
strict energy interval, stored eligibility and finite dx; its x-matched subset
adds `abs(dx) <= occupancy_dx_max_m`. It adds no LE/ToT, calibration, spatial,
y, ellipse or configurable time gate. The baseline setting changes the plotted
comparison, not the definition of the three cached occupancy histograms.

Projection's reference requires strict configured ECal energy and ADC-time
bounds, positive finite cluster count and finite x/y. It does not require any
pulse. A zero-pulse event can therefore add an opportunity but no response.
Complete response means any pulse on that pixel; selected response means at
least one final-selected pulse. Both count once per pixel/event. The separate
`hECalX` histogram counts events with any final-selected pulse; it is not forced
through the reference selection. These details match the original macro.

Changing `analysis.ecal_time_min/max` changes the independent reference only;
it cannot rewrite stored `ecal_eligible`. Likewise, changing projection geometry
or x windows does not recompute stored `spatial_pass`. The layer view is retained
as stored; physical layer identification here uses ID/1344, as in the original.

## Validation and limits

Validated locally with ROOT 6.40 on 2026-10-19. The original macro and copy ran
in the same ROOT process. The real sample was the first 1,000 Run 6077 entries
with segment range 0..0 (the helper selected the base file and its rollover
part). Only the first file's entries were needed for this sample. Every stored
Double_t value was checked bit-for-bit against a second independent reader,
including 1,070,420 NaN values across the fields.

| Real-sample measurement | Result |
| --- | --- |
| Cached events / complete pulses | 1,000 / 488,252 |
| Mean / maximum pulse multiplicity | 488.252 / 563 |
| Zero-pulse events in real sample | 0; tested separately with synthetic data |
| Estimated vector allocation | 108.767 MiB |
| Load time | 0.277 s |
| Original full analysis including I/O, plots/export | 0.707 s |
| Cached full analysis, first / repeat | 0.371 / 0.394 s |
| Baseline study / control / detector selected pulses | 4 / 32 / 4,413 |
| Baseline independent ECal reference events | 284 |
| Edited-cut study / control / detector selected pulses | 1 / 10 / 975 |
| Edited-cut independent ECal reference events | 247 |

These are single-run wall times after compilation, with warm filesystem effects
and other work on the machine; they do not establish a general speedup. The
allocation estimate includes vector capacities and nested-vector objects, not
allocator metadata, maps, ROOT buffers, histograms or process overhead. At this
sample's multiplicity, linear extrapolation is roughly **106 GiB per million
events**, before those overheads. A full-run cache was not attempted. The source
configuration's large event limit remains unchanged; the examples deliberately
load 1,000 entries first.

The regression driver compares all entries, bin contents, bin errors and axis
edges in **114 histograms with ellipse off**, and **117 with ellipse on**,
including all 1D/2D underflow and overflow cells. Both fraction graphs and their
y confidence intervals also match. Repeated full-suite and timing-only calls
pass the same comparisons and do not increase the canvas count. Timing-only
calls leave the other routine objects intact. Editing LE maximum to 45 ns,
energy minimum to 3.2 GeV, enabling the ToT window and setting ellipse radius 2
refills from the same cache generation and matches a new original-macro pass.
Changing the event limit is rejected until explicit reload; reloading 500 events
replaces the 1,000-event cache. An intentionally missing input directory fails
and clears the cache. Fit/drawing helpers are retained from the original; fit
parameter equality and full-run statistical performance were not separate tests.

The synthetic fixture adds a zero-pulse event, repeated pulses on the same
pixel, strict ECal cut endpoints, inclusive LE/ToT endpoints, populated flow
bins, NaNs and contradictory stored flags to expose accidental extra cuts.
It adds Float_t arrays and ULong64_t scalar values above 2^53 to check discovery
and precision. It is split into a base file and rollover part to exercise a
chain transition. Its reference includes all seven qualifying ECal events,
including the zero-pulse event; duplicate pulses contribute only once to pixel
occupancy. A second segment with a missing extra field is rejected for schema
mismatch. The fixture and real-data driver both compare the full original/copy
suite and exercise reload/failure behavior.

To reproduce, use absolute paths and a scratch output directory. Run the two
drivers in **separate fresh ROOT sessions**; each includes both macros:

```cpp
// Session 1, from gocdetscripts:
gSystem->SetBuildDir("/private/tmp/cdet_vectors_validation_build", true);
.L Validate_CDet_PulsePixelTiming_Vectors.C+
Validate_CDet_PulsePixelTiming_Vectors(
  "/Users/spaubt/CDet_replay/git-repo/sbs_devel/SBS-replay/scripts/cdet/CDet_run6077_pulse_pixel_timing_vectors.conf",
  "/Users/spaubt/CDet_replay/sbs/data",
  "/Users/spaubt/CDet_replay/git-repo/sbs_devel/SBS-replay/DB/db_earm.cdet.dat",
  "/private/tmp/cdet_vectors_validation", 1000);

// Session 2:
gSystem->SetBuildDir("/private/tmp/cdet_vectors_fixture_build", true);
.L Validate_CDet_PulseVectors_Fixture.C+
Validate_CDet_PulseVectors_Fixture(
  "/Users/spaubt/CDet_replay/git-repo/sbs_devel/SBS-replay/scripts/cdet/CDet_run6077_pulse_pixel_timing_vectors.conf",
  "/Users/spaubt/CDet_replay/git-repo/sbs_devel/SBS-replay/DB/db_earm.cdet.dat",
  "/private/tmp/cdet_vectors_fixture");
```

The first driver is intended for a Run 6077 sample with at least the requested
entry count and uses segment 0. Both drivers write scratch configurations and
plots, including original-macro outputs, only under their output directory.
The synthetic driver creates/replaces its named fixture ROOT files there.
Normal analysis only reads replay ROOT files.

Original-file SHA256 checks before and after the implementation matched:

```text
Plot_CDet_PulsePixelTiming.C
c5f282950db4af41c886c1067311c85110ee8163c67ad7355b8b6648b96fe55e
CDet_run6077_pulse_pixel_timing.conf
a23e8240c0a1f8eb8f0a7837b1dcf8dfdd6e12b965ac64279dc9c8fbfec4a48d
```
