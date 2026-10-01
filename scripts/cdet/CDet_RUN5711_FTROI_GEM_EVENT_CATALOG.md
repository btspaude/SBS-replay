# Run 5711 FTROI/GEM event catalogue

This study turns the 358 Run 5711 events containing at least one reconstructed
front-tracker GEM track into a reproducible event-by-event review queue. It is
based on the 18 topology-aware replay files containing 261,010 events.

The generated catalogue is
[`CDet_run5711_FTROI_GEM_event_catalog.csv`](CDet_run5711_FTROI_GEM_event_catalog.csv).
It contains one row per GEM-track event and identifies the event number, global
and source-file-local entries, segment, chain tree, and exact source ROOT file.

## Cohorts

| Cohort | Events | Meaning |
|---|---:|---|
| GEM track and at least one FTROI CDet hypothesis | 117 | Primary cross-arm/elastic-closure validation sample |
| GEM track with no FTROI CDet hypothesis | 241 | CDet acceptance and candidate-formation control sample |
| **Total** | **358** | All events with `sbs.gemFT.track.ntrack > 0` |

All 358 events have valid HEEP output, calibrated CDet timing status 2, and a
unique event number. The 241 GEM-only events all have CDet ROI status 0: none
has a detector-local pair or exclusive single-layer candidate. They are
therefore not instances of candidates being lost during FTROI hypothesis
export.

The GEM multiplicity is:

| GEM tracks in event | Events |
|---:|---:|
| 1 | 345 |
| 2 | 11 |
| 3 | 2 |

The six segment populations are reasonably uniform: 57, 63, 57, 62, 57, and
62 catalogue events for segments 0 through 5.

## Review tiers

The `review_tier` column supplies an initial deterministic queue. Lower tier
numbers should be inspected first; this is triage, not a physics selection.

| Tier | Events | First event numbers | Reason |
|---|---:|---|---|
| `1_overlap_interior_best_z` | 5 | 38340, 83370, 140892, 683628, 757116 | Only overlap events whose globally best constrained-x association is not at a scan boundary |
| `2_overlap_single` | 3 | 473571, 540852, 712965 | Contains an exclusive single-layer FTROI hypothesis |
| `3_overlap_seam` | 36 | 30855, 44091, 68763, 78552, 80124 | Contains an opposite-side seam hypothesis |
| `4_overlap_multi_track` | 2 | 139083, 659493 | Multiple GEM tracks after earlier overlap tiers are removed |
| `5_overlap_standard` | 71 | 147, 627, 31185, 35553, 51615 | Remaining one-track overlap events |
| `6_gem_only_multi_track` | 10 | 14217, 230049, 313122, 351753, 401040 | Multiple GEM tracks and no CDet candidate |
| `7_gem_only_no_candidate` | 231 | 1287, 8331, 8898, 10419, 23160 | One GEM track and no CDet candidate |

The tiers are mutually exclusive and sum to 358. The free-form `review_tags`
column additionally records overlapping properties such as `seam`, `single`,
`multi_track`, and `boundary_best_z`.

Among the 117 overlap events, 40 contain at least one seam hypothesis, three
contain a single-layer hypothesis, and 112 have their globally smallest
constrained-x chi-square at the first or last target-z scan bin. The latter is
consistent with the known target-coordinate problem and must not be treated as
physical vertex reconstruction. At the catalogue-selected best association,
102 events pass topology-aware y compatibility and 15 do not. The selected
best topology is same-side pair in 97 events, seam pair in 17, and single layer
in three.

## Catalogue contents

Each row contains:

- stable identity: catalogue index, event number, global and local entries,
  chain tree, segment, and source filename;
- cohort, review tier, and review tags;
- GEM multiplicity plus the replay-selected GEM track's hit counts, fit
  quality, intercepts, slopes, and time;
- HEEP validity and elastic diagnostics (`dpe`, `dpp`, `dphi`, acoplanarity,
  ECal residuals, and ECal-CDet ADC time difference);
- ECal energy and position;
- CDet timing/ROI status and pulse, pair-candidate, and single-candidate
  counts;
- FTROI hypothesis counts by pair, same-side, seam, and single-layer topology;
- per-event y-compatible and total target associations;
- the globally minimum constrained-x association, including hypothesis,
  topology, target bin, target z, chi-square/NDF, y decision, and boundary flag.

`sbs.gemFT` describes the proton-arm front tracker, whereas CDet and FTROI
describe the electron arm. Their slopes are **not directly comparable**. Any
cross-arm ranking must use elastic closure quantities such as the HEEP
residuals or a reviewed two-arm kinematic construction.

## Reproducing the catalogue

Start in this directory with the intended Podd 1.7.12 environment and run:

```text
root [0] Catalog_CDet_FTROI_GEMEvents.C(
             "/Users/brash/CDet_replay/sbs/Rootfiles/FTROI_step5/rootfiles/gep5_replayed_5711_stream0_2_seg*_firstevent0_nevent100000*.root",
             "CDet_run5711_FTROI_GEM_event_catalog.csv");
```

The macro must report 18 input files, 358 rows, and an overlap/GEM-only split
of 117/241.

## Inspecting one event

`Inspect_CDet_FTROI_GEMEvent.C` accepts a Run 5711 event number and prints the
exact source file, global entry, every GEM track, HEEP diagnostics, CDet
candidate state, and every FTROI hypothesis with its best target association.
For the first tier-1 event:

```text
root [0] Inspect_CDet_FTROI_GEMEvent.C(38340)
```

Event 38340 contains one GEM track and 24 pair hypotheses. One of those
hypotheses has an interior best bin at target z = -0.1825 m with constrained-x
chi-square/NDF approximately 0.00686; many alternative hypotheses still select
scan boundaries. This illustrates why the event-level view is more useful than
retaining only one aggregate best value.

Two additional smoke-test examples are:

```text
root [0] Inspect_CDet_FTROI_GEMEvent.C(1287)    // GEM-only, ROI status 0
root [1] Inspect_CDet_FTROI_GEMEvent.C(139083)  // Three GEM tracks, one FTROI pair
```

The existing graphical CDet browser is coupled to the older calibration
analysis vectors and does not yet consume this FTROI catalogue. The catalogue
and textual inspector establish stable event identity and review ordering; a
dedicated FTROI/GEM graphical display can now be added without redefining the
sample.
