# CDet APS/DNP 2026 talk

This folder contains the self-contained first draft of Benjamin Spaude's
APS/DNP presentation, **The Coordinate Detector in GEp-V**.

## Build

From this directory, run:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error CDet_APS_DNP_2026.tex
```

The generated presentation is `CDet_APS_DNP_2026.pdf`.

All figures used by the Beamer source are stored in `assets/`; the build does
not depend on files elsewhere in the repository.

## Current scope

The draft contains 11 slides for a 10-minute talk followed by 2 minutes of
questions. Comments in the TeX source give suggested pacing. The ROI material
is intentionally described as validation and diagnostic work: CDet hypotheses
do not yet constrain production GEM pattern recognition.

