#!/usr/bin/env python3
"""Compare LH2 timing shifts with their cross-target era fallbacks."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

import cdet_matplotlib_style  # noqa: F401


HERE = Path(__file__).resolve().parent
HISTOGRAM = "hCDetSelectedPairMeanCorrectedLE"
RESIDUAL_HISTOGRAM = "hCDetSelectedPairMeanECalResidual"
TARGET_PAIR_TIME_NS = 30.0
OUTPUT_STEM = HERE / "CDet_lh2_shift_relative_to_era"
OUTPUT_TABLE = OUTPUT_STEM.with_suffix(".tsv")

# The replay shift is the value actually applied when each ROOT file was
# produced. Runs 5711 and 6077 were deliberately replayed from the same Era-4
# cross-target fallback as Runs 5727 and 5886 for a like-for-like comparison.
# Their older independently obtained shifts are retained as historical context.
RUNS = {
    3649: {"era": 2, "era_shift": -9.359456, "replay_shift": -9.359456,
           "historical_shift": None},
    4346: {"era": 3, "era_shift": -11.060945, "replay_shift": -11.060945,
           "historical_shift": None},
    4724: {"era": 3, "era_shift": -11.060945, "replay_shift": -11.060945,
           "historical_shift": None},
    5295: {"era": 3, "era_shift": -11.060945, "replay_shift": -11.060945,
           "historical_shift": None},
    5711: {"era": 4, "era_shift": 1.095437, "replay_shift": 1.095437,
           "historical_shift": 3.422000},
    5727: {"era": 4, "era_shift": 1.095437, "replay_shift": 1.095437,
           "historical_shift": None},
    5886: {"era": 4, "era_shift": 1.095437, "replay_shift": 1.095437,
           "historical_shift": None},
    6077: {"era": 4, "era_shift": 1.095437, "replay_shift": 1.095437,
           "historical_shift": 2.913449},
}

ERA_COLORS = {2: "#008837", 3: "#c51b7d", 4: "#0571b0"}


def histogram_summary(histogram) -> tuple[float, float, float]:
    probability = np.array([0.5], dtype="d")
    quantile = np.zeros(1, dtype="d")
    histogram.GetQuantiles(1, quantile, probability)
    return float(histogram.GetMean()), float(quantile[0]), float(histogram.GetRMS())


def read_histogram(run: int) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    path = HERE / f"CDet_run{run}_good_pulse_tdc" / "CDetGoodPulse_AllTDC.root"
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    histogram = root_file.Get(HISTOGRAM)
    residual_histogram = root_file.Get(RESIDUAL_HISTOGRAM)
    if not histogram or not residual_histogram:
        raise RuntimeError(f"Required selected-pair histograms are missing from {path}")

    bins = histogram.GetNbinsX()
    centers = np.array(
        [histogram.GetBinCenter(index) for index in range(1, bins + 1)]
    )
    counts = np.array(
        [histogram.GetBinContent(index) for index in range(1, bins + 1)]
    )
    mean, median, rms = histogram_summary(histogram)
    residual_mean, residual_median, residual_rms = histogram_summary(
        residual_histogram
    )
    summary = {
        "entries": float(histogram.GetEntries()),
        "mean": mean,
        "median": median,
        "rms": rms,
        "residual_mean": residual_mean,
        "residual_median": residual_median,
        "residual_rms": residual_rms,
    }
    root_file.Close()
    return centers, counts, summary


def main() -> None:
    spectra: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    results: dict[int, dict[str, float]] = {}
    for run, settings in RUNS.items():
        centers, counts, summary = read_histogram(run)
        spectra[run] = (centers, counts)
        results[run] = {**settings, **summary}

    for result in results.values():
        correction = TARGET_PAIR_TIME_NS - result["median"]
        result["target_correction"] = correction
        result["candidate_shift"] = result["replay_shift"] + correction
        result["increment_above_era"] = (
            result["candidate_shift"] - result["era_shift"]
        )

    with OUTPUT_TABLE.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "run",
            "era",
            "entries",
            "era_shift_ns",
            "replay_shift_ns",
            "historical_independent_shift_ns",
            "pair_time_mean_ns",
            "pair_time_median_ns",
            "pair_time_rms_ns",
            "correction_to_30_ns",
            "candidate_shift_ns",
            "increment_above_era_ns",
            "ecal_minus_pair_mean_ns",
            "ecal_minus_pair_median_ns",
            "ecal_minus_pair_rms_ns",
        ]
        writer = csv.DictWriter(
            stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for run, result in results.items():
            writer.writerow(
                {
                    "run": run,
                    "era": int(result["era"]),
                    "entries": int(result["entries"]),
                    "era_shift_ns": f'{result["era_shift"]:.6f}',
                    "replay_shift_ns": f'{result["replay_shift"]:.6f}',
                    "historical_independent_shift_ns": (
                        "" if result["historical_shift"] is None
                        else f'{result["historical_shift"]:.6f}'
                    ),
                    "pair_time_mean_ns": f'{result["mean"]:.6f}',
                    "pair_time_median_ns": f'{result["median"]:.6f}',
                    "pair_time_rms_ns": f'{result["rms"]:.6f}',
                    "correction_to_30_ns": f'{result["target_correction"]:.6f}',
                    "candidate_shift_ns": f'{result["candidate_shift"]:.6f}',
                    "increment_above_era_ns":
                        f'{result["increment_above_era"]:.6f}',
                    "ecal_minus_pair_mean_ns": f'{result["residual_mean"]:.6f}',
                    "ecal_minus_pair_median_ns": f'{result["residual_median"]:.6f}',
                    "ecal_minus_pair_rms_ns": f'{result["residual_rms"]:.6f}',
                }
            )

    fig = plt.figure(figsize=(14.5, 9.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.2, 1.0))
    spectra_axis = fig.add_subplot(grid[0, :])
    shift_axis = fig.add_subplot(grid[1, 0])
    increment_axis = fig.add_subplot(grid[1, 1])

    linestyles = ["-", "--", "-.", ":", "-", "--", ":", "-."]
    for (run, result), linestyle in zip(results.items(), linestyles):
        centers, counts = spectra[run]
        normalized = counts / counts.sum()
        color = ERA_COLORS[int(result["era"])]
        width = 2.5 if run in (5711, 6077) else 1.8
        spectra_axis.step(
            centers,
            normalized,
            where="mid",
            color=color,
            linewidth=width,
            linestyle=linestyle,
            label=(
                f'Run {run}, Era {int(result["era"])}: '
                f'median {result["median"]:.2f} ns'
            ),
        )
    spectra_axis.axvline(30.0, color="#555555", linestyle=":", linewidth=1.5)
    spectra_axis.text(30.25, 0.97, "30 ns reference", rotation=90,
                      transform=spectra_axis.get_xaxis_transform(), va="top",
                      color="#555555")
    spectra_axis.set_xlim(8.0, 48.0)
    spectra_axis.set_xlabel(r"Absolute corrected selected-pair mean time (ns)")
    spectra_axis.set_ylabel("Normalized pairs per 1 ns")
    spectra_axis.set_title(
        r"Narrow diagnostic: $-5<t_{ECal}<5$ ns and trajectory-time radius 1"
    )
    spectra_axis.legend(ncol=2, fontsize=10.2, framealpha=0.95)
    spectra_axis.grid(alpha=0.25)

    runs = np.array(list(results))
    era_shifts = np.array([results[run]["era_shift"] for run in runs])
    candidate_shifts = np.array([results[run]["candidate_shift"] for run in runs])
    increments = np.array([results[run]["increment_above_era"] for run in runs])
    colors = [ERA_COLORS[int(results[run]["era"])] for run in runs]

    for run, baseline, candidate, color in zip(
        runs, era_shifts, candidate_shifts, colors
    ):
        shift_axis.plot([run, run], [baseline, candidate], color=color,
                        linewidth=2.3, alpha=0.8)
    shift_axis.scatter(runs, era_shifts, marker="s", s=60, color=colors,
                       edgecolor="white", linewidth=0.7,
                       label="Cross-target era fallback")
    shift_axis.scatter(runs, candidate_shifts, marker="o", s=75, color=colors,
                       edgecolor="black", linewidth=0.6,
                       label="LH2 candidate from 30 ns target")
    shift_axis.set_xlabel("Run number")
    shift_axis.set_ylabel(r"$shift_{ns}$ (ns)")
    shift_axis.set_title("Cross-target fallback and candidate LH2 shift")
    shift_axis.grid(alpha=0.25)
    shift_axis.legend(fontsize=9.5)

    increment_axis.bar([str(run) for run in runs], increments, color=colors,
                       edgecolor="white", linewidth=0.8)
    increment_axis.axhline(0.0, color="black", linewidth=0.8)
    for index, value in enumerate(increments):
        increment_axis.text(index, value + 0.08, f"{value:.2f}",
                            ha="center", va="bottom", fontsize=10)
    increment_axis.set_xlabel("LH2 survey run")
    increment_axis.set_ylabel("Additional shift above cross-target fallback (ns)")
    increment_axis.set_title("LH2 increment required to place the median at 30 ns")
    increment_axis.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "CDet LH2 timing relative to cross-target era fallbacks",
        fontsize=17,
        fontweight="semibold",
    )
    metadata = {
        "Title": "CDet LH2 timing relative to cross-target era fallbacks",
        "Subject": "Fallback-start spectra and candidate LH2 timing increments",
        "Creator": "plot_cdet_lh2_shift_relative_to_era.py",
    }
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=200, metadata=metadata)
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"), metadata=metadata)


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    main()
