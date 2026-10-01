#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "CDetRunDataset.h"

namespace {

struct HydrogenPair {
  int pulse1;
  int pulse2;
  double meanTime;
  double dt;
  double dx;
  double dy;
  double ecalResidual;
  double trajectoryResidual;
  double ecalScore;
};

struct HydrogenSingle {
  int pulse;
  int layer;
  double time;
  double ecalResidual;
  double xResidual;
  double score;
};

void StyleComparison(TH1D &oldHist, TH1D &newHist) {
  oldHist.SetLineColor(kBlue + 1);
  oldHist.SetLineWidth(3);
  newHist.SetLineColor(kRed + 1);
  newHist.SetLineStyle(2);
  newHist.SetLineWidth(3);
  oldHist.SetStats(false);
  newHist.SetStats(false);
}

void DrawComparison(TH1D &oldHist, TH1D &newHist, TLegend &legend) {
  StyleComparison(oldHist, newHist);
  const double maximum = std::max(oldHist.GetMaximum(), newHist.GetMaximum());
  oldHist.SetMaximum(maximum > 0.0 ? 1.12 * maximum : 1.0);
  oldHist.Draw("hist");
  newHist.Draw("hist same");
  legend.DrawClone();
}

void DrawTimingPeak(TH1D &hist, const char *fitName, bool narrowCore = false) {
  hist.SetLineColor(kBlack);
  hist.SetLineWidth(2);
  hist.SetStats(true);
  hist.Draw("hist");
  if (hist.GetEntries() < 20)
    return;
  double fitMin = -55.0;
  double fitMax = -10.0;
  int peakBin = hist.FindBin(-40.0);
  for (int bin = peakBin + 1; bin <= hist.FindBin(-15.0); ++bin)
    if (hist.GetBinContent(bin) > hist.GetBinContent(peakBin))
      peakBin = bin;
  const double peak = hist.GetBinCenter(peakBin);
  if (narrowCore) {
    fitMin = peak - 5.0;
    fitMax = peak + 5.0;
  }
  const double background =
      0.5 * (hist.GetBinContent(hist.FindBin(fitMin)) +
             hist.GetBinContent(hist.FindBin(fitMax)));
  TF1 fit(fitName, narrowCore ? "gaus" : "gaus(0)+pol1(3)", fitMin,
          fitMax);
  if (narrowCore)
    fit.SetParameters(hist.GetBinContent(peakBin), peak, 3.0);
  else
    fit.SetParameters(std::max(1.0, hist.GetBinContent(peakBin) - background),
                      peak, 3.0, background, 0.0);
  fit.SetLineColor(kRed + 1);
  hist.Fit(&fit, "RQN");
  fit.DrawCopy("same");
}

template <class T>
bool SameValue(T oldValue, T newValue, double tolerance) {
  return std::isfinite(double(oldValue)) && std::isfinite(double(newValue)) &&
         std::fabs(double(oldValue) - double(newValue)) <= tolerance;
}

} // namespace

// Hydrogen-policy parity study. The old path independently reconstructs the
// validated pair/seam and exclusive single-layer selections from pulse.*; the
// new path reads pair_candidate.*, single_candidate.*, and roi.status only.
void Compare_CDet_HydrogenCandidateBranches(
    int runNumber = 6077, const char *inputDirectory = nullptr,
    const char *outputDirectory = nullptr) {
  constexpr double tolerance = 1.0e-10;
  constexpr double dtMax = 15.0;
  constexpr double dxMax = 0.15;
  constexpr double dyMax = 0.08;
  constexpr double oppositeDYCenter = 0.51;
  constexpr double oppositeDYTolerance = 0.08;
  constexpr double oppositeProjectedYMax = 0.17;
  constexpr double pairXScale = 0.020;
  constexpr double pairTimeCenter = -26.0;
  constexpr double pairTimeScale = 5.0;
  constexpr double pairRadius = 2.0;
  constexpr double singleXScale = 0.040;
  constexpr double singleTimeCenter = -26.0;
  constexpr double singleTimeScale = 5.0;
  constexpr double singleRadius = 2.0;
  constexpr double yResidualOffset = 0.10;
  constexpr double ecalTimeMin = -10.0;
  constexpr double ecalTimeMax = 10.0;
  constexpr double ecalEnergyMin = 3.0;
  constexpr double ecalEnergyMax = 4.5;

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[hydrogen candidate comparison] inputDirectory or OUT_DIR "
                 "is required.\n";
    return;
  }
  const TString defaultOutput =
      TString::Format("CDet_run%d_hydrogen_candidate_comparison", runNumber);
  const char *output = outputDirectory && outputDirectory[0]
                           ? outputDirectory
                           : defaultOutput.Data();

  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[hydrogen candidate comparison] No events found for run "
              << runNumber << " in " << input << ".\n";
    return;
  }
  const char *required[] = {
      "earm.cdet.pulse.pmtnum", "earm.cdet.pulse.tdc_le_corr",
      "earm.cdet.pulse.tdc_tot_ns", "earm.ecal.adctime", "earm.ecal.e",
      "earm.cdet.pulse.ecal_residual", "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass",
      "earm.cdet.pulse.ecal_eligible", "earm.cdet.pulse.spatial_pass",
      "earm.cdet.pulse.x_corr", "earm.cdet.pulse.y",
      "earm.cdet.pulse.ecal_x_proj", "earm.cdet.pulse.ecal_y_proj",
      "earm.cdet.pair_candidate.pulse_index_l1",
      "earm.cdet.pair_candidate.pulse_index_l2",
      "earm.cdet.pair_candidate.time_mean",
      "earm.cdet.pair_candidate.dt", "earm.cdet.pair_candidate.dx",
      "earm.cdet.pair_candidate.dy",
      "earm.cdet.pair_candidate.ecal_residual",
      "earm.cdet.pair_candidate.trajectory_residual",
      "earm.cdet.pair_candidate.ecal_score",
      "earm.cdet.single_candidate.pulse_index",
      "earm.cdet.single_candidate.layer", "earm.cdet.single_candidate.time",
      "earm.cdet.single_candidate.ecal_residual",
      "earm.cdet.single_candidate.x_residual",
      "earm.cdet.single_candidate.score", "earm.cdet.roi.status"};
  chain.LoadTree(0);
  for (const char *branch : required) {
    if (!chain.GetBranch(branch)) {
      std::cerr << "[hydrogen candidate comparison] Missing branch " << branch
                << ".\n";
      return;
    }
  }

  gSystem->mkdir(output, true);
  gStyle->SetOptStat(0);

  TH1D oldPairN("hOldPairN", "Pair hypotheses per event;hypotheses;events",
                31, -0.5, 30.5);
  TH1D newPairN("hNewPairN", "Pair hypotheses per event;hypotheses;events",
                31, -0.5, 30.5);
  TH1D oldPairTime("hOldPairTime", "Pair mean corrected LE;time (ns);pairs",
                   120, 0, 60);
  TH1D newPairTime("hNewPairTime", "Pair mean corrected LE;time (ns);pairs",
                   120, 0, 60);
  TH1D oldPairResidual("hOldPairResidual",
                       "Pair ECal timing residual;t_{ECal}-<t_{CDet}> (ns);pairs",
                       120, -40, 20);
  TH1D newPairResidual("hNewPairResidual",
                       "Pair ECal timing residual;t_{ECal}-<t_{CDet}> (ns);pairs",
                       120, -40, 20);
  TH1D oldTrajectory("hOldPairTrajectory",
                     "Pair trajectory residual;residual (m);pairs", 120,
                     -0.12, 0.12);
  TH1D newTrajectory("hNewPairTrajectory",
                     "Pair trajectory residual;residual (m);pairs", 120,
                     -0.12, 0.12);
  TH1D oldSingleN("hOldSingleN",
                  "Exclusive single-layer candidates per event;candidates;events",
                  16, -0.5, 15.5);
  TH1D newSingleN("hNewSingleN",
                  "Exclusive single-layer candidates per event;candidates;events",
                  16, -0.5, 15.5);
  TH1D oldSingleResidual("hOldSingleResidual",
      "Single-layer ECal timing residual;t_{ECal}-t_{CDet} (ns);candidates",
      120, -40, 20);
  TH1D newSingleResidual("hNewSingleResidual",
      "Single-layer ECal timing residual;t_{ECal}-t_{CDet} (ns);candidates",
      120, -40, 20);
  TH2D oldPair2D("hOldHydrogenPair2D",
      "Old reconstruction;trajectory residual (m);t_{ECal}-<t_{CDet}> (ns)",
      120, -0.12, 0.12, 120, -40, 20);
  TH2D newPair2D("hNewHydrogenPair2D",
      "SBSCDet pair_candidate.*;trajectory residual (m);t_{ECal}-<t_{CDet}> (ns)",
      120, -0.12, 0.12, 120, -40, 20);
  TH2D oldSingle2D("hOldHydrogenSingle2D",
      "Old reconstruction;single x residual (m);t_{ECal}-t_{CDet} (ns)",
      120, -0.16, 0.16, 120, -40, 20);
  TH2D newSingle2D("hNewHydrogenSingle2D",
      "SBSCDet single_candidate.*;single x residual (m);t_{ECal}-t_{CDet} (ns)",
      120, -0.16, 0.16, 120, -40, 20);
  TH1D oldDetectorECalCut("hOldDetectorECalCut",
      "Old: ECal energy cut, events with an accepted CDet pair;"
      "t_{ECal}-t_{CDet,corr} (ns);Selected hits", 180, -60, 30);
  TH1D newDetectorECalCut("hNewDetectorECalCut",
      "Native: ECal energy cut, events with an accepted CDet pair;"
      "t_{ECal}-t_{CDet,corr} (ns);Selected hits", 180, -60, 30);
  TH1D oldDetectorProjected("hOldDetectorProjected",
      "Old: projection-matched pulses, events with an accepted CDet pair;"
      "t_{ECal}-t_{CDet,corr} (ns);Trajectory-matched hits", 180, -60, 30);
  TH1D newDetectorProjected("hNewDetectorProjected",
      "Native: projection-matched pulses, events with an accepted CDet pair;"
      "t_{ECal}-t_{CDet,corr} (ns);Trajectory-matched hits", 180, -60, 30);
  TH1D oldAmalgamatedPairResidual("hOldAmalgamatedPairResidual",
      "Old: trajectory-time selected pair mean, all CDet;"
      "t_{ECal}-<t_{CDet,corr}>_{pair} (ns);Selected pairs", 180, -60, 30);
  TH1D newAmalgamatedPairResidual("hNewAmalgamatedPairResidual",
      "Native: trajectory-time selected pair mean, all CDet;"
      "t_{ECal}-<t_{CDet,corr}>_{pair} (ns);Selected pairs", 180, -60, 30);
  TH2D oldAmalgamatedPairVsToT("hOldAmalgamatedPairResidualVsMeanToT",
      "Old: trajectory-time selected pair mean, all CDet;"
      "<CDet ToT>_{pair} (ns);t_{ECal}-<t_{CDet,corr}>_{pair} (ns)",
      40, 0, 40, 180, -60, 30);
  TH2D newAmalgamatedPairVsToT("hNewAmalgamatedPairResidualVsMeanToT",
      "Native: trajectory-time selected pair mean, all CDet;"
      "<CDet ToT>_{pair} (ns);t_{ECal}-<t_{CDet,corr}>_{pair} (ns)",
      40, 0, 40, 180, -60, 30);

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pmt(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> tot(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderArray<Double_t> ecalResidual(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> calib(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> quality(reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> eligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatial(reader, "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<Double_t> x(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> y(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<Double_t> projX(reader, "earm.cdet.pulse.ecal_x_proj");
  TTreeReaderArray<Double_t> projY(reader, "earm.cdet.pulse.ecal_y_proj");
  TTreeReaderArray<Double_t> pairP1(reader, "earm.cdet.pair_candidate.pulse_index_l1");
  TTreeReaderArray<Double_t> pairP2(reader, "earm.cdet.pair_candidate.pulse_index_l2");
  TTreeReaderArray<Double_t> pairMean(reader, "earm.cdet.pair_candidate.time_mean");
  TTreeReaderArray<Double_t> pairDT(reader, "earm.cdet.pair_candidate.dt");
  TTreeReaderArray<Double_t> pairDX(reader, "earm.cdet.pair_candidate.dx");
  TTreeReaderArray<Double_t> pairDY(reader, "earm.cdet.pair_candidate.dy");
  TTreeReaderArray<Double_t> pairER(reader, "earm.cdet.pair_candidate.ecal_residual");
  TTreeReaderArray<Double_t> pairTR(reader, "earm.cdet.pair_candidate.trajectory_residual");
  TTreeReaderArray<Double_t> pairES(reader, "earm.cdet.pair_candidate.ecal_score");
  TTreeReaderArray<Double_t> singlePulse(reader, "earm.cdet.single_candidate.pulse_index");
  TTreeReaderArray<Double_t> singleLayer(reader, "earm.cdet.single_candidate.layer");
  TTreeReaderArray<Double_t> singleTime(reader, "earm.cdet.single_candidate.time");
  TTreeReaderArray<Double_t> singleER(reader, "earm.cdet.single_candidate.ecal_residual");
  TTreeReaderArray<Double_t> singleXR(reader, "earm.cdet.single_candidate.x_residual");
  TTreeReaderArray<Double_t> singleScore(reader, "earm.cdet.single_candidate.score");
  TTreeReaderValue<Double_t> roiStatus(reader, "earm.cdet.roi.status");

  Long64_t events = 0, oldPairs = 0, newPairs = 0, oldSingles = 0, newSingles = 0;
  Long64_t pairCountMismatch = 0, pairIdentityMismatch = 0, pairValueMismatch = 0;
  Long64_t singleCountMismatch = 0, singleIdentityMismatch = 0;
  Long64_t singleValueMismatch = 0, roiMismatch = 0, malformed = 0;

  while (reader.Next()) {
    ++events;
    const size_t n = std::min({pmt.GetSize(), le.GetSize(), tot.GetSize(), ecalResidual.GetSize(),
        calib.GetSize(), quality.GetSize(), eligible.GetSize(), spatial.GetSize(),
        x.GetSize(), y.GetSize(), projX.GetSize(), projY.GetSize()});
    if (n != pmt.GetSize() || n != le.GetSize() || n != tot.GetSize() || n != ecalResidual.GetSize() ||
        n != calib.GetSize() || n != quality.GetSize() || n != eligible.GetSize() ||
        n != spatial.GetSize() || n != x.GetSize() || n != y.GetSize() ||
        n != projX.GetSize() || n != projY.GetSize())
      ++malformed;

    std::vector<int> layer1, layer2;
    for (size_t i = 0; i < n; ++i) {
      if (!(calib[i] > 0.5 && quality[i] > 0.5 && eligible[i] > 0.5 &&
            spatial[i] > 0.5))
        continue;
      const int pixel = int(std::lround(pmt[i]));
      if (pixel >= 0 && pixel < 1344) layer1.push_back(int(i));
      if (pixel >= 1344 && pixel < 2688) layer2.push_back(int(i));
    }

    std::vector<HydrogenPair> oldPair;
    for (int p1 : layer1) for (int p2 : layer2) {
      const double dt = le[p2] - le[p1];
      const double dx = x[p2] - x[p1];
      const double dy = y[p2] - y[p1];
      const bool sameSide = std::fabs(dy) <= dyMax;
      const double projectedY = 0.5 * (projY[p1] + projY[p2]) + yResidualOffset;
      const bool oppositeSide =
          std::fabs(std::fabs(dy) - oppositeDYCenter) <= oppositeDYTolerance &&
          std::fabs(projectedY) <= oppositeProjectedYMax;
      if (std::fabs(dt) > dtMax || std::fabs(dx) > dxMax ||
          (!sameSide && !oppositeSide))
        continue;
      const double tr = dx - (projX[p2] - projX[p1]);
      const double er = 0.5 * (ecalResidual[p1] + ecalResidual[p2]);
      const double es = (tr / pairXScale) * (tr / pairXScale) +
          ((er - pairTimeCenter) / pairTimeScale) *
          ((er - pairTimeCenter) / pairTimeScale);
      if (es > pairRadius * pairRadius)
        continue;
      oldPair.push_back({p1, p2, 0.5 * (le[p1] + le[p2]), dt, dx, dy, er, tr, es});
    }
    std::stable_sort(oldPair.begin(), oldPair.end(),
        [](const HydrogenPair &a, const HydrogenPair &b) {
          return a.ecalScore < b.ecalScore;
        });

    std::vector<HydrogenSingle> oldSingle;
    if (layer1.empty() != layer2.empty()) {
      const std::vector<int> &populated = layer1.empty() ? layer2 : layer1;
      for (int pulse : populated) {
        const int layer = int(std::lround(pmt[pulse])) / 1344;
        const double xr = x[pulse] - projX[pulse];
        const double er = ecalResidual[pulse];
        const double score = (xr / singleXScale) * (xr / singleXScale) +
            ((er - singleTimeCenter) / singleTimeScale) *
            ((er - singleTimeCenter) / singleTimeScale);
        if (score <= singleRadius * singleRadius)
          oldSingle.push_back({pulse, layer, le[pulse], er, xr, score});
      }
    }

    const size_t np = std::min({pairP1.GetSize(), pairP2.GetSize(), pairMean.GetSize(),
        pairDT.GetSize(), pairDX.GetSize(), pairDY.GetSize(), pairER.GetSize(),
        pairTR.GetSize(), pairES.GetSize()});
    const size_t ns = std::min({singlePulse.GetSize(), singleLayer.GetSize(),
        singleTime.GetSize(), singleER.GetSize(), singleXR.GetSize(),
        singleScore.GetSize()});
    oldPairN.Fill(oldPair.size()); newPairN.Fill(np);
    oldSingleN.Fill(oldSingle.size()); newSingleN.Fill(ns);
    oldPairs += oldPair.size(); newPairs += np;
    oldSingles += oldSingle.size(); newSingles += ns;
    if (oldPair.size() != np) ++pairCountMismatch;
    if (oldSingle.size() != ns) ++singleCountMismatch;

    const bool passesECalEnergy = std::isfinite(*ecalEnergy) &&
        *ecalEnergy >= ecalEnergyMin && *ecalEnergy <= ecalEnergyMax;
    const bool passesECalTime = std::isfinite(*ecalTime) &&
        *ecalTime >= ecalTimeMin && *ecalTime <= ecalTimeMax;
    if (passesECalEnergy) {
      if (!oldPair.empty()) {
        for (size_t i = 0; i < n; ++i) {
          if (!(calib[i] > 0.5 && eligible[i] > 0.5) ||
              !std::isfinite(ecalResidual[i]))
            continue;
          oldDetectorECalCut.Fill(ecalResidual[i]);
          if (spatial[i] > 0.5)
            oldDetectorProjected.Fill(ecalResidual[i]);
        }
      }
      if (np > 0) {
        for (size_t i = 0; i < n; ++i) {
          if (!(calib[i] > 0.5 && eligible[i] > 0.5) ||
              !std::isfinite(ecalResidual[i]))
            continue;
          newDetectorECalCut.Fill(ecalResidual[i]);
          if (spatial[i] > 0.5)
            newDetectorProjected.Fill(ecalResidual[i]);
        }
      }
    }

    std::map<std::pair<int,int>, HydrogenPair> pairMap;
    for (const auto &pair : oldPair) {
      pairMap[{pair.pulse1, pair.pulse2}] = pair;
      oldPairTime.Fill(pair.meanTime); oldPairResidual.Fill(pair.ecalResidual);
      oldTrajectory.Fill(pair.trajectoryResidual);
      oldPair2D.Fill(pair.trajectoryResidual, pair.ecalResidual);
      if (passesECalEnergy && passesECalTime) {
        oldAmalgamatedPairResidual.Fill(pair.ecalResidual);
        const double meanToT = 0.5 * (tot[pair.pulse1] + tot[pair.pulse2]);
        if (std::isfinite(meanToT))
          oldAmalgamatedPairVsToT.Fill(meanToT, pair.ecalResidual);
      }
    }
    bool pairIdentityBad = false;
    for (size_t i = 0; i < np; ++i) {
      const std::pair<int,int> key = {int(std::lround(pairP1[i])),
                                      int(std::lround(pairP2[i]))};
      newPairTime.Fill(pairMean[i]); newPairResidual.Fill(pairER[i]);
      newTrajectory.Fill(pairTR[i]); newPair2D.Fill(pairTR[i], pairER[i]);
      if (passesECalEnergy && passesECalTime) {
        newAmalgamatedPairResidual.Fill(pairER[i]);
        const int source1 = int(std::lround(pairP1[i]));
        const int source2 = int(std::lround(pairP2[i]));
        if (source1 >= 0 && source2 >= 0 && source1 < int(n) &&
            source2 < int(n)) {
          const double meanToT = 0.5 * (tot[source1] + tot[source2]);
          if (std::isfinite(meanToT))
            newAmalgamatedPairVsToT.Fill(meanToT, pairER[i]);
        }
      }
      auto found = pairMap.find(key);
      if (found == pairMap.end()) { pairIdentityBad = true; continue; }
      const auto &old = found->second;
      const double oldValues[] = {old.meanTime, old.dt, old.dx, old.dy,
                                  old.ecalResidual, old.trajectoryResidual,
                                  old.ecalScore};
      const double newValues[] = {pairMean[i], pairDT[i], pairDX[i], pairDY[i],
                                  pairER[i], pairTR[i], pairES[i]};
      for (int k = 0; k < 7; ++k)
        if (!SameValue(oldValues[k], newValues[k], tolerance))
          ++pairValueMismatch;
      pairMap.erase(found);
    }
    if (!pairMap.empty()) pairIdentityBad = true;
    if (pairIdentityBad) ++pairIdentityMismatch;

    std::map<int, HydrogenSingle> singleMap;
    for (const auto &single : oldSingle) {
      singleMap[single.pulse] = single;
      oldSingleResidual.Fill(single.ecalResidual);
      oldSingle2D.Fill(single.xResidual, single.ecalResidual);
    }
    bool singleIdentityBad = false;
    for (size_t i = 0; i < ns; ++i) {
      const int pulse = int(std::lround(singlePulse[i]));
      newSingleResidual.Fill(singleER[i]); newSingle2D.Fill(singleXR[i], singleER[i]);
      auto found = singleMap.find(pulse);
      if (found == singleMap.end()) { singleIdentityBad = true; continue; }
      const auto &old = found->second;
      if (old.layer != int(std::lround(singleLayer[i]))) ++singleValueMismatch;
      const double oldValues[] = {old.time, old.ecalResidual, old.xResidual, old.score};
      const double newValues[] = {singleTime[i], singleER[i], singleXR[i], singleScore[i]};
      for (int k = 0; k < 4; ++k)
        if (!SameValue(oldValues[k], newValues[k], tolerance))
          ++singleValueMismatch;
      singleMap.erase(found);
    }
    if (!singleMap.empty()) singleIdentityBad = true;
    if (singleIdentityBad) ++singleIdentityMismatch;

    int expectedStatus = 0;
    if (!oldPair.empty()) expectedStatus = 1;
    else if (!oldSingle.empty()) expectedStatus = oldSingle.front().layer + 2;
    if (int(std::lround(*roiStatus)) != expectedStatus) ++roiMismatch;
  }

  StyleComparison(oldPairResidual, newPairResidual);
  TLegend legend(0.52, 0.74, 0.88, 0.88);
  legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.AddEntry(&oldPairResidual, "Old: reconstructed from pulse.*", "l");
  legend.AddEntry(&newPairResidual, "New: SBSCDet candidates", "l");

  TCanvas spectra("cHydrogenCandidateComparison",
                  "Hydrogen candidate comparison", 1800, 1050);
  spectra.Divide(3, 2);
  spectra.cd(1); DrawComparison(oldPairN, newPairN, legend);
  spectra.cd(2); DrawComparison(oldPairTime, newPairTime, legend);
  spectra.cd(3); DrawComparison(oldPairResidual, newPairResidual, legend);
  spectra.cd(4); DrawComparison(oldTrajectory, newTrajectory, legend);
  spectra.cd(5); DrawComparison(oldSingleN, newSingleN, legend);
  spectra.cd(6); DrawComparison(oldSingleResidual, newSingleResidual, legend);
  spectra.SaveAs(Form("%s/CDet_Run%d_HydrogenCandidateComparison.pdf", output,
                      runNumber));
  spectra.SaveAs(Form("%s/CDet_Run%d_HydrogenCandidateComparison.png", output,
                      runNumber));

  TCanvas maps("cHydrogenCandidateComparison2D",
               "Hydrogen candidate comparison 2D", 1700, 1000);
  maps.Divide(2, 2);
  maps.cd(1); oldPair2D.Draw("colz"); maps.cd(2); newPair2D.Draw("colz");
  maps.cd(3); oldSingle2D.Draw("colz"); maps.cd(4); newSingle2D.Draw("colz");
  maps.SaveAs(Form("%s/CDet_Run%d_HydrogenCandidateComparison_2D.pdf", output,
                   runNumber));
  maps.SaveAs(Form("%s/CDet_Run%d_HydrogenCandidateComparison_2D.png", output,
                   runNumber));

  gStyle->SetOptStat(1110);
  TCanvas amalgamated("cHydrogenCandidateAmalgamatedComparison",
      "Hydrogen detector-amalgamated comparison", 1800, 1050);
  amalgamated.Divide(4, 2);
  amalgamated.cd(1); DrawTimingPeak(oldDetectorECalCut, "fOldDetectorECalCut");
  amalgamated.cd(2); DrawTimingPeak(oldDetectorProjected, "fOldDetectorProjected");
  amalgamated.cd(3); DrawTimingPeak(oldAmalgamatedPairResidual, "fOldPairMean", true);
  amalgamated.cd(4); oldAmalgamatedPairVsToT.SetStats(false);
  oldAmalgamatedPairVsToT.Draw("colz");
  amalgamated.cd(5); DrawTimingPeak(newDetectorECalCut, "fNewDetectorECalCut");
  amalgamated.cd(6); DrawTimingPeak(newDetectorProjected, "fNewDetectorProjected");
  amalgamated.cd(7); DrawTimingPeak(newAmalgamatedPairResidual, "fNewPairMean", true);
  amalgamated.cd(8); newAmalgamatedPairVsToT.SetStats(false);
  newAmalgamatedPairVsToT.Draw("colz");
  amalgamated.SaveAs(Form(
      "%s/CDet_Run%d_HydrogenCandidateComparison_Amalgamated.pdf", output,
      runNumber));
  amalgamated.SaveAs(Form(
      "%s/CDet_Run%d_HydrogenCandidateComparison_Amalgamated.png", output,
      runNumber));

  TFile rootOutput(Form("%s/CDet_Run%d_HydrogenCandidateComparison.root",
                        output, runNumber), "RECREATE");
  oldPairN.Write(); newPairN.Write(); oldPairTime.Write(); newPairTime.Write();
  oldPairResidual.Write(); newPairResidual.Write();
  oldTrajectory.Write(); newTrajectory.Write();
  oldSingleN.Write(); newSingleN.Write();
  oldSingleResidual.Write(); newSingleResidual.Write();
  oldPair2D.Write(); newPair2D.Write(); oldSingle2D.Write(); newSingle2D.Write();
  oldDetectorECalCut.Write(); newDetectorECalCut.Write();
  oldDetectorProjected.Write(); newDetectorProjected.Write();
  oldAmalgamatedPairResidual.Write(); newAmalgamatedPairResidual.Write();
  oldAmalgamatedPairVsToT.Write(); newAmalgamatedPairVsToT.Write();
  rootOutput.Close();

  std::ofstream summary(Form("%s/CDet_Run%d_HydrogenCandidateComparison.txt",
                             output, runNumber));
  summary << "run " << runNumber << "\nfiles " << filesAdded << "\nevents "
          << events << "\nold_pairs " << oldPairs << "\nnew_pairs " << newPairs
          << "\nold_singles " << oldSingles << "\nnew_singles " << newSingles
          << "\nmalformed_events " << malformed
          << "\npair_count_mismatch_events " << pairCountMismatch
          << "\npair_identity_mismatch_events " << pairIdentityMismatch
          << "\npair_value_mismatches " << pairValueMismatch
          << "\nsingle_count_mismatch_events " << singleCountMismatch
          << "\nsingle_identity_mismatch_events " << singleIdentityMismatch
          << "\nsingle_value_mismatches " << singleValueMismatch
          << "\nroi_status_mismatch_events " << roiMismatch
          << "\nold_amalgamated_pair_entries "
          << static_cast<Long64_t>(oldAmalgamatedPairResidual.GetEntries())
          << "\nnew_amalgamated_pair_entries "
          << static_cast<Long64_t>(newAmalgamatedPairResidual.GetEntries())
          << "\ntolerance " << tolerance << "\n";
  summary.close();

  std::cout << "\n[hydrogen candidate comparison] files/events: " << filesAdded
            << "/" << events << "\nold/new pairs: " << oldPairs << "/" << newPairs
            << "\nold/new singles: " << oldSingles << "/" << newSingles
            << "\npair count/identity/value mismatches: " << pairCountMismatch
            << "/" << pairIdentityMismatch << "/" << pairValueMismatch
            << "\nsingle count/identity/value mismatches: " << singleCountMismatch
            << "/" << singleIdentityMismatch << "/" << singleValueMismatch
            << "\nROI-status mismatches: " << roiMismatch
            << "\nmalformed pulse-array events: " << malformed
            << "\noutput: " << output << std::endl;
}
