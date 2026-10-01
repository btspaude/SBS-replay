#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>

#include "CDetRunDataset.h"
#include "CDetGoodPulseConfig.h"

namespace {

const int kCDetNPixels = 2688;
const int kCDetNBars = 168;
const int kPixelsPerBar = 16;

bool HasRequiredBranches(TChain &chain) {
  const char *required[] = {
      "earm.cdet.pulse.pmtnum", "earm.cdet.pulse.tdc_le_corr",
      "earm.cdet.pulse.tdc_te_corr", "earm.cdet.pulse.tdc_tot_ns",
      "earm.cdet.pulse.ecal_residual",
      "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass",
      "earm.cdet.pulse.ecal_eligible", "earm.cdet.pulse.spatial_pass",
      "earm.cdet.pulse.x_corr", "earm.cdet.pulse.y",
      "earm.cdet.pulse.z"};

  chain.LoadTree(0);
  for (const char *name : required) {
    if (!chain.GetBranch(name)) {
      std::cerr << "[good-pulse TDC] Missing required branch: " << name
                << std::endl;
      return false;
    }
  }
  const bool hasNativePairs =
      chain.GetBranch("earm.cdet.pair_candidate.pulse_index_l1") &&
      chain.GetBranch("earm.cdet.pair_candidate.pulse_index_l2");
  const bool hasLegacyPairs = chain.GetBranch("earm.cdet.pair.pulse_index_l1") &&
                              chain.GetBranch("earm.cdet.pair.pulse_index_l2");
  if (!hasNativePairs && !hasLegacyPairs) {
    std::cerr << "[good-pulse TDC] Missing both pair_candidate.* and pair.* "
                 "source branches."
              << std::endl;
    return false;
  }
  return true;
}

void DrawBarPage(TCanvas *canvas, const std::vector<TH1D *> &histograms,
                 int firstBar, int lastBar) {
  canvas->Clear();
  canvas->Divide(7, 6, 0.001, 0.001);
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    canvas->cd(bar - firstBar + 1);
    TH1D *hist = histograms[bar];
    hist->SetLineColor(kBlue + 1);
    hist->SetStats(false);
    hist->Draw();

    if (hist->GetEntries() >= 20 && hist->GetRMS() > 0.0) {
      const double fitLow = std::max(hist->GetXaxis()->GetXmin(),
                                     hist->GetMean() - 2.0 * hist->GetRMS());
      const double fitHigh = std::min(hist->GetXaxis()->GetXmax(),
                                      hist->GetMean() + 2.0 * hist->GetRMS());
      if (fitHigh > fitLow) {
        TF1 fit(Form("fGoodPulseBar%d", bar), "gaus", fitLow, fitHigh);
        fit.SetLineColor(kRed + 1);
        hist->Fit(&fit, "QNR");
        fit.DrawCopy("same");
      }
    }
  }
  canvas->Modified();
  canvas->Update();
}

} // namespace

// Reproduce the plotAllTDC timing views using only calibrated good pulse
// candidates.  A pulse is accepted when all four analyzer flags are true:
// calib_valid, broad_quality_pass, ecal_eligible, and spatial_pass.  An event
// contributes to the legacy-comparison timing views only if its accepted
// candidates include both CDet layers.  A separate diagnostic studies
// exclusive one-layer recovery in events without a selected pair. New trees
// use analyzer-produced pair_candidate.* and single_candidate.* collections;
// older trees fall back to the historical macro-level reconstruction.
//
// Example:
// root -l -b -q 'Plot_CDet_GoodPulseCandidates_AllTDC.C+(5710,"/path/to/Rootfiles","CDet_run5710_good_pulse_tdc")'
void Plot_CDet_GoodPulseCandidates_AllTDC(
    int runNumber = 5710, const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_good_pulse_tdc",
    double binWidthNs = 1.0, double leMinNs = 0.0, double leMaxNs = 60.0,
    double totMinNs = 0.0, double totMaxNs = 40.0,
    bool recoveredOnly = false, double pairResidualCenterM = 0.0,
    double pairTimingCenterNs = -26.0, double pairResidualScaleM = 0.020,
    double pairTimingScaleNs = 5.0, double pairCutRadius = 2.0,
    double singleResidualCenterM = 0.0,
    double singleTimingCenterNs = -26.0,
    double singleResidualScaleM = 0.040,
    double singleTimingScaleNs = 5.0, double singleCutRadius = 2.0,
    double ecalTimeMinNs = -10.0, double ecalTimeMaxNs = 10.0,
    double ecalEnergyMinGeV = 3.0, double ecalEnergyMaxGeV = 4.5,
    bool oppositeSideEnabled = false, double oppositeDYCenterM = 0.51,
    double oppositeDYToleranceM = 0.08,
    double oppositeProjectedYCenterM = 0.0,
    double oppositeProjectedYMaxM = 0.17) {
  TString resolvedInputDirectory(inputDirectory ? inputDirectory : "");
  if (resolvedInputDirectory.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      resolvedInputDirectory = outDir;
  }
  if (resolvedInputDirectory.IsNull()) {
    std::cerr << "[good-pulse TDC] An input directory or OUT_DIR is required."
              << std::endl;
    return;
  }
  if (binWidthNs <= 0.0 || leMaxNs <= leMinNs || totMaxNs <= totMinNs ||
      pairResidualScaleM <= 0.0 || pairTimingScaleNs <= 0.0 ||
      pairCutRadius <= 0.0 || singleResidualScaleM <= 0.0 ||
      singleTimingScaleNs <= 0.0 || singleCutRadius <= 0.0 ||
      ecalTimeMaxNs <= ecalTimeMinNs ||
      ecalEnergyMaxGeV <= ecalEnergyMinGeV ||
      (oppositeSideEnabled &&
       (oppositeDYCenterM <= 0.0 || oppositeDYToleranceM <= 0.0 ||
        oppositeProjectedYMaxM <= 0.0))) {
    std::cerr << "[good-pulse TDC] Invalid histogram limits." << std::endl;
    return;
  }

  TChain chain("T");
  const int filesAdded = CDetRunDataset::AddToChain(
      &chain, runNumber, resolvedInputDirectory.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[good-pulse TDC] No events found for run " << runNumber
              << " in " << resolvedInputDirectory << std::endl;
    return;
  }
  if (!HasRequiredBranches(chain))
    return;

  const int nLEBins = std::max(1, int(std::ceil((leMaxNs - leMinNs) /
                                                binWidthNs)));
  const double teMaxNs = leMaxNs + totMaxNs;
  const int nTEBins = std::max(1, int(std::ceil((teMaxNs - leMinNs) /
                                                binWidthNs)));
  const int nToTBins = std::max(1, int(std::ceil((totMaxNs - totMinNs) /
                                                 binWidthNs)));
  const int bar30FirstPixel = 30 * kPixelsPerBar;
  const double bar30DeltaTMinNs = -60.0;
  const double bar30DeltaTMaxNs = 30.0;
  const int nBar30DeltaTBins = std::max(
      1, int(std::ceil((bar30DeltaTMaxNs - bar30DeltaTMinNs) / binWidthNs)));

  const char *populationTitle = recoveredOnly
      ? "Recovered CDet pulses absent from legacy hAllGoodLe"
      : "Good CDet pulse candidates, two-layer events";
  TH1D hGoodLE("hCDetGoodPulseLE",
               Form("%s;Corrected LE time (ns);Pulses", populationTitle),
               nLEBins, leMinNs, leMaxNs);
  TH1D hGoodTE("hCDetGoodPulseTE",
               Form("%s;Corrected TE time (ns);Pulses", populationTitle),
               nTEBins, leMinNs, teMaxNs);
  TH1D hGoodToT("hCDetGoodPulseToT",
                Form("%s;Corrected ToT (ns);Pulses", populationTitle),
                nToTBins, totMinNs, totMaxNs);
  TH1D hGoodPixel("hCDetGoodPulsePixel",
                  "Good CDet pulse candidates;Pixel ID;Pulses",
                  kCDetNPixels, -0.5, kCDetNPixels - 0.5);
  TH1D hGoodBar("hCDetGoodPulseBar",
                "Good CDet pulse candidates;Bar ID;Pulses", kCDetNBars,
                -0.5, kCDetNBars - 0.5);
  TH1D hGoodMultiplicity(
      "hCDetGoodPulseMultiplicity",
      "Good CDet pulse candidates per event;Candidate pulses;Events", 101,
      -0.5, 100.5);
  TH1D hBar30ECalCut(
      "hCDetBar30ECalMinusCDet_ECalCut",
      "After ECal energy cut, CDet bar 30 (all instrumented pixels);"
      "t_{ECal} - t_{CDet,corr} (ns);Selected hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hBar30Projected(
      "hCDetBar30ECalMinusCDet_Projected",
      "ECal projection in CDet bar 30;"
      "t_{ECal} - t_{CDet,corr} (ns);Trajectory-matched hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hBar30ProjectedQuality(
      "hCDetBar30ECalMinusCDet_ProjectedQuality",
      "ECal projection + pulse-quality selection, CDet bar 30;"
      "t_{ECal} - t_{CDet,corr} (ns);Selected hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hBar30ProjectedQualityVsToT(
      "hCDetBar30ECalMinusCDetVsToT_ProjectedQuality",
      "ECal projection + pulse-quality selection, CDet bar 30;"
      "CDet ToT (ns);t_{ECal} - t_{CDet,corr} (ns)",
      nToTBins, totMinNs, totMaxNs, nBar30DeltaTBins,
      bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hDetectorECalCut(
      "hCDetDetectorECalMinusCDet_ECalCut",
      "ECal energy cut, events with an accepted CDet pair;"
      "t_{ECal} - t_{CDet,corr} (ns);Selected hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hDetectorProjected(
      "hCDetDetectorECalMinusCDet_Projected",
      "Projection-matched pulses, events with an accepted CDet pair;"
      "t_{ECal} - t_{CDet,corr} (ns);Trajectory-matched hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hDetectorProjectedQuality(
      "hCDetDetectorECalMinusCDet_ProjectedQuality",
      "Trajectory-time selected Layer-1/Layer-2 pair members, all CDet;"
      "t_{ECal} - t_{CDet,corr} (ns);Paired hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hDetectorProjectedQualityVsToT(
      "hCDetDetectorECalMinusCDetVsToT_ProjectedQuality",
      "Trajectory-time selected Layer-1/Layer-2 pair members, all CDet;"
      "CDet ToT (ns);t_{ECal} - t_{CDet,corr} (ns)",
      nToTBins, totMinNs, totMaxNs, nBar30DeltaTBins,
      bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hSelectedPairMeanResidual(
      "hCDetSelectedPairMeanECalResidual",
      "Trajectory-time selected pair mean, all CDet;"
      "t_{ECal} - <t_{CDet,corr}>_{pair} (ns);Selected pairs",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hSelectedPairMeanCorrectedLE(
      "hCDetSelectedPairMeanCorrectedLE",
      "Trajectory-time selected pair mean, all CDet;"
      "<t_{CDet,corr}>_{pair} (ns);Selected pairs",
      nLEBins, leMinNs, leMaxNs);
  TH2D hSelectedPairMeanResidualVsMeanToT(
      "hCDetSelectedPairMeanECalResidualVsMeanToT",
      "Trajectory-time selected pair mean, all CDet;"
      "<CDet ToT>_{pair} (ns);t_{ECal} - <t_{CDet,corr}>_{pair} (ns)",
      nToTBins, totMinNs, totMaxNs, nBar30DeltaTBins,
      bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hSelectedPairXCorrelation(
      "hCDetSelectedPairXVsProjectedECalX",
      "Trajectory-time selected pairs;"
      "x_{ECal} projected to pair mean z (m);"
      "<x_{CDet,corr}>_{pair} (m)",
      160, -1.6, 1.6, 160, -1.6, 1.6);
  TH1D hSelectedPairXResidual(
      "hCDetSelectedPairXResidual",
      "Trajectory-time selected pairs;"
      "<x_{CDet,corr}>_{pair} - x_{ECal projected} (m);Selected pairs",
      160, -0.20, 0.20);
  TH1D hPairTrajectoryResidual(
      "hCDetPairTrajectoryResidual",
      "All accepted pairs;#Delta x_{pair} - (x_{ECal}/z_{ECal})#Delta z (m);Pairs",
      160, -0.08, 0.08);
  TH1D hBestPairTrajectoryResidual(
      "hCDetBestPairTrajectoryResidual",
      "Best trajectory-matched pair per event;#Delta x_{pair} - (x_{ECal}/z_{ECal})#Delta z (m);Events",
      160, -0.08, 0.08);
  TH2D hPairTimingVsTrajectoryResidual(
      "hCDetPairTimingVsTrajectoryResidual",
      "Accepted pairs;#Delta x_{pair} - (x_{ECal}/z_{ECal})#Delta z (m);"
      "t_{ECal} - <t_{CDet}>_{pair} (ns)",
      160, -0.08, 0.08, nBar30DeltaTBins, bar30DeltaTMinNs,
      bar30DeltaTMaxNs);
  TH2D hSingleLayer1TimingVsXResidual(
      "hCDetSingleLayer1TimingVsXResidual",
      "Pairless one-layer events, Layer 1 good pulses;"
      "x_{CDet,corr} - x_{ECal projected} (m);"
      "t_{ECal} - t_{CDet,corr} (ns)",
      160, -0.16, 0.16, nBar30DeltaTBins, bar30DeltaTMinNs,
      bar30DeltaTMaxNs);
  TH2D hSingleLayer2TimingVsXResidual(
      "hCDetSingleLayer2TimingVsXResidual",
      "Pairless one-layer events, Layer 2 good pulses;"
      "x_{CDet,corr} - x_{ECal projected} (m);"
      "t_{ECal} - t_{CDet,corr} (ns)",
      160, -0.16, 0.16, nBar30DeltaTBins, bar30DeltaTMinNs,
      bar30DeltaTMaxNs);
  TH1D hSelectedSingleLayer1Timing(
      "hCDetSelectedSingleLayer1Timing",
      "Selected single-layer recovery, Layer 1;"
      "t_{ECal} - t_{CDet,corr} (ns);Selected pulses",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hSelectedSingleLayer2Timing(
      "hCDetSelectedSingleLayer2Timing",
      "Selected single-layer recovery, Layer 2;"
      "t_{ECal} - t_{CDet,corr} (ns);Selected pulses",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hECalAdmittedEventOutcome(
      "hCDetECalAdmittedEventOutcome",
      "CDet outcome for ECal-energy-and-time admitted events;"
      "Exclusive event category;Events",
      4, 0.5, 4.5);
  hECalAdmittedEventOutcome.GetXaxis()->SetBinLabel(1, "Selected pair");
  hECalAdmittedEventOutcome.GetXaxis()->SetBinLabel(
      2, "Good both; pair fails");
  hECalAdmittedEventOutcome.GetXaxis()->SetBinLabel(
      3, "Calib both; quality/spatial loss");
  hECalAdmittedEventOutcome.GetXaxis()->SetBinLabel(
      4, "Missing calibrated layer");
  TH2D hFailedPairBestTimingVsTrajectory(
      "hCDetFailedPairBestTimingVsTrajectory",
      "Best stored pair in good-both events failing final ellipse;"
      "trajectory residual (m);t_{ECal} - <t_{CDet,corr}>_{pair} (ns)",
      160, -0.16, 0.16, nBar30DeltaTBins, bar30DeltaTMinNs,
      bar30DeltaTMaxNs);
  TH1D hFailedPairMinimumRadius(
      "hCDetFailedPairMinimumNormalizedRadius",
      "Good-both events failing final ellipse;minimum normalized pair radius;"
      "Events",
      120, 0.0, 12.0);
  TH1D hNoStoredPairFailureReason(
      "hCDetNoStoredPairFailureReason",
      "Nearest cross-layer combination when no stored pair;failed analyzer gate;"
      "Events",
      7, 0.5, 7.5);
  const char *noStoredFailureLabels[7] = {
      "#Deltat", "#Deltax", "#Deltay", "#Deltat+#Deltax",
      "#Deltat+#Deltay", "#Deltax+#Deltay", "all three"};
  for (int bin = 1; bin <= 7; ++bin)
    hNoStoredPairFailureReason.GetXaxis()->SetBinLabel(
        bin, noStoredFailureLabels[bin - 1]);
  TH1D hAllCombinationRecoveryOutcome(
      "hCDetAllCombinationRecoveryOutcome",
      "Good-both events rejected by stored-pair ellipse;all-combination outcome;"
      "Events",
      3, 0.5, 3.5);
  hAllCombinationRecoveryOutcome.GetXaxis()->SetBinLabel(
      1, "No combination in ellipse");
  hAllCombinationRecoveryOutcome.GetXaxis()->SetBinLabel(
      2, "Ellipse pair fails hard gates");
  hAllCombinationRecoveryOutcome.GetXaxis()->SetBinLabel(
      3, "Greedy pairing lost ellipse pair");
  TH1D hEllipsePairHardGateFailureReason(
      "hCDetEllipsePairHardGateFailureReason",
      "Best ellipse combination rejected before pairing;failed analyzer gate;"
      "Events",
      7, 0.5, 7.5);
  for (int bin = 1; bin <= 7; ++bin)
    hEllipsePairHardGateFailureReason.GetXaxis()->SetBinLabel(
        bin, noStoredFailureLabels[bin - 1]);
  TH1D hAlternativePairMultiplicity(
      "hCDetAlternativePairMultiplicity",
      "Recovered failed events;hard-gate-eligible combinations inside ellipse;"
      "Compatible combinations;Events",
      51, -0.5, 50.5);
  TH2D hRecoveredGreedyBestTimingVsTrajectory(
      "hCDetRecoveredGreedyBestTimingVsTrajectory",
      "Best ECal-informed pair recovered after greedy loss;"
      "trajectory residual (m);t_{ECal} - <t_{CDet,corr}>_{pair} (ns)",
      160, -0.08, 0.08, nBar30DeltaTBins, bar30DeltaTMinNs,
      bar30DeltaTMaxNs);
  TH1D hRecoveredGreedyBestMeanResidual(
      "hCDetRecoveredGreedyBestMeanECalResidual",
      "Best ECal-informed pair recovered after greedy loss;"
      "t_{ECal} - <t_{CDet,corr}>_{pair} (ns);Recovered events",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hRecoveredGreedyBestXCorrelation(
      "hCDetRecoveredGreedyBestXVsProjectedECalX",
      "Best ECal-informed pair recovered after greedy loss;"
      "x_{ECal} projected to pair mean z (m);"
      "<x_{CDet,corr}>_{pair} (m)",
      160, -1.6, 1.6, 160, -1.6, 1.6);
  TH1D hRecoveredGreedyBestXResidual(
      "hCDetRecoveredGreedyBestXResidual",
      "Best ECal-informed pair recovered after greedy loss;"
      "<x_{CDet,corr}>_{pair} - x_{ECal projected} (m);Recovered events",
      160, -0.20, 0.20);

  std::vector<TH1D *> pixelLE(kCDetNPixels, nullptr);
  for (int pixel = 0; pixel < kCDetNPixels; ++pixel) {
    pixelLE[pixel] =
        new TH1D(Form("hCDetGoodPulseLE_pixel%04d", pixel),
                 Form("Good pulse LE, pixel %d;Corrected LE time (ns);Pulses",
                      pixel),
                 nLEBins, leMinNs, leMaxNs);
  }

  std::vector<TH1D *> barLE(kCDetNBars, nullptr);
  for (int bar = 0; bar < kCDetNBars; ++bar) {
    barLE[bar] =
        new TH1D(Form("hCDetGoodPulseLE_bar%03d", bar),
                 Form("Good pulse LE, bar %d;Corrected LE time (ns);Pulses",
                      bar),
                 nLEBins, leMinNs, leMaxNs);
  }

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> te(reader, "earm.cdet.pulse.tdc_te_corr");
  TTreeReaderArray<Double_t> tot(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<Double_t> ecalResidual(
      reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> calibValid(reader,
                                        "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> ecalEligible(reader,
                                          "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatialPass(reader,
                                         "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<Double_t> correctedX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> pulseY(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  const bool useAnalyzerPairCandidates =
      chain.GetBranch("earm.cdet.pair_candidate.pulse_index_l1") &&
      chain.GetBranch("earm.cdet.pair_candidate.pulse_index_l2");
  const bool useAnalyzerSingleCandidates =
      chain.GetBranch("earm.cdet.single_candidate.pulse_index") &&
      chain.GetBranch("earm.cdet.single_candidate.layer") &&
      chain.GetBranch("earm.cdet.single_candidate.ecal_residual") &&
      chain.GetBranch("earm.cdet.single_candidate.x_residual");
  const bool hasROIStatus = chain.GetBranch("earm.cdet.roi.status");
  const char *pairIndexL1Branch = useAnalyzerPairCandidates
      ? "earm.cdet.pair_candidate.pulse_index_l1"
      : "earm.cdet.pair.pulse_index_l1";
  const char *pairIndexL2Branch = useAnalyzerPairCandidates
      ? "earm.cdet.pair_candidate.pulse_index_l2"
      : "earm.cdet.pair.pulse_index_l2";
  TTreeReaderArray<Double_t> pairPulseIndexL1(reader, pairIndexL1Branch);
  TTreeReaderArray<Double_t> pairPulseIndexL2(reader, pairIndexL2Branch);
  std::unique_ptr<TTreeReaderArray<Double_t>> singlePulseIndex;
  std::unique_ptr<TTreeReaderArray<Double_t>> singleLayer;
  std::unique_ptr<TTreeReaderArray<Double_t>> singleECalResidual;
  std::unique_ptr<TTreeReaderArray<Double_t>> singleXResidual;
  std::unique_ptr<TTreeReaderArray<Double_t>> candidateGreedySelected;
  std::unique_ptr<TTreeReaderArray<Double_t>> candidateSelectedPairIndex;
  std::unique_ptr<TTreeReaderArray<Double_t>> selectedPairPulseIndexL1;
  std::unique_ptr<TTreeReaderArray<Double_t>> selectedPairPulseIndexL2;
  std::unique_ptr<TTreeReaderValue<Double_t>> roiStatus;
  if (useAnalyzerPairCandidates &&
      chain.GetBranch("earm.cdet.pair_candidate.greedy_selected") &&
      chain.GetBranch("earm.cdet.pair_candidate.selected_pair_index") &&
      chain.GetBranch("earm.cdet.pair.pulse_index_l1") &&
      chain.GetBranch("earm.cdet.pair.pulse_index_l2")) {
    candidateGreedySelected.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.pair_candidate.greedy_selected"));
    candidateSelectedPairIndex.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.pair_candidate.selected_pair_index"));
    selectedPairPulseIndexL1.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.pair.pulse_index_l1"));
    selectedPairPulseIndexL2.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.pair.pulse_index_l2"));
  }
  if (useAnalyzerSingleCandidates) {
    singlePulseIndex.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.single_candidate.pulse_index"));
    singleLayer.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.single_candidate.layer"));
    singleECalResidual.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.single_candidate.ecal_residual"));
    singleXResidual.reset(new TTreeReaderArray<Double_t>(
        reader, "earm.cdet.single_candidate.x_residual"));
  }
  if (hasROIStatus)
    roiStatus.reset(
        new TTreeReaderValue<Double_t>(reader, "earm.cdet.roi.status"));
  TTreeReaderArray<Double_t> legacyPMT(reader, "earm.cdet.hit.pmtnum");
  TTreeReaderArray<Double_t> legacyLE(reader, "earm.cdet.hit.tdc_le");
  TTreeReaderArray<Double_t> legacyToT(reader, "earm.cdet.hit.tdc_tot");
  TTreeReaderArray<Double_t> legacyX(reader, "earm.cdet.hit.xhit");
  TTreeReaderArray<Double_t> legacyY(reader, "earm.cdet.hit.yhit");
  TTreeReaderArray<Double_t> legacyZ(reader, "earm.cdet.hit.zhit");
  TTreeReaderArray<Double_t> legacyMultiplicity(reader, "earm.cdet.tdc_mult");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");

  Long64_t eventCount = 0;
  Long64_t pulseCount = 0;
  Long64_t goodPulseCount = 0;
  Long64_t twoLayerEventCount = 0;
  Long64_t plottedEventCount = 0;
  Long64_t malformedEvents = 0;
  Long64_t ecalEnergyEventCount = 0;
  Long64_t detectorBaselineEventCount = 0;
  Long64_t detectorProjectedEventCount = 0;
  Long64_t detectorProjectedQualityEventCount = 0;
  Long64_t acceptedPairCount = 0;
  Long64_t trajectoryTimeSelectedPairCount = 0;
  Long64_t pairlessSingleLayer1EventCount = 0;
  Long64_t pairlessSingleLayer2EventCount = 0;
  Long64_t recoveredSingleLayer1EventCount = 0;
  Long64_t recoveredSingleLayer2EventCount = 0;
  Long64_t recoveredSingleLayer1PulseCount = 0;
  Long64_t recoveredSingleLayer2PulseCount = 0;
  Long64_t ecalAdmittedEventCount = 0;
  Long64_t selectedPairEventCount = 0;
  Long64_t goodBothPairFailureEventCount = 0;
  Long64_t goodBothNoStoredPairEventCount = 0;
  Long64_t goodBothEllipseFailureEventCount = 0;
  Long64_t calibratedBothQualitySpatialLossEventCount = 0;
  Long64_t goodLayer1OnlyEventCount = 0;
  Long64_t goodLayer2OnlyEventCount = 0;
  Long64_t noGoodLayerEventCount = 0;
  Long64_t missingCalibratedLayerEventCount = 0;
  Long64_t calibratedLayer1OnlyEventCount = 0;
  Long64_t calibratedLayer2OnlyEventCount = 0;
  Long64_t noCalibratedPulseEventCount = 0;
  Long64_t failedPairTimingOnlyEventCount = 0;
  Long64_t failedPairTrajectoryOnlyEventCount = 0;
  Long64_t failedPairBothAxesEventCount = 0;
  Long64_t failedPairEllipseCornerEventCount = 0;
  Long64_t noStoredPairUnexpectedEventCount = 0;
  Long64_t noAllGoodEllipseCombinationEventCount = 0;
  Long64_t ellipseCombinationFailsHardGatesEventCount = 0;
  Long64_t greedyPairingLostEllipseCombinationEventCount = 0;
  Long64_t allCombinationRecoveredPairCount = 0;
  Long64_t analyzerSingleCandidateEventCount = 0;
  Long64_t pairClassificationCountMismatch = 0;
  Long64_t pairClassificationIdentityMismatch = 0;
  Long64_t singleCandidateCountMismatch = 0;
  Long64_t singleCandidateIdentityMismatch = 0;
  Long64_t roiStatusMismatch = 0;
  constexpr double kECalZFromTargetM = 6.144;
  constexpr double kAnalyzerPairDeltaTimeMaxNs = 15.0;
  constexpr double kAnalyzerPairDeltaXMaxM = 0.15;
  constexpr double kAnalyzerPairDeltaYMaxM = 0.08;

  while (reader.Next()) {
    ++eventCount;
    const std::vector<size_t> sizes = {
        pixel.GetSize(),        le.GetSize(),          te.GetSize(),
        tot.GetSize(),          ecalResidual.GetSize(), calibValid.GetSize(),
        broadQuality.GetSize(), ecalEligible.GetSize(), spatialPass.GetSize(),
        correctedX.GetSize(), pulseY.GetSize(), pulseZ.GetSize()};
    const size_t nPulses = *std::min_element(sizes.begin(), sizes.end());
    const size_t maxPulses = *std::max_element(sizes.begin(), sizes.end());
    if (nPulses != maxPulses)
      ++malformedEvents;
    pulseCount += nPulses;
    const bool passesECalTime =
        std::isfinite(*ecalTime) && *ecalTime >= ecalTimeMinNs &&
        *ecalTime <= ecalTimeMaxNs;

    // Analyzer-native reproduction of the historical Bar 30 timing canvas.
    // The saved residual has the desired historical sign, ECal minus CDet.
    if (std::isfinite(*ecalEnergy) && *ecalEnergy >= ecalEnergyMinGeV &&
        *ecalEnergy <= ecalEnergyMaxGeV) {
      ++ecalEnergyEventCount;
      std::unordered_set<size_t> pairedPulseIndices;
      std::unordered_set<size_t> trajectoryTimeSelectedPulseIndices;
      const size_t nPairs =
          std::min(pairPulseIndexL1.GetSize(), pairPulseIndexL2.GetSize());
      if (candidateGreedySelected) {
        const size_t nClassified = std::min(
            {nPairs, candidateGreedySelected->GetSize(),
             candidateSelectedPairIndex->GetSize()});
        size_t selectedCandidateCount = 0;
        for (size_t candidate = 0; candidate < nClassified; ++candidate) {
          if ((*candidateGreedySelected)[candidate] <= 0.5)
            continue;
          ++selectedCandidateCount;
          const Long64_t selectedIndex =
              std::llround((*candidateSelectedPairIndex)[candidate]);
          if (selectedIndex < 0 ||
              selectedIndex >=
                  static_cast<Long64_t>(selectedPairPulseIndexL1->GetSize()) ||
              selectedIndex >=
                  static_cast<Long64_t>(selectedPairPulseIndexL2->GetSize()) ||
              std::llround(pairPulseIndexL1[candidate]) !=
                  std::llround((*selectedPairPulseIndexL1)[selectedIndex]) ||
              std::llround(pairPulseIndexL2[candidate]) !=
                  std::llround((*selectedPairPulseIndexL2)[selectedIndex]))
            ++pairClassificationIdentityMismatch;
        }
        if (nClassified != nPairs ||
            selectedCandidateCount != selectedPairPulseIndexL1->GetSize() ||
            selectedCandidateCount != selectedPairPulseIndexL2->GetSize())
          ++pairClassificationCountMismatch;
      }
      double bestTrajectoryAbsResidual = std::numeric_limits<double>::infinity();
      double bestTrajectoryResidual = std::numeric_limits<double>::quiet_NaN();
      for (size_t pair = 0; pair < nPairs; ++pair) {
        if (!std::isfinite(pairPulseIndexL1[pair]) ||
            !std::isfinite(pairPulseIndexL2[pair]))
          continue;
        const Long64_t indexL1 = std::llround(pairPulseIndexL1[pair]);
        const Long64_t indexL2 = std::llround(pairPulseIndexL2[pair]);
        if (indexL1 < 0 || indexL2 < 0 ||
            indexL1 >= static_cast<Long64_t>(nPulses) ||
            indexL2 >= static_cast<Long64_t>(nPulses))
          continue;
        pairedPulseIndices.insert(static_cast<size_t>(indexL1));
        pairedPulseIndices.insert(static_cast<size_t>(indexL2));
        const double deltaZ = pulseZ[indexL2] - pulseZ[indexL1];
        const double trajectoryResidual =
            (correctedX[indexL2] - correctedX[indexL1]) -
            (*ecalX / kECalZFromTargetM) * deltaZ;
        const double pairTimingResidual =
            0.5 * (ecalResidual[indexL1] + ecalResidual[indexL2]);
        if (std::isfinite(trajectoryResidual) &&
            std::isfinite(pairTimingResidual)) {
          hPairTrajectoryResidual.Fill(trajectoryResidual);
          hPairTimingVsTrajectoryResidual.Fill(trajectoryResidual,
                                               pairTimingResidual);
          if (std::fabs(trajectoryResidual) < bestTrajectoryAbsResidual) {
            bestTrajectoryAbsResidual = std::fabs(trajectoryResidual);
            bestTrajectoryResidual = trajectoryResidual;
          }
          const double normalizedResidual =
              (trajectoryResidual - pairResidualCenterM) /
              pairResidualScaleM;
          const double normalizedTiming =
              (pairTimingResidual - pairTimingCenterNs) /
              pairTimingScaleNs;
          if (passesECalTime &&
              normalizedResidual * normalizedResidual +
                  normalizedTiming * normalizedTiming <=
              pairCutRadius * pairCutRadius) {
            trajectoryTimeSelectedPulseIndices.insert(
                static_cast<size_t>(indexL1));
            trajectoryTimeSelectedPulseIndices.insert(
                static_cast<size_t>(indexL2));
            hSelectedPairMeanResidual.Fill(pairTimingResidual);
            const double pairMeanCorrectedLE =
                0.5 * (le[indexL1] + le[indexL2]);
            if (std::isfinite(pairMeanCorrectedLE))
              hSelectedPairMeanCorrectedLE.Fill(pairMeanCorrectedLE);
            const double pairMeanToT =
                0.5 * (tot[indexL1] + tot[indexL2]);
            if (std::isfinite(pairMeanToT))
              hSelectedPairMeanResidualVsMeanToT.Fill(pairMeanToT,
                                                      pairTimingResidual);
            const double pairMeanZ =
                0.5 * (pulseZ[indexL1] + pulseZ[indexL2]);
            const double pairMeanX =
                0.5 * (correctedX[indexL1] + correctedX[indexL2]);
            const double projectedECalX =
                *ecalX * pairMeanZ / kECalZFromTargetM;
            if (std::isfinite(pairMeanX) && std::isfinite(projectedECalX)) {
              hSelectedPairXCorrelation.Fill(projectedECalX, pairMeanX);
              hSelectedPairXResidual.Fill(pairMeanX - projectedECalX);
            }
            ++trajectoryTimeSelectedPairCount;
          }
        }
        ++acceptedPairCount;
      }
      if (std::isfinite(bestTrajectoryResidual))
        hBestPairTrajectoryResidual.Fill(bestTrajectoryResidual);
      const bool hasTrajectoryTimeSelectedPair =
          !trajectoryTimeSelectedPulseIndices.empty();

      // Give every event admitted by the configured ECal energy and timing
      // cuts one mutually exclusive CDet outcome.  "Calibrated" means a
      // complete analyzer pulse with valid calibration; "good" additionally
      // requires broad quality, ECal eligibility, and spatial compatibility.
      if (passesECalTime) {
        ++ecalAdmittedEventCount;
        bool calibratedLayer1 = false;
        bool calibratedLayer2 = false;
        bool goodLayer1 = false;
        bool goodLayer2 = false;
        for (size_t i = 0; i < nPulses; ++i) {
          if (!std::isfinite(pixel[i]))
            continue;
          const int pixelID = int(std::lround(pixel[i]));
          if (pixelID < 0 || pixelID >= kCDetNPixels)
            continue;
          const bool isLayer1 = pixelID < kCDetNPixels / 2;
          const bool calibrated = calibValid[i] > 0.5;
          const bool good = calibrated && broadQuality[i] > 0.5 &&
              ecalEligible[i] > 0.5 && spatialPass[i] > 0.5;
          if (isLayer1) {
            calibratedLayer1 |= calibrated;
            goodLayer1 |= good;
          } else {
            calibratedLayer2 |= calibrated;
            goodLayer2 |= good;
          }
        }

        if (hasTrajectoryTimeSelectedPair) {
          ++selectedPairEventCount;
          hECalAdmittedEventOutcome.Fill(1.0);
        } else if (goodLayer1 && goodLayer2) {
          ++goodBothPairFailureEventCount;
          std::vector<size_t> goodLayer1Indices;
          std::vector<size_t> goodLayer2Indices;
          for (size_t i = 0; i < nPulses; ++i) {
            if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
                  ecalEligible[i] > 0.5 && spatialPass[i] > 0.5) ||
                !std::isfinite(pixel[i]) || !std::isfinite(le[i]) ||
                !std::isfinite(ecalResidual[i]) ||
                !std::isfinite(correctedX[i]) ||
                !std::isfinite(pulseY[i]) || !std::isfinite(pulseZ[i]))
              continue;
            const int pixelID = int(std::lround(pixel[i]));
            if (pixelID < 0 || pixelID >= kCDetNPixels)
              continue;
            if (pixelID < kCDetNPixels / 2)
              goodLayer1Indices.push_back(i);
            else
              goodLayer2Indices.push_back(i);
          }

          bool hasAllGoodEllipseCombination = false;
          int hardGateEllipseCombinationCount = 0;
          double bestEllipseRadiusSquared =
              std::numeric_limits<double>::infinity();
          int bestEllipseFailureMask = 0;
          double bestHardGateRadiusSquared =
              std::numeric_limits<double>::infinity();
          size_t bestHardGateIndexL1 = 0;
          size_t bestHardGateIndexL2 = 0;
          double bestHardGateTrajectoryResidual =
              std::numeric_limits<double>::quiet_NaN();
          double bestHardGateTimingResidual =
              std::numeric_limits<double>::quiet_NaN();
          for (size_t indexL1 : goodLayer1Indices) {
            for (size_t indexL2 : goodLayer2Indices) {
              const double dt = le[indexL2] - le[indexL1];
              const double dx = correctedX[indexL2] - correctedX[indexL1];
              const double dy = pulseY[indexL2] - pulseY[indexL1];
              const double alignedProjectedY = *ecalY *
                  0.5 * (pulseZ[indexL1] + pulseZ[indexL2]) /
                  kECalZFromTargetM + 0.10;
              const bool sameSide =
                  std::fabs(dy) <= kAnalyzerPairDeltaYMaxM;
              const bool oppositeSide = oppositeSideEnabled &&
                  std::fabs(std::fabs(dy) - oppositeDYCenterM) <=
                      oppositeDYToleranceM &&
                  std::fabs(alignedProjectedY -
                            oppositeProjectedYCenterM) <=
                      oppositeProjectedYMaxM;
              const bool passesYTopology = sameSide || oppositeSide;
              const double deltaZ = pulseZ[indexL2] - pulseZ[indexL1];
              const double trajectoryResidual =
                  dx - (*ecalX / kECalZFromTargetM) * deltaZ;
              const double pairTimingResidual =
                  0.5 * (ecalResidual[indexL1] + ecalResidual[indexL2]);
              const double normalizedResidual =
                  (trajectoryResidual - pairResidualCenterM) /
                  pairResidualScaleM;
              const double normalizedTiming =
                  (pairTimingResidual - pairTimingCenterNs) /
                  pairTimingScaleNs;
              const bool insideEllipse =
                  normalizedResidual * normalizedResidual +
                      normalizedTiming * normalizedTiming <=
                  pairCutRadius * pairCutRadius;
              if (!insideEllipse)
                continue;
              hasAllGoodEllipseCombination = true;
              const double radiusSquared =
                  normalizedResidual * normalizedResidual +
                  normalizedTiming * normalizedTiming;
              const int failureMask =
                  (std::fabs(dt) > kAnalyzerPairDeltaTimeMaxNs ? 1 : 0) |
                  (std::fabs(dx) > kAnalyzerPairDeltaXMaxM ? 2 : 0) |
                  (!passesYTopology ? 4 : 0);
              if (radiusSquared < bestEllipseRadiusSquared) {
                bestEllipseRadiusSquared = radiusSquared;
                bestEllipseFailureMask = failureMask;
              }
              const bool passesAnalyzerHardGates =
                  std::fabs(dt) <= kAnalyzerPairDeltaTimeMaxNs &&
                  std::fabs(dx) <= kAnalyzerPairDeltaXMaxM &&
                  passesYTopology;
              if (passesAnalyzerHardGates) {
                ++hardGateEllipseCombinationCount;
                if (radiusSquared < bestHardGateRadiusSquared) {
                  bestHardGateRadiusSquared = radiusSquared;
                  bestHardGateIndexL1 = indexL1;
                  bestHardGateIndexL2 = indexL2;
                  bestHardGateTrajectoryResidual = trajectoryResidual;
                  bestHardGateTimingResidual = pairTimingResidual;
                }
              }
            }
          }

          if (!hasAllGoodEllipseCombination) {
            ++noAllGoodEllipseCombinationEventCount;
            hAllCombinationRecoveryOutcome.Fill(1.0);
          } else if (hardGateEllipseCombinationCount == 0) {
            ++ellipseCombinationFailsHardGatesEventCount;
            hAllCombinationRecoveryOutcome.Fill(2.0);
            int failureBin = 0;
            switch (bestEllipseFailureMask) {
            case 1: failureBin = 1; break;
            case 2: failureBin = 2; break;
            case 4: failureBin = 3; break;
            case 3: failureBin = 4; break;
            case 5: failureBin = 5; break;
            case 6: failureBin = 6; break;
            case 7: failureBin = 7; break;
            default: break;
            }
            if (failureBin > 0)
              hEllipsePairHardGateFailureReason.Fill(failureBin);
          } else {
            ++greedyPairingLostEllipseCombinationEventCount;
            allCombinationRecoveredPairCount += hardGateEllipseCombinationCount;
            hAlternativePairMultiplicity.Fill(hardGateEllipseCombinationCount);
            hAllCombinationRecoveryOutcome.Fill(3.0);
            hRecoveredGreedyBestTimingVsTrajectory.Fill(
                bestHardGateTrajectoryResidual, bestHardGateTimingResidual);
            hRecoveredGreedyBestMeanResidual.Fill(
                bestHardGateTimingResidual);
            const double pairMeanZ =
                0.5 * (pulseZ[bestHardGateIndexL1] +
                       pulseZ[bestHardGateIndexL2]);
            const double pairMeanX =
                0.5 * (correctedX[bestHardGateIndexL1] +
                       correctedX[bestHardGateIndexL2]);
            const double projectedECalX =
                *ecalX * pairMeanZ / kECalZFromTargetM;
            if (std::isfinite(pairMeanX) && std::isfinite(projectedECalX)) {
              hRecoveredGreedyBestXCorrelation.Fill(projectedECalX,
                                                     pairMeanX);
              hRecoveredGreedyBestXResidual.Fill(pairMeanX - projectedECalX);
            }
          }

          if (pairedPulseIndices.empty()) {
            ++goodBothNoStoredPairEventCount;
            double bestGateMetric = std::numeric_limits<double>::infinity();
            int bestFailureMask = 0;
            for (size_t indexL1 : goodLayer1Indices) {
              for (size_t indexL2 : goodLayer2Indices) {
                const double dt = le[indexL2] - le[indexL1];
                const double dx = correctedX[indexL2] - correctedX[indexL1];
                const double dy = pulseY[indexL2] - pulseY[indexL1];
                const double dtRatio =
                    std::fabs(dt) / kAnalyzerPairDeltaTimeMaxNs;
                const double dxRatio =
                    std::fabs(dx) / kAnalyzerPairDeltaXMaxM;
                const double dyRatio =
                    std::fabs(dy) / kAnalyzerPairDeltaYMaxM;
                const double gateMetric = std::max({dtRatio, dxRatio, dyRatio});
                if (gateMetric < bestGateMetric) {
                  bestGateMetric = gateMetric;
                  bestFailureMask = (dtRatio > 1.0 ? 1 : 0) |
                      (dxRatio > 1.0 ? 2 : 0) |
                      (dyRatio > 1.0 ? 4 : 0);
                }
              }
            }
            int failureBin = 0;
            switch (bestFailureMask) {
            case 1: failureBin = 1; break;
            case 2: failureBin = 2; break;
            case 4: failureBin = 3; break;
            case 3: failureBin = 4; break;
            case 5: failureBin = 5; break;
            case 6: failureBin = 6; break;
            case 7: failureBin = 7; break;
            default: break;
            }
            if (failureBin > 0)
              hNoStoredPairFailureReason.Fill(failureBin);
            else
              ++noStoredPairUnexpectedEventCount;
          } else {
            ++goodBothEllipseFailureEventCount;
            double bestRadiusSquared = std::numeric_limits<double>::infinity();
            double bestResidual = std::numeric_limits<double>::quiet_NaN();
            double bestTiming = std::numeric_limits<double>::quiet_NaN();
            for (size_t pair = 0; pair < nPairs; ++pair) {
              if (!std::isfinite(pairPulseIndexL1[pair]) ||
                  !std::isfinite(pairPulseIndexL2[pair]))
                continue;
              const Long64_t indexL1 = std::llround(pairPulseIndexL1[pair]);
              const Long64_t indexL2 = std::llround(pairPulseIndexL2[pair]);
              if (indexL1 < 0 || indexL2 < 0 ||
                  indexL1 >= static_cast<Long64_t>(nPulses) ||
                  indexL2 >= static_cast<Long64_t>(nPulses))
                continue;
              const double deltaZ = pulseZ[indexL2] - pulseZ[indexL1];
              const double residual =
                  (correctedX[indexL2] - correctedX[indexL1]) -
                  (*ecalX / kECalZFromTargetM) * deltaZ;
              const double timing =
                  0.5 * (ecalResidual[indexL1] + ecalResidual[indexL2]);
              const double nx =
                  (residual - pairResidualCenterM) / pairResidualScaleM;
              const double nt =
                  (timing - pairTimingCenterNs) / pairTimingScaleNs;
              const double radiusSquared = nx * nx + nt * nt;
              if (std::isfinite(radiusSquared) &&
                  radiusSquared < bestRadiusSquared) {
                bestRadiusSquared = radiusSquared;
                bestResidual = residual;
                bestTiming = timing;
              }
            }
            if (std::isfinite(bestRadiusSquared)) {
              const double minimumRadius = std::sqrt(bestRadiusSquared);
              hFailedPairMinimumRadius.Fill(minimumRadius);
              hFailedPairBestTimingVsTrajectory.Fill(bestResidual, bestTiming);
              const bool trajectoryAxisFails =
                  std::fabs((bestResidual - pairResidualCenterM) /
                            pairResidualScaleM) > pairCutRadius;
              const bool timingAxisFails =
                  std::fabs((bestTiming - pairTimingCenterNs) /
                            pairTimingScaleNs) > pairCutRadius;
              if (timingAxisFails && !trajectoryAxisFails)
                ++failedPairTimingOnlyEventCount;
              else if (trajectoryAxisFails && !timingAxisFails)
                ++failedPairTrajectoryOnlyEventCount;
              else if (trajectoryAxisFails && timingAxisFails)
                ++failedPairBothAxesEventCount;
              else
                ++failedPairEllipseCornerEventCount;
            }
          }
          hECalAdmittedEventOutcome.Fill(2.0);
        } else if (calibratedLayer1 && calibratedLayer2) {
          ++calibratedBothQualitySpatialLossEventCount;
          if (goodLayer1)
            ++goodLayer1OnlyEventCount;
          else if (goodLayer2)
            ++goodLayer2OnlyEventCount;
          else
            ++noGoodLayerEventCount;
          hECalAdmittedEventOutcome.Fill(3.0);
        } else {
          ++missingCalibratedLayerEventCount;
          if (calibratedLayer1)
            ++calibratedLayer1OnlyEventCount;
          else if (calibratedLayer2)
            ++calibratedLayer2OnlyEventCount;
          else
            ++noCalibratedPulseEventCount;
          hECalAdmittedEventOutcome.Fill(4.0);
        }
      }

      // Reconstruct the historical single-layer definition independently for
      // diagnostics and regression testing. On new trees, the selected sample
      // itself comes from SBSCDet's single_candidate.* collection.
      if (!hasTrajectoryTimeSelectedPair) {
        std::vector<size_t> goodLayer1Indices;
        std::vector<size_t> goodLayer2Indices;
        for (size_t i = 0; i < nPulses; ++i) {
          if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
                ecalEligible[i] > 0.5 && spatialPass[i] > 0.5) ||
              !std::isfinite(pixel[i]) || !std::isfinite(ecalResidual[i]) ||
              !std::isfinite(correctedX[i]) || !std::isfinite(pulseZ[i]))
            continue;
          const int pixelID = int(std::lround(pixel[i]));
          if (pixelID < 0 || pixelID >= kCDetNPixels)
            continue;
          if (pixelID < kCDetNPixels / 2)
            goodLayer1Indices.push_back(i);
          else
            goodLayer2Indices.push_back(i);
        }

        const bool layer1Only =
            !goodLayer1Indices.empty() && goodLayer2Indices.empty();
        const bool layer2Only =
            goodLayer1Indices.empty() && !goodLayer2Indices.empty();
        if (layer1Only || layer2Only) {
          const std::vector<size_t> &indices =
              layer1Only ? goodLayer1Indices : goodLayer2Indices;
          TH2D &diagnostic = layer1Only ? hSingleLayer1TimingVsXResidual
                                       : hSingleLayer2TimingVsXResidual;
          TH1D &selectedTiming = layer1Only ? hSelectedSingleLayer1Timing
                                            : hSelectedSingleLayer2Timing;
          if (layer1Only)
            ++pairlessSingleLayer1EventCount;
          else
            ++pairlessSingleLayer2EventCount;

          std::vector<size_t> reconstructedSelectedIndices;
          for (size_t i : indices) {
            const double projectedECalX =
                *ecalX * pulseZ[i] / kECalZFromTargetM;
            const double xResidual = correctedX[i] - projectedECalX;
            const double timingResidual = ecalResidual[i];
            if (!std::isfinite(xResidual) || !std::isfinite(timingResidual))
              continue;
            diagnostic.Fill(xResidual, timingResidual);
            const double normalizedX =
                (xResidual - singleResidualCenterM) / singleResidualScaleM;
            const double normalizedTiming =
                (timingResidual - singleTimingCenterNs) /
                singleTimingScaleNs;
            if (normalizedX * normalizedX +
                    normalizedTiming * normalizedTiming <=
                singleCutRadius * singleCutRadius) {
              reconstructedSelectedIndices.push_back(i);
            }
          }

          std::vector<size_t> selectedIndices = reconstructedSelectedIndices;
          if (useAnalyzerSingleCandidates) {
            selectedIndices.clear();
            const size_t nNativeSingles = std::min(
                {singlePulseIndex->GetSize(), singleLayer->GetSize(),
                 singleECalResidual->GetSize(), singleXResidual->GetSize()});
            if (nNativeSingles > 0)
              ++analyzerSingleCandidateEventCount;
            for (size_t candidate = 0; candidate < nNativeSingles;
                 ++candidate) {
              if (!std::isfinite((*singlePulseIndex)[candidate]) ||
                  !std::isfinite((*singleLayer)[candidate]))
                continue;
              const Long64_t sourceIndex =
                  std::llround((*singlePulseIndex)[candidate]);
              const int nativeLayer =
                  int(std::lround((*singleLayer)[candidate]));
              if (sourceIndex < 0 ||
                  sourceIndex >= static_cast<Long64_t>(nPulses) ||
                  nativeLayer != (layer1Only ? 0 : 1))
                continue;
              selectedIndices.push_back(static_cast<size_t>(sourceIndex));
            }
            std::sort(selectedIndices.begin(), selectedIndices.end());
            std::sort(reconstructedSelectedIndices.begin(),
                      reconstructedSelectedIndices.end());
            if (selectedIndices.size() != reconstructedSelectedIndices.size())
              ++singleCandidateCountMismatch;
            else if (selectedIndices != reconstructedSelectedIndices)
              ++singleCandidateIdentityMismatch;
          }

          for (size_t i : selectedIndices) {
            selectedTiming.Fill(ecalResidual[i]);
            if (layer1Only)
              ++recoveredSingleLayer1PulseCount;
            else
              ++recoveredSingleLayer2PulseCount;
          }
          if (!selectedIndices.empty()) {
            if (layer1Only)
              ++recoveredSingleLayer1EventCount;
            else
              ++recoveredSingleLayer2EventCount;
          }
        }
      }

      if (roiStatus) {
        const int expectedStatus = nPairs > 0 ? 1 :
            (useAnalyzerSingleCandidates && singlePulseIndex->GetSize() > 0
                 ? int(std::lround((*singleLayer)[0])) + 2
                 : 0);
        if (int(std::lround(**roiStatus)) != expectedStatus)
          ++roiStatusMismatch;
      }
      const bool hasAcceptedPair = !pairedPulseIndices.empty();
      bool hasDetectorBaseline = false;
      bool hasDetectorProjected = false;
      for (size_t i = 0; i < nPulses; ++i) {
        if (!(calibValid[i] > 0.5 && ecalEligible[i] > 0.5) ||
            !std::isfinite(pixel[i]) || !std::isfinite(ecalResidual[i]) ||
            !std::isfinite(tot[i]))
          continue;
        const int pixelID = int(std::lround(pixel[i]));
        if (pixelID < 0 || pixelID >= kCDetNPixels)
          continue;
        if (hasAcceptedPair) {
          hDetectorECalCut.Fill(ecalResidual[i]);
          hasDetectorBaseline = true;
          if (spatialPass[i] > 0.5) {
            hDetectorProjected.Fill(ecalResidual[i]);
            hasDetectorProjected = true;
          }
          if (trajectoryTimeSelectedPulseIndices.count(i) != 0U) {
            hDetectorProjectedQuality.Fill(ecalResidual[i]);
            hDetectorProjectedQualityVsToT.Fill(tot[i], ecalResidual[i]);
          }
        }
        if (pixelID < bar30FirstPixel ||
            pixelID >= bar30FirstPixel + kPixelsPerBar)
          continue;
        hBar30ECalCut.Fill(ecalResidual[i]);
        if (!(spatialPass[i] > 0.5))
          continue;
        hBar30Projected.Fill(ecalResidual[i]);
        if (!(broadQuality[i] > 0.5))
          continue;
        hBar30ProjectedQuality.Fill(ecalResidual[i]);
        hBar30ProjectedQualityVsToT.Fill(tot[i], ecalResidual[i]);
      }
      if (hasDetectorBaseline)
        ++detectorBaselineEventCount;
      if (hasDetectorProjected)
        ++detectorProjectedEventCount;
      if (!trajectoryTimeSelectedPulseIndices.empty())
        ++detectorProjectedQualityEventCount;
    }

    std::vector<size_t> acceptedIndices;
    bool hasLayer1 = false;
    bool hasLayer2 = false;
    for (size_t i = 0; i < nPulses; ++i) {
      if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
            ecalEligible[i] > 0.5 && spatialPass[i] > 0.5))
        continue;
      if (!std::isfinite(pixel[i]) || !std::isfinite(le[i]) ||
          !std::isfinite(te[i]) || !std::isfinite(tot[i]) ||
          !std::isfinite(ecalResidual[i]))
        continue;

      const int pixelID = int(std::lround(pixel[i]));
      if (pixelID < 0 || pixelID >= kCDetNPixels)
        continue;
      acceptedIndices.push_back(i);
      hasLayer1 |= pixelID < kCDetNPixels / 2;
      hasLayer2 |= pixelID >= kCDetNPixels / 2;
    }

    if (!(hasLayer1 && hasLayer2)) {
      hGoodMultiplicity.Fill(0);
      continue;
    }

    std::unordered_set<int> legacyAcceptedIDs;
    if (recoveredOnly) {
      const size_t nLegacy = std::min(
          {legacyPMT.GetSize(), legacyLE.GetSize(), legacyToT.GetSize(),
           legacyX.GetSize(), legacyY.GetSize(), legacyZ.GetSize(),
           legacyMultiplicity.GetSize()});
      std::vector<bool> legacyBasic(nLegacy, false);
      int legacyLayer1 = 0;
      int legacyLayer2 = 0;
      for (size_t i = 0; i < nLegacy; ++i) {
        const int id = static_cast<int>(legacyPMT[i]);
        const int layer = id < kCDetNPixels / 2 ? 0 : 1;
        const double correctedX = legacyX[i] * 1.08 - 0.03;
        const bool pass = *ecalY > -1.2 && *ecalY < 1.2 &&
            *ecalX > -1.5 && *ecalX < 1.5 && *ecalX != 0.0 &&
            *ecalY != 0.0 && legacyLE[i] * 0.01 >= 0.02 &&
            legacyLE[i] * 0.01 <= 60.0 && legacyToT[i] * 0.01 >= 4.0 &&
            legacyToT[i] * 0.01 <= 30.0 && legacyMultiplicity[i] < 100.0 &&
            std::fabs(correctedX - (*ecalX) * legacyZ[i] / 6.144) <= 0.08 &&
            std::fabs(legacyY[i] - (*ecalY) * legacyZ[i] / 6.144 - 0.10) <=
                0.36 &&
            *ecalTime > 10.0 && *ecalTime < 35.0;
        legacyBasic[i] = pass;
        if (pass) {
          if (layer == 0)
            ++legacyLayer1;
          else
            ++legacyLayer2;
        }
      }
      if (legacyLayer1 >= 1 && legacyLayer1 <= 100 && legacyLayer2 >= 1 &&
          legacyLayer2 <= 100) {
        for (size_t i = 0; i < nLegacy; ++i) {
          if (legacyBasic[i])
            legacyAcceptedIDs.insert(static_cast<int>(legacyPMT[i]));
        }
      }
    }

    ++twoLayerEventCount;
    int plottedThisEvent = 0;
    for (size_t i : acceptedIndices) {
      const int pixelID = int(std::lround(pixel[i]));
      if (recoveredOnly && legacyAcceptedIDs.count(pixelID))
        continue;
      const int barID = pixelID / kPixelsPerBar;
      ++plottedThisEvent;
      ++goodPulseCount;
      hGoodLE.Fill(le[i]);
      hGoodTE.Fill(te[i]);
      hGoodToT.Fill(tot[i]);
      hGoodPixel.Fill(pixelID);
      hGoodBar.Fill(barID);
      pixelLE[pixelID]->Fill(le[i]);
      barLE[barID]->Fill(le[i]);
    }
    if (plottedThisEvent > 0)
      ++plottedEventCount;
    hGoodMultiplicity.Fill(plottedThisEvent);
  }

  if (gSystem->mkdir(outputDirectory, true) != 0 &&
      gSystem->AccessPathName(outputDirectory)) {
    std::cerr << "[good-pulse TDC] Cannot create output directory "
              << outputDirectory << std::endl;
    return;
  }

  gStyle->SetOptStat(1110);
  TCanvas cAll("cCDetGoodPulseAllTDC", "Good pulse candidate timing", 1500,
               900);
  cAll.Divide(2, 2);
  cAll.cd(1);
  hGoodLE.Draw();
  cAll.cd(2);
  hGoodTE.Draw();
  cAll.cd(3);
  hGoodToT.Draw();
  cAll.cd(4);
  hGoodMultiplicity.Draw();
  cAll.SaveAs(Form("%s/CDetGoodPulse_AllTDC.pdf", outputDirectory));

  TCanvas cChannels("cCDetGoodPulseChannels", "Good pulse channels", 1500,
                    700);
  cChannels.Divide(2, 1);
  cChannels.cd(1);
  hGoodPixel.Draw();
  cChannels.cd(2);
  hGoodBar.Draw();
  cChannels.SaveAs(Form("%s/CDetGoodPulse_AllChannels.pdf", outputDirectory));

  auto drawBar30Peak = [](TH1D &hist, const char *fitName) {
    hist.SetLineColor(kBlack);
    hist.SetLineWidth(2);
    hist.Draw();
    if (hist.GetEntries() < 20)
      return;
    const double fitMin = -55.0;
    const double fitMax = -10.0;
    int peakBin = hist.FindBin(-40.0);
    for (int bin = peakBin + 1; bin <= hist.FindBin(-15.0); ++bin) {
      if (hist.GetBinContent(bin) > hist.GetBinContent(peakBin))
        peakBin = bin;
    }
    const double peak = hist.GetBinCenter(peakBin);
    const double background =
        0.5 * (hist.GetBinContent(hist.FindBin(fitMin)) +
               hist.GetBinContent(hist.FindBin(fitMax)));
    TF1 fit(fitName, "gaus(0)+pol1(3)", fitMin, fitMax);
    fit.SetParameters(std::max(1.0, hist.GetBinContent(peakBin) - background),
                      peak, 3.0, background, 0.0);
    fit.SetLineColor(kRed + 1);
    hist.Fit(&fit, "RQN");
    fit.DrawCopy("same");
  };

  TCanvas cBar30("cCDetGoodPulseBar30Amalgamated",
                 "CDet bar 30 amalgamated ECal-CDet timing", 1200, 800);
  cBar30.Divide(2, 2);
  cBar30.cd(1);
  drawBar30Peak(hBar30ECalCut, "fCDetBar30ECalCut");
  cBar30.cd(2);
  drawBar30Peak(hBar30Projected, "fCDetBar30Projected");
  cBar30.cd(3);
  drawBar30Peak(hBar30ProjectedQuality, "fCDetBar30ProjectedQuality");
  cBar30.cd(4);
  hBar30ProjectedQualityVsToT.SetStats(false);
  hBar30ProjectedQualityVsToT.Draw("COLZ");
  cBar30.SaveAs(Form("%s/CDetGoodPulse_Bar030_Amalgamated.pdf",
                     outputDirectory));
  cBar30.SaveAs(Form("%s/CDetGoodPulse_Bar030_Amalgamated.png",
                     outputDirectory));

  TCanvas cDetector("cCDetGoodPulseDetectorAmalgamated",
                    "Whole-detector amalgamated ECal-CDet timing", 1200, 800);
  cDetector.Divide(2, 2);
  cDetector.cd(1);
  drawBar30Peak(hDetectorECalCut, "fCDetDetectorECalCut");
  cDetector.cd(2);
  drawBar30Peak(hDetectorProjected, "fCDetDetectorProjected");
  cDetector.cd(3);
  drawBar30Peak(hSelectedPairMeanResidual,
                "fCDetSelectedPairMeanResidual");
  cDetector.cd(4);
  hSelectedPairMeanResidualVsMeanToT.SetStats(false);
  hSelectedPairMeanResidualVsMeanToT.Draw("COLZ");
  cDetector.SaveAs(Form("%s/CDetGoodPulse_Detector_Amalgamated.pdf",
                        outputDirectory));
  cDetector.SaveAs(Form("%s/CDetGoodPulse_Detector_Amalgamated.png",
                        outputDirectory));

  // Unlike the ECal-minus-CDet residual, the absolute corrected pair-mean
  // time isolates movement of the CDet timing origin from movement of the
  // ECal centroid. Fit a narrow local core around the modal bin and display
  // the conventional 30 ns timing reference. This trajectory-selected sample
  // is not identical to the projected-half-bar shift-calibration sample, so
  // compare like-for-like reference runs before changing shift_ns.
  int pairMeanModeBin = hSelectedPairMeanCorrectedLE.GetMaximumBin();
  const double pairMeanMode =
      hSelectedPairMeanCorrectedLE.GetBinCenter(pairMeanModeBin);
  const double pairMeanFitMin = std::max(leMinNs, pairMeanMode - 5.0);
  const double pairMeanFitMax = std::min(leMaxNs, pairMeanMode + 5.0);
  const double pairMeanBackground = 0.5 *
      (hSelectedPairMeanCorrectedLE.GetBinContent(
           hSelectedPairMeanCorrectedLE.FindBin(pairMeanFitMin)) +
       hSelectedPairMeanCorrectedLE.GetBinContent(
           hSelectedPairMeanCorrectedLE.FindBin(pairMeanFitMax)));
  TF1 selectedPairMeanFit("fCDetSelectedPairMeanCorrectedLE",
                          "gaus(0)+pol1(3)", pairMeanFitMin,
                          pairMeanFitMax);
  selectedPairMeanFit.SetParameters(
      std::max(1.0, hSelectedPairMeanCorrectedLE.GetMaximum() -
                        pairMeanBackground),
      pairMeanMode, 3.0, pairMeanBackground, 0.0);
  selectedPairMeanFit.SetParLimits(2, 0.2, 10.0);
  selectedPairMeanFit.SetLineColor(kRed + 1);
  int selectedPairMeanFitStatus = -1;
  if (hSelectedPairMeanCorrectedLE.GetEntries() >= 20)
    selectedPairMeanFitStatus = static_cast<int>(
        hSelectedPairMeanCorrectedLE.Fit(&selectedPairMeanFit, "RQ0"));

  TCanvas cSelectedPairMeanTime(
      "cCDetGoodPulseSelectedPairMeanTime",
      "Absolute corrected CDet selected-pair mean time", 1400, 650);
  cSelectedPairMeanTime.Divide(2, 1);
  for (int pad = 1; pad <= 2; ++pad) {
    cSelectedPairMeanTime.cd(pad);
    if (pad == 2) gPad->SetLogy();
    hSelectedPairMeanCorrectedLE.SetLineColor(kBlack);
    hSelectedPairMeanCorrectedLE.SetLineWidth(2);
    hSelectedPairMeanCorrectedLE.Draw("HIST");
    if (selectedPairMeanFitStatus == 0) selectedPairMeanFit.Draw("same");
    gPad->Update();
    TLine targetLine(30.0, gPad->GetUymin(), 30.0, gPad->GetUymax());
    targetLine.SetLineColor(kBlue + 1);
    targetLine.SetLineStyle(2);
    targetLine.SetLineWidth(3);
    targetLine.DrawClone("same");
  }
  cSelectedPairMeanTime.SaveAs(Form(
      "%s/CDetGoodPulse_SelectedPairMeanCorrectedLE.pdf", outputDirectory));
  cSelectedPairMeanTime.SaveAs(Form(
      "%s/CDetGoodPulse_SelectedPairMeanCorrectedLE.png", outputDirectory));

  TCanvas cPairSlope("cCDetGoodPulsePairSlopeDiagnostics",
                     "CDet pair trajectory diagnostics", 1200, 500);
  cPairSlope.Divide(3, 1);
  cPairSlope.cd(1);
  hPairTrajectoryResidual.Draw();
  cPairSlope.cd(2);
  hBestPairTrajectoryResidual.Draw();
  cPairSlope.cd(3);
  hPairTimingVsTrajectoryResidual.SetStats(false);
  hPairTimingVsTrajectoryResidual.Draw("COLZ");
  TEllipse pairSelectionEllipse(pairResidualCenterM, pairTimingCenterNs,
                                pairResidualScaleM * pairCutRadius,
                                pairTimingScaleNs * pairCutRadius);
  pairSelectionEllipse.SetFillStyle(0);
  pairSelectionEllipse.SetLineColor(kRed + 1);
  pairSelectionEllipse.SetLineWidth(3);
  pairSelectionEllipse.Draw("same");
  cPairSlope.SaveAs(Form("%s/CDetGoodPulse_PairTrajectoryDiagnostics.pdf",
                         outputDirectory));
  cPairSlope.SaveAs(Form("%s/CDetGoodPulse_PairTrajectoryDiagnostics.png",
                         outputDirectory));

  TCanvas cSingleLayerRecovery(
      "cCDetGoodPulseSingleLayerRecovery",
      "Pairless single-layer good-pulse recovery", 1500, 1100);
  cSingleLayerRecovery.Divide(2, 2);
  cSingleLayerRecovery.cd(1);
  gPad->SetRightMargin(0.14);
  hSingleLayer1TimingVsXResidual.SetStats(false);
  hSingleLayer1TimingVsXResidual.Draw("COLZ");
  TEllipse singleLayer1Ellipse(singleResidualCenterM,
                               singleTimingCenterNs,
                               singleResidualScaleM * singleCutRadius,
                               singleTimingScaleNs * singleCutRadius);
  singleLayer1Ellipse.SetFillStyle(0);
  singleLayer1Ellipse.SetLineColor(kRed + 1);
  singleLayer1Ellipse.SetLineWidth(3);
  singleLayer1Ellipse.Draw("same");
  cSingleLayerRecovery.cd(2);
  gPad->SetRightMargin(0.14);
  hSingleLayer2TimingVsXResidual.SetStats(false);
  hSingleLayer2TimingVsXResidual.Draw("COLZ");
  TEllipse singleLayer2Ellipse(singleResidualCenterM,
                               singleTimingCenterNs,
                               singleResidualScaleM * singleCutRadius,
                               singleTimingScaleNs * singleCutRadius);
  singleLayer2Ellipse.SetFillStyle(0);
  singleLayer2Ellipse.SetLineColor(kRed + 1);
  singleLayer2Ellipse.SetLineWidth(3);
  singleLayer2Ellipse.Draw("same");
  cSingleLayerRecovery.cd(3);
  hSelectedSingleLayer1Timing.SetLineWidth(2);
  hSelectedSingleLayer1Timing.Draw("HIST");
  cSingleLayerRecovery.cd(4);
  hSelectedSingleLayer2Timing.SetLineWidth(2);
  hSelectedSingleLayer2Timing.Draw("HIST");
  cSingleLayerRecovery.SaveAs(
      Form("%s/CDetGoodPulse_SingleLayerRecovery.pdf", outputDirectory));
  cSingleLayerRecovery.SaveAs(
      Form("%s/CDetGoodPulse_SingleLayerRecovery.png", outputDirectory));

  TCanvas cEventOutcome("cCDetGoodPulseECalAdmittedEventOutcome",
                        "CDet event-selection outcome", 1200, 750);
  cEventOutcome.SetBottomMargin(0.20);
  hECalAdmittedEventOutcome.SetStats(false);
  hECalAdmittedEventOutcome.SetFillColor(kAzure - 4);
  hECalAdmittedEventOutcome.SetLineColor(kBlue + 2);
  hECalAdmittedEventOutcome.SetLineWidth(2);
  hECalAdmittedEventOutcome.GetXaxis()->SetLabelSize(0.040);
  hECalAdmittedEventOutcome.Draw("HIST TEXT0");
  cEventOutcome.SaveAs(
      Form("%s/CDetGoodPulse_ECalAdmittedEventOutcome.pdf", outputDirectory));
  cEventOutcome.SaveAs(
      Form("%s/CDetGoodPulse_ECalAdmittedEventOutcome.png", outputDirectory));

  TCanvas cPairFailure("cCDetGoodPulsePairFailureDiagnostics",
                       "Why good two-layer events fail pair selection",
                       1800, 650);
  cPairFailure.Divide(3, 1);
  cPairFailure.cd(1);
  gPad->SetRightMargin(0.14);
  hFailedPairBestTimingVsTrajectory.SetStats(false);
  hFailedPairBestTimingVsTrajectory.Draw("COLZ");
  TEllipse failedPairEllipse(pairResidualCenterM, pairTimingCenterNs,
                             pairResidualScaleM * pairCutRadius,
                             pairTimingScaleNs * pairCutRadius);
  failedPairEllipse.SetFillStyle(0);
  failedPairEllipse.SetLineColor(kRed + 1);
  failedPairEllipse.SetLineWidth(3);
  failedPairEllipse.Draw("same");
  cPairFailure.cd(2);
  hFailedPairMinimumRadius.SetLineWidth(2);
  hFailedPairMinimumRadius.Draw("HIST");
  TLine pairRadiusCut(pairCutRadius, 0.0, pairCutRadius,
                      1.05 * hFailedPairMinimumRadius.GetMaximum());
  pairRadiusCut.SetLineColor(kRed + 1);
  pairRadiusCut.SetLineWidth(3);
  pairRadiusCut.Draw("same");
  cPairFailure.cd(3);
  gPad->SetBottomMargin(0.20);
  hNoStoredPairFailureReason.SetStats(false);
  hNoStoredPairFailureReason.SetFillColor(kOrange - 3);
  hNoStoredPairFailureReason.SetLineColor(kOrange + 7);
  hNoStoredPairFailureReason.SetLineWidth(2);
  hNoStoredPairFailureReason.GetXaxis()->SetLabelSize(0.045);
  hNoStoredPairFailureReason.Draw("HIST TEXT0");
  cPairFailure.SaveAs(
      Form("%s/CDetGoodPulse_PairFailureDiagnostics.pdf", outputDirectory));
  cPairFailure.SaveAs(
      Form("%s/CDetGoodPulse_PairFailureDiagnostics.png", outputDirectory));

  TCanvas cAllCombinationRecovery(
      "cCDetGoodPulseAllCombinationRecovery",
      "All-combination recovery before greedy pairing", 1900, 650);
  cAllCombinationRecovery.Divide(3, 1);
  cAllCombinationRecovery.cd(1);
  gPad->SetBottomMargin(0.32);
  hAllCombinationRecoveryOutcome.SetStats(false);
  hAllCombinationRecoveryOutcome.SetMinimum(0.0);
  hAllCombinationRecoveryOutcome.SetFillColor(kGreen - 6);
  hAllCombinationRecoveryOutcome.SetLineColor(kGreen + 3);
  hAllCombinationRecoveryOutcome.SetLineWidth(2);
  hAllCombinationRecoveryOutcome.GetXaxis()->SetLabelSize(0.043);
  hAllCombinationRecoveryOutcome.GetXaxis()->LabelsOption("v");
  hAllCombinationRecoveryOutcome.Draw("HIST TEXT0");
  cAllCombinationRecovery.cd(2);
  gPad->SetBottomMargin(0.20);
  hEllipsePairHardGateFailureReason.SetStats(false);
  hEllipsePairHardGateFailureReason.SetMinimum(0.0);
  hEllipsePairHardGateFailureReason.SetFillColor(kOrange - 3);
  hEllipsePairHardGateFailureReason.SetLineColor(kOrange + 7);
  hEllipsePairHardGateFailureReason.SetLineWidth(2);
  hEllipsePairHardGateFailureReason.GetXaxis()->SetLabelSize(0.043);
  hEllipsePairHardGateFailureReason.Draw("HIST TEXT0");
  cAllCombinationRecovery.cd(3);
  hAlternativePairMultiplicity.SetLineWidth(2);
  hAlternativePairMultiplicity.Draw("HIST");
  cAllCombinationRecovery.SaveAs(Form(
      "%s/CDetGoodPulse_AllCombinationRecovery.pdf", outputDirectory));
  cAllCombinationRecovery.SaveAs(Form(
      "%s/CDetGoodPulse_AllCombinationRecovery.png", outputDirectory));

  TCanvas cRecoveredGreedy(
      "cCDetGoodPulseRecoveredGreedyPairs",
      "Best ECal-informed pairs recovered after greedy loss", 1500, 1100);
  cRecoveredGreedy.Divide(2, 2);
  cRecoveredGreedy.cd(1);
  gPad->SetRightMargin(0.14);
  hRecoveredGreedyBestTimingVsTrajectory.SetStats(false);
  hRecoveredGreedyBestTimingVsTrajectory.Draw("COLZ");
  TEllipse recoveredGreedyEllipse(
      pairResidualCenterM, pairTimingCenterNs,
      pairResidualScaleM * pairCutRadius,
      pairTimingScaleNs * pairCutRadius);
  recoveredGreedyEllipse.SetFillStyle(0);
  recoveredGreedyEllipse.SetLineColor(kRed + 1);
  recoveredGreedyEllipse.SetLineWidth(3);
  recoveredGreedyEllipse.Draw("same");
  cRecoveredGreedy.cd(2);
  hRecoveredGreedyBestMeanResidual.SetLineWidth(2);
  hRecoveredGreedyBestMeanResidual.Draw("HIST");
  cRecoveredGreedy.cd(3);
  gPad->SetRightMargin(0.14);
  hRecoveredGreedyBestXCorrelation.SetStats(false);
  hRecoveredGreedyBestXCorrelation.Draw("COLZ");
  TLine recoveredGreedyXDiagonal(-1.6, -1.6, 1.6, 1.6);
  recoveredGreedyXDiagonal.SetLineColor(kRed + 1);
  recoveredGreedyXDiagonal.SetLineWidth(3);
  recoveredGreedyXDiagonal.Draw("same");
  cRecoveredGreedy.cd(4);
  hRecoveredGreedyBestXResidual.SetLineWidth(2);
  hRecoveredGreedyBestXResidual.Draw("HIST");
  TLine recoveredGreedyXZero(
      0.0, 0.0, 0.0, 1.05 * hRecoveredGreedyBestXResidual.GetMaximum());
  recoveredGreedyXZero.SetLineColor(kRed + 1);
  recoveredGreedyXZero.SetLineWidth(3);
  recoveredGreedyXZero.Draw("same");
  cRecoveredGreedy.SaveAs(Form(
      "%s/CDetGoodPulse_RecoveredGreedyPairs.pdf", outputDirectory));
  cRecoveredGreedy.SaveAs(Form(
      "%s/CDetGoodPulse_RecoveredGreedyPairs.png", outputDirectory));

  TCanvas cPairX("cCDetGoodPulseSelectedPairXCorrelation",
                 "Selected CDet pair x versus projected ECal x", 800, 700);
  cPairX.SetRightMargin(0.14);
  hSelectedPairXCorrelation.SetStats(false);
  hSelectedPairXCorrelation.Draw("COLZ");
  TLine pairXDiagonal(-1.6, -1.6, 1.6, 1.6);
  pairXDiagonal.SetLineColor(kRed + 1);
  pairXDiagonal.SetLineWidth(3);
  pairXDiagonal.Draw("same");
  cPairX.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXCorrelation.pdf",
                     outputDirectory));
  cPairX.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXCorrelation.png",
                     outputDirectory));

  TCanvas cPairXResidual("cCDetGoodPulseSelectedPairXResidual",
                         "Selected CDet pair x residual", 800, 700);
  hSelectedPairXResidual.SetLineWidth(2);
  hSelectedPairXResidual.Draw("HIST");
  TLine pairXResidualZero(0.0, 0.0, 0.0,
                          1.05 * hSelectedPairXResidual.GetMaximum());
  pairXResidualZero.SetLineColor(kRed + 1);
  pairXResidualZero.SetLineWidth(3);
  pairXResidualZero.Draw("same");
  cPairXResidual.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXResidual.pdf",
                             outputDirectory));
  cPairXResidual.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXResidual.png",
                             outputDirectory));

  TCanvas cPairXDiagnostics("cCDetGoodPulseSelectedPairXDiagnostics",
                            "Selected CDet pair x diagnostics", 1500, 650);
  cPairXDiagnostics.Divide(2, 1);
  cPairXDiagnostics.cd(1);
  gPad->SetRightMargin(0.14);
  hSelectedPairXCorrelation.Draw("COLZ");
  pairXDiagonal.Draw("same");
  cPairXDiagnostics.cd(2);
  hSelectedPairXResidual.Draw("HIST");
  pairXResidualZero.Draw("same");
  cPairXDiagnostics.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXDiagnostics.pdf",
                                outputDirectory));
  cPairXDiagnostics.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXDiagnostics.png",
                                outputDirectory));

  const char *pageNames[4] = {"Layer1_Left", "Layer1_Right",
                              "Layer2_Left", "Layer2_Right"};
  TCanvas cBars("cCDetGoodPulseBars", "Good pulse bar timing", 1800, 1200);
  for (int page = 0; page < 4; ++page) {
    const int firstBar = page * 42;
    DrawBarPage(&cBars, barLE, firstBar, firstBar + 41);
    cBars.SaveAs(Form("%s/CDetGoodPulse_Bars_%s.pdf", outputDirectory,
                      pageNames[page]));
  }

  TFile output(Form("%s/CDetGoodPulse_AllTDC.root", outputDirectory),
               "RECREATE");
  hGoodLE.Write();
  hGoodTE.Write();
  hGoodToT.Write();
  hGoodPixel.Write();
  hGoodBar.Write();
  hGoodMultiplicity.Write();
  hBar30ECalCut.Write();
  hBar30Projected.Write();
  hBar30ProjectedQuality.Write();
  hBar30ProjectedQualityVsToT.Write();
  hDetectorECalCut.Write();
  hDetectorProjected.Write();
  hDetectorProjectedQuality.Write();
  hDetectorProjectedQualityVsToT.Write();
  hSelectedPairMeanResidual.Write();
  hSelectedPairMeanCorrectedLE.Write();
  hSelectedPairMeanResidualVsMeanToT.Write();
  hSelectedPairXCorrelation.Write();
  hSelectedPairXResidual.Write();
  hPairTrajectoryResidual.Write();
  hBestPairTrajectoryResidual.Write();
  hPairTimingVsTrajectoryResidual.Write();
  hSingleLayer1TimingVsXResidual.Write();
  hSingleLayer2TimingVsXResidual.Write();
  hSelectedSingleLayer1Timing.Write();
  hSelectedSingleLayer2Timing.Write();
  hECalAdmittedEventOutcome.Write();
  hFailedPairBestTimingVsTrajectory.Write();
  hFailedPairMinimumRadius.Write();
  hNoStoredPairFailureReason.Write();
  hAllCombinationRecoveryOutcome.Write();
  hEllipsePairHardGateFailureReason.Write();
  hAlternativePairMultiplicity.Write();
  hRecoveredGreedyBestTimingVsTrajectory.Write();
  hRecoveredGreedyBestMeanResidual.Write();
  hRecoveredGreedyBestXCorrelation.Write();
  hRecoveredGreedyBestXResidual.Write();
  for (TH1D *hist : pixelLE)
    hist->Write();
  for (TH1D *hist : barLE)
    hist->Write();
  output.Close();

  std::cout << "\n[good-pulse TDC] Files/events/pulses: " << filesAdded << "/"
            << eventCount << "/" << pulseCount << std::endl;
  std::cout << "[good-pulse TDC] Candidate source: pairs="
            << (useAnalyzerPairCandidates ? "analyzer pair_candidate.*"
                                          : "legacy pair.* fallback")
            << ", singles="
            << (useAnalyzerSingleCandidates ? "analyzer single_candidate.*"
                                             : "macro reconstruction fallback")
            << std::endl;
  std::cout << "[good-pulse TDC] Accepted calibrated good pulse candidates: "
            << goodPulseCount << std::endl;
  std::cout << "[good-pulse TDC] Population mode: "
            << (recoveredOnly ? "recovered-only" : "all new candidates")
            << std::endl;
  std::cout << "[good-pulse TDC] Events with accepted candidates in both layers: "
            << twoLayerEventCount << std::endl;
  std::cout << "[good-pulse TDC] Events contributing to plotted population: "
            << plottedEventCount << std::endl;
  std::cout << "[good-pulse TDC] Events with inconsistent pulse-array sizes: "
            << malformedEvents << std::endl;
  std::cout << "[good-pulse TDC] ECal-energy-selected events: "
            << ecalEnergyEventCount << std::endl;
  std::cout << "[good-pulse TDC] ECal-energy events with >=1 accepted Layer-1/Layer-2 pair: "
            << detectorBaselineEventCount << std::endl;
  std::cout << "[good-pulse TDC] Pair-containing events with >=1 projected CDet pulse: "
            << detectorProjectedEventCount << std::endl;
  std::cout << "[good-pulse TDC] ECal-energy events contributing trajectory-time selected pair members: "
            << detectorProjectedQualityEventCount << std::endl;
  std::cout << "[good-pulse TDC] "
            << (useAnalyzerPairCandidates
                    ? "Accepted Layer-1/Layer-2 pair hypotheses"
                    : "Accepted one-to-one Layer-1/Layer-2 pairs")
            << " in ECal-energy events: "
            << acceptedPairCount << std::endl;
  if (candidateGreedySelected) {
    std::cout << "[good-pulse TDC] Analyzer pair-classification regression: "
              << "count-mismatch events=" << pairClassificationCountMismatch
              << ", identity mismatches="
              << pairClassificationIdentityMismatch << std::endl;
  }
  std::cout << "[good-pulse TDC] Trajectory-time ellipse: center=("
            << pairResidualCenterM << " m, " << pairTimingCenterNs
            << " ns), scales=(" << pairResidualScaleM << " m, "
            << pairTimingScaleNs << " ns), radius=" << pairCutRadius
            << std::endl;
  std::cout << "[good-pulse TDC] Pairs inside trajectory-time ellipse: "
            << trajectoryTimeSelectedPairCount << std::endl;
  std::cout << "[good-pulse TDC] Selected pair-mean timing entries: "
            << static_cast<Long64_t>(hSelectedPairMeanResidual.GetEntries())
            << std::endl;
  std::cout << "[good-pulse TDC] Absolute selected-pair mean CDet time: "
            << "mode=" << pairMeanMode << " ns, local-fit status="
            << selectedPairMeanFitStatus;
  if (selectedPairMeanFitStatus == 0) {
    std::cout << ", peak=" << selectedPairMeanFit.GetParameter(1)
              << " +/- " << selectedPairMeanFit.GetParError(1)
              << " ns, sigma=" << selectedPairMeanFit.GetParameter(2)
              << " +/- " << selectedPairMeanFit.GetParError(2) << " ns"
              << ", chi2/NDF=" << selectedPairMeanFit.GetChisquare() << "/"
              << selectedPairMeanFit.GetNDF();
  }
  std::cout << std::endl;
  std::cout << "[good-pulse TDC] Single-layer recovery ellipse: center=("
            << singleResidualCenterM << " m, " << singleTimingCenterNs
            << " ns), scales=(" << singleResidualScaleM << " m, "
            << singleTimingScaleNs << " ns), radius=" << singleCutRadius
            << std::endl;
  std::cout << "[good-pulse TDC] Pairless Layer-1-only events / recovered events / selected pulses: "
            << pairlessSingleLayer1EventCount << " / "
            << recoveredSingleLayer1EventCount << " / "
            << recoveredSingleLayer1PulseCount << std::endl;
  std::cout << "[good-pulse TDC] Pairless Layer-2-only events / recovered events / selected pulses: "
            << pairlessSingleLayer2EventCount << " / "
            << recoveredSingleLayer2EventCount << " / "
            << recoveredSingleLayer2PulseCount << std::endl;
  if (useAnalyzerSingleCandidates) {
    if (analyzerSingleCandidateEventCount > 0) {
      std::cout << "[good-pulse TDC] Analyzer single-candidate regression: events="
                << analyzerSingleCandidateEventCount
                << ", count mismatches=" << singleCandidateCountMismatch
                << ", identity mismatches=" << singleCandidateIdentityMismatch
                << std::endl;
    } else {
      std::cout << "[good-pulse TDC] Analyzer single_candidate.* is present "
                   "but empty; regression is not applicable (selection may "
                   "be disabled by database policy)."
                << std::endl;
    }
  }
  if (roiStatus)
    std::cout << "[good-pulse TDC] ROI-status consistency mismatches: "
              << roiStatusMismatch << std::endl;
  std::cout << "[good-pulse TDC] ECal-admitted event outcome denominator ("
            << ecalTimeMinNs << " < t_ECal < " << ecalTimeMaxNs
            << " ns): " << ecalAdmittedEventCount << std::endl;
  std::cout << "[good-pulse TDC]   selected pair: "
            << selectedPairEventCount << std::endl;
  std::cout << "[good-pulse TDC]   good pulses in both layers, pair selection failed: "
            << goodBothPairFailureEventCount << " (no stored pair "
            << goodBothNoStoredPairEventCount << ", stored pair outside ellipse "
            << goodBothEllipseFailureEventCount << ")" << std::endl;
  std::cout << "[good-pulse TDC]     best stored-pair failure: timing axis only "
            << failedPairTimingOnlyEventCount << ", trajectory axis only "
            << failedPairTrajectoryOnlyEventCount << ", both axes "
            << failedPairBothAxesEventCount << ", ellipse corner "
            << failedPairEllipseCornerEventCount << std::endl;
  std::cout << "[good-pulse TDC]     no-stored-pair nearest-combination gate failures: dt "
            << static_cast<Long64_t>(hNoStoredPairFailureReason.GetBinContent(1))
            << ", dx "
            << static_cast<Long64_t>(hNoStoredPairFailureReason.GetBinContent(2))
            << ", dy "
            << static_cast<Long64_t>(hNoStoredPairFailureReason.GetBinContent(3))
            << ", dt+dx "
            << static_cast<Long64_t>(hNoStoredPairFailureReason.GetBinContent(4))
            << ", dt+dy "
            << static_cast<Long64_t>(hNoStoredPairFailureReason.GetBinContent(5))
            << ", dx+dy "
            << static_cast<Long64_t>(hNoStoredPairFailureReason.GetBinContent(6))
            << ", all "
            << static_cast<Long64_t>(hNoStoredPairFailureReason.GetBinContent(7))
            << ", unexpected " << noStoredPairUnexpectedEventCount
            << std::endl;
  std::cout << "[good-pulse TDC]     all-combination decomposition: no ellipse combination "
            << noAllGoodEllipseCombinationEventCount
            << ", ellipse combination fails analyzer hard gates "
            << ellipseCombinationFailsHardGatesEventCount
            << ", hard-gate ellipse combination lost by greedy pairing "
            << greedyPairingLostEllipseCombinationEventCount << std::endl;
  std::cout << "[good-pulse TDC]     recovered hard-gate ellipse combinations: "
            << allCombinationRecoveredPairCount << std::endl;
  std::cout << "[good-pulse TDC]     best greedy-recovered pair timing mean/RMS: "
            << hRecoveredGreedyBestMeanResidual.GetMean() << " / "
            << hRecoveredGreedyBestMeanResidual.GetRMS() << " ns"
            << std::endl;
  std::cout << "[good-pulse TDC]     best greedy-recovered pair x residual mean/RMS: "
            << hRecoveredGreedyBestXResidual.GetMean() << " / "
            << hRecoveredGreedyBestXResidual.GetRMS() << " m"
            << std::endl;
  std::cout << "[good-pulse TDC]     best ellipse-combination hard-gate failures: dt "
            << static_cast<Long64_t>(hEllipsePairHardGateFailureReason.GetBinContent(1))
            << ", dx "
            << static_cast<Long64_t>(hEllipsePairHardGateFailureReason.GetBinContent(2))
            << ", dy "
            << static_cast<Long64_t>(hEllipsePairHardGateFailureReason.GetBinContent(3))
            << ", dt+dx "
            << static_cast<Long64_t>(hEllipsePairHardGateFailureReason.GetBinContent(4))
            << ", dt+dy "
            << static_cast<Long64_t>(hEllipsePairHardGateFailureReason.GetBinContent(5))
            << ", dx+dy "
            << static_cast<Long64_t>(hEllipsePairHardGateFailureReason.GetBinContent(6))
            << ", all "
            << static_cast<Long64_t>(hEllipsePairHardGateFailureReason.GetBinContent(7))
            << std::endl;
  std::cout << "[good-pulse TDC]   calibrated pulses in both layers, quality/spatial coverage lost: "
            << calibratedBothQualitySpatialLossEventCount
            << " (good Layer 1 only " << goodLayer1OnlyEventCount
            << ", good Layer 2 only " << goodLayer2OnlyEventCount
            << ", neither layer good " << noGoodLayerEventCount << ")"
            << std::endl;
  std::cout << "[good-pulse TDC]   missing calibrated pulse coverage in >=1 layer: "
            << missingCalibratedLayerEventCount << " (Layer 1 only "
            << calibratedLayer1OnlyEventCount << ", Layer 2 only "
            << calibratedLayer2OnlyEventCount << ", neither "
            << noCalibratedPulseEventCount << ")" << std::endl;
  std::cout << "[good-pulse TDC] Output directory: " << outputDirectory
            << std::endl;

  for (TH1D *hist : pixelLE)
    delete hist;
  for (TH1D *hist : barLE)
    delete hist;
}

// Configuration-file interface. The original positional interface above is
// retained unchanged (apart from appended optional energy limits), so existing
// ROOT commands remain valid.
void Plot_CDet_GoodPulseCandidates_AllTDC(
    const char *configFile, const char *inputDirectory = nullptr,
    const char *outputDirectory = nullptr) {
  CDetGoodPulseConfig::Values config;
  if (!CDetGoodPulseConfig::Load(
          configFile, config, "good-pulse TDC configuration"))
    return;

  const TString defaultOutput = TString::Format(
      "CDet_run%d_good_pulse_tdc", config.runNumber);
  const char *resolvedOutput =
      outputDirectory && outputDirectory[0] ? outputDirectory : defaultOutput.Data();
  Plot_CDet_GoodPulseCandidates_AllTDC(
      config.runNumber, inputDirectory, resolvedOutput, config.binWidthNs,
      config.leMinNs, config.leMaxNs, config.totMinNs, config.totMaxNs,
      config.recoveredOnly, config.pairResidualCenterM,
      config.pairTimingCenterNs, config.pairResidualScaleM,
      config.pairTimingScaleNs, config.pairCutRadius,
      config.singleResidualCenterM, config.singleTimingCenterNs,
      config.singleResidualScaleM, config.singleTimingScaleNs,
      config.singleCutRadius, config.ecalTimeMinNs, config.ecalTimeMaxNs,
      config.ecalEnergyMinGeV, config.ecalEnergyMaxGeV,
      config.oppositeSideEnabled, config.oppositeDYCenterM,
      config.oppositeDYToleranceM, config.oppositeProjectedYCenterM,
      config.oppositeProjectedYMaxM);
}
