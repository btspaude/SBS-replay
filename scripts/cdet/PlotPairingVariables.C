#include <TCanvas.h>
#include <TEnv.h>
#include <THashList.h>
#include <TH2D.h>
#include <TEllipse.h>
#include <set>
#include <stdexcept>
#include <string>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <array>
#include <fstream>
#include <memory>
#include <algorithm>
#include <vector>
#include <limits>
#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <iostream>

#include "CDetRunDataset.h"
#include "PairingPlotHelpers.h"

// Interactive ROOT:
// .L PlotPairingVariables.C+
// PlotPairingVariables(6077, "/path/to/Rootfiles", 10000);
// maxEvents = -1 reads the full selected dataset. Input defaults to OUT_DIR.
// Only pulse and pair branches are read. Corrected times below are already ns.
// Preferred configuration interface:
// PlotPairingVariables("PlotPairingVariables.conf");
// Optionally override only its input directory:
// PlotPairingVariables("PlotPairingVariables.conf", "/path/to/Rootfiles");
// Positional arguments after maxEvents: L1 bar, bin width, LE min/max,
// layer-dt min/max, ECal-dt min/max, member ToT min/max, maximum |L2-L1|,
// maximum ECal ellipse radius, savePlots, outputDirectory, minEntriesPerBar,
// minimum and maximum segment (both -1 means all segments), then geometry
// display limits for x1-x2, x1, x residual, and angle, followed by the ECal
// event adctime and energy cut limits.
// -1 disables either additional pair cut. Default selects stored pairs as-is.
// One event pass produces all eight canvases, including focused/all-detector
// geometry and the all-bar scan. The upper row is the existing all-pair
// population; the lower row contains one best pair per event. savePlots
// defaults to false. When true, save all eight as PNG AND PDF plus
// the per-bar CSV, under outputDirectory with run/bar-specific file names.
// These tighten the STORED pair population; they do not rerun assignment.
// The best-pair row chooses the accepted pair with the smallest
// pair.ecal_score in each event (falling back to the explicitly recomputed
// ellipse radius when that diagnostic score is nonfinite).
// Standard deviations and their ROOT moment-based error estimates include
// all finite selected values, including histogram tails. Tail counts are kept
// in the CSV/terminal diagnostics but are not drawn on the thesis canvases.
// No Gaussian fit is used. Pair/event correlations are not modeled in errors.
// Geometry uses transport x/z for the out-of-plane angle and assumes the
// trajectory starts at (y,z)=(0,0); the angle is a rough diagnostic only.
void PlotPairingVariables(int runNumber = 6077,
                         const char *inputDirectory = nullptr,
                         Long64_t maxEvents = -1, int selectedLayer1Bar = 30,
                         double binWidthNs = 0.5, double leMinNs = 0, double leMaxNs = 80,
                         double dtMinNs = -30, double dtMaxNs = 30,
                         double ecalDTMinNs = -100, double ecalDTMaxNs = 50,
                         double memberToTMinNs = 0, double memberToTMaxNs = 1e9,
                         double layerDTMaxNs = -1, double pairRadiusMax = -1,
                         bool savePlots = false, const char *outputDirectory = "pairing_plots",
                         int minEntriesPerBar = 30, int segmentMin = -1,
                         int segmentMax = -1, double xDiffMinM = -0.5,
                         double xDiffMaxM = 0.5, double x1MinM = -1.6,
                         double x1MaxM = 1.6, double xResidualMinM = -0.2,
                         double xResidualMaxM = 0.2, double angleMinDeg = -40,
                         double angleMaxDeg = 40, double angleBinWidthDeg = 0.1,
                         double ecalTimeMinNs = -10,
                         double ecalTimeMaxNs = 4, double ecalEnergyMinGeV = 3.0,
                         double ecalEnergyMaxGeV = 4.5) {
  using namespace CDetPairingPlots;
  const int nLE = Bins(binWidthNs, leMinNs, leMaxNs);
  const int nDT = Bins(binWidthNs, dtMinNs, dtMaxNs);
  const int nECalDT = Bins(binWidthNs, ecalDTMinNs, ecalDTMaxNs);
  const int nAngleBins = Bins(angleBinWidthDeg, angleMinDeg, angleMaxDeg);
  if ((savePlots && (!outputDirectory || !outputDirectory[0])) || minEntriesPerBar < 2 || !nLE || !nDT || !nECalDT ||
      !std::isfinite(xDiffMinM) || !std::isfinite(xDiffMaxM) || xDiffMaxM <= xDiffMinM ||
      !std::isfinite(x1MinM) || !std::isfinite(x1MaxM) || x1MaxM <= x1MinM ||
      !std::isfinite(xResidualMinM) || !std::isfinite(xResidualMaxM) || xResidualMaxM <= xResidualMinM ||
      !std::isfinite(angleMinDeg) || !std::isfinite(angleMaxDeg) || angleMaxDeg <= angleMinDeg ||
      !nAngleBins ||
      !std::isfinite(ecalTimeMinNs) || !std::isfinite(ecalTimeMaxNs) || ecalTimeMaxNs <= ecalTimeMinNs ||
      !std::isfinite(ecalEnergyMinGeV) || !std::isfinite(ecalEnergyMaxGeV) || ecalEnergyMaxGeV <= ecalEnergyMinGeV ||
      selectedLayer1Bar < 0 || selectedLayer1Bar >= 84 ||
      !ValidCuts(memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax)) {
    std::cerr << "Invalid histogram bounds, Layer-1 bar, or cut settings.\n";
    return;
  }
  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull() || maxEvents < -1 || segmentMin < -1 || segmentMax < -1) {
    std::cerr << "Supply an input directory (or OUT_DIR), maxEvents >= -1, and valid segment limits.\n";
    return;
  }

  TChain chain("T");
  if (CDetRunDataset::AddToChain(&chain, runNumber, input.Data(), segmentMin, segmentMax) <= 0 ||
      chain.GetEntries() == 0) {
    std::cerr << "No replay events found for run " << runNumber << ".\n";
    return;
  }
  const Long64_t totalEntries = chain.GetEntries();

  TTreeReader reader(&chain);

  // Pulse identity: array position identifies a pulse within this event.
  // These replay arrays use Double_t even for IDs, indices, and flags.
  TTreeReaderArray<Double_t> pulsePixelID(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> pulseIndexInChannel(reader, "earm.cdet.pulse.index");
  TTreeReaderArray<Double_t> pulseLEIndex(reader, "earm.cdet.pulse.le_index");
  TTreeReaderArray<Double_t> pulseTEIndex(reader, "earm.cdet.pulse.te_index");

  // All complete pulses, including those that were not accepted into a pair.
  TTreeReaderArray<Double_t> pulseLE(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> pulseTE(reader, "earm.cdet.pulse.tdc_te_corr");
  TTreeReaderArray<Double_t> pulseToT(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<Double_t> pulseECalDT(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> pulseCalibrationValid(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> pulseBroadQualityPass(reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> pulseECalEligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> pulseSpatialPass(reader, "earm.cdet.pulse.spatial_pass");

  // Positions and projection residuals are in meters.
  TTreeReaderArray<Double_t> pulseX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> pulseY(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderArray<Double_t> pulseECalXProjected(reader, "earm.cdet.pulse.ecal_x_proj");
  TTreeReaderArray<Double_t> pulseECalYProjected(reader, "earm.cdet.pulse.ecal_y_proj");
  TTreeReaderArray<Double_t> pulseECalXResidual(reader, "earm.cdet.pulse.ecal_x_residual");
  TTreeReaderArray<Double_t> pulseECalYResidual(reader, "earm.cdet.pulse.ecal_y_residual");

  // Accepted pairs: all pair arrays share the same pair index.
  // Source indices below refer to pulse ARRAY POSITIONS, not pixel IDs or
  // pulse.index (which is the accepted-pulse index within one channel).
  TTreeReaderArray<Double_t> pairPulseIndexL1(reader, "earm.cdet.pair.pulse_index_l1");
  TTreeReaderArray<Double_t> pairPulseIndexL2(reader, "earm.cdet.pair.pulse_index_l2");
  TTreeReaderArray<Double_t> pairPixelL1(reader, "earm.cdet.pair.pmtnum_l1");
  TTreeReaderArray<Double_t> pairPixelL2(reader, "earm.cdet.pair.pmtnum_l2");
  TTreeReaderArray<Double_t> pairLEL1(reader, "earm.cdet.pair.time_l1");
  TTreeReaderArray<Double_t> pairLEL2(reader, "earm.cdet.pair.time_l2");
  TTreeReaderArray<Double_t> pairMeanLE(reader, "earm.cdet.pair.time_mean");
  TTreeReaderArray<Double_t> pairLayerDT(reader, "earm.cdet.pair.dt"); // L2 - L1
  TTreeReaderArray<Double_t> pairECalDT(reader, "earm.cdet.pair.ecal_residual"); // ECal - mean LE
  TTreeReaderArray<Double_t> pairDX(reader, "earm.cdet.pair.dx");
  TTreeReaderArray<Double_t> pairDY(reader, "earm.cdet.pair.dy");
  TTreeReaderArray<Double_t> pairTrajectoryResidual(reader, "earm.cdet.pair.trajectory_residual");
  TTreeReaderArray<Double_t> pairECalScore(reader, "earm.cdet.pair.ecal_score");
  TTreeReaderArray<Double_t> pairCDetScore(reader, "earm.cdet.pair.score");
  TTreeReaderArray<Double_t> pairYTopology(reader, "earm.cdet.pair.y_topology");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");


  // 1. Define histograms before the event loop. Display bounds are not cuts.
  // Unique names permit comparisons across repeated interactive invocations.
  static unsigned invocation = 0;
  const TString tag = TString::Format("run%d_bar%d_%u", runNumber, selectedLayer1Bar, ++invocation);
  const TString title = TString::Format("Run %d, L1 bar %d", runNumber, selectedLayer1Bar);
  TH1D hPairMeanLE("hPairMeanLE_"+tag, title+";Corrected pair-mean LE (ns);Pairs", nLE, leMinNs, leMaxNs);
  TH1D hLayerDT("hLayerDT_"+tag, title+";t_{L2} - t_{L1} (ns);Pairs", nDT, dtMinNs, dtMaxNs);
  TH1D hECalPairDT("hECalPairDT_"+tag, title+";t_{ECal} - (t_{L1}+t_{L2})/2 (ns);Pairs", nECalDT, ecalDTMinNs, ecalDTMaxNs);
  TH1D hLEL1("hLEL1_"+tag, title+";Corrected member LE (ns);Pairs", nLE, leMinNs, leMaxNs);
  TH1D hLEL2("hLEL2_"+tag, title+";Corrected member LE (ns);Pairs", nLE, leMinNs, leMaxNs);
  constexpr double kECalZFromTargetM = 6.144;
  const int nGeometryBins = 160;
  TH2D hPairXDiffVsX1("hPairXDiffVsX1_"+tag,
      title+";x_{1} - x_{2} (m);x_{1} (m)", nGeometryBins, xDiffMinM,
      xDiffMaxM, nGeometryBins, x1MinM, x1MaxM);
  TH1D hPairXResidual("hPairXResidual_"+tag,
      title+";<x>_{CDet,pair} - x_{ECal projected} (m);Pairs",
      nGeometryBins, xResidualMinM, xResidualMaxM);
  TH1D hPairX1MinusXw("hPairX1MinusXw_"+tag,
      title+";x_{1} - x_{w} (m);Pairs", nGeometryBins,
      xResidualMinM, xResidualMaxM);
  TH1D hOutOfPlaneAngle("hOutOfPlaneAngle_"+tag,
      title+";Out-of-plane angle (degrees);Pairs", nAngleBins,
      angleMinDeg, angleMaxDeg);
  TH2D hBestPairXDiffVsX1("hBestPairXDiffVsX1_"+tag,
      title+";x_{1} - x_{2} (m);x_{1} (m)", nGeometryBins, xDiffMinM,
      xDiffMaxM, nGeometryBins, x1MinM, x1MaxM);
  TH1D hBestPairXResidual("hBestPairXResidual_"+tag,
      title+";<x>_{CDet,best pair} - x_{ECal projected} (m);Best pairs",
      nGeometryBins, xResidualMinM, xResidualMaxM);
  TH1D hBestPairX1MinusXw("hBestPairX1MinusXw_"+tag,
      title+";x_{1} - x_{w} (m);Best pairs", nGeometryBins,
      xResidualMinM, xResidualMaxM);
  TH1D hBestOutOfPlaneAngle("hBestOutOfPlaneAngle_"+tag,
      title+";Out-of-plane angle (degrees);Best pairs", nAngleBins,
      angleMinDeg, angleMaxDeg);
  TH1D hBestPairMeanLE("hBestPairMeanLE_"+tag,
      title+";Corrected best-pair mean LE (ns);Best pair per event", nLE, leMinNs, leMaxNs);
  TH1D hBestLayerDT("hBestLayerDT_"+tag,
      title+";t_{L2} - t_{L1} (ns);Best pair per event", nDT, dtMinNs, dtMaxNs);
  TH1D hBestECalPairDT("hBestECalPairDT_"+tag,
      title+";t_{ECal} - (t_{L1}+t_{L2})/2 (ns);Best pair per event", nECalDT, ecalDTMinNs, ecalDTMaxNs);
  TH1D hBestLEL1("hBestLEL1_"+tag,
      title+";Corrected best-pair member LE (ns);Best pairs", nLE, leMinNs, leMaxNs);
  TH1D hBestLEL2("hBestLEL2_"+tag,
      title+";Corrected best-pair member LE (ns);Best pairs", nLE, leMinNs, leMaxNs);
  const TString allTitle = TString::Format("Run %d, all bars", runNumber);
  TH2D hAllPairXDiffVsX1("hAllPairXDiffVsX1_"+tag,
      allTitle+";x_{1} - x_{2} (m);x_{1} (m)", nGeometryBins, xDiffMinM,
      xDiffMaxM, nGeometryBins, x1MinM, x1MaxM);
  TH1D hAllPairXResidual("hAllPairXResidual_"+tag,
      allTitle+";<x>_{CDet,pair} - x_{ECal projected} (m);Pairs",
      nGeometryBins, xResidualMinM, xResidualMaxM);
  TH1D hAllPairX1MinusXw("hAllPairX1MinusXw_"+tag,
      allTitle+";x_{1} - x_{w} (m);Pairs", nGeometryBins,
      xResidualMinM, xResidualMaxM);
  TH1D hAllOutOfPlaneAngle("hAllOutOfPlaneAngle_"+tag,
      allTitle+";Out-of-plane angle (degrees);Pairs", nAngleBins,
      angleMinDeg, angleMaxDeg);
  TH2D hBestAllPairXDiffVsX1("hBestAllPairXDiffVsX1_"+tag,
      allTitle+";x_{1} - x_{2} (m);x_{1} (m)", nGeometryBins, xDiffMinM,
      xDiffMaxM, nGeometryBins, x1MinM, x1MaxM);
  TH1D hBestAllPairXResidual("hBestAllPairXResidual_"+tag,
      allTitle+";<x>_{CDet,best pair} - x_{ECal projected} (m);Best pairs",
      nGeometryBins, xResidualMinM, xResidualMaxM);
  TH1D hBestAllPairX1MinusXw("hBestAllPairX1MinusXw_"+tag,
      allTitle+";x_{1} - x_{w} (m);Best pairs", nGeometryBins,
      xResidualMinM, xResidualMaxM);
  TH1D hBestAllOutOfPlaneAngle("hBestAllOutOfPlaneAngle_"+tag,
      allTitle+";Out-of-plane angle (degrees);Best pairs", nAngleBins,
      angleMinDeg, angleMaxDeg);
  TH2D hPairEllipseBefore("hPairEllipseBefore_"+tag,
      allTitle+";r_{x} (m);#Delta t_{pair} = t_{ECal} - <t_{CDet}>_{pair} (ns)",
      160, -0.16, 0.16, nECalDT, ecalDTMinNs, ecalDTMaxNs);
  TH2D hPairEllipseBest("hPairEllipseBest_"+tag,
      allTitle+";r_{x} (m);#Delta t_{pair} = t_{ECal} - <t_{CDet}>_{pair} (ns)",
      160, -0.16, 0.16, nECalDT, ecalDTMinNs, ecalDTMaxNs);
  hPairXDiffVsX1.SetDirectory(nullptr);
  hPairXDiffVsX1.SetStats(false);
  hBestPairXDiffVsX1.SetDirectory(nullptr);
  hBestPairXDiffVsX1.SetStats(false);
  hAllPairXDiffVsX1.SetDirectory(nullptr);
  hAllPairXDiffVsX1.SetStats(false);
  hPairEllipseBefore.SetDirectory(nullptr);
  hPairEllipseBefore.SetStats(false);
  hPairEllipseBest.SetDirectory(nullptr);
  hPairEllipseBest.SetStats(false);
  for (auto *h : {&hPairMeanLE, &hLayerDT, &hECalPairDT, &hLEL1, &hLEL2,
                  &hBestPairMeanLE, &hBestLayerDT, &hBestECalPairDT,
                  &hBestLEL1, &hBestLEL2,
                  &hPairXResidual, &hPairX1MinusXw, &hOutOfPlaneAngle,
                  &hBestPairXResidual, &hBestPairX1MinusXw,
                  &hBestOutOfPlaneAngle,
                  &hAllPairXResidual, &hAllPairX1MinusXw, &hAllOutOfPlaneAngle,
                  &hBestAllPairXResidual, &hBestAllPairX1MinusXw,
                  &hBestAllOutOfPlaneAngle})
    Prepare(*h);
  hBestAllPairXDiffVsX1.SetDirectory(nullptr);
  hBestAllPairXDiffVsX1.SetStats(false);

  auto fillGeometry = [&](TH2D &xDiffVsX1, TH1D &xResidual,
                          TH1D &x1MinusXw, TH1D &outOfPlaneAngle,
                          size_t i1, size_t i2) {
    if (!std::isfinite(pulseX[i1]) || !std::isfinite(pulseX[i2]) ||
        !std::isfinite(pulseZ[i1]) || !std::isfinite(pulseZ[i2]) ||
        !std::isfinite(*ecalX))
      return;
    const double pairMeanZ = 0.5 * (pulseZ[i1] + pulseZ[i2]);
    const double pairMeanX = 0.5 * (pulseX[i1] + pulseX[i2]);
    const double projectedECalX = *ecalX * pairMeanZ / kECalZFromTargetM;
    const double deltaZ = pulseZ[i1] - pulseZ[i2];
    const double xw = pulseX[i2] + (*ecalX / kECalZFromTargetM) * deltaZ;
    xDiffVsX1.Fill(pulseX[i1] - pulseX[i2], pulseX[i1]);
    xResidual.Fill(pairMeanX - projectedECalX);
    x1MinusXw.Fill(pulseX[i1] - xw);
    const double z1 = pulseZ[i1], z2 = pulseZ[i2];
    const double numerator = z1 * pulseX[i1] + z2 * pulseX[i2] +
                             kECalZFromTargetM * (*ecalX);
    const double denominator = z1*z1 + z2*z2 +
                               kECalZFromTargetM*kECalZFromTargetM;
    if (denominator > 0.0)
      outOfPlaneAngle.Fill(std::atan(numerator / denominator) * 180.0 / 3.14159265358979323846);
  };

  // All 168 paired-member histograms are filled before the selected-bar cut.
  std::array<std::unique_ptr<TH1D>, 168> histograms;
  std::array<std::unique_ptr<TH1D>, 168> bestHistograms;
  for (int bar = 0; bar < 168; ++bar) {
    histograms[bar].reset(new TH1D(TString::Format("hPairMemberDT_%s_bar%d", tag.Data(), bar), TString::Format("Run %d, global bar %d;t_{ECal} - t_{CDet,member} (ns);Paired members", runNumber, bar), nECalDT, ecalDTMinNs, ecalDTMaxNs));
    Prepare(*histograms[bar]);
    bestHistograms[bar].reset(new TH1D(TString::Format("hBestPairMemberDT_%s_bar%d", tag.Data(), bar), TString::Format("Run %d, global bar %d; t_{ECal} - t_{CDet,member} (ns);Best pairs", runNumber, bar), nECalDT, ecalDTMinNs, ecalDTMaxNs));
    Prepare(*bestHistograms[bar]);
  }

  Long64_t eventsRead = 0, goodPairEventCount = 0, malformedPairs = 0;
  const TString cuts = CutLabel(memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax);
  std::cout << cuts << "\nECal adctime event cut: [" << ecalTimeMinNs
            << ", " << ecalTimeMaxNs << "] ns; ECal energy cut: ["
            << ecalEnergyMinGeV << ", " << ecalEnergyMaxGeV << "] GeV.\n";
  while ((maxEvents < 0 || eventsRead < maxEvents) && reader.Next()) {
    ++eventsRead;
    bool eventHasGoodPair = false;
    if (!std::isfinite(*ecalTime) || *ecalTime < ecalTimeMinNs ||
        *ecalTime > ecalTimeMaxNs || !std::isfinite(*ecalEnergy) ||
        *ecalEnergy < ecalEnergyMinGeV || *ecalEnergy > ecalEnergyMaxGeV)
      continue;
    const size_t nPairs = pairPixelL1.GetSize();
    if (pulsePixelID.GetSize() != pulseToT.GetSize() ||
        pulsePixelID.GetSize() != pulseECalDT.GetSize() ||
        pulsePixelID.GetSize() != pulseX.GetSize() ||
        pulsePixelID.GetSize() != pulseY.GetSize() ||
        pulsePixelID.GetSize() != pulseZ.GetSize() ||
        pairPixelL2.GetSize() != nPairs || pairPulseIndexL1.GetSize() != nPairs || pairPulseIndexL2.GetSize() != nPairs ||
        pairLEL1.GetSize() != nPairs || pairLEL2.GetSize() != nPairs ||
        pairMeanLE.GetSize() != nPairs || pairLayerDT.GetSize() != nPairs ||
        pairECalDT.GetSize() != nPairs || pairTrajectoryResidual.GetSize() != nPairs ||
        pairECalScore.GetSize() != nPairs) {
      std::cerr << "Mismatched pair arrays at entry " << reader.GetCurrentEntry() << '\n';
      return;
    }
    struct PairCandidate {
      size_t pair = 0, i1 = 0, i2 = 0, id1 = 0, id2 = 0;
      double rank = std::numeric_limits<double>::infinity();
    };
    std::vector<PairCandidate> preEllipseCandidates;
    std::vector<PairCandidate> acceptedCandidates;
    for (size_t pair = 0; pair < nPairs; ++pair) {
      // 2. Validate member references and apply common cuts to all bars.
      size_t i1, i2, id1, id2;
      if (!Index(pairPulseIndexL1[pair], pulseToT.GetSize(), i1) ||
          !Index(pairPulseIndexL2[pair], pulseToT.GetSize(), i2) ||
          !Index(pulsePixelID[i1], 2688, id1) || !Index(pulsePixelID[i2], 2688, id2) ||
          id1 >= 1344 || id2 < 1344 || pairPixelL1[pair] != id1 || pairPixelL2[pair] != id2 ||
          !std::isfinite(pulseECalDT[i1]) || !std::isfinite(pulseECalDT[i2])) {
        ++malformedPairs;
        continue;
      }
      // All five histograms use exactly the same finite selected pairs.
      if (!std::isfinite(pairLEL1[pair]) || !std::isfinite(pairLEL2[pair]) ||
          !std::isfinite(pairMeanLE[pair]) || !std::isfinite(pairECalDT[pair]) ||
          !std::isfinite(pairTrajectoryResidual[pair])) {
        ++malformedPairs;
        continue;
      }
      // Fill the candidate population before applying the ECal trajectory-time
      // ellipse. ToT and optional inter-layer timing cuts are still applied.
      const bool passesPreEllipse = PassCuts(
          pulseToT[i1], pulseToT[i2], pairLayerDT[pair],
          pairTrajectoryResidual[pair], pairECalDT[pair],
          memberToTMinNs, memberToTMaxNs, layerDTMaxNs, -1);
      if (passesPreEllipse) {
        hPairEllipseBefore.Fill(pairTrajectoryResidual[pair], pairECalDT[pair]);
        const double radiusSquared = PairEllipseRadiusSquared(
            pairTrajectoryResidual[pair], pairECalDT[pair]);
        const double rank = std::isfinite(pairECalScore[pair])
                                ? pairECalScore[pair] : radiusSquared;
        preEllipseCandidates.push_back({pair, i1, i2, id1, id2, rank});
      }
      if (!PassCuts(pulseToT[i1], pulseToT[i2], pairLayerDT[pair],
                    pairTrajectoryResidual[pair], pairECalDT[pair],
                    memberToTMinNs, memberToTMaxNs, layerDTMaxNs,
                    pairRadiusMax))
        continue;
      eventHasGoodPair = true;
      const double radiusSquared = PairEllipseRadiusSquared(
          pairTrajectoryResidual[pair], pairECalDT[pair]);
      const double rank = std::isfinite(pairECalScore[pair])
                              ? pairECalScore[pair] : radiusSquared;
      acceptedCandidates.push_back({pair, i1, i2, id1, id2, rank});
      // 3. Detector-wide paired-member timing, under the same study cuts.
      // This is ECal minus each member's LE, NOT ECal minus pair-mean time.
      histograms[id1/16]->Fill(pulseECalDT[i1]);
      histograms[id2/16]->Fill(pulseECalDT[i2]);
      // Geometry diagnostics are filled for every accepted pair before the
      // selected-bar restriction, so the all-detector canvas is independent
      // of the focused bar.
      fillGeometry(hAllPairXDiffVsX1, hAllPairXResidual, hAllPairX1MinusXw,
                   hAllOutOfPlaneAngle, i1, i2);

      // Only the first two canvases restrict the Layer-1 bar.
      // Keep the actual matched L2 partner regardless of its bar number.
      if (id1/16 != static_cast<size_t>(selectedLayer1Bar))
        continue;
      hPairMeanLE.Fill(pairMeanLE[pair]);
      hLayerDT.Fill(pairLayerDT[pair]);
      hECalPairDT.Fill(pairECalDT[pair]);
      hLEL1.Fill(pairLEL1[pair]);
      hLEL2.Fill(pairLEL2[pair]);
      fillGeometry(hPairXDiffVsX1, hPairXResidual, hPairX1MinusXw,
                   hOutOfPlaneAngle, i1, i2);
    }
    auto bestCandidate = [](const std::vector<PairCandidate>& candidates,
                            bool selectedBarOnly, int selectedBar)
        -> const PairCandidate* {
      const PairCandidate *best = nullptr;
      for (const PairCandidate& candidate : candidates) {
        if (selectedBarOnly && candidate.id1 / 16 !=
                                  static_cast<size_t>(selectedBar))
          continue;
        if (!best || candidate.rank < best->rank ||
            (candidate.rank == best->rank && candidate.pair < best->pair))
          best = &candidate;
      }
      return best;
    };
    const PairCandidate *bestPreEllipse =
        bestCandidate(preEllipseCandidates, false, selectedLayer1Bar);
    if (bestPreEllipse)
      hPairEllipseBest.Fill(pairTrajectoryResidual[bestPreEllipse->pair],
                            pairECalDT[bestPreEllipse->pair]);

    const PairCandidate *bestGlobal =
        bestCandidate(acceptedCandidates, false, selectedLayer1Bar);
    if (bestGlobal) {
      const size_t pair = bestGlobal->pair;
      bestHistograms[bestGlobal->id1 / 16]->Fill(
          pulseECalDT[bestGlobal->i1]);
      bestHistograms[bestGlobal->id2 / 16]->Fill(
          pulseECalDT[bestGlobal->i2]);
      fillGeometry(hBestAllPairXDiffVsX1, hBestAllPairXResidual,
                   hBestAllPairX1MinusXw, hBestAllOutOfPlaneAngle,
                   bestGlobal->i1, bestGlobal->i2);
    }
    const PairCandidate *bestSelected =
        bestCandidate(acceptedCandidates, true, selectedLayer1Bar);
    if (bestSelected) {
      const size_t pair = bestSelected->pair;
      hBestPairMeanLE.Fill(pairMeanLE[pair]);
      hBestLayerDT.Fill(pairLayerDT[pair]);
      hBestECalPairDT.Fill(pairECalDT[pair]);
      hBestLEL1.Fill(pairLEL1[pair]);
      hBestLEL2.Fill(pairLEL2[pair]);
      fillGeometry(hBestPairXDiffVsX1, hBestPairXResidual,
                   hBestPairX1MinusXw, hBestOutOfPlaneAngle,
                   bestSelected->i1, bestSelected->i2);
    }
    if (eventHasGoodPair)
      ++goodPairEventCount;
    if (eventsRead % 1000 == 0)
      std::cout << "[pairing] event " << eventsRead << " of " << totalEntries
                << "; events with a good pair: " << goodPairEventCount << '\n';
  }
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd &&
      (maxEvents < 0 || eventsRead < maxEvents)) {
    std::cerr << "Tree reading stopped early; check branch availability/types.\n";
    return;
  }
  std::cout << "Read " << eventsRead << " of " << totalEntries
            << " events; events with a good pair: " << goodPairEventCount
            << "; skipped malformed pairs: " << malformedPairs << '\n';

  // Saving is opt-in; savePlots=false creates no directory or output files.
  TString outputPrefix;
  if (savePlots) {
    TString directory(outputDirectory);
    gSystem->ExpandPathName(directory);
    if (gSystem->AccessPathName(directory) && gSystem->mkdir(directory, true) != 0) {
      std::cerr << "Cannot create output directory: " << directory << '\n';
      return;
    }
    if (gSystem->AccessPathName(directory, kWritePermission)) {
      std::cerr << "Output directory is not writable: " << directory << '\n';
      return;
    }
    outputPrefix = directory + TString::Format("/CDet_run%d_bar%03d", runNumber, selectedLayer1Bar);
  }

  // 4. Sigma here is full-sample standard deviation, NOT a Gaussian-core fit.
  for (auto *h : {&hPairMeanLE, &hLayerDT, &hECalPairDT, &hLEL1, &hLEL2,
                  &hBestPairMeanLE, &hBestLayerDT, &hBestECalPairDT,
                  &hBestLEL1, &hBestLEL2})
    Report(*h);
  if (hLayerDT.GetEntries() >= 2)
    std::cout << "Equal-independent-layer estimate SD(dt)/sqrt(2) = "
              << hLayerDT.GetStdDev()/std::sqrt(2.0) << " ns (assumption-dependent).\n";
  std::cout << "Widths are conditioned by stored-pair and study cuts; ECal residuals also include reference timing.\n";

  // 5. Draw copies so histograms remain on interactive canvases after return.
  auto *timingCanvas = new TCanvas("cPairTiming_"+tag, title+" | "+cuts, 1500, 900);
  timingCanvas->Divide(3, 2);
  int pad = 0;
  for (auto *h : {&hPairMeanLE, &hLayerDT, &hECalPairDT}) {
    timingCanvas->cd(++pad);
    h->DrawCopy("HIST");
    Annotate(*h);
  }
  for (auto *h : {&hBestPairMeanLE, &hBestLayerDT, &hBestECalPairDT}) {
    timingCanvas->cd(++pad);
    h->DrawCopy("HIST");
    Annotate(*h);
  }
  timingCanvas->Update();

  // 6. Matched layer spectra, with separate standard deviations.
  auto *layerCanvas = new TCanvas("cPairLayers_"+tag, title+" | "+cuts, 900, 900);
  layerCanvas->Divide(1, 2);
  hLEL1.SetLineColor(kBlue+1);
  hLEL2.SetLineColor(kRed+1);
  layerCanvas->cd(1);
  hLEL1.SetMaximum(1.3*std::max(1.0, std::max(hLEL1.GetMaximum(), hLEL2.GetMaximum())));
  auto *drawL1 = hLEL1.DrawCopy("HIST");
  auto *drawL2 = hLEL2.DrawCopy("HIST SAME");
  auto *legend = new TLegend(0.40, 0.73, 0.89, 0.89);
  legend->SetBorderSize(0);
  legend->AddEntry(drawL1, "L1: "+WidthLabel(hLEL1), "l");
  legend->AddEntry(drawL2, "L2 partner: "+WidthLabel(hLEL2), "l");
  legend->Draw();
  layerCanvas->cd(2);
  hBestLEL1.SetLineColor(kBlue+1);
  hBestLEL2.SetLineColor(kRed+1);
  hBestLEL1.SetMaximum(1.3*std::max(1.0, std::max(hBestLEL1.GetMaximum(), hBestLEL2.GetMaximum())));
  auto *drawBestL1 = hBestLEL1.DrawCopy("HIST");
  auto *drawBestL2 = hBestLEL2.DrawCopy("HIST SAME");
  auto *bestLegend = new TLegend(0.40, 0.73, 0.89, 0.89);
  bestLegend->SetBorderSize(0);
  bestLegend->AddEntry(drawBestL1, "Best L1: "+WidthLabel(hBestLEL1), "l");
  bestLegend->AddEntry(drawBestL2, "Best L2 partner: "+WidthLabel(hBestLEL2), "l");
  bestLegend->Draw();
  layerCanvas->Update();
  if (savePlots) {
    timingCanvas->SaveAs(outputPrefix+"_timing.pdf");
    timingCanvas->SaveAs(outputPrefix+"_timing.png");
    layerCanvas->SaveAs(outputPrefix+"_layers.pdf");
    layerCanvas->SaveAs(outputPrefix+"_layers.png");
  }

  // 7. Pair geometry and rough transport-coordinate angle diagnostics.
  Report(hPairXResidual, "m");
  Report(hOutOfPlaneAngle, "degrees");
  Report(hBestPairXResidual, "m");
  Report(hBestOutOfPlaneAngle, "degrees");
  if (hPairXResidual.GetEntries() >= 2)
    std::cout << "Rough CDet x position resolution estimate = "
              << hPairXResidual.GetStdDev() * 1000.0 << " mm (pair residual SD).\n";
  if (hBestPairXResidual.GetEntries() >= 2)
    std::cout << "Rough best-pair CDet x position estimate = "
              << hBestPairXResidual.GetStdDev() * 1000.0 << " mm (best-pair residual SD).\n";
  auto *geometryCanvas = new TCanvas("cPairGeometry_"+tag,
      title+" | pair geometry", 1500, 900);
  geometryCanvas->Divide(3, 2);
  geometryCanvas->cd(1);
  hPairXDiffVsX1.DrawCopy("COLZ");
  geometryCanvas->cd(2);
  hPairXResidual.DrawCopy("HIST");
  Annotate(hPairXResidual, "m");
  geometryCanvas->cd(3);
  hOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hOutOfPlaneAngle, "degrees");
  geometryCanvas->cd(4);
  hBestPairXDiffVsX1.DrawCopy("COLZ");
  geometryCanvas->cd(5);
  hBestPairXResidual.DrawCopy("HIST");
  Annotate(hBestPairXResidual, "m");
  geometryCanvas->cd(6);
  hBestOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hBestOutOfPlaneAngle, "degrees");
  geometryCanvas->Update();
  if (savePlots) {
    geometryCanvas->SaveAs(outputPrefix+"_geometry.pdf");
    geometryCanvas->SaveAs(outputPrefix+"_geometry.png");
  }
  auto *allGeometryCanvas = new TCanvas("cAllPairGeometry_"+tag,
      allTitle+" | pair geometry", 1500, 900);
  allGeometryCanvas->Divide(3, 2);
  allGeometryCanvas->cd(1);
  hAllPairXDiffVsX1.DrawCopy("COLZ");
  allGeometryCanvas->cd(2);
  hAllPairXResidual.DrawCopy("HIST");
  Annotate(hAllPairXResidual, "m");
  allGeometryCanvas->cd(3);
  hAllOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hAllOutOfPlaneAngle, "degrees");
  allGeometryCanvas->cd(4);
  hBestAllPairXDiffVsX1.DrawCopy("COLZ");
  allGeometryCanvas->cd(5);
  hBestAllPairXResidual.DrawCopy("HIST");
  Annotate(hBestAllPairXResidual, "m");
  allGeometryCanvas->cd(6);
  hBestAllOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hBestAllOutOfPlaneAngle, "degrees");
  allGeometryCanvas->Update();
  if (savePlots) {
    allGeometryCanvas->SaveAs(outputPrefix+"_geometry_all.pdf");
    allGeometryCanvas->SaveAs(outputPrefix+"_geometry_all.png");
  }

  // Alternative geometry canvases: use the ECal-guided Layer-2-to-Layer-1
  // projection x_w, which is the trajectory residual up to a sign.
  Report(hPairX1MinusXw, "m");
  Report(hBestPairX1MinusXw, "m");
  if (hPairX1MinusXw.GetEntries() >= 2)
    std::cout << "Rough Layer-1 x resolution estimate from x1-xw = "
              << hPairX1MinusXw.GetStdDev() * 1000.0 << " mm.\n";
  if (hBestPairX1MinusXw.GetEntries() >= 2)
    std::cout << "Rough best-pair Layer-1 x estimate from x1-xw = "
              << hBestPairX1MinusXw.GetStdDev() * 1000.0 << " mm.\n";
  auto *xwGeometryCanvas = new TCanvas("cPairGeometryXw_"+tag,
      title+" | x1-xw geometry", 1500, 900);
  xwGeometryCanvas->Divide(3, 2);
  xwGeometryCanvas->cd(1);
  hPairXDiffVsX1.DrawCopy("COLZ");
  xwGeometryCanvas->cd(2);
  hPairX1MinusXw.DrawCopy("HIST");
  Annotate(hPairX1MinusXw, "m");
  xwGeometryCanvas->cd(3);
  hOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hOutOfPlaneAngle, "degrees");
  xwGeometryCanvas->cd(4);
  hBestPairXDiffVsX1.DrawCopy("COLZ");
  xwGeometryCanvas->cd(5);
  hBestPairX1MinusXw.DrawCopy("HIST");
  Annotate(hBestPairX1MinusXw, "m");
  xwGeometryCanvas->cd(6);
  hBestOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hBestOutOfPlaneAngle, "degrees");
  xwGeometryCanvas->Update();
  if (savePlots) {
    xwGeometryCanvas->SaveAs(outputPrefix+"_geometry_xw.pdf");
    xwGeometryCanvas->SaveAs(outputPrefix+"_geometry_xw.png");
  }
  auto *allXwGeometryCanvas = new TCanvas("cAllPairGeometryXw_"+tag,
      allTitle+" | x1-xw geometry", 1500, 900);
  allXwGeometryCanvas->Divide(3, 2);
  allXwGeometryCanvas->cd(1);
  hAllPairXDiffVsX1.DrawCopy("COLZ");
  allXwGeometryCanvas->cd(2);
  hAllPairX1MinusXw.DrawCopy("HIST");
  Annotate(hAllPairX1MinusXw, "m");
  allXwGeometryCanvas->cd(3);
  hAllOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hAllOutOfPlaneAngle, "degrees");
  allXwGeometryCanvas->cd(4);
  hBestAllPairXDiffVsX1.DrawCopy("COLZ");
  allXwGeometryCanvas->cd(5);
  hBestAllPairX1MinusXw.DrawCopy("HIST");
  Annotate(hBestAllPairX1MinusXw, "m");
  allXwGeometryCanvas->cd(6);
  hBestAllOutOfPlaneAngle.DrawCopy("HIST");
  Annotate(hBestAllOutOfPlaneAngle, "degrees");
  allXwGeometryCanvas->Update();
  if (savePlots) {
    allXwGeometryCanvas->SaveAs(outputPrefix+"_geometry_xw_all.pdf");
    allXwGeometryCanvas->SaveAs(outputPrefix+"_geometry_xw_all.png");
  }

  // 11. Candidate pair population before the ECal trajectory-time ellipse.
  auto *ellipseCanvas = new TCanvas("cPairEllipse_"+tag,
      title+" | pre-ellipse pair candidates", 900, 900);
  ellipseCanvas->Divide(1, 2);
  ellipseCanvas->cd(1);
  ellipseCanvas->SetRightMargin(0.14);
  hPairEllipseBefore.DrawCopy("COLZ");
  if (pairRadiusMax > 0) {
    TEllipse ellipse(0.0, CDetPairingPlots::kPairEllipseTimingCenterNs,
        pairRadiusMax * CDetPairingPlots::kPairEllipseTrajectoryScaleM,
        pairRadiusMax * CDetPairingPlots::kPairEllipseTimingScaleNs);
    ellipse.SetFillStyle(0);
    ellipse.SetLineColor(kRed+1);
    ellipse.SetLineWidth(3);
    ellipse.Draw("SAME");
  }
  ellipseCanvas->cd(2);
  hPairEllipseBest.DrawCopy("COLZ");
  if (pairRadiusMax > 0) {
    TEllipse ellipse(0.0, CDetPairingPlots::kPairEllipseTimingCenterNs,
        pairRadiusMax * CDetPairingPlots::kPairEllipseTrajectoryScaleM,
        pairRadiusMax * CDetPairingPlots::kPairEllipseTimingScaleNs);
    ellipse.SetFillStyle(0);
    ellipse.SetLineColor(kRed+1);
    ellipse.SetLineWidth(3);
    ellipse.Draw("SAME");
  }
  ellipseCanvas->Update();
  if (savePlots) {
    ellipseCanvas->SaveAs(outputPrefix+"_ellipse.pdf");
    ellipseCanvas->SaveAs(outputPrefix+"_ellipse.png");
  }

  // 12. Third canvas: sigma versus bar number within each detector layer.
  // Graphs omit low-statistics bars; the optional CSV retains every bar/status.
  auto *g1 = new TGraphErrors();
  auto *g2 = new TGraphErrors();
  auto *gBest1 = new TGraphErrors();
  auto *gBest2 = new TGraphErrors();
  g1->SetName("gPairSigmaL1_"+tag);
  g2->SetName("gPairSigmaL2_"+tag);
  gBest1->SetName("gBestPairSigmaL1_"+tag);
  gBest2->SetName("gBestPairSigmaL2_"+tag);
  std::ofstream csv;
  if (savePlots) {
    csv.open((TString(outputPrefix)+"_sigma_vs_bar.csv").Data());
    if (!csv) std::cerr << "Could not open CSV output.\n";
    else {
      csv.precision(12);
      csv << "global_bar,layer,local_bar,entries,mean_ns,stddev_ns,stddev_error_ns,underflow,overflow,status\n";
    }
  }
  for (int bar = 0; bar < 168; ++bar) {
    const auto &h = *histograms[bar];
    const bool valid = h.GetEntries() >= minEntriesPerBar;
    if (valid) {
      auto *g = bar < 84 ? g1 : g2;
      const int point = g->GetN();
      g->SetPoint(point, bar%84, h.GetStdDev());
      g->SetPointError(point, 0, h.GetStdDevError());
    }
    const auto &bestH = *bestHistograms[bar];
    if (bestH.GetEntries() >= minEntriesPerBar) {
      auto *g = bar < 84 ? gBest1 : gBest2;
      const int point = g->GetN();
      g->SetPoint(point, bar%84, bestH.GetStdDev());
      g->SetPointError(point, 0, bestH.GetStdDevError());
    }
    // if (csv) {
    //   csv << bar << ',' << bar/84+1 << ',' << bar%84 << ',' << h.GetEntries() << ',';
    //   if (h.GetEntries() >= 2)
    //     csv << h.GetMean() << ',' << h.GetStdDev() << ',' << h.GetStdDevError();
    //   else csv << ",,";
    //   csv << ',' << h.GetBinContent(0) << ',' << h.GetBinContent(nECalDT+1)
    //       << ',' << (valid ? "accepted" : "insufficient_entries") << '\n';
    // }
  }
  auto *canvas = new TCanvas("cPairSigmaVsBar_"+tag, TString::Format("Run %d | ", runNumber)+cuts, 1200, 900);
  canvas->Divide(1, 2);
  double ymax = 1;
  for (auto *g : {g1, g2})
    for (int i = 0; i < g->GetN(); ++i)
      ymax = std::max(ymax, g->GetPointY(i)+g->GetErrorY(i));
  canvas->cd(1);
  canvas->DrawFrame(-0.5, 0, 83.5, 1.2*ymax, TString::Format("Run %d: all paired-member timing widths;Bar number within layer;SD(t_{ECal}-t_{CDet,member}) (ns)", runNumber));
  g1->SetMarkerStyle(20); g1->SetMarkerColor(kBlue+1); g1->SetLineColor(kBlue+1);
  g2->SetMarkerStyle(22); g2->SetMarkerColor(kRed+1); g2->SetLineColor(kRed+1);
  if (g1->GetN()) g1->Draw("P SAME");
  if (g2->GetN()) g2->Draw("P SAME");
  auto *barLegend = new TLegend(0.65, 0.76, 0.89, 0.89);
  barLegend->AddEntry(g1, "Layer 1", "p");
  barLegend->AddEntry(g2, "Layer 2", "p");
  barLegend->Draw();
  double bestYmax = 1;
  for (auto *g : {gBest1, gBest2})
    for (int i = 0; i < g->GetN(); ++i)
      bestYmax = std::max(bestYmax, g->GetPointY(i)+g->GetErrorY(i));
  canvas->cd(2);
  canvas->DrawFrame(-0.5, 0, 83.5, 1.2*bestYmax, TString::Format("Run %d: best-pair timing widths;Bar number within layer;SD(t_{ECal}-t_{CDet,member}) (ns)", runNumber));
  gBest1->SetMarkerStyle(20); gBest1->SetMarkerColor(kBlue+1); gBest1->SetLineColor(kBlue+1);
  gBest2->SetMarkerStyle(22); gBest2->SetMarkerColor(kRed+1); gBest2->SetLineColor(kRed+1);
  if (gBest1->GetN()) gBest1->Draw("P SAME");
  if (gBest2->GetN()) gBest2->Draw("P SAME");
  auto *bestBarLegend = new TLegend(0.65, 0.76, 0.89, 0.89);
  bestBarLegend->AddEntry(gBest1, "Layer 1", "p");
  bestBarLegend->AddEntry(gBest2, "Layer 2", "p");
  bestBarLegend->Draw();
  canvas->Update();
  if (savePlots) {
    canvas->SaveAs(outputPrefix+"_sigma_vs_bar.pdf");
    canvas->SaveAs(outputPrefix+"_sigma_vs_bar.png");
  }
  std::cout << "Widths include overflow and selected-sample tails; they are not Gaussian fit sigmas or intrinsic resolutions.\n";
}

// TEnv-style configuration, kept separate from the event loop. Unknown keys
// and malformed numbers are errors so a mistyped cut cannot silently default.
void PlotPairingVariables(const char *configFile,
                         const char *inputDirectoryOverride = nullptr) {
  if (!configFile || !configFile[0]) {
    std::cerr << "Supply a pairing configuration file.\n";
    return;
  }
  TEnv config;
  if (config.ReadFile(configFile, kEnvLocal) != 0) {
    std::cerr << "Cannot read configuration: " << configFile << '\n';
    return;
  }
  const std::set<std::string> keys = {
    "config.version", "analysis.run_number", "analysis.input_directory",
    "analysis.events", "analysis.min_segment", "analysis.max_segment",
    "analysis.layer1_bar", "cuts.member_tot_min_ns",
    "cuts.member_tot_max_ns", "cuts.layer_dt_max_ns", "cuts.pair_radius_max",
    "cuts.ecal_time_min_ns", "cuts.ecal_time_max_ns",
    "cuts.ecal_energy_min_gev", "cuts.ecal_energy_max_gev",
    "plots.bin_width_ns", "plots.le_min_ns", "plots.le_max_ns",
    "plots.layer_dt_min_ns", "plots.layer_dt_max_ns", "plots.ecal_dt_min_ns",
    "plots.ecal_dt_max_ns", "plots.min_entries_per_bar", "plots.x_diff_min_m",
    "plots.x_diff_max_m", "plots.x1_min_m", "plots.x1_max_m",
    "plots.x_residual_min_m", "plots.x_residual_max_m", "plots.angle_min_deg",
    "plots.angle_max_deg", "plots.angle_bin_width_deg", "output.save_plots",
    "output.directory"
  };
  TIter next(config.GetTable());
  while (auto *record = next()) {
    if (!keys.count(record->GetName())) {
      std::cerr << "Unknown pairing configuration key: " << record->GetName() << '\n';
      return;
    }
  }
  try {
    auto number = [&](const char *key, double fallback) {
      if (!config.Defined(key)) return fallback;
      std::string text = config.GetValue(key, "");
      size_t used = 0;
      const double value = std::stod(text, &used);
      if (text.find_first_not_of(" \t\r\n", used) != std::string::npos || !std::isfinite(value))
        throw std::runtime_error(std::string("Invalid number for ")+key);
      return value;
    };
    auto integer = [&](const char *key, Long64_t fallback) {
      if (!config.Defined(key)) return fallback;
      std::string text = config.GetValue(key, "");
      size_t used = 0;
      const Long64_t value = std::stoll(text, &used);
      if (text.find_first_not_of(" \t\r\n", used) != std::string::npos)
        throw std::runtime_error(std::string("Invalid integer for ")+key);
      return value;
    };
    if (integer("config.version", -1) != 1)
      throw std::runtime_error("config.version must be 1");
    const Long64_t run = integer("analysis.run_number", 6077);
    const Long64_t bar = integer("analysis.layer1_bar", 30);
    const Long64_t segmentMin = integer("analysis.min_segment", -1);
    const Long64_t segmentMax = integer("analysis.max_segment", -1);
    const Long64_t minimum = integer("plots.min_entries_per_bar", 30);
    const Long64_t save = integer("output.save_plots", 0);
    if (run <= 0 || run > 2147483647 || bar < 0 || bar >= 84 ||
        segmentMin < -1 || segmentMax < -1 ||
        minimum < 2 || minimum > 2147483647 || (save != 0 && save != 1))
      throw std::runtime_error("Invalid run, bar, minimum entries, or save flag (use 0 or 1)");
    const TString input = inputDirectoryOverride ? inputDirectoryOverride :
        config.GetValue("analysis.input_directory", "");
    std::cout << "Pairing configuration: " << configFile << '\n';
    PlotPairingVariables(static_cast<int>(run), input.Data(), integer("analysis.events", -1),
        static_cast<int>(bar), number("plots.bin_width_ns", 0.5),
        number("plots.le_min_ns", 0), number("plots.le_max_ns", 80),
        number("plots.layer_dt_min_ns", -30), number("plots.layer_dt_max_ns", 30),
        number("plots.ecal_dt_min_ns", -100), number("plots.ecal_dt_max_ns", 50),
        number("cuts.member_tot_min_ns", 0), number("cuts.member_tot_max_ns", 1e9),
        number("cuts.layer_dt_max_ns", -1), number("cuts.pair_radius_max", -1),
        save != 0, config.GetValue("output.directory", "pairing_plots"),
        static_cast<int>(minimum), static_cast<int>(segmentMin),
        static_cast<int>(segmentMax), number("plots.x_diff_min_m", -0.5),
        number("plots.x_diff_max_m", 0.5), number("plots.x1_min_m", -1.6),
        number("plots.x1_max_m", 1.6), number("plots.x_residual_min_m", -0.2),
        number("plots.x_residual_max_m", 0.2), number("plots.angle_min_deg", -40),
        number("plots.angle_max_deg", 40), number("plots.angle_bin_width_deg", 0.1),
        number("cuts.ecal_time_min_ns", -10),
        number("cuts.ecal_time_max_ns", 4), number("cuts.ecal_energy_min_gev", 3.0),
        number("cuts.ecal_energy_max_gev", 4.5));
  } catch (const std::exception &error) {
    std::cerr << "Invalid pairing configuration: " << error.what() << '\n';
  }
}
