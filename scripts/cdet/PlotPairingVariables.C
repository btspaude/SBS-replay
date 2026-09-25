#include <TCanvas.h>
#include <TLegend.h>
#include <algorithm>
#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <iostream>

#include "CDetRunDataset.h"
#include "PairingPlotHelpers.h"

// Interactive ROOT:
// .L PlotPairingVariables.C+
// PlotPairingVariables(6077, "/path/to/Rootfiles", 10000);
// maxEvents = -1 reads the full selected dataset. Input defaults to OUT_DIR.
// Only pulse and pair branches are read. Corrected times below are already ns.
// Arguments after maxEvents: L1 bar, bin width, LE display min/max,
// layer-dt display min/max, ECal-dt display min/max, member ToT min/max,
// maximum |L2-L1|, maximum ECal ellipse radius, optional output prefix.
// -1 disables either additional pair cut. Default selects stored pairs as-is.
// Example with both members ToT 12..35 ns, |dt| <= 10 ns, and radius <= 2:
// PlotPairingVariables(6077, "/path/to/Rootfiles", -1, 30,
//     0.5, 0, 80, -30, 30, -100, 50, 12, 35, 10, 2, "bar30_tight");
// Outputs (when requested): <prefix>_timing.pdf and <prefix>_layers.pdf.
// These tighten the STORED pair population; they do not rerun assignment.
// Standard deviations and their ROOT moment-based error estimates include
// all finite selected values, including histogram underflow/overflow.
// No Gaussian fit is used. Pair/event correlations are not modeled in errors.
void PlotPairingVariables(int runNumber = 6077,
                         const char *inputDirectory = nullptr,
                         Long64_t maxEvents = -1, int selectedLayer1Bar = 30,
                         double binWidthNs = 0.5, double leMinNs = 0, double leMaxNs = 80,
                         double dtMinNs = -30, double dtMaxNs = 30,
                         double ecalDTMinNs = -100, double ecalDTMaxNs = 50,
                         double memberToTMinNs = 0, double memberToTMaxNs = 1e9,
                         double layerDTMaxNs = -1, double pairRadiusMax = -1,
                         const char *outputPrefix = "") {
  using namespace CDetPairingPlots;
  const int nLE = Bins(binWidthNs, leMinNs, leMaxNs);
  const int nDT = Bins(binWidthNs, dtMinNs, dtMaxNs);
  const int nECalDT = Bins(binWidthNs, ecalDTMinNs, ecalDTMaxNs);
  if (!nLE || !nDT || !nECalDT || selectedLayer1Bar < 0 || selectedLayer1Bar >= 84 ||
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
  if (input.IsNull() || maxEvents < -1) {
    std::cerr << "Supply an input directory (or OUT_DIR) and maxEvents >= -1.\n";
    return;
  }

  TChain chain("T");
  if (CDetRunDataset::AddToChain(&chain, runNumber, input.Data()) <= 0 ||
      chain.GetEntries() == 0) {
    std::cerr << "No replay events found for run " << runNumber << ".\n";
    return;
  }

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
  TTreeReaderArray<Double_t> pairCDetScore(reader, "earm.cdet.pair.score");
  TTreeReaderArray<Double_t> pairECalScore(reader, "earm.cdet.pair.ecal_score"); // radius squared
  TTreeReaderArray<Double_t> pairYTopology(reader, "earm.cdet.pair.y_topology");


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
  for (auto *h : {&hPairMeanLE, &hLayerDT, &hECalPairDT, &hLEL1, &hLEL2})
    Prepare(*h);

  Long64_t eventsRead = 0, malformedPairs = 0;
  const TString cuts = CutLabel(memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax);
  std::cout << cuts << "\nNo additional ECal event cut is imposed.\n";
  while ((maxEvents < 0 || eventsRead < maxEvents) && reader.Next()) {
    ++eventsRead;
    const size_t nPairs = pairPixelL1.GetSize();
    if (pairPulseIndexL1.GetSize() != nPairs || pairPulseIndexL2.GetSize() != nPairs ||
        pairLEL1.GetSize() != nPairs || pairLEL2.GetSize() != nPairs ||
        pairMeanLE.GetSize() != nPairs || pairLayerDT.GetSize() != nPairs ||
        pairECalDT.GetSize() != nPairs || pairECalScore.GetSize() != nPairs) {
      std::cerr << "Mismatched pair arrays at entry " << reader.GetCurrentEntry() << '\n';
      return;
    }
    for (size_t pair = 0; pair < nPairs; ++pair) {
      // 2. Select the L1 bar and keep its actual matched L2 partner.
      if (!(pairPixelL1[pair] >= 16*selectedLayer1Bar &&
            pairPixelL1[pair] < 16*(selectedLayer1Bar+1)))
        continue;
      size_t i1, i2;
      if (!Index(pairPulseIndexL1[pair], pulseToT.GetSize(), i1) ||
          !Index(pairPulseIndexL2[pair], pulseToT.GetSize(), i2)) {
        ++malformedPairs;
        continue;
      }
      if (!PassCuts(pulseToT[i1], pulseToT[i2], pairLayerDT[pair], pairECalScore[pair],
                    memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax))
        continue;
      // All five histograms use exactly the same finite selected pairs.
      if (!std::isfinite(pairLEL1[pair]) || !std::isfinite(pairLEL2[pair]) ||
          !std::isfinite(pairMeanLE[pair]) || !std::isfinite(pairECalDT[pair])) {
        ++malformedPairs;
        continue;
      }
      // 3. Add additional study cuts above this common filling block.
      hPairMeanLE.Fill(pairMeanLE[pair]);
      hLayerDT.Fill(pairLayerDT[pair]);
      hECalPairDT.Fill(pairECalDT[pair]);
      hLEL1.Fill(pairLEL1[pair]);
      hLEL2.Fill(pairLEL2[pair]);
    }
  }
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd &&
      (maxEvents < 0 || eventsRead < maxEvents)) {
    std::cerr << "Tree reading stopped early; check branch availability/types.\n";
    return;
  }
  std::cout << "Read " << eventsRead << " events; skipped malformed pairs: " << malformedPairs << '\n';

  // 4. Sigma here is full-sample standard deviation, NOT a Gaussian-core fit.
  for (auto *h : {&hPairMeanLE, &hLayerDT, &hECalPairDT, &hLEL1, &hLEL2})
    Report(*h);
  if (hLayerDT.GetEntries() >= 2)
    std::cout << "Equal-independent-layer estimate SD(dt)/sqrt(2) = "
              << hLayerDT.GetStdDev()/std::sqrt(2.0) << " ns (assumption-dependent).\n";
  std::cout << "Widths are conditioned by stored-pair and study cuts; ECal residuals also include reference timing.\n";

  // 5. Draw copies so histograms remain on interactive canvases after return.
  auto *timingCanvas = new TCanvas("cPairTiming_"+tag, title+" | "+cuts, 1500, 500);
  timingCanvas->Divide(3, 1);
  int pad = 0;
  for (auto *h : {&hPairMeanLE, &hLayerDT, &hECalPairDT}) {
    timingCanvas->cd(++pad);
    h->DrawCopy("HIST");
    Annotate(*h);
  }
  timingCanvas->Update();

  // 6. Matched layer spectra, with separate standard deviations.
  auto *layerCanvas = new TCanvas("cPairLayers_"+tag, title+" | "+cuts, 900, 650);
  hLEL1.SetLineColor(kBlue+1);
  hLEL2.SetLineColor(kRed+1);
  hLEL1.SetMaximum(1.3*std::max(1.0, std::max(hLEL1.GetMaximum(), hLEL2.GetMaximum())));
  auto *drawL1 = hLEL1.DrawCopy("HIST");
  auto *drawL2 = hLEL2.DrawCopy("HIST SAME");
  auto *legend = new TLegend(0.40, 0.73, 0.89, 0.89);
  legend->SetBorderSize(0);
  legend->AddEntry(drawL1, "L1: "+WidthLabel(hLEL1), "l");
  legend->AddEntry(drawL2, "L2 partner: "+WidthLabel(hLEL2), "l");
  legend->Draw();
  layerCanvas->Update();
  if (outputPrefix && outputPrefix[0]) {
    timingCanvas->SaveAs(TString(outputPrefix)+"_timing.pdf");
    layerCanvas->SaveAs(TString(outputPrefix)+"_layers.pdf");
  }
}
