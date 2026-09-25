#include <TCanvas.h>
#include <TEnv.h>
#include <THashList.h>
#include <set>
#include <stdexcept>
#include <string>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <array>
#include <fstream>
#include <memory>
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
// Preferred configuration interface:
// PlotPairingVariables("PlotPairingVariables.conf");
// Optionally override only its input directory:
// PlotPairingVariables("PlotPairingVariables.conf", "/path/to/Rootfiles");
// Positional arguments after maxEvents: L1 bar, bin width, LE min/max,
// layer-dt min/max, ECal-dt min/max, member ToT min/max, maximum |L2-L1|,
// maximum ECal ellipse radius, savePlots, outputDirectory, minEntriesPerBar.
// -1 disables either additional pair cut. Default selects stored pairs as-is.
// One event pass produces all three canvases, including the all-bar scan.
// savePlots defaults to false. When true, save all three as PNG AND PDF plus
// the per-bar CSV, under outputDirectory with run/bar-specific file names.
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
                         bool savePlots = false, const char *outputDirectory = "pairing_plots",
                         int minEntriesPerBar = 30) {
  using namespace CDetPairingPlots;
  const int nLE = Bins(binWidthNs, leMinNs, leMaxNs);
  const int nDT = Bins(binWidthNs, dtMinNs, dtMaxNs);
  const int nECalDT = Bins(binWidthNs, ecalDTMinNs, ecalDTMaxNs);
  if ((savePlots && (!outputDirectory || !outputDirectory[0])) || minEntriesPerBar < 2 || !nLE || !nDT || !nECalDT || selectedLayer1Bar < 0 || selectedLayer1Bar >= 84 ||
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

  // All 168 paired-member histograms are filled before the selected-bar cut.
  std::array<std::unique_ptr<TH1D>, 168> histograms;
  for (int bar = 0; bar < 168; ++bar) {
    histograms[bar].reset(new TH1D(TString::Format("hPairMemberDT_%s_bar%d", tag.Data(), bar), TString::Format("Run %d, global half-bar %d;t_{ECal} - t_{CDet,member} (ns);Paired members", runNumber, bar), nECalDT, ecalDTMinNs, ecalDTMaxNs));
    Prepare(*histograms[bar]);
  }

  Long64_t eventsRead = 0, malformedPairs = 0;
  const TString cuts = CutLabel(memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax);
  std::cout << cuts << "\nNo additional ECal event cut is imposed.\n";
  while ((maxEvents < 0 || eventsRead < maxEvents) && reader.Next()) {
    ++eventsRead;
    const size_t nPairs = pairPixelL1.GetSize();
    if (pulsePixelID.GetSize() != pulseToT.GetSize() ||
        pulsePixelID.GetSize() != pulseECalDT.GetSize() ||
        pairPixelL2.GetSize() != nPairs || pairPulseIndexL1.GetSize() != nPairs || pairPulseIndexL2.GetSize() != nPairs ||
        pairLEL1.GetSize() != nPairs || pairLEL2.GetSize() != nPairs ||
        pairMeanLE.GetSize() != nPairs || pairLayerDT.GetSize() != nPairs ||
        pairECalDT.GetSize() != nPairs || pairECalScore.GetSize() != nPairs) {
      std::cerr << "Mismatched pair arrays at entry " << reader.GetCurrentEntry() << '\n';
      return;
    }
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
      if (!PassCuts(pulseToT[i1], pulseToT[i2], pairLayerDT[pair], pairECalScore[pair],
                    memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax))
        continue;
      // All five histograms use exactly the same finite selected pairs.
      if (!std::isfinite(pairLEL1[pair]) || !std::isfinite(pairLEL2[pair]) ||
          !std::isfinite(pairMeanLE[pair]) || !std::isfinite(pairECalDT[pair])) {
        ++malformedPairs;
        continue;
      }
      // 3. Detector-wide paired-member timing, under the same study cuts.
      // This is ECal minus each member's LE, NOT ECal minus pair-mean time.
      histograms[id1/16]->Fill(pulseECalDT[i1]);
      histograms[id2/16]->Fill(pulseECalDT[i2]);

      // Only the first two canvases restrict the Layer-1 bar.
      // Keep the actual matched L2 partner regardless of its bar number.
      if (id1/16 != static_cast<size_t>(selectedLayer1Bar))
        continue;
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
  if (savePlots) {
    timingCanvas->SaveAs(outputPrefix+"_timing.pdf");
    timingCanvas->SaveAs(outputPrefix+"_timing.png");
    layerCanvas->SaveAs(outputPrefix+"_layers.pdf");
    layerCanvas->SaveAs(outputPrefix+"_layers.png");
  }

  // 7. Third canvas: sigma versus local half-bar, both layers.
  // Graphs omit low-statistics bars; the optional CSV retains every bar/status.
  auto *g1 = new TGraphErrors();
  auto *g2 = new TGraphErrors();
  g1->SetName("gPairSigmaL1_"+tag);
  g2->SetName("gPairSigmaL2_"+tag);
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
    Report(h);
    if (valid) {
      auto *g = bar < 84 ? g1 : g2;
      const int point = g->GetN();
      g->SetPoint(point, bar%84, h.GetStdDev());
      g->SetPointError(point, 0, h.GetStdDevError());
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
  auto *canvas = new TCanvas("cPairSigmaVsBar_"+tag, TString::Format("Run %d | ", runNumber)+cuts, 1200, 700);
  double ymax = 1;
  for (auto *g : {g1, g2})
    for (int i = 0; i < g->GetN(); ++i)
      ymax = std::max(ymax, g->GetPointY(i)+g->GetErrorY(i));
  canvas->DrawFrame(-0.5, 0, 83.5, 1.2*ymax, TString::Format("Run %d: paired-member timing widths;Half-bar within layer;SD(t_{ECal}-t_{CDet,member}) (ns)", runNumber));
  g1->SetMarkerStyle(20); g1->SetMarkerColor(kBlue+1); g1->SetLineColor(kBlue+1);
  g2->SetMarkerStyle(22); g2->SetMarkerColor(kRed+1); g2->SetLineColor(kRed+1);
  if (g1->GetN()) g1->Draw("P SAME");
  if (g2->GetN()) g2->Draw("P SAME");
  auto *barLegend = new TLegend(0.65, 0.76, 0.89, 0.89);
  barLegend->AddEntry(g1, "Layer 1", "p");
  barLegend->AddEntry(g2, "Layer 2", "p");
  barLegend->Draw();
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
    "analysis.events", "analysis.layer1_bar", "cuts.member_tot_min_ns",
    "cuts.member_tot_max_ns", "cuts.layer_dt_max_ns", "cuts.pair_radius_max",
    "plots.bin_width_ns", "plots.le_min_ns", "plots.le_max_ns",
    "plots.layer_dt_min_ns", "plots.layer_dt_max_ns", "plots.ecal_dt_min_ns",
    "plots.ecal_dt_max_ns", "plots.min_entries_per_bar", "output.save_plots",
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
    const Long64_t minimum = integer("plots.min_entries_per_bar", 30);
    const Long64_t save = integer("output.save_plots", 0);
    if (run <= 0 || run > 2147483647 || bar < 0 || bar >= 84 ||
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
        save != 0, config.GetValue("output.directory", "pairing_plots"), static_cast<int>(minimum));
  } catch (const std::exception &error) {
    std::cerr << "Invalid pairing configuration: " << error.what() << '\n';
  }
}
