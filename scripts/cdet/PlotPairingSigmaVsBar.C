#include <TCanvas.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <array>
#include <fstream>
#include <memory>
#include "CDetRunDataset.h"
#include "PairingPlotHelpers.h"

// .L PlotPairingSigmaVsBar.C+
// PlotPairingSigmaVsBar(6077, "/path/to/Rootfiles", -1);
// Use matching ToT/dt/radius settings when comparing with PlotPairingVariables.
// Example: PlotPairingSigmaVsBar(6077, "/path/to/Rootfiles", -1, 30,
//     0.5, -100, 50, 12, 35, 10, 2, "bars_tight");
// Optional outputs: <prefix>_sigma_vs_bar.pdf and _sigma_vs_bar.csv.
// SD/error use ROOT full-sample moments, including underflow/overflow;
// error bars do not account for correlations between pairs in one event.
// Each member contributes ECal-minus-its-own-LE to its electronic half-bar.
// This is NOT pair-mean timing assigned to both bars. Only stored pairs enter.
void PlotPairingSigmaVsBar(int runNumber = 6077, const char *inputDirectory = nullptr,
                          Long64_t maxEvents = -1, int minEntries = 30,
                          double binWidthNs = 0.5, double dtMinNs = -100, double dtMaxNs = 50,
                          double memberToTMinNs = 0, double memberToTMaxNs = 1e9,
                          double layerDTMaxNs = -1, double pairRadiusMax = -1,
                          const char *outputPrefix = "") {
  using namespace CDetPairingPlots;
  const int nBins = Bins(binWidthNs, dtMinNs, dtMaxNs);
  if (!nBins || minEntries < 2 || maxEvents < -1 ||
      !ValidCuts(memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax)) {
    std::cerr << "Invalid histogram limits, event limit, minimum entries, or cuts.\n";
    return;
  }
  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull() && gSystem->Getenv("OUT_DIR"))
    input = gSystem->Getenv("OUT_DIR");
  if (input.IsNull()) {
    std::cerr << "Supply an input directory or OUT_DIR.\n";
    return;
  }
  TChain chain("T");
  if (CDetRunDataset::AddToChain(&chain, runNumber, input.Data()) <= 0 || chain.GetEntries() == 0)
    return;
  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> tot(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<Double_t> residual(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> index1(reader, "earm.cdet.pair.pulse_index_l1");
  TTreeReaderArray<Double_t> index2(reader, "earm.cdet.pair.pulse_index_l2");
  TTreeReaderArray<Double_t> layerDT(reader, "earm.cdet.pair.dt");
  TTreeReaderArray<Double_t> score(reader, "earm.cdet.pair.ecal_score");

  static unsigned invocation = 0;
  const TString tag = TString::Format("run%d_%u", runNumber, ++invocation);
  std::array<std::unique_ptr<TH1D>, 168> histograms;
  for (int bar = 0; bar < 168; ++bar) {
    histograms[bar].reset(new TH1D(TString::Format("hPairMemberDT_%s_bar%d", tag.Data(), bar), TString::Format("Run %d, global half-bar %d;t_{ECal} - t_{CDet,member} (ns);Paired members", runNumber, bar), nBins, dtMinNs, dtMaxNs));
    Prepare(*histograms[bar]);
  }
  const TString cuts = CutLabel(memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax);
  std::cout << cuts << "\nNo additional ECal event cut is imposed.\n";
  Long64_t eventsRead = 0, malformed = 0;
  while ((maxEvents < 0 || eventsRead < maxEvents) && reader.Next()) {
    ++eventsRead;
    if (pixel.GetSize() != tot.GetSize() || pixel.GetSize() != residual.GetSize() ||
        index1.GetSize() != index2.GetSize() || index1.GetSize() != layerDT.GetSize() ||
        index1.GetSize() != score.GetSize()) {
      std::cerr << "Mismatched pulse/pair arrays at entry " << reader.GetCurrentEntry() << '\n';
      return;
    }
    for (size_t pair = 0; pair < index1.GetSize(); ++pair) {
      size_t i1, i2, id1, id2;
      if (!Index(index1[pair], pixel.GetSize(), i1) || !Index(index2[pair], pixel.GetSize(), i2) ||
          !Index(pixel[i1], 2688, id1) || !Index(pixel[i2], 2688, id2) ||
          id1 >= 1344 || id2 < 1344 ||
          !std::isfinite(residual[i1]) || !std::isfinite(residual[i2])) {
        ++malformed;
        continue;
      }
      if (!PassCuts(tot[i1], tot[i2], layerDT[pair], score[pair],
                    memberToTMinNs, memberToTMaxNs, layerDTMaxNs, pairRadiusMax))
        continue;
      histograms[id1/16]->Fill(residual[i1]);
      histograms[id2/16]->Fill(residual[i2]);
    }
  }
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd &&
      (maxEvents < 0 || eventsRead < maxEvents)) {
    std::cerr << "Tree reading stopped early; no completed results produced.\n";
    return;
  }
  std::cout << "Read " << eventsRead << " events; skipped malformed pairs: " << malformed << '\n';

  // Graphs omit low-statistics bars; the optional CSV retains every bar/status.
  auto *g1 = new TGraphErrors();
  auto *g2 = new TGraphErrors();
  g1->SetName("gPairSigmaL1_"+tag);
  g2->SetName("gPairSigmaL2_"+tag);
  std::ofstream csv;
  if (outputPrefix && outputPrefix[0]) {
    csv.open((TString(outputPrefix)+"_sigma_vs_bar.csv").Data());
    if (!csv) std::cerr << "Could not open CSV output.\n";
    else {
      csv.precision(12);
      csv << "global_bar,layer,local_bar,entries,mean_ns,stddev_ns,stddev_error_ns,underflow,overflow,status\n";
    }
  }
  for (int bar = 0; bar < 168; ++bar) {
    const auto &h = *histograms[bar];
    const bool valid = h.GetEntries() >= minEntries;
    Report(h);
    if (valid) {
      auto *g = bar < 84 ? g1 : g2;
      const int point = g->GetN();
      g->SetPoint(point, bar%84, h.GetStdDev());
      g->SetPointError(point, 0, h.GetStdDevError());
    }
    if (csv) {
      csv << bar << ',' << bar/84+1 << ',' << bar%84 << ',' << h.GetEntries() << ',';
      if (h.GetEntries() >= 2)
        csv << h.GetMean() << ',' << h.GetStdDev() << ',' << h.GetStdDevError();
      else csv << ",,";
      csv << ',' << h.GetBinContent(0) << ',' << h.GetBinContent(nBins+1)
          << ',' << (valid ? "accepted" : "insufficient_entries") << '\n';
    }
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
  auto *legend = new TLegend(0.65, 0.76, 0.89, 0.89);
  legend->AddEntry(g1, "Layer 1", "p");
  legend->AddEntry(g2, "Layer 2", "p");
  legend->Draw();
  canvas->Update();
  if (outputPrefix && outputPrefix[0])
    canvas->SaveAs(TString(outputPrefix)+"_sigma_vs_bar.pdf");
  std::cout << "Widths include overflow and selected-sample tails; they are not Gaussian fit sigmas or intrinsic resolutions.\n";
}
