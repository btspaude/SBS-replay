// Standalone pulse timing comparison of an investigation bar and a control bar.
// From scripts/cdet, with OUT_DIR set to the replay ROOT-file directory:
//   root -l
//   .L Plot_CDet_PulsePixelTiming.C+
//   Plot_CDet_PulsePixelTiming("CDet_run6077_pulse_pixel_timing.conf");
// Override the investigation bar, input directory, or event limit if desired:
//   Plot_CDet_PulsePixelTiming("CDet_run6077_pulse_pixel_timing.conf", 19, "/path/to/Rootfiles", 10000);
// Numeric shorthand reads CDet_run<run>_pulse_pixel_timing.conf:
//   Plot_CDet_PulsePixelTiming(6077, 19);
//
// Global bar IDs are 0..167; bar 30 contains logical pixels 480..495.
// pulse.ecal_residual is already t_ECal - t_CDet,corr in ns, including the
// calibration used during replay. A zero pixel offset still gives a value.
// Each bar has a 4x4 pixel canvas and a summed-pulse timing canvas. Fits use
// gaus(0)+pol1(3); their centroids are residual timing, not absolute offsets.
// The configuration supplies the run, event/segment limits, ECal energy cut,
// bar choices, and fit settings. Pulse-quality/projection flags remain those
// stored during replay. Manual polygons are not applied. No calibration files
// are written and no offsets are activated by this macro.

#include <TCanvas.h>
#include <TChain.h>
#include <TEnv.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TList.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <iostream>

#include "CDetRunDataset.h"

// Available at the ROOT prompt after the function returns.
// Indices 0..15 refer to the pixels of the most recently initialized bar.
TH1D* hPixelDt[16] = {};
TH1D* hControlPixelDt[16] = {};
TH1D* hBarDt = nullptr;
TH1D* hControlBarDt = nullptr;

namespace CDetPulsePixelTiming {

// Keep an accepted fit attached to its histogram for inspection from ROOT.
// Failed or boundary-limited fits are reported rather than drawn as results.
TF1* FitSpectrum(TH1D* hist, double fitMin, double fitMax,
                 int minFitEntries, TString& fitNote)
{
  const int firstBin = std::max(1, hist->FindFixBin(fitMin));
  const int lastBin = std::min(hist->GetNbinsX(), hist->FindFixBin(std::nextafter(fitMax, fitMin)));
  const double fitEntries = hist->Integral(firstBin, lastBin);
  if (fitEntries < minFitEntries) {
    fitNote = "Fit: too few entries in range";
    std::cout << "[CDet timing fit] " << hist->GetName() << ": "
              << fitEntries << " entries in fit window; need " << minFitEntries << '\n';
    return nullptr;
  }

  int peakBin = firstBin;
  for (int bin = firstBin + 1; bin <= lastBin; ++bin)
    if (hist->GetBinContent(bin) > hist->GetBinContent(peakBin)) peakBin = bin;
  const double peak = hist->GetBinContent(peakBin);
  const double background = 0.5 * (hist->GetBinContent(firstBin) + hist->GetBinContent(lastBin));
  const double sigmaMin = 0.5 * hist->GetBinWidth(1);
  const double sigmaMax = 0.5 * (fitMax - fitMin);
  const TString name = TString(hist->GetName()) + "_gaus_linear";
  TF1 model(name, "gaus(0)+pol1(3)", fitMin, fitMax, TF1::EAddToList::kNo);
  model.SetParameters(std::max(1.0, peak - background), hist->GetBinCenter(peakBin),
                      std::max(sigmaMin, std::min(3.0, sigmaMax)), background, 0.0);
  model.SetParLimits(0, 0.0, 2.0 * std::max(1.0, peak));
  model.SetParLimits(1, fitMin, fitMax);
  model.SetParLimits(2, sigmaMin, sigmaMax);
  const TFitResultPtr result = hist->Fit(&model, "RQSN0");
  const int status = int(result);
  const double mean = model.GetParameter(1), sigma = model.GetParameter(2);
  const double meanError = model.GetParError(1), sigmaError = model.GetParError(2);
  const double tolerance = 1.0e-3 * (fitMax - fitMin);
  if (status != 0 || !result.Get() || !result->IsValid() || model.GetNDF() <= 0 ||
      !std::isfinite(mean) || !std::isfinite(sigma) ||
      !std::isfinite(meanError) || !std::isfinite(sigmaError) ||
      meanError <= 0.0 || sigmaError <= 0.0) {
    fitNote = "Fit: failed";
    std::cout << "[CDet timing fit] " << hist->GetName() << ": failed, status " << status << '\n';
    return nullptr;
  }
  if (mean <= fitMin + tolerance || mean >= fitMax - tolerance ||
      sigma <= sigmaMin + tolerance || sigma >= sigmaMax - tolerance) {
    fitNote = "Fit: at parameter limit";
    std::cout << "[CDet timing fit] " << hist->GetName() << ": at parameter limit\n";
    return nullptr;
  }

  TF1* fit = static_cast<TF1*>(model.Clone());
  fit->AddToGlobalList(false);
  fit->SetLineColor(kRed + 1);
  fit->SetLineWidth(2);
  hist->GetListOfFunctions()->Add(fit);
  std::cout << "[CDet timing fit] " << hist->GetName() << ": mu = " << mean
            << " +/- " << meanError << " ns, sigma = " << sigma << " +/- "
            << sigmaError << " ns, chi2/ndf = " << model.GetChisquare()
            << '/' << model.GetNDF() << '\n';
  return fit;
}

void DrawSpectrum(TH1D* hist, double fitMin, double fitMax, int minFitEntries)
{
  TString fitNote;
  TF1* fit = FitSpectrum(hist, fitMin, fitMax, minFitEntries, fitNote);
  hist->Draw("HIST");
  if (fit) fit->Draw("SAME");

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(33);
  label.SetTextFont(42);
  label.SetTextSize(0.035);
  label.DrawLatex(0.87, 0.85, TString::Format("N = %.0f pulses", hist->GetEntries()));
  if (fit) {
    label.DrawLatex(0.87, 0.79, TString::Format("Fit #mu = %.2f #pm %.2f ns", fit->GetParameter(1), fit->GetParError(1)));
    label.DrawLatex(0.87, 0.73, TString::Format("Fit #sigma = %.2f #pm %.2f ns", fit->GetParameter(2), fit->GetParError(2)));
  } else {
    label.DrawLatex(0.87, 0.79, fitNote);
  }
}

} // namespace CDetPulsePixelTiming

void Plot_CDet_PulsePixelTiming(
    const char *configFile = "CDet_run6077_pulse_pixel_timing.conf", int bar = -1,
    const char *inputDirectory = nullptr, Long64_t eventsOverride = -2,
    double ecalDTBinWidthNs = 1.0,
    double ecalDTMinNs = -60.0, double ecalDTMaxNs = 30.0)
{
  using namespace CDetPulsePixelTiming;
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0) {
    std::cerr << "[CDet pulse pixel timing] Cannot read the configuration file.\n";
    return;
  }
  if (config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pulse pixel timing] Expected config.version: 1.\n";
    return;
  }
  const int runNumber = config.GetValue("analysis.run_number", 0);
  // -2 uses the configuration; -1 reads all entries; 0 or more limits entries.
  const Long64_t events = eventsOverride == -2
      ? config.GetValue("analysis.events", -1) : eventsOverride;
  const int minSegment = config.GetValue("analysis.min_segment", -1);
  const int maxSegment = config.GetValue("analysis.max_segment", -1);
  const double energyMin = config.GetValue("analysis.ecal_energy_min", NAN);
  const double energyMax = config.GetValue("analysis.ecal_energy_max", NAN);
  if (bar == -1) bar = config.GetValue("pulse_pixel_timing.study_bar", 19);
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double fitMin = config.GetValue("pulse_pixel_timing.fit_min_ns", -55.0);
  const double fitMax = config.GetValue("pulse_pixel_timing.fit_max_ns", -10.0);
  const int minFitEntries = config.GetValue("pulse_pixel_timing.min_fit_entries", 35);
  const int pixelsPerBar = 16;
  const int nBars = 168;
  if (!std::isfinite(energyMin) || !std::isfinite(energyMax) || energyMin >= energyMax) {
    std::cerr << "[CDet pulse pixel timing] Set analysis.ecal_energy_min/max to a valid GeV interval.\n";
    return;
  }
  if (runNumber <= 0 || bar < 0 || bar >= nBars ||
      controlBar < 0 || controlBar >= nBars || events < -1 ||
      !std::isfinite(ecalDTBinWidthNs) || ecalDTBinWidthNs <= 0.0 ||
      !std::isfinite(ecalDTMinNs) || !std::isfinite(ecalDTMaxNs) ||
      ecalDTMinNs >= ecalDTMaxNs || !std::isfinite(fitMin) || !std::isfinite(fitMax) ||
      fitMin < ecalDTMinNs || fitMax > ecalDTMaxNs ||
      fitMax - fitMin <= 2.0 * ecalDTBinWidthNs || minFitEntries < 6) {
    std::cerr << "[CDet pulse pixel timing] Invalid run, bar, event limit, binning, or fit settings.\n";
    return;
  }
  const int firstPixel = pixelsPerBar * bar;
  const int firstControlPixel = pixelsPerBar * controlBar;
  const int nECalDT = int(std::ceil((ecalDTMaxNs - ecalDTMinNs) / ecalDTBinWidthNs));

  // 1. Load one run, including its ROOT rollover files.
  const char *directory = inputDirectory && inputDirectory[0]
      ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0]) {
    std::cerr << "[CDet pulse pixel timing] Set OUT_DIR or pass inputDirectory.\n";
    return;
  }
  TChain chain("T");
  if (CDetRunDataset::AddToChain(&chain, runNumber, directory, minSegment, maxSegment) <= 0) return;
  const Long64_t total = chain.GetEntries();
  if (total <= 0) {
    std::cerr << "[CDet pulse pixel timing] No events found.\n";
    return;
  }
  // -1 reads all entries; a nonnegative value limits the number processed.
  const Long64_t limit = events >= 0 ? std::min(events, total) : total;

  // 2. Bind the initial pulse branches. No pair.* branches are needed.

  TTreeReader reader(&chain);
  TTreeReaderArray<double> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<double> pulseDt(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<double> calibValid(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderValue<double> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<double> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderArray<double> cdetPulseX(reader, "earm.cdet.pulse.x");
  TTreeReaderArray<double> cdetPulseXCorr(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<double> broadQuality(reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<double> ecalEligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<double> spatialPass(reader, "earm.cdet.pulse.spatial_pass");


  // 3. Define the investigation/control pixel histograms and bar sums.
  // hPixelDt[i] corresponds to logical pixel firstPixel + i.
  for (int i = 0; i < pixelsPerBar; ++i) {
    const int pixelID = firstPixel + i;
    const TString name = TString::Format("hPixelDt_run%d_pixel%d", runNumber, pixelID);
    const TString label = TString::Format("t_{ECal} - t_{CDet,%d}", pixelID);
    const TString title = TString::Format("Bar %d: %s;%s (ns);Counts", bar, label.Data(), label.Data());
    delete hPixelDt[i]; // Replace the previous histogram when rerunning.
    hPixelDt[i] = new TH1D(name, title, nECalDT, ecalDTMinNs, ecalDTMaxNs);
    hPixelDt[i]->SetDirectory(nullptr);
    hPixelDt[i]->SetStats(false);

    const int controlPixelID = firstControlPixel + i;
    const TString controlName = TString::Format("hControlPixelDt_run%d_pixel%d", runNumber, controlPixelID);
    const TString controlLabel = TString::Format("t_{ECal} - t_{CDet,%d}", controlPixelID);
    const TString controlTitle = TString::Format("Control bar %d: %s;%s (ns);Counts", controlBar, controlLabel.Data(), controlLabel.Data());
    delete hControlPixelDt[i];
    hControlPixelDt[i] = new TH1D(controlName, controlTitle, nECalDT, ecalDTMinNs, ecalDTMaxNs);
    hControlPixelDt[i]->SetDirectory(nullptr);
    hControlPixelDt[i]->SetStats(false);
  }
  delete hBarDt;
  delete hControlBarDt;
  hBarDt = new TH1D(TString::Format("hBarDt_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: all 16 pixels;t_{ECal} - t_{CDet} (ns);Counts", runNumber, bar), nECalDT, ecalDTMinNs, ecalDTMaxNs);
  hControlBarDt = new TH1D(TString::Format("hControlBarDt_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: all 16 pixels;t_{ECal} - t_{CDet} (ns);Counts", runNumber, controlBar), nECalDT, ecalDTMinNs, ecalDTMaxNs);
  for (TH1D* hist : {hBarDt, hControlBarDt}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

  std::cout << "[CDet pulse pixel timing] Configuration: " << configFile << '\n'
            << "  Run " << runNumber << ", investigation bar " << bar
            << ", pixels " << firstPixel << '-' << firstPixel + pixelsPerBar - 1
            << "; control bar " << controlBar << ", pixels " << firstControlPixel
            << '-' << firstControlPixel + pixelsPerBar - 1 << ".\n"
            << "  Reading " << limit << " entries; ECal energy: " << energyMin
            << " < E < " << energyMax << " GeV.\n"
            << "  DT binning: " << nECalDT << " bins from "
            << ecalDTMinNs << " to " << ecalDTMaxNs << " ns.\n"
            << "  Gaussian + linear fit window: " << fitMin << " to " << fitMax
            << " ns; minimum entries: " << minFitEntries << ".\n";

  // 4. Loop over events and individual pulses in the requested bar.
  Long64_t processed = 0, energyEvents = 0, selectedPulses = 0, controlPulses = 0;
  while (processed < limit && reader.Next()) {
    ++processed;
    if (std::isfinite(*ecalEnergy) && *ecalEnergy > energyMin && *ecalEnergy < energyMax) {
      ++energyEvents;
      for (size_t i = 0; i < pixel.GetSize(); ++i) {
        if (!std::isfinite(pixel[i]) || !std::isfinite(pulseDt[i])) continue;
        const int pixelID = int(pixel[i]);
        const int iPixel = pixelID - firstPixel;
        const int iControlPixel = pixelID - firstControlPixel;
        const bool inBar = iPixel >= 0 && iPixel < pixelsPerBar;
        const bool inControlBar = iControlPixel >= 0 && iControlPixel < pixelsPerBar;
        if (!inBar && !inControlBar) continue;

        // Apply the same stored flags to both bars. No pairing or polygons.
        if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
              ecalEligible[i] > 0.5 && spatialPass[i] > 0.5)) continue;
        if (inBar) {
          ++selectedPulses;
          hPixelDt[iPixel]->Fill(pulseDt[i]);
        }
        if (inControlBar) {
          ++controlPulses;
          hControlPixelDt[iControlPixel]->Fill(pulseDt[i]);
        }
      }
    }

    if (processed % 1000 == 0)
      std::cout << "[CDet pulse pixel timing] " << processed << '/' << limit
                << " events; ECal energy passed: " << energyEvents
                << "; study/control pulses: " << selectedPulses << '/' << controlPulses << '\n';
  } //end event loop
  if (processed != limit) {
    std::cerr << "[CDet pulse pixel timing] Reader stopped after " << processed
              << " of " << limit << " requested events; status "
              << reader.GetEntryStatus() << '\n';
    return;
  }
  // 5. Sum all 16 pixel spectra, including their underflow/overflow counts.
  // These are pulse distributions, not an average time per event or bar.
  for (int i = 0; i < pixelsPerBar; ++i) {
    hBarDt->Add(hPixelDt[i]);
    hControlBarDt->Add(hControlPixelDt[i]);
  }

  // 6. Fit and draw both bars: one pixel canvas and one summed canvas each.
  TCanvas *cPixelDt = new TCanvas(TString::Format("cPixelDt_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: pulse timing", runNumber, bar), 1400, 1000);
  cPixelDt->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cPixelDt->cd(i + 1);
    DrawSpectrum(hPixelDt[i], fitMin, fitMax, minFitEntries);
  }
  cPixelDt->Modified();
  cPixelDt->Update();

  TCanvas *cBarDt = new TCanvas(TString::Format("cBarDt_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: summed pulse timing", runNumber, bar), 900, 650);
  DrawSpectrum(hBarDt, fitMin, fitMax, minFitEntries);
  cBarDt->Modified();
  cBarDt->Update();

  TCanvas *cControlPixelDt = new TCanvas(TString::Format("cControlPixelDt_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: pulse timing", runNumber, controlBar), 1400, 1000);
  cControlPixelDt->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cControlPixelDt->cd(i + 1);
    DrawSpectrum(hControlPixelDt[i], fitMin, fitMax, minFitEntries);
  }
  cControlPixelDt->Modified();
  cControlPixelDt->Update();

  TCanvas *cControlBarDt = new TCanvas(TString::Format("cControlBarDt_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: summed pulse timing", runNumber, controlBar), 900, 650);
  DrawSpectrum(hControlBarDt, fitMin, fitMax, minFitEntries);
  cControlBarDt->Modified();
  cControlBarDt->Update();

  std::cout << "[CDet pulse pixel timing] Read " << processed << " events; "
            << energyEvents << " passed the ECal energy cut; " << selectedPulses
            << " pulses in bar " << bar << " and " << controlPulses
            << " in control bar " << controlBar << ". Drew four canvases.\n";
}

// Retain the run-number shorthand using the separate pulse-timing config.
void Plot_CDet_PulsePixelTiming(
    int runNumber, int bar = -1, const char *inputDirectory = nullptr,
    Long64_t eventsOverride = -2, double ecalDTBinWidthNs = 1.0,
    double ecalDTMinNs = -60.0, double ecalDTMaxNs = 30.0)
{
  const TString configFile = TString::Format("CDet_run%d_pulse_pixel_timing.conf", runNumber);
  Plot_CDet_PulsePixelTiming(configFile.Data(), bar, inputDirectory, eventsOverride,
                           ecalDTBinWidthNs, ecalDTMinNs, ecalDTMaxNs);
}
