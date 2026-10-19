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
// Each bar also has a dt-versus-dx canvas before the optional pulse ellipse,
// and an after panel only when good_pulse.single_cut_radius > 0. Radius -1
// disables the ellipse. The 1D spectra use the selected pulses when enabled.
// This is a pulse-level analysis cut, not exclusive offline single-layer recovery.
// LE and LE-versus-ToT plots use the same selected pulses, before fitted pixel
// offsets/time walk. A detector-wide dt-versus-dx canvas covers both layers.
// The configuration supplies explicit LE and optional ToT cuts; broad_quality_pass
// is NOT required, so its replay ToT window can be widened or disabled here.
// Calibration, ECal eligibility, and spatial flags remain those stored in replay.
// Manual polygons are not applied. No calibration files or offsets are written.

#include <TCanvas.h>
#include <TChain.h>
#include <TEnv.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TList.h>
#include <TPad.h>
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
TH2D* hBarDtVsDxBefore = nullptr;
TH2D* hBarDtVsDxAfter = nullptr;
TH2D* hControlBarDtVsDxBefore = nullptr;
TH2D* hControlBarDtVsDxAfter = nullptr;
TH2D* hDetectorDtVsDxBefore = nullptr;
TH2D* hDetectorDtVsDxAfter = nullptr;
TH1D* hPixelLe[16] = {};
TH1D* hControlPixelLe[16] = {};
TH2D* hPixelLeVsToT[16] = {};
TH2D* hControlPixelLeVsToT[16] = {};
TH1D* hBarLe = nullptr;
TH1D* hControlBarLe = nullptr;
TH2D* hBarLeVsToT = nullptr;
TH2D* hControlBarLeVsToT = nullptr;

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
  const double outside = hist->GetEntries() - hist->Integral();
  if (outside > 0.0)
    label.DrawLatex(0.87, 0.67, TString::Format("Outside axes: %.0f", outside));
}

void DrawPulseDistribution(TH1* hist, bool twoDimensional)
{
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(twoDimensional ? 0.18 : 0.05);
  if (twoDimensional) {
    hist->SetMinimum(0.0);
    hist->SetMaximum(std::max(1.0, hist->GetMaximum()));
    hist->GetZaxis()->SetLabelSize(0.03);
    hist->GetZaxis()->SetTitleSize(0.035);
    hist->GetZaxis()->SetTitleOffset(1.5);
  }
  hist->Draw(twoDimensional ? "COLZ" : "HIST");
  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.035);
  if (!twoDimensional) label.SetTextAlign(33);
  const double labelX = twoDimensional ? 0.17 : 0.92;
  label.DrawLatex(labelX, 0.85, TString::Format("N = %.0f pulses", hist->GetEntries()));
  const double outside = hist->GetEntries() - hist->Integral();
  if (outside > 0.0)
    label.DrawLatex(labelX, 0.79, TString::Format("Outside axes: %.0f", outside));
}

void DrawDtVsDx(TH2D* hist, bool useEllipse, double xCenter, double timeCenter,
                double xScale, double timeScale, double radius)
{
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.18);
  hist->GetZaxis()->SetLabelSize(0.03);
  hist->GetZaxis()->SetTitleSize(0.035);
  hist->GetZaxis()->SetTitleOffset(1.5);
  hist->Draw("COLZ");
  if (useEllipse) {
    TEllipse boundary(xCenter, timeCenter, radius * xScale, radius * timeScale);
    boundary.SetFillStyle(0);
    boundary.SetLineColor(kRed + 1);
    boundary.SetLineWidth(2);
    boundary.DrawClone("SAME");
  }
  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.03);
  label.DrawLatex(0.17, 0.85, TString::Format("N = %.0f pulses", hist->GetEntries()));
  const double outside = hist->GetEntries() - hist->Integral();
  if (outside > 0.0)
    label.DrawLatex(0.17, 0.79, TString::Format("Outside axes: %.0f", outside));
}

} // namespace CDetPulsePixelTiming

void Plot_CDet_PulsePixelTiming(
    const char *configFile = "CDet_run6077_pulse_pixel_timing.conf", int bar = -1,
    const char *inputDirectory = nullptr, Long64_t eventsOverride = -2,
    double ecalDTBinWidthNs = 1.0,
    double ecalDTMinNs = -150.0, double ecalDTMaxNs = 150.0,
    double ecalDXBinWidthM = 0.005, double ecalDXMinM = -0.10, double ecalDXMaxM = 0.10,
    double dtVsDxMinNs = -150.0, double dtVsDxMaxNs = 150.0,
    double leBinWidthNs = 1.0, double leMinNs = -20.0, double leMaxNs = 100.0,
    double totBinWidthNs = 1.0, double totMinNs = 0.0, double totMaxNs = 100.0)
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
  const double leCutMinNs = config.GetValue("analysis.le_min", 0.02);
  const double leCutMaxNs = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_pixel_timing.apply_tot_cut", 0) != 0;
  const double totCutMinNs = config.GetValue("analysis.tot_min", 4.0);
  const double totCutMaxNs = config.GetValue("analysis.tot_max", 30.0);
  if (bar == -1) bar = config.GetValue("pulse_pixel_timing.study_bar", 19);
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double fitMin = config.GetValue("pulse_pixel_timing.fit_min_ns", -55.0);
  const double fitMax = config.GetValue("pulse_pixel_timing.fit_max_ns", -10.0);
  const int minFitEntries = config.GetValue("pulse_pixel_timing.min_fit_entries", 35);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const bool useEllipse = ellipseRadius > 0.0;
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const int pixelsPerBar = 16;
  const int nBars = 168;
  if (!std::isfinite(energyMin) || !std::isfinite(energyMax) || energyMin >= energyMax) {
    std::cerr << "[CDet pulse pixel timing] Set analysis.ecal_energy_min/max to a valid GeV interval.\n";
    return;
  }
  if (!std::isfinite(leCutMinNs) || !std::isfinite(leCutMaxNs) || leCutMinNs >= leCutMaxNs ||
      (applyToTCut && (!std::isfinite(totCutMinNs) || !std::isfinite(totCutMaxNs) ||
       totCutMinNs < 0.0 || totCutMinNs >= totCutMaxNs))) {
    std::cerr << "[CDet pulse pixel timing] Invalid analysis.le_min/max or enabled analysis.tot_min/max (ns).\n";
    return;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse pixel timing] Ellipse radius must be -1 (off) or positive; enabled ellipses require finite centers and positive finite scales.\n";
    return;
  }
  if (runNumber <= 0 || bar < 0 || bar >= nBars ||
      controlBar < 0 || controlBar >= nBars || events < -1 ||
      !std::isfinite(ecalDTBinWidthNs) || ecalDTBinWidthNs <= 0.0 ||
      !std::isfinite(ecalDTMinNs) || !std::isfinite(ecalDTMaxNs) ||
      ecalDTMinNs >= ecalDTMaxNs || !std::isfinite(fitMin) || !std::isfinite(fitMax) ||
      fitMin < ecalDTMinNs || fitMax > ecalDTMaxNs ||
      fitMax - fitMin <= 2.0 * ecalDTBinWidthNs || minFitEntries < 6 ||
      !std::isfinite(ecalDXBinWidthM) || ecalDXBinWidthM <= 0.0 ||
      !std::isfinite(ecalDXMinM) || !std::isfinite(ecalDXMaxM) || ecalDXMinM >= ecalDXMaxM ||
      !std::isfinite(dtVsDxMinNs) || !std::isfinite(dtVsDxMaxNs) || dtVsDxMinNs >= dtVsDxMaxNs ||
      !std::isfinite(leBinWidthNs) || leBinWidthNs <= 0.0 ||
      !std::isfinite(leMinNs) || !std::isfinite(leMaxNs) || leMinNs >= leMaxNs ||
      !std::isfinite(totBinWidthNs) || totBinWidthNs <= 0.0 ||
      !std::isfinite(totMinNs) || !std::isfinite(totMaxNs) || totMinNs >= totMaxNs) {
    std::cerr << "[CDet pulse pixel timing] Invalid run, bar, event limit, binning, or fit settings.\n";
    return;
  }
  const int firstPixel = pixelsPerBar * bar;
  const int firstControlPixel = pixelsPerBar * controlBar;
  const int nECalDT = int(std::ceil((ecalDTMaxNs - ecalDTMinNs) / ecalDTBinWidthNs));
  const int nECalDX = int(std::ceil((ecalDXMaxM - ecalDXMinM) / ecalDXBinWidthM));
  const int nDtVsDx = int(std::ceil((dtVsDxMaxNs - dtVsDxMinNs) / ecalDTBinWidthNs));
  const int nLe = int(std::ceil((leMaxNs - leMinNs) / leBinWidthNs));
  const int nToT = int(std::ceil((totMaxNs - totMinNs) / totBinWidthNs));
  const double tdcToNs = 0.01; // Matches timing.tdc_to_ns in this replay DB.
  const char* selectionLabel = useEllipse ? "ellipse selected" : "ellipse off";

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
  // Before fitted pixel offsets/time walk; multiply by 0.01 for ns in this replay.
  TTreeReaderArray<double> pulseLe(reader, "earm.cdet.pulse.tdc_le");
  TTreeReaderArray<double> pulseToTNs(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<double> pulseDt(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<double> pulseDx(reader, "earm.cdet.pulse.ecal_x_residual");
  TTreeReaderArray<double> calibValid(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderValue<double> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<double> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderArray<double> cdetPulseX(reader, "earm.cdet.pulse.x");
  TTreeReaderArray<double> cdetPulseXCorr(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<double> ecalEligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<double> spatialPass(reader, "earm.cdet.pulse.spatial_pass");


  // 3. Define the investigation/control pixel histograms and bar sums.
  // hPixelDt[i] corresponds to logical pixel firstPixel + i.
  for (int i = 0; i < pixelsPerBar; ++i) {
    const int pixelID = firstPixel + i;
    const TString name = TString::Format("hPixelDt_run%d_pixel%d", runNumber, pixelID);
    const TString label = TString::Format("t_{ECal} - t_{CDet,%d}", pixelID);
    const TString title = TString::Format("Bar %d: %s (%s);%s (ns);Counts", bar, label.Data(), selectionLabel, label.Data());
    delete hPixelDt[i]; // Replace the previous histogram when rerunning.
    hPixelDt[i] = new TH1D(name, title, nECalDT, ecalDTMinNs, ecalDTMaxNs);
    hPixelDt[i]->SetDirectory(nullptr);
    hPixelDt[i]->SetStats(false);

    const int controlPixelID = firstControlPixel + i;
    const TString controlName = TString::Format("hControlPixelDt_run%d_pixel%d", runNumber, controlPixelID);
    const TString controlLabel = TString::Format("t_{ECal} - t_{CDet,%d}", controlPixelID);
    const TString controlTitle = TString::Format("Control bar %d: %s (%s);%s (ns);Counts", controlBar, controlLabel.Data(), selectionLabel, controlLabel.Data());
    delete hControlPixelDt[i];
    hControlPixelDt[i] = new TH1D(controlName, controlTitle, nECalDT, ecalDTMinNs, ecalDTMaxNs);
    hControlPixelDt[i]->SetDirectory(nullptr);
    hControlPixelDt[i]->SetStats(false);

    delete hPixelLe[i];
    delete hControlPixelLe[i];
    delete hPixelLeVsToT[i];
    delete hControlPixelLeVsToT[i];
    hPixelLe[i] = new TH1D(TString::Format("hPixelLe_run%d_pixel%d", runNumber, pixelID), TString::Format("Bar %d, pixel %d: LE (%s);t_{LE} (ns);Counts", bar, pixelID, selectionLabel), nLe, leMinNs, leMaxNs);
    hControlPixelLe[i] = new TH1D(TString::Format("hControlPixelLe_run%d_pixel%d", runNumber, controlPixelID), TString::Format("Control bar %d, pixel %d: LE (%s);t_{LE} (ns);Counts", controlBar, controlPixelID, selectionLabel), nLe, leMinNs, leMaxNs);
    hPixelLeVsToT[i] = new TH2D(TString::Format("hPixelLeVsToT_run%d_pixel%d", runNumber, pixelID), TString::Format("Bar %d, pixel %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", bar, pixelID, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
    hControlPixelLeVsToT[i] = new TH2D(TString::Format("hControlPixelLeVsToT_run%d_pixel%d", runNumber, controlPixelID), TString::Format("Control bar %d, pixel %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", controlBar, controlPixelID, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
    TH1* pixelPulsePlots[] = {hPixelLe[i], hControlPixelLe[i], hPixelLeVsToT[i], hControlPixelLeVsToT[i]};
    for (TH1* hist : pixelPulsePlots) {
      hist->SetDirectory(nullptr);
      hist->SetStats(false);
    }
  }
  delete hBarDt;
  delete hControlBarDt;
  hBarDt = new TH1D(TString::Format("hBarDt_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: all 16 pixels (%s);t_{ECal} - t_{CDet} (ns);Counts", runNumber, bar, selectionLabel), nECalDT, ecalDTMinNs, ecalDTMaxNs);
  hControlBarDt = new TH1D(TString::Format("hControlBarDt_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: all 16 pixels (%s);t_{ECal} - t_{CDet} (ns);Counts", runNumber, controlBar, selectionLabel), nECalDT, ecalDTMinNs, ecalDTMaxNs);
  for (TH1D* hist : {hBarDt, hControlBarDt}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

  delete hBarLe;
  delete hControlBarLe;
  delete hBarLeVsToT;
  delete hControlBarLeVsToT;
  hBarLe = new TH1D(TString::Format("hBarLe_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: LE (%s);t_{LE} (ns);Counts", runNumber, bar, selectionLabel), nLe, leMinNs, leMaxNs);
  hControlBarLe = new TH1D(TString::Format("hControlBarLe_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: LE (%s);t_{LE} (ns);Counts", runNumber, controlBar, selectionLabel), nLe, leMinNs, leMaxNs);
  hBarLeVsToT = new TH2D(TString::Format("hBarLeVsToT_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", runNumber, bar, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
  hControlBarLeVsToT = new TH2D(TString::Format("hControlBarLeVsToT_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", runNumber, controlBar, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
  TH1* barPulsePlots[] = {hBarLe, hControlBarLe, hBarLeVsToT, hControlBarLeVsToT};
  for (TH1* hist : barPulsePlots) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

  // Before means after energy, explicit LE/ToT, and replay eligibility/spatial cuts.
  // dx is the stored aligned CDet minus projected ECal x, minus the DB offset.
  // Plot ranges do not reject pulses; underflow/overflow counts are retained.
  delete hBarDtVsDxBefore;
  delete hControlBarDtVsDxBefore;
  delete hBarDtVsDxAfter;
  delete hControlBarDtVsDxAfter;
  delete hDetectorDtVsDxBefore;
  delete hDetectorDtVsDxAfter;
  hBarDtVsDxAfter = nullptr;
  hControlBarDtVsDxAfter = nullptr;
  hDetectorDtVsDxAfter = nullptr;
  const char* axes = ";#Deltax_{CDet-ECal} (m);t_{ECal} - t_{CDet,corr} (ns);Pulses";
  hBarDtVsDxBefore = new TH2D(TString::Format("hBarDtVsDxBefore_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: before ellipse%s", runNumber, bar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  hControlBarDtVsDxBefore = new TH2D(TString::Format("hControlBarDtVsDxBefore_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: before ellipse%s", runNumber, controlBar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  hDetectorDtVsDxBefore = new TH2D(TString::Format("hDetectorDtVsDxBefore_run%d", runNumber), TString::Format("Run %d, all CDet pixels: before ellipse%s", runNumber, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  if (useEllipse) {
    hBarDtVsDxAfter = new TH2D(TString::Format("hBarDtVsDxAfter_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: after ellipse%s", runNumber, bar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
    hControlBarDtVsDxAfter = new TH2D(TString::Format("hControlBarDtVsDxAfter_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: after ellipse%s", runNumber, controlBar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
    hDetectorDtVsDxAfter = new TH2D(TString::Format("hDetectorDtVsDxAfter_run%d", runNumber), TString::Format("Run %d, all CDet pixels: after ellipse%s", runNumber, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  }
  for (TH2D* hist : {hBarDtVsDxBefore, hControlBarDtVsDxBefore, hBarDtVsDxAfter, hControlBarDtVsDxAfter, hDetectorDtVsDxBefore, hDetectorDtVsDxAfter}) {
    if (!hist) continue;
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
  std::cout << "  Pre-offset LE cut: " << leCutMinNs << " <= LE <= " << leCutMaxNs << " ns.\n";
  if (applyToTCut)
    std::cout << "  ToT cut: " << totCutMinNs << " <= ToT <= " << totCutMaxNs << " ns.\n";
  else
    std::cout << "  ToT window disabled; finite positive pulse ToT still required.\n";
  std::cout << "  Replay broad_quality_pass is not required; explicit LE/ToT cuts replace it.\n"
            << "  LE plots use 0.01 * pulse.tdc_le, before fitted pixel offsets/time walk.\n";
  if (useEllipse) {
    std::cout << "  Pulse ellipse: center (" << ellipseXCenter << " m, " << ellipseTimeCenter
              << " ns), scales (" << ellipseXScale << " m, " << ellipseTimeScale
              << " ns), radius " << ellipseRadius << ".\n"
              << "  1D spectra/fits use ellipse-selected pulses; timing selection can bias offset fits.\n";
  } else {
    std::cout << "  Pulse ellipse disabled (radius -1); no after-ellipse panel.\n";
  }

  // 4. Loop over all pulses; bar plots and detector-wide plots share the cuts.
  Long64_t processed = 0, energyEvents = 0, selectedPulses = 0, controlPulses = 0, detectorPulses = 0;
  while (processed < limit && reader.Next()) {
    ++processed;
    if (std::isfinite(*ecalEnergy) && *ecalEnergy > energyMin && *ecalEnergy < energyMax) {
      ++energyEvents;
      for (size_t i = 0; i < pixel.GetSize(); ++i) {
        if (!std::isfinite(pixel[i]) || !std::isfinite(pulseDt[i]) || !std::isfinite(pulseDx[i])) continue;
        const double leNs = tdcToNs * pulseLe[i];
        const double totNs = pulseToTNs[i];
        if (!std::isfinite(leNs) || !std::isfinite(totNs) || totNs <= 0.0) continue;
        const int pixelID = int(pixel[i]);
        const int iPixel = pixelID - firstPixel;
        const int iControlPixel = pixelID - firstControlPixel;
        const bool inBar = iPixel >= 0 && iPixel < pixelsPerBar;
        const bool inControlBar = iControlPixel >= 0 && iControlPixel < pixelsPerBar;
        // Do not require broad_quality_pass: it would silently restore 4--30 ns.
        // Spatial/ECal flags are independent of that replay quality decision.
        if (!(calibValid[i] > 0.5 && ecalEligible[i] > 0.5 && spatialPass[i] > 0.5)) continue;
        if (leNs < leCutMinNs || leNs > leCutMaxNs) continue;
        if (applyToTCut && (totNs < totCutMinNs || totNs > totCutMaxNs)) continue;
        hDetectorDtVsDxBefore->Fill(pulseDx[i], pulseDt[i]);
        if (inBar) hBarDtVsDxBefore->Fill(pulseDx[i], pulseDt[i]);
        if (inControlBar) hControlBarDtVsDxBefore->Fill(pulseDx[i], pulseDt[i]);
        if (useEllipse) {
          const double normalizedX = (pulseDx[i] - ellipseXCenter) / ellipseXScale;
          const double normalizedTime = (pulseDt[i] - ellipseTimeCenter) / ellipseTimeScale;
          if (normalizedX * normalizedX + normalizedTime * normalizedTime > ellipseRadius * ellipseRadius)
            continue;
          hDetectorDtVsDxAfter->Fill(pulseDx[i], pulseDt[i]);
          if (inBar) hBarDtVsDxAfter->Fill(pulseDx[i], pulseDt[i]);
          if (inControlBar) hControlBarDtVsDxAfter->Fill(pulseDx[i], pulseDt[i]);
        }
        ++detectorPulses;
        if (inBar) {
          ++selectedPulses;
          hPixelDt[iPixel]->Fill(pulseDt[i]);
          hPixelLe[iPixel]->Fill(leNs);
          hPixelLeVsToT[iPixel]->Fill(totNs, leNs);
        }
        if (inControlBar) {
          ++controlPulses;
          hControlPixelDt[iControlPixel]->Fill(pulseDt[i]);
          hControlPixelLe[iControlPixel]->Fill(leNs);
          hControlPixelLeVsToT[iControlPixel]->Fill(totNs, leNs);
        }
      }
    }

    if (processed % 1000 == 0)
      std::cout << "[CDet pulse pixel timing] " << processed << '/' << limit
                << " events; ECal energy passed: " << energyEvents
                << "; study/control/detector pulses: " << selectedPulses << '/' << controlPulses << '/' << detectorPulses << '\n';
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
    hBarLe->Add(hPixelLe[i]);
    hControlBarLe->Add(hControlPixelLe[i]);
    hBarLeVsToT->Add(hPixelLeVsToT[i]);
    hControlBarLeVsToT->Add(hControlPixelLeVsToT[i]);
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

  // 7. Before/after panels use identical axes and color scales within each bar.
  // Reusing canvas names replaces any previous after panels when rerunning off.
  TCanvas *cBarDtVsDx = new TCanvas(TString::Format("cBarDtVsDx_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: dt vs dx", runNumber, bar), useEllipse ? 1400 : 800, 650);
  if (useEllipse) cBarDtVsDx->Divide(2, 1);
  cBarDtVsDx->cd(useEllipse ? 1 : 0);
  hBarDtVsDxBefore->SetMinimum(0.0);
  hBarDtVsDxBefore->SetMaximum(std::max(1.0, hBarDtVsDxBefore->GetMaximum()));
  DrawDtVsDx(hBarDtVsDxBefore, useEllipse, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale, ellipseRadius);
  if (useEllipse) {
    cBarDtVsDx->cd(2);
    hBarDtVsDxAfter->SetMinimum(0.0);
    hBarDtVsDxAfter->SetMaximum(hBarDtVsDxBefore->GetMaximum());
    DrawDtVsDx(hBarDtVsDxAfter, true, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale, ellipseRadius);
  }
  cBarDtVsDx->Modified();
  cBarDtVsDx->Update();

  TCanvas *cControlBarDtVsDx = new TCanvas(TString::Format("cControlBarDtVsDx_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: dt vs dx", runNumber, controlBar), useEllipse ? 1400 : 800, 650);
  if (useEllipse) cControlBarDtVsDx->Divide(2, 1);
  cControlBarDtVsDx->cd(useEllipse ? 1 : 0);
  hControlBarDtVsDxBefore->SetMinimum(0.0);
  hControlBarDtVsDxBefore->SetMaximum(std::max(1.0, hControlBarDtVsDxBefore->GetMaximum()));
  DrawDtVsDx(hControlBarDtVsDxBefore, useEllipse, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale, ellipseRadius);
  if (useEllipse) {
    cControlBarDtVsDx->cd(2);
    hControlBarDtVsDxAfter->SetMinimum(0.0);
    hControlBarDtVsDxAfter->SetMaximum(hControlBarDtVsDxBefore->GetMaximum());
    DrawDtVsDx(hControlBarDtVsDxAfter, true, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale, ellipseRadius);
  }
  cControlBarDtVsDx->Modified();
  cControlBarDtVsDx->Update();

  TCanvas *cDetectorDtVsDx = new TCanvas(TString::Format("cDetectorDtVsDx_run%d", runNumber), TString::Format("CDet run %d, entire detector: dt vs dx", runNumber), useEllipse ? 1400 : 800, 650);
  if (useEllipse) cDetectorDtVsDx->Divide(2, 1);
  cDetectorDtVsDx->cd(useEllipse ? 1 : 0);
  hDetectorDtVsDxBefore->SetMinimum(0.0);
  hDetectorDtVsDxBefore->SetMaximum(std::max(1.0, hDetectorDtVsDxBefore->GetMaximum()));
  DrawDtVsDx(hDetectorDtVsDxBefore, useEllipse, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale, ellipseRadius);
  if (useEllipse) {
    cDetectorDtVsDx->cd(2);
    hDetectorDtVsDxAfter->SetMinimum(0.0);
    hDetectorDtVsDxAfter->SetMaximum(hDetectorDtVsDxBefore->GetMaximum());
    DrawDtVsDx(hDetectorDtVsDxAfter, true, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale, ellipseRadius);
  }
  cDetectorDtVsDx->Modified();
  cDetectorDtVsDx->Update();

  // 8. Pre-offset LE and LE-versus-ToT for the exact pulses in the 1D dt plots.
  TCanvas *cPixelLe = new TCanvas(TString::Format("cPixelLe_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: LE before pixel offsets", runNumber, bar), 1400, 1000);
  cPixelLe->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cPixelLe->cd(i + 1);
    DrawPulseDistribution(hPixelLe[i], false);
  }
  cPixelLe->Modified();
  cPixelLe->Update();

  TCanvas *cControlPixelLe = new TCanvas(TString::Format("cControlPixelLe_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: LE before pixel offsets", runNumber, controlBar), 1400, 1000);
  cControlPixelLe->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cControlPixelLe->cd(i + 1);
    DrawPulseDistribution(hControlPixelLe[i], false);
  }
  cControlPixelLe->Modified();
  cControlPixelLe->Update();

  TCanvas *cPixelLeVsToT = new TCanvas(TString::Format("cPixelLeVsToT_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: LE vs ToT", runNumber, bar), 1400, 1000);
  cPixelLeVsToT->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cPixelLeVsToT->cd(i + 1);
    DrawPulseDistribution(hPixelLeVsToT[i], true);
  }
  cPixelLeVsToT->Modified();
  cPixelLeVsToT->Update();

  TCanvas *cControlPixelLeVsToT = new TCanvas(TString::Format("cControlPixelLeVsToT_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: LE vs ToT", runNumber, controlBar), 1400, 1000);
  cControlPixelLeVsToT->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cControlPixelLeVsToT->cd(i + 1);
    DrawPulseDistribution(hControlPixelLeVsToT[i], true);
  }
  cControlPixelLeVsToT->Modified();
  cControlPixelLeVsToT->Update();

  TCanvas *cBarLeToT = new TCanvas(TString::Format("cBarLeToT_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: summed LE and LE vs ToT", runNumber, bar), 1400, 650);
  cBarLeToT->Divide(2, 1);
  cBarLeToT->cd(1);
  DrawPulseDistribution(hBarLe, false);
  cBarLeToT->cd(2);
  DrawPulseDistribution(hBarLeVsToT, true);
  cBarLeToT->Modified();
  cBarLeToT->Update();

  TCanvas *cControlBarLeToT = new TCanvas(TString::Format("cControlBarLeToT_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: summed LE and LE vs ToT", runNumber, controlBar), 1400, 650);
  cControlBarLeToT->Divide(2, 1);
  cControlBarLeToT->cd(1);
  DrawPulseDistribution(hControlBarLe, false);
  cControlBarLeToT->cd(2);
  DrawPulseDistribution(hControlBarLeVsToT, true);
  cControlBarLeToT->Modified();
  cControlBarLeToT->Update();

  std::cout << "[CDet pulse pixel timing] Read " << processed << " events; "
            << energyEvents << " passed the ECal energy cut; " << selectedPulses
            << " pulses in bar " << bar << " and " << controlPulses
            << " in control bar " << controlBar << "; " << detectorPulses
            << " across the detector. Drew thirteen canvases.\n";
  if (useEllipse)
    std::cout << "  Before -> after ellipse: study " << hBarDtVsDxBefore->GetEntries()
              << " -> " << selectedPulses << "; control " << hControlBarDtVsDxBefore->GetEntries()
              << " -> " << controlPulses << "; detector " << hDetectorDtVsDxBefore->GetEntries()
              << " -> " << detectorPulses << " pulses.\n";
}

// Retain the run-number shorthand using the separate pulse-timing config.
void Plot_CDet_PulsePixelTiming(
    int runNumber, int bar = -1, const char *inputDirectory = nullptr,
    Long64_t eventsOverride = -2, double ecalDTBinWidthNs = 1.0,
    double ecalDTMinNs = -150.0, double ecalDTMaxNs = 150.0,
    double ecalDXBinWidthM = 0.005, double ecalDXMinM = -0.10, double ecalDXMaxM = 0.10,
    double dtVsDxMinNs = -150.0, double dtVsDxMaxNs = 150.0,
    double leBinWidthNs = 1.0, double leMinNs = -20.0, double leMaxNs = 100.0,
    double totBinWidthNs = 1.0, double totMinNs = 0.0, double totMaxNs = 100.0)
{
  const TString configFile = TString::Format("CDet_run%d_pulse_pixel_timing.conf", runNumber);
  Plot_CDet_PulsePixelTiming(configFile.Data(), bar, inputDirectory, eventsOverride,
                           ecalDTBinWidthNs, ecalDTMinNs, ecalDTMaxNs,
                           ecalDXBinWidthM, ecalDXMinM, ecalDXMaxM, dtVsDxMinNs, dtVsDxMaxNs,
                           leBinWidthNs, leMinNs, leMaxNs, totBinWidthNs, totMinNs, totMaxNs);
}
