// Standalone pulse timing comparison of an investigation bar and a control bar.
// From scripts/cdet, with OUT_DIR set to the replay ROOT-file directory:
//   root -l
//   .L Plot_CDet_PulsePixelTiming.C+
//   Plot_CDet_PulsePixelTiming("CDet_run6077_pulse_pixel_timing.conf");
// Override the investigation bar, input directory, or event limit if desired:
//   Plot_CDet_PulsePixelTiming("CDet_run6077_pulse_pixel_timing.conf", 19, "/path/to/Rootfiles", 10000);
// Numeric shorthand reads CDet_run<run>_pulse_pixel_timing.conf:
//   Plot_CDet_PulsePixelTiming(6077, 19);
// All seventeen canvases save as PDFs under PulsePixelTiming/<tag>/, relative
// to the working directory. Set pulse_pixel_timing.output_tag in the config,
// or override it with the convenient tagged call (remaining arguments optional):
//   Plot_CDet_PulsePixelTiming(6077, "bar19_open_tot", 19);
//   Plot_CDet_PulsePixelTiming("CDet_run6077_pulse_pixel_timing.conf", "bar19_open_tot", 19, "/path/to/Rootfiles", 10000);
// An empty tag defaults to run<run>_bar<study>_control<control>. Reusing a tag
// replaces matching PDFs. Tags use letters, numbers, underscores, and hyphens.
// display.ecal_cdet_dt_min/max set the dt axes in ns for 1D and dt-versus-dx
// plots. Explicit axis arguments override the config; NAN uses the defaults.
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
// A separate event-occupancy comparison counts each pixel once per event,
// before and after an x-only match, without LE/ToT windows, y cuts, or an ellipse.
// Projection opportunities use ECal energy/time from this config, independently
// of CDet activity. DB geometry includes nonresponding instrumented paddles.
// Response fractions in the broad spatial window are not intrinsic efficiencies.
// Manual polygons are not applied. No calibration files or offsets are written.

#include <TCanvas.h>
#include <TBox.h>
#include <TChain.h>
#include <TEnv.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TList.h>
#include <TPad.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

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
TH1D* hPixelOccupancyRaw = nullptr;
TH1D* hPixelOccupancyECal = nullptr;
TH1D* hPixelOccupancyXMatched = nullptr;

TH1D* hECalX = nullptr;
TH1D* hECalXReference = nullptr;
TH1D* hProjectionOpportunities = nullptr;
TH1D* hProjectionComplete = nullptr;
TH1D* hProjectionSelected = nullptr;
TH1D* hProjectionInstrumented = nullptr;
TGraphAsymmErrors* gProjectionCompleteFraction = nullptr;
TGraphAsymmErrors* gProjectionSelectedFraction = nullptr;

namespace CDetPulsePixelTiming {

// Read only projection geometry/alignment keys, inheriting the latest value at
// or before the explicit run timestamp. Do not blindly take the last zpos block:
// the database contains both Kin1 and Kin3 geometries.
bool ReadProjectionDatabase(const char* path, const char* timestamp,
                            std::map<std::string, std::vector<double>>& values)
{
  std::ifstream input(path);
  std::tm requestedDate{};
  std::istringstream dateInput(timestamp);
  dateInput >> std::get_time(&requestedDate, "%Y-%m-%d %H:%M:%S");
  if (!input || dateInput.fail() || std::string(timestamp).size() != 19) {
    std::cerr << "[CDet projection] Need a readable database and timestamp YYYY-MM-DD HH:MM:SS.\n";
    return false;
  }
  const char* keys[] = {"xpos", "ypos", "zpos", "selection.x_scale_l1", "selection.x_scale_l2",
      "selection.x_alignment_l1", "selection.x_alignment_l2", "selection.x_residual_offset",
      "selection.y_residual_offset", "selection.x_residual_max", "selection.y_residual_max"};
  for (const char* key : keys) values[key] = {};
  std::map<std::string, std::string> selectedDates;
  std::string line, currentDate, readingKey;
  auto trim = [](const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    return first == std::string::npos ? std::string() : s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
  };
  while (std::getline(input, line)) {
    std::string content = trim(line.substr(0, line.find('#')));
    if (content.empty()) continue;
    const auto bracket = content.find('[');
    if (bracket != std::string::npos && content.find(']', bracket) != std::string::npos) {
      currentDate = trim(content.substr(bracket + 1, content.find(']', bracket) - bracket - 1));
      readingKey.clear();
      continue;
    }
    const auto equals = content.find('=');
    if (equals != std::string::npos) {
      readingKey.clear();
      const std::string prefix = "earm.cdet.";
      const std::string fullKey = trim(content.substr(0, equals));
      if (fullKey.compare(0, prefix.size(), prefix) != 0) continue;
      const std::string key = fullKey.substr(prefix.size());
      if (!values.count(key) || currentDate > timestamp || currentDate < selectedDates[key]) continue;
      readingKey = key;
      selectedDates[key] = currentDate;
      values[key].clear();
      content = content.substr(equals + 1);
    }
    if (readingKey.empty()) continue;
    std::istringstream numbers(content);
    double number;
    while (numbers >> number) values[readingKey].push_back(number);
  }
  for (const char* key : keys) {
    const bool positions = std::string(key).find("pos") != std::string::npos;
    if ((positions && values[key].size() < 2688) || (!positions && values[key].size() != 1)) {
      std::cerr << "[CDet projection] Missing or malformed DB key at " << timestamp << ": " << key << '\n';
      return false;
    }
  }
  return true;
}

void DrawProjectionResponse(int firstPixel, int lastPixel, const TString& title, bool fractions)
{
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  const TString frameName = TString::Format("projectionFrame_%d_%d", firstPixel, int(fractions));
  TH1D frame(frameName, title, lastPixel - firstPixel + 1, firstPixel - 0.5, lastPixel + 0.5);
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitle("Logical pixel number");
  frame.GetYaxis()->SetTitle(fractions ? "Conditional response fraction" : "Events");
  double maximum = 1.0;
  for (int pixel = firstPixel; pixel <= lastPixel; ++pixel)
    maximum = std::max(maximum, hProjectionOpportunities->GetBinContent(pixel + 1));
  frame.SetMinimum(0.0);
  frame.SetMaximum(fractions ? 1.25 : 1.4 * maximum);
  frame.DrawCopy("AXIS");
  for (int pixel = firstPixel; pixel <= lastPixel; ++pixel) {
    if (hProjectionInstrumented->GetBinContent(pixel + 1) > 0.5) continue;
    TBox* unused = new TBox(pixel - 0.5, 0.0, pixel + 0.5, frame.GetMaximum());
    unused->SetFillColor(kGray);
    unused->SetFillStyle(3004);
    unused->Draw();
  }
  TLegend* legend = new TLegend(0.15, 0.73, 0.89, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.031);
  if (fractions) {
    auto* complete = gProjectionCompleteFraction->DrawClone("P SAME");
    auto* selected = gProjectionSelectedFraction->DrawClone("P SAME");
    legend->AddEntry(complete, "Complete pulse / opportunity", "pe");
    legend->AddEntry(selected, "Final-selected pulse / opportunity", "pe");
    legend->AddEntry((TObject*)nullptr, "No opportunity: no point; shaded: unused slot", "");
  } else {
    for (TH1D* hist : {hProjectionOpportunities, hProjectionComplete, hProjectionSelected}) {
      hist->GetXaxis()->SetRange(firstPixel + 1, lastPixel + 1);
      auto* drawn = hist->DrawCopy("HIST SAME");
      legend->AddEntry(drawn, hist == hProjectionOpportunities ? "ECal projection opportunities" :
          hist == hProjectionComplete ? "With a complete pulse" : "With a final-selected pulse", "l");
      hist->GetXaxis()->SetRange(0, 0);
    }
  }
  legend->Draw();
  gPad->RedrawAxis();
  TLatex panelTitle;
  panelTitle.SetNDC();
  panelTitle.SetTextAlign(23);
  panelTitle.SetTextSize(0.042);
  panelTitle.DrawLatex(0.5, 0.98, title);
}

void DrawPixelOccupancy(int firstPixel, int lastPixel, const TString& title,
                        bool useECalBaseline, double dxMaxM)
{
  TH1D* baseline = useECalBaseline ? hPixelOccupancyECal : hPixelOccupancyRaw;
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.13);
  baseline->GetXaxis()->SetRange(firstPixel + 1, lastPixel + 1);
  baseline->SetTitle(title);
  baseline->SetMinimum(0.0);
  baseline->SetMaximum(1.4 * std::max(1.0, baseline->GetMaximum()));
  TH1* before = baseline->DrawCopy("HIST");
  hPixelOccupancyXMatched->GetXaxis()->SetRange(firstPixel + 1, lastPixel + 1);
  TH1* after = hPixelOccupancyXMatched->DrawCopy("HIST SAME");
  TLegend* legend = new TLegend(0.14, 0.76, 0.88, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.033);
  legend->AddEntry(before, useECalBaseline ? "ECal selected, before x match" : "All stored complete pulses", "l");
  legend->AddEntry(after, TString::Format("ECal selected + |#Deltax| #leq %.3f m", dxMaxM), "l");
  legend->Draw();
  baseline->GetXaxis()->SetRange(0, 0);
  baseline->SetMaximum(-1111);
  hPixelOccupancyXMatched->GetXaxis()->SetRange(0, 0);
}

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
    double ecalDTMinNs = NAN, double ecalDTMaxNs = NAN,
    double ecalDXBinWidthM = 0.005, double ecalDXMinM = -0.10, double ecalDXMaxM = 0.10,
    double dtVsDxMinNs = NAN, double dtVsDxMaxNs = NAN,
    double leBinWidthNs = 1.0, double leMinNs = -20.0, double leMaxNs = 100.0,
    double totBinWidthNs = 1.0, double totMinNs = 0.0, double totMaxNs = 100.0,
    const char *outputTag = nullptr)
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
  // Resolve optional axis overrides before validating ranges and making bins.
  if (std::isnan(ecalDTMinNs)) ecalDTMinNs = config.GetValue("display.ecal_cdet_dt_min", -150.0);
  if (std::isnan(ecalDTMaxNs)) ecalDTMaxNs = config.GetValue("display.ecal_cdet_dt_max", 150.0);
  if (std::isnan(dtVsDxMinNs)) dtVsDxMinNs = ecalDTMinNs;
  if (std::isnan(dtVsDxMaxNs)) dtVsDxMaxNs = ecalDTMaxNs;
  const double energyMin = config.GetValue("analysis.ecal_energy_min", NAN);
  const double energyMax = config.GetValue("analysis.ecal_energy_max", NAN);
  const double ecalTimeMinNs = config.GetValue("analysis.ecal_time_min", NAN);
  const double ecalTimeMaxNs = config.GetValue("analysis.ecal_time_max", NAN);
  const double leCutMinNs = config.GetValue("analysis.le_min", 0.02);
  const double leCutMaxNs = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_pixel_timing.apply_tot_cut", 0) != 0;
  const bool occupancyUseECalBaseline = config.GetValue("pulse_pixel_timing.occupancy_use_ecal_baseline", 1) != 0;
  const double occupancyDxMaxM = config.GetValue("pulse_pixel_timing.occupancy_dx_max_m", 0.08);
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
  const int nPixels = nBars * pixelsPerBar;
  const double projectionECalZM = config.GetValue("projection.ecal_z_m", NAN);
  const TString geometryTimestamp = config.GetValue("projection.database_timestamp", "");
  TString geometryFile = config.GetValue("projection.database_file", "../../DB/db_earm.cdet.dat");
  if (!gSystem->IsAbsoluteFileName(geometryFile))
    geometryFile = TString::Format("%s/%s", gSystem->DirName(configFile), geometryFile.Data());
  if (!std::isfinite(ecalTimeMinNs) || !std::isfinite(ecalTimeMaxNs) || ecalTimeMinNs >= ecalTimeMaxNs ||
      !std::isfinite(projectionECalZM) || projectionECalZM <= 0.0) {
    std::cerr << "[CDet projection] Set valid analysis.ecal_time_min/max and projection.ecal_z_m.\n";
    return;
  }
  std::map<std::string, std::vector<double>> geometryDB;
  if (!ReadProjectionDatabase(geometryFile, geometryTimestamp, geometryDB)) return;
  for (const auto& item : geometryDB) {
    for (double value : item.second) {
      if (!std::isfinite(value)) {
        std::cerr << "[CDet projection] Nonfinite database value: " << item.first << '\n';
        return;
      }
    }
  }
  const double projectionDxMaxM = config.GetValue("projection.dx_max_m", geometryDB["selection.x_residual_max"][0]);
  const double projectionDyMaxM = config.GetValue("projection.dy_max_m", geometryDB["selection.y_residual_max"][0]);
  const double projectionXOffsetM = geometryDB["selection.x_residual_offset"][0];
  const double projectionYOffsetM = geometryDB["selection.y_residual_offset"][0];
  if (!std::isfinite(projectionDxMaxM) || projectionDxMaxM <= 0.0 ||
      !std::isfinite(projectionDyMaxM) || projectionDyMaxM <= 0.0) {
    std::cerr << "[CDet projection] Projection x/y half-windows must be positive and finite.\n";
    return;
  }
  std::array<double, 2688> geometryX{}, geometryY{}, geometryZ{};
  std::array<bool, 2688> instrumented{};
  std::array<int, 2> instrumentedPerLayer{};
  for (int id = 0; id < nPixels; ++id) {
    const double x = geometryDB["xpos"][id], y = geometryDB["ypos"][id], z = geometryDB["zpos"][id];
    // Sentinel channels are unused. Low-occupancy but physical channels stay in.
    if (std::fabs(x) >= 900.0 || std::fabs(y) >= 900.0 || std::fabs(z) >= 900.0) continue;
    const int layer = id / 1344;
    const double xScale = geometryDB[layer == 0 ? "selection.x_scale_l1" : "selection.x_scale_l2"][0];
    const double xAlignment = geometryDB[layer == 0 ? "selection.x_alignment_l1" : "selection.x_alignment_l2"][0];
    if (xScale <= 0.0 || z <= 0.0 || z >= projectionECalZM) {
      std::cerr << "[CDet projection] Invalid alignment or layer z; check DB timestamp and ECal distance.\n";
      return;
    }
    geometryX[id] = xScale * x - xAlignment;
    geometryY[id] = y;
    geometryZ[id] = z;
    instrumented[id] = true;
    ++instrumentedPerLayer[layer];
  }
  if (instrumentedPerLayer[0] == 0 || instrumentedPerLayer[1] == 0) {
    std::cerr << "[CDet projection] No instrumented geometry in one or both layers.\n";
    return;
  }
  if (!std::isfinite(occupancyDxMaxM) || occupancyDxMaxM <= 0.0) {
    std::cerr << "[CDet pulse pixel timing] occupancy_dx_max_m must be positive and finite.\n";
    return;
  }
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
  const double ecalXMinM = -1.5, ecalXMaxM = 1.5;
  const double ecalXBinWidthM = 0.0425; // Target width; 71 equal bins give about 4.23 cm.
  const int nECalX = int(std::ceil((ecalXMaxM - ecalXMinM) / ecalXBinWidthM));
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

  // Prepare a tagged output directory before processing the events.
  TString tag = outputTag ? outputTag : config.GetValue("pulse_pixel_timing.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d_bar%d_control%d", runNumber, bar, controlBar);
  const TString allowedTagCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowedTagCharacters.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet pulse pixel timing] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/PulsePixelTiming/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet pulse pixel timing] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }
  std::cout << "[CDet pulse pixel timing] PDF output: " << outputDirectory << '\n';

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
  TTreeReaderValue<double> ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<double> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<double> ecalNclus(reader, "earm.ecal.nclus");
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

  delete hECalX;
  delete hECalXReference;
  hECalX = new TH1D(TString::Format("hECalX_run%d", runNumber), TString::Format("Run %d: event with a final-selected CDet pulse;ECal x (m);Events", runNumber), nECalX, ecalXMinM, ecalXMaxM);
  hECalXReference = new TH1D(TString::Format("hECalXReference_run%d", runNumber), TString::Format("Run %d: ECal energy/time reference, no CDet requirement;ECal x (m);Events", runNumber), nECalX, ecalXMinM, ecalXMaxM);
  for (TH1D* hist : {hECalX, hECalXReference}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

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

  delete hPixelOccupancyRaw;
  delete hPixelOccupancyECal;
  delete hPixelOccupancyXMatched;
  hPixelOccupancyRaw = new TH1D(TString::Format("hPixelOccupancyRaw_run%d", runNumber), "All stored complete pulses;Logical pixel number;Events with #geq1 pulse", nPixels, -0.5, nPixels - 0.5);
  hPixelOccupancyECal = new TH1D(TString::Format("hPixelOccupancyECal_run%d", runNumber), "ECal-selected pulses before x match;Logical pixel number;Events with #geq1 pulse", nPixels, -0.5, nPixels - 0.5);
  hPixelOccupancyXMatched = new TH1D(TString::Format("hPixelOccupancyXMatched_run%d", runNumber), "ECal-selected pulses with x match;Logical pixel number;Events with #geq1 pulse", nPixels, -0.5, nPixels - 0.5);
  for (TH1D* hist : {hPixelOccupancyRaw, hPixelOccupancyECal, hPixelOccupancyXMatched}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue + 1);
  }
  hPixelOccupancyXMatched->SetLineColor(kRed + 1);

  delete hProjectionOpportunities;
  delete hProjectionComplete;
  delete hProjectionSelected;
  delete hProjectionInstrumented;
  hProjectionOpportunities = new TH1D(TString::Format("hProjectionOpportunities_run%d", runNumber), "ECal projection opportunities;Logical pixel number;Events", nPixels, -0.5, nPixels - 0.5);
  hProjectionComplete = new TH1D(TString::Format("hProjectionComplete_run%d", runNumber), "Projection opportunity with a complete pulse;Logical pixel number;Events", nPixels, -0.5, nPixels - 0.5);
  hProjectionSelected = new TH1D(TString::Format("hProjectionSelected_run%d", runNumber), "Projection opportunity with a final-selected pulse;Logical pixel number;Events", nPixels, -0.5, nPixels - 0.5);
  hProjectionInstrumented = new TH1D(TString::Format("hProjectionInstrumented_run%d", runNumber), "DB instrumented pixel mask;Logical pixel number;Instrumented", nPixels, -0.5, nPixels - 0.5);
  for (TH1D* hist : {hProjectionOpportunities, hProjectionComplete, hProjectionSelected, hProjectionInstrumented}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
    hist->SetLineWidth(2);
  }
  hProjectionOpportunities->SetLineColor(kBlack);
  hProjectionComplete->SetLineColor(kBlue + 1);
  hProjectionSelected->SetLineColor(kRed + 1);
  for (int id = 0; id < nPixels; ++id) hProjectionInstrumented->SetBinContent(id + 1, instrumented[id]);

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

  std::cout << "  Occupancy baseline: " << (occupancyUseECalBaseline ? "ECal selected" : "all stored complete pulses")
            << "; x match |pulse.ecal_x_residual| <= " << occupancyDxMaxM << " m.\n"
            << "  Occupancy counts each event once per pixel. ECal selection uses energy,\n"
            << "  replay eligibility, and finite stored dx; no LE/ToT window, y match, or ellipse.\n";

  // 4. Define individual cuts, combine them into named selections, then fill.
  // Occupancy and timing spectra deliberately use different selections.
  Long64_t processed = 0, energyEvents = 0, selectedPulses = 0, controlPulses = 0, detectorPulses = 0;
  Long64_t referenceEvents = 0, referenceWithoutPulses = 0, referenceWithoutSelectedPulses = 0;
  std::cout << "[CDet projection] DB: " << geometryFile << "; timestamp: " << geometryTimestamp
            << "; instrumented L1/L2: " << instrumentedPerLayer[0] << '/' << instrumentedPerLayer[1] << '\n'
            << "  ECal-only reference: " << energyMin << " < E < " << energyMax << " GeV, "
            << ecalTimeMinNs << " < ADC time < " << ecalTimeMaxNs << " ns; valid cluster and finite x/y.\n"
            << "  Target-origin rays; ECal z = " << projectionECalZM << " m; |dx| <= "
            << projectionDxMaxM << " m, |dy| <= " << projectionDyMaxM << " m after DB alignment.\n"
            << "  Complete response: any stored pulse. Selected response: existing good_final cuts.\n";
  while (processed < limit && reader.Next()) {
    ++processed;
    bool eventHasGoodPulse = false;
    // Event-level cut: evaluated once, before looking at individual pulses.
    const bool good_ECal_energy = std::isfinite(*ecalEnergy) &&
        *ecalEnergy > energyMin && *ecalEnergy < energyMax;
    const bool good_ECal_time = std::isfinite(*ecalTime) &&
        *ecalTime > ecalTimeMinNs && *ecalTime < ecalTimeMaxNs;
    const bool good_ECal_geometry = std::isfinite(*ecalNclus) && *ecalNclus > 0.0 &&
        std::isfinite(*ecalX) && std::isfinite(*ecalY);
    // This denominator does not use ANY pulse branch, replay CDet flag, or pair.
    const bool good_ECal_reference = good_ECal_energy && good_ECal_time && good_ECal_geometry;
    if (good_ECal_energy) ++energyEvents;

    std::array<bool, 2688> seenRaw{}, seenECal{}, seenXMatched{};
    std::array<bool, 2688> seenSelected{};
    for (size_t i = 0; i < pixel.GetSize(); ++i) {
      const bool good_pixel_value = std::isfinite(pixel[i]);
      const int pixelID = good_pixel_value ? int(pixel[i]) : -1;
      const int iPixel = pixelID - firstPixel;
      const int iControlPixel = pixelID - firstControlPixel;
      const bool inBar = iPixel >= 0 && iPixel < pixelsPerBar;
      const bool inControlBar = iControlPixel >= 0 && iControlPixel < pixelsPerBar;
      const double leNs = tdcToNs * pulseLe[i];
      const double totNs = pulseToTNs[i];

      // Individual replay decisions. These read stored flags, not conf windows:
      // - ECal eligibility includes a valid cluster, ECal time (LH2: -10..10 ns),
      //   ECal x/y acceptance, and valid projection geometry from the replay DB.
      // - Spatial pass includes BOTH CDet-ECal x and y residual cuts.
      // - Calibration validity does not imply a measured nonzero pixel offset.
      const bool good_ECal_replay = ecalEligible[i] > 0.5;
      const bool good_spatial_replay = spatialPass[i] > 0.5;
      const bool good_calibration = calibValid[i] > 0.5;

      // Individual cuts evaluated here. No replay broad_quality_pass gate:
      // it would reimpose the replay's LE/ToT windows instead of these conf cuts.
      const bool good_dx_value = std::isfinite(pulseDx[i]);
      const bool good_timing_values = good_pixel_value && good_dx_value &&
          std::isfinite(pulseDt[i]) && std::isfinite(leNs) &&
          std::isfinite(totNs) && totNs > 0.0;
      const bool good_le_time = leNs >= leCutMinNs && leNs <= leCutMaxNs;
      const bool good_tot_window = !applyToTCut ||
          (totNs >= totCutMinNs && totNs <= totCutMaxNs);
      const bool good_occupancy_x = good_dx_value && std::fabs(pulseDx[i]) <= occupancyDxMaxM;

      bool good_ellipse = true; // Disabled ellipse accepts every otherwise-selected pulse.
      if (useEllipse) {
        const double normalizedX = (pulseDx[i] - ellipseXCenter) / ellipseXScale;
        const double normalizedTime = (pulseDt[i] - ellipseTimeCenter) / ellipseTimeScale;
        good_ellipse = normalizedX * normalizedX + normalizedTime * normalizedTime <= ellipseRadius * ellipseRadius;
      }

      // Combined selections: this is the cut list for each histogram group.
      // Finite stored dx inherits replay's ability to form a calibrated pulse/
      // ECal projection; occupancy makes no cut on the corrected CDet time.
      const bool good_raw_occupancy = good_pixel_value && pixelID >= 0 && pixelID < nPixels;
      const bool good_ECal_occupancy = good_raw_occupancy && good_ECal_energy &&
          good_ECal_replay && good_dx_value;
      const bool good_x_occupancy = good_ECal_occupancy && good_occupancy_x;
      const bool good_before_ellipse = good_ECal_energy && good_timing_values &&
          good_calibration && good_ECal_replay && good_spatial_replay &&
          good_le_time && good_tot_window;
      const bool good_final = good_before_ellipse && good_ellipse;

      // Occupancy: each event contributes at most once per pixel per selection.
      if (good_raw_occupancy && !seenRaw[pixelID]) {
        hPixelOccupancyRaw->Fill(pixelID);
        seenRaw[pixelID] = true;
      }
      if (good_ECal_occupancy && !seenECal[pixelID]) {
        hPixelOccupancyECal->Fill(pixelID);
        seenECal[pixelID] = true;
      }
      if (good_x_occupancy && !seenXMatched[pixelID]) {
        hPixelOccupancyXMatched->Fill(pixelID);
        seenXMatched[pixelID] = true;
      }

      // Timing diagnostic before ellipse: all other timing-spectrum cuts applied.
      if (good_before_ellipse) {
        hDetectorDtVsDxBefore->Fill(pulseDx[i], pulseDt[i]);
        if (inBar) hBarDtVsDxBefore->Fill(pulseDx[i], pulseDt[i]);
        if (inControlBar) hControlBarDtVsDxBefore->Fill(pulseDx[i], pulseDt[i]);
      }

      // Final timing selection: shared by dt, LE, and LE-versus-ToT spectra.
      if (good_final) {
        eventHasGoodPulse = true;
        if (good_raw_occupancy) seenSelected[pixelID] = true;
        if (useEllipse) {
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
    if (eventHasGoodPulse) {
      hECalX->Fill(*ecalX);
    }
    // Count an opportunity even when the pulse arrays are empty. A broad window
    // can admit several candidate paddles; it does not predict they all fire.
    if (good_ECal_reference) {
      ++referenceEvents;
      hECalXReference->Fill(*ecalX);
      if (pixel.GetSize() == 0) ++referenceWithoutPulses;
      if (!eventHasGoodPulse) ++referenceWithoutSelectedPulses;
      for (int id = 0; id < nPixels; ++id) {
        if (!instrumented[id]) continue;
        const double projectedX = *ecalX * geometryZ[id] / projectionECalZM;
        const double projectedY = *ecalY * geometryZ[id] / projectionECalZM;
        const double dx = geometryX[id] - projectedX - projectionXOffsetM;
        const double dy = geometryY[id] - projectedY - projectionYOffsetM;
        const bool good_projection_x = std::fabs(dx) <= projectionDxMaxM;
        const bool good_projection_y = std::fabs(dy) <= projectionDyMaxM;
        const bool good_opportunity = good_projection_x && good_projection_y;
        if (good_opportunity) {
          hProjectionOpportunities->Fill(id);
          if (seenRaw[id]) hProjectionComplete->Fill(id);
          if (seenSelected[id]) hProjectionSelected->Fill(id);
        }
      }
    }
    if (processed % 1000 == 0)
      std::cout << "[CDet pulse pixel timing] " << processed << '/' << limit
                << " events; ECal energy passed: " << energyEvents
                << "; study/control/detector pulses: " << selectedPulses << '/' << controlPulses << '/' << detectorPulses
                << "; ECal-only reference events: " << referenceEvents << '\n';
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

  TCanvas *cPixelOccupancy = new TCanvas(TString::Format("cPixelOccupancy_run%d_bar%d_control%d", runNumber, bar, controlBar), "CDet event occupancy before and after ECal x match", 1600, 1000);
  cPixelOccupancy->Divide(2, 2);
  cPixelOccupancy->cd(1);
  DrawPixelOccupancy(0, 1343, TString::Format("Run %d, Layer 1;Logical pixel number;Events with #geq1 pulse", runNumber), occupancyUseECalBaseline, occupancyDxMaxM);
  cPixelOccupancy->cd(2);
  DrawPixelOccupancy(1344, 2687, TString::Format("Run %d, Layer 2;Logical pixel number;Events with #geq1 pulse", runNumber), occupancyUseECalBaseline, occupancyDxMaxM);
  cPixelOccupancy->cd(3);
  DrawPixelOccupancy(firstPixel, firstPixel + pixelsPerBar - 1, TString::Format("Run %d, study bar %d;Logical pixel number;Events with #geq1 pulse", runNumber, bar), occupancyUseECalBaseline, occupancyDxMaxM);
  cPixelOccupancy->cd(4);
  DrawPixelOccupancy(firstControlPixel, firstControlPixel + pixelsPerBar - 1, TString::Format("Run %d, control bar %d;Logical pixel number;Events with #geq1 pulse", runNumber, controlBar), occupancyUseECalBaseline, occupancyDxMaxM);
  cPixelOccupancy->Modified();
  cPixelOccupancy->Update();

  // Binomial intervals for each conditional fraction; bins with no denominator
  // are omitted from the graphs, rather than displayed as zero response.
  delete gProjectionCompleteFraction;
  delete gProjectionSelectedFraction;
  gProjectionCompleteFraction = new TGraphAsymmErrors;
  gProjectionSelectedFraction = new TGraphAsymmErrors;
  gProjectionCompleteFraction->SetName(TString::Format("gProjectionCompleteFraction_run%d", runNumber));
  gProjectionSelectedFraction->SetName(TString::Format("gProjectionSelectedFraction_run%d", runNumber));
  gProjectionCompleteFraction->Divide(hProjectionComplete, hProjectionOpportunities, "cl=0.683 cp");
  gProjectionSelectedFraction->Divide(hProjectionSelected, hProjectionOpportunities, "cl=0.683 cp");
  for (TGraphAsymmErrors* graph : {gProjectionCompleteFraction, gProjectionSelectedFraction}) {
    const int color = graph == gProjectionCompleteFraction ? kBlue + 1 : kRed + 1;
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(graph == gProjectionCompleteFraction ? 20 : 24);
    graph->SetMarkerSize(0.5);
  }
  TCanvas *cProjectionCounts = new TCanvas(TString::Format("cProjectionCounts_run%d_bar%d_control%d", runNumber, bar, controlBar), "ECal projection opportunities and CDet responses", 1600, 1000);
  TCanvas *cProjectionFractions = new TCanvas(TString::Format("cProjectionFractions_run%d_bar%d_control%d", runNumber, bar, controlBar), "Conditional CDet response fractions", 1600, 1000);
  const int firstProjectionPixel[] = {0, 1344, firstPixel, firstControlPixel};
  const int lastProjectionPixel[] = {1343, 2687, firstPixel + 15, firstControlPixel + 15};
  const TString projectionTitles[] = {TString::Format("Run %d, Layer 1", runNumber), TString::Format("Run %d, Layer 2", runNumber), TString::Format("Run %d, study bar %d", runNumber, bar), TString::Format("Run %d, control bar %d", runNumber, controlBar)};
  for (TCanvas* canvas : {cProjectionCounts, cProjectionFractions}) {
    canvas->Divide(2, 2);
    for (int panel = 0; panel < 4; ++panel) {
      canvas->cd(panel + 1);
      DrawProjectionResponse(firstProjectionPixel[panel], lastProjectionPixel[panel], projectionTitles[panel], canvas == cProjectionFractions);
    }
    canvas->Modified();
    canvas->Update();
  }
  TCanvas *cECalXSelection = new TCanvas(TString::Format("cECalXSelection_run%d", runNumber), "ECal x: independent reference and final CDet selection", 1400, 650);
  cECalXSelection->Divide(2, 1);
  cECalXSelection->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.13);
  hECalXReference->Draw("HIST");
  cECalXSelection->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.13);
  hECalX->Draw("HIST");
  cECalXSelection->Modified();
  cECalXSelection->Update();

  const TString projectionCSV = TString::Format("%s/projection_response_run%d.csv", outputDirectory.Data(), runNumber);
  std::ofstream responseTable(projectionCSV.Data());
  responseTable << "# DB=" << geometryFile << "; timestamp=" << geometryTimestamp << '\n'
                << "# ECal E (GeV): " << energyMin << " < E < " << energyMax
                << "; ADC time (ns): " << ecalTimeMinNs << " < t < " << ecalTimeMaxNs << '\n'
                << "# ECal z (m)=" << projectionECalZM << "; x/y half-windows (m)=" << projectionDxMaxM << '/' << projectionDyMaxM
                << "; origin projection; final-selected response uses existing good_final cuts\n"
                << "# good_final LE (ns)=" << leCutMinNs << ".." << leCutMaxNs
                << "; ToT window enabled=" << applyToTCut << ", bounds (ns)=" << totCutMinNs << ".." << totCutMaxNs
                << "; ellipse radius=" << ellipseRadius << "; replay calibration/ECal/spatial flags required\n"
                << "pixel,bar,layer,instrumented,x_aligned_m,y_m,z_m,opportunities,complete_responses,selected_responses,complete_fraction,selected_fraction\n";
  for (int id = 0; id < nPixels; ++id) {
    const double opportunities = hProjectionOpportunities->GetBinContent(id + 1);
    const double complete = hProjectionComplete->GetBinContent(id + 1);
    const double selected = hProjectionSelected->GetBinContent(id + 1);
    responseTable << id << ',' << id / 16 << ',' << id / 1344 + 1 << ',' << instrumented[id] << ',';
    if (instrumented[id]) responseTable << geometryX[id] << ',' << geometryY[id] << ',' << geometryZ[id];
    else responseTable << ",,";
    responseTable << ',' << opportunities << ',' << complete << ',' << selected << ',';
    if (opportunities > 0.0) responseTable << complete / opportunities << ',' << selected / opportunities;
    else responseTable << ',';
    responseTable << '\n';
  }
  responseTable.close();
  if (!responseTable) std::cerr << "[CDet projection] Could not write complete response table: " << projectionCSV << '\n';

  // Save only canvases from this invocation, preserving the displayed panels.
  TCanvas* outputCanvases[] = {cPixelDt, cBarDt, cControlPixelDt, cControlBarDt,
      cBarDtVsDx, cControlBarDtVsDx, cDetectorDtVsDx, cPixelLe, cControlPixelLe,
      cPixelLeVsToT, cControlPixelLeVsToT, cBarLeToT, cControlBarLeToT, cPixelOccupancy,
      cProjectionCounts, cProjectionFractions, cECalXSelection};
  for (TCanvas* canvas : outputCanvases) {
    canvas->SaveAs(TString::Format("%s/%s.pdf", outputDirectory.Data(), canvas->GetName()));
  }

  std::cout << "[CDet pulse pixel timing] Read " << processed << " events; "
            << energyEvents << " passed the ECal energy cut; " << selectedPulses
            << " pulses in bar " << bar << " and " << controlPulses
            << " in control bar " << controlBar << "; " << detectorPulses
            << " across the detector. Drew seventeen canvases.\n"
            << "  PDF output directory: " << outputDirectory << '\n';
  std::cout << "[CDet projection] " << referenceEvents << " ECal-only reference events; "
            << referenceWithoutPulses << " have no stored CDet pulses; " << referenceWithoutSelectedPulses
            << " have no final-selected CDet pulse.\n";
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
    double ecalDTMinNs = NAN, double ecalDTMaxNs = NAN,
    double ecalDXBinWidthM = 0.005, double ecalDXMinM = -0.10, double ecalDXMaxM = 0.10,
    double dtVsDxMinNs = NAN, double dtVsDxMaxNs = NAN,
    double leBinWidthNs = 1.0, double leMinNs = -20.0, double leMaxNs = 100.0,
    double totBinWidthNs = 1.0, double totMinNs = 0.0, double totMaxNs = 100.0,
    const char *outputTag = nullptr)
{
  const TString configFile = TString::Format("CDet_run%d_pulse_pixel_timing.conf", runNumber);
  Plot_CDet_PulsePixelTiming(configFile.Data(), bar, inputDirectory, eventsOverride,
                           ecalDTBinWidthNs, ecalDTMinNs, ecalDTMaxNs,
                           ecalDXBinWidthM, ecalDXMinM, ecalDXMaxM, dtVsDxMinNs, dtVsDxMaxNs,
                           leBinWidthNs, leMinNs, leMaxNs, totBinWidthNs, totMinNs, totMaxNs,
                           outputTag);
}

// Convenient tagged calls use the configured dt axes and existing call forms.
void Plot_CDet_PulsePixelTiming(
    const char *configFile, const char *outputTag, int bar = -1,
    const char *inputDirectory = nullptr, Long64_t eventsOverride = -2)
{
  Plot_CDet_PulsePixelTiming(configFile, bar, inputDirectory, eventsOverride,
                           1.0, NAN, NAN, 0.005, -0.10, 0.10, NAN, NAN,
                           1.0, -20.0, 100.0, 1.0, 0.0, 100.0, outputTag);
}

void Plot_CDet_PulsePixelTiming(
    int runNumber, const char *outputTag, int bar = -1,
    const char *inputDirectory = nullptr, Long64_t eventsOverride = -2)
{
  const TString configFile = TString::Format("CDet_run%d_pulse_pixel_timing.conf", runNumber);
  Plot_CDet_PulsePixelTiming(configFile.Data(), outputTag, bar, inputDirectory, eventsOverride);
}
