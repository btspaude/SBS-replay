// CDet occupancy and ECal-projected hit opportunities, independent of timing plots.
// .L Plot_CDet_OccupancyResponse.C+
// Plot_CDet_OccupancyResponse(6077, 19, nullptr, 1000);
// Reads CDet_run<run>_occupancy_response.conf; OUT_DIR supplies replay files.
// Four PDFs and a per-pixel CSV go to OccupancyResponse/<tag>/.
// Counts are once per event/pixel. Broad-window fractions are conditional
// responses, not intrinsic efficiencies. No calibration or replay files are changed.
#include <TCanvas.h>
#include <TBox.h>
#include <TChain.h>
#include <TEnv.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
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

// Available at the ROOT prompt after this macro returns.
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

namespace CDetOccupancyResponse {

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

} // namespace CDetOccupancyResponse

void Plot_CDet_OccupancyResponse(
    const char* configFile = "CDet_run6077_occupancy_response.conf", int bar = -1,
    const char* inputDirectory = nullptr, Long64_t eventsOverride = -2,
    double ecalXBinWidthM = 0.0425, double ecalXMinM = -1.5, double ecalXMaxM = 1.5,
    const char* outputTag = nullptr)
{
  using namespace CDetOccupancyResponse;
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0) {
    std::cerr << "[CDet occupancy response] Cannot read the configuration file.\n";
    return;
  }
  if (config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet occupancy response] Expected config.version: 1.\n";
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
  const double ecalTimeMinNs = config.GetValue("analysis.ecal_time_min", NAN);
  const double ecalTimeMaxNs = config.GetValue("analysis.ecal_time_max", NAN);
  const double leCutMinNs = config.GetValue("analysis.le_min", 0.02);
  const double leCutMaxNs = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("occupancy_response.apply_tot_cut", 0) != 0;
  const bool occupancyUseECalBaseline = config.GetValue("occupancy_response.occupancy_use_ecal_baseline", 1) != 0;
  const double occupancyDxMaxM = config.GetValue("occupancy_response.occupancy_dx_max_m", 0.08);
  const double totCutMinNs = config.GetValue("analysis.tot_min", 4.0);
  const double totCutMaxNs = config.GetValue("analysis.tot_max", 30.0);
  if (bar == -1) bar = config.GetValue("occupancy_response.study_bar", 19);
  const int controlBar = config.GetValue("occupancy_response.control_bar", 30);
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
    std::cerr << "[CDet occupancy response] occupancy_dx_max_m must be positive and finite.\n";
    return;
  }
  if (!std::isfinite(energyMin) || !std::isfinite(energyMax) || energyMin >= energyMax) {
    std::cerr << "[CDet occupancy response] Set analysis.ecal_energy_min/max to a valid GeV interval.\n";
    return;
  }
  if (!std::isfinite(leCutMinNs) || !std::isfinite(leCutMaxNs) || leCutMinNs >= leCutMaxNs ||
      (applyToTCut && (!std::isfinite(totCutMinNs) || !std::isfinite(totCutMaxNs) ||
       totCutMinNs < 0.0 || totCutMinNs >= totCutMaxNs))) {
    std::cerr << "[CDet occupancy response] Invalid analysis.le_min/max or enabled analysis.tot_min/max (ns).\n";
    return;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet occupancy response] Ellipse radius must be -1 (off) or positive; enabled ellipses require finite centers and positive finite scales.\n";
    return;
  }
  if (runNumber <= 0 || bar < 0 || bar >= nBars || controlBar < 0 ||
      controlBar >= nBars || events < -1 || !std::isfinite(ecalXBinWidthM) ||
      ecalXBinWidthM <= 0.0 || !std::isfinite(ecalXMinM) ||
      !std::isfinite(ecalXMaxM) || ecalXMinM >= ecalXMaxM) {
    std::cerr << "[CDet occupancy response] Invalid run, bar, event limit or ECal x binning.\n";
    return;
  }
  const int firstPixel = pixelsPerBar * bar;
  const int firstControlPixel = pixelsPerBar * controlBar;
  const int nECalX = int(std::ceil((ecalXMaxM - ecalXMinM) / ecalXBinWidthM));
  const double tdcToNs = 0.01; // Legacy LE units; corrected dt and ToT are already ns.

  // 1. Load one run, including its ROOT rollover files.
  const char *directory = inputDirectory && inputDirectory[0]
      ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0]) {
    std::cerr << "[CDet occupancy response] Set OUT_DIR or pass inputDirectory.\n";
    return;
  }
  TChain chain("T");
  if (CDetRunDataset::AddToChain(&chain, runNumber, directory, minSegment, maxSegment) <= 0) return;
  const Long64_t total = chain.GetEntries();
  if (total <= 0) {
    std::cerr << "[CDet occupancy response] No events found.\n";
    return;
  }
  // -1 reads all entries; a nonnegative value limits the number processed.
  const Long64_t limit = events >= 0 ? std::min(events, total) : total;

  // Prepare a tagged output directory before processing the events.
  TString tag = outputTag ? outputTag : config.GetValue("occupancy_response.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d_bar%d_control%d", runNumber, bar, controlBar);
  const TString allowedTagCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowedTagCharacters.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet occupancy response] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/OccupancyResponse/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet occupancy response] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }
  std::cout << "[CDet occupancy response] PDF output: " << outputDirectory << '\n';

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
  TTreeReaderArray<double> ecalEligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<double> spatialPass(reader, "earm.cdet.pulse.spatial_pass");


  // 3. Define occupancy, reference and response histograms.
  delete hECalX;
  delete hECalXReference;
  hECalX = new TH1D(TString::Format("hECalX_run%d", runNumber), TString::Format("Run %d: event with a final-selected CDet pulse;ECal x (m);Events", runNumber), nECalX, ecalXMinM, ecalXMaxM);
  hECalXReference = new TH1D(TString::Format("hECalXReference_run%d", runNumber), TString::Format("Run %d: ECal energy/time reference, no CDet requirement;ECal x (m);Events", runNumber), nECalX, ecalXMinM, ecalXMaxM);
  for (TH1D* hist : {hECalX, hECalXReference}) {
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

  std::cout << "  Occupancy baseline: " << (occupancyUseECalBaseline ? "ECal selected" : "all stored complete pulses")
            << "; x match |pulse.ecal_x_residual| <= " << occupancyDxMaxM << " m.\n"
            << "  Occupancy counts each event once per pixel. ECal selection uses energy,\n"
            << "  replay eligibility, and finite stored dx; no LE/ToT window, y match, or ellipse.\n";

  // 4. Evaluate event/pulse selections and fill once per event/pixel.
  Long64_t processed = 0, energyEvents = 0;
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

      if (good_final) {
        eventHasGoodPulse = true;
        if (good_raw_occupancy) seenSelected[pixelID] = true;
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
      std::cout << "[CDet occupancy response] " << processed << '/' << limit
                << " events; ECal energy passed: " << energyEvents
                << "; ECal-only reference events: " << referenceEvents << '\n';
  } //end event loop

  if (processed != limit) {
    std::cerr << "[CDet occupancy response] Reader stopped after " << processed
              << " of " << limit << " requested events; status "
              << reader.GetEntryStatus() << '\n';
    return;
  }
  // 5. Draw occupancy, opportunity counts, fractions and ECal reference.
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

  TCanvas* outputCanvases[] = {cPixelOccupancy, cProjectionCounts,
      cProjectionFractions, cECalXSelection};
  for (TCanvas* canvas : outputCanvases)
    canvas->SaveAs(TString::Format("%s/%s.pdf", outputDirectory.Data(), canvas->GetName()));
  std::cout << "[CDet occupancy response] Read " << processed << " events; "
            << energyEvents << " passed the ECal energy cut. Drew four canvases.\n"
            << "  PDF/CSV output directory: " << outputDirectory << '\n';
  std::cout << "[CDet projection] " << referenceEvents << " ECal-only reference events; "
            << referenceWithoutPulses << " have no stored CDet pulses; " << referenceWithoutSelectedPulses
            << " have no final-selected CDet pulse.\n";
}

void Plot_CDet_OccupancyResponse(
    int runNumber, int bar = -1, const char* inputDirectory = nullptr,
    Long64_t eventsOverride = -2, double ecalXBinWidthM = 0.0425,
    double ecalXMinM = -1.5, double ecalXMaxM = 1.5, const char* outputTag = nullptr)
{
  const TString configFile = TString::Format("CDet_run%d_occupancy_response.conf", runNumber);
  Plot_CDet_OccupancyResponse(configFile.Data(), bar, inputDirectory, eventsOverride,
                            ecalXBinWidthM, ecalXMinM, ecalXMaxM, outputTag);
}

void Plot_CDet_OccupancyResponse(
    const char* configFile, const char* outputTag, int bar = -1,
    const char* inputDirectory = nullptr, Long64_t eventsOverride = -2)
{
  Plot_CDet_OccupancyResponse(configFile, bar, inputDirectory, eventsOverride,
                            0.0425, -1.5, 1.5, outputTag);
}

void Plot_CDet_OccupancyResponse(
    int runNumber, const char* outputTag, int bar = -1,
    const char* inputDirectory = nullptr, Long64_t eventsOverride = -2)
{
  const TString configFile = TString::Format("CDet_run%d_occupancy_response.conf", runNumber);
  Plot_CDet_OccupancyResponse(configFile.Data(), outputTag, bar, inputDirectory, eventsOverride);
}
