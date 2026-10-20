// Read-only pulse-based offset scan. Run from scripts/cdet with OUT_DIR set:
// root -l -b -q 'CDet_CrossTargetPulseOffsetScan.C+()'
// Optional list: one "run [projection.conf]" per line. Run 5710 runs first.
// Keeps the PulseStudy architecture: load aligned event/hit vectors once,
// then fill histograms and fit them separately. Only one run stays in memory.
// No dependency on PlotElastic or PulseStudy; no calibration .dat is written.

#include <TChain.h>
#include <TEnv.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "CDetRunDataset.h"

namespace CDetPulseOffsetScan {
constexpr int nPixels = 2688, nGroups = 336, nBars = 168;
// All vectors share [event][hit]; clear them before loading the next run.
std::vector<std::vector<int>> vPixel;
std::vector<std::vector<double>> vDt;
TH1D* hPixelDt[nPixels] = {};
TH1D* hHalfDt[nGroups] = {};
TH1D* hBarDt[nBars] = {};

// Fit result slots, shared by the pixel/half/bar arrays below. No custom structs.
enum { Mean, MeanError, Sigma, Chi2Ndf, Significance, FitStatus,
       WindowEntries, FitLow, FitHigh, Valid, nFitValues };
using FitValues = std::array<double, nFitValues>;

std::string CSV(const std::string& value) {
  std::string result = "\"";
  for (char c : value) result += c == '"' ? "\"\"" : std::string(1, c);
  return result + '"';
}
void Number(std::ostream& out, double value) { if (std::isfinite(value)) out << value; }
double Median(std::vector<double> values) {
  if (values.empty()) return NAN;
  std::sort(values.begin(), values.end());
  const size_t middle = values.size() / 2;
  return values.size() % 2 ? values[middle] : 0.5 * (values[middle-1] + values[middle]);
}

// The reviewed table supplies the known unused-slot mask and the historical
// comparison. Parse numeric tokens with stod: its rows legitimately contain nan.
bool ReadReviewed(const char* path, std::vector<double>& offsets,
                  std::vector<long long>& entries, std::vector<std::string>& sources,
                  std::vector<bool>& unused) {
  offsets.assign(nPixels, NAN); entries.assign(nPixels, 0);
  sources.assign(nPixels, "missing"); unused.assign(nPixels, false);
  std::ifstream input(path);
  if (!input) { std::cerr << "[Pulse offset scan] Cannot read reviewed reference: " << path << '\n'; return false; }
  std::vector<bool> seen(nPixels, false);
  std::string line;
  int count = 0;
  while (std::getline(input, line)) {
    std::istringstream row(line);
    std::vector<std::string> tokens;
    std::string token;
    while (row >> token) tokens.push_back(token);
    if (tokens.empty() || tokens[0][0] == '#') continue;
    try {
      if (tokens.size() < 14) throw std::runtime_error("short row");
      const int pixel = std::stoi(tokens[0]);
      if (pixel < 0 || pixel >= nPixels || seen[pixel]) throw std::runtime_error("invalid/duplicate pixel");
      seen[pixel] = true; ++count;
      entries[pixel] = std::stoll(tokens[3]); sources[pixel] = tokens[5];
      unused[pixel] = sources[pixel] == "unused";
      const double offset = std::stod(tokens[12]);
      if (!unused[pixel] && std::isfinite(offset) &&
          (sources[pixel] != "retained_existing" || offset != 0.0)) offsets[pixel] = offset;
    } catch (const std::exception& error) {
      std::cerr << "[Pulse offset scan] Malformed reviewed table: " << error.what() << '\n'; return false;
    }
  }
  if (count != nPixels) { std::cerr << "[Pulse offset scan] Reviewed table must contain all 2688 slots.\n"; return false; }
  return true;
}

bool ReadShift(const TString& path, double& shift) {
  shift = 0.0;
  std::ifstream input(path.Data());
  if (!input) return false;
  std::string line, section;
  while (std::getline(input, line)) {
    line = line.substr(0, line.find('#'));
    std::istringstream row(line);
    std::string key;
    if (!(row >> key)) continue;
    if (key[0] == '[') { section = key; continue; }
    if (section != "[GlobalTiming]" || key != "shift_ns") continue;
    row >> std::ws;
    if (row.peek() == '=') row.get();
    return bool(row >> shift) && std::isfinite(shift);
  }
  return false;
}

// Store only the fields needed for the scan; selection is explicit here.
bool LoadRun(int run, const char* configFile, const char* inputDirectory,
             Long64_t events, const std::vector<bool>& unused,
             Long64_t& processed, Long64_t& ecalEvents, std::string& cuts,
             std::string& status) {
  vPixel.clear(); vDt.clear(); processed = 0; ecalEvents = 0;
  TEnv config;
  if (config.ReadFile(configFile, kEnvLocal) != 0 || config.GetValue("config.version", 0) != 1 ||
      config.GetValue("analysis.run_number", 0) != run) { status = "invalid_config"; return false; }
  const double ECalMinE = config.GetValue("analysis.ecal_energy_min", 1.0);
  const double ECalMaxE = config.GetValue("analysis.ecal_energy_max", 12.0);
  const double ECalMinT = config.GetValue("analysis.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("analysis.ecal_time_max", NAN);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const int minSeg = config.GetValue("analysis.min_segment", -1);
  const int maxSeg = config.GetValue("analysis.max_segment", -1);
  if (!std::isfinite(ECalMinE) || !std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE ||
      !std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      !std::isfinite(TotMin) || !std::isfinite(TotMax) || TotMin < 0 || TotMin >= TotMax ||
      minSeg < -1 || maxSeg < -1 || (minSeg >= 0 && maxSeg >= 0 && minSeg > maxSeg)) {
    status = "invalid_cuts"; return false;
  }
  std::ostringstream description;
  description << "ECal energy " << ECalMinE << " < E < " << ECalMaxE
      << " GeV; ECal ADC time " << ECalMinT << " < t < " << ECalMaxT
      << " ns; LE " << LeMin << " <= t <= " << LeMax
      << " ns; ToT " << TotMin << " < ToT < " << TotMax
      << " ns; replay calib_valid/ecal_eligible/spatial_pass; no ellipse or polygons";
  cuts = description.str();
  TChain chain("T");
  if (CDetRunDataset::AddToChain(&chain, run, inputDirectory, minSeg, maxSeg) <= 0) { status = "missing_input"; return false; }
  const Long64_t entries = chain.GetEntries();
  if (entries <= 0) { status = "empty_input"; return false; }
  const Long64_t limit = events < 0 ? entries : std::min(events, entries);
  TTreeReader reader(&chain);
  TTreeReaderValue<double> ECalE(reader, "earm.ecal.e");
  TTreeReaderValue<double> ECalTime(reader, "earm.ecal.adctime");
  TTreeReaderArray<double> PulseID(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<double> PulseLe(reader, "earm.cdet.pulse.tdc_le");
  TTreeReaderArray<double> PulseToT(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<double> PulseCalib(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<double> PulseECal(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<double> PulseSpatial(reader, "earm.cdet.pulse.spatial_pass");
  std::cout << "[Pulse offset scan] Run " << run << ": " << limit << " events; " << cuts << std::endl;
  while (processed < limit && reader.Next()) {
    ++processed;
    if (processed % 1000 == 0) std::cout << "[Pulse offset scan] Run " << run << ' ' << processed << '/' << limit << " events" << std::endl;
    const bool good_ECal_energy = std::isfinite(*ECalE) && *ECalE > ECalMinE && *ECalE < ECalMaxE;
    const bool good_ECal_time = std::isfinite(*ECalTime) && *ECalTime > ECalMinT && *ECalTime < ECalMaxT;
    const bool good_ECal_event = good_ECal_energy && good_ECal_time;
    if (!good_ECal_event) continue;
    ++ecalEvents;
    std::vector<int> thisEventPixel;
    std::vector<double> thisEventDt;
    for (size_t hit = 0; hit < PulseID.GetSize(); ++hit) {
      const double id = PulseID[hit], le = 0.01 * PulseLe[hit], tot = PulseToT[hit];
      const bool good_id = std::isfinite(id) && id >= 0 && id < nPixels && std::floor(id) == id;
      if (!good_id || unused[int(id)]) continue;
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(tot) && tot > TotMin && tot < TotMax;
      const bool good_calibration = PulseCalib[hit] > 0.5;
      const bool good_ecal_eligible = PulseECal[hit] > 0.5;
      const bool good_spatial = PulseSpatial[hit] > 0.5;
      const bool good_hit = good_le && good_tot && good_calibration && good_ecal_eligible && good_spatial;
      if (!good_hit) continue;
      thisEventPixel.push_back(int(id));
      // Low-level calibrated LE, before pixel/ECal/timewalk/global corrections.
      thisEventDt.push_back(*ECalTime - le);
    }
    if (thisEventPixel.empty()) continue;
    vPixel.push_back(std::move(thisEventPixel));
    vDt.push_back(std::move(thisEventDt));
  }
  if (processed != limit) { status = "reader_stopped"; return false; }
  status = "loaded"; return true;
}

double EntriesIn(TH1D* hist, double low, double high) {
  // Count the bins whose centers ROOT includes in this fit interval.
  double entries = 0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
    if (hist->GetBinCenter(bin) >= low && hist->GetBinCenter(bin) <= high) entries += hist->GetBinContent(bin);
  return entries;
}

FitValues FitPeak(TH1D* hist, double low, double high, int minEntries,
                  double seed = NAN, double seedSigma = 3.0) {
  FitValues result; result.fill(NAN);
  result[FitStatus] = -1; result[Valid] = 0;
  result[FitLow] = low; result[FitHigh] = high;
  result[WindowEntries] = EntriesIn(hist, low, high);
  if (result[WindowEntries] < minEntries) return result;
  int peakBin = 0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    if (hist->GetBinCenter(bin) < low || hist->GetBinCenter(bin) > high) continue;
    if (!peakBin || hist->GetBinContent(bin) > hist->GetBinContent(peakBin)) peakBin = bin;
  }
  if (!peakBin) return result;
  if (!std::isfinite(seed)) seed = hist->GetBinCenter(peakBin);
  TF1 fit(TString::Format("fit_%s", hist->GetName()), "gaus(0)+pol1(3)", low, high);
  fit.SetParameters(std::max(1.0, hist->GetBinContent(hist->FindBin(seed))), seed, std::min(3.0, seedSigma), 0.0, 0.0);
  fit.SetParLimits(0, 0.0, 1.5 * std::max(1.0, hist->GetMaximum()));
  fit.SetParLimits(1, low, high); fit.SetParLimits(2, 0.5, 8.0);
  const int status = hist->Fit(&fit, "RQSN0");
  result[FitStatus] = status; result[Mean] = fit.GetParameter(1);
  result[MeanError] = fit.GetParError(1); result[Sigma] = std::fabs(fit.GetParameter(2));
  result[Chi2Ndf] = fit.GetNDF() > 0 ? fit.GetChisquare() / fit.GetNDF() : NAN;
  result[Significance] = fit.GetParError(0) > 0 ? fit.GetParameter(0) / fit.GetParError(0) : NAN;
  const double margin = 0.5 * hist->GetBinWidth(1);
  result[Valid] = status == 0 && std::isfinite(result[Mean]) &&
      result[Mean] > low + margin && result[Mean] < high - margin &&
      std::isfinite(result[MeanError]) && result[MeanError] > 0 && result[MeanError] <= 2.0 &&
      result[Sigma] >= 0.5 && result[Sigma] <= 8.0 &&
      std::isfinite(result[Chi2Ndf]) && result[Chi2Ndf] <= 15.0 &&
      std::isfinite(result[Significance]) && result[Significance] >= 1.5 && fit.GetParameter(0) > 0;
  // Retain the attempted fit for inspection, including failed quality checks.
  hist->GetListOfFunctions()->Add(fit.Clone());
  return result;
}

void FillHistograms(int run, double Width, double HistMin, double HistMax,
                    std::vector<long long>& pixelEvents) {
  const int bins = int(std::ceil((HistMax - HistMin) / Width));
  pixelEvents.assign(nPixels, 0);
  for (int pixel = 0; pixel < nPixels; ++pixel) {
    delete hPixelDt[pixel];
    hPixelDt[pixel] = new TH1D(TString::Format("hPixelDt_%d", pixel), TString::Format("Run %d pixel %d;t_{ECal}-0.01 t_{CDet,LE} (ns);Pulses", run, pixel), bins, HistMin, HistMax);
    hPixelDt[pixel]->SetDirectory(nullptr);
  }
  for (int group = 0; group < nGroups; ++group) {
    delete hHalfDt[group];
    hHalfDt[group] = new TH1D(TString::Format("hHalfDt_%d", group), TString::Format("Run %d bar %d half %d;t_{ECal}-0.01 t_{CDet,LE} (ns);Pulses", run, group/2, group%2), bins, HistMin, HistMax);
    hHalfDt[group]->SetDirectory(nullptr);
  }
  for (int bar = 0; bar < nBars; ++bar) {
    delete hBarDt[bar];
    hBarDt[bar] = new TH1D(TString::Format("hBarDt_%d", bar), TString::Format("Run %d bar %d;t_{ECal}-0.01 t_{CDet,LE} (ns);Pulses", run, bar), bins, HistMin, HistMax);
    hBarDt[bar]->SetDirectory(nullptr);
  }
  for (size_t event = 0; event < vPixel.size(); ++event) {
    std::array<bool, nPixels> seen{};
    for (size_t hit = 0; hit < vPixel[event].size(); ++hit) {
      const int id = vPixel[event][hit];
      const double dt = vDt[event][hit];
      hPixelDt[id]->Fill(dt); hHalfDt[id/8]->Fill(dt); hBarDt[id/16]->Fill(dt);
      if (!seen[id]) { ++pixelEvents[id]; seen[id] = true; }
    }
  }
}
} // namespace CDetPulseOffsetScan

void CDet_CrossTargetPulseOffsetScan(
    const char* runList = "", const char* outputDirectory = "PulseOffsetScan",
    const char* inputDirectory = nullptr, Long64_t events = -1,
    const char* reviewed5710 = "CDet_pixel_timing_fit_results_run5710_reviewed_102_cut_pass.dat",
    double Width = 1.0, double HistMin = -100.0, double HistMax = 100.0,
    double BroadFitMin = -30.0, double BroadFitMax = 10.0,
    int minPixelEntries = 35, int minHalfEntries = 100, int minBarEntries = 100) {
  using namespace CDetPulseOffsetScan;
  const char* directory = inputDirectory && inputDirectory[0] ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0] || !outputDirectory || !outputDirectory[0] || events < -1 ||
      !std::isfinite(Width) || Width <= 0 || !std::isfinite(HistMin) || !std::isfinite(HistMax) || HistMin >= HistMax ||
      !std::isfinite(BroadFitMin) || !std::isfinite(BroadFitMax) || BroadFitMin >= BroadFitMax ||
      (HistMax-HistMin)/Width > 100000 || minPixelEntries < 1 || minHalfEntries < 1 || minBarEntries < 1) {
    std::cerr << "[Pulse offset scan] Set OUT_DIR/inputDirectory and valid output, axes, limits, thresholds.\n"; return;
  }
  std::vector<double> reviewed;
  std::vector<long long> reviewedEntries;
  std::vector<std::string> reviewedSources;
  std::vector<bool> unused;
  if (!ReadReviewed(reviewed5710, reviewed, reviewedEntries, reviewedSources, unused)) return;
  // Approved inventory, not crossRuns.txt (which still contains 5976/5979).
  const int defaults[] = {5710,3573,3575,3602,3603,3604,3648,3682,3685,3686,3694,
      3757,3758,3786,3788,3844,3845,3846,4003,4006,4344,4345,4400,4722,4735,
      4926,4929,4985,5066,5068,5290,5292,5294,5423,5722,5723,5724,5726,5794,5795,5992};
  std::vector<int> runs{5710};
  std::vector<std::string> configs{"CDet_run5710_projection.conf"};
  if (runList && runList[0]) {
    std::ifstream input(runList);
    if (!input) { std::cerr << "[Pulse offset scan] Cannot read run list.\n"; return; }
    std::string line;
    while (std::getline(input, line)) {
      std::istringstream row(line.substr(0, line.find('#')));
      int run; std::string config;
      if (!(row >> run)) continue;
      if (run <= 0) { std::cerr << "[Pulse offset scan] Invalid run number.\n"; return; }
      row >> config;
      if (config.empty()) config = TString::Format("CDet_run%d_projection.conf", run).Data();
      if (run == 5710) { configs[0] = config; continue; }
      if (std::find(runs.begin(), runs.end(), run) == runs.end()) { runs.push_back(run); configs.push_back(config); }
    }
  } else {
    for (int run : defaults) if (run != 5710) {
      runs.push_back(run); configs.push_back(TString::Format("CDet_run%d_projection.conf", run).Data());
    }
  }
  // New directory only: a repeated invocation cannot silently replace a scan.
  if (!gSystem->AccessPathName(outputDirectory) || gSystem->mkdir(outputDirectory, true) != 0) {
    std::cerr << "[Pulse offset scan] Choose a NEW output directory: " << outputDirectory << '\n'; return;
  }
  std::ofstream table(TString::Format("%s/pixel_offsets.csv", outputDirectory));
  std::ofstream summary(TString::Format("%s/run_summary.csv", outputDirectory));
  if (!table || !summary) { std::cerr << "[Pulse offset scan] Cannot create CSV outputs.\n"; return; }
  std::ofstream settings(TString::Format("%s/scan_settings.txt", outputDirectory));
  settings << std::setprecision(17)
      << "input_directory " << directory << "\nreviewed5710 " << reviewed5710
      << "\nevents_per_run " << events << "\nhistogram_width_min_max_ns " << Width << ' ' << HistMin << ' ' << HistMax
      << "\nlegacy_fit_min_max_ns " << BroadFitMin << ' ' << BroadFitMax
      << "\nmin_pixel_half_bar_entries " << minPixelEntries << ' ' << minHalfEntries << ' ' << minBarEntries
      << "\nfit_model gaus(0)+pol1(3)\nsigma_bounds_ns 0.5 8\nmax_chi2_ndf 15\nmax_mean_error_ns 2\nmin_amplitude_significance 1.5\n"
      << "dt_definition ECal ADC time minus 0.01*pulse.tdc_le; no full timing corrections\n"
      << "offset_definition fitted mean minus per-run detector reference median\n"
      << "fit_window_translation add legacy shift_ns to fit bounds in raw-dt coordinates\n"
      << "fit_mean_error_is_conditional excludes reference uncertainty and fallback systematics\n";
  for (size_t i = 0; i < runs.size(); ++i) settings << "run_config " << runs[i] << ' ' << configs[i] << '\n';
  settings.close();
  if (!settings) { std::cerr << "[Pulse offset scan] Cannot write scan settings.\n"; return; }
  table << std::setprecision(10); summary << std::setprecision(10);
  table << "run,pixel,layer,bar,half,run_status,source,pixel_events,pixel_pulses,pixel_fit_entries,half_pulses,half_fit_entries,bar_pulses,bar_fit_entries,used_fit_entries,pixel_underflow,pixel_overflow,fit_status,mean_ns,mean_err_ns,sigma_ns,chi2_ndf,signal_significance,fit_low_ns,fit_high_ns,detector_reference_ns,reference_source,reference_count,offset_ns,run5710_offset_ns,minus_run5710_ns,run5710_source,run5710_pixel_fit_entries,run5710_used_fit_entries,reviewed5710_offset_ns,minus_reviewed5710_ns,reviewed5710_entries,reviewed5710_source,pixel_fit_status,half_fit_status,bar_fit_status,pixel_attempt_entries,pixel_fit_valid,half_fit_valid,bar_fit_valid,reason\n";
  summary << "run,status,config,processed_events,ecal_events,pulse_events,selected_pulses,pixel_fits,half_fallbacks,bar_fallbacks,unavailable,unused,detector_reference_ns,reference_source,reference_count,legacy_run_shift_ns,shift_loaded,fit_min_ns,fit_max_ns,cuts\n";
  std::vector<double> referenceOffsets(nPixels, NAN), referencePixelEntries(nPixels, NAN), referenceUsedEntries(nPixels, NAN);
  std::vector<std::string> referenceSources(nPixels, "unavailable");
  for (size_t r = 0; r < runs.size(); ++r) {
    const int run = runs[r];
    std::cout << "[Pulse offset scan] Starting " << r+1 << '/' << runs.size() << ": run " << run << std::endl;
    double runShift = 0;
    const TString timingFile = TString::Format("%s/CDet_run%d.dat", gSystem->DirName(configs[r].c_str()), run);
    const bool shiftLoaded = ReadShift(timingFile, runShift);
    if (!shiftLoaded) { runShift = 0; std::cerr << "[Pulse offset scan] No valid run shift in " << timingFile << "; using 0 for fit-window translation.\n"; }
    // The legacy stage-0 scan adds shift_ns to LE. Our histograms keep raw dt;
    // translate its search window instead. Relative offsets are unchanged by
    // a uniform shift when the same peaks are fit.
    const double fitMin = BroadFitMin + runShift, fitMax = BroadFitMax + runShift;
    Long64_t processed = 0, ecalEvents = 0;
    std::string cuts, status;
    bool loaded = false;
    if (fitMin < HistMin || fitMax > HistMax) status = "fit_window_outside_histogram";
    else loaded = LoadRun(run, configs[r].c_str(), directory, events, unused, processed, ecalEvents, cuts, status);
    std::vector<long long> pixelEvents(nPixels, 0);
    std::vector<FitValues> pixelFits(nPixels), halfFits(nGroups), barFits(nBars);
    FitValues invalid; invalid.fill(NAN); invalid[Valid] = 0; invalid[FitStatus] = -1;
    std::fill(pixelFits.begin(), pixelFits.end(), invalid);
    std::fill(halfFits.begin(), halfFits.end(), invalid);
    std::fill(barFits.begin(), barFits.end(), invalid);
    if (loaded) {
      FillHistograms(run, Width, HistMin, HistMax, pixelEvents);
      std::cout << "[Pulse offset scan] Run " << run << ": fitting pixels, half-NINOs, and bars." << std::endl;
      for (int bar = 0; bar < nBars; ++bar) barFits[bar] = FitPeak(hBarDt[bar], fitMin, fitMax, minBarEntries);
      for (int half = 0; half < nGroups; ++half) halfFits[half] = FitPeak(hHalfDt[half], fitMin, fitMax, minHalfEntries);
      for (int pixel = 0; pixel < nPixels; ++pixel) {
        if (unused[pixel]) continue;
        // As in the earlier hierarchical extractor, an accepted half fit guides
        // the pixel peak away from cross-talk. Otherwise use a broad pixel fit.
        const auto& half = halfFits[pixel/8];
        if (half[Valid]) {
          const double width = std::min(8.0, std::max(4.0, 2.0*half[Sigma]));
          pixelFits[pixel] = FitPeak(hPixelDt[pixel], std::max(fitMin, half[Mean]-width),
              std::min(fitMax, half[Mean]+width), minPixelEntries, half[Mean], half[Sigma]);
        } else pixelFits[pixel] = FitPeak(hPixelDt[pixel], fitMin, fitMax, minPixelEntries);
      }
    }
    // Reference counts independent fits, never repeated inherited pixel rows.
    std::vector<double> centroids;
    std::string referenceSource = "pixel_median";
    for (const auto& fit : pixelFits) if (fit[Valid]) centroids.push_back(fit[Mean]);
    if (centroids.empty()) {
      referenceSource = "half_median";
      for (const auto& fit : halfFits) if (fit[Valid]) centroids.push_back(fit[Mean]);
    }
    if (centroids.empty()) {
      referenceSource = "bar_median";
      for (const auto& fit : barFits) if (fit[Valid]) centroids.push_back(fit[Mean]);
    }
    const double detectorReference = Median(centroids);
    if (centroids.empty()) referenceSource = "unavailable";
    if (loaded) status = std::isfinite(detectorReference) ? "ok" : "no_valid_reference";
    int pixelCount = 0, halfCount = 0, barCount = 0, missingCount = 0, unusedCount = 0;
    double selectedPulses = 0;
    for (int pixel = 0; pixel < nPixels; ++pixel) {
      const auto& p = pixelFits[pixel]; const auto& h = halfFits[pixel/8]; const auto& b = barFits[pixel/16];
      FitValues used = invalid;
      std::string source = "unavailable", reason = status;
      if (unused[pixel]) { source = "unused"; reason = "known unused electronic slot"; ++unusedCount; }
      else if (loaded) {
        if (p[Valid]) { used = p; source = "pixel"; reason = "accepted pixel fit"; ++pixelCount; }
        else if (h[Valid]) { used = h; source = "half_nino"; reason = "pixel insufficient statistics or failed fit quality"; ++halfCount; }
        else if (b[Valid]) { used = b; source = "bar"; reason = "pixel and half insufficient statistics or failed fit quality"; ++barCount; }
        else { reason = "no acceptable pixel half or bar fit"; ++missingCount; }
      } else ++missingCount;
      const double offset = used[Valid] && std::isfinite(detectorReference) ? used[Mean]-detectorReference : NAN;
      const double pixelFitEntries = loaded ? EntriesIn(hPixelDt[pixel], fitMin, fitMax) : NAN;
      if (run == 5710) {
        referenceOffsets[pixel] = offset; referencePixelEntries[pixel] = pixelFitEntries;
        referenceUsedEntries[pixel] = used[WindowEntries]; referenceSources[pixel] = source;
      }
      const double pulses = loaded ? hPixelDt[pixel]->GetEntries() : NAN;
      if (loaded) selectedPulses += pulses;
      table << run << ',' << pixel << ',' << pixel/1344+1 << ',' << pixel/16 << ',' << pixel%16/8
            << ',' << status << ',' << source << ',';
      Number(table, loaded ? double(pixelEvents[pixel]) : NAN); table << ',';
      for (double value : {pulses, pixelFitEntries,
          loaded ? hHalfDt[pixel/8]->GetEntries() : NAN,
          loaded ? EntriesIn(hHalfDt[pixel/8], fitMin, fitMax) : NAN,
          loaded ? hBarDt[pixel/16]->GetEntries() : NAN,
          loaded ? EntriesIn(hBarDt[pixel/16], fitMin, fitMax) : NAN,
          used[WindowEntries], loaded ? hPixelDt[pixel]->GetBinContent(0) : NAN,
          loaded ? hPixelDt[pixel]->GetBinContent(hPixelDt[pixel]->GetNbinsX()+1) : NAN,
          used[FitStatus], used[Mean], used[MeanError], used[Sigma], used[Chi2Ndf], used[Significance],
          used[FitLow], used[FitHigh], detectorReference}) { Number(table, value); table << ','; }
      table << referenceSource << ',' << centroids.size() << ',';
      for (double value : {offset, referenceOffsets[pixel], offset-referenceOffsets[pixel]}) { Number(table, value); table << ','; }
      table << referenceSources[pixel] << ',';
      for (double value : {referencePixelEntries[pixel], referenceUsedEntries[pixel], reviewed[pixel], offset-reviewed[pixel]}) { Number(table, value); table << ','; }
      table << reviewedEntries[pixel] << ',' << reviewedSources[pixel] << ',' << p[FitStatus] << ',' << h[FitStatus] << ',' << b[FitStatus] << ',';
      Number(table, p[WindowEntries]);
      table << ',' << p[Valid] << ',' << h[Valid] << ',' << b[Valid] << ',' << CSV(reason) << '\n';
    }
    summary << run << ',' << status << ',' << CSV(configs[r]) << ',' << processed << ',' << ecalEvents << ','
            << (loaded ? vPixel.size() : 0) << ',' << selectedPulses << ',' << pixelCount << ',' << halfCount << ','
            << barCount << ',' << missingCount << ',' << unusedCount << ',';
    Number(summary, detectorReference);
    summary << ',' << referenceSource << ',' << centroids.size() << ',' << runShift << ',' << shiftLoaded << ',' << fitMin << ',' << fitMax << ',' << CSV(cuts) << '\n';
    table.flush(); summary.flush();
    if (!table || !summary) { std::cerr << "[Pulse offset scan] CSV write failed; stopping.\n"; return; }
    if (loaded) {
      TFile output(TString::Format("%s/run%d_dt.root", outputDirectory, run), "RECREATE");
      if (output.IsZombie()) { std::cerr << "[Pulse offset scan] Cannot save ROOT diagnostics; CSV is available.\n"; return; }
      for (TH1D* hist : hPixelDt) hist->Write();
      for (TH1D* hist : hHalfDt) hist->Write();
      for (TH1D* hist : hBarDt) hist->Write();
      TNamed("selection", cuts.c_str()).Write();
      TNamed("timing_definition", "dt = earm.ecal.adctime - 0.01*earm.cdet.pulse.tdc_le; offset = fitted mean - reference median").Write();
      TParameter<double>("detector_reference_ns", detectorReference).Write();
      TParameter<double>("legacy_run_shift_ns", runShift).Write();
      output.Close();
    }
    std::cout << "[Pulse offset scan] Run " << run << ": " << status << "; pixel/half/bar = "
              << pixelCount << '/' << halfCount << '/' << barCount << "; unavailable = " << missingCount << std::endl;
    vPixel.clear(); vDt.clear();
  }
  std::cout << "[Pulse offset scan] Finished. CSV tables and per-run dt histograms: " << outputDirectory
            << "\n  No calibration files were changed.\n";
}
