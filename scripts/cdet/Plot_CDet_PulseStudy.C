// Read the run once; keep vectors available for later plotting, like PlotElastic.
// From gocdetscripts, with OUT_DIR set:
//   .L Plot_CDet_PulseStudy.C+
//   Plot_CDet_PulseStudy(6077, 1000);
//   plotPulseStudy("CDet_run6077_pulse_pixel_timing.conf", "baseline");
// Or override display settings and plot only bar 19:
//   plotPulseStudy("CDet_run6077_pulse_pixel_timing.conf", "bar19", 1., -150., 150., 19, 19);
// Saves four PDFs per bar under PusleStudy/<tag>/ without opening plot windows.
// Interactive study/control plots (eight visible canvases per bar):
//   plotPulseStudyBars("CDet_run6077_pulse_pixel_timing.conf", "first_pass");
// Position/rx overlays (one visible canvas per bar):
//   plotPulseStudyPositions("CDet_run6077_pulse_pixel_timing.conf", "positions", 19, 30);
// Background plots, including matched rx canvases for multiple/single pulses in bar 30:
//   plotPulseStudyBackground("CDet_run6077_pulse_pixel_timing.conf", "multiplicity", 19, 30);
// Shares CDet_run<run>_pulse_pixel_timing.conf with the pulse-timing macro.
// Cut arguments override that file for this call only; NAN means use the file.
// Loading always applies TotMin/TotMax. The secondary plotting ToT switch is
// pulse_study.apply_tot_cut; it does not disable the loading cut.
// Add your pulse readers, vectors and hit cuts in the marked sections below.

#include <TChain.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TList.h>
#include <TROOT.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "CDetRunDataset.h"

// One element per ECal-selected event with at least one selected CDet pulse.
// All pulse vectors share [event][hit].

//Global constants
const int pixelsPerBar = 16;
const int numBars = 168;
constexpr double kECalZFromTargetM = 6.144;

// Histograms from the most recent plotPulseStudy() call, available in ROOT.
TH1D* hStudyBarDt[numBars] = {};
TH1D* hStudyPixelDt[numBars][pixelsPerBar] = {};
TH2D* hStudyBarLeVsToT[numBars] = {};
TH2D* hStudyPixelLeVsToT[numBars][pixelsPerBar] = {};

// Separate summary histograms: index 0 = study bar, 1 = control bar.
TH1D* hPulseSummaryBarDt[2] = {};
TH1D* hPulseSummaryBarLe[2] = {};
TH1D* hPulseSummaryPixelDt[2][pixelsPerBar] = {};
TH1D* hPulseSummaryPixelLe[2][pixelsPerBar] = {};
TH2D* hPulseSummaryBarLeVsToT[2] = {};
TH2D* hPulseSummaryPixelLeVsToT[2][pixelsPerBar] = {};
TH2D* hPulseBarDtVsRx[2] = {nullptr, nullptr};
TH2D* hPulsePixelDtVsRx[2][pixelsPerBar] = {};
// Separate position-overlay histograms: index 0 = study bar, 1 = control bar.
TH1D* hPulsePositionBarX[2] = {};
TH1D* hPulsePositionBarRx[2] = {};
// Multiplicity is counted per event in the background routine's selected bar.
TH1D* hPulseBackgroundRxMultiple = nullptr;
TH1D* hPulseBackgroundRxSingle = nullptr;
TH2D* hPulseBackgroundDtVsRxMultiple = nullptr;
TH2D* hPulseBackgroundDtVsRxSingle = nullptr;

//-------ECal Vectors ---------
std::vector<double> vECalE;
std::vector<double> vECalAdcTime;
std::vector<double> vECalX;
std::vector<double> vECalY;
std::vector<std::vector<double>> vECalXProjection;
std::vector<std::vector<double>> vECalYProjection;


//-------CDet Vectors ---------
std::vector<std::vector<double>> vPulseID; // Example pulse vector, indexed [event][pulse].
std::vector<std::vector<double>> vPulseLe; // Raw leading edge in ns (0.01 conversion only; no timing corrections).
std::vector<std::vector<double>> vPulseLeCorr; //calibrated leading edge
std::vector<std::vector<double>> vPulseToT;
std::vector<std::vector<double>> vPulseX; //uncalibrated x position
std::vector<std::vector<double>> vPulseXCorr; //calibrated x position
std::vector<std::vector<double>> vPulseY; //calibrated y position
std::vector<std::vector<double>> vPulseZ; //calibrated z position
std::vector<std::vector<bool>> vPulseGoodCalib;

//------- residuals ---------
std::vector<std::vector<double>> vPulseXResiduals;
std::vector<std::vector<double>> vPulseYResiduals;
std::vector<std::vector<double>> vPulseDt; //Ecal t - Cdet t

// Add each new vector to ClearPulseStudyVectors() and the reserve/fill sections.

void ClearPulseStudyVectors()
{
  // ECal
  vECalE.clear();
  vECalAdcTime.clear();
  vECalX.clear();
  vECalY.clear();
  vECalXProjection.clear();
  vECalYProjection.clear();

  // CDet pulses: clear() also removes the inner pulse vectors.
  vPulseID.clear();
  vPulseLe.clear();
  vPulseLeCorr.clear();
  vPulseToT.clear();
  vPulseX.clear();
  vPulseXCorr.clear();
  vPulseY.clear();
  vPulseZ.clear();
  vPulseGoodCalib.clear();

  // Residuals
  vPulseXResiduals.clear();
  vPulseYResiduals.clear();
  vPulseDt.clear();
}

void Plot_CDet_PulseStudy(
    const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
    Long64_t eventsOverride = -2, const char* inputDirectory = nullptr,
    double LeMin = NAN, double LeMax = NAN, double TotMin = NAN, double TotMax = NAN,
    int ApplyToTCut = -1, double ECalMinE = NAN, double ECalMaxE = NAN,
    double EllipseRadius = NAN, double ECalMinT = NAN, double ECalMaxT = NAN)
{
  // A second main-macro call replaces the previous sample.
  ClearPulseStudyVectors();

  // 1. Read the same analysis.* keys used by the previous configurations.
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0) {
    std::cerr << "[CDet pulse study] Cannot read configuration.\n";
    return;
  }
  if (config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pulse study] Expected config.version: 1.\n";
    return;
  }
  const int RunNumber = config.GetValue("analysis.run_number", 0);
  const int minSeg = config.GetValue("analysis.min_segment", -1);
  const int maxSeg = config.GetValue("analysis.max_segment", -1);
  // -2 uses the config; -1 reads all selected entries; 0 reads none.
  const Long64_t nevents = eventsOverride == -2
      ? config.GetValue("analysis.events", -1) : eventsOverride;
  // Same cut keys/defaults as Plot_CDet_PulsePixelTiming.C.
  // An explicit argument wins; otherwise read the value from the config.
  if (std::isnan(LeMin)) LeMin = config.GetValue("analysis.le_min", 0.02);
  if (std::isnan(LeMax)) LeMax = config.GetValue("analysis.le_max", 60.0);
  if (std::isnan(TotMin)) TotMin = config.GetValue("analysis.tot_min", 4.0);
  if (std::isnan(TotMax)) TotMax = config.GetValue("analysis.tot_max", 30.0);
  if (ApplyToTCut == -1) ApplyToTCut = config.GetValue("pulse_study.apply_tot_cut", 0) != 0;
  if (std::isnan(ECalMinE)) ECalMinE = config.GetValue("analysis.ecal_energy_min", NAN);
  if (std::isnan(ECalMaxE)) ECalMaxE = config.GetValue("analysis.ecal_energy_max", NAN);
  if (std::isnan(ECalMinT)) ECalMinT = config.GetValue("analysis.ecal_time_min", NAN);
  if (std::isnan(ECalMaxT)) ECalMaxT = config.GetValue("analysis.ecal_time_max", NAN);
  if (std::isnan(EllipseRadius)) EllipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool applyToTCut = ApplyToTCut != 0;
  const bool useEllipse = EllipseRadius > 0.0;
  // Reject events outside the configured ECal ADC-time window before the hit loop.
  // The per-pulse replay eligibility flag still supplies its cluster/position checks.
  // Optional tighter ADC-time cuts belong to the secondary plotting routines.
  // Display ranges and fit bounds in the config are NOT selection cuts.

  if (RunNumber <= 0 || nevents < -1 || minSeg < -1 || maxSeg < -1 ||
      !std::isfinite(ECalMinE) || !std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE ||
      !std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT) {
    std::cerr << "[CDet pulse study] Invalid run, event/segment limit or ECal energy/time bounds.\n";
    return;
  }
  if ((ApplyToTCut != 0 && ApplyToTCut != 1) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || !std::isfinite(TotMax) ||
       TotMin < 0.0 || TotMin >= TotMax))) {
    std::cerr << "[CDet pulse study] Invalid LE/ToT bounds or ToT switch.\n";
    return;
  }
  if (!std::isfinite(EllipseRadius) || (!useEllipse && EllipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse study] Invalid ellipse settings; radius -1 disables it.\n";
    return;
  }

  // 2. Chain the run's selected segment groups and their ROOT rollover files.
  const char* directory = inputDirectory && inputDirectory[0]
      ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0]) {
    std::cerr << "[CDet pulse study] Set OUT_DIR or pass inputDirectory.\n";
    return;
  }
  TChain T("T");
  if (CDetRunDataset::AddToChain(&T, RunNumber, directory, minSeg, maxSeg) <= 0) return;
  const Long64_t entries = T.GetEntries();
  if (entries <= 0) {
    std::cerr << "[CDet pulse study] No input events found.\n";
    return;
  }
  const Long64_t limit = nevents >= 0 ? std::min(nevents, entries) : entries;

  // 3. SET YOUR TTREEREADERS HERE, before the event loop.
  TTreeReader reader(&T);
  //-------ECal readers--------
  TTreeReaderValue<double> ECalE(reader, "earm.ecal.e");
  TTreeReaderValue<double> ECalAdcTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<double> ECalX(reader, "earm.ecal.x");
  TTreeReaderValue<double> ECalY(reader, "earm.ecal.y");
  TTreeReaderArray<double> ECalXProjection(reader, "earm.cdet.pulse.ecal_x_proj");
  TTreeReaderArray<double> ECalYProjection(reader, "earm.cdet.pulse.ecal_y_proj");
  //------CDet pulse readers--------
  TTreeReaderArray<double> PulseID(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<double> PulseLe(reader, "earm.cdet.pulse.tdc_le");
  TTreeReaderArray<double> PulseLeCorr(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<double> PulseToT(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<double> PulseX(reader, "earm.cdet.pulse.x");
  TTreeReaderArray<double> PulseXCorr(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<double> PulseY(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<double> PulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderArray<double> PulseGoodCalib(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<double> PulseXResiduals(reader, "earm.cdet.pulse.ecal_x_residual");
  TTreeReaderArray<double> PulseYResiduals(reader, "earm.cdet.pulse.ecal_y_residual");
  TTreeReaderArray<double> PulseDt(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<double> PulseECalEligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<double> PulseSpatialPass(reader, "earm.cdet.pulse.spatial_pass");

  vECalE.reserve(limit);
  vECalAdcTime.reserve(limit);
  vECalX.reserve(limit);
  vECalY.reserve(limit);
  vECalXProjection.reserve(limit);
  vECalYProjection.reserve(limit);
  vPulseID.reserve(limit);
  vPulseLe.reserve(limit);
  vPulseLeCorr.reserve(limit);
  vPulseToT.reserve(limit);
  vPulseX.reserve(limit);
  vPulseXCorr.reserve(limit);
  vPulseY.reserve(limit);
  vPulseZ.reserve(limit);
  vPulseGoodCalib.reserve(limit);
  vPulseXResiduals.reserve(limit);
  vPulseYResiduals.reserve(limit);
  vPulseDt.reserve(limit);

  std::cout << "[CDet pulse study] Run " << RunNumber << ": reading " << limit << " events.\n"
            << "  Configuration: " << configFile << '\n'
            << "  ECal energy: " << ECalMinE << " < E < " << ECalMaxE << " GeV\n"
            << "  ECal ADC time: " << ECalMinT << " < t < " << ECalMaxT << " ns (before hit loop).\n"
            << "  LE: " << LeMin << " <= LE <= " << LeMax << " ns\n"
            << "  Loading ToT window: " << TotMin << " < ToT < " << TotMax << " ns.\n"
            << "  Ellipse radius: " << EllipseRadius << " (-1 = off).\n"
            << "  Loading requires ECal energy/time, replay eligibility/spatial/calibration flags, and ToT.\n";

  // 4. ONE EVENT LOOP. Save values here; do histogram filling later.
  Long64_t EventCounter = 0;
  while (EventCounter < limit && reader.Next()) {
    ++EventCounter;
    if (EventCounter % 1000 == 0)
      std::cout << "[CDet pulse study] " << EventCounter << '/' << limit << " events\n";

    // ADD YOUR EVENT-LEVEL CUT BOOLEANS HERE.
    //ECal bools
    bool good_ECal_energy = std::isfinite(*ECalE) && *ECalE > ECalMinE && *ECalE < ECalMaxE;
    bool good_ECal_time = std::isfinite(*ECalAdcTime) && *ECalAdcTime > ECalMinT && *ECalAdcTime < ECalMaxT;
    bool good_ECal_event = good_ECal_energy && good_ECal_time;
    if (!good_ECal_event) continue;

    // Fresh pulse vectors for this event. Projections are per pulse too.
    std::vector<double> thisEventECalXProjection;
    std::vector<double> thisEventECalYProjection;
    std::vector<double> thisEventPulseID;
    std::vector<double> thisEventPulseLe;
    std::vector<double> thisEventPulseLeCorr;
    std::vector<double> thisEventPulseToT;
    std::vector<double> thisEventPulseX;
    std::vector<double> thisEventPulseXCorr;
    std::vector<double> thisEventPulseY;
    std::vector<double> thisEventPulseZ;
    std::vector<bool> thisEventPulseGoodCalib;
    std::vector<double> thisEventPulseXResiduals;
    std::vector<double> thisEventPulseYResiduals;
    std::vector<double> thisEventPulseDt;

    for (size_t hit = 0; hit < PulseID.GetSize(); ++hit) {
      // ADD YOUR CDET HIT-CUT BOOLEANS HERE, then combine them in good_hit.
      bool good_ecal_eligible = PulseECalEligible[hit] > 0.5;
      bool good_spatial_pass = PulseSpatialPass[hit] > 0.5;
      bool good_calibration = PulseGoodCalib[hit] > 0.5;
      bool good_pulse_tot = std::isfinite(PulseToT[hit]) && PulseToT[hit] > TotMin && PulseToT[hit] < TotMax;
      const bool good_hit = good_ecal_eligible && good_spatial_pass && good_calibration && good_pulse_tot;
      if (!good_hit) continue;

      // Copy all fields for the SAME passing pulse; convert raw LE to ns.
      thisEventECalXProjection.push_back(ECalXProjection[hit]);
      thisEventECalYProjection.push_back(ECalYProjection[hit]);
      thisEventPulseID.push_back(PulseID[hit]);
      thisEventPulseLe.push_back(0.01 * PulseLe[hit]); // Convert raw TDC units to ns once, when saving.
      thisEventPulseLeCorr.push_back(PulseLeCorr[hit]);
      thisEventPulseToT.push_back(PulseToT[hit]);
      thisEventPulseX.push_back(PulseX[hit]);
      thisEventPulseXCorr.push_back(PulseXCorr[hit]);
      thisEventPulseY.push_back(PulseY[hit]);
      thisEventPulseZ.push_back(PulseZ[hit]);
      thisEventPulseGoodCalib.push_back(PulseGoodCalib[hit] > 0.5);
      thisEventPulseXResiduals.push_back(PulseXResiduals[hit]);
      thisEventPulseYResiduals.push_back(PulseYResiduals[hit]);
      thisEventPulseDt.push_back(PulseDt[hit]);
    }
    if (thisEventPulseID.empty()) continue; // Skip events with no selected pulses.

    // Append once AFTER the hit loop
    vECalE.push_back(*ECalE);
    vECalAdcTime.push_back(*ECalAdcTime);
    vECalX.push_back(*ECalX);
    vECalY.push_back(*ECalY);
    vECalXProjection.push_back(thisEventECalXProjection);
    vECalYProjection.push_back(thisEventECalYProjection);
    vPulseID.push_back(thisEventPulseID);
    vPulseLe.push_back(thisEventPulseLe);
    vPulseLeCorr.push_back(thisEventPulseLeCorr);
    vPulseToT.push_back(thisEventPulseToT);
    vPulseX.push_back(thisEventPulseX);
    vPulseXCorr.push_back(thisEventPulseXCorr);
    vPulseY.push_back(thisEventPulseY);
    vPulseZ.push_back(thisEventPulseZ);
    vPulseGoodCalib.push_back(thisEventPulseGoodCalib);
    vPulseXResiduals.push_back(thisEventPulseXResiduals);
    vPulseYResiduals.push_back(thisEventPulseYResiduals);
    vPulseDt.push_back(thisEventPulseDt);
  }
  if (EventCounter != limit) {
    std::cerr << "[CDet pulse study] Reader stopped after " << EventCounter << '/' << limit
              << " events; status " << reader.GetEntryStatus()
              << ". Saved vectors contain only ECal-selected events processed so far.\n";
    return;
  }
  std::cout << "[CDet pulse study] Read " << EventCounter << " events; saved "
            << vECalE.size() << " ECal-selected events.\n";
  // Call your plotting routine from the ROOT prompt after this returns.
}

void Plot_CDet_PulseStudy(
    int RunNumber, Long64_t nevents = -2, const char* inputDirectory = nullptr,
    double LeMin = NAN, double LeMax = NAN, double TotMin = NAN, double TotMax = NAN,
    int ApplyToTCut = -1, double ECalMinE = NAN, double ECalMaxE = NAN,
    double EllipseRadius = NAN, double ECalMinT = NAN, double ECalMaxT = NAN)
{
  const TString configFile = TString::Format("CDet_run%d_pulse_pixel_timing.conf", RunNumber);
  Plot_CDet_PulseStudy(configFile.Data(), nevents, inputDirectory,
                      LeMin, LeMax, TotMin, TotMax, ApplyToTCut, ECalMinE, ECalMaxE,
                      EllipseRadius, ECalMinT, ECalMaxT);
}

// 5. Plot stored branch dt and raw LE versus ToT: summed bar and 4x4 pixels for each.
// Bar IDs are global and zero-based. Display bounds do not impose hit cuts.
// Each call rereads the config. Cuts can only select from the already-saved sample.
// These are config cuts, not the replay eligibility/spatial flags used by PulsePixelTiming.
void plotPulseStudy(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                    const char* outputTag = nullptr,
                    double BinWidth = NAN, double DtMin = NAN, double DtMax = NAN,
                    int firstBar = -1, int lastBar = -1,
                    double LeBinWidth = NAN, double LeDisplayMin = NAN, double LeDisplayMax = NAN,
                    double TotBinWidth = NAN, double TotDisplayMin = NAN, double TotDisplayMax = NAN)
{
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pulse study] Cannot read configuration with config.version: 1.\n";
    return;
  }
  if (std::isnan(BinWidth)) BinWidth = config.GetValue("pulse_study.dt_bin_width_ns", 1.0);
  if (std::isnan(DtMin)) DtMin = config.GetValue("display.ecal_cdet_dt_min", -150.0);
  if (std::isnan(DtMax)) DtMax = config.GetValue("display.ecal_cdet_dt_max", 150.0);
  if (firstBar == -1) firstBar = config.GetValue("pulse_study.first_bar", 0);
  if (lastBar == -1) lastBar = config.GetValue("pulse_study.last_bar", numBars - 1);
  if (std::isnan(LeBinWidth)) LeBinWidth = config.GetValue("pulse_study.le_bin_width_ns", 1.0);
  if (std::isnan(LeDisplayMin)) LeDisplayMin = config.GetValue("display.le_min", 0.0);
  if (std::isnan(LeDisplayMax)) LeDisplayMax = config.GetValue("display.le_max", 60.0);
  if (std::isnan(TotBinWidth)) TotBinWidth = config.GetValue("pulse_study.tot_bin_width_ns", 1.0);
  if (std::isnan(TotDisplayMin)) TotDisplayMin = config.GetValue("display.tot_min", 0.0);
  if (std::isnan(TotDisplayMax)) TotDisplayMax = config.GetValue("display.tot_max", 60.0);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 0) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 0) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if ((applyECalEnergyCut && (!std::isfinite(ECalMinE) || !std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || !std::isfinite(TotMax) || TotMin < 0 || TotMin >= TotMax))) {
    std::cerr << "[CDet pulse study] Invalid configured ECal, LE, or ToT cuts.\n";
    return;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse study] Invalid ellipse settings; radius -1 disables it.\n";
    return;
  }
  if (!std::isfinite(BinWidth) || BinWidth <= 0.0 || !std::isfinite(DtMin) ||
      !std::isfinite(DtMax) || DtMin >= DtMax || firstBar < 0 ||
      lastBar >= numBars || firstBar > lastBar ||
      !std::isfinite((DtMax - DtMin) / BinWidth) ||
      std::ceil((DtMax - DtMin) / BinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse study] Invalid dt binning or bar range (0..167).\n";
    return;
  }

  if (!std::isfinite(LeBinWidth) || LeBinWidth <= 0 || !std::isfinite(LeDisplayMin) ||
      !std::isfinite(LeDisplayMax) || LeDisplayMin >= LeDisplayMax ||
      !std::isfinite(TotBinWidth) || TotBinWidth <= 0 || !std::isfinite(TotDisplayMin) ||
      !std::isfinite(TotDisplayMax) || TotDisplayMin >= TotDisplayMax ||
      !std::isfinite((LeDisplayMax - LeDisplayMin) / LeBinWidth) ||
      !std::isfinite((TotDisplayMax - TotDisplayMin) / TotBinWidth) ||
      std::ceil((LeDisplayMax - LeDisplayMin) / LeBinWidth) > std::numeric_limits<int>::max() - 2 ||
      std::ceil((TotDisplayMax - TotDisplayMin) / TotBinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse study] Invalid LE/ToT display settings.\n";
    return;
  }
  TString tag = outputTag ? outputTag : config.GetValue("pulse_study.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d", config.GetValue("analysis.run_number", 0));
  const TString allowedTagCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowedTagCharacters.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet pulse study] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/PusleStudy/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet pulse study] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }
  const int nBins = int(std::ceil((DtMax - DtMin) / BinWidth));
  const int nLeBins = int(std::ceil((LeDisplayMax - LeDisplayMin) / LeBinWidth));
  const int nTotBins = int(std::ceil((TotDisplayMax - TotDisplayMin) / TotBinWidth));
  if (vPulseID.empty()) {
    std::cerr << "[CDet pulse study] No saved pulse events; run the main macro first.\n";
    return;
  }

  // Only this routine is batch-only; restore the caller's ROOT mode afterward.
  const bool wasBatch = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  // Replace this routine's previous canvases/histograms when redrawing.
  for (int bar = 0; bar < numBars; ++bar) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cStudyBarDt_bar%d", bar));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cStudyPixelDt_bar%d", bar));
    delete hStudyBarDt[bar];
    hStudyBarDt[bar] = nullptr;
    delete hStudyBarLeVsToT[bar];
    hStudyBarLeVsToT[bar] = nullptr;
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      delete hStudyPixelDt[bar][pixel];
      hStudyPixelDt[bar][pixel] = nullptr;
      delete hStudyPixelLeVsToT[bar][pixel];
      hStudyPixelLeVsToT[bar][pixel] = nullptr;
    }
  }

  // Define the bar and pixel histograms together before filling.
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    hStudyBarDt[bar] = new TH1D(TString::Format("hStudyBarDt_bar%d", bar), TString::Format("Bar %d: branch dt, all 16 pixels;Stored pulse.ecal_residual (ns);Pulses", bar), nBins, DtMin, DtMax);
    hStudyBarDt[bar]->SetDirectory(nullptr);
    hStudyBarLeVsToT[bar] = new TH2D(TString::Format("hStudyBarLeVsToT_bar%d", bar), TString::Format("Bar %d: raw LE versus ToT, all 16 pixels;ToT (ns);Raw CDet LE (ns);Pulses", bar), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
    hStudyBarLeVsToT[bar]->SetDirectory(nullptr);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      const int pixelID = bar * pixelsPerBar + pixel;
      hStudyPixelDt[bar][pixel] = new TH1D(TString::Format("hStudyPixelDt_pixel%d", pixelID), TString::Format("Bar %d, pixel %d (ID %d): branch dt;Stored pulse.ecal_residual (ns);Pulses", bar, pixel, pixelID), nBins, DtMin, DtMax);
      hStudyPixelDt[bar][pixel]->SetDirectory(nullptr);
      hStudyPixelLeVsToT[bar][pixel] = new TH2D(TString::Format("hStudyPixelLeVsToT_pixel%d", pixelID), TString::Format("Bar %d, pixel %d (ID %d): raw LE versus ToT;ToT (ns);Raw CDet LE (ns);Pulses", bar, pixel, pixelID), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
      hStudyPixelLeVsToT[bar][pixel]->SetDirectory(nullptr);
    }
  }

  // One pass through the aligned saved vectors. dt is already in ns.
  Long64_t filledPulses = 0;
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    const bool good_ECal_energy = std::isfinite(vECalE[event]) &&
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && vECalE[event] < ECalMaxE));
    const bool good_ECal_time = std::isfinite(vECalAdcTime[event]) &&
        (!applyECalTimeCut || (vECalAdcTime[event] > ECalMinT && vECalAdcTime[event] < ECalMaxT));
    if (!good_ECal_energy || !good_ECal_time) continue;
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit];
      const double dt = vPulseDt[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars * pixelsPerBar ||
          std::floor(id) != id || !std::isfinite(dt)) continue;
      const int bar = int(id) / pixelsPerBar;
      const int pixel = int(id) % pixelsPerBar;
      if (bar < firstBar || bar > lastBar) continue;
      // All saved timing vectors used here are in ns.
      const double le = vPulseLe[event][hit];
      const double totNs = vPulseToT[event][hit];
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(totNs) && totNs > 0.0 &&
          (!applyToTCut || (totNs >= TotMin && totNs <= TotMax));
      bool good_ellipse = true;
      if (useEllipse) {
        // Select once using the branch dt, so both plots contain identical pulses.
        const double rx = vPulseXResiduals[event][hit];
        const double normalizedX = (rx - ellipseXCenter) / ellipseXScale;
        const double normalizedTime = (dt - ellipseTimeCenter) / ellipseTimeScale;
        good_ellipse = std::isfinite(rx) &&
            normalizedX * normalizedX + normalizedTime * normalizedTime <= ellipseRadius * ellipseRadius;
      }
      if (!good_le || !good_tot || !good_ellipse) continue;
      hStudyBarDt[bar]->Fill(dt);
      hStudyPixelDt[bar][pixel]->Fill(dt);
      hStudyBarLeVsToT[bar]->Fill(totNs, le);
      hStudyPixelLeVsToT[bar][pixel]->Fill(totNs, le);
      ++filledPulses;
    }
  }

  // Save four canvases per bar; pixel pads run 0..15 in row-major order.
  // Delete each canvas after saving so no windows remain when batch mode ends.
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    TCanvas* cBar = new TCanvas(TString::Format("cStudyBarDt_bar%d", bar), TString::Format("CDet bar %d: summed dt", bar), 900, 650);
    hStudyBarDt[bar]->Draw("HIST");
    cBar->Modified();
    cBar->Update();
    cBar->SaveAs(TString::Format("%s/bar%03d_dt.pdf", outputDirectory.Data(), bar));
    delete cBar;

    TCanvas* cPixels = new TCanvas(TString::Format("cStudyPixelDt_bar%d", bar), TString::Format("CDet bar %d: pixel dt", bar), 1400, 1000);
    cPixels->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixels->cd(pixel + 1);
      hStudyPixelDt[bar][pixel]->Draw("HIST");
    }
    cPixels->Modified();
    cPixels->Update();
    cPixels->SaveAs(TString::Format("%s/bar%03d_pixels_dt.pdf", outputDirectory.Data(), bar));
    delete cPixels;

    TCanvas* cLeToTBar = new TCanvas(TString::Format("cStudyBarLeVsToT_bar%d", bar), TString::Format("CDet bar %d: summed raw LE versus ToT", bar), 900, 650);
    hStudyBarLeVsToT[bar]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hStudyBarLeVsToT[bar]->Draw("COLZ");
    cLeToTBar->SaveAs(TString::Format("%s/bar%03d_le_vs_tot.pdf", outputDirectory.Data(), bar));
    delete cLeToTBar;

    TCanvas* cLeToTPixels = new TCanvas(TString::Format("cStudyPixelLeVsToT_bar%d", bar), TString::Format("CDet bar %d: pixel raw LE versus ToT", bar), 1400, 1000);
    cLeToTPixels->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cLeToTPixels->cd(pixel + 1);
      hStudyPixelLeVsToT[bar][pixel]->SetStats(false);
      gPad->SetRightMargin(0.15);
      hStudyPixelLeVsToT[bar][pixel]->Draw("COLZ");
    }
    cLeToTPixels->SaveAs(TString::Format("%s/bar%03d_pixels_le_vs_tot.pdf", outputDirectory.Data(), bar));
    delete cLeToTPixels;
  }
  gROOT->SetBatch(wasBatch);
  std::cout << "[CDet pulse study] Saved " << 4 * (lastBar - firstBar + 1)
            << " PDFs for bars " << firstBar << ".." << lastBar
            << "; filled " << filledPulses << " pulses passing configured cuts.\n"
            << "  Secondary ECal energy cut: " << (applyECalEnergyCut ? "on" : "off") << '\n'
            << "  Secondary ECal ADC-time cut: " << (applyECalTimeCut ? "on" : "off") << '\n'
            << "  Configuration: " << configFile << '\n'
            << "  PDF output: " << outputDirectory << '\n'
            << "  dt uses vPulseDt; LE versus ToT uses uncorrected vPulseLe in ns. Both use identical pulses.\n"
            << "  Out-of-range values remain in histogram underflow/overflow bins.\n";
}

// 6. Interactive plots for the study/control bars: one spectrum or 4x4 grid per canvas.
// The optional ellipse uses branch dt; otherwise the configured x-residual cut applies.
// dt uses vPulseDt. Raw LE is already saved in ns and has no extra timing corrections.
// Every dt plot uses the stored pulse.ecal_residual branch, including dt versus rx.
// Leaves canvases open and does not change ROOT's batch setting. Use an interactive session.
void plotPulseStudyBars(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                        const char* outputTag = nullptr, int studyBar = -1, int controlBar = -1,
                        double BinWidth = NAN, double DtMin = NAN, double DtMax = NAN,
                        double LeBinWidth = NAN, double LeDisplayMin = NAN, double LeDisplayMax = NAN,
                        double TotBinWidth = NAN, double TotDisplayMin = NAN, double TotDisplayMax = NAN)
{
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pulse bar summary] Cannot read configuration with config.version: 1.\n";
    return;
  }
  if (std::isnan(BinWidth)) BinWidth = config.GetValue("pulse_study.dt_bin_width_ns", 1.0);
  if (std::isnan(DtMin)) DtMin = config.GetValue("display.ecal_cdet_dt_min", -150.0);
  if (std::isnan(DtMax)) DtMax = config.GetValue("display.ecal_cdet_dt_max", 150.0);
  if (studyBar == -1) studyBar = config.GetValue("pulse_pixel_timing.study_bar", 19);
  if (controlBar == -1) controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  if (std::isnan(LeBinWidth)) LeBinWidth = config.GetValue("pulse_study.le_bin_width_ns", 1.0);
  if (std::isnan(LeDisplayMin)) LeDisplayMin = config.GetValue("display.le_min", 0.0);
  if (std::isnan(LeDisplayMax)) LeDisplayMax = config.GetValue("display.le_max", 60.0);
  if (std::isnan(TotBinWidth)) TotBinWidth = config.GetValue("pulse_study.tot_bin_width_ns", 1.0);
  if (std::isnan(TotDisplayMin)) TotDisplayMin = config.GetValue("display.tot_min", 0.0);
  if (std::isnan(TotDisplayMax)) TotDisplayMax = config.GetValue("display.tot_max", 60.0);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 0) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 0) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if ((applyECalEnergyCut && (!std::isfinite(ECalMinE) || !std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || !std::isfinite(TotMax) || TotMin < 0 || TotMin >= TotMax))) {
    std::cerr << "[CDet pulse bar summary] Invalid configured ECal, LE, or ToT cuts.\n";
    return;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse bar summary] Invalid ellipse settings; radius -1 disables it.\n";
    return;
  }
  if (studyBar < 0 || studyBar >= numBars || controlBar < 0 || controlBar >= numBars ||
      !std::isfinite(BinWidth) || BinWidth <= 0 || !std::isfinite(DtMin) ||
      !std::isfinite(DtMax) || DtMin >= DtMax ||
      !std::isfinite((DtMax - DtMin) / BinWidth) ||
      std::ceil((DtMax - DtMin) / BinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse bar summary] Invalid bar IDs (0..167) or dt display settings.\n";
    return;
  }
  if (vPulseID.empty()) {
    std::cerr << "[CDet pulse bar summary] No saved pulse events; run the main macro first.\n";
    return;
  }
  const int bars[2] = {studyBar, controlBar};
  const char* roles[2] = {"Study", "Control"};
  if (!std::isfinite(LeBinWidth) || LeBinWidth <= 0 || !std::isfinite(LeDisplayMin) ||
      !std::isfinite(LeDisplayMax) || LeDisplayMin >= LeDisplayMax ||
      !std::isfinite(TotBinWidth) || TotBinWidth <= 0 || !std::isfinite(TotDisplayMin) ||
      !std::isfinite(TotDisplayMax) || TotDisplayMin >= TotDisplayMax ||
      !std::isfinite((LeDisplayMax - LeDisplayMin) / LeBinWidth) ||
      !std::isfinite((TotDisplayMax - TotDisplayMin) / TotBinWidth) ||
      std::ceil((LeDisplayMax - LeDisplayMin) / LeBinWidth) > std::numeric_limits<int>::max() - 2 ||
      std::ceil((TotDisplayMax - TotDisplayMin) / TotBinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse study] Invalid LE/ToT display settings.\n";
    return;
  }
  const int nDtBins = int(std::ceil((DtMax - DtMin) / BinWidth));
  const int nLeBins = int(std::ceil((LeDisplayMax - LeDisplayMin) / LeBinWidth));
  const int nTotBins = int(std::ceil((TotDisplayMax - TotDisplayMin) / TotBinWidth));
  TString tag = outputTag ? outputTag : config.GetValue("pulse_study.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d", config.GetValue("analysis.run_number", 0));
  const TString allowedTagCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowedTagCharacters.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet pulse bar summary] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/PusleStudy/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet pulse bar summary] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }

  // Replace only this routine's previous canvases, including when the bars change.
  delete gROOT->GetListOfCanvases()->FindObject("cControlBar");
  for (int selected = 0; selected < 2; ++selected) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%s", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarDt", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarLe", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelDt", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelLe", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarLeVsToT", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelLeVsToT", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarDtVsRx", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelDtVsRx", roles[selected]));
    delete hPulseSummaryBarDt[selected];
    delete hPulseSummaryBarLe[selected];
    delete hPulseSummaryBarLeVsToT[selected];
    delete hPulseBarDtVsRx[selected];
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      delete hPulseSummaryPixelDt[selected][pixel];
      delete hPulseSummaryPixelLe[selected][pixel];
      delete hPulseSummaryPixelLeVsToT[selected][pixel];
      delete hPulsePixelDtVsRx[selected][pixel];
    }
  }

  // Define all histograms before filling. Display limits are not additional cuts.
  for (int selected = 0; selected < 2; ++selected) {
    const int bar = bars[selected];
    hPulseSummaryBarDt[selected] = new TH1D(TString::Format("hPulseSummary%sBarDt", roles[selected]), TString::Format("%s bar %d: branch dt, all pixels;Stored pulse.ecal_residual (ns);Pulses", roles[selected], bar), nDtBins, DtMin, DtMax);
    hPulseSummaryBarLe[selected] = new TH1D(TString::Format("hPulseSummary%sBarLe", roles[selected]), TString::Format("%s bar %d: raw LE, all pixels;Raw CDet LE (ns);Pulses", roles[selected], bar), nLeBins, LeDisplayMin, LeDisplayMax);
    hPulseSummaryBarDt[selected]->SetDirectory(nullptr);
    hPulseSummaryBarLe[selected]->SetDirectory(nullptr);
    hPulseSummaryBarLeVsToT[selected] = new TH2D(TString::Format("hPulseSummary%sBarLeVsToT", roles[selected]), TString::Format("%s bar %d: raw LE versus ToT, all pixels;ToT (ns);Raw CDet LE (ns);Pulses", roles[selected], bar), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
    hPulseSummaryBarLeVsToT[selected]->SetDirectory(nullptr);
    hPulseBarDtVsRx[selected] = new TH2D(TString::Format("hPulseSummary%sBarDtVsRx", roles[selected]), TString::Format("%s bar %d: branch dt versus rx, all pixels;rx (m);Stored pulse.ecal_residual (ns);Pulses", roles[selected], bar), 100, -xRes, xRes, nDtBins, DtMin, DtMax);
    hPulseBarDtVsRx[selected]->SetDirectory(nullptr);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      const int id = bar * pixelsPerBar + pixel;
      hPulseSummaryPixelDt[selected][pixel] = new TH1D(TString::Format("hPulseSummary%sPixel%dDt", roles[selected], pixel), TString::Format("Pixel %d (ID %d): branch dt;Stored pulse.ecal_residual (ns);Pulses", pixel, id), nDtBins, DtMin, DtMax);
      hPulseSummaryPixelLe[selected][pixel] = new TH1D(TString::Format("hPulseSummary%sPixel%dLe", roles[selected], pixel), TString::Format("Pixel %d (ID %d): raw LE;Raw CDet LE (ns);Pulses", pixel, id), nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelDt[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLe[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLeVsToT[selected][pixel] = new TH2D(TString::Format("hPulseSummary%sPixel%dLeVsToT", roles[selected], pixel), TString::Format("Pixel %d (ID %d): raw LE versus ToT;ToT (ns);Raw CDet LE (ns);Pulses", pixel, id), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelLeVsToT[selected][pixel]->SetDirectory(nullptr);
      hPulsePixelDtVsRx[selected][pixel] = new TH2D(TString::Format("hPulseSummary%sPixel%dDtVsRx", roles[selected], pixel), TString::Format("Pixel %d (ID %d): branch dt versus rx;rx (m);Stored pulse.ecal_residual (ns);Pulses", pixel, id), 100, -xRes, xRes, nDtBins, DtMin, DtMax);
      hPulsePixelDtVsRx[selected][pixel]->SetDirectory(nullptr);
    }
  }

  // One pass through the aligned saved vectors. dt is already in ns.
  Long64_t filledPulses = 0;
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    const bool good_ECal_energy = std::isfinite(vECalE[event]) &&
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && vECalE[event] < ECalMaxE));
    const bool good_ECal_time = std::isfinite(vECalAdcTime[event]) &&
        (!applyECalTimeCut || (vECalAdcTime[event] > ECalMinT && vECalAdcTime[event] < ECalMaxT));
    if (!good_ECal_energy || !good_ECal_time) continue;
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit];
      const double dt = vPulseDt[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars * pixelsPerBar ||
          std::floor(id) != id || !std::isfinite(dt)) continue;
      const int bar = int(id) / pixelsPerBar;
      const int pixel = int(id) % pixelsPerBar;

      if (bar != studyBar && bar != controlBar) continue;
      // All saved timing vectors used here are in ns.
      const double le = vPulseLe[event][hit];
      const double totNs = vPulseToT[event][hit];
      const double rx = vPulseXResiduals[event][hit];
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(totNs) && totNs > 0.0 &&
          (!applyToTCut || (totNs >= TotMin && totNs <= TotMax));
      const bool good_x_res = std::isfinite(rx) && std::abs(rx) <= xRes;

      bool good_ellipse = true;
      if (useEllipse) {
        // Select once using the branch dt, so both plots contain identical pulses.
        const double normalizedX = (rx - ellipseXCenter) / ellipseXScale;
        const double normalizedTime = (dt - ellipseTimeCenter) / ellipseTimeScale;
        good_ellipse = std::isfinite(rx) &&
            normalizedX * normalizedX + normalizedTime * normalizedTime <= ellipseRadius * ellipseRadius;
      }
      else {
        // If the ellipse is not used, still apply the x residual cut.
        good_ellipse = good_x_res;
      }
      if (!good_le || !good_tot || !good_ellipse) continue;
      for (int selected = 0; selected < 2; ++selected) {
        if (bar != bars[selected]) continue;
        hPulseSummaryBarDt[selected]->Fill(dt);
        hPulseSummaryBarLe[selected]->Fill(le);
        hPulseSummaryPixelDt[selected][pixel]->Fill(dt);
        hPulseSummaryPixelLe[selected][pixel]->Fill(le);
        hPulseSummaryBarLeVsToT[selected]->Fill(totNs, le);
        hPulseSummaryPixelLeVsToT[selected][pixel]->Fill(totNs, le);
        hPulseBarDtVsRx[selected]->Fill(rx, dt);
        hPulsePixelDtVsRx[selected][pixel]->Fill(rx, dt);
      }
      ++filledPulses;
    }
  }
  // Eight independent canvases per bar; each pixel canvas keeps its 4x4 grid.
  for (int selected = 0; selected < 2; ++selected) {
    TCanvas* cBarDt = new TCanvas(TString::Format("cPulseSummary%sBarDt", roles[selected]), TString::Format("%s bar %d: summed branch dt", roles[selected], bars[selected]), 900, 650);
    hPulseSummaryBarDt[selected]->Draw("HIST");
    cBarDt->Modified();
    cBarDt->Update();
    cBarDt->SaveAs(TString::Format("%s/bar%03d_%s_dt.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cBarLe = new TCanvas(TString::Format("cPulseSummary%sBarLe", roles[selected]), TString::Format("%s bar %d: summed raw LE", roles[selected], bars[selected]), 900, 650);
    hPulseSummaryBarLe[selected]->Draw("HIST");
    cBarLe->Modified();
    cBarLe->Update();
    cBarLe->SaveAs(TString::Format("%s/bar%03d_%s_raw_le.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cPixelDt = new TCanvas(TString::Format("cPulseSummary%sPixelDt", roles[selected]), TString::Format("%s bar %d: pixel branch dt", roles[selected], bars[selected]), 1400, 1000);
    cPixelDt->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelDt->cd(pixel + 1);
      hPulseSummaryPixelDt[selected][pixel]->Draw("HIST");
    }
    cPixelDt->Modified();
    cPixelDt->Update();
    cPixelDt->SaveAs(TString::Format("%s/bar%03d_%s_pixels_dt.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cPixelLe = new TCanvas(TString::Format("cPulseSummary%sPixelLe", roles[selected]), TString::Format("%s bar %d: pixel raw LE", roles[selected], bars[selected]), 1400, 1000);
    cPixelLe->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelLe->cd(pixel + 1);
      hPulseSummaryPixelLe[selected][pixel]->Draw("HIST");
    }
    cPixelLe->Modified();
    cPixelLe->Update();
    cPixelLe->SaveAs(TString::Format("%s/bar%03d_%s_pixels_raw_le.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cBarLeVsToT = new TCanvas(TString::Format("cPulseSummary%sBarLeVsToT", roles[selected]), TString::Format("%s bar %d: summed raw LE versus ToT", roles[selected], bars[selected]), 900, 650);
    hPulseSummaryBarLeVsToT[selected]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hPulseSummaryBarLeVsToT[selected]->Draw("COLZ");
    cBarLeVsToT->Modified();
    cBarLeVsToT->Update();
    cBarLeVsToT->SaveAs(TString::Format("%s/bar%03d_%s_le_vs_tot.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cPixelLeVsToT = new TCanvas(TString::Format("cPulseSummary%sPixelLeVsToT", roles[selected]), TString::Format("%s bar %d: pixel raw LE versus ToT", roles[selected], bars[selected]), 1400, 1000);
    cPixelLeVsToT->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelLeVsToT->cd(pixel + 1);
      hPulseSummaryPixelLeVsToT[selected][pixel]->SetStats(false);
      gPad->SetRightMargin(0.15);
      hPulseSummaryPixelLeVsToT[selected][pixel]->Draw("COLZ");
    }
    cPixelLeVsToT->Modified();
    cPixelLeVsToT->Update();
    cPixelLeVsToT->SaveAs(TString::Format("%s/bar%03d_%s_pixels_le_vs_tot.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cBarDtVsRx = new TCanvas(TString::Format("cPulseSummary%sBarDtVsRx", roles[selected]), TString::Format("%s bar %d: summed branch dt versus rx", roles[selected], bars[selected]), 900, 650);
    hPulseBarDtVsRx[selected]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hPulseBarDtVsRx[selected]->Draw("COLZ");
    cBarDtVsRx->Modified();
    cBarDtVsRx->Update();
    cBarDtVsRx->SaveAs(TString::Format("%s/bar%03d_%s_dt_vs_rx.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cPixelDtVsRx = new TCanvas(TString::Format("cPulseSummary%sPixelDtVsRx", roles[selected]), TString::Format("%s bar %d: pixel branch dt versus rx", roles[selected], bars[selected]), 1400, 1000);
    cPixelDtVsRx->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelDtVsRx->cd(pixel + 1);
      hPulsePixelDtVsRx[selected][pixel]->SetStats(false);
      gPad->SetRightMargin(0.15);
      hPulsePixelDtVsRx[selected][pixel]->Draw("COLZ");
    }
    cPixelDtVsRx->Modified();
    cPixelDtVsRx->Update();
    cPixelDtVsRx->SaveAs(TString::Format("%s/bar%03d_%s_pixels_dt_vs_rx.pdf", outputDirectory.Data(), bars[selected], roles[selected]));
  }
  std::cout << "[CDet pulse bar summary] Drew study bar " << studyBar << " and control bar " << controlBar
            << "; selected " << filledPulses << " pulses. Sixteen canvases remain open.\n"
            << "  dt uses vPulseDt; LE versus ToT uses uncorrected vPulseLe. All times are in ns.\n"
            << "  dt versus rx uses vPulseDt against stored x residual.\n"
            << "  PDF output: " << outputDirectory << '\n';
}

// 7. Interactive x/rx overlays for the study/control bars, with the same cuts as
// plotPulseStudyBars. UseCorrectedX: -1 reads config, 0 channel x, 1 aligned x.
// Both histograms contain the same pulses; display bounds do not select pulses.
void plotPulseStudyPositions(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                            const char* outputTag = nullptr, int studyBar = -1, int controlBar = -1,
                            double PositionBinWidth = NAN, double PositionMin = NAN, double PositionMax = NAN,
                            int UseCorrectedX = -1)
{
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pulse positions] Cannot read configuration with config.version: 1.\n";
    return;
  }
  if (studyBar == -1) studyBar = config.GetValue("pulse_pixel_timing.study_bar", 19);
  if (controlBar == -1) controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  if (std::isnan(PositionBinWidth)) PositionBinWidth = config.GetValue("pulse_study.position_bin_width_m", 0.005);
  if (std::isnan(PositionMin)) PositionMin = config.GetValue("pulse_study.position_min_m", -1.6);
  if (std::isnan(PositionMax)) PositionMax = config.GetValue("pulse_study.position_max_m", 1.6);
  if (UseCorrectedX == -1) UseCorrectedX = config.GetValue("pulse_study.position_use_corrected_x", 1);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 0) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 0) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if ((applyECalEnergyCut && (!std::isfinite(ECalMinE) || !std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || !std::isfinite(TotMax) || TotMin < 0 || TotMin >= TotMax))) {
    std::cerr << "[CDet pulse positions] Invalid configured ECal, LE, or ToT cuts.\n";
    return;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse positions] Invalid ellipse settings; radius -1 disables it.\n";
    return;
  }
  if (studyBar < 0 || studyBar >= numBars || controlBar < 0 || controlBar >= numBars ||
      (UseCorrectedX != 0 && UseCorrectedX != 1) ||
      !std::isfinite(PositionBinWidth) || PositionBinWidth <= 0 ||
      !std::isfinite(PositionMin) || !std::isfinite(PositionMax) || PositionMin >= PositionMax ||
      !std::isfinite((PositionMax - PositionMin) / PositionBinWidth) ||
      std::ceil((PositionMax - PositionMin) / PositionBinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse positions] Invalid bar IDs, x-source switch, or position display settings.\n";
    return;
  }
  if (vPulseID.empty()) {
    std::cerr << "[CDet pulse positions] No saved pulse events; run the main macro first.\n";
    return;
  }
  const int nPositionBins = int(std::ceil((PositionMax - PositionMin) / PositionBinWidth));
  const int bars[2] = {studyBar, controlBar};
  const char* roles[2] = {"Study", "Control"};
  const char* xLabel = UseCorrectedX ? "Aligned CDet x (pulse.x_corr)" : "Channel CDet x (pulse.x)";
  TString tag = outputTag ? outputTag : config.GetValue("pulse_study.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d", config.GetValue("analysis.run_number", 0));
  const TString allowedTagCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowedTagCharacters.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet pulse positions] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/PusleStudy/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet pulse positions] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }

  // Replace this routine's own canvases/histograms before defining new ones.
  for (int selected = 0; selected < 2; ++selected) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulsePosition%s", roles[selected]));
    delete hPulsePositionBarX[selected];
    delete hPulsePositionBarRx[selected];
  }

  // Define all histograms before the event/hit loop, using shared position axes.
  for (int selected = 0; selected < 2; ++selected) {
    hPulsePositionBarX[selected] = new TH1D(TString::Format("hPulsePosition%sX", roles[selected]), TString::Format("%s bar %d: %s and rx;Position / residual (m);Pulses", roles[selected], bars[selected], UseCorrectedX ? "aligned x" : "channel x"), nPositionBins, PositionMin, PositionMax);
    hPulsePositionBarRx[selected] = new TH1D(TString::Format("hPulsePosition%sRx", roles[selected]), TString::Format("%s bar %d: stored rx;Position / residual (m);Pulses", roles[selected], bars[selected]), nPositionBins, PositionMin, PositionMax);
    hPulsePositionBarX[selected]->SetDirectory(nullptr);
    hPulsePositionBarRx[selected]->SetDirectory(nullptr);
  }

  // One pass through the aligned saved vectors. dt is already in ns.
  Long64_t filledPulses = 0;
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    const bool good_ECal_energy = std::isfinite(vECalE[event]) &&
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && vECalE[event] < ECalMaxE));
    const bool good_ECal_time = std::isfinite(vECalAdcTime[event]) &&
        (!applyECalTimeCut || (vECalAdcTime[event] > ECalMinT && vECalAdcTime[event] < ECalMaxT));
    if (!good_ECal_energy || !good_ECal_time) continue;
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit];
      const double dt = vPulseDt[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars * pixelsPerBar ||
          std::floor(id) != id || !std::isfinite(dt)) continue;
      const int bar = int(id) / pixelsPerBar;

      if (bar != studyBar && bar != controlBar) continue;
      // All saved timing vectors used here are in ns.
      const double le = vPulseLe[event][hit];
      const double totNs = vPulseToT[event][hit];
      const double rx = vPulseXResiduals[event][hit];
      const double x = UseCorrectedX ? vPulseXCorr[event][hit] : vPulseX[event][hit];
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(totNs) && totNs > 0.0 &&
          (!applyToTCut || (totNs >= TotMin && totNs <= TotMax));
      const bool good_x_res = std::isfinite(rx) && std::abs(rx) <= xRes;

      bool good_ellipse = true;
      if (useEllipse) {
        // Select once using the branch dt, so both plots contain identical pulses.
        const double normalizedX = (rx - ellipseXCenter) / ellipseXScale;
        const double normalizedTime = (dt - ellipseTimeCenter) / ellipseTimeScale;
        good_ellipse = std::isfinite(rx) &&
            normalizedX * normalizedX + normalizedTime * normalizedTime <= ellipseRadius * ellipseRadius;
      }
      else {
        // If the ellipse is not used, still apply the x residual cut.
        good_ellipse = good_x_res;
      }
      if (!good_le || !good_tot || !good_ellipse || !std::isfinite(x)) continue;
      for (int selected = 0; selected < 2; ++selected) {
        if (bar != bars[selected]) continue;
        hPulsePositionBarX[selected]->Fill(x);
        hPulsePositionBarRx[selected]->Fill(rx);
      }
      ++filledPulses;
    }
  }
  // Draw one overlay per bar after filling; retain canvases in the caller's mode.
  for (int selected = 0; selected < 2; ++selected) {
    TCanvas* canvas = new TCanvas(TString::Format("cPulsePosition%s", roles[selected]), TString::Format("%s bar %d: x and rx", roles[selected], bars[selected]), 1000, 650);
    hPulsePositionBarX[selected]->SetStats(false);
    hPulsePositionBarRx[selected]->SetStats(false);
    hPulsePositionBarX[selected]->SetLineColor(kBlue + 1);
    hPulsePositionBarRx[selected]->SetLineColor(kRed + 1);
    hPulsePositionBarX[selected]->SetLineWidth(2);
    hPulsePositionBarRx[selected]->SetLineWidth(2);
    hPulsePositionBarX[selected]->SetMinimum(0.0);
    hPulsePositionBarX[selected]->SetMaximum(1.35 * std::max(1.0, std::max(hPulsePositionBarX[selected]->GetMaximum(), hPulsePositionBarRx[selected]->GetMaximum())));
    hPulsePositionBarX[selected]->Draw("HIST");
    hPulsePositionBarRx[selected]->Draw("HIST SAME");
    TLegend* legend = new TLegend(0.52, 0.76, 0.89, 0.89);
    legend->AddEntry(hPulsePositionBarX[selected], xLabel, "l");
    legend->AddEntry(hPulsePositionBarRx[selected], "rx (pulse.ecal_x_residual)", "l");
    legend->Draw();
    legend->SetBit(TObject::kCanDelete);
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(TString::Format("%s/bar%03d_%s_x_and_rx.pdf", outputDirectory.Data(), bars[selected], roles[selected]));
  }
  std::cout << "[CDet pulse positions] Drew study bar " << studyBar << " and control bar " << controlBar
            << "; selected " << filledPulses << " pulses. Two position canvases remain open.\n"
            << "  Blue: " << xLabel << "; red: stored rx. Positions are in metres.\n"
            << "  Values outside the display bounds remain in underflow/overflow bins.\n"
            << "  PDF output: " << outputDirectory << '\n';
}

// calibrationFile must contain the pixel constants used to replay the saved events.
// Results are proposals only; this routine does not modify calibration files.
// CSV writing is opt-in. Reusing an output filename replaces its previous contents.
// CSV: all 16 pixels (ns). Successful donors retain their own fitted mu;
// other pixels inherit the weighted donor mu from the full bar or their half.
// c = mu - mu0; final offset = this pixel's existing offset + c.
// Empty fields mean unresolved. Errors exclude uncertainty in the inheritance assumption.
void extractBarOffsets(bool viewOnly = true, int myBar = 0, const std::vector<int> &goodPixels = {}, bool fullNINO = true, const char* configFile = "CDet_run6077_pulse_pixel_timing.conf", const char* calibrationFile = "CDet_calibration_dt.dat", bool writeCSV = false, const char* csvFile = nullptr)
{
  const int fitMode = fullNINO ? 2 : 1; // 2: full bar; 1: separate eight-pixel halves.

  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pulse offsets] Cannot read configuration with config.version: 1.\n";
    return;
  }

  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  if (myBar < 0 || myBar >= numBars || controlBar < 0 || controlBar >= numBars) {
    std::cerr << "[CDet pulse offsets] Bar IDs must be 0..167.\n";
    return;
  }
  for (int pixel : goodPixels) {
    if (pixel < 0 || pixel >= pixelsPerBar) {
      std::cerr << "[CDet pulse offsets] Selected local pixel IDs must be 0..15.\n";
      return;
    }
  }

  // Master-file rows are global pixel ID, offset (ns), and hit count.
  // Missing offsets must not silently become zero.
  std::vector<double> existingOffsets(numBars * pixelsPerBar, NAN);
  if (!viewOnly) {
    std::ifstream calibration(calibrationFile ? calibrationFile : "");
    if (!calibration) {
      std::cerr << "[CDet pulse offsets] Cannot read calibration file.\n";
      return;
    }
    bool inPixelOffsets = false;
    std::string line;
    while (std::getline(calibration, line)) {
      line = line.substr(0, line.find('#'));
      const size_t first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos) continue;
      line = line.substr(first, line.find_last_not_of(" \t\r") - first + 1);
      if (line[0] == '[') {
        inPixelOffsets = line == "[PixelOffsets]";
        continue;
      }
      if (!inPixelOffsets) continue;
      std::istringstream row(line);
      int pixelID;
      double offset;
      Long64_t nhits;
      if (!(row >> pixelID >> offset >> nhits) || pixelID < 0 ||
          pixelID >= numBars * pixelsPerBar || !std::isfinite(offset) ||
          std::isfinite(existingOffsets[pixelID])) {
        std::cerr << "[CDet pulse offsets] Invalid or duplicate [PixelOffsets] row: " << line << '\n';
        return;
      }
      existingOffsets[pixelID] = offset;
    }
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      if (!std::isfinite(existingOffsets[myBar * pixelsPerBar + pixel])) {
        std::cerr << "[CDet pulse offsets] Missing existing offset for pixel "
                  << myBar * pixelsPerBar + pixel << ".\n";
        return;
      }
    }
    std::cout << "[CDet pulse offsets] Existing constants: " << calibrationFile
              << " (must match the ROOT-file replay calibration).\n";
  }

  // define potential cuts
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 0) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double controlFitMin = config.GetValue("pulse_study.control_fit_min", -40.0);
  const double controlFitMax = config.GetValue("pulse_study.control_fit_max", 0.0);
  const double FitMin = config.GetValue("pulse_study.fit_min", -50.0);
  const double FitMax = config.GetValue("pulse_study.fit_max", 0.0);
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 0) != 0;
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08);

  if (!viewOnly && (!std::isfinite(controlFitMin) || !std::isfinite(controlFitMax) || controlFitMin >= controlFitMax ||
      !std::isfinite(FitMin) || !std::isfinite(FitMax) || FitMin >= FitMax)) {
    std::cerr << "[CDet pulse offsets] Invalid fit windows.\n";
    return;
  }
  if (vPulseID.empty()) {
    std::cerr << "[CDet pulse offsets] No saved events; run the main macro first.\n";
    return;
  }

  // Always define all 16 spectra so empty/unselected pixels have their own pad.
  TH1D* hBarDt[pixelsPerBar] = {};
  TH1D* hControlBarDt[pixelsPerBar] = {};
  for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
    hBarDt[pixel] = new TH1D(TString::Format("hBar%dPixel%dDt", myBar, pixel), TString::Format("Bar %d pixel %d;dt [ns];Pulses", myBar, pixel), 100, -50, 50);
    hBarDt[pixel]->SetDirectory(nullptr);
    hControlBarDt[pixel] = new TH1D(TString::Format("hControlBar%dPixel%dDt", controlBar, pixel), TString::Format("Control bar %d pixel %d;dt [ns];Pulses", controlBar, pixel), 100, -50, 50);
    hControlBarDt[pixel]->SetDirectory(nullptr);
  }

  Int_t nevents = vECalAdcTime.size();
  for (Int_t ev = 0; ev < nevents; ev++){
    Int_t nhits = vPulseID[ev].size();
    for (Int_t hit = 0; hit < nhits; hit++){
      Int_t bar = Int_t(vPulseID[ev][hit])/pixelsPerBar; //determine which bar
      Int_t pixel = Int_t(vPulseID[ev][hit])%pixelsPerBar; //convert pixelID into pixel within bar
      double rx = vPulseXResiduals[ev][hit];
      if (!std::isfinite(rx) || std::abs(rx) > xRes || !std::isfinite(vPulseDt[ev][hit])) continue; //apply x residual cut
      if (bar == myBar) hBarDt[pixel]->Fill(vPulseDt[ev][hit]);
      if (bar == controlBar) hControlBarDt[pixel]->Fill(vPulseDt[ev][hit]);
    } //end hit loop
  } //end event loop

  // Same global unused-channel IDs as PlotElastic.C (not inferred from occupancy).
  const std::vector<int> unusedPixels = {
    3, 13, 28, 31, 41, 42, 57, 59, 65, 79, 83, 95, 109, 111, 115, 127,
    140, 143, 145, 156, 172, 175, 176, 188, 195, 199, 213, 220, 236, 239, 244, 255,
    268, 271, 284, 287, 300, 303, 307, 319, 332, 335, 339, 351, 354, 364, 371, 381,
    384, 396, 401, 410, 419, 423, 435, 436, 451, 461, 465, 479, 480, 483, 508, 511,
    512, 515, 540, 543, 546, 559, 563, 573, 576, 589, 596, 605, 609, 610, 627, 638,
    643, 655, 656, 665, 674, 675, 696, 703, 707, 709, 725, 729, 738, 748, 752, 766,
    777, 780, 784, 791, 800, 812, 818, 828, 844, 847, 850, 860, 867, 868, 884, 885,
    900, 904, 912, 927, 940, 943, 945, 947, 967, 971, 986, 991, 1005, 1007, 1011, 1023,
    1027, 1028, 1043, 1050, 1066, 1068, 1072, 1075, 1088, 1102, 1106, 1119, 1121, 1135, 1148, 1151,
    1162, 1166, 1178, 1182, 1184, 1186, 1203, 1215, 1228, 1231, 1235, 1247, 1249, 1252, 1267, 1274,
    1282, 1285, 1299, 1310, 1317, 1321, 1340, 1341, 1349, 1359, 1372, 1375, 1376, 1391, 1392, 1405,
    1409, 1420, 1428, 1439, 1443, 1455, 1468, 1471, 1486, 1487, 1500, 1503, 1516, 1519, 1520, 1523,
    1536, 1551, 1557, 1567, 1568, 1583, 1584, 1597, 1603, 1615, 1617, 1629, 1646, 1647, 1648, 1654,
    1676, 1679, 1692, 1695, 1708, 1711, 1715, 1725, 1732, 1743, 1744, 1757, 1761, 1770, 1778, 1786,
    1804, 1807, 1820, 1823, 1836, 1839, 1854, 1855, 1856, 1868, 1877, 1887, 1902, 1903, 1916, 1919,
    1934, 1935, 1942, 1951, 1964, 1967, 1973, 1983, 1988, 1999, 2000, 2013, 2028, 2031, 2034, 2047,
    2048, 2051, 2064, 2067, 2080, 2085, 2099, 2104, 2112, 2122, 2131, 2143, 2144, 2147, 2160, 2163,
    2177, 2188, 2202, 2207, 2208, 2221, 2227, 2239, 2243, 2254, 2259, 2271, 2279, 2283, 2300, 2303,
    2307, 2316, 2320, 2334, 2339, 2348, 2355, 2367, 2369, 2383, 2384, 2395, 2405, 2409, 2416, 2422,
    2432, 2435, 2448, 2451, 2464, 2479, 2483, 2493, 2499, 2508, 2512, 2513, 2531, 2537, 2544, 2547,
    2563, 2570, 2576, 2591, 2592, 2607, 2611, 2621, 2633, 2636, 2643, 2650, 2656, 2657, 2675, 2679
  };
  if (gSystem->AccessPathName("PusleStudy") && gSystem->mkdir("PusleStudy", true) != 0) {
    std::cerr << "[CDet pulse offsets] Cannot create PusleStudy output directory.\n";
    return;
  }
  const int run = config.GetValue("analysis.run_number", 0);
  const char* plotMode = viewOnly ? "view" : "fit";
  const TString studyPDF = TString::Format("PusleStudy/run%d_bar%03d_dt_pixels_%s.pdf", run, myBar, plotMode);
  const TString controlPDF = TString::Format("PusleStudy/run%d_bar%03d_control%03d_dt_pixels_%s.pdf", run, myBar, controlBar, plotMode);
  if (TObject* old = gROOT->GetListOfCanvases()->FindObject("cPulseOffsetStudy")) delete old;
  if (TObject* old = gROOT->GetListOfCanvases()->FindObject("cPulseOffsetControl")) delete old;
  TCanvas* cBarDt = new TCanvas("cPulseOffsetStudy", TString::Format("Bar %d: pixel dt spectra (%s)", myBar, plotMode), 1500, 1100);
  TCanvas* cControlDt = new TCanvas("cPulseOffsetControl", TString::Format("Control bar %d: pixel dt spectra (%s)", controlBar, plotMode), 1500, 1100);
  TCanvas* canvases[2] = {cBarDt, cControlDt};
  const int plotBars[2] = {myBar, controlBar};
  for (int which = 0; which < 2; ++which) {
    canvases[which]->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      canvases[which]->cd(pixel + 1);
      TH1D* hist = which == 0 ? hBarDt[pixel] : hControlBarDt[pixel];
      hist->SetBit(kCanDelete); // Canvas owns its displayed histogram across repeat calls.
      hist->SetStats(false); // Draw a consistent statistics snapshot, including empty spectra.
      hist->SetLineColor(kBlue + 1);
      hist->SetLineWidth(2);
      hist->GetXaxis()->SetTitleSize(0.05);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->GetXaxis()->SetLabelSize(0.04);
      hist->GetYaxis()->SetLabelSize(0.04);
      gPad->SetBottomMargin(0.14);
      gPad->SetLeftMargin(0.14);
      hist->Draw();
      gPad->Update();
      TPaveStats* stats = new TPaveStats(0.63, 0.68, 0.98, 0.91, "NDC");
      stats->SetName("pixelStats");
      stats->SetBit(kCanDelete);
      stats->SetOptStat(1110);
      stats->SetOptFit(0);
      stats->SetFillColor(kWhite);
      stats->SetBorderSize(1);
      stats->SetTextFont(42);
      stats->SetTextSize(0.035);
      stats->SetTextAlign(12);
      stats->AddText(TString::Format("Entries = %.0f", hist->GetEntries()));
      stats->AddText(TString::Format("Mean = %.4g", hist->GetMean()));
      stats->AddText(TString::Format("Std Dev = %.4g", hist->GetStdDev()));
      stats->Draw("same");
      const int pixelID = plotBars[which] * pixelsPerBar + pixel;
      if (std::find(unusedPixels.begin(), unusedPixels.end(), pixelID) != unusedPixels.end()) {
        TPaveText* mask = new TPaveText(0.63, 0.68, 0.98, 0.91, "NDC");
        mask->SetName("unusedPixelMask");
        mask->SetBit(kCanDelete);
        mask->SetFillColor(kBlack);
        mask->SetFillStyle(1001);
        mask->SetLineColor(kBlack);
        mask->SetBorderSize(0);
        mask->Draw("same");
      }
      gPad->Modified();
      gPad->Update();
    }
    canvases[which]->Modified();
    canvases[which]->Update();
  }
  // Keep histogram-only PDFs even if a later control fit cannot establish mu0.
  cBarDt->SaveAs(studyPDF);
  cControlDt->SaveAs(controlPDF);
  if (viewOnly) {
    std::cout << "[CDet pulse offsets] View only: no fits, offset calculation, or CSV writing.\n";
    return;
  }

  //fit control bar histograms to get mu0
  TF1* fControlBarFit[pixelsPerBar] = {};
  double muControl[pixelsPerBar]; //values of mu for control bar
  double muControlError[pixelsPerBar]; //errors on mu for control bar
  double mu0 = 0.0;
  double mu0Error = 0.0;

  for (Int_t pixel = 0; pixel < pixelsPerBar; pixel++){
    if (hControlBarDt[pixel]->GetEntries() > 0){
      fControlBarFit[pixel] = new TF1(TString::Format("fControlBar%dPixel%dFit", controlBar, pixel), "gaus(0)+pol1(3)", controlFitMin, controlFitMax);
      hControlBarDt[pixel]->GetXaxis()->SetRangeUser(controlFitMin, controlFitMax);
      const double peak = hControlBarDt[pixel]->GetBinCenter(hControlBarDt[pixel]->GetMaximumBin());
      const double background = hControlBarDt[pixel]->GetMinimum();
      fControlBarFit[pixel]->SetParameters(std::max(1.0, hControlBarDt[pixel]->GetMaximum() - background), peak, (controlFitMax - controlFitMin) / 10.0, background, 0.0);
      hControlBarDt[pixel]->GetXaxis()->SetRange(0, 0);
      TFitResultPtr fitResult = hControlBarDt[pixel]->Fit(fControlBarFit[pixel], "QRS0");
      muControl[pixel] = fControlBarFit[pixel]->GetParameter(1);
      muControlError[pixel] = fControlBarFit[pixel]->GetParError(1);
      if (int(fitResult) != 0 || !fitResult.Get() || fitResult->CovMatrixStatus() < 2 ||
          !std::isfinite(muControl[pixel]) || !std::isfinite(muControlError[pixel]) || muControlError[pixel] <= 0.0 ||
          muControl[pixel] <= controlFitMin || muControl[pixel] >= controlFitMax ||
          fControlBarFit[pixel]->GetParameter(0) <= 0.0 ||
          !std::isfinite(fControlBarFit[pixel]->GetParameter(2)) || std::abs(fControlBarFit[pixel]->GetParameter(2)) <= 0.0) {
        std::cout << "  Skipping control pixel " << pixel << ": unusable fit.\n";
        continue;
      }
      TF1* drawnFit = hControlBarDt[pixel]->GetFunction(fControlBarFit[pixel]->GetName());
      if (drawnFit) {
        drawnFit->ResetBit(TF1::kNotDraw);
        drawnFit->SetLineColor(kRed + 1);
        drawnFit->SetLineWidth(2);
      }
      mu0 += muControl[pixel] / (muControlError[pixel] * muControlError[pixel]);
      mu0Error += 1.0 / (muControlError[pixel] * muControlError[pixel]);
    }
  }

  if (!std::isfinite(mu0Error) || mu0Error <= 0.0 || !std::isfinite(mu0)) {
    std::cerr << "[CDet pulse offsets] No usable weighted control reference; no offsets proposed.\n";
    return;
  }
  mu0 /= mu0Error;
  mu0Error = std::sqrt(1.0 / mu0Error);

  std::cout << "mu0 = " << mu0 << " +/- " << mu0Error << std::endl;
  if (!viewOnly){
    // Fit selected pixels individually, grouped for the later offset calculation.
    // Mode 2: one group (0-15). Mode 1: two groups (0-7 and 8-15).
    TF1* fBarFit[pixelsPerBar] = {};
    double muBar[pixelsPerBar]; //values of mu for bar of interest
    double muBarError[pixelsPerBar]; //errors on mu for bar of interest
    double barOffset[pixelsPerBar]; //proposed per-pixel final constants, not applied to the DB
    double barOffsetError[pixelsPerBar]; //fit/reference errors only; excludes sharing/model uncertainty
    bool hasIndividualFit[pixelsPerBar] = {};
    std::fill_n(muBar, pixelsPerBar, NAN);
    std::fill_n(muBarError, pixelsPerBar, NAN);
    std::fill_n(barOffset, pixelsPerBar, NAN);
    std::fill_n(barOffsetError, pixelsPerBar, NAN);

    const int nGroups = fitMode == 2 ? 1 : 2;
    const int pixelsPerGroup = pixelsPerBar / nGroups;
    for (Int_t group = 0; group < nGroups; group++){
      const int firstPixel = group * pixelsPerGroup;
      const int lastPixel = firstPixel + pixelsPerGroup;
      double weightedMuSum = 0.0;
      double weightSum = 0.0;
      int donors = 0;
      for (Int_t pixel = firstPixel; pixel < lastPixel; pixel++){
        // An empty list considers all pixels in this group.
        if (!goodPixels.empty() &&
            std::find(goodPixels.begin(), goodPixels.end(), pixel) == goodPixels.end()) continue;
        if (hBarDt[pixel]->GetEntries() > 0){
          fBarFit[pixel] = new TF1(TString::Format("fBar%dPixel%dFit", myBar, pixel), "gaus(0)+pol1(3)", FitMin, FitMax);
          hBarDt[pixel]->GetXaxis()->SetRangeUser(FitMin, FitMax);
          const double peak = hBarDt[pixel]->GetBinCenter(hBarDt[pixel]->GetMaximumBin());
          const double background = hBarDt[pixel]->GetMinimum();
          fBarFit[pixel]->SetParameters(std::max(1.0, hBarDt[pixel]->GetMaximum() - background), peak, (FitMax - FitMin) / 10.0, background, 0.0);
          hBarDt[pixel]->GetXaxis()->SetRange(0, 0);
          TFitResultPtr fitResult = hBarDt[pixel]->Fit(fBarFit[pixel], "QRS0");
          muBar[pixel] = fBarFit[pixel]->GetParameter(1);
          muBarError[pixel] = fBarFit[pixel]->GetParError(1);
          if (int(fitResult) != 0 || !fitResult.Get() || fitResult->CovMatrixStatus() < 2 ||
              !std::isfinite(muBar[pixel]) || !std::isfinite(muBarError[pixel]) || muBarError[pixel] <= 0.0 ||
              muBar[pixel] <= FitMin || muBar[pixel] >= FitMax || fBarFit[pixel]->GetParameter(0) <= 0.0 ||
              !std::isfinite(fBarFit[pixel]->GetParameter(2)) || std::abs(fBarFit[pixel]->GetParameter(2)) <= 0.0) {
            std::cout << "  Skipping donor pixel " << pixel << ": unusable fit.\n";
            continue;
          }
          const int pixelID = myBar * pixelsPerBar + pixel;
          const double proposedOffset = existingOffsets[pixelID] + muBar[pixel] - mu0;
          const double weight = 1.0 / (muBarError[pixel] * muBarError[pixel]);
          TF1* drawnFit = hBarDt[pixel]->GetFunction(fBarFit[pixel]->GetName());
          if (drawnFit) {
            drawnFit->ResetBit(TF1::kNotDraw);
            drawnFit->SetLineColor(kRed + 1);
            drawnFit->SetLineWidth(2);
          }
          hasIndividualFit[pixel] = true;
          barOffset[pixel] = proposedOffset;
          barOffsetError[pixel] = myBar == controlBar ? NAN : std::hypot(muBarError[pixel], mu0Error);
          weightedMuSum += weight * muBar[pixel];
          weightSum += weight;
          ++donors;
          std::cout << "  Donor local pixel " << pixel << " (ID " << pixelID << "): old="
                    << existingOffsets[pixelID] << ", mu=" << muBar[pixel]
                    << ", proposed=" << proposedOffset << " ns\n";
        }
      }

      if (!std::isfinite(weightSum) || weightSum <= 0.0 || !std::isfinite(weightedMuSum)) {
        std::cout << "[CDet pulse offsets] Bar " << myBar << " pixels " << firstPixel << "-" << lastPixel - 1
                  << ": unresolved (no usable donor weights).\n";
        continue;
      }
      const double groupMu = weightedMuSum / weightSum;
      const double groupMuError = std::sqrt(1.0 / weightSum);
      const double groupCorrection = groupMu - mu0;
      // The common reference uncertainty is added once, not reduced by donor count.
      // This assumes donor fits and control fits are independent.
      const double groupError = myBar == controlBar ? NAN : std::sqrt(1.0 / weightSum + mu0Error * mu0Error);
      int inheritedPixels = 0;
      for (Int_t pixel = firstPixel; pixel < lastPixel; pixel++) {
        if (hasIndividualFit[pixel]) continue; // Never replace a successful individual fit.
        muBar[pixel] = groupMu;
        muBarError[pixel] = groupMuError;
        barOffset[pixel] = existingOffsets[myBar * pixelsPerBar + pixel] + groupCorrection;
        barOffsetError[pixel] = groupError;
        ++inheritedPixels;
        std::cout << "  Inherited local pixel " << pixel << ": mu=" << groupMu
                  << ", c=" << groupCorrection << ", proposed=" << barOffset[pixel] << " ns\n";
      }
      std::cout << "[CDet pulse offsets] Bar " << myBar << " pixels " << firstPixel << "-" << lastPixel - 1
                << ": weighted donor mu = " << groupMu << " +/- " << groupMuError
                << " ns; inherited c = " << groupCorrection << " +/- " << groupError << " ns.\n"
                << "  Kept " << donors << " individual fit(s); assigned " << inheritedPixels << " other pixels.\n";
      if (myBar == controlBar) std::cout << "  Error not estimated: donors overlap the control reference.\n";
    }
    std::cout << "[CDet pulse offsets] Proposals only: inspect peaks before assigning to unresolved pixels;\n"
              << "  preserve established offsets. Numerical fit checks do not establish a signal peak.\n";

    for (int which = 0; which < 2; ++which) {
      for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
        canvases[which]->cd(pixel + 1);
        gPad->Modified();
        gPad->Update();
      }
      canvases[which]->Modified();
      canvases[which]->Update();
    }
    cBarDt->SaveAs(studyPDF);
    cControlDt->SaveAs(controlPDF);

    if (writeCSV) {
      std::ostringstream csvRows;
      csvRows << std::setprecision(17);
      int resolvedPixels = 0;
      for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
        csvRows << "Bar " << myBar << " pixel " << pixel << ',';
        if (std::isfinite(barOffset[pixel])) {
          const double correction = muBar[pixel] - mu0;
          csvRows << muBar[pixel] << ',' << muBarError[pixel] << ',' << correction << ',';
          // Overlapping control/donor fits need covariance treatment; error stays blank.
          if (std::isfinite(barOffsetError[pixel])) csvRows << barOffsetError[pixel];
          ++resolvedPixels;
        } else {
          csvRows << ",,,";
        }
        csvRows << '\n';
      }
      TString csvPath = csvFile ? csvFile : "";
      if (csvPath.IsNull()) {
        if (gSystem->AccessPathName("PusleStudy") && gSystem->mkdir("PusleStudy", true) != 0) {
          std::cerr << "[CDet pulse offsets] Cannot create PusleStudy output directory.\n";
          return;
        }
        csvPath = TString::Format("PusleStudy/run%d_bar%03d_%s_offsets.csv",
                                  config.GetValue("analysis.run_number", 0), myBar, fullNINO ? "full" : "halves");
      }
      if (!csvPath.EndsWith(".csv")) {
        std::cerr << "[CDet pulse offsets] Output filename must end in .csv.\n";
        return;
      }
      std::ofstream csv(csvPath.Data());
      if (!csv) {
        std::cerr << "[CDet pulse offsets] Cannot write CSV: " << csvPath << '\n';
        return;
      }
      csv << "ID,mu,mu error,c,c error\n" << csvRows.str();
      csv.close();
      if (!csv) {
        std::cerr << "[CDet pulse offsets] CSV write failed: " << csvPath << '\n';
        return;
      }
      std::cout << "[CDet pulse offsets] Saved " << pixelsPerBar << " pixel rows to " << csvPath
                << " (" << resolvedPixels << " resolved; remaining rows have blank values).\n"
                << "  Successful donors retain fitted mu; other resolved pixels inherit their group's weighted mu.\n"
                << "  c = mu - mu0; proposed final offset = existing pixel offset + c (ns).\n";
    }
  } //end viewOnly check
} //end extractBarOffsets

// 8. Background study. Uses the saved vectors; run the main macro first.
// ECal time is filled once per event passing the energy/time selections below.
// The selected bar is controlBar (fourth argument, default from the config).
// The rx and dt-versus-rx distributions split events AFTER this routine's corrected-LE/ToT
// cuts: N >= 2 versus N == 1 selected pulses in that bar. Repeated pulses on the
// same pixel count separately; pulses on other bars do not affect N.
void plotPulseStudyBackground(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                             const char* outputTag = nullptr, int studyBar = -1, int controlBar = -1, int controlPixel = -1,
                             double BinWidth = NAN, double DtMin = NAN, double DtMax = NAN,
                             double ECalBinWidth = 0.5, double ECalDisplayMin = NAN, double ECalDisplayMax = NAN,
                             double PositionBinWidth = NAN, double PositionMin = NAN, double PositionMax = NAN,
                             double RxBinWidth = 0.002, double RxMin = NAN, double RxMax = NAN)
{
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pulse background] Cannot read configuration with config.version: 1.\n";
    return;
  }
  if (studyBar == -1) studyBar = config.GetValue("pulse_pixel_timing.study_bar", 19);
  if (controlBar == -1) controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  if (controlPixel == -1) controlPixel = config.GetValue("pulse_pixel_timing.control_pixel", 9);
  if (std::isnan(BinWidth)) BinWidth = config.GetValue("pulse_study.dt_bin_width_ns", 1.0);
  if (std::isnan(DtMin)) DtMin = config.GetValue("display.ecal_cdet_dt_min", -90.0);
  if (std::isnan(DtMax)) DtMax = config.GetValue("display.ecal_cdet_dt_max", 0.0);

  if (std::isnan(ECalDisplayMin)) ECalDisplayMin = config.GetValue("display.ecal_time_min", -15.0);
  if (std::isnan(ECalDisplayMax)) ECalDisplayMax = config.GetValue("display.ecal_time_max", 15.0);

  // Existing secondary-cut settings.
  // Saved vectors already have the main macro's loading cuts applied.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 0) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 0) != 0;
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08); // Existing rx-cut key.

  if (std::isnan(PositionBinWidth)) PositionBinWidth = config.GetValue("pulse_study.position_bin_width_m", 0.005);
  if (std::isnan(PositionMin)) PositionMin = config.GetValue("pulse_study.position_min_m", -1.6);
  if (std::isnan(PositionMax)) PositionMax = config.GetValue("pulse_study.position_max_m", 1.6);
  if (std::isnan(RxMin)) RxMin = -xRes;
  if (std::isnan(RxMax)) RxMax = xRes;
  if (!std::isfinite(PositionBinWidth) || PositionBinWidth <= 0 ||
      !std::isfinite(PositionMin) || !std::isfinite(PositionMax) || PositionMin >= PositionMax ||
      !std::isfinite((PositionMax - PositionMin) / PositionBinWidth) ||
      std::ceil((PositionMax - PositionMin) / PositionBinWidth) > std::numeric_limits<int>::max() - 2 ||
      !std::isfinite(RxBinWidth) || RxBinWidth <= 0 ||
      !std::isfinite(RxMin) || !std::isfinite(RxMax) || RxMin >= RxMax ||
      !std::isfinite((RxMax - RxMin) / RxBinWidth) ||
      std::ceil((RxMax - RxMin) / RxBinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse background] Invalid position or rx display settings.\n";
    return;
  }

  if (studyBar < 0 || studyBar >= numBars || controlBar < 0 || controlBar >= numBars ||
      !std::isfinite(BinWidth) || BinWidth <= 0 || !std::isfinite(DtMin) ||
      !std::isfinite(DtMax) || DtMin >= DtMax ||
      !std::isfinite((DtMax - DtMin) / BinWidth) ||
      std::ceil((DtMax - DtMin) / BinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse background] Invalid bar IDs or dt display settings.\n";
    return;
  }
  if (!std::isfinite(ECalBinWidth) || ECalBinWidth <= 0 ||
      !std::isfinite(ECalDisplayMin) || !std::isfinite(ECalDisplayMax) || ECalDisplayMin >= ECalDisplayMax ||
      !std::isfinite((ECalDisplayMax - ECalDisplayMin) / ECalBinWidth) ||
      std::ceil((ECalDisplayMax - ECalDisplayMin) / ECalBinWidth) > std::numeric_limits<int>::max() - 2 ||
      !std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT ||
      !std::isfinite(ECalMinE) || !std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE) {
    std::cerr << "[CDet pulse background] Invalid ECal display settings or cuts.\n";
    return;
  }
  if (vPulseID.empty()) {
    std::cerr << "[CDet pulse background] No saved pulse events; run the main macro first.\n";
    return;
  }
  const int nECalTimeBins = int(std::ceil((ECalDisplayMax - ECalDisplayMin) / ECalBinWidth));
  const int nDtBins = int(std::ceil((DtMax - DtMin) / BinWidth));
  const int nPositionBins = int(std::ceil((PositionMax - PositionMin) / PositionBinWidth));
  const int nRxBins = int(std::ceil((RxMax - RxMin) / RxBinWidth));
  const int nCDetTotBins = int(std::ceil((TotMax - TotMin) / BinWidth));
  TString tag = outputTag ? outputTag : config.GetValue("pulse_study.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d", config.GetValue("analysis.run_number", 0));
  const TString allowedTagCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowedTagCharacters.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet pulse background] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/PusleStudy/%s", gSystem->WorkingDirectory(), tag.Data());

  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet pulse background] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }

  // Replace this routine's canvases/histograms on repeated calls.
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundECalTime");
  static TH1D* hECalT = nullptr;
  delete hECalT;

  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundCDetTotVsRx");
  static TH2D* hCDetTotVsRx = nullptr;
  delete hCDetTotVsRx;
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundDtVsRx");
  static TH2D* hDtVsRx = nullptr;
  delete hDtVsRx;
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundDt");
  static TH1D* hDt = nullptr;
  delete hDt;
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundCDetX");
  static TH1D* hCDetX = nullptr;
  delete hCDetX;
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundCDetXCorr");
  static TH1D* hCDetXCorr = nullptr;
  delete hCDetXCorr;
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundCDetXCorrVsX");
  static TH2D* hCDetXCorrVsX = nullptr;
  delete hCDetXCorrVsX;
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundRxMultiple");
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundRxSingle");
  delete hPulseBackgroundRxMultiple;
  delete hPulseBackgroundRxSingle;
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundDtVsRxMultiple");
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBackgroundDtVsRxSingle");
  delete hPulseBackgroundDtVsRxMultiple;
  delete hPulseBackgroundDtVsRxSingle;

  // Define histograms before the event/hit loop.
  hECalT = new TH1D("hPulseBackgroundECalTime", TString::Format("ECal ADC time: %.1f #leq t #leq %.1f ns;ECal ADC time (ns);Events", ECalMinT, ECalMaxT), nECalTimeBins, ECalDisplayMin, ECalDisplayMax);
  hECalT->SetDirectory(nullptr);
  hDtVsRx = new TH2D("hPulseBackgroundDtVsRx", TString::Format("Bar %d: dt versus rx;rx (m);Stored pulse.ecal_residual (ns);Pulses", controlBar), nRxBins, RxMin, RxMax, nDtBins, DtMin, DtMax);
  hDt = new TH1D("hPulseBackgroundDt", TString::Format("Bar %d: dt;Stored pulse.ecal_residual (ns);Pulses", controlBar), nDtBins, DtMin, DtMax);
  hCDetX = new TH1D("hPulseBackgroundCDetX", TString::Format("Bar %d: channel x;Channel CDet x (m);Pulses", controlBar), nPositionBins, PositionMin, PositionMax);
  hCDetXCorr = new TH1D("hPulseBackgroundCDetXCorr", TString::Format("Bar %d: aligned x;Aligned CDet x (m);Pulses", controlBar), nPositionBins, PositionMin, PositionMax);
  hCDetXCorrVsX = new TH2D("hPulseBackgroundCDetXCorrVsX", TString::Format("Bar %d: aligned x versus channel x;Channel CDet x (m);Aligned CDet x (m);Pulses", controlBar), nPositionBins, PositionMin, PositionMax, nPositionBins, PositionMin, PositionMax);
  hCDetTotVsRx = new TH2D("hPulseBackgroundCDetTotVsRx", TString::Format("Bar %d Pixel %d: CDet Tot versus rx;rx (m);CDet Tot (ns);Pulses", controlBar, controlPixel), nRxBins, RxMin, RxMax, nCDetTotBins, TotMin, TotMax);
  hCDetTotVsRx->SetDirectory(nullptr);
  hDtVsRx->SetDirectory(nullptr);
  hDt->SetDirectory(nullptr);
  hCDetX->SetDirectory(nullptr);
  hCDetXCorr->SetDirectory(nullptr);
  hCDetXCorrVsX->SetDirectory(nullptr);
  hPulseBackgroundRxMultiple = new TH1D("hPulseBackgroundRxMultiple", TString::Format("Bar %d: at least two selected pulses/event;rx (m);Pulses", controlBar), nRxBins, RxMin, RxMax);
  hPulseBackgroundRxSingle = new TH1D("hPulseBackgroundRxSingle", TString::Format("Bar %d: exactly one selected pulse/event;rx (m);Pulses", controlBar), nRxBins, RxMin, RxMax);
  hPulseBackgroundRxMultiple->SetDirectory(nullptr);
  hPulseBackgroundRxSingle->SetDirectory(nullptr);
  hPulseBackgroundDtVsRxMultiple = new TH2D("hPulseBackgroundDtVsRxMultiple", TString::Format("Bar %d: at least two selected pulses/event;rx (m);Stored pulse.ecal_residual (ns);Pulses", controlBar), nRxBins, RxMin, RxMax, nDtBins, DtMin, DtMax);
  hPulseBackgroundDtVsRxSingle = new TH2D("hPulseBackgroundDtVsRxSingle", TString::Format("Bar %d: exactly one selected pulse/event;rx (m);Stored pulse.ecal_residual (ns);Pulses", controlBar), nRxBins, RxMin, RxMax, nDtBins, DtMin, DtMax);
  hPulseBackgroundDtVsRxMultiple->SetDirectory(nullptr);
  hPulseBackgroundDtVsRxSingle->SetDirectory(nullptr);

  Long64_t multiplePulseEvents = 0, singlePulseEvents = 0;
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    const double e_ECal = vECalE[event];
    const double t_ECal = vECalAdcTime[event];
    const double x_ECal = vECalX[event];
    const double y_ECal = vECalY[event];

    bool good_ecal_time = t_ECal >= ECalMinT && t_ECal <= ECalMaxT;
    bool good_ecal_energy = e_ECal >= ECalMinE && e_ECal <= ECalMaxE;
    const bool good_ecal_event = good_ecal_time && good_ecal_energy;
    if (!good_ecal_event) continue;
    hECalT->Fill(t_ECal);

    // Collect this event's selected bar pulses before deciding its multiplicity.
    std::vector<double> eventBarRx;
    std::vector<double> eventBarDt;
    std::vector<double> eventPixel;
    std::vector<double> eventTot;
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars * pixelsPerBar || std::floor(id) != id) continue;
      const int bar = int(id) / pixelsPerBar;
      const int pixel = int(id) % pixelsPerBar;
      if (bar != controlBar) continue;
      const double dt = vPulseDt[event][hit]; // pulse.ecal_residual, ns.
      const double le = vPulseLe[event][hit]; // Uncorrected LE, already in ns.
      const double leCorr = vPulseLeCorr[event][hit];
      const double tot = vPulseToT[event][hit];
      const double rx = vPulseXResiduals[event][hit];
      const double ry = vPulseYResiduals[event][hit];
      const double x = vPulseX[event][hit];
      const double xCorr = vPulseXCorr[event][hit];
      const double y = vPulseY[event][hit];
      const double z = vPulseZ[event][hit];
      const double xProj = vECalXProjection[event][hit];
      const double yProj = vECalYProjection[event][hit];
      const bool goodCalib = vPulseGoodCalib[event][hit];

      bool good_cdet_le = leCorr >= LeMin && leCorr <= LeMax;
      bool good_cdet_tot = tot >= TotMin && tot <= TotMax;
      const bool good_cdet_hit = good_cdet_le && good_cdet_tot; //good cdet hit condition
      if (!good_cdet_hit) continue;
      hDtVsRx->Fill(rx,dt);
      hDt->Fill(dt);
      hCDetX->Fill(x);
      hCDetXCorr->Fill(xCorr);
      hCDetXCorrVsX->Fill(x,xCorr);
      eventBarRx.push_back(rx);
      eventBarDt.push_back(dt);
      eventPixel.push_back(pixel);
      eventTot.push_back(tot);
    }

    // Only this event determines its category; display bounds do not change N.
    TH1D* rxHistogram = nullptr;
    TH2D* dtVsRxHistogram = nullptr;
    if (eventBarRx.size() >= 2) {
      ++multiplePulseEvents;
      rxHistogram = hPulseBackgroundRxMultiple;
      dtVsRxHistogram = hPulseBackgroundDtVsRxMultiple;
    } else if (eventBarRx.size() == 1) {
      ++singlePulseEvents;
      rxHistogram = hPulseBackgroundRxSingle;
      dtVsRxHistogram = hPulseBackgroundDtVsRxSingle;
    }
    if (rxHistogram) {
      for (size_t pulse = 0; pulse < eventBarRx.size(); ++pulse) {
        if (eventPixel[pulse] != controlPixel) continue; // Only the selected pixel.
        const double rx = eventBarRx[pulse];
        const double dt = eventBarDt[pulse];
        if (std::isfinite(rx)) rxHistogram->Fill(rx);
        if (std::isfinite(rx) && std::isfinite(dt)) dtVsRxHistogram->Fill(rx, dt);
        if (std::isfinite(rx) && std::isfinite(eventTot[pulse])) hCDetTotVsRx->Fill(rx, eventTot[pulse]);
      }
    }
  }

  // Draw after filling; keep the caller's ROOT mode so canvases stay interactive.
  TCanvas* cCDetTotVsRx = new TCanvas("cPulseBackgroundCDetTotVsRx", TString::Format("Bar %d Pixel %d: CDet Tot versus rx", controlBar, controlPixel), 900, 650);
  gPad->SetRightMargin(0.15);
  hCDetTotVsRx->Draw("COLZ");
  cCDetTotVsRx->Modified();
  cCDetTotVsRx->Update();
  cCDetTotVsRx->SaveAs(TString::Format("%s/bar%03d_pixel%02d_background_cdet_tot_vs_rx.pdf", outputDirectory.Data(), controlBar, controlPixel));

  TCanvas* cECalT = new TCanvas("cPulseBackgroundECalTime", TString::Format("ECal ADC time: %.1f <= t <= %.1f ns", ECalMinT, ECalMaxT), 900, 650);
  hECalT->Draw("HIST");
  cECalT->Modified();
  cECalT->Update();
  cECalT->SaveAs(TString::Format("%s/background_ecal_adc_time.pdf", outputDirectory.Data()));

  TCanvas* cDtVsRx = new TCanvas("cPulseBackgroundDtVsRx", TString::Format("Bar %d: dt versus rx", controlBar), 900, 650);
  hDtVsRx->SetStats(false);
  gPad->SetRightMargin(0.15);
  hDtVsRx->Draw("COLZ");
  cDtVsRx->Modified();
  cDtVsRx->Update();
  cDtVsRx->SaveAs(TString::Format("%s/bar%03d_background_dt_vs_rx.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cDt = new TCanvas("cPulseBackgroundDt", TString::Format("Bar %d: dt", controlBar), 900, 650);
  hDt->Draw("HIST");
  cDt->Modified();
  cDt->Update();
  cDt->SaveAs(TString::Format("%s/bar%03d_background_dt.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cCDetX = new TCanvas("cPulseBackgroundCDetX", TString::Format("Bar %d: channel x", controlBar), 900, 650);
  hCDetX->Draw("HIST");
  cCDetX->Modified();
  cCDetX->Update();
  cCDetX->SaveAs(TString::Format("%s/bar%03d_background_x.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cCDetXCorr = new TCanvas("cPulseBackgroundCDetXCorr", TString::Format("Bar %d: aligned x", controlBar), 900, 650);
  hCDetXCorr->Draw("HIST");
  cCDetXCorr->Modified();
  cCDetXCorr->Update();
  cCDetXCorr->SaveAs(TString::Format("%s/bar%03d_background_x_corr.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cCDetXCorrVsX = new TCanvas("cPulseBackgroundCDetXCorrVsX", TString::Format("Bar %d: aligned x versus channel x", controlBar), 900, 650);
  hCDetXCorrVsX->SetStats(false);
  gPad->SetRightMargin(0.15);
  hCDetXCorrVsX->Draw("COLZ");
  cCDetXCorrVsX->Modified();
  cCDetXCorrVsX->Update();
  cCDetXCorrVsX->SaveAs(TString::Format("%s/bar%03d_background_x_corr_vs_x.pdf", outputDirectory.Data(), controlBar));

  // Identical bins, axes, count scale, and styling; separate windows side by side.
  TH1D* rxHistograms[] = {hPulseBackgroundRxMultiple, hPulseBackgroundRxSingle};
  const char* categories[] = {"Multiple", "Single"};
  const Long64_t eventCounts[] = {multiplePulseEvents, singlePulseEvents};
  const double rxMaximum = 1.2 * std::max(1.0, std::max(rxHistograms[0]->GetMaximum(), rxHistograms[1]->GetMaximum()));
  for (int category = 0; category < 2; ++category) {
    TCanvas* canvas = new TCanvas(TString::Format("cPulseBackgroundRx%s", categories[category]), TString::Format("Bar %d: %s-pulse events", controlBar, categories[category]), 40 + 720 * category, 80, 700, 600);
    TH1D* histogram = rxHistograms[category];
    histogram->SetStats(false);
    histogram->SetLineColor(kBlue + 1);
    histogram->SetLineWidth(2);
    histogram->SetMinimum(0.0);
    histogram->SetMaximum(rxMaximum);
    histogram->Draw("HIST");
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.032);
    label.DrawLatex(0.14, 0.84, TString::Format("Events: %lld; finite-rx pulses: %.0f", eventCounts[category], histogram->GetEntries()));
    label.DrawLatex(0.14, 0.79, TString::Format("Outside axes: %.0f pulses", histogram->GetEntries() - histogram->Integral()));
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(TString::Format("%s/bar%03d_background_rx_%s.pdf", outputDirectory.Data(), controlBar, categories[category]));
  }

  // Compare raw pulse counts with identical axes and a shared colour scale.
  TH2D* dtVsRxHistograms[] = {hPulseBackgroundDtVsRxMultiple, hPulseBackgroundDtVsRxSingle};
  const double dtVsRxMaximum = std::max(1.0, std::max(dtVsRxHistograms[0]->GetMaximum(), dtVsRxHistograms[1]->GetMaximum()));
  for (int category = 0; category < 2; ++category) {
    TCanvas* canvas = new TCanvas(TString::Format("cPulseBackgroundDtVsRx%s", categories[category]), TString::Format("Bar %d: dt versus rx, %s-pulse events", controlBar, categories[category]), 40 + 920 * category, 120, 900, 650);
    TH2D* histogram = dtVsRxHistograms[category];
    histogram->SetStats(false);
    histogram->SetMinimum(0.0);
    histogram->SetMaximum(dtVsRxMaximum);
    canvas->SetRightMargin(0.15);
    canvas->SetTopMargin(0.20);
    histogram->Draw("COLZ");
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.028);
    label.DrawLatex(0.14, 0.87, TString::Format("Events: %lld; finite dt/rx pulses: %.0f", eventCounts[category], histogram->GetEntries()));
    label.DrawLatex(0.14, 0.82, TString::Format("Outside axes: %.0f pulses", histogram->GetEntries() - histogram->Integral()));
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(TString::Format("%s/bar%03d_background_dt_vs_rx_%s.pdf", outputDirectory.Data(), controlBar, categories[category]));
  }

  std::cout << "[CDet pulse background] Filled ECal time for " << hECalT->GetEntries()
            << " selected events and CDet histograms for " << hDt->GetEntries()
            << " selected pulses in bar " << controlBar << ". Ten canvases remain open.\n"
            << "  rx split after corrected-LE/ToT cuts: " << multiplePulseEvents
            << " events with >=2 pulses, " << singlePulseEvents << " events with exactly 1 pulse.\n"
            << "  rx histograms count pulses, not events; repeated pulses on one pixel count separately.\n"
            << "  PDF output: " << outputDirectory << '\n';
}
