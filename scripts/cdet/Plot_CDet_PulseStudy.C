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
// Shares CDet_run<run>_pulse_pixel_timing.conf with the pulse-timing macro.
// Cut arguments override that file for this call only; NAN means use the file.
// Loading always applies TotMin/TotMax. The secondary plotting ToT switch is
// pulse_study.apply_tot_cut; it does not disable the loading cut.
// Add your pulse readers, vectors and hit cuts in the marked sections below.

#include <TChain.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TROOT.h>
#include <TPad.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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
// The dt-versus-rx diagnostic instead uses ECal ADC time minus raw LE.
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
    hPulseBarDtVsRx[selected] = new TH2D(TString::Format("hPulseSummary%sBarDtVsRx", roles[selected]), TString::Format("%s bar %d: manual dt versus rx, all pixels;rx (m);t_{ECal,ADC} - t_{CDet,raw} (ns);Pulses", roles[selected], bar), 100, -xRes, xRes, nDtBins, DtMin, DtMax);
    hPulseBarDtVsRx[selected]->SetDirectory(nullptr);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      const int id = bar * pixelsPerBar + pixel;
      hPulseSummaryPixelDt[selected][pixel] = new TH1D(TString::Format("hPulseSummary%sPixel%dDt", roles[selected], pixel), TString::Format("Pixel %d (ID %d): branch dt;Stored pulse.ecal_residual (ns);Pulses", pixel, id), nDtBins, DtMin, DtMax);
      hPulseSummaryPixelLe[selected][pixel] = new TH1D(TString::Format("hPulseSummary%sPixel%dLe", roles[selected], pixel), TString::Format("Pixel %d (ID %d): raw LE;Raw CDet LE (ns);Pulses", pixel, id), nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelDt[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLe[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLeVsToT[selected][pixel] = new TH2D(TString::Format("hPulseSummary%sPixel%dLeVsToT", roles[selected], pixel), TString::Format("Pixel %d (ID %d): raw LE versus ToT;ToT (ns);Raw CDet LE (ns);Pulses", pixel, id), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelLeVsToT[selected][pixel]->SetDirectory(nullptr);
      hPulsePixelDtVsRx[selected][pixel] = new TH2D(TString::Format("hPulseSummary%sPixel%dDtVsRx", roles[selected], pixel), TString::Format("Pixel %d (ID %d): manual dt versus rx;rx (m);t_{ECal,ADC} - t_{CDet,raw} (ns);Pulses", pixel, id), 100, -xRes, xRes, nDtBins, DtMin, DtMax);
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
    double t_ECal = vECalAdcTime[event];
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
      const double dt_manual = t_ECal - le;
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
        hPulseBarDtVsRx[selected]->Fill(rx, dt_manual);
        hPulsePixelDtVsRx[selected][pixel]->Fill(rx, dt_manual);
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

    TCanvas* cBarDtVsRx = new TCanvas(TString::Format("cPulseSummary%sBarDtVsRx", roles[selected]), TString::Format("%s bar %d: summed manual dt versus rx", roles[selected], bars[selected]), 900, 650);
    hPulseBarDtVsRx[selected]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hPulseBarDtVsRx[selected]->Draw("COLZ");
    cBarDtVsRx->Modified();
    cBarDtVsRx->Update();
    cBarDtVsRx->SaveAs(TString::Format("%s/bar%03d_%s_dt_manual_vs_rx.pdf", outputDirectory.Data(), bars[selected], roles[selected]));

    TCanvas* cPixelDtVsRx = new TCanvas(TString::Format("cPulseSummary%sPixelDtVsRx", roles[selected]), TString::Format("%s bar %d: pixel manual dt versus rx", roles[selected], bars[selected]), 1400, 1000);
    cPixelDtVsRx->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelDtVsRx->cd(pixel + 1);
      hPulsePixelDtVsRx[selected][pixel]->SetStats(false);
      gPad->SetRightMargin(0.15);
      hPulsePixelDtVsRx[selected][pixel]->Draw("COLZ");
    }
    cPixelDtVsRx->Modified();
    cPixelDtVsRx->Update();
    cPixelDtVsRx->SaveAs(TString::Format("%s/bar%03d_%s_pixels_dt_manual_vs_rx.pdf", outputDirectory.Data(), bars[selected], roles[selected]));
  }
  std::cout << "[CDet pulse bar summary] Drew study bar " << studyBar << " and control bar " << controlBar
            << "; selected " << filledPulses << " pulses. Sixteen canvases remain open.\n"
            << "  dt uses vPulseDt; LE versus ToT uses uncorrected vPulseLe. All times are in ns.\n"
            << "  dt versus rx uses ECal ADC time - vPulseLe against stored x residual.\n"
            << "  PDF output: " << outputDirectory << '\n';
}
