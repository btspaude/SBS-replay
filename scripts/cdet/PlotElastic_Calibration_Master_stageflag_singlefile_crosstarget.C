#pragma once

#include "CDetPlotStyle.h"

#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <THashList.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>      // for sscanf
#include <cstdlib>     // for strtod
#include <algorithm>   // for std::sort
#include <limits>
#include <TMath.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "CDetRunDataset.h"
#include <TCanvas.h>
#include <TControlBar.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TText.h>
#include <TProfile.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TParameter.h>
#include <TCutG.h>
#include <TEnv.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>

// Defaults used by the LH2 bar-timing diagnostics.  Loading a TEnv
// configuration updates these values; explicit function arguments still win.
double gCDetDiagnosticAcceptedTotMin = 4.0;
double gCDetDiagnosticAcceptedTotMax = 30.0;

static const std::unordered_set<std::string>& CDetConfigurationKeys()
{
  static const std::unordered_set<std::string> keys = {
    "config.version",
    "analysis.run_number", "analysis.events", "analysis.calibration_stage",
    "analysis.group_index", "analysis.min_segment", "analysis.max_segment",
    "analysis.le_min", "analysis.le_max", "analysis.tot_min", "analysis.tot_max",
    "analysis.ecal_time_min", "analysis.ecal_time_max",
    "analysis.ecal_energy_min", "analysis.ecal_energy_max",
    "analysis.layer1_hits_min", "analysis.layer1_hits_max",
    "analysis.layer2_hits_min", "analysis.layer2_hits_max",
    "analysis.x_difference_max", "analysis.x_offset", "analysis.y_offset",
    "analysis.layer_choice", "analysis.suppress_bad", "analysis.number_of_runs",
    "analysis.max_stream", "analysis.first_event", "analysis.use_reference_timing",
    "diagnostics.accepted_tot_min", "diagnostics.accepted_tot_max",
    "diagnostics.bar_ecal_fit_min", "diagnostics.bar_ecal_fit_max",
    "diagnostics.bar_ecal_peak_seed_min",
    "diagnostics.bar_ecal_peak_seed_max",
    "good_pulse.bin_width_ns", "good_pulse.le_min_ns",
    "good_pulse.le_max_ns", "good_pulse.tot_min_ns",
    "good_pulse.tot_max_ns", "good_pulse.recovered_only",
    "good_pulse.pair_residual_center_m",
    "good_pulse.pair_timing_center_ns",
    "good_pulse.pair_residual_scale_m",
    "good_pulse.pair_timing_scale_ns", "good_pulse.pair_cut_radius",
    "good_pulse.single_residual_center_m",
    "good_pulse.single_timing_center_ns",
    "good_pulse.single_residual_scale_m",
    "good_pulse.single_timing_scale_ns", "good_pulse.single_cut_radius",
    "good_pulse.opposite_side_enable", "good_pulse.opposite_dy_center_m",
    "good_pulse.opposite_dy_tolerance_m",
    "good_pulse.opposite_projected_y_center_m",
    "good_pulse.opposite_projected_y_max_m",
    "pair_scan.min_half_width_ns", "pair_scan.max_half_width_ns",
    "pair_scan.step_ns",
    "display.overwrite", "display.pixel", "display.histogram_width",
    "display.layer_dt_min", "display.layer_dt_max",
    "display.layer_dx_min", "display.layer_dx_max",
    "display.le_min", "display.le_max", "display.tot_min", "display.tot_max",
    "display.dt_hist_min", "display.dt_hist_max",
    "display.cdet_time_min", "display.cdet_time_max",
    "display.cdet_tot_min", "display.cdet_tot_max",
    "display.ecal_time_min", "display.ecal_time_max",
    "display.hcal_time_min", "display.hcal_time_max",
    "display.ecal_cdet_dt_min", "display.ecal_cdet_dt_max",
    "display.y_correction_refractive_index",
    "display.hit_mode", "display.best_hit_peak_mean",
    "display.best_hit_peak_sigma", "display.best_hit_nsigma",
    "display.allow_multiple_pairs", "display.x_bin_width", "display.x_min", "display.x_max",
    "display.z_bin_width", "display.z_min", "display.z_max"
  };
  return keys;
}

static bool LoadCDetConfiguration(TEnv& env, const char *configFile,
                                  const char *caller)
{
  if (!configFile || !configFile[0]) {
    std::cerr << "[" << caller << "] ERROR: configuration filename is empty.\n";
    return false;
  }
  if (env.ReadFile(configFile, kEnvLocal) != 0) {
    std::cerr << "[" << caller << "] ERROR: cannot read configuration file "
              << configFile << ".\n";
    return false;
  }
  if (env.GetValue("config.version", 0) != 1) {
    std::cerr << "[" << caller << "] ERROR: " << configFile
              << " must set config.version = 1.\n";
    return false;
  }

  TIter next(env.GetTable());
  while (TObject *record = next()) {
    const std::string key = record->GetName();
    if (CDetConfigurationKeys().count(key) == 0) {
      std::cerr << "[" << caller << "] ERROR: unknown configuration key '"
                << key << "' in " << configFile << ".\n";
      return false;
    }
  }

  const double acceptedTotMin = env.GetValue("diagnostics.accepted_tot_min", 4.0);
  const double acceptedTotMax = env.GetValue("diagnostics.accepted_tot_max", 30.0);
  if (acceptedTotMin >= acceptedTotMax) {
    std::cerr << "[" << caller << "] ERROR: diagnostics.accepted_tot_min must be less than "
              << "diagnostics.accepted_tot_max in " << configFile << ".\n";
    return false;
  }
  const double barECalFitMin =
      env.GetValue("diagnostics.bar_ecal_fit_min", -45.0);
  const double barECalFitMax =
      env.GetValue("diagnostics.bar_ecal_fit_max", -10.0);
  const double barECalPeakSeedMin =
      env.GetValue("diagnostics.bar_ecal_peak_seed_min", -40.0);
  const double barECalPeakSeedMax =
      env.GetValue("diagnostics.bar_ecal_peak_seed_max", -15.0);
  if (barECalFitMin >= barECalFitMax ||
      barECalPeakSeedMin >= barECalPeakSeedMax ||
      barECalPeakSeedMin < barECalFitMin ||
      barECalPeakSeedMax > barECalFitMax) {
    std::cerr << "[" << caller << "] ERROR: the bar ECal-CDet peak-seed "
                 "interval must be ordered and contained within the fit "
                 "interval in " << configFile << ".\n";
    return false;
  }
  gCDetDiagnosticAcceptedTotMin = acceptedTotMin;
  gCDetDiagnosticAcceptedTotMax = acceptedTotMax;
  std::cout << "[" << caller << "] Loaded configuration " << configFile << ".\n";
  std::cout << "[" << caller << "] Diagnostic accepted-ToT window: ("
            << gCDetDiagnosticAcceptedTotMin << ", "
            << gCDetDiagnosticAcceptedTotMax
            << ") ns for polygon-aware pixel-calibration diagnostics. "
               "Run analysis and plotCDetLayersTimeComp do not apply saved polygons.\n";
  return true;
}

std::vector<TCanvas*> canvas_vector;

static const int TDCmult_cut = 100;
static const double xcut = 998.0;
//static const int nhitcutlow = 2;
//static const int nhitcuthigh = 20;
static const double TDC_calib_to_ns = 0.01;
static const double HotChannelRatio = .01;
static const double RawSinglesWindowSeconds = 60.0e-9;
static const double ECalClusterWindowSeconds = 250.0e-9;
static const double ECalClusterEnergyBinWidthGeV = 0.01;

static const int NumPaddles = 16;
static const int NumBars = 14;
static const int NumHalfBarsPerBank = NumBars/2;
static const int NumLayers = 2;
static const int NumSides = 2;
static const int NumModules = 3;
static const int NumHalfModules = NumModules*NumSides*NumLayers;

static const int NumCDetPaddles = NumHalfModules*NumBars*NumPaddles; //2688
static const int nRef = 4;
static const int NumRefPaddles = 4;
static const int nTdc = NumCDetPaddles+NumRefPaddles; //2704

static const int NumSidesTotal = NumSides*NumLayers;
static const int NumCDetPaddlesPerSide = NumCDetPaddles/NumSidesTotal; //672
static const int NumLogicalPaddlesPerSide = NumCDetPaddlesPerSide+nRef; //676

static const int nBarsADC = 0;
static const double ADCCUT = 150.;   //100.0

static const double ECal_dist = 6.144; // from db_run.dat as of 2026-08-31
static const double CDet_y_half_length = 0.30;

static const double XCorr1 = 1.08; // Layer 1 CDet x scale relative to the ECal projection
static const double XCorr2 = 1.08; // Layer 2 CDet x scale relative to the ECal projection
static const double CDetXOffset1 = 0.03; // subtract from scaled Layer 1 x (m)
static const double CDetXOffset2 = 0.03; // subtract from scaled Layer 2 x (m)

static double CDetXCorrForLayer(int layer)
{
  return layer == 0 ? XCorr1 : XCorr2;
}

static double CDetXOffsetForLayer(int layer)
{
  return layer == 0 ? CDetXOffset1 : CDetXOffset2;
}

static double CorrectCDetX(double rawX, int layer)
{
  return rawX * CDetXCorrForLayer(layer) - CDetXOffsetForLayer(layer);
}

int NXDiffBins;
double XDiffLow;
double XDiffHigh;

// For generating sorting cdet layer 1 and 2 hits into pairs
struct Cand {
  int i1;       // index into layer-1 per-event vectors
  int j2;       // index into layer-2 per-event vectors
  double dt;    // t2 - t1
  double dx;    // x2 - x1
  double score; // ranking metric (lower = better)
};

// For storing pairs for cdet layer 1 and layer 2 hits, eventually need to add y,z, tot
struct PairHit {
  double t1, tot1, x1, y1, z1;
  int    id1;

  double t2, tot2, x2, y2, z2;
  int    id2;

  double dt, dx;
  double score;
};

struct CDetYPairedSample {
  int halfBar1, halfBar2;
  double t1, t2;
};

std::vector<std::vector<PairHit>> pairs_CDet;
std::vector<TH1D*> gCDetPairedLeSpectra;
std::vector<TH2D*> gCDetBarECalTimingSpectra;
std::vector<std::vector<std::pair<double,double>>> gCDetHalfBarECalTimingSamples;
struct CDetFrozenECalTimingSample {
  size_t eventIndex;
  size_t pairIndex;
  int layer;
  int pixel;
  int halfBar;
  double ecalTime;
  double cdetTime;
  double tot;
};
std::vector<CDetFrozenECalTimingSample> gCDetFrozenECalTimingSamples;
std::vector<std::vector<std::pair<double,double>>> gCDetHalfBarECalYSamples;
std::vector<CDetYPairedSample> gCDetYPairedSamples;
TH2D *gCDetPairedMeanTimeVsECal = nullptr;
std::vector<double> gCDetAcceptedPairMeanTimes;
std::vector<double> gCDetProjectedHalfBarPairMeanTimes;

struct CDetDisplayHit {
  int id, layer, side, module, bar, pixel;
  double x, y, z, le, te, tot;
};

struct CDetDisplayEvent {
  size_t savedEventIndex;
  Long64_t treeEntry;
  int runNumber;
  int selectedBar;
  double ecalX, ecalY, ecalZ, ecalEnergy, ecalTime;
  std::vector<CDetDisplayHit> hits;
  std::vector<PairHit> selectedPairs;
};

struct CDetPixelReviewCandidate {
  int pixel;
  int entries;
  double score;
  double broadGroupDifference;
  double usedGroupDifference;
  double usedSigma;
  std::string source;
  std::string reasons;
};

std::vector<CDetDisplayEvent> gCDetDisplayEvents;
Long64_t gCDetDisplayIndex = -1;
TCanvas *gCDetEventCanvas = nullptr;
TCanvas *gCDetProjectedHalfBarTimingCanvas = nullptr;
TControlBar *gCDetEventControl = nullptr;
double gCDetDisplayXMin = -1.5;
double gCDetDisplayXMax = 1.5;
double gCDetDisplayZMin = 0.0;
double gCDetDisplayZMax = 7.0;
bool gCDetDisplayBestHits = false;
bool gCDetDisplayBestHitTimingValid = false;
double gCDetDisplayBestHitPeakMean = std::numeric_limits<double>::quiet_NaN();
double gCDetDisplayBestHitPeakSigma = std::numeric_limits<double>::quiet_NaN();
double gCDetDisplayBestHitNSigma = 3.0;
std::vector<CDetPixelReviewCandidate> gCDetPixelReviewQueue;
Long64_t gCDetPixelReviewIndex = -1;
TString gCDetPixelReviewCutFile = "CDet_pixel_quality_cuts.root";

void BuildCDetEventDisplay(int selectedBar, double diffMinCut, double diffMaxCut,
                           double xdiffMinCut, double xdiffMaxCut,
                           double tdiffECalCDetMin, double tdiffECalCDetMax,
                           double xMin, double xMax, double zMin, double zMax,
                           bool initialBestHits = false,
                           double bestHitPeakMean = std::numeric_limits<double>::quiet_NaN(),
                           double bestHitPeakSigma = std::numeric_limits<double>::quiet_NaN(),
                           double bestHitNSigma = 3.0);
void ShowCDetEvent(Long64_t displayIndex = 0);
void NextCDetEvent();
void PreviousCDetEvent();
void ShowAllCDetHits();
void ShowBestCDetHits();
void PrintCDetEvent();
void SaveCDetEvent();
void writeAllCDetPairedLeSpectra(TString outputFile = "CDet_all_bars_paired_le.root");
void showCDetPairedLeBar(int bar = 29,
                         TString inputFile = "CDet_all_bars_paired_le.root");
void writeAllCDetBarECalTimingDiagnostics(
    TString outputFile = "CDet_bar_ecal_timing.root",
    TString summaryFile = "CDet_bar_ecal_timing_summary.dat");
void showCDetBarECalTiming(
    int bar = 29, TString inputFile = "CDet_bar_ecal_timing.root");
void plotAllCDetBarECalTiming(
    TString outputFile = "CDet_bar_ecal_timing_overview.root",
    TString summaryFile = "CDet_bar_ecal_timing_overview_summary.dat");
void calibrateCDetHalfBarIntercepts(
    bool overwrite = false, double referenceECalTime = 22.0,
    int minEntries = 100, double maxInterceptError = 1.0,
    TString outputRoot = "CDet_halfbar_intercept_diagnostics.root",
    TString outputSummary = "CDet_halfbar_intercept_corrections.dat");
void reportCDetPairedTimeResolution(bool draw = true,
                                    double gaussianFitMin = 25.0,
                                    double gaussianFitMax = 35.0);
void extractCDetYPositionCalibration(
    TString outputFile = "CDet_y_position_calibration.dat",
    int minEntries = 100);
void plotCDetYPositionResolution(
    TString calibrationFile = "CDet_y_position_calibration.dat");
void ResetCalibrationGlobals();

// Interactive, persistent LE-versus-TOT selections for pixels whose physical
// population cannot be isolated reliably by the automatic one-dimensional fit.
void plotCDetPixelLeVsTot(int logicalPixelID = 468,
                          TString cutFile = "CDet_pixel_quality_cuts.root",
                          double totMin = 0.0, double totMax = 60.0,
                          double leMin = 0.0, double leMax = 60.0,
                          double binWidth = 0.5);
void editCDetPixelLeTotCut(int logicalPixelID = 468,
                           TString cutFile = "CDet_pixel_quality_cuts.root",
                           double totMin = 0.0, double totMax = 60.0,
                           double leMin = 0.0, double leMax = 60.0,
                           double binWidth = 0.5);
void buildCDetPixelLeTotReviewQueue(
    TString resultsFile = "CDet_pixel_timing_fit_results_manual_le_tot.dat",
    TString cutFile = "CDet_pixel_quality_cuts.root",
    TString queueOutput = "CDet_pixel_le_tot_review_queue.dat",
    int minEntries = 35, double centroidDifference = 2.0,
    double sigmaThreshold = 3.0);
void printCDetPixelLeTotReviewQueue(int maximumRows = 25);
void reviewCurrentCDetPixelLeTotCandidate();
void reviewNextCDetPixelLeTotCandidate();
void reviewPreviousCDetPixelLeTotCandidate();

// List of x-positions (or bins) for unused pixels ----- 1/19 verified correct
static std::vector<double> missingPixelBins = {
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
2563, 2570, 2576, 2591, 2592, 2607, 2611, 2621, 2633, 2636, 2643, 2650, 2656, 2657, 2675, 2679};

static const std::unordered_set<int> kUnusedCDetPixels = {
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
2563, 2570, 2576, 2591, 2592, 2607, 2611, 2621, 2633, 2636, 2643, 2650, 2656, 2657, 2675, 2679};

inline bool IsUnusedPixel(int elID) {
    return kUnusedCDetPixels.count(elID) != 0;
}


//const TString REPLAYED_DIR = TString(gSystem->Getenv("OUT_DIR")) + "/wrongdbRootfiles";
const TString REPLAYED_DIR = TString(gSystem->Getenv("OUT_DIR"));

TString GetAnalysedDir() {
  const char *analysed = gSystem->Getenv("ANALYSED_DIR");
  if (analysed && *analysed) return TString(analysed);

  const char *out = gSystem->Getenv("OUT_DIR");
  if (out && *out) return TString(out) + "/cdetFiles";

  return TString();
}

const TString ANALYSED_DIR = GetAnalysedDir();

// Parse the "segX_Y" part: returns true and fills firstSeg/lastSeg if found.
// -------- 4/21/2026 B.Spaude modified so we can use this for either a single run or put multiple runs in chain
bool GetSegRange(const TString& fname, int& firstSeg, int& lastSeg) {
  Ssiz_t pos = fname.Index("_seg");
  if (pos == kNPOS) return false;

  TString tail = fname(pos + 4, fname.Length() - (pos + 4));

  int a = -1, b = -1;
  if (sscanf(tail.Data(), "%d_%d", &a, &b) == 2) {
    firstSeg = a;
    lastSeg  = b;
    return true;
  }
  return false;
}

void AddRunFilesToChain(TChain *chain, const char *dir, int runnum,
                        int segMin = -1, int segMax = -1) {
  CDetRunDataset::AddToChain(chain, runnum, dir, segMin, segMax);
}
// 5/18/2026 B. Spaude: modified to allow for chaining different run groups
std::vector<int> ReadRunList(const char *runListFile, int groupIndex = 0) {
  std::vector<int> runs;
  std::unordered_set<int> seenRuns;
  std::ifstream infile(runListFile);

  if (!infile.is_open()) {
    std::cerr << "ERROR: Could not open run list file: "
              << runListFile << "\n";
    return runs;
  }

  if (groupIndex < 0) {
    std::cerr << "ERROR: groupIndex must be >= 0. Got "
              << groupIndex << "\n";
    return runs;
  }

  std::string targetHeader = TString::Format("[Group %d]", groupIndex).Data();

  bool inTargetGroup = false;
  std::string line;

  while (std::getline(infile, line)) {
    // Trim leading whitespace
    size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) continue;

    line = line.substr(first);

    // Skip comments
    if (line[0] == '#') continue;

    // Header line
    if (line[0] == '[') {
      if (line.find(targetHeader) == 0) {
        inTargetGroup = true;
      } else if (inTargetGroup) {
        // We reached the next group, so stop reading
        break;
      } else {
        inTargetGroup = false;
      }

      continue;
    }

    // Run lines inside requested group
    if (inTargetGroup) {
      std::stringstream ss(line);
      int runnum;

      if (ss >> runnum) {
        if (seenRuns.insert(runnum).second) {
          runs.push_back(runnum);
        } else {
          std::cerr << "WARNING: Duplicate run " << runnum
                    << " ignored in " << runListFile << "\n";
        }
      }
    }
  }

  if (runs.empty()) {
    std::cerr << "WARNING: No runs found for "
              << targetHeader
              << " in file "
              << runListFile << "\n";
  }

  return runs;
}

bool IsCrossTargetRun(int runNumber, const char *runListFile = "runs.txt") {
  std::ifstream infile(runListFile);
  if (!infile.is_open()) {
    infile.clear();
    infile.open(TString::Format("scripts/cdet/%s", runListFile).Data());
  }

  int listedRun = 0;
  while (infile >> listedRun) {
    if (listedRun == runNumber) return true;
  }
  return false;
}

void AddFitResultsToStatsBox(TH1 *hist, double fitMean, double fitMeanError, double fitSigma, double fitSigmaError, double fitChi2Ndf) {
  if (!hist || !gPad) return;
  gPad->Update();
  TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");
  if (!stats) return;
  stats->SetName(TString::Format("stats_%s", hist->GetName()));
  TLatex *meanLine = new TLatex(0.0, 0.0, TString::Format("#mu_{fit} = %.3f #pm %.3f ns", fitMean, fitMeanError));
  TLatex *sigmaLine = new TLatex(0.0, 0.0, TString::Format("#sigma_{fit} = %.3f #pm %.3f ns", fitSigma, fitSigmaError));
  TLatex *chi2Line = new TLatex(0.0, 0.0, TString::Format("#chi^{2}/NDF = %.2f", fitChi2Ndf));
  TText *defaultLine = (TText*)stats->GetListOfLines()->At(1);
  if (defaultLine) {
    meanLine->SetTextFont(defaultLine->GetTextFont());
    meanLine->SetTextSize(defaultLine->GetTextSize());
    meanLine->SetTextColor(defaultLine->GetTextColor());
    meanLine->SetTextAlign(defaultLine->GetTextAlign());
    sigmaLine->SetTextFont(defaultLine->GetTextFont());
    sigmaLine->SetTextSize(defaultLine->GetTextSize());
    sigmaLine->SetTextColor(defaultLine->GetTextColor());
    sigmaLine->SetTextAlign(defaultLine->GetTextAlign());
    chi2Line->SetTextFont(defaultLine->GetTextFont());
    chi2Line->SetTextSize(defaultLine->GetTextSize());
    chi2Line->SetTextColor(defaultLine->GetTextColor());
    chi2Line->SetTextAlign(defaultLine->GetTextAlign());
  }
  stats->GetListOfLines()->Add(meanLine);
  stats->GetListOfLines()->Add(sigmaLine);
  stats->GetListOfLines()->Add(chi2Line);
  hist->SetStats(kFALSE);
  stats->Draw("SAME");
  gPad->Modified();
  gPad->Update();
}

static double CDetTimingBackgroundRejectLow = 0.0;
static double CDetTimingBackgroundRejectHigh = 0.0;

double CDetTimingBackgroundGaussianReject(double *x, double *par) {
  if (CDetTimingBackgroundRejectLow <= x[0] && x[0] <= CDetTimingBackgroundRejectHigh) {
    TF1::RejectPoint();
    return 0.0;
  }
  return par[0]*TMath::Gaus(x[0], par[1], par[2], false);
}

// 4/21/2026 B. Spaude: Use this to get the run list from crossRuns.txt
void AddRunListFilesToChain(TChain *chain, const char *dir,
                            const char *runListFile,
                            int segMin = -1, int segMax = -1, int groupIndex = 0) {
  std::vector<int> runs = ReadRunList(runListFile, groupIndex);

  for (int runnum : runs) {
    AddRunFilesToChain(chain, dir, runnum, segMin, segMax);

  }
}

/* Create globals for vectors */
// Scalars (1D vectors)
std::vector<double> vheep_dpp;
std::vector<double> vheep_dt_ADC;
std::vector<double> vheep_ECalo;
std::vector<double> vheep_eprime_eth;
std::vector<double> vheep_dxECAL;
std::vector<double> vearm_ECal_x;
std::vector<double> vsbs_gemFPP_track_ntrack;
std::vector<double> vheep_dyECAL;

// Arrays (2D vectors)
std::vector<std::vector<double>> vsbs_tr_vz;
std::vector<std::vector<double>> vsbs_gemFPP_track_sclose;
std::vector<std::vector<double>> vsbs_gemFT_track_nhits;
std::vector<std::vector<double>> vsbs_gemFT_track_ngoodhits;

/* CDet & ECal Vectors */
//1D vectors
std::vector<double> vRefRawLe;
std::vector<double> vGoodRefRawLe; // ref LE time aligned to GOOD-event vectors (same indexing as vGoodLe)
std::vector<double> vRefRawTe;
std::vector<double> vRefRawTot;
std::vector<int>    vRefRawPMT;

std::vector<double> vRefGoodLe;
std::vector<double> vRefGoodTe;
std::vector<double> vRefGoodTot;
std::vector<int>    vRefGoodPMT;

std::vector<double> vAllRawLe;
std::vector<double> vAllRawTe;
std::vector<double> vAllRawTot;
std::vector<int> vAllRawPMT;
std::vector<int> vAllRawBar;

std::vector<double> vAllGoodLe;
std::vector<double> vAllGoodECalT;                 // per-hit ECal ADC time aligned with vAllGoodLe/Te
std::vector<std::vector<double>> vBarGoodLeECalT;        // per-bar per-hit ECal time aligned with vBarGoodLe

// --- ECal-time linear correction: remove correlation t_CDet vs t_ECal, then apply global shift
//double gECalFitP0 = -40.303;    // p0 from fit: <t_CDet> = p0 + p1*t_ECal
//double gECalFitP1 =  0.83357;   // p1 from fit
// After Stage 3
//double gECalFitP0 = -59.6812;    // p0 from fit: <t_CDet> = p0 + p1*t_ECal
//double gECalFitP1 =  0.816631;   // p1 from fit
// After Stage 4a
double gECalFitP0 = -43.1529;    // p0 from fit: <t_CDet> = p0 + p1*t_ECal
double gECalFitP1 =  0.9866;   // p1 from fit
double gTargetMeanLE = 30.0;    // desired mean corrected LE (ns)
double gECalDeltaShift = 0.0;   // fixed detector-wide shift applied after removing correlation
bool   gECalDeltaLoaded = false;
bool   gUseECalTimeCorr = true; // enable/disable ECal-time correction

// --- TOT time-walk correction (layer-dependent)
// Correction applied after pixel-offset and ECal-time corrections:
//   t_corr = t - p1*(1/sqrt(TOT) - 1/sqrt(TOT_ref))
// so that the correction is zero at TOT_ref.  Update the p1 values below
// with those obtained from plotGoodLeVsTotByLayer(...).
double gTimeWalkP1_L1 = 19.736;      // ns*sqrt(ns)
double gTimeWalkP1_L2 = 18.037;      // ns*sqrt(ns)
double gTimeWalkTotRef_L1 = 12.0; // ns
double gTimeWalkTotRef_L2 = 12.0; // ns
double gTimeWalkFitP0_L1 = 0.0;   // stored from plotting fit
double gTimeWalkFitP0_L2 = 0.0;   // stored from plotting fit
double gTimeWalkFitP1_L1 = 0.0;   // stored from plotting fit
double gTimeWalkFitP1_L2 = 0.0;   // stored from plotting fit
bool   gUseTimeWalkCorr = true;   // enable/disable time-walk correction
double gTimeWalkTotMin = 5.0;  // lower edge of the calibration fit domain
double gTimeWalkTotMax = 25.0; // upper edge of the calibration fit domain

inline int GetLayerFromID(int elID) {
  if (elID >= 0 && elID <= 1343) return 1;
  if (elID >= 1344 && elID <= 2687) return 2;
  return 0;
}

inline double GetTimeWalkCorrection(int layer, double tot) {
  if (!gUseTimeWalkCorr) return 0.0;
  if (!std::isfinite(tot) || tot <= 0.0) return 0.0;

  if (layer == 1 && std::isfinite(gTimeWalkTotRef_L1) && gTimeWalkTotRef_L1 > 0.0) {
    return gTimeWalkP1_L1 * (1.0/std::sqrt(tot) - 1.0/std::sqrt(gTimeWalkTotRef_L1));
  }
  if (layer == 2 && std::isfinite(gTimeWalkTotRef_L2) && gTimeWalkTotRef_L2 > 0.0) {
    return gTimeWalkP1_L2 * (1.0/std::sqrt(tot) - 1.0/std::sqrt(gTimeWalkTotRef_L2));
  }
  return 0.0;
}

std::vector<std::vector<double>> vBarGoodLe;
std::vector<std::vector<double>> vPaddleGoodLe;
std::vector<std::vector<double>> vPaddleGoodTot;
std::vector<std::vector<double>> vPaddleMatchHCalTime;
std::vector<TH1F*> hBarGoodLe; // one histogram per bar (PMT group), built in plotAllTDC()
std::vector<TH1F*> hPaddleGoodLe; // one histogram per paddle, built in plotAllTDC()
std::vector<TH1F*> hPaddleGoodTot;
std::vector<TH2F*> hPaddleLEvsTOT;
static const int NumPMTs = NumHalfModules*NumBars; // 168 (does not include 4 ref paddles)

std::vector<int> gPixelToffsetNhits;      // store Nhits for each physical logical pixel when computing toffsets
std::vector<double> gPixelToffsetCorr;    // size NumCDetPaddles, correction to ADD to LE/TE: (mean_all - mean_pixel)
// Immutable baseline corresponding to the currently populated analysis vectors.
// Offset extraction must add residuals to this snapshot, not to constants that
// may already have been updated by an earlier extraction in the same session.
std::vector<int> gAnalysisPixelToffsetNhits;
std::vector<double> gAnalysisPixelToffsetCorr;
bool gAnalysisPixelToffsetSnapshotValid = false;
bool gPixelToffsetLoaded = false;         // true if offsets were read from file
bool gECalParamsLoaded = false;
bool gTimeWalkParamsLoaded = false;
bool gCalibrationLoaded = false;
std::string gCalibrationFile = "CDet_calibration_dt.dat";
Int_t gNumEventsInRun = 0;
Int_t gRunNumber = 0;
Int_t gCalibrationStage = 0;
bool gLastCalibrationStageSucceeded = false;
bool gLastCalibrationFitSucceeded = false;
bool gLastCalibrationSequenceSucceeded = false;
bool gDisableRunTimingConstants = false;
double gLastRunMeanGoodLe = std::numeric_limits<double>::quiet_NaN();
double gLastRunMeanECalAdcTime = std::numeric_limits<double>::quiet_NaN();
std::vector<double> gLastRunGroupMeanGoodLe(3,
    std::numeric_limits<double>::quiet_NaN());
std::size_t gLastRunGoodLeCount = 0;
std::size_t gLastRunGoodECalEventCount = 0;
// Most recent combined within-half-bar ECal timing fit.  Exposed so a
// run-specific p1 calibrator can write CDet_run<run>.dat without modifying the
// detector-wide master calibration.
double gLastECalFixedEffectsSlope = std::numeric_limits<double>::quiet_NaN();
double gLastECalFixedEffectsSlopeError = std::numeric_limits<double>::quiet_NaN();
bool gLastECalFixedEffectsValid = false;
// Most recent accepted-pair Gaussian-core timing fit.  These values make the
// final timing-origin definition available to reproducible run-shift drivers.
double gLastPairedCoreMean = std::numeric_limits<double>::quiet_NaN();
double gLastPairedCoreMeanError = std::numeric_limits<double>::quiet_NaN();
double gLastPairedCoreSigma = std::numeric_limits<double>::quiet_NaN();
double gLastPairedCoreSigmaError = std::numeric_limits<double>::quiet_NaN();
double gLastPairedTimeEntries = 0.0;
bool gLastPairedCoreFitValid = false;

// Per-run timing overrides, kept separate from the main detector calibration file.
// Expected file name: CDet_run<RunNumber>.dat
double gGlobalTimingShift = 0.0;   // ns, additive final timing shift
bool   gGlobalTimingLoaded = false;
bool   gRunECalP0Loaded = false;
bool   gRunECalP1Loaded = false;
double gRunECalTimeMin = 10.0;
double gRunECalTimeMax = 35.0;
bool   gRunECalTimeWindowLoaded = false;
bool   gConfigurationSelectionAuthoritative = false;
std::string gRunTimingFile = "";
// When enabled, absolute pair-timing cuts are evaluated before applying the
// run origin shift.  The reported times still include that shift.  This makes
// s_run a pure coordinate translation for its calibration/closure workflow.
bool gShiftInvariantPairTimingCuts = false;

bool LoadCalibrationConstants(const std::string& fname);
bool LoadRunTimingConstants(const std::string& fname);
bool WriteCalibrationConstants(const std::string& fname);

inline double GetPixelToffsetCorr(int elID) {
  if (0 <= elID && elID < NumCDetPaddles && (int)gPixelToffsetCorr.size() == NumCDetPaddles) return gPixelToffsetCorr[elID];
  return 0.0;
}


enum ECalibrationStage {
  kStage0_NoCalibration            = 0,
  kStage1_FitPixelOffsets          = 1,
  kStage2_ApplyPixelOffsets        = 2,
  kStage3_FitECalTiming            = 3,
  kStage4_ApplyPixelOffsetsECal    = 4,
  kStage5_InspectTimeWalk          = 5,
  kStage6_FitTimeWalk              = 6,
  kStage7_ApplyAllCorrections      = 7,
  kStage8_FitECalTiming_TimeWalk   = 8
};

inline bool StageUsesPixelOffsets(int stage) {
  return stage >= kStage2_ApplyPixelOffsets;
}

inline bool StageUsesECalCorrection(int stage) {
  return stage >= kStage4_ApplyPixelOffsetsECal && stage != kStage8_FitECalTiming_TimeWalk;
}

inline bool StageUsesTimeWalkCorrection(int stage) {
  return stage >= kStage7_ApplyAllCorrections || stage == kStage8_FitECalTiming_TimeWalk;
}

inline const char* CalibrationStageDescription(int stage) {
  switch (stage) {
    case kStage0_NoCalibration:       return "0 = no calibration";
    case kStage1_FitPixelOffsets:       return "1 = fit pixel time offsets from CDet LE";
    case kStage2_ApplyPixelOffsets:     return "2 = apply pixel time offsets";
    case kStage3_FitECalTiming:         return "3 = fit CDet/ECal timing using pixel offsets";
    case kStage4_ApplyPixelOffsetsECal: return "4 = apply pixel offsets + ECal timing correction";
    case kStage5_InspectTimeWalk:       return "5 = inspect time-walk using pixel offsets + ECal timing correction";
    case kStage6_FitTimeWalk:           return "6 = fit time-walk using pixel offsets + ECal timing correction";
    case kStage7_ApplyAllCorrections:   return "7 = apply pixel offsets + ECal timing correction + time-walk correction";
    case kStage8_FitECalTiming_TimeWalk: return "8 = fit CDet/ECal timing using pixel offsets + time-walk correction (no ECal correction)";
    default:                          return "unknown stage";
  }
}

inline void ConfigureCalibrationStage(int stage) {
  LoadCalibrationConstants(gCalibrationFile);

  const bool usePixelOffsets =
      StageUsesPixelOffsets(stage) || (stage == kStage1_FitPixelOffsets && gPixelToffsetLoaded);
  const bool useECalCorrection =
      StageUsesECalCorrection(stage) || (stage == kStage3_FitECalTiming && gECalParamsLoaded);
  const bool useTimeWalkCorrection =
      StageUsesTimeWalkCorrection(stage) || (stage == kStage6_FitTimeWalk && gTimeWalkParamsLoaded);

  gUseECalTimeCorr = useECalCorrection;
  gUseTimeWalkCorr = useTimeWalkCorrection;

  if (!usePixelOffsets) {
    gPixelToffsetCorr.assign(NumCDetPaddles, 0.0);
    gPixelToffsetNhits.assign(NumCDetPaddles, 0);
    std::cout << "[CDet] Stage " << stage
              << ": pixel offsets disabled; proceeding with zero pixel offsets.\n";
  }

  std::cout << "[CDet] Calibration stage " << stage << ": "
            << CalibrationStageDescription(stage) << "\n"
            << "        applyPixelOffsets=" << (usePixelOffsets ? "true" : "false")
            << "  applyECalCorr=" << (gUseECalTimeCorr ? "true" : "false")
            << "  applyTimeWalk=" << (gUseTimeWalkCorr ? "true" : "false")
            << std::endl;
}


bool LoadRunTimingConstants(const std::string& fname) {
  gGlobalTimingShift = 0.0;
  gGlobalTimingLoaded = false;
  gRunECalP0Loaded = false;
  gRunECalP1Loaded = false;
  gRunECalTimeMin = 10.0;
  gRunECalTimeMax = 35.0;
  gRunECalTimeWindowLoaded = false;
  bool runECalTimeMinSeen = false;
  bool runECalTimeMaxSeen = false;

  std::ifstream fin(fname.c_str());
  if (!fin) {
    std::cout << "[CDet] No run-timing file '" << fname
              << "' found; using global timing shift = 0 ns.\n";
    return false;
  }

  std::string line, section;
  while (std::getline(fin, line)) {
    const std::string::size_type first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos) continue;
    line.erase(0, first);
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    if (line[0] == '[') {
      section = line;
      continue;
    }

    std::istringstream iss(line);

    if (section == "[GlobalTiming]" || section == "[ECalTiming]" ||
        section == "[ECalSelection]") {
      std::string key;
      if (!(iss >> key)) continue;

      if (iss >> std::ws && iss.peek() == '=') {
        iss.get(); // accept optional '='
      }
      double val = 0.0;
      if (!(iss >> val)) continue;

      if (section == "[GlobalTiming]" && key == "shift_ns") {
          gGlobalTimingShift = val;
          gGlobalTimingLoaded = true;
      } else if (section == "[ECalTiming]" && key == "p0") {
        gECalFitP0 = val;
        gRunECalP0Loaded = true;
      } else if (section == "[ECalTiming]" && key == "p1") {
        gECalFitP1 = val;
        gRunECalP1Loaded = true;
      } else if (section == "[ECalSelection]" && key == "time_min") {
        gRunECalTimeMin = val;
        runECalTimeMinSeen = true;
      } else if (section == "[ECalSelection]" && key == "time_max") {
        gRunECalTimeMax = val;
        runECalTimeMaxSeen = true;
      }
    }
  }

  gRunECalTimeWindowLoaded = runECalTimeMinSeen && runECalTimeMaxSeen &&
      std::isfinite(gRunECalTimeMin) &&
      std::isfinite(gRunECalTimeMax) && gRunECalTimeMin < gRunECalTimeMax;
  if (runECalTimeMinSeen != runECalTimeMaxSeen ||
      ((runECalTimeMinSeen && runECalTimeMaxSeen) && !gRunECalTimeWindowLoaded)) {
    std::cerr << "[CDet] WARNING: invalid or incomplete [ECalSelection] in '"
              << fname << "'; retaining the caller's ECal timing window.\n";
  }

  if (gGlobalTimingLoaded) {
    std::cout << "[CDet] Loaded run timing shift from '" << fname
              << "': shift_ns = " << gGlobalTimingShift << "\n";
  }
  if (gRunECalP0Loaded || gRunECalP1Loaded) {
    std::cout << "[CDet] Loaded run-specific ECal timing override from '" << fname
              << "': p0=" << gECalFitP0 << " p1=" << gECalFitP1 << "\n";
  }
  if (gRunECalTimeWindowLoaded && !gConfigurationSelectionAuthoritative) {
    std::cout << "[CDet] Loaded run-specific ECal selection window from '"
              << fname << "': " << gRunECalTimeMin << " < time < "
              << gRunECalTimeMax << " ns\n";
  }

  return gGlobalTimingLoaded || gRunECalP0Loaded || gRunECalP1Loaded ||
      gRunECalTimeWindowLoaded;
}

bool LoadCalibrationConstants(const std::string& fname) {
  gCalibrationLoaded = false;
  gPixelToffsetCorr.assign(NumCDetPaddles, 0.0);
  gPixelToffsetNhits.assign(NumCDetPaddles, 0);
  gPixelToffsetLoaded = false;
  gECalParamsLoaded = false;
  gECalDeltaLoaded = false;
  gTimeWalkParamsLoaded = false;

  std::ifstream fin(fname.c_str());
  if (!fin) {
    std::cout << "[CDet] No calibration file '" << fname << "' found; proceeding with defaults.\n";
    return false;
  }

  std::string line, section;
  std::vector<bool> pixelSeen(NumCDetPaddles, false);
  std::vector<double> pixelCorr(NumCDetPaddles, 0.0);
  std::vector<int> pixelNhits(NumCDetPaddles, 0);
  double ecalP0 = gECalFitP0, ecalP1 = gECalFitP1;
  double ecalDelta = gECalDeltaShift;
  double twP1L1 = gTimeWalkP1_L1, twP1L2 = gTimeWalkP1_L2;
  double twRefL1 = gTimeWalkTotRef_L1, twRefL2 = gTimeWalkTotRef_L2;
  double twMin = gTimeWalkTotMin, twMax = gTimeWalkTotMax;
  bool ecalP0Seen = false, ecalP1Seen = false, ecalDeltaSeen = false;
  bool twP1L1Seen = false, twP1L2Seen = false;
  bool twRefL1Seen = false, twRefL2Seen = false;
  bool twMinSeen = false, twMaxSeen = false;
  bool legacyBarOffsetsSeen = false;

  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    if (line[0] == '[') {
      section = line;
      if (section == "[BarOffsets]") legacyBarOffsetsSeen = true;
      continue;
    }

    std::istringstream iss(line);

    if (section == "[PixelOffsets]") {
      int pixelID = -1;
      double dt = 0.0;
      int nhits = 0;
      if (!(iss >> pixelID >> dt)) continue;
      if (!(iss >> nhits)) nhits = 0;
      if (pixelID >= 0 && pixelID < NumCDetPaddles) {
        pixelCorr[pixelID] = dt;
        pixelNhits[pixelID] = nhits;
        pixelSeen[pixelID] = true;
      }
    } else if (section == "[ECalTiming]") {
      std::string key;
      double val = 0.0;
      if (!(iss >> key >> val)) continue;
      if (key == "p0") { ecalP0 = val; ecalP0Seen = true; }
      else if (key == "p1") { ecalP1 = val; ecalP1Seen = true; }
      else if (key == "delta") { ecalDelta = val; ecalDeltaSeen = true; }
    } else if (section == "[TimeWalk]") {
      std::string key;
      double val = 0.0;
      if (!(iss >> key >> val)) continue;
      if (key == "p1_L1") { twP1L1 = val; twP1L1Seen = true; }
      else if (key == "p1_L2") { twP1L2 = val; twP1L2Seen = true; }
      else if (key == "totref_L1") { twRefL1 = val; twRefL1Seen = true; }
      else if (key == "totref_L2") { twRefL2 = val; twRefL2Seen = true; }
      else if (key == "totmin") { twMin = val; twMinSeen = true; }
      else if (key == "totmax") { twMax = val; twMaxSeen = true; }
    }
  }

  const int nPixelRead = std::count(pixelSeen.begin(), pixelSeen.end(), true);
  gPixelToffsetLoaded = (nPixelRead == NumCDetPaddles);
  gECalParamsLoaded = ecalP0Seen && ecalP1Seen;
  gECalDeltaLoaded = ecalDeltaSeen;
  gTimeWalkParamsLoaded = twP1L1Seen && twP1L2Seen && twRefL1Seen &&
                          twRefL2Seen && twMinSeen && twMaxSeen;
  if (gPixelToffsetLoaded) {
    gPixelToffsetCorr.swap(pixelCorr);
    gPixelToffsetNhits.swap(pixelNhits);
  }
  if (gECalParamsLoaded) {
    gECalFitP0 = ecalP0;
    gECalFitP1 = ecalP1;
  }
  if (gECalDeltaLoaded) gECalDeltaShift = ecalDelta;
  if (gTimeWalkParamsLoaded) {
    gTimeWalkP1_L1 = twP1L1;
    gTimeWalkP1_L2 = twP1L2;
    gTimeWalkTotRef_L1 = twRefL1;
    gTimeWalkTotRef_L2 = twRefL2;
    gTimeWalkTotMin = twMin;
    gTimeWalkTotMax = twMax;
  }
  gCalibrationLoaded = gPixelToffsetLoaded || gECalParamsLoaded || gTimeWalkParamsLoaded;

  if (nPixelRead > 0 && !gPixelToffsetLoaded)
    std::cerr << "[CDet] WARNING: ignoring incomplete [PixelOffsets] section ("
              << nPixelRead << "/" << NumCDetPaddles << ") in '" << fname << "'.\n";
  if (legacyBarOffsetsSeen)
    std::cerr << "[CDet] WARNING: ignoring legacy [BarOffsets] section in '" << fname << "'; regenerate this calibration to produce [PixelOffsets].\n";
  if ((ecalP0Seen || ecalP1Seen) && !gECalParamsLoaded)
    std::cerr << "[CDet] WARNING: ignoring incomplete [ECalTiming] section in '"
              << fname << "'.\n";
  if ((twP1L1Seen || twP1L2Seen || twRefL1Seen || twRefL2Seen || twMinSeen || twMaxSeen) &&
      !gTimeWalkParamsLoaded)
    std::cerr << "[CDet] WARNING: ignoring incomplete [TimeWalk] section in '"
              << fname << "'.\n";

  if (gPixelToffsetLoaded) {
    std::cout << "[CDet] Loaded " << nPixelRead << " pixel offsets from '" << fname << "'.\n";
  }
  if (gECalParamsLoaded) {
    std::cout << "[CDet] Loaded ECal timing parameters from '" << fname
              << "': p0=" << gECalFitP0 << " p1=" << gECalFitP1;
    if (gECalDeltaLoaded)
      std::cout << " delta=" << gECalDeltaShift << " ns (fixed)";
    else
      std::cout << "\n[CDet] WARNING: no ECalTiming delta is stored; this legacy file "
                   "will use sample-dependent recentering until it is rewritten";
    std::cout << "\n";
  }
  if (gTimeWalkParamsLoaded) {
    std::cout << "[CDet] Loaded time-walk parameters from '" << fname
              << "': p1_L1=" << gTimeWalkP1_L1
              << " p1_L2=" << gTimeWalkP1_L2 << "\n";
  }

  return gCalibrationLoaded;
}

bool WriteCalibrationConstants(const std::string& fname) {
  if (gPixelToffsetCorr.size() != NumCDetPaddles || gPixelToffsetNhits.size() != NumCDetPaddles) {
    std::cerr << "[CDet] ERROR: refusing to write incomplete pixel-offset vectors.\n";
    return false;
  }
  const std::string tmpname = fname + ".tmp";
  std::ofstream fout(tmpname.c_str(), std::ios::out | std::ios::trunc);
  if (!fout) {
    std::cout << "[CDet] ERROR: could not write calibration file '" << fname << "'\n";
    return false;
  }

  fout << "# CDet master calibration constants\n"
       << "# run " << gRunNumber << "\n"
       << "# calibration_stage " << gCalibrationStage << "\n"
       << "# events_processed " << gNumEventsInRun << "\n\n";

  fout << "[PixelOffsets]\n";
  fout.setf(std::ios::fixed);
  fout.precision(6);
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    fout << pixelID << " " << gPixelToffsetCorr[pixelID] << " " << gPixelToffsetNhits[pixelID] << "\n";
  }
  fout << "\n";

  fout << "[ECalTiming]\n";
  fout << "p0 " << gECalFitP0 << "\n";
  fout << "p1 " << gECalFitP1 << "\n";
  fout << "delta " << gECalDeltaShift << "\n\n";

  fout << "[TimeWalk]\n";
  fout << "p1_L1 " << gTimeWalkP1_L1 << "\n";
  fout << "p1_L2 " << gTimeWalkP1_L2 << "\n";
  fout << "totref_L1 " << gTimeWalkTotRef_L1 << "\n";
  fout << "totref_L2 " << gTimeWalkTotRef_L2 << "\n";
  fout << "totmin " << gTimeWalkTotMin << "\n";
  fout << "totmax " << gTimeWalkTotMax << "\n";

  fout.flush();
  if (!fout) {
    std::cerr << "[CDet] ERROR: failed while writing '" << tmpname << "'.\n";
    fout.close();
    std::remove(tmpname.c_str());
    return false;
  }
  fout.close();
  if (std::rename(tmpname.c_str(), fname.c_str()) != 0) {
    std::cerr << "[CDet] ERROR: could not replace calibration file '" << fname << "'.\n";
    std::remove(tmpname.c_str());
    return false;
  }
  std::cout << "[CDet] Wrote calibration constants atomically to '" << fname << "'.\n";
  return true;
}

std::vector<double> vAllGoodTe;
std::vector<double> vAllGoodTot;
std::vector<int> vAllGoodPMT;
std::vector<int> vAllGoodBar;

std::vector<int> vhitCDetPMT;
std::vector<int> vRow;
std::vector<std::vector<int>> vGoodCol;
std::vector<std::vector<int>> vGoodLayer;

std::vector<std::vector<double>> vCDetX;
std::vector<std::vector<double>> vCDetY;
std::vector<std::vector<double>> vCDetZ;

std::vector<std::vector<double>> vCDetGoodX;
std::vector<std::vector<double>> vCDetGoodY;
std::vector<std::vector<double>> vCDetGoodZ;
std::vector<Long64_t> vTreeEntry; // maps saved passing-event index -> TTree entry

std::vector<int> vRowLayer1Side1;
std::vector<int> vRowLayer2Side1;
std::vector<int> vRowLayer1Side2;
std::vector<int> vRowLayer2Side2;

std::vector<int> vnhits1;
std::vector<int> vnpaddles;
std::vector<int> vngoodpaddles;
std::vector<int> vngoodTDCpaddles;
std::vector<int> vngoodhits1;
std::vector<int> vngoodTDChits1;
std::vector<int> vnhits2;
std::vector<int> vngoodhits2;
std::vector<int> vngoodTDChits2;

std::vector<std::vector<double>> vCDetPaddleRawTot;
std::vector<std::vector<double>> vCDetPaddleCutTot;

//2D vectors
std::vector<std::vector<double>> vRawLe;
std::vector<std::vector<double>> vRawTe;
std::vector<std::vector<double>> vRawTot;
std::vector<std::vector<int>> vRawID;

std::vector<std::vector<double>> vGoodLe;
std::vector<std::vector<double>> vGoodTe;
std::vector<std::vector<double>> vGoodTot;
std::vector<std::vector<int>> vGoodID;

std::vector<std::vector<double>> vTestLe;
std::vector<std::vector<double>> vTestTe;
std::vector<std::vector<double>> vTestTot;
std::vector<std::vector<int>> vTestID;

// per-event 1D
std::vector<double> vCDetMultAll;          // like hMultiplicity
std::vector<int>    vCDetMultAllPMT;       // ids belonging to those mult values

// 2D vectors for all ECal cluster arrays
std::vector<std::vector<double>> v_ECal_clus_adctime;
std::vector<std::vector<double>> v_ECal_clus_again;
std::vector<std::vector<double>> v_ECal_clus_atimeblk;
std::vector<std::vector<double>> v_ECal_clus_col;
std::vector<std::vector<double>> v_ECal_clus_col_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_e;
std::vector<std::vector<double>> v_ECal_clus_e_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_eblk;
std::vector<std::vector<double>> v_ECal_clus_eblk_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_id;
std::vector<std::vector<double>> v_ECal_clus_id_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_nblk;
std::vector<std::vector<double>> v_ECal_clus_nblk_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_row;
std::vector<std::vector<double>> v_ECal_clus_row_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_tdctime;
std::vector<std::vector<double>> v_ECal_clus_tdctime_tw;
std::vector<std::vector<double>> v_ECal_clus_tdctimeblk;
std::vector<std::vector<double>> v_ECal_clus_tdctimeblk_tw;
std::vector<std::vector<double>> v_ECal_clus_x;
std::vector<std::vector<double>> v_ECal_clus_y;
std::vector<std::vector<double>> v_ECal_a_p;
std::vector<std::vector<double>> v_ECal_a_amp_p;
std::vector<std::vector<double>> v_ECal_a_time;
std::vector<std::vector<double>> v_ECal_adcxpos;

// 1D (event-wise) vectors ECal for scalars:
std::vector<double> v_ECal_nclus;
std::vector<double> v_ECalX;
std::vector<double> v_ECalY;
std::vector<double> v_ECalE;
std::vector<double> v_ECalAdcTime;
//1D Good ECal Events
std::vector<double> v_GoodECalX;
std::vector<double> v_GoodECalY;
std::vector<double> v_GoodECalE;
std::vector<double> v_GoodECalAdcTime;

//1D Good HCal Events
std::vector<double> v_GoodHCalAdcTime;
std::vector<double> v_GoodHCalE;

struct CDetHit {
  int    id;     // pixelID
  double le_ns;  // LE
  double tot_ns; // TOT
  double te_ns;  // TE
};

struct AdjPair {
  int event;      // event number/index
  int id1, id2;   // id1 < id2
  double le1, te1, tot1;
  double le2, te2, tot2;
  int i1, i2;     // hit indices within that pixel for this event (for dedupe)
};

static std::vector<AdjPair> vAdjPairs;

// using Hit = std::pair<int,double>; // (pixelID, tot_ns)
std::vector<std::vector<CDetHit>> vEventHits; // [event][hit]
std::vector<std::vector<CDetHit>> vGoodEventHits;

std::vector<int> rawHitCount(2688, 0);
int occupancyEventCount = 0;
std::vector<double> rawHitOccupancy(2688, 0.0);
std::vector<double> rawSinglesRateHz(2688, 0.0);
std::vector<double> rawSinglesRateErrorHz(2688, 0.0);
TH1D* hRawSinglesRateVsID = nullptr;
std::vector<double> ecalClusterEnergiesGeV;
std::vector<double> ecalClusterEnergiesGeVInTime;
int ecalClusterProcessedEventCount = 0;
int ecalClusterCountMismatchEvents = 0;
int ecalClusterNonFiniteEnergyCount = 0;
int ecalClusterNonFiniteEnergyCountInTime = 0;
int ecalClusterGroupIndex = 0;
TH1D* hECalClusterEnergySpectrum = nullptr;
TH1D* hECalClusterEnergySpectrumInTime = nullptr;
std::vector<int> totCutHitCount(2688, 0);
std::vector<double> totCutHitOccupancy(2688, 0.0);
std::vector<double> ave_tot(2688,0);
std::vector<int> vNumRawAdjacentHits;
std::vector<int> vNumGoodAdjacentHits;

//copy a TTreeReaderArray<double> into a std::vector<double>, makes it easier to fill the 2D vector
inline std::vector<double> copyArray(const TTreeReaderArray<double>& arr) {
  std::vector<double> v;
  v.reserve(arr.GetSize());
  for (size_t i = 0; i < arr.GetSize(); ++i) v.push_back(arr[i]);
  return v;
}


// per-event vectors grouped by PMT index (2D)
std::vector<std::vector<double>> vCDetMultPerPMT(nTdc);


/*
namespace TCDet {
  Int_t NdataMult;
  Double_t TDCmult[nTdc*2];

  Int_t NdataRawElID;
  Double_t RawElID[nTdc*2];
  Int_t NdataRawElLE;
  Double_t RawElLE[nTdc*2];
  Int_t NdataRawElTE;
  Double_t RawElTE[nTdc*2];
  Int_t NdataRawElTot;
  Double_t RawElTot[nTdc*2];

  Int_t NdataGoodRow;
  Double_t GoodRow[nTdc*2];
  Int_t NdataGoodCol;
  Double_t GoodCol[nTdc*2];
  Int_t NdataGoodLayer;
  Double_t GoodLayer[nTdc*2];

  Int_t NdataGoodElID;
  Double_t GoodElID[nTdc*2];
  Int_t NdataGoodElLE;
  Double_t GoodElLE[nTdc*2];
  Int_t NdataGoodElTE;
  Double_t GoodElTE[nTdc*2];
  Int_t NdataGoodElTot;
  Double_t GoodElTot[nTdc*2];

  Int_t NdataGoodX;
  Double_t GoodX[nTdc*2];
  Int_t NdataGoodY;
  Double_t GoodY[nTdc*2];
  Int_t NdataGoodZ;
  Double_t GoodZ[nTdc*2];

  Double_t ECalX;
  Double_t ECalY;
  Double_t ECalE;
  Double_t nhits;
  Double_t ngoodhits;
  Double_t ngoodTDChits;



};*/

//===================================================== Globals for paddle hits
Double_t nhits_paddles[nTdc*2];
Double_t ngoodhits_paddles[nTdc*2];
Double_t ngoodTDChits_paddles[nTdc*2];
Double_t npaddles;
Double_t ngoodpaddles;
Double_t ngoodTDCpaddles;

TChain *T = 0;

//===================================================== Histogram Declarations
// number of histo bins
const int NTotBins = 200;
const double TotBinLow = 1.;
const double TotBinHigh = 51.;
const int RefNTotBins = 800;
const double RefTotBinLow = 1.;
const double RefTotBinHigh = 201.;

double TDCBinLow;
double TDCBinHigh;
int NTDCBins;
double RefTDCBinLow;
double RefTDCBinHigh;
int RefNTDCBins;


//const int num_bad = 0;
//

const int num_bad = 5;

const int bad_channels[] = {
	1161, 1472, 1670, 1696, 1896
};

//const int num_bad = 24;
//
//const int bad_channels[] = {
//	61,
//	64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
//	408,417,1286,
//	2124, 2404,2406,2414
//};

//const int num_bad = 63;

//const int bad_channels[] = {
//	35,40,42,45,
//	80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,
//	136,169,183,184,192,193,194,198,199,
//	314,
//	390,391,392,395,397,433,506,507,
//	678,800,902,933,949,
//	1177,1205,1215,1255,
//	1296,1297,1299,1302,1303,1304,1307,1311,
//	1912,2140,
//	2420,2422,2430,
//	2629,2630,2662
//};


/*
const int bad_channels[] = {
		45,   224, 226, 234, 241, 284, 240, 604, 690, 698, 702,
		864, 1174,1437,1488,1726,1912,2131,2132,2213,2214,
		2411,2417,2425,2426,2427,2434,2436,2437,2438,2448, 2457,
		2466,2467,2468,2469,2470,2471,2472,2473,2474,2475,
		2476,2477,2478,2479,2502,2517,2518,2519,2520,2522,
		2592,2607,2613,2614,2616,2620,2624,2628,2632,2633,
		2634,2656,2661,2664,2665,2666,2667,2673,2675,2681,
		2683,2687};
*/

// Raw hits ie all hits
TH1F *hRawLe[nTdc];
TH1F *hRawTe[nTdc];
TH1F *hRawTot[nTdc];
TH1F *hGoodLe[nTdc];
TH1F *hGoodTe[nTdc];
TH1F *hGoodTot[nTdc];

TH1F *hAllRawLe;
TH1F *hAllRawTe;
TH1F *hAllRawTot;
TH1F *hAllRawPMT;
TH1F *hAllRawBar;
TH1F *hAllGoodLe;
TH1F *hAllGoodTe;
TH1F *hAllGoodTot;
TH1F *hAllGoodPMT;
TH1F *hAllGoodBar;

TH2D* hGoodLeVsTot_L1;
TH2D* hGoodLeVsTot_L2;

TH2F *h2AllGoodLe;
TH2F *h2AllGoodTe;
TH2F *h2AllGoodTot;

TH2F *h2TDCTOTvsLE;
TH2F *h2CDetX1vsX2;

TH2F *h2TOTvsXDiff1;
TH2F *h2TOTvsXDiff2;
TH2F *h2LEvsXDiff1;
TH2F *h2LEvsXDiff2;

TH2F *hBarRateHV;

TH1F *hRefRawLe;
TH1F *hRefRawTe;
TH1F *hRefRawTot;
TH1F *hRefRawPMT;
TH1F *hRefGoodLe;
TH1F *hRefGoodTe;
TH1F *hRefGoodTot;
TH1F *hRefGoodPMT;

TH1F *hMultiplicityL[nTdc];
TH1F *hMultiplicity;


// hit channel id
TH1F *hHitPMT;
TH1F *hRow;
TH1F *hRowLayer1Side1;
TH1F *hRowLayer1Side2;
TH1F *hRowLayer2Side1;
TH1F *hRowLayer2Side2;
TH1F *hLayer;
TH1F *hCol;

TH1F *hnhits1;
TH1F *hngoodhits1;
TH1F *hngoodTDChits1;
TH1F *hnhits2;
TH1F *hngoodhits2;
TH1F *hngoodTDChits2;
TH1F *hnhits_ev;
TH1F *hngoodhits_ev;
TH1F *hngoodTDChits_ev;

TH1F *hnpaddles;
TH1F *hngoodpaddles;
TH1F *hngoodTDCpaddles;

TH1F *hHitX;
TH1F *hHitY;
TH1F *hHitZ;

TH2F *hHitXY1;
TH2F *hHitXY2;

TH1F *hXECal;
TH1F *hYECal;
TH1F *hEECal;

TH2F *hXECalCDet1;
TH2F *hXECalCDet2;
TH2F *hXECalCDet1_min;
TH2F *hXECalCDet2_min;
TH2F *hYECalCDet1;
TH2F *hYECalCDet2;
TH2F *hEECalCDet1;
TH2F *hEECalCDet2;

TH1F *hXDiffECalCDet1;
TH1F *hXPlusECalCDet1;
TH1F *hXDiffECalCDet2;
TH1F *hXPlusECalCDet2;
TProfile *pXDiffECalCDet1VsXCDet1;
TProfile *pXDiffECalCDet2VsXCDet2;

TH2F *hXCDet1CDet2;

// 2D histograms
TH2F* h2d_RawLE;
TH2F* h2d_RawTE;
TH2F* h2d_RawTot;

TH2F* h2d_GoodLE;
TH2F* h2d_GoodTE;
TH2F* h2d_GoodTot;

TH2F* h2d_Mult;

using namespace std;

bool check_bad(int pmt, bool suppress_bad) {
	bool flag = false;
	if (!suppress_bad) return flag;
	for (int i=0;i<num_bad;i++) {
		if (pmt == bad_channels[i])  flag = true;
	}
	return flag;
}

std::vector<double> extractBinContents(const TH1* hist) {
    if (!hist) {
        throw std::invalid_argument("Null histogram pointer passed.");
    }

    int nBins = hist->GetNbinsX();
    std::vector<double> contents;

    // Loop over all *visible* bins (skip underflow bin 0 and overflow bin nBins+1)
    for (int i = 1; i <= nBins; ++i) {
        contents.push_back(hist->GetBinContent(i));
    }

    return contents;
}


vector<vector<double>> readDataFromFiles(const vector<string>& filenames) {
    const int NUM_VALUES = 42;
    vector<vector<double>> allData;

    for (const auto& filename : filenames) {
        ifstream infile(filename);
        if (!infile) {
            cerr << "Error opening file: " << filename << endl;
            continue; // Skip this file and move on
        }

        vector<double> fileData;
        double value;
        for (int i = 0; i < NUM_VALUES; ++i) {
            infile >> value;
            if (!infile) {
                cerr << "Error reading value " << i << " from file: " << filename << endl;
                break; // Stop reading this file
            }
            fileData.push_back(value);
        }

        if (fileData.size() == NUM_VALUES) {
            allData.push_back(fileData);
        } else {
            cerr << "Incomplete data in file: " << filename << endl;
        }

        infile.close();
    }

    return allData;
}
void pixelsFromBar(int bar)
{
    if (bar < 0) {
        std::cerr << "Error: bar number must be >= 0" << std::endl;
        return;
    }

    const int pixelsPerBar = 16;

    int firstPixel = bar * pixelsPerBar;
    int lastPixel  = firstPixel + pixelsPerBar - 1;

    std::cout << "Bar/PMT " << bar << " has pixel slots: "
              << firstPixel << " to " << lastPixel << std::endl;

    std::cout << "Pixels: ";
    for (int pix = firstPixel; pix <= lastPixel; pix++) {
        std::cout << pix;
        if (pix < lastPixel) std::cout << ", ";
    }
    std::cout << std::endl;
}
int getPixelID(int layerNum, int sideNum, int submoduleNum, int pmtNum, int pixelNum){
  // Calculate paddle number, note that missing pixels are included here
  // Validate inputs
  if (layerNum < 1 || layerNum > 2 ||
    sideNum < 1 || sideNum > 2 ||
    submoduleNum < 1 || submoduleNum > 3 ||
    pmtNum < 1 || pmtNum > 14 ||
    pixelNum < 1 || pixelNum > 16) {
    std::cerr << "Error: Invalid input values.\n"
              << "  layerNum must be 1 or 2\n"
              << "  sideNum must be 1 or 2\n"
              << "  submoduleNum must be 1 to 3\n"
              << "  pmtNum must be 1 to 14\n"
              << "  pixelNum must be 1 to 16\n";
    return -1;  // Error code
}

  int pixel = (layerNum - 1) * 1344 + //1344 pixels per layer, 0-1343 in layer 1, 1344-2687 in layer 2
  (submoduleNum - 1) * 224 + //224 pixels per side of a module
  (sideNum - 1) * 672 + //672 pixels per side
  (pmtNum-1) * 16 + //16 pixels per pmt
  pixelNum;
  return pixel - 1;
}

std::vector<int> getLocation(int pixelID) {
  // Check valid range
  if (pixelID < 0 || pixelID > 2687) {
      std::cerr << "Error: pixelID must be in the range 0 to 2687.\n";
      return {};  // return empty vector to signal error
  }

  int layerNum      = pixelID / 1344; //1344 pixels per layer
  pixelID               %= 1344;

  int sideNum       = pixelID / 672;  //672 pixels per side
  pixelID               %= 672;

  int submoduleNum  = pixelID / 224; //224 pixels per side of module
  pixelID               %= 224;

  int pmtNum        = pixelID / 16; //16 pixels per bar
  pixelID               %= 16;
  int pixelNum      = pixelID % 16;

  return {layerNum, sideNum, submoduleNum, pmtNum, pixelNum};
}
//Used to fill 2D arrays
template<typename T>
std::vector<T> fill2D(const TTreeReaderArray<T>& arr) {
  std::vector<T> tmp;
  int n = arr.GetSize();
  tmp.reserve(n);
  for (int i=0; i<n; i++) tmp.push_back(arr[i]);
  return tmp;
}

void PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(Int_t RunNumber1=5811, Int_t nevents=50000, Int_t calibStage = 7,
  Int_t groupIndex = 0, Int_t minSeg = -1, Int_t maxSeg = -1,
	Double_t LeMin = 0.02, Double_t LeMax = 100,
	Double_t TotMin = 0.02, Double_t TotMax = 150.0,
  Double_t ECalMinT = 10.0, Double_t ECalMaxT = 35.0,
	Int_t nhitcutlow1 = 1, Int_t nhitcuthigh1 = 100,
	Int_t nhitcutlow2 = 0, Int_t nhitcuthigh2 = 100,
	Double_t XDiffCut = 0.05, Double_t XOffset = 0.0, Double_t YOffset = 0.1,
        Int_t layer_choice=3,
	bool suppress_bad = false,
	Int_t nruns=30, Int_t maxstream = 2, Int_t firstevent = 1,
        bool useReferenceTiming = false,
        Double_t ECalMinE = -1.0e9, Double_t ECalMaxE = 1.0e9)
{
  // Each invocation owns a fresh analysis sample. Without this reset, calling
  // the main macro twice in one ROOT process appends duplicate events to the
  // global vectors and corrupts fit counts and uncertainties.
  ResetCalibrationGlobals();
  gLastCalibrationStageSucceeded = false;
  gRunNumber = RunNumber1;
  gCalibrationStage = calibStage;
  gNumEventsInRun = 0;
  gLastRunMeanGoodLe = std::numeric_limits<double>::quiet_NaN();
  gLastRunMeanECalAdcTime = std::numeric_limits<double>::quiet_NaN();
  gLastRunGroupMeanGoodLe.assign(3,
      std::numeric_limits<double>::quiet_NaN());
  gLastRunGoodLeCount = 0;
  gLastRunGoodECalEventCount = 0;
  ecalClusterGroupIndex = groupIndex;
  ecalClusterEnergiesGeV.clear();
  ecalClusterProcessedEventCount = 0;
  ecalClusterCountMismatchEvents = 0;
  ecalClusterNonFiniteEnergyCount = 0;
  (void)firstevent; // currently unused
  gCalibrationFile = RunNumber1 == 0
      ? TString::Format("CDet_calibration_dt_group%d.dat", groupIndex).Data()
      : "CDet_calibration_dt.dat";
  const char *runListFile = "crossRuns.txt"; //file with run numbers to analyze
  Int_t nseg = nruns/(maxstream+1);
	Double_t RefLeMin = 0.02;
	Double_t RefLeMax = 251.0;
	RefNTDCBins = (RefLeMax-RefLeMin)/4;
	Double_t RefTotMin = 0.02;
	Double_t RefTotMax = 251.0;

	NTDCBins = 2*(LeMax-LeMin)/.0160167; // 4 ns is the trigger time, 0.018 ns is the expected time resolution, if we use a reference TDC ?
					// 4 ns resolution is the best we can hope for, I think, using only the module trigger time.

	NXDiffBins = (int)((2*XDiffCut)/0.0073);
	XDiffLow = XOffset-XDiffCut;
	XDiffHigh = XOffset+XDiffCut;

  // InFile is the input file without absolute path and without .root suffix
  // nevents is how many events to analyse, -1 for all

  // To execute
  // root -l
  // .L PlotRawTDC2D.C+
  // PlotRawTDC2D("filename", -1)
  TDCBinLow = LeMin;
  TDCBinHigh = LeMax;
  RefTDCBinLow = RefLeMin;
  RefTDCBinHigh = RefLeMax;


// Configure which calibrations are active for this run.
ConfigureCalibrationStage(calibStage);
  gAnalysisPixelToffsetCorr = gPixelToffsetCorr;
  gAnalysisPixelToffsetNhits = gPixelToffsetNhits;
  gAnalysisPixelToffsetSnapshotValid =
      (int)gAnalysisPixelToffsetCorr.size() == NumCDetPaddles;

// A single shift is meaningful only for single-run chains. Applying a
// CDet_run0.dat shift to a mixed-run group would silently bias every run.
if (RunNumber1 != 0 && !gDisableRunTimingConstants) {
  gRunTimingFile = Form("CDet_run%d.dat", (int)RunNumber1);
  LoadRunTimingConstants(gRunTimingFile);
  if (gRunECalTimeWindowLoaded) {
    ECalMinT = gRunECalTimeMin;
    ECalMaxT = gRunECalTimeMax;
  }
} else {
  gRunTimingFile.clear();
  gGlobalTimingShift = 0.0;
  gGlobalTimingLoaded = false;
  gRunECalTimeMin = 10.0;
  gRunECalTimeMax = 35.0;
  gRunECalTimeWindowLoaded = false;
  std::cout << "[CDet] Per-run timing constants are disabled"
            << (RunNumber1 == 0 ? " in group mode" : " for this analysis")
            << ".\n";
}

std::cout << "[CDet] Reference timing subtraction is "
          << (useReferenceTiming ? "enabled" : "disabled") << ".\n";


  // hit channel id

  hHitPMT = new TH1F("hHitPMT","hHitPMT",nTdc,0,nTdc);

  hnhits1 = new TH1F("hnhits1","hnhits1",150,1,151);
  hngoodhits1 = new TH1F("hngoothits1","hngoodhits1",75,1,76);
  hngoodTDChits1 = new TH1F("hngoodTDChits1","hngoodTDChits1",75,1,76);

  hnhits2 = new TH1F("hnhits2","hnhits2",150,1,151);
  hngoodhits2 = new TH1F("hngoothits2","hngoodhits2",75,1,76);
  hngoodTDChits2 = new TH1F("hngoodTDChits2","hngoodTDChits2",75,1,76);

  hnhits_ev = new TH1F("hnhits_ev","hnhits_ev",500,0,50000);
  hngoodhits_ev = new TH1F("hngoothits_ev","hngoodhits_ev",500,0,50000);
  hngoodTDChits_ev = new TH1F("hngoodTDChits_ev","hngoodTDChits_ev",500,0,50000);

  hnpaddles = new TH1F("hnpaddles","hnpaddles",200,1,201);
  hngoodpaddles = new TH1F("hngoodpaddles","hnTDCpaddles",75,1,71);
  hngoodTDCpaddles = new TH1F("hngoodTDCpaddles","hngoodTDCpaddles",75,1,76);

  hRow = new TH1F("RowNumber","RowNumber",680, 0, 680);
  hRowLayer1Side1 = new TH1F("RowNumberL1S1","RowNumberL1S1",680, 0, 680);
  hRowLayer1Side2 = new TH1F("RowNumberL1S2","RowNumberL1S2",680, 0, 680);
  hRowLayer2Side1 = new TH1F("RowNumberL2S1","RowNumberL2S1",680, 0, 680);
  hRowLayer2Side2 = new TH1F("RowNumberL2S2","RowNumberL2S2",680, 0, 680);
  hLayer = new TH1F("LayerNumber","LayerNumber",3, 0, 3);
  hCol = new TH1F("ColNumber","ColNumber",3, 0, 3);

  hHitX = new TH1F("HitXposition","HitXPosition",1000,-2.0,2.0);
  hHitY = new TH1F("HitYposition","HitYPosition",200,-0.5,0.5);
  hHitZ = new TH1F("HitZposition","HitZPosition",200,7.5,8.0);

  hHitXY1 = new TH2F("HitXY1position","HitXY1Position",200,-0.5,0.5,800,-2.0,2.0);
  hHitXY2 = new TH2F("HitXY2position","HitXY2Position",200,-0.5,0.5,800,-2.0,2.0);

  hXECal = new TH1F("XECal","XECal",200,-1.5,1.5);
  hYECal = new TH1F("YECal","YECal",200,-1.0,1.0);
  hEECal = new TH1F("EECal","EECal",200,0.0,20.0);

  hXECalCDet1 = new TH2F("XECalCDet1","XECalCDet1",100,-2.0,2.0,100,-2.0,2.0);
  hXECalCDet2 = new TH2F("XECalCDet2","XECalCDet2",100,-2.0,2.0,100,-2.0,2.0);
  hXECalCDet1_min = new TH2F("XECalCDet1_min","XECalCDet1_min (min |x_{CDet}-x_{ECal->CDet}| per event)",100,-2.0,2.0,100,-2.0,2.0);
  hXECalCDet2_min = new TH2F("XECalCDet2_min","XECalCDet2_min (min |x_{CDet}-x_{ECal->CDet}| per event)",100,-2.0,2.0,100,-2.0,2.0);
  hYECalCDet1 = new TH2F("YECalCDet1","YECalCDet1",100,-1.0,1.0,9,-1.0,1.0);
  hYECalCDet2 = new TH2F("YECalCDet2","YECalCDet2",100,-1.0,1.0,9,-1.0,1.0);
  hEECalCDet1 = new TH2F("EECalCDet1","EECalCDet1",100,0.0,20.0,100,-2.0,2.0);
  hEECalCDet2 = new TH2F("EECalCDet2","EECalCDet2",100,0.0,20.0,100,-2.0,2.0);

  hXDiffECalCDet1 = new TH1F("XDiffECalCDet1","XDiffECalCDet1",NXDiffBins,XDiffLow,XDiffHigh);
  hXPlusECalCDet1 = new TH1F("XPlusECalCDet1","XPlusECalCDet1",NXDiffBins,XDiffLow,XDiffHigh);
  hXDiffECalCDet2 = new TH1F("XDiffECalCDet2","XDiffECalCDet2",NXDiffBins,XDiffLow,XDiffHigh);
  hXPlusECalCDet2 = new TH1F("XPlusECalCDet2","XPlusECalCDet2",NXDiffBins,XDiffLow,XDiffHigh);
  pXDiffECalCDet1VsXCDet1 = new TProfile(
      "pXDiffECalCDet1VsXCDet1",
      "Layer 1 mean ECal-CDet x residual vs corrected CDet x;corrected CDet Layer 1 x (m);#LT x_{CDet}-x_{ECal#rightarrowCDet} #GT (m)",
      200, -2.0, 2.0, XDiffLow, XDiffHigh);
  pXDiffECalCDet2VsXCDet2 = new TProfile(
      "pXDiffECalCDet2VsXCDet2",
      "Layer 2 mean ECal-CDet x residual vs corrected CDet x;corrected CDet Layer 2 x (m);#LT x_{CDet}-x_{ECal#rightarrowCDet} #GT (m)",
      200, -2.0, 2.0, XDiffLow, XDiffHigh);

  hXCDet1CDet2 = new TH2F("XCDet1CDet2","XCDet1CDet2",200,-0.5,0.5,200,-0.5,0.5);

  // 2D histograms
  h2d_RawLE  = new TH2F("h2d_RawLE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTE  = new TH2F("h2d_RawTE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTot = new TH2F("h2d_RawTot","", NTotBins,TotBinLow,TotBinHigh,nTdc+1,0,nTdc+1);
  h2d_Mult   = new TH2F("h2d_Mult","", 100,0,100,nTdc+1,0,nTdc+1);

  hMultiplicity = new TH1F("hMultiplicity","hMultiplicity",20,0,20);

  for(Int_t tdc=0; tdc<nTdc; tdc++){
    hMultiplicityL[tdc] =  new TH1F(TString::Format("hMultiplicity_Bar%d",tdc),
		      TString::Format("hMultiplicity_Bar%d",tdc),
		      10, 0, 10);
  }// element loop

  hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllRawTe = new TH1F(TString::Format("hRawTe"),
            TString::Format("hRawTe"),
            NTDCBins, TDCBinLow, TDCBinHigh+TotBinHigh);
  hAllRawTot = new TH1F(TString::Format("hRawTot"),
            TString::Format("hRawTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllRawPMT = new TH1F(TString::Format("hRawPMT"),
            TString::Format("hRawPMT"),
            nTdc, 0, nTdc);
  hAllRawBar = new TH1F(TString::Format("hRawBar"),
            TString::Format("hRawBar"),
            168, 0, 168);
  hAllGoodLe = new TH1F(TString::Format("hAllGoodLe"),
            TString::Format("hAllGoodLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllGoodTe = new TH1F(TString::Format("hAllGoodTe"),
            TString::Format("hAllGoodTe"),
            NTDCBins, TDCBinLow, TDCBinHigh+TotBinHigh);
  hAllGoodTot = new TH1F(TString::Format("hAllGoodTot"),
            TString::Format("hAllGoodTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllGoodPMT = new TH1F(TString::Format("hAllGoodPMT"),
            TString::Format("hAllGoodPMT"),
            nTdc, 0, nTdc);
  hAllGoodBar = new TH1F(TString::Format("hAllGoodBar"),
            TString::Format("hAllGoodBar"),
            168, 0, 168);
  h2AllGoodLe = new TH2F(TString::Format("h2AllGoodLe"),
            TString::Format("h2AllGoodLe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTe = new TH2F(TString::Format("h2AllGoodTe"),
            TString::Format("h2AllGoodTe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTot = new TH2F(TString::Format("h2AllGoodTot"),
            TString::Format("h2AllGoodTot"),nTdc,0,nTdc,
            NTotBins, TotBinLow, TotBinHigh);
  hBarRateHV = new TH2F(TString::Format("hBarRateHV"),
	    TString::Format("hBarRateHV"), 250,1,250000,
	    50,600,800);

  h2TDCTOTvsLE = new TH2F(TString::Format("h2TDCTOTvsLE"),
            TString::Format("h2TDCTOTvsLE"),NTotBins,TotBinLow,TotBinHigh,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2CDetX1vsX2 = new TH2F(TString::Format("h2CDetX1vsX2"),
            TString::Format("h2CDetX1vsX2"),1000, -2.0, 2.0,
            1000, -2.0, 2.0);
  h2TOTvsXDiff1 = new TH2F(TString::Format("h2TOTvsXDiff1"),
            TString::Format("h2TOTvsXDiff1"),NTotBins,TotBinLow,TotBinHigh,
            1000, -0.3, 0.3);
  h2TOTvsXDiff2 = new TH2F(TString::Format("h2TOTvsXDiff2"),
            TString::Format("h2TOTvsXDiff2"),NTotBins,TotBinLow,TotBinHigh,
            1000, -0.3, 0.3);
  h2LEvsXDiff1 = new TH2F(TString::Format("h2LEvsXDiff1"),
            TString::Format("h2LEvsXDiff1"),NTDCBins,TDCBinLow,TDCBinHigh,
            1000, -0.3, 0.3);
  h2LEvsXDiff2 = new TH2F(TString::Format("h2LEvsXDiff2"),
            TString::Format("h2LEvsXDiff2"),NTDCBins,TDCBinLow,TDCBinHigh,
            1000, -0.3, 0.3);
  // hRefRawLe = new TH1F(TString::Format("hRefRawLe"),
  //           TString::Format("hRefRawLe"),
  //           RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  // hRefRawTe = new TH1F(TString::Format("hRefRawTe"),
  //           TString::Format("hRefRawTe"),
  //           RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  // hRefRawTot = new TH1F(TString::Format("hRefRawTot"),
  //           TString::Format("hRefRawTot"),
  //           RefNTotBins, RefTotBinLow, RefTotBinHigh);
  // hRefRawPMT = new TH1F(TString::Format("hRefRawPMT"),
  //           TString::Format("hRefRawPMT"),
  //           32, 2688, 2720);
  hRefGoodLe = new TH1F(TString::Format("hRefGoodLe"),
            TString::Format("hRefGoodLe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTe = new TH1F(TString::Format("hRefGoodTe"),
            TString::Format("hRefGoodTe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTot = new TH1F(TString::Format("hRefGoodTot"),
            TString::Format("hRefGoodTot"),
            RefNTotBins, RefTotBinLow, RefTotBinHigh);
  hRefGoodPMT = new TH1F(TString::Format("hRefGoodPMT"),
            TString::Format("hRefGoodPMT"),
            2720, 0, 2720);


  /* The per-channel raw/good histogram fills are disabled below, so avoid
     allocating six large, empty histograms for every mapped TDC channel.
  for(Int_t bar=0; bar<(nTdc); bar++){
    hRawLe[bar] = new TH1F(TString::Format("hRawLe_Bar%d",bar), TString::Format("hRawLe_Bar%d",bar), NTDCBins, TDCBinLow, TDCBinHigh);
    hRawTe[bar] = new TH1F(TString::Format("hRawTe_Bar%d",bar), TString::Format("hRawTe_Bar%d",bar), NTDCBins, TDCBinLow, TDCBinHigh);
    hRawTot[bar] = new TH1F(TString::Format("hRawTot_Bar%d",bar), TString::Format("hRawTot_Bar%d",bar), NTotBins, TotBinLow, TotBinHigh);
    hGoodLe[bar] = new TH1F(TString::Format("hGoodLe_Bar%d",bar), TString::Format("hGoodLe_Bar%d",bar), NTDCBins, TDCBinLow, TDCBinHigh);
    hGoodTe[bar] = new TH1F(TString::Format("hGoodTe_Bar%d",bar), TString::Format("hGoodTe_Bar%d",bar), NTDCBins, TDCBinLow, TDCBinHigh);
    hGoodTot[bar] = new TH1F(TString::Format("hGoodTot_Bar%d",bar), TString::Format("hGoodTot_Bar%d",bar), NTotBins, TotBinLow, TotBinHigh);
  }
  */



  //========================================================= Get data from tree
  if (!T) {
  T = new TChain("T");

  vCDetPaddleRawTot.assign(2688, std::vector<double>{});
  vCDetPaddleCutTot.assign(2688, std::vector<double>{});
  // Per-bar storage for good leading-edge times (bar = GoodElID/16)
  vBarGoodLe.assign(NumPMTs, std::vector<double>{});
  vPaddleGoodLe.assign(2688, std::vector<double>{});
  vPaddleGoodTot.assign(2688, std::vector<double>{});
  vPaddleMatchHCalTime.assign(2688, std::vector<double>{});
  vBarGoodLeECalT.assign(NumPMTs, std::vector<double>{});
  vAllGoodECalT.clear();
  //int onlySegment = -1; // set to >=0 to pick just one

  if (RunNumber1 == 0) {
    std::cout << "RunNumber1 == 0, using hardcoded run list file: "
              << runListFile << "\n";
    AddRunListFilesToChain(T, REPLAYED_DIR.Data(), runListFile, minSeg, maxSeg, groupIndex);
  } else {
    std::cout << "Using single run mode for run " << RunNumber1 << "\n";
    AddRunFilesToChain(T, REPLAYED_DIR.Data(), RunNumber1, minSeg, maxSeg);
  }
}

  if (!T || T->GetEntries() <= 0) {
    std::cerr << "[CDet] ERROR: no input events were found; aborting calibration stage.\n";
    delete T;
    T = nullptr;
    return;
  }

  TTreeReader reader(T);

  //Set TTreeReaders
  /* ----- Earm ----- */

  // ******CDet******
  // ----- CDet arrays -----
  TTreeReaderArray<double> TDCmult(reader, "earm.cdet.tdc_mult");

  TTreeReaderArray<double> RawElID   (reader, "earm.cdet.hits.TDCelemID");
  TTreeReaderArray<double> RawElLE   (reader, "earm.cdet.hits.t");
  TTreeReaderArray<double> RawElTE   (reader, "earm.cdet.hits.t_te");
  TTreeReaderArray<double> RawElTot  (reader, "earm.cdet.hits.t_tot");

  TTreeReaderArray<double> GoodElID  (reader, "earm.cdet.hit.pmtnum");
  TTreeReaderArray<double> GoodElLE  (reader, "earm.cdet.hit.tdc_le");
  TTreeReaderArray<double> GoodElTE  (reader, "earm.cdet.hit.tdc_te");
  TTreeReaderArray<double> GoodElTot (reader, "earm.cdet.hit.tdc_tot");

  TTreeReaderArray<double> GoodX     (reader, "earm.cdet.hit.xhit");
  TTreeReaderArray<double> GoodY     (reader, "earm.cdet.hit.yhit");
  TTreeReaderArray<double> GoodZ     (reader, "earm.cdet.hit.zhit");
  TTreeReaderArray<double> GoodCol   (reader, "earm.cdet.hit.row");
  TTreeReaderArray<double> GoodRow   (reader, "earm.cdet.hit.col");
  TTreeReaderArray<double> GoodLayer (reader, "earm.cdet.hit.layer");

  // ----- cdet scalars -----
  TTreeReaderValue<double> nhits        (reader, "earm.cdet.nhits");
  TTreeReaderValue<double> ngoodhits    (reader, "earm.cdet.ngoodhits");
  TTreeReaderValue<double> ngoodTDChits (reader, "earm.cdet.ngoodTDChits");

  //------ECal-------
  // Cluster arrays
  // TTreeReaderArray<double> ECal_clus_adctime      (reader, "earm.ecal.clus.adctime");
  // // TTreeReaderArray<double> ECal_clus_again        (reader, "earm.ecal.clus.again");
  // TTreeReaderArray<double> ECal_clus_atimeblk     (reader, "earm.ecal.clus.atimeblk");
  // // TTreeReaderArray<double> ECal_clus_col          (reader, "earm.ecal.clus.col");
  TTreeReaderArray<double> ECal_clus_e            (reader, "earm.ecal.clus.e");
  // TTreeReaderArray<double> ECal_clus_eblk         (reader, "earm.ecal.clus.eblk");
  // TTreeReaderArray<double> ECal_clus_id           (reader, "earm.ecal.clus.id");
  // TTreeReaderArray<double> ECal_clus_nblk         (reader, "earm.ecal.clus.nblk");
  // TTreeReaderArray<double> ECal_clus_row          (reader, "earm.ecal.clus.row");
  // TTreeReaderArray<double> ECal_clus_x            (reader, "earm.ecal.clus.x");
  // TTreeReaderArray<double> ECal_clus_y            (reader, "earm.ecal.clus.y");
  TTreeReaderArray<double> ECal_a_time      (reader, "earm.ecal.a_time");
  TTreeReaderArray<double> ECal_a_p         (reader, "earm.ecal.a_p");
  TTreeReaderArray<double> ECal_a_amp_p         (reader, "earm.ecal.a_amp_p");
  TTreeReaderArray<double> ECal_adcxpos         (reader, "earm.ecal.adcxpos");

  // Cluster count (scalar)
  TTreeReaderValue<double> ECal_nclus(reader, "earm.ecal.nclus");

  //event-level ECal branches
  TTreeReaderValue<double> ECalX       (reader, "earm.ecal.x");
  TTreeReaderValue<double> ECalY       (reader, "earm.ecal.y");
  TTreeReaderValue<double> ECalE       (reader, "earm.ecal.e");
  TTreeReaderValue<double> ECalAdcTime (reader, "earm.ecal.adctime");

  // ----- SBS branches -----
  //event-level HCal branches
  TTreeReaderValue<double> HCalAdcTime (reader, "sbs.hcal.adctime");
  TTreeReaderValue<double> HCalE       (reader, "sbs.hcal.e");
  /*  ------- comment out for now ---------
  // Scalars
  TTreeReaderValue<double> heep_dpp(reader, "heep.dpp");
  TTreeReaderValue<double> heep_dt_ADC(reader, "heep.dt_ADC");
  TTreeReaderValue<double> heep_ECalo(reader, "heep.ECalo");
  TTreeReaderValue<double> heep_eprime_eth(reader, "heep.eprime_eth");
  TTreeReaderValue<double> heep_dxECAL(reader, "heep.dxECAL");
  TTreeReaderValue<double> earm_ECal_x(reader, "earm.ECal.x");
  TTreeReaderValue<double> sbs_gemFPP_track_ntrack(reader, "sbs.gemFPP.track.ntrack");
  TTreeReaderValue<double> heep_dyECAL(reader, "heep.dyECAL");

  // Arrays
  TTreeReaderArray<double> sbs_tr_vz(reader, "sbs.tr.vz");
  TTreeReaderArray<double> sbs_gemFPP_track_sclose(reader, "sbs.gemFPP.track.sclose");
  TTreeReaderArray<double> sbs_gemFT_track_nhits(reader, "sbs.gemFT.track.nhits");
  TTreeReaderArray<double> sbs_gemFT_track_ngoodhits(reader, "sbs.gemFT.track.ngoodhits");
  */
  //========================================================= Check no of events


  //Should likely just have Nev = T->GetEntries();, so it just grabs all events ------- NOT WORKING --------
  Int_t Nev = T->GetEntries();
  cout << "N entries in tree is " << Nev << endl;
  Int_t NEventsAnalysis;// = Nev;
  if(nevents==-1) NEventsAnalysis = Nev;
  else NEventsAnalysis = nevents;
  cout << "Running analysis for " << NEventsAnalysis << " events" << endl;


  /*
  //==================================================== Create output root file
  // root file for viewing fits
  TString subfile;
  subfile = TString::Format("gep5_replayed_nogems_%d_50k_events.root",RunNumber1);
  TString outrootfile = ANALYSED_DIR + "/RawTDC_" + subfile;
  TFile *f = new TFile(outrootfile, "RECREATE");
  */



  //================================================================= Event Loop
  // variables outside event loop
  Int_t EventCounter = 0;

  // DEBUG controls
  const bool DBG = false;        // master on/off
  const long DBG_ENTRY = -1;    // set to a specific tree entry, or -1 for all

  cout << "Starting Event Loop" << endl;

  int eff_denominator = 0;
  int eff_numerator_layer1 = 0;
  int eff_numerator_layer2 = 0;
  int eff_numerator = 0;
  int nh0_counter = 0;
  int goodEvCount = 0;
  // event loop start
  Int_t event = 0;
  while(reader.Next()){
    event++;
    event = event - 1;
    EventCounter++;
    // Only stop early if nevents > 0
    if (nevents > 0 && EventCounter > nevents) {
        break;
    }
    gNumEventsInRun++;
    ecalClusterProcessedEventCount++;

    const std::size_t storedECalClusterCount = ECal_clus_e.GetSize();
    const double reportedECalClusterCount = *ECal_nclus;
    const bool ecalClusterCountMatches =
        std::isfinite(reportedECalClusterCount) &&
        reportedECalClusterCount >= 0.0 &&
        std::floor(reportedECalClusterCount) == reportedECalClusterCount &&
        static_cast<std::size_t>(reportedECalClusterCount) == storedECalClusterCount;
    if (!ecalClusterCountMatches) {
      ecalClusterCountMismatchEvents++;
      if (ecalClusterCountMismatchEvents <= 10) {
        std::cerr << "[ECal rate] Cluster-count mismatch at tree entry "
                  << reader.GetCurrentEntry() << ": earm.ecal.nclus="
                  << reportedECalClusterCount << ", clus.e size="
                  << storedECalClusterCount << std::endl;
      }
    }
    for (std::size_t iclus = 0; iclus < storedECalClusterCount; ++iclus) {
      const double energyGeV = ECal_clus_e[iclus];
      if (std::isfinite(energyGeV)) {
        ecalClusterEnergiesGeV.push_back(energyGeV);

        bool inTime = *ECalAdcTime > ECalMinT && *ECalAdcTime < ECalMaxT;
        if (inTime) {
          ecalClusterEnergiesGeVInTime.push_back(energyGeV);
        }
      } else {
        ecalClusterNonFiniteEnergyCount++;
      }
    }

    Int_t nh = *nhits;
    Int_t ngh = *ngoodhits;
    Int_t ngth = *ngoodTDChits;

    if (nh == 0){
      nh0_counter++;
    }
    if (EventCounter % 10000 == 0) {
	cout << EventCounter << "/" << NEventsAnalysis << "/ Nhits = " << (Int_t)nh << endl;
    	for (Int_t nfill=0; nfill<nh; nfill++) {hnhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<ngh; nfill++) {hngoodhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<ngth; nfill++) {hngoodTDChits_ev->Fill(EventCounter);}
    }

    /* Fill ECal cluster vectors */ /////-------------- These need to get moved into the if statement after second pass, rawEventCounter>=1
    // ---- Per-event filling ----
    //v_ECal_clus_adctime.push_back(copyArray(ECal_clus_adctime));
    //v_ECal_clus_again.push_back(copyArray(ECal_clus_again));
    // v_ECal_clus_atimeblk.push_back(copyArray(ECal_clus_atimeblk));

    //v_ECal_clus_col.push_back(copyArray(ECal_clus_col));
    //v_ECal_clus_row.push_back(copyArray(ECal_clus_row));

    //v_ECal_clus_e.push_back(copyArray(ECal_clus_e));
    //v_ECal_clus_eblk.push_back(copyArray(ECal_clus_eblk));

    // v_ECal_clus_id.push_back(copyArray(ECal_clus_id));

    // v_ECal_clus_nblk.push_back(copyArray(ECal_clus_nblk));

    // v_ECal_clus_x.push_back(copyArray(ECal_clus_x));
    // v_ECal_clus_y.push_back(copyArray(ECal_clus_y));

    // Event-level scalars
    //v_ECal_nclus.push_back(*ECal_nclus);



    bool good_elastic = false;
    good_elastic = true; //abs(*heep_dt_ADC)<10 && abs(sbs_tr_vz[0]+0.1)<0.18 && *heep_ECalo/(*heep_eprime_eth) > 0.7 && abs(*heep_dxECAL - 0.01 + 0.025 * (*earm_ECal_x)) < 0.05 && *sbs_gemFPP_track_ntrack > 0 && abs(*heep_dyECAL - 0.01) < 0.06 && sbs_gemFPP_track_sclose[0] < 0.01 && (sbs_gemFT_track_nhits[0] > 4 || sbs_gemFT_track_ngoodhits[0] > 2);
    //currently not using full replays, so elastic cuts dont work, just assume all elastic
    // else if (elastic == 1) good_elastic = true; //incase one does not want to use the elastic cut
    if (good_elastic){

      //fill vectors we wish to make cuts on for selecting elastics
      // Scalars — push_back the dereferenced values
      /* Comment out heep and gems vectors, do not use root files witht them yet
      vheep_dpp.push_back(*heep_dpp);
      vheep_dt_ADC.push_back(*heep_dt_ADC);
      vheep_ECalo.push_back(*heep_ECalo);
      vheep_eprime_eth.push_back(*heep_eprime_eth);
      vheep_dxECAL.push_back(*heep_dxECAL);
      vearm_ECal_x.push_back(*earm_ECal_x);
      vsbs_gemFPP_track_ntrack.push_back(*sbs_gemFPP_track_ntrack);
      vheep_dyECAL.push_back(*heep_dyECAL);

      //2D arrays
      vsbs_tr_vz.push_back(fill2D(sbs_tr_vz));
      vsbs_gemFPP_track_sclose.push_back(fill2D(sbs_gemFPP_track_sclose));
      vsbs_gemFT_track_nhits.push_back(fill2D(sbs_gemFT_track_nhits));
      vsbs_gemFT_track_ngoodhits.push_back(fill2D(sbs_gemFT_track_ngoodhits));
      */
            // Per-entry reference time (ns) for aligning with GOOD-event vectors
      double thisEvent_refRawLe_ns = std::nan("");

// First pass through hits:  purpose is to get reference LE TDC Value for this event

      double event_ref_tdc = 0.0;
      double ref_int = 0;
      double ref_corr = 0;
      for (std::size_t el = 0; el < RawElID.GetSize(); ++el) {
        if ((Int_t)RawElID[el] == 2696) {  // only look at ref PMT
          bool good_ref_le_time = RawElLE[el] > 0.0/TDC_calib_to_ns && RawElLE[el] <= 100.0/TDC_calib_to_ns;
          bool good_ref_tot = RawElTot[el] >= 0.0/TDC_calib_to_ns && RawElTot[el] <= 200.0/TDC_calib_to_ns;
          bool good_ref_event = good_ref_le_time && good_ref_tot;
          if ( good_ref_event ) {

            //if (RawElID[el] > 2687) {
            //	cout << "el = " << el << " Raw ID = " << RawElID[el] << " raw le = " <<
          //	RawElLE[el] << " raw te = " << RawElTE[el] << " raw tot = " <<
          //	RawElTot[el] << " corrected CDet X = " << correctedX << " ECal X = " << ECalX << endl;
            //}
            if ( !check_bad(RawElID[el],suppress_bad) ) {
            //cout << " el = " << el << endl;
            //cout << " tdc = " << RawElLE[el]*TDC_calib_to_ns << endl;
              if ( (Int_t)RawElID[el] == 2696 && (Int_t)RawElLE[el]>0 && (Int_t)RawElTot[el]>0 ) {
                //cout << " Ref  ID = " << (Int_t)RawElID[el] << " el = " << el << "    LE = " << RawElLE[el]*TDC_calib_to_ns
                //		<< "    TE = " << RawElTE[el]*TDC_calib_to_ns << "    ToT = " << RawElTot[el]*TDC_calib_to_ns << endl;

                thisEvent_refRawLe_ns = RawElLE[el] * TDC_calib_to_ns;
                vRefRawLe.push_back(thisEvent_refRawLe_ns);
                vRefRawTe.push_back(RawElTE[el] * TDC_calib_to_ns);
                vRefRawTot.push_back(RawElTot[el] * TDC_calib_to_ns);
                vRefRawPMT.push_back((int)RawElID[el]);
                //hRefRawLe->Fill(RawElLE[el]*TDC_calib_to_ns);
                //hRefRawTe->Fill(RawElTE[el]*TDC_calib_to_ns);
                //hRefRawTot->Fill(RawElTot[el]*TDC_calib_to_ns);
                //hRefRawPMT->Fill(RawElID[el]);

                event_ref_tdc = RawElLE[el]*TDC_calib_to_ns;
                ref_int = std::floor(event_ref_tdc);
                ref_corr = event_ref_tdc - ref_int;
                if (!useReferenceTiming) event_ref_tdc = 0.0;

              }
            }
          }
        }
      }// end ref TDC loop

      // second pass: fill raw CDet TDC histos
      std::vector<double> thisEvent_LE;
      std::vector<double> thisEvent_TE;
      std::vector<double> thisEvent_TOT;
      std::vector<int> thisEvent_ID;

      std::vector<double> thisEvent_TestLE;
      std::vector<double> thisEvent_TestTE;
      std::vector<double> thisEvent_TestTOT;
      std::vector<int> thisEvent_TestID;

      std::vector<double> thisEvent_CDetX;
      std::vector<double> thisEvent_CDetY;
      std::vector<double> thisEvent_CDetZ;
      std::vector<CDetHit> eventHits;
      occupancyEventCount++;

      // Build lookup from PMT id -> index in Good* arrays for this TTree entry
      std::unordered_map<int,int> goodIdx;
      goodIdx.reserve(GoodElID.GetSize());
      for (int ig = 0; ig < (int)GoodElID.GetSize(); ig++) {
        goodIdx[(int)GoodElID[ig]] = ig;
      }

      // Raw hits and good-hit multiplicities have different indexing. Count
      // raw hits by channel so the raw-hit cut uses the collection it filters.
      std::unordered_map<int,int> rawMultiplicity;
      rawMultiplicity.reserve(RawElID.GetSize());
      for (std::size_t iraw = 0; iraw < RawElID.GetSize(); ++iraw) {
        ++rawMultiplicity[(int)RawElID[iraw]];
      }

      int rawEventCounter = 0;
      for (std::size_t el = 0; el < RawElID.GetSize(); ++el){

      const int raw_pmt = (int)RawElID[el];
      auto itGood = goodIdx.find(raw_pmt);
      const bool hasGood = (itGood != goodIdx.end());
      const int ig = hasGood ? itGood->second : -1;
      const double gx = hasGood ? GoodX[ig] : 1.0e9;
      const double gy = hasGood ? GoodY[ig] : 1.0e9;
      const double gz = hasGood ? GoodZ[ig] : 1.0e9;

      // Singles occupancy/rate uses every accepted raw leading edge on a
      // physical CDet channel. Do not apply analysis-level cuts here, and do
      // not deduplicate repeated hits from the same channel in one event.
      const int raw_channel = static_cast<int>(RawElID[el]);
      if (0 <= raw_channel && raw_channel < NumCDetPaddles) {
        rawHitCount[raw_channel]++;
      }

      bool good_raw_le_time = RawElLE[el] >= LeMin/TDC_calib_to_ns && RawElLE[el] <= LeMax/TDC_calib_to_ns;
      bool good_raw_tot = RawElTot[el] >= TotMin/TDC_calib_to_ns && RawElTot[el] <= TotMax/TDC_calib_to_ns;
      bool good_mult = rawMultiplicity[raw_pmt] < TDCmult_cut;
      bool good_CDet_X = hasGood && (fabs(gx) < xcut);
      bool good_ECal_atime = *ECalAdcTime > ECalMinT && *ECalAdcTime < ECalMaxT;
      // Apply the appropriate layer-specific CDet x correction before comparing to ECal.
      // bool good_ECal_diff_y = (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) <= 1.2*CDet_y_half_length &&
      //     (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) >= -1.2*CDet_y_half_length;

      bool good_raw_event = good_raw_le_time && good_raw_tot && good_mult && good_CDet_X && good_ECal_atime;

      //if ((Int_t)RawElID[el] > 1000) cout << "el = " << el << " Hit ID = " << (Int_t)RawElID[el] << "    TDC = " << RawElLE[el]*TDC_calib_to_ns << endl;
      //cout << "Raw ID = " << RawElID[el] << " raw le = " << RawElLE[el] << " raw te = " << RawElTE[el] << " raw tot = " << RawElTot[el] << endl;
      if ( good_raw_event ) {
        rawEventCounter++;
        int idx = RawElID[el];
        if (0 <= idx && idx < 2688) {
          double le_ns = RawElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr(idx);
          double tot_ns = RawElTot[el]*TDC_calib_to_ns;
          double te_ns = RawElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr(idx);
          vCDetPaddleRawTot[idx].push_back(tot_ns);
          eventHits.push_back({idx, le_ns, tot_ns, te_ns});
        } // collecting hit counts and TOT values for pixel occupancy
        if ( !check_bad(RawElID[el],suppress_bad) ) {
        //cout << " el = " << el << endl;
        //cout << " tdc = " << RawElLE[el]*TDC_calib_to_ns << endl;
          if ( (Int_t)RawElID[el] < 2688 ) {
            if (RawElLE[el]*TDC_calib_to_ns - event_ref_tdc < 0) {
              thisEvent_TestLE.push_back(RawElLE[el]*TDC_calib_to_ns);
              thisEvent_TestTE.push_back(RawElTE[el]*TDC_calib_to_ns);
              thisEvent_TestTOT.push_back(RawElTot[el]*TDC_calib_to_ns);
              thisEvent_TestID.push_back((Int_t)RawElID[el]);
            }
            //if ((Int_t)RawElID[el] > nTdc) cout << " CDet ID = " << (Int_t)RawElID[el] << "    TDC = " << RawElLE[el]*TDC_calib_to_ns << endl;

            //fill this events vectors
            thisEvent_LE.push_back(RawElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)RawElID[el]));
            thisEvent_TE.push_back(RawElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)RawElID[el]));
            thisEvent_TOT.push_back(RawElTot[el]*TDC_calib_to_ns); //- event_ref_tdc);
            thisEvent_ID.push_back((Int_t)RawElID[el]);

            //fill all hits vectors
            vAllRawLe.push_back(RawElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)RawElID[el]));
            vAllRawTe.push_back(RawElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)RawElID[el]));
            vAllRawTot.push_back(RawElTot[el]*TDC_calib_to_ns);
            vAllRawPMT.push_back(RawElID[el]);
            vAllRawBar.push_back((Int_t)(RawElID[el]/16));

            thisEvent_CDetX.push_back(gx);
            thisEvent_CDetY.push_back(gy);
            thisEvent_CDetZ.push_back(gz);
            //if (fabs(gx) == 999 && GoodZ[el] != 999){
            if (DBG && (DBG_ENTRY < 0 || reader.GetCurrentEntry() == DBG_ENTRY) && rawEventCounter<20) {
            std::cout << "event = " << rawEventCounter << " " << "cdetX = " << gx << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetY = " << gy << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetZ = " << gz << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetID = " << (Int_t)RawElID[el] << std::endl;
            std::cout << " " <<std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetLE = " << RawElLE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTE = " << RawElTE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTot = " << RawElTot[el]*TDC_calib_to_ns << std::endl;
            std::cout << "-------------------- " <<std::endl;
            }
            /*}
            if (fabs(gx) == 999 && GoodZ[el] != -999){
            std::cout << "event = " << rawEventCounter << " " << "cdetX = " << gx << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetZ = " << gz << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetID = " << (Int_t)RawElID[el] << std::endl;
            std::cout << " " <<std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetLE = " << RawElLE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTE = " << RawElTE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTot = " << RawElTot[el]*TDC_calib_to_ns << std::endl;
            }*/
            /* Comment out histograms for now
            hRawLe[(Int_t)RawElID[el]]->Fill(RawElLE[el]*TDC_calib_to_ns-event_ref_tdc);
            hRawTe[(Int_t)RawElID[el]]->Fill(RawElTE[el]*TDC_calib_to_ns-event_ref_tdc);
            hRawTot[(Int_t)RawElID[el]]->Fill(RawElTot[el]*TDC_calib_to_ns);
            hAllRawLe->Fill(RawElLE[el]*TDC_calib_to_ns-event_ref_tdc);
            hAllRawTe->Fill(RawElTE[el]*TDC_calib_to_ns-event_ref_tdc);
            hAllRawTot->Fill(RawElTot[el]*TDC_calib_to_ns);
            hAllRawPMT->Fill(RawElID[el]);
            hAllRawBar->Fill((Int_t)(RawElID[el]/16));



            h2d_RawLE->Fill(RawElLE[el]*TDC_calib_to_ns-event_ref_tdc, (Int_t)RawElID[el]);
            h2d_RawTE->Fill(RawElTE[el]*TDC_calib_to_ns-event_ref_tdc, (Int_t)RawElID[el]);
            h2d_RawTot->Fill(RawElTot[el]*TDC_calib_to_ns, (Int_t)RawElID[el]);
            */
          }
        }
      }
    }// all raw tdc hit loop
    if (rawEventCounter >= 1){
      vRawLe.push_back(thisEvent_LE);
      vRawTe.push_back(thisEvent_TE);
      vRawTot.push_back(thisEvent_TOT);
      vRawID.push_back(thisEvent_ID);
      vCDetX.push_back(thisEvent_CDetX);
      vCDetY.push_back(thisEvent_CDetY);
      vCDetZ.push_back(thisEvent_CDetZ);
      v_ECalX.push_back(*ECalX);
      v_ECalY.push_back(*ECalY);
      v_ECalE.push_back(*ECalE);
      v_ECalAdcTime.push_back(*ECalAdcTime);
      vEventHits.push_back(eventHits);
      vTestLe.push_back(thisEvent_TestLE);
      vTestTe.push_back(thisEvent_TestTE);
      vTestTot.push_back(thisEvent_TestTOT);
      vTestID.push_back(thisEvent_TestID);
    }
    //check nadjacent pairs for each event
    auto is_unused = [&](int id) -> bool {
      return kUnusedCDetPixels.count(id) != 0;
    };

    auto is_adjacent_with_skip = [&](int a, int b) -> bool {
      if (b == a + 1) {
        // adjacent normally, but only if that neighbor isn't unused
        return !is_unused(b);
      }
      if (b == a + 2) {
        // treat as adjacent if the in-between pixel is unused
        return is_unused(a + 1) && !is_unused(b);
      }
      return false;
    };
    if (rawEventCounter >= 1) {
      std::vector<int> ids;
      const auto& currentEvent = vEventHits.back();
      ids.reserve(currentEvent.size());

      for (const auto& hit : currentEvent) {
        int id = hit.id;
        // if (260 <= id && id <= 270) ids.push_back(id);
        ids.push_back(id);
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

      int nAdjacentHits = 0;
      for (size_t i = 0; i + 1 < ids.size(); i++) {
        if (is_adjacent_with_skip(ids[i], ids[i + 1])) {
          nAdjacentHits++;
        }
      }

      vNumRawAdjacentHits.push_back(nAdjacentHits);
    }

    // Third pass:  Get layer occupancies

    int nhitsc1 = 0;
    int nhitsc2 = 0;
    int ngoodhitsc1 = 0;
    int ngoodhitsc2 = 0;
    int ngoodTDChitsc1 = 0;
    int ngoodTDChitsc2 = 0;
    for (int j=0; j<nTdc; j++) {
      nhits_paddles[j]=0;
      ngoodhits_paddles[j]=0;
      ngoodTDChits_paddles[j]=0;
    }
    npaddles=0;
    ngoodpaddles=0;
    ngoodTDCpaddles=0;

    for (std::size_t el = 0; el < GoodElID.GetSize(); ++el){
      int sbselemid = (Int_t)GoodElID[el];
      int sbsrown = sbselemid%672;
      int sbscoln = sbselemid/672;
      //int sbsrown = (Int_t)GoodRow[el];
      //int sbscoln = (Int_t)GoodCol[el];
      int mylayern = sbscoln/2;
      int mypaddlen = sbscoln*672 + sbsrown;

      if (mylayern == 0) {
        nhitsc1++;
      } else {
        nhitsc2++;
      }
      nhits_paddles[mypaddlen]++;

      bool good_ECal_reconstruction = *ECalY > -1.2 && *ECalY < 1.2 &&
                                      *ECalX > -1.5 && *ECalX < 1.5 &&
                                      *ECalX != 0.00 && *ECalY != 0.00 ;
      bool good_le_time = GoodElLE[el] >= LeMin/TDC_calib_to_ns && GoodElLE[el] <= LeMax/TDC_calib_to_ns;
      bool good_tot = GoodElTot[el] >= TotMin/TDC_calib_to_ns && GoodElTot[el] <= TotMax/TDC_calib_to_ns;
      bool good_hit_mult = TDCmult[el] < TDCmult_cut;
      const double correctedX = CorrectCDetX(GoodX[el], mylayern);
      bool good_CDet_X = correctedX < xcut;
      bool good_ECal_diff_x = (correctedX-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) <= XDiffCut &&
          (correctedX-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) >= -1.0*XDiffCut;
      bool good_ECal_diff_y = (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) <= 1.2*CDet_y_half_length &&
          (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) >= -1.2*CDet_y_half_length;
      bool good_ECal_atime = *ECalAdcTime > ECalMinT && *ECalAdcTime < ECalMaxT;
      bool good_ECal_energy = *ECalE > ECalMinE && *ECalE < ECalMaxE;


      bool good_CDet_event = good_ECal_reconstruction && good_ECal_diff_x && good_ECal_diff_y && good_le_time && good_tot && good_hit_mult && good_CDet_X && good_ECal_atime && good_ECal_energy;


      if (good_CDet_event) {
        if ( !check_bad(GoodElID[el], suppress_bad) ) {
          if ( (Int_t)GoodElID[el]%NumSidesTotal < NumCDetPaddlesPerSide )  {

            //cout << "Hit number " << el << ":    Paddle = " << mypaddlen << " Row = " << sbsrown  << " Col = " << sbscoln  << " hits = " << ngoodTDChits_paddles[mypaddlen] << endl;
            //cout << "el = " << el << " Good ID = " << GoodElID[el] << " Good le = " <<
        //	GoodElLE[el] << " Good te = " << GoodElTE[el] << " Good tot = " <<
        //	GoodElTot[el] << " corrected CDet X = " << correctedX << " ECal X = " << ECalX << endl;
            if (mylayern == 0) {
                ngoodhitsc1++;
            } else {
                ngoodhitsc2++;
            }
            ngoodhits_paddles[mypaddlen]++;
          }
        }
      }
    }
    for (int j=0; j<nTdc; j++) {
      if (nhits_paddles[j] > 0) {
        npaddles++;
        //cout << "Paddle = " << j <<  "  nhits = " << ngoodTDChits_paddles[j] << endl;
      }
    }
    for (int j=0; j<nTdc; j++) {
      if (ngoodhits_paddles[j] > 0) {
        ngoodpaddles++;
        //cout << "Paddle = " << j <<  "  nhits = " << ngoodTDChits_paddles[j] << endl;
      }
    }
        //cout << "event " << event << endl;
        //cout << "Number of good layer 1 hits: " << ngoodTDChitsc1 << endl;
        //cout << "Number of good layer 2 hits: " << ngoodTDChitsc2 << endl;
        //cout << "Layer 1 Hit Cut " << nhitcutlow1 << " " << nhitcuthigh1 << endl;
        //cout << "Layer 2 Hit Cut " << nhitcutlow2 << " " << nhitcuthigh2 << endl;
    vnpaddles.push_back(npaddles);
    vngoodpaddles.push_back(ngoodpaddles);
    //hnpaddles->Fill(npaddles);
    //hngoodpaddles->Fill(ngoodpaddles);

    // Fourth pass:  use layer occupancies to apply additional cuts

    //before loop temp vectors to fill into vGoodLE, etc.
    std::vector<double> thisEvent_GoodLE;
    std::vector<double> thisEvent_GoodTE;
    std::vector<double> thisEvent_GoodTOT;
    std::vector<int> thisEvent_GoodID;

    std::vector<double> thisEvent_GoodX;
    std::vector<double> thisEvent_GoodY;
    std::vector<double> thisEvent_GoodZ;
    std::vector<int> thisEvent_GoodLayer;
    std::vector<int> thisEvent_GoodCol;
    std::vector<CDetHit> goodEventHits;

    std::vector<double> thisEvent_ECal_a_p;
    std::vector<double> thisEvent_ECal_a_amp_p;
    std::vector<double> thisEvent_ECal_a_time;
    std::vector<double> thisEvent_ECal_adcxpos;

    // --- NEW: per-event best (smallest |x-diff|) hit in each layer
    double bestAbsXDiff[2] = {1e99, 1e99};
    double bestXCDet[2]    = {0.0, 0.0};
    double bestXECalProj[2]= {0.0, 0.0};
    bool   foundBest[2]    = {false, false};

    int CDetPassedBoolCount = 0;

    for (std::size_t el = 0; el < GoodElID.GetSize(); ++el){
      const int hitLayer = ((Int_t)GoodElID[el] / NumCDetPaddlesPerSide) / NumSides;
      const double correctedHitX = CorrectCDetX(GoodX[el], hitLayer);
      bool goodhit_ECal_reconstruction = *ECalY > -1.2 && *ECalY < 1.2 &&
                                         *ECalX > -1.5 && *ECalX < 1.5 &&
                                         *ECalX != 0.00 && *ECalY != 0.00;
      bool goodhit_le_time = GoodElLE[el] >= LeMin/TDC_calib_to_ns && GoodElLE[el] <= LeMax/TDC_calib_to_ns;
      bool goodhit_tot = GoodElTot[el] >= TotMin/TDC_calib_to_ns && GoodElTot[el] <= TotMax/TDC_calib_to_ns;
      bool goodhit_hit_mult = TDCmult[el] < TDCmult_cut;
      bool goodhit_CDet_X = correctedHitX < xcut;
      bool goodhit_low = ngoodhitsc1 >= nhitcutlow1  && ngoodhitsc2 >= nhitcutlow2;
      bool goodhit_high  = ngoodhitsc1 <= nhitcuthigh1 && ngoodhitsc2 <= nhitcuthigh2;
      bool goodhit_ECal_atime = *ECalAdcTime > ECalMinT && *ECalAdcTime < ECalMaxT;
      bool goodhit_ECal_energy = *ECalE > ECalMinE && *ECalE < ECalMaxE;
      bool goodhit_ECal_diff_x = (correctedHitX-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) <= XDiffCut &&
          (correctedHitX-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) >= -1.0*XDiffCut;
      bool goodhit_ECal_diff_y = (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) <= 1.2*CDet_y_half_length &&
           (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) >= -1.2*CDet_y_half_length;
      bool goodhit_CDet_event = goodhit_ECal_reconstruction && goodhit_ECal_diff_x && goodhit_ECal_diff_y && goodhit_le_time && goodhit_tot
        && goodhit_hit_mult && goodhit_CDet_X && goodhit_low && goodhit_high && goodhit_ECal_atime && goodhit_ECal_energy;

      if (goodhit_CDet_event) {
        // correctedHitX-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset
        // std::cout << " gx = " << correctedHitX << " & ECalX_Proj = " << (*ECalX)*GoodZ[el]/ECal_dist - XOffset <<std::endl;
        CDetPassedBoolCount++;
        int idx = GoodElID[el];
        if (0 <= idx && idx < 2688) {
          double le_ns = GoodElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr(idx);
          double tot_ns = GoodElTot[el]*TDC_calib_to_ns;
          double te_ns = GoodElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr(idx);
          goodEventHits.push_back({idx, le_ns, tot_ns, te_ns});
        } // collecting good-hit information for pixels
        if ( !check_bad(GoodElID[el], suppress_bad) ) {
          if ( (Int_t)GoodElID[el]%NumSidesTotal < NumCDetPaddlesPerSide )  {
                //cout << "event " << event << endl;
            //cout << "el = " << el << " Good ID = " << GoodElID[el] << " Good le = " <<
          //GoodElLE[el] << " Good te = " << GoodElTE[el] << " Good tot = " <<
          //GoodElTot[el] << " corrected CDet X = " << correctedHitX << " ECal X = " << ECalX << endl;

            //cout << "Filling good timing histos ... " << ngoodTDChitsc1 << " " << endl;

            //std::cout << "Layer = " << (Int_t)GoodLayer[el] << " Side = " << (Int_t)GoodCol[el] << std::endl;

            int sbselem = (Int_t)GoodElID[el];
            int sbsrow = sbselem%672;
            int sbscol = sbselem/672;
            //int sbscol = (Int_t)GoodCol[el];
            //int sbsrow = (Int_t)GoodRow[el];
            int myside = sbscol%2;
            int mylayer = sbscol/2;
            int mypaddle = sbscol*672 + sbsrow;

            if (mylayer == 0) {
              ngoodTDChitsc1++;
              eff_numerator_layer1++;
            }
            else {
              ngoodTDChitsc2++;
              eff_numerator_layer2++;
            }

            ngoodTDChits_paddles[mypaddle]++;

            if ( (layer_choice == 1 && mylayer == 0) || (layer_choice == 2 && mylayer == 1) ||
            (layer_choice == 3 && ngoodhitsc1>=1 && ngoodhitsc2 >= 1) ) {
              eff_numerator++;

              thisEvent_GoodLE.push_back(GoodElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)GoodElID[el]));
              thisEvent_GoodTE.push_back(GoodElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)GoodElID[el]));
              thisEvent_GoodTOT.push_back(GoodElTot[el]*TDC_calib_to_ns);
              thisEvent_GoodID.push_back((Int_t)GoodElID[el]);

              // hGoodLe[(Int_t)GoodElID[el]]->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hGoodTe[(Int_t)GoodElID[el]]->Fill(GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hGoodTot[(Int_t)GoodElID[el]]->Fill(GoodElTot[el]*TDC_calib_to_ns);

              double t_ECal_event = *ECalAdcTime;
              double t_HCal_event = *HCalAdcTime;
              double tLE_bar = GoodElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)GoodElID[el]);
              double tTE_bar = GoodElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetPixelToffsetCorr((Int_t)GoodElID[el]);
              vAllGoodLe.push_back(tLE_bar);
              vAllGoodTe.push_back(tTE_bar);
              vAllGoodTot.push_back(GoodElTot[el]*TDC_calib_to_ns);
              vAllGoodPMT.push_back(GoodElID[el]);
              vAllGoodBar.push_back((Int_t)(GoodElID[el]/16));
              vAllGoodECalT.push_back(t_ECal_event);
              {
                const int bar = (int)(GoodElID[el] / 16);
                if (0 <= bar && bar < NumPMTs) {
                  vBarGoodLe[bar].push_back(tLE_bar);
                  vBarGoodLeECalT[bar].push_back(t_ECal_event);
                }
                const int paddle = GoodElID[el];
                if (0 <= paddle && paddle < 2688){
                  vPaddleGoodLe[paddle].push_back(tLE_bar);
                  vPaddleGoodTot[paddle].push_back(GoodElTot[el]*TDC_calib_to_ns);
                  vPaddleMatchHCalTime[paddle].push_back(t_HCal_event);
                }
              }

              // hAllGoodLe->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hAllGoodTe->Fill(GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hAllGoodTot->Fill(GoodElTot[el]*TDC_calib_to_ns);
              // hAllGoodPMT->Fill(GoodElID[el]);
              // hAllGoodBar->Fill((Int_t)(GoodElID[el]/16));

              // h2AllGoodLe->Fill(GoodElID[el],GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // h2AllGoodTe->Fill(GoodElID[el],GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // h2AllGoodTot->Fill(GoodElID[el],GoodElTot[el]*TDC_calib_to_ns);

              // h2TDCTOTvsLE->Fill(GoodElTot[el]*TDC_calib_to_ns,GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);

              vhitCDetPMT.push_back((Int_t)GoodElID[el]);
              //hHitPMT->Fill((Int_t)GoodElID[el]);
              vRow.push_back((Int_t)GoodRow[el]);
              //hRow->Fill((Int_t)GoodRow[el]);
            }

            if (myside == 0) {
              if (mylayer == 0) {
                vRowLayer1Side1.push_back((Int_t)GoodRow[el]);
                //hRowLayer1Side1->Fill((Int_t)GoodRow[el]);
              }
              else {
                vRowLayer2Side1.push_back((Int_t)GoodRow[el]);
                //hRowLayer2Side1->Fill((Int_t)GoodRow[el]);
              }
            }
            else {
              if(mylayer == 0) {
                vRowLayer1Side2.push_back((Int_t)GoodRow[el]);
                //hRowLayer1Side2->Fill((Int_t)GoodRow[el]);
              }
              else {
                vRowLayer2Side2.push_back((Int_t)GoodRow[el]);
                //hRowLayer2Side2->Fill((Int_t)GoodRow[el]);
              }
            }
            thisEvent_GoodCol.push_back(myside);
            //hCol->Fill(myside);
            thisEvent_GoodLayer.push_back(mylayer);
            //hLayer->Fill(mylayer);

            thisEvent_GoodX.push_back(correctedHitX);
            thisEvent_GoodY.push_back(GoodY[el]);
            thisEvent_GoodZ.push_back(GoodZ[el]);

            // --- NEW: compute projected ECal X at this hit's Z, and update per-layer best if this is smallest |x-diff|
            if (*ECalX != 0.00) {
                const double xECalProj = (*ECalX) * (GoodZ[el]) / ECal_dist;
                const double xdiff     = correctedHitX - xECalProj;
                const double axdiff    = fabs(xdiff);

                if (axdiff < bestAbsXDiff[mylayer]) {
                    bestAbsXDiff[mylayer] = axdiff;
                    bestXCDet[mylayer]    = correctedHitX;
                    bestXECalProj[mylayer]= xECalProj;
                    foundBest[mylayer]    = true;
                }
            }

//------------------------------------------------------- replace hist below
             if (mylayer==0) { //layer 1 "good" histograms & higher level
               //i think we can remove these histograms from here, and put them in their own plot routine, they just need vectors for GoodX positions from CDet and ECal
               h2TOTvsXDiff1->Fill(GoodElTot[el]*TDC_calib_to_ns,correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               h2LEvsXDiff1->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0,correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               hHitXY1->Fill(GoodY[el],correctedHitX);
               hXECalCDet1->Fill(correctedHitX,(*ECalX)*(GoodZ[el])/ECal_dist);
               hYECalCDet1->Fill(GoodY[el],(*ECalY)*(GoodZ[el])/ECal_dist);
               hXDiffECalCDet1->Fill(correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               pXDiffECalCDet1VsXCDet1->Fill(
                   correctedHitX,
                   correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               hXPlusECalCDet1->Fill(correctedHitX+(*ECalX)*(GoodZ[el])/ECal_dist);
               hEECalCDet1->Fill(*ECalE,correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
             }
             else { //layer 2
               h2TOTvsXDiff2->Fill(GoodElTot[el]*TDC_calib_to_ns,correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               h2LEvsXDiff2->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0,correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               hHitXY2->Fill(GoodY[el],correctedHitX);
               hXECalCDet2->Fill(correctedHitX,(*ECalX)*(GoodZ[el])/ECal_dist);
               hYECalCDet2->Fill(GoodY[el],(*ECalY)*(GoodZ[el])/ECal_dist);
               hXDiffECalCDet2->Fill(correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               pXDiffECalCDet2VsXCDet2->Fill(
                   correctedHitX,
                   correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
               hXPlusECalCDet2->Fill(correctedHitX+(*ECalX)*(GoodZ[el])/ECal_dist);
               hEECalCDet2->Fill(*ECalE,correctedHitX-(*ECalX)*(GoodZ[el])/ECal_dist);
             }


          }
          else {
            if (GoodElID[el]==2696){
              vRefGoodLe.push_back(GoodElLE[el]*TDC_calib_to_ns);
              vRefGoodTe.push_back(GoodElTE[el]*TDC_calib_to_ns);
              vRefGoodTot.push_back(GoodElTot[el]*TDC_calib_to_ns);
              vRefGoodPMT.push_back(GoodElID[el]);
            }

            // hRefGoodLe->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
            // hRefGoodTe->Fill(GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
            // hRefGoodTot->Fill(GoodElTot[el]*TDC_calib_to_ns);
            // hRefGoodPMT->Fill(GoodElID[el]);
          }
        }
      }
    }// all good tdc hit loop

    // --- NEW: fill ONLY the single best-matching hit per layer for this event
    if (foundBest[0]) hXECalCDet1_min->Fill(bestXCDet[0], bestXECalProj[0]);
    if (foundBest[1]) hXECalCDet2_min->Fill(bestXCDet[1], bestXECalProj[1]);

    if (CDetPassedBoolCount >= 1){
      const auto nECalBlocks = std::min(
          std::min(ECal_a_p.GetSize(), ECal_a_amp_p.GetSize()),
          std::min(ECal_a_time.GetSize(), ECal_adcxpos.GetSize()));
      for (std::size_t ihit = 0; ihit < nECalBlocks; ++ihit) {
        if (ECal_a_p[ihit] > 0.0) {
          thisEvent_ECal_a_amp_p.push_back(ECal_a_amp_p[ihit]);
          thisEvent_ECal_a_p.push_back(ECal_a_p[ihit]);
          thisEvent_ECal_a_time.push_back(ECal_a_time[ihit]);
          thisEvent_ECal_adcxpos.push_back(ECal_adcxpos[ihit]);
        }
      }
      goodEvCount++;
      v_GoodECalX.push_back(*ECalX);
      v_GoodECalY.push_back(*ECalY);
      v_GoodECalE.push_back(*ECalE);
      v_GoodECalAdcTime.push_back(*ECalAdcTime);
      v_GoodHCalAdcTime.push_back(*HCalAdcTime);
      v_GoodHCalE.push_back(*HCalE);
      v_ECal_a_p.push_back(thisEvent_ECal_a_p);
      v_ECal_a_amp_p.push_back(thisEvent_ECal_a_amp_p);
      v_ECal_a_time.push_back(thisEvent_ECal_a_time);
      v_ECal_adcxpos.push_back(thisEvent_ECal_adcxpos);

      vGoodCol.push_back(thisEvent_GoodCol);
      //hCol->Fill(myside);
      vGoodLayer.push_back(thisEvent_GoodLayer);
      //hLayer->Fill(mylayer);

      vCDetGoodX.push_back(thisEvent_GoodX);
      vCDetGoodY.push_back(thisEvent_GoodY);
      vCDetGoodZ.push_back(thisEvent_GoodZ);
      vTreeEntry.push_back(reader.GetCurrentEntry());
      vGoodRefRawLe.push_back(thisEvent_refRawLe_ns);

      vGoodLe.push_back(thisEvent_GoodLE);
      vGoodTe.push_back(thisEvent_GoodTE);
      vGoodTot.push_back(thisEvent_GoodTOT);
      vGoodID.push_back(thisEvent_GoodID);

      vGoodEventHits.push_back(goodEventHits);
    }
    if (CDetPassedBoolCount >= 1) {
      std::vector<int> ids;
      const auto& currentEvent = vGoodEventHits.back();
      ids.reserve(currentEvent.size());

      for (const auto& hit : currentEvent) {
        int id = hit.id;
        // if (260 <= id && id <= 270) ids.push_back(id);
        ids.push_back(id);
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

      int nAdjacentHits = 0;
      for (size_t i = 0; i + 1 < ids.size(); i++) {
        if (is_adjacent_with_skip(ids[i], ids[i + 1])) {
          nAdjacentHits++;
        }
      }

      vNumGoodAdjacentHits.push_back(nAdjacentHits);
    }//end adjacent check

    if (*ECalX != 0.00 && *ECalY != 0.00) {//double check this later, probably want to fill vectors with ECal hit position
      eff_denominator++;
      hXECal->Fill(*ECalX);
      hYECal->Fill(*ECalY);
      hEECal->Fill(*ECalE);
    };

    // vnhits1.push_back(nhitsc1);
    // vngoodhits1.push_back(ngoodhitsc1);
    // vngoodTDChits1.push_back(ngoodTDChitsc1);
    // vnhits2.push_back(nhitsc2);
    // vngoodhits2.push_back(ngoodhitsc2);
    // vngoodTDChits2.push_back(ngoodTDChitsc2);
    //
    hnhits1->Fill(nhitsc1);
    hngoodhits1->Fill(ngoodhitsc1);
    hngoodTDChits1->Fill(ngoodTDChitsc1);
    hnhits2->Fill(nhitsc2);
    hngoodhits2->Fill(ngoodhitsc2);
    hngoodTDChits2->Fill(ngoodTDChitsc2);

    for (int j=0; j<nTdc; j++) {
      if (ngoodTDChits_paddles[j] > 0) {
        ngoodTDCpaddles++;
        //cout << "Paddle = " << j <<  "  nhits = " << ngoodTDChits_paddles[j] << endl;
      }
    }
    // vngoodTDCpaddles.push_back(ngoodTDCpaddles);
        hngoodTDCpaddles->Fill(ngoodTDCpaddles);

        //cout << "Element loop: " << NdataMult << endl;
    const auto nGoodMultiplicity = std::min(TDCmult.GetSize(), GoodElID.GetSize());
    for (std::size_t tdc = 0; tdc < nGoodMultiplicity; ++tdc){
      if (!check_bad(GoodElID[tdc],suppress_bad)) {
        hMultiplicity->Fill(TDCmult[tdc]);
        hMultiplicityL[(Int_t)GoodElID[tdc]]->Fill(TDCmult[tdc]);
        if( TDCmult[tdc] != 0 ){
          h2d_Mult->Fill(TDCmult[tdc], (Int_t)GoodElID[tdc] );
        }
      }
    }// element loop
    }//good elastic bool
  }//event loop

  std::cout << "events used for occupancy = " << occupancyEventCount << std::endl;
  const double observedTimeSeconds =
      static_cast<double>(occupancyEventCount) * RawSinglesWindowSeconds;

  if (hRawSinglesRateVsID) {
    delete hRawSinglesRateVsID;
    hRawSinglesRateVsID = nullptr;
  }
  hRawSinglesRateVsID = new TH1D(
      "hRawSinglesRateVsID",
      "CDet Raw Singles Rate vs Channel;Channel ID;Raw singles rate [kHz]",
      NumCDetPaddles, 0, NumCDetPaddles);
  hRawSinglesRateVsID->SetDirectory(nullptr);
  hRawSinglesRateVsID->SetStats(0);

  for (int i = 0; i < 2688; i++){
    rawHitOccupancy[i] = occupancyEventCount > 0
        ? static_cast<double>(rawHitCount[i]) / occupancyEventCount
        : 0.0;
    rawSinglesRateHz[i] = observedTimeSeconds > 0.0
        ? static_cast<double>(rawHitCount[i]) / observedTimeSeconds
        : 0.0;
    rawSinglesRateErrorHz[i] = observedTimeSeconds > 0.0
        ? std::sqrt(static_cast<double>(rawHitCount[i])) / observedTimeSeconds
        : 0.0;

    const int bin = i + 1;
    hRawSinglesRateVsID->SetBinContent(bin, rawSinglesRateHz[i] / 1000); //values in units of kHz for histogram
    hRawSinglesRateVsID->SetBinError(bin, rawSinglesRateErrorHz[i] / 1000);
  }

  hRawSinglesRateVsID->GetListOfFunctions()->Add(
      new TParameter<double>("raw_singles_window_s", RawSinglesWindowSeconds));
  hRawSinglesRateVsID->GetListOfFunctions()->Add(
      new TParameter<int>("processed_event_count", occupancyEventCount));

  const double validationRateHz =
      300.0 / (50000.0 * RawSinglesWindowSeconds);
  std::cout << "Raw singles-rate validation: 300 hits / (50000 events x 60 ns) = "
            << validationRateHz << " hits/s"
            << (std::fabs(validationRateHz - 100000.0) < 1.0e-9
                    ? " (expected 100000 hits/s)"
                    : " (VALIDATION FAILED)")
            << std::endl;

  //Second Pass over all events for tot_ave calc

  for (Int_t idx = 0; idx < 2688; idx++){
    Int_t ihits = vCDetPaddleRawTot[idx].size();
    double sumTot = 0;
    for (Int_t i = 0; i < ihits; i++){
      sumTot += vCDetPaddleRawTot[idx][i];
    }
    ave_tot[idx] = ihits > 0 ? sumTot / ihits : 0.0;
  }

  for (Int_t idx = 0; idx < 2688; idx++){
    Int_t ihits = vCDetPaddleRawTot[idx].size();
    for (Int_t i = 0; i < ihits; i++){
      double hit_tot = vCDetPaddleRawTot[idx][i];
      if (hit_tot > ave_tot[idx]){
        vCDetPaddleCutTot[idx].push_back(hit_tot);
        totCutHitCount[idx]++;
      }
    }
    totCutHitOccupancy[idx] = occupancyEventCount > 0
        ? static_cast<double>(totCutHitCount[idx]) / occupancyEventCount
        : 0.0;
  }

  std::cout << "Candidate Events = " << eff_denominator << std::endl;
  std::cout << "Layer 1 Events = " << eff_numerator_layer1 << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator_layer1/eff_denominator  << std::endl;
  std::cout << "Layer 2 Events = " << eff_numerator_layer2 << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator_layer2/eff_denominator <<  std::endl;
  std::cout << "One Good Layer Events = " << eff_numerator << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator/eff_denominator <<  std::endl;

/*
  for (Int_t b=0; b<NumCDetPaddles; b++) {
	//if (hRawLe[b]->GetEntries() > EventCounter/HotChannelRatio) {
    if (hRawLe[b]->GetEntries() > EventCounter/NumCDetPaddles*2*1000) {
      int myhotlayer = b/1344 + 1;
      int myhotside = (b%1344)/672 + 1;
      int myhotmodule = (b%672)/244 + 1;
      int myhotbar = (b%672)%224/16 + 1;
      int myhotpaddle = ((b%672)%224)%16 + 1;
      int mycable = b/16;
      //std::cout << "Hot PMT!! ID = " << b << "  layer = " << myhotlayer <<
      //"   side = " << myhotside << "   module = " << myhotmodule <<
      //"   bar = " << myhotbar << "   paddle_PMT = " << myhotpaddle << " CDet Cable = " << mycable << "   Entries = " << hRawLe[b]->GetEntries() << std::endl;
    }
  }

    /// Get rid of this whole chunk?
    //========================================================== Write histos
  for(Int_t b=0; b<nTdc; b++){
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawLe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawLe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawLe[b]->Write();
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawTe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTe[b]->Write();
    // hRawTe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTot[b]->GetXaxis()->SetTitle("tot (ns)");
    // hRawTot[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTot[b]->Write();
    // hMultiplicityL[b]->GetXaxis()->SetLabelSize(0.06);
    hMultiplicityL[b]->GetXaxis()->SetTitle("tdc ref hit mult");
    hMultiplicityL[b]->Write();
  }
  // hMultiplicity->GetXaxis()->SetLabelSize(0.06);
  hMultiplicity->GetXaxis()->SetTitle("tdc hit multiplicity");
  hMultiplicity->Write();
  // hHitPMT->GetXaxis()->SetLabelSize(0.06);
  hHitPMT->GetXaxis()->SetTitle("Bar ID of Left PMT Hit");
  hHitPMT->SetTitle("");
  hHitPMT->Write();

  // 2D histograms
  h2d_RawLE->GetXaxis()->SetTitle("TDC Leading Edge Time [ns]");
  h2d_RawLE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawLE->SetTitle("");
  h2d_RawLE->Write();
  h2d_RawTE->GetXaxis()->SetTitle("TDC Trailing Edge Time [ns]");
  h2d_RawTE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTE->SetTitle("");
  h2d_RawTE->Write();

  h2d_RawTot->GetXaxis()->SetTitle("TDC Time-over-threshold [ns]");
  h2d_RawTot->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTot->SetTitle("");
  h2d_RawTot->Write();

  h2d_Mult->GetXaxis()->SetTitle("TDC Multiplicity [ns]");
  h2d_Mult->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_Mult->SetTitle("");
  h2d_Mult->Write();
  */
  // Get HV values
  vector<string> HVfilenames = {"l1Left.dat", "l1Right.dat", "l2Left.dat", "l2Right.dat"};

  vector<vector<double>> data = readDataFromFiles(HVfilenames);
  std::vector<double> barRateContents = extractBinContents(hAllRawBar);
  Double_t xval,yval;
  if (data.size() != HVfilenames.size()) {
    std::cerr << "[CDet] WARNING: skipping HV/rate correlation because one or more HV files are incomplete.\n";
  }
  for (std::size_t ii=0; ii<data.size() && ii<HVfilenames.size(); ++ii) {
    for (int jj=0;jj<42;jj++) {
      xval = data[ii][jj];
            yval = barRateContents[ii*42+jj];
      //cout << "Contents:  " << xval << " "  << yval << endl;
      hBarRateHV->Fill(yval,-xval);
    }
  }

  //========================================================== Close output file
  //f->Close();




  //==================================================== ECal-time correction
  // Remove linear correlation between CDet time and ECal ADC time:
  //   <t_CDet> = p0 + p1*t_ECal
  // New calibration files carry a fixed detector-wide delta derived from the
  // calibration run.  Only legacy files lacking delta use the historical
  // sample-dependent recentering below, so existing files remain readable.
  if (gUseECalTimeCorr) {
    const size_t NevCorr = std::min(vGoodLe.size(), v_GoodECalAdcTime.size());
    if (!gECalDeltaLoaded) {
      double sumResLE = 0.0;
      long long nResLE = 0;
      for (size_t ev = 0; ev < NevCorr; ++ev) {
        const double tE = v_GoodECalAdcTime[ev];
        const size_t Nh = std::min(vGoodLe[ev].size(), vGoodID[ev].size());
        for (size_t ih = 0; ih < Nh; ++ih) {
          const int bar = vGoodID[ev][ih] / 16;
          if (bar < 0 || bar >= NumPMTs) continue;
          const double tLE = vGoodLe[ev][ih]; // already includes pixel offset correction
          sumResLE += tLE - (gECalFitP0 + gECalFitP1*tE);
          ++nResLE;
        }
      }
      const double muResLE = nResLE > 0 ? sumResLE/(double)nResLE : 0.0;
      gECalDeltaShift = gTargetMeanLE - muResLE;
      std::cout << "[CDet] WARNING: legacy sample-dependent ECal recentering: muResLE="
                << muResLE << " ns => delta=" << gECalDeltaShift << " ns\n";
    } else {
      std::cout << "[CDet] ECal-time corr enabled with fixed detector calibration: p0="
                << gECalFitP0 << " p1=" << gECalFitP1
                << " delta=" << gECalDeltaShift << " ns\n";
    }

    // --- Apply to ALL stored CDet LE/TE times (good-hit level vectors) ---
    // Convention: corrected time = (t_barcorr - (p0+p1*tE)) + delta
    for (size_t i = 0; i < vAllGoodLe.size() && i < vAllGoodECalT.size() && i < vAllGoodBar.size(); ++i) {
      const int bar = vAllGoodBar[i];
      if (bar < 0 || bar >= NumPMTs) continue;
      const double tE = vAllGoodECalT[i];
      vAllGoodLe[i] = (vAllGoodLe[i] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
      vAllGoodTe[i] = (vAllGoodTe[i] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
    }

    // per-bar LE storage used for bar histograms
    if (vBarGoodLe.size() == (size_t)NumPMTs && vBarGoodLeECalT.size() == (size_t)NumPMTs) {
      for (int bar = 0; bar < NumPMTs; ++bar) {
        const size_t Nh = std::min(vBarGoodLe[bar].size(), vBarGoodLeECalT[bar].size());
        for (size_t j = 0; j < Nh; ++j) {
          const double tE = vBarGoodLeECalT[bar][j];
          vBarGoodLe[bar][j] = (vBarGoodLe[bar][j] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
        }
      }
    }

    // event-level good-hit vectors (used by many timing comparison plots)
    for (size_t ev = 0; ev < NevCorr; ++ev) {
      const double tE = v_GoodECalAdcTime[ev];
      const size_t Nh = std::min(vGoodLe[ev].size(), vGoodID[ev].size());
      for (size_t ih = 0; ih < Nh; ++ih) {
        vGoodLe[ev][ih] = (vGoodLe[ev][ih] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
        vGoodTe[ev][ih] = (vGoodTe[ev][ih] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
      }
    }
  }

  //==================================================== TOT time-walk correction
  // Apply a layer-dependent correction using:
  //   t_corr = t - p1*(1/sqrt(TOT) - 1/sqrt(TOT_ref))
  // This keeps the correction zero at TOT_ref and preserves the overall timing
  // scale while flattening LE vs TOT.
  if (gUseTimeWalkCorr) {
    std::cout << "[CDet] Time-walk corr enabled."
              << "  p1(L1)=" << gTimeWalkP1_L1
              << "  p1(L2)=" << gTimeWalkP1_L2
              << "  TOTref(L1)=" << gTimeWalkTotRef_L1
              << "  TOTref(L2)=" << gTimeWalkTotRef_L2
              << "  calibration fit range=[" << gTimeWalkTotMin << ", "
              << gTimeWalkTotMax << "] ns; applied to every accepted good hit";

    // Flat good-hit vectors used in global timing histograms
    const size_t Nall = std::min(vAllGoodLe.size(),
                         std::min(vAllGoodTe.size(),
                         std::min(vAllGoodTot.size(), vAllGoodPMT.size())));
    for (size_t i = 0; i < Nall; ++i) {
      const int layer = GetLayerFromID(vAllGoodPMT[i]);
      const double twcorr = GetTimeWalkCorrection(layer, vAllGoodTot[i]);
      vAllGoodLe[i] -= twcorr;
      vAllGoodTe[i] -= twcorr;
    }

    // Per-bar LE storage used for bar histograms
    if (vBarGoodLe.size() == (size_t)NumPMTs) {
      for (int bar = 0; bar < NumPMTs; ++bar) {
        const int layer = (bar < 84) ? 1 : 2;
        const size_t Nh = std::min(vBarGoodLe[bar].size(), vGoodTot.size()); // dummy upper bound not used directly
        for (size_t j = 0; j < vBarGoodLe[bar].size(); ++j) {
          // Need the corresponding TOT for this per-bar hit; if not stored, skip here.
          // The authoritative event-level vectors below are corrected and are what the
          // time-walk study/plots should use. The bar histograms remain ECal-corrected.
        }
      }
    }

    // Event-level good-hit vectors used by downstream plots
    const size_t NevTW = std::min(vGoodLe.size(),
                          std::min(vGoodTe.size(),
                          std::min(vGoodTot.size(), vGoodID.size())));
    for (size_t ev = 0; ev < NevTW; ++ev) {
      const size_t Nh = std::min(vGoodLe[ev].size(),
                        std::min(vGoodTe[ev].size(),
                        std::min(vGoodTot[ev].size(), vGoodID[ev].size())));
      for (size_t ih = 0; ih < Nh; ++ih) {
        const int layer = GetLayerFromID(vGoodID[ev][ih]);
        const double twcorr = GetTimeWalkCorrection(layer, vGoodTot[ev][ih]);
        vGoodLe[ev][ih] -= twcorr;
        vGoodTe[ev][ih] -= twcorr;
      }
    }

    // Rebuild the per-bar good-LE vectors from the corrected event-level vectors so the
    // bar-by-bar timing plots stay consistent with the final corrected times.
    vBarGoodLe.assign(NumPMTs, std::vector<double>{});
    for (size_t ev = 0; ev < NevTW; ++ev) {
      const size_t Nh = std::min(vGoodLe[ev].size(), vGoodID[ev].size());
      for (size_t ih = 0; ih < Nh; ++ih) {
        const int bar = vGoodID[ev][ih] / 16;
        if (bar >= 0 && bar < NumPMTs) vBarGoodLe[bar].push_back(vGoodLe[ev][ih]);
      }
    }
  }
  //==================================================== Run-dependent global timing shift
  // Apply a simple additive shift read from CDet_run<run>.dat:
  //   [GlobalTiming]
  //   shift_ns <value>
  if (std::isfinite(gGlobalTimingShift) && std::fabs(gGlobalTimingShift) > 0.0) {
    std::cout << "[CDet] Applying run global timing shift: "
              << gGlobalTimingShift << " ns"
              << (gGlobalTimingLoaded ? " (from file)" : " (default)") << "\n";

    for (std::size_t i = 0; i < vAllGoodLe.size() && i < vAllGoodTe.size(); ++i) {
      vAllGoodLe[i] += gGlobalTimingShift;
      vAllGoodTe[i] += gGlobalTimingShift;
    }

    if (vBarGoodLe.size() == (std::size_t)NumPMTs) {
      for (int bar = 0; bar < NumPMTs; ++bar) {
        for (std::size_t j = 0; j < vBarGoodLe[bar].size(); ++j) {
          vBarGoodLe[bar][j] += gGlobalTimingShift;
        }
      }
    }

    for (std::size_t ev = 0; ev < vGoodLe.size() && ev < vGoodTe.size(); ++ev) {
      const std::size_t Nh = std::min(vGoodLe[ev].size(), vGoodTe[ev].size());
      for (std::size_t ih = 0; ih < Nh; ++ih) {
        vGoodLe[ev][ih] += gGlobalTimingShift;
        vGoodTe[ev][ih] += gGlobalTimingShift;
      }
    }
  }
  // Expose timing-summary values after any requested analysis stage. This lets
  // lightweight run-list drivers compare either raw or calibrated good-hit LE.
  {
    double leSum = 0.0;
    for (double le : vAllGoodLe) leSum += le;
    const double meanLe = vAllGoodLe.empty() ? -999.0 : leSum/vAllGoodLe.size();

    double ecalAdcTimeSum = 0.0;
    for (double time : v_GoodECalAdcTime) ecalAdcTimeSum += time;
    const double meanECalAdcTime = v_GoodECalAdcTime.empty()
        ? -999.0 : ecalAdcTimeSum/v_GoodECalAdcTime.size();

    const int groupMinID[3] = {1472, 1200, 480};
    const int groupMaxID[3] = {1487, 1215, 495};
    std::vector<double> sumEventMeanLe(3, 0.0);
    std::vector<int> eventsWithHit(3, 0);

    const size_t nEvents = std::min(vGoodLe.size(), vGoodID.size());
    for (size_t event = 0; event < nEvents; ++event) {
      std::vector<double> eventLeSum(3, 0.0);
      std::vector<int> eventLeCount(3, 0);
      const size_t nHits = std::min(vGoodLe[event].size(), vGoodID[event].size());

      for (size_t hit = 0; hit < nHits; ++hit) {
        for (int group = 0; group < 3; ++group) {
          if (vGoodID[event][hit] >= groupMinID[group] &&
              vGoodID[event][hit] <= groupMaxID[group]) {
            eventLeSum[group] += vGoodLe[event][hit];
            ++eventLeCount[group];
          }
        }
      }

      for (int group = 0; group < 3; ++group) {
        if (eventLeCount[group] > 0) {
          sumEventMeanLe[group] += eventLeSum[group]/eventLeCount[group];
          ++eventsWithHit[group];
        }
      }
    }

    std::vector<double> groupMeanLe(3, -999.0);
    for (int group = 0; group < 3; ++group) {
      if (eventsWithHit[group] > 0)
        groupMeanLe[group] = sumEventMeanLe[group]/eventsWithHit[group];
    }

    gLastRunMeanGoodLe = meanLe;
    gLastRunMeanECalAdcTime = meanECalAdcTime;
    gLastRunGroupMeanGoodLe = groupMeanLe;
    gLastRunGoodLeCount = vAllGoodLe.size();
    gLastRunGoodECalEventCount = v_GoodECalAdcTime.size();
  }
  std::cout << "nevents with 0 hits = " << nh0_counter << std::endl;
  std::cout << "nGoodEvents = " << goodEvCount << std::endl;
  gLastCalibrationStageSucceeded = gNumEventsInRun > 0;
  if (!gLastCalibrationStageSucceeded)
    std::cerr << "[CDet] ERROR: calibration stage processed no events.\n";
  /* 7/6/2026 BS : debugging stuff, dont need normally
  TH1F*hForMean = new TH1F("hForMean", "hForMean", 61, -1, 60);
  TH1F*hForMeanRef = new TH1F("hForMeanRef", "hForMeanRef", 100, 0, 100);
  for (double x : vAllRawLe) hForMean->Fill(x);
  for (double x : vRefGoodLe) hForMeanRef->Fill(x);
  double le_mean = hForMean->GetMean();
  double le_mean_ref = hForMeanRef->GetMean();

  TCanvas* cMean = new TCanvas("cMean", "cMean", 800, 600);
  hForMean->Draw();

  std::cout << "le_mean = " << le_mean << std::endl;
  std::cout << "le_mean_ref = " << le_mean_ref << std::endl;

  int testBin = LeMax - LeMin;
  int testTotBin = TotMax - TotMin;
  TH1F* hTestLe = new TH1F("hTestLe", "hTestLe", testBin, LeMin, LeMax);
  TH1F* hTestTot = new TH1F("hTestTot", "hTestTot", testTotBin, TotMin, TotMax);
  TH1F* hTestTe = new TH1F("hTestTe", "hTestTe", testBin+10, LeMin, LeMax+10);
  TH2F* hTestLeVsTe = new TH2F("hTestLeVsTe", "hTestLeVsTe", testBin, LeMin, LeMax, testBin+10, LeMin, LeMax+10);
  TH1F* hTestId = new TH1F("hTestId", "hTestId", 2688, 0, 2688);

  for (size_t ev = 0; ev < vTestLe.size(); ++ev) {
    const size_t Nh = vTestLe[ev].size();
    for (size_t ih = 0; ih < Nh; ++ih) {
      hTestLe->Fill(vTestLe[ev][ih]);
      hTestTe->Fill(vTestTe[ev][ih]);
      hTestTot->Fill(vTestTot[ev][ih]);
      hTestLeVsTe->Fill(vTestLe[ev][ih], vTestTe[ev][ih]);
      hTestId->Fill(vTestID[ev][ih]);

    }
  }
  TCanvas* cTest = new TCanvas("cTest", "cTest", 1200, 400);
  cTest->Divide(2,2);
  cTest->cd(1);
  hTestLe->Draw();
  cTest->cd(2);
  hTestTe->Draw();
  cTest->cd(3);
  hTestTot->Draw();
  cTest->cd(4);
  hTestId->Draw();
  */
  //================================================================== End Macro
}// end main

void testFunction(){

  //define histograms

  //loop through hits

  //fill histograms

  //draw histograms/canvas
  //save histograms??

}

void calculateECalClusterRate(double energyThresholdGeV, double energyBinWidthGeV = ECalClusterEnergyBinWidthGeV, double binMinGeV = 0.0, double binMaxGeV = 12.0)
{
  if (ecalClusterProcessedEventCount <= 0) {
    std::cerr << "[ECal rate] No processed-event data are available. Run the main analysis first." << std::endl;
    return;
  }
  if (!std::isfinite(energyThresholdGeV)) {
    std::cerr << "[ECal rate] Energy threshold must be finite." << std::endl;
    return;
  }
  if (!std::isfinite(energyBinWidthGeV) || energyBinWidthGeV <= 0.0) {
    std::cerr << "[ECal rate] Energy-bin width must be positive and finite." << std::endl;
    return;
  }
  if (!std::isfinite(binMinGeV) || !std::isfinite(binMaxGeV) || binMaxGeV <= binMinGeV) {
    std::cerr << "[ECal rate] Energy bounds must be finite, with maximum greater than minimum." << std::endl;
    return;
  }

  if (hECalClusterEnergySpectrum) {
    delete hECalClusterEnergySpectrum;
    hECalClusterEnergySpectrum = nullptr;
  }
  if (hECalClusterEnergySpectrumInTime) {
    delete hECalClusterEnergySpectrumInTime;
    hECalClusterEnergySpectrumInTime = nullptr;
  }
  const int ecalEnergyBinCount = std::max(1, static_cast<int>(std::ceil((binMaxGeV - binMinGeV) / energyBinWidthGeV)));
  hECalClusterEnergySpectrum = new TH1D("hECalClusterEnergySpectrum", "ECal reconstructed-cluster energy;Cluster energy [GeV];Clusters / bin", ecalEnergyBinCount, binMinGeV, binMaxGeV);
  hECalClusterEnergySpectrum->SetDirectory(nullptr);
  hECalClusterEnergySpectrum->SetStats(0);
  hECalClusterEnergySpectrumInTime = new TH1D("hECalClusterEnergySpectrumInTime", "ECal reconstructed-cluster energy (in time);Cluster energy [GeV];Clusters / bin", ecalEnergyBinCount, binMinGeV, binMaxGeV);
  hECalClusterEnergySpectrumInTime->SetDirectory(nullptr);
  hECalClusterEnergySpectrumInTime->SetStats(0);
  for (const double energyGeV : ecalClusterEnergiesGeV) hECalClusterEnergySpectrum->Fill(energyGeV);
  for (const double energyGeVInTime : ecalClusterEnergiesGeVInTime) hECalClusterEnergySpectrumInTime->Fill(energyGeVInTime); //does not currently get rates

  const double ecalObservedTimeSeconds = static_cast<double>(ecalClusterProcessedEventCount) * ECalClusterWindowSeconds;
  const Long64_t passingClusterCount = static_cast<Long64_t>(std::count_if(ecalClusterEnergiesGeV.begin(), ecalClusterEnergiesGeV.end(), [energyThresholdGeV](double energyGeV) { return energyGeV >= energyThresholdGeV; }));
  const double rateHz = static_cast<double>(passingClusterCount) / ecalObservedTimeSeconds;
  const double rateErrorHz = std::sqrt(static_cast<double>(passingClusterCount)) / ecalObservedTimeSeconds;

  hECalClusterEnergySpectrum->GetListOfFunctions()->Add(new TParameter<double>("ecal_cluster_window_s", ECalClusterWindowSeconds));
  hECalClusterEnergySpectrum->GetListOfFunctions()->Add(new TParameter<double>("energy_threshold_gev", energyThresholdGeV));
  hECalClusterEnergySpectrum->GetListOfFunctions()->Add(new TParameter<double>("energy_bin_width_gev", energyBinWidthGeV));
  hECalClusterEnergySpectrum->GetListOfFunctions()->Add(new TParameter<double>("bin_min_gev", binMinGeV));
  hECalClusterEnergySpectrum->GetListOfFunctions()->Add(new TParameter<double>("bin_max_gev", binMaxGeV));
  hECalClusterEnergySpectrum->GetListOfFunctions()->Add(new TParameter<int>("processed_event_count", ecalClusterProcessedEventCount));

  TCanvas *cECalClusterEnergySpectrum = new TCanvas("cECalClusterEnergySpectrum", "ECal reconstructed-cluster energy", 1000, 700);
  hECalClusterEnergySpectrum->Draw("HIST");
  const double lineMaximum = std::max(1.0, 1.05 * hECalClusterEnergySpectrum->GetMaximum());
  TLine *energyThresholdLine = new TLine(energyThresholdGeV, 0.0, energyThresholdGeV, lineMaximum);
  energyThresholdLine->SetLineColor(kRed + 1);
  energyThresholdLine->SetLineWidth(3);
  energyThresholdLine->SetLineStyle(2);
  energyThresholdLine->Draw("SAME");
  TLegend *thresholdLegend = new TLegend(0.58, 0.76, 0.88, 0.88);
  thresholdLegend->AddEntry(energyThresholdLine, TString::Format("E_{thr} = %.3f GeV", energyThresholdGeV), "l");
  thresholdLegend->AddEntry((TObject*)nullptr, TString::Format("Rate = %.3g #pm %.2g Hz", rateHz, rateErrorHz), "");
  thresholdLegend->Draw();
  cECalClusterEnergySpectrum->Update();

  TCanvas *cECalClusterEnergySpectrumInTime = new TCanvas("cECalClusterEnergySpectrumInTime", "ECal reconstructed-cluster energy (in time)", 1000, 700);
  hECalClusterEnergySpectrumInTime->Draw("HIST");
  TLine *energyThresholdLineInTime = new TLine(energyThresholdGeV, 0.0, energyThresholdGeV, lineMaximum);
  energyThresholdLineInTime->SetLineColor(kRed + 1);
  energyThresholdLineInTime->SetLineWidth(3);
  energyThresholdLineInTime->SetLineStyle(2);
  energyThresholdLineInTime->Draw("SAME");
  TLegend *thresholdLegendInTime = new TLegend(0.58, 0.76, 0.88, 0.88);
  thresholdLegendInTime->AddEntry(energyThresholdLineInTime, TString::Format("E_{thr} = %.3f GeV", energyThresholdGeV), "l");
  thresholdLegendInTime->Draw();
  cECalClusterEnergySpectrumInTime->Update();

  std::cout << "[ECal rate] Processed events: " << ecalClusterProcessedEventCount
            << ", cluster entries: " << ecalClusterEnergiesGeV.size() + ecalClusterNonFiniteEnergyCount
            << ", nclus/array mismatches: " << ecalClusterCountMismatchEvents
            << ", nonfinite energies: " << ecalClusterNonFiniteEnergyCount << std::endl;
  std::cout << "[ECal rate] E >= " << energyThresholdGeV << " GeV: " << passingClusterCount
            << " clusters, rate = " << rateHz << " +/- " << rateErrorHz << " Hz" << std::endl;
}

void ResetCalibrationGlobals()
{
    gShiftInvariantPairTimingCuts = false;
    canvas_vector.clear();

    vRefRawLe.clear();
    vGoodRefRawLe.clear();
    vRefRawTe.clear();
    vRefRawTot.clear();
    vRefRawPMT.clear();

    vRefGoodLe.clear();
    vRefGoodTe.clear();
    vRefGoodTot.clear();
    vRefGoodPMT.clear();

    vAllRawLe.clear();
    vAllRawTe.clear();
    vAllRawTot.clear();
    vAllRawPMT.clear();
    vAllRawBar.clear();

    vAllGoodLe.clear();
    vAllGoodTe.clear();
    vAllGoodTot.clear();
    vAllGoodPMT.clear();
    vAllGoodBar.clear();
    vAllGoodECalT.clear();

    vBarGoodLe.clear();
    vBarGoodLeECalT.clear();

    vRawLe.clear();
    vRawTe.clear();
    vRawTot.clear();
    vRawID.clear();

    vGoodLe.clear();
    vGoodTe.clear();
    vGoodTot.clear();
    vGoodID.clear();

    vCDetX.clear();
    vCDetY.clear();
    vCDetZ.clear();

    vCDetGoodX.clear();
    vCDetGoodY.clear();
    vCDetGoodZ.clear();

    vTreeEntry.clear();

    v_GoodECalX.clear();
    v_GoodECalY.clear();
    v_GoodECalE.clear();
    v_GoodECalAdcTime.clear();
    v_GoodHCalAdcTime.clear();
    v_GoodHCalE.clear();

    v_ECal_a_p.clear();
    v_ECal_a_amp_p.clear();
    v_ECal_a_time.clear();
    v_ECal_adcxpos.clear();

    v_ECalX.clear();
    v_ECalY.clear();
    v_ECalE.clear();
    v_ECalAdcTime.clear();

    vTestLe.clear();
    vTestTe.clear();
    vTestTot.clear();
    vTestID.clear();

    vhitCDetPMT.clear();
    vRow.clear();
    vGoodCol.clear();
    vGoodLayer.clear();
    vRowLayer1Side1.clear();
    vRowLayer2Side1.clear();
    vRowLayer1Side2.clear();
    vRowLayer2Side2.clear();
    vnpaddles.clear();
    vngoodpaddles.clear();
    vngoodTDCpaddles.clear();

    vEventHits.clear();
    vGoodEventHits.clear();
    pairs_CDet.clear();
    gCDetDisplayEvents.clear();
    gCDetDisplayIndex = -1;

    vNumRawAdjacentHits.clear();
    vNumGoodAdjacentHits.clear();

    rawHitCount.assign(2688, 0);
    rawHitOccupancy.assign(2688, 0.0);
    rawSinglesRateHz.assign(2688, 0.0);
    rawSinglesRateErrorHz.assign(2688, 0.0);
    totCutHitCount.assign(2688, 0);
    totCutHitOccupancy.assign(2688, 0.0);
    ave_tot.assign(2688, 0.0);

    occupancyEventCount = 0;
    ecalClusterEnergiesGeV.clear();
    ecalClusterProcessedEventCount = 0;
    ecalClusterCountMismatchEvents = 0;
    ecalClusterNonFiniteEnergyCount = 0;

    if (hRawSinglesRateVsID) {
        delete hRawSinglesRateVsID;
        hRawSinglesRateVsID = nullptr;
    }
    if (hECalClusterEnergySpectrum) {
        delete hECalClusterEnergySpectrum;
        hECalClusterEnergySpectrum = nullptr;
    }

    if (T) {
        delete T;
        T = 0;
    }
}

void PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
    const char *configFile, Int_t calibrationStageOverride = -1,
    Int_t eventsOverride = std::numeric_limits<Int_t>::min())
{
  TEnv env;
  const char *caller = "CDet cross-target analysis";
  if (!LoadCDetConfiguration(env, configFile, caller)) return;

  const Int_t runNumber = env.GetValue("analysis.run_number", 5811);
  const Int_t events = eventsOverride != std::numeric_limits<Int_t>::min()
                           ? eventsOverride
                           : env.GetValue("analysis.events", 50000);
  const Int_t calibrationStage = calibrationStageOverride >= 0
                                     ? calibrationStageOverride
                                     : env.GetValue("analysis.calibration_stage", 7);
  const Int_t groupIndex = env.GetValue("analysis.group_index", 0);
  const Int_t minSegment = env.GetValue("analysis.min_segment", -1);
  const Int_t maxSegment = env.GetValue("analysis.max_segment", -1);
  const Double_t leMin = env.GetValue("analysis.le_min", 0.02);
  const Double_t leMax = env.GetValue("analysis.le_max", 100.0);
  const Double_t totMin = env.GetValue("analysis.tot_min", 0.02);
  const Double_t totMax = env.GetValue("analysis.tot_max", 150.0);
  const Double_t ecalTimeMin = env.GetValue("analysis.ecal_time_min", 10.0);
  const Double_t ecalTimeMax = env.GetValue("analysis.ecal_time_max", 35.0);
  const Double_t ecalEnergyMin = env.GetValue("analysis.ecal_energy_min", -1.0e9);
  const Double_t ecalEnergyMax = env.GetValue("analysis.ecal_energy_max", 1.0e9);
  const Int_t layer1HitsMin = env.GetValue("analysis.layer1_hits_min", 1);
  const Int_t layer1HitsMax = env.GetValue("analysis.layer1_hits_max", 100);
  const Int_t layer2HitsMin = env.GetValue("analysis.layer2_hits_min", 0);
  const Int_t layer2HitsMax = env.GetValue("analysis.layer2_hits_max", 100);
  const Double_t xDifferenceMax = env.GetValue("analysis.x_difference_max", 0.05);
  const Double_t xOffset = env.GetValue("analysis.x_offset", 0.0);
  const Double_t yOffset = env.GetValue("analysis.y_offset", 0.1);
  const Int_t layerChoice = env.GetValue("analysis.layer_choice", 3);
  const Bool_t suppressBad = env.GetValue("analysis.suppress_bad", false);
  const Int_t numberOfRuns = env.GetValue("analysis.number_of_runs", 30);
  const Int_t maxStream = env.GetValue("analysis.max_stream", 2);
  const Int_t firstEvent = env.GetValue("analysis.first_event", 1);
  const Bool_t useReferenceTiming = env.GetValue("analysis.use_reference_timing", false);

  if (leMin >= leMax || totMin >= totMax || ecalTimeMin >= ecalTimeMax ||
      ecalEnergyMin >= ecalEnergyMax ||
      layer1HitsMin < 0 || layer1HitsMin > layer1HitsMax ||
      layer2HitsMin < 0 || layer2HitsMin > layer2HitsMax || xDifferenceMax <= 0.0) {
    std::cerr << "[" << caller << "] ERROR: invalid min/max range or occupancy cut in "
              << configFile << ".\n";
    return;
  }

  std::cout << "[" << caller << "] Authoritative configuration: "
            << configFile << "; stage " << calibrationStage
            << "; events " << events << ".\n"
            << "[" << caller << "] Effective cuts: ECal ADC time "
            << ecalTimeMin << " to " << ecalTimeMax << " ns; ECal energy "
            << ecalEnergyMin << " to " << ecalEnergyMax << " GeV; Layer 1 hits "
            << layer1HitsMin << " to " << layer1HitsMax << "; Layer 2 hits "
            << layer2HitsMin << " to " << layer2HitsMax << ".\n";

  const bool previousConfigAuthority = gConfigurationSelectionAuthoritative;
  gConfigurationSelectionAuthoritative = true;
  PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
      runNumber, events, calibrationStage, groupIndex, minSegment, maxSegment,
      leMin, leMax, totMin, totMax, ecalTimeMin, ecalTimeMax,
      layer1HitsMin, layer1HitsMax, layer2HitsMin, layer2HitsMax,
      xDifferenceMax, xOffset, yOffset, layerChoice, suppressBad,
      numberOfRuns, maxStream, firstEvent, useReferenceTiming,
      ecalEnergyMin, ecalEnergyMax);
  gConfigurationSelectionAuthoritative = previousConfigAuthority;
}

void runStats(){
    std::ofstream outfile("RunStats.csv", std::ios::app);

    outfile << gRunNumber << "," << gNumEventsInRun << "\n";

    outfile.close();
}

void plotECalHCalTimeComp(double threshold = 200, int binSizeAmp = 100, int ampMin = 0, int ampMax = 1000,
                         int paddleA = 485, int paddleB = 489, double cdet_shift = 82.0,
                         double leMinA = 10, double leMaxA = 30, double leMinB = 25, double leMaxB = 45,
                         double hcalCutMin = 90, double hcalCutMax = 100,
                         double tMinECal = -40, double tMaxECal = 40,
                         double tMinHCal = 0, double tMaxHCal = 250,
                         double tMinCDet = 0, double tMaxCDet = 60)
{
  TH1::AddDirectory(kFALSE);
  int nbinsECal = tMaxECal - tMinECal;
  int nbinsHCal = tMaxHCal - tMinHCal;
  int nbinsCDet = tMaxCDet - tMinCDet;
  int nbinsAmp = (ampMax - ampMin)/binSizeAmp;
  // Define histograms ---- 3pm July 6 2026 need to add hcal
  TH1D* hECalNoCut = new TH1D("hECalNoCut", "ECal ADC Time (no cut);Time (ns);Counts", nbinsECal, tMinECal, tMaxECal);
  TH1D* hECalWithPaddleA = new TH1D("hECalWithPaddleA", "ECal ADC Time (paddle A cut);Time (ns);Counts", nbinsECal, tMinECal, tMaxECal);
  TH1D* hECalWithPaddleB = new TH1D("hECalWithPaddleB", "ECal ADC Time (paddle B cut);Time (ns);Counts", nbinsECal, tMinECal, tMaxECal);
  TH1D* hHCalNoCut = new TH1D("hHCalNoCut", "HCal ADC Time (no cut);Time (ns);Counts", nbinsHCal, tMinHCal, tMaxHCal);
  TH1D* hHCalWithPaddleA = new TH1D("hHCalWithPaddleA", "HCal ADC Time (with paddle A);Time (ns);Counts", nbinsHCal, tMinHCal, tMaxHCal);
  TH1D* hHCalWithPaddleB = new TH1D("hHCalWithPaddleB", "HCal ADC Time (with paddle B);Time (ns);Counts", nbinsHCal, tMinHCal, tMaxHCal);
  TH1D* hCDetWithPaddleA = new TH1D("hCDetWithPaddleA", "CDet TDC Time (paddle A);Time (ns);Counts", nbinsECal, tMinECal, tMaxECal);
  TH1D* hCDetWithPaddleB = new TH1D("hCDetWithPaddleB", "CDet TDC Time (paddle B);Time (ns);Counts", nbinsCDet, tMinCDet, tMaxCDet);
  TH1D* hECalTimeWithHCalCut = new TH1D("hECalTimeWithHCalCut", "ECal ADC Time (with HCal cut);Time (ns);Counts", nbinsECal, tMinECal, tMaxECal);
  TH2D* hECalVsCDetShifted = new TH2D("hECalVsCDetShifted", "ECal ADC Time vs CDet Shifted Time;CDet Shifted Time (ns);ECal ADC Time (ns)", nbinsECal, tMinECal, tMaxECal, nbinsECal, tMinECal, tMaxECal);
  TH2D* hECalAmpVsTime = new TH2D("hECalAmpVsTime", "ECal Amplitude vs Time;ECal ADC Time (ns);Amplitude (ADC counts)", nbinsECal, tMinECal, tMaxECal, nbinsAmp, ampMin, ampMax);
  TH2D* hECalXVsCDetLE = new TH2D("hECalXVsCDetLE", "ECal X vs CDet Paddle 485 LE;CDet Paddle 485 LE (ns);ECalX (m)",nbinsCDet, tMinCDet, tMaxCDet, 200,-1.5,1.5);
  TH2D* hECalAAmpPVsTime = new TH2D("hECalAAmpVsTime", "ECal A_Amp_P vs ADC Time;ECal ADC Time (ns);ECal A_Amp_P;",nbinsECal, tMinECal, tMaxECal, nbinsAmp, ampMin, ampMax);
  // Extract times
  for (size_t ev = 0; ev < vGoodLe.size(); ++ev) {
    double ecal_t = v_GoodECalAdcTime[ev];
    double hcal_t = v_GoodHCalAdcTime[ev];
    if (v_GoodHCalAdcTime[ev] > 0){
      hHCalNoCut->Fill(hcal_t);
    }

    hECalNoCut->Fill(ecal_t);

    bool hasPaddleA = false;
    bool hasPaddleB = false;
    for (size_t ihE = 0; ihE < v_ECal_a_p[ev].size(); ++ihE) {
      double ecal_amp = v_ECal_a_p[ev][ihE];
      double ecal_time = v_ECal_a_time[ev][ihE];
      hECalAmpVsTime->Fill(ecal_time, ecal_amp);
      hECalAAmpPVsTime->Fill(ecal_time, ecal_amp);
    }
    // Find CDet hits belonging to the selected paddle
    std::vector<double> paddleLE;
    for (size_t ih = 0; ih < vGoodID[ev].size(); ++ih) {
      int id = vGoodID[ev][ih];
      double cdet_t = vGoodLe[ev][ih];
      if (id == paddleA && cdet_t >= leMinA && cdet_t <= leMaxA) {
        paddleLE.push_back(cdet_t);
        hasPaddleA = true;
        hECalVsCDetShifted->Fill(cdet_t+cdet_shift, ecal_t);
        hCDetWithPaddleA->Fill(cdet_t+cdet_shift); // cdet_shift is the offset between CDet and ECal times
      }
      if (id == paddleB && cdet_t >= leMinB && cdet_t <= leMaxB) {
        hasPaddleB = true;
        hCDetWithPaddleB->Fill(cdet_t+cdet_shift); // cdet_shift is the offset between CDet and ECal times
      }
    }

    if (paddleLE.size() != 1) continue;
    const double cdetLE = paddleLE[0];
    // Plot every ECal block position in that event against that one CDet time
    for (size_t ihit = 0; ihit < v_ECal_adcxpos[ev].size(); ihit++) {
      double ecalX = v_ECal_adcxpos[ev][ihit];
      double ecalAP = v_ECal_a_p[ev][ihit];
      if (ecalAP > threshold && hasPaddleA) hECalXVsCDetLE->Fill(cdetLE, ecalX);
    }

    if (hasPaddleA) {
      hECalWithPaddleA->Fill(ecal_t);
      hHCalWithPaddleA->Fill(hcal_t);
    }
    if (hasPaddleB) {
      hECalWithPaddleB->Fill(ecal_t);
      hHCalWithPaddleB->Fill(hcal_t);
    }
  }//all events done

  //looking at ecal time using hcal time cut
  for (size_t ev = 0; ev < v_GoodECalAdcTime.size(); ++ev) {
    double ecal_t = v_GoodECalAdcTime[ev];
    double hcal_t = v_GoodHCalAdcTime[ev];
    if (hcal_t >= hcalCutMin && hcal_t <= hcalCutMax){
      hECalTimeWithHCalCut->Fill(ecal_t);
    }
  }

  TCanvas* cECalTimeComp = new TCanvas("cECalTimeComp", "ECal Time Comparison", 1000, 620);
  cECalTimeComp->SetLogy();
  hECalNoCut->SetLineColor(kBlack);
  hECalWithPaddleA->SetLineColor(kBlue);
  hECalWithPaddleB->SetLineColor(kRed);
  hECalNoCut->Draw();
  hECalWithPaddleA->Draw("SAME");
  hECalWithPaddleB->Draw("SAME");
  cECalTimeComp->BuildLegend(0.7, 0.7, 0.9, 0.9);
  cECalTimeComp->SaveAs("ECalTimeComparison.pdf");

  TCanvas* cHCalTimeComp = new TCanvas("cHCalTimeComp", "HCal Time Comparison", 1000, 620);
  cHCalTimeComp->SetLogy();
  hHCalNoCut->SetLineColor(kBlack);
  hHCalWithPaddleA->SetLineColor(kBlue);
  hHCalWithPaddleB->SetLineColor(kRed);
  hHCalNoCut->Draw();
  hHCalWithPaddleA->Draw("SAME");
  hHCalWithPaddleB->Draw("SAME");
  cHCalTimeComp->BuildLegend(0.7, 0.7, 0.9, 0.9);
  cHCalTimeComp->SaveAs("HCalTimeComparison.pdf");

  TCanvas* cHCalTimeCompNoLog = new TCanvas("cHCalTimeCompNoLog", "HCal Time Comparison (No Log)", 1000, 620);
  hHCalNoCut->SetLineColor(kBlack);
  hHCalWithPaddleA->SetLineColor(kBlue);
  hHCalWithPaddleB->SetLineColor(kRed);
  hHCalNoCut->Draw();
  hHCalWithPaddleA->Draw("SAME");
  hHCalWithPaddleB->Draw("SAME");
  cHCalTimeCompNoLog->BuildLegend(0.7, 0.7, 0.9, 0.9);
  cHCalTimeCompNoLog->SaveAs("HCalTimeComparisonNoLog.pdf");

  TCanvas* cECalTimeCompWithHCalCut = new TCanvas("cECalTimeCompWithHCalCut", "ECal Time Comparison with HCal Cut", 1000, 620);
  cECalTimeCompWithHCalCut->SetLogy();
  hECalNoCut->SetLineColor(kBlack);
  hECalTimeWithHCalCut->SetLineColor(kGreen+2);
  hECalTimeWithHCalCut->Draw();
  hECalNoCut->Draw("SAME");
  cECalTimeCompWithHCalCut->BuildLegend(0.7, 0.7, 0.9, 0.9);
  cECalTimeCompWithHCalCut->SaveAs("ECalTimeComparisonWithHCalCut.pdf");

  TCanvas* cCDetECalOverlay = new TCanvas("cCDetECalOverlay", "CDet and ECal Time Overlay", 1000, 620);
  cCDetECalOverlay->SetLogy();
  hCDetWithPaddleA->SetLineColor(kBlue);
  hECalWithPaddleA->SetLineColor(kRed);
  hECalWithPaddleA->Draw();
  hCDetWithPaddleA->Draw("SAME");
  cCDetECalOverlay->BuildLegend(0.7, 0.7, 0.9, 0.9);
  cCDetECalOverlay->SaveAs("CDetECalTimeOverlay.pdf");

  TCanvas* cECalVsCDetShifted = new TCanvas("cECalVsCDetShifted", "ECal vs CDet Shifted Time", 1000, 620);
  hECalVsCDetShifted->Draw("COLZ");
  cECalVsCDetShifted->SaveAs("ECalVsCDetShifted.pdf");

  TCanvas* cECalAmpVsTime = new TCanvas("cECalAmpVsTime", "ECal Amplitude vs ADC Time", 1000, 620);
  hECalAmpVsTime->Draw("COLZ");
  cECalAmpVsTime->SaveAs("ECalAmpVsTime.pdf");

  TCanvas* cECalXvsCDetLE = new TCanvas("cECalXvsCDetLE", "ECal X vs CDet LE Paddle 485", 1000,620);
  hECalXVsCDetLE->Draw("COLZ");
  cECalXvsCDetLE->SaveAs("ECalXvsCDetLE.pdf");

  // TCanvas* cCDetTimeComp = new TCanvas("cCDetTimeComp", "CDet Time Comparison", 1000, 620);
  // // cCDetTimeComp->SetLogy();
  // hCDetWithPaddleA->SetLineColor(kBlue);
  // hCDetWithPaddleB->SetLineColor(kRed);
  // hCDetWithPaddleA->Draw();
  // hCDetWithPaddleB->Draw("SAME");
  // cCDetTimeComp->BuildLegend(0.7, 0.7, 0.9, 0.9);
  // cCDetTimeComp->SaveAs("CDetTimeComparison.pdf");

}

void plotGoodLeVsTotByLayer(bool overwrite = false,
                            double leMin = 0.0, double leMax = 60.0,
                            double totMin = 0.0, double totMax = 80.0,
                            double leBinWidth = 0.25, double totBinWidth = 0.5,
                            bool drawProfiles = true,
                            bool fitTimeWalk = true,
                            double fitTotMin = 5.0,
                            double fitTotMax = 25.0)
{
  gLastCalibrationFitSucceeded = !fitTimeWalk;
  TH1::AddDirectory(kFALSE);

  int nLeBins  = std::max(1, (int)((leMax  - leMin ) / leBinWidth));
  int nTotBins = std::max(1, (int)((totMax - totMin) / totBinWidth));

  hGoodLeVsTot_L1 = new TH2D("hGoodLeVsTot_L1",
      "Layer 1 Good Hits: LE vs TOT;TOT (ns);LE time (ns)",
      nTotBins, totMin, totMax,
      nLeBins,  leMin,  leMax);

  hGoodLeVsTot_L2 = new TH2D("hGoodLeVsTot_L2",
      "Layer 2 Good Hits: LE vs TOT;TOT (ns);LE time (ns)",
      nTotBins, totMin, totMax,
      nLeBins,  leMin,  leMax);

  const size_t Nev = std::min(vGoodLe.size(),
                     std::min(vGoodTot.size(), vGoodID.size()));

  for (size_t ev = 0; ev < Nev; ++ev) {
    const size_t Nh = std::min(vGoodLe[ev].size(),
                      std::min(vGoodTot[ev].size(), vGoodID[ev].size()));

    for (size_t ih = 0; ih < Nh; ++ih) {
      const double le  = vGoodLe[ev][ih];
      const double tot = vGoodTot[ev][ih];
      const int    id  = vGoodID[ev][ih];

      if (!std::isfinite(le) || !std::isfinite(tot)) continue;
      if (tot <= 0.0) continue;
      if (tot < totMin || tot > totMax) continue;
      if (le  < leMin  || le  > leMax ) continue;

      if (id >= 0 && id <= 1343) {
        hGoodLeVsTot_L1->Fill(tot, le);
      } else if (id >= 1344 && id <= 2687) {
        hGoodLeVsTot_L2->Fill(tot, le);
      }
    }
  }

  TCanvas* cLeVsTot = new TCanvas("cGoodLeVsTotByLayer",
                                  "Good LE vs TOT by Layer", 1200, 600);
  cLeVsTot->Divide(2,1);

  cLeVsTot->cd(1);
  gPad->SetRightMargin(0.12);
  hGoodLeVsTot_L1->Draw("COLZ");

  TProfile* pL1 = nullptr;
  TF1* fTW_L1 = nullptr;
  bool fitL1OK = false;
  if (drawProfiles) {
    pL1 = hGoodLeVsTot_L1->ProfileX("pGoodLeVsTot_L1");
    pL1->SetMarkerStyle(20);
    pL1->SetMarkerSize(0.7);
    pL1->SetLineWidth(2);
    pL1->Draw("SAME");

    if (fitTimeWalk) {
      const double fitMin = std::max(fitTotMin, totMin + 1.0e-6);
      const double fitMax = std::min(fitTotMax, totMax);
      if (fitMax > fitMin) {
        fTW_L1 = new TF1("fTW_L1", "[0] + [1]/sqrt(x)", fitMin, fitMax);
        const double yMean = pL1->GetMean(2);
        fTW_L1->SetParameters(std::isfinite(yMean) ? yMean : 30.0, 5.0);
        const int fitStatus = pL1->Fit(fTW_L1, "QRS");
        fTW_L1->SetLineColor(kRed);
        fTW_L1->SetLineWidth(2);
        fTW_L1->Draw("SAME");

        const double p0 = fTW_L1->GetParameter(0);
        const double p1 = fTW_L1->GetParameter(1);
        const double e0 = fTW_L1->GetParError(0);
        const double e1 = fTW_L1->GetParError(1);
        fitL1OK = fitStatus == 0 && std::isfinite(p0) && std::isfinite(p1);
        if (fitL1OK) {
          gTimeWalkFitP0_L1 = p0;
          gTimeWalkFitP1_L1 = p1;
        }
        std::cout << "\n[plotGoodLeVsTotByLayer] Layer 1 time-walk fit:\n"
                  << "  LE(TOT) = p0 + p1/sqrt(TOT)\n"
                  << "  fit range: " << fitMin << " to " << fitMax << " ns\n"
                  << "  p0 = " << p0 << " +/- " << e0 << " ns\n"
                  << "  p1 = " << p1 << " +/- " << e1 << " ns*sqrt(ns)\n";
        TPaveText* pt1 = new TPaveText(0.14, 0.78, 0.55, 0.92, "NDC");
        pt1->SetFillColor(0);
        pt1->SetTextAlign(12);
        pt1->AddText("LE = p0 + p1/#sqrt{TOT}");
        pt1->AddText(Form("p0 = %.3f #pm %.3f", p0, e0));
        pt1->AddText(Form("p1 = %.3f #pm %.3f", p1, e1));
        pt1->Draw("SAME");
      }
    }
  }

  cLeVsTot->cd(2);
  gPad->SetRightMargin(0.12);
  hGoodLeVsTot_L2->Draw("COLZ");

  TProfile* pL2 = nullptr;
  TF1* fTW_L2 = nullptr;
  bool fitL2OK = false;
  if (drawProfiles) {
    pL2 = hGoodLeVsTot_L2->ProfileX("pGoodLeVsTot_L2");
    pL2->SetMarkerStyle(20);
    pL2->SetMarkerSize(0.7);
    pL2->SetLineWidth(2);
    pL2->Draw("SAME");

    if (fitTimeWalk) {
      const double fitMin = std::max(fitTotMin, totMin + 1.0e-6);
      const double fitMax = std::min(fitTotMax, totMax);
      if (fitMax > fitMin) {
        fTW_L2 = new TF1("fTW_L2", "[0] + [1]/sqrt(x)", fitMin, fitMax);
        const double yMean = pL2->GetMean(2);
        fTW_L2->SetParameters(std::isfinite(yMean) ? yMean : 30.0, 5.0);
        const int fitStatus = pL2->Fit(fTW_L2, "QRS");
        fTW_L2->SetLineColor(kRed);
        fTW_L2->SetLineWidth(2);
        fTW_L2->Draw("SAME");

        const double p0 = fTW_L2->GetParameter(0);
        const double p1 = fTW_L2->GetParameter(1);
        const double e0 = fTW_L2->GetParError(0);
        const double e1 = fTW_L2->GetParError(1);
        fitL2OK = fitStatus == 0 && std::isfinite(p0) && std::isfinite(p1);
        if (fitL2OK) {
          gTimeWalkFitP0_L2 = p0;
          gTimeWalkFitP1_L2 = p1;
        }
        std::cout << "\n[plotGoodLeVsTotByLayer] Layer 2 time-walk fit:\n"
                  << "  LE(TOT) = p0 + p1/sqrt(TOT)\n"
                  << "  fit range: " << fitMin << " to " << fitMax << " ns\n"
                  << "  p0 = " << p0 << " +/- " << e0 << " ns\n"
                  << "  p1 = " << p1 << " +/- " << e1 << " ns*sqrt(ns)\n";
        TPaveText* pt2 = new TPaveText(0.14, 0.78, 0.55, 0.92, "NDC");
        pt2->SetFillColor(0);
        pt2->SetTextAlign(12);
        pt2->AddText("LE = p0 + p1/#sqrt{TOT}");
        pt2->AddText(Form("p0 = %.3f #pm %.3f", p0, e0));
        pt2->AddText(Form("p1 = %.3f #pm %.3f", p1, e1));
        pt2->Draw("SAME");
      }
    }
  }

  if (overwrite && fitTimeWalk) {
    if (!fitL1OK || !fitL2OK) {
      std::cerr << "[CDet] ERROR: time-walk fit failed for one or more layers; constants were not written.\n";
      return;
    }
    bool updated = false;
    if (fitL1OK) {
      gTimeWalkP1_L1 += gTimeWalkFitP1_L1;
      updated = true;
    }
    if (fitL2OK) {
      gTimeWalkP1_L2 += gTimeWalkFitP1_L2;
      updated = true;
    }
    if (updated) {
      // These bounds document the domain used to determine the coefficients.
      // Correction application is controlled by the upstream good-hit ToT cut.
      gTimeWalkTotMin = std::max(fitTotMin, totMin + 1.0e-6);
      gTimeWalkTotMax = std::min(fitTotMax, totMax);
      gTimeWalkParamsLoaded = true;
      gLastCalibrationFitSucceeded = WriteCalibrationConstants(gCalibrationFile);
      std::cout << "[CDet] Updated time-walk parameters in calibration file: "
                << "p1_L1=" << gTimeWalkP1_L1
                << "  p1_L2=" << gTimeWalkP1_L2 << "\n";
    }
  } else if (fitTimeWalk) {
    gLastCalibrationFitSucceeded = fitL1OK && fitL2OK;
  }
}


void plotNumAdjacent(int nbins = 50){
  TH1::AddDirectory(kFALSE);
  TH1D* hNumRawAdjacentHits = new TH1D("hNumRawAdjacentHits", "Number Raw Hits in Adjacent Pixels", nbins, 0, nbins);
  TH1D* hNumGoodAdjacentHits = new TH1D("hNumGoodAdjacentHits", "Number Raw Hits in Adjacent Pixels", nbins, 0, nbins);

  for (const auto& hit : vNumRawAdjacentHits){
    if (hit > 0) hNumRawAdjacentHits->Fill(hit);
  }
  for (const auto& hit : vNumGoodAdjacentHits){
    if (hit > 0) hNumGoodAdjacentHits->Fill(hit);
  }

  TCanvas* cNumAdjacentHits = new TCanvas("cNumAdjacentHits", "Number of Adjacent Hits", 900,700);
  cNumAdjacentHits->Divide(1,2);
  cNumAdjacentHits->cd(1);
  hNumRawAdjacentHits->Draw();
  cNumAdjacentHits->cd(2);
  hNumGoodAdjacentHits->Draw();

}

void plotAveTotPerPixel() {
  const int nPixels = ave_tot.size();

  TH1::AddDirectory(kFALSE);

  TH1D* hAveTot = new TH1D("hAveTot","Average TOT per CDet pixel;Pixel ID;Average TOT (ns)",nPixels, -0.5, nPixels - 0.5);

  for (int p = 0; p < nPixels; ++p) {
    if (ave_tot[p] > 0) {   // optional guard
      hAveTot->SetBinContent(p + 1, ave_tot[p]);
    }
  }

  TCanvas* c = new TCanvas("cAveTot", "Average TOT per Pixel", 1200, 500);
  hAveTot->Draw("HIST");
}

void plotSingleTot(int pixel_base = 0, bool raw = true, double width = 1, double totMin=1, double totMax=80){
  TH1::AddDirectory(kFALSE);
  if (pixel_base % 16 != 0) {
    Error("plotSingleTot", "pixel_base = %d is not a multiple of 16", pixel_base);
    return;
  }

  const int nPlots = 16;
  int TDCBinNum = (int)((totMax-totMin)/width);

  // Decode pixel → {layer, side, submodule, pmt, pixel}
  auto info = getLocation(pixel_base);

  int layer     = info[0] + 1;  // display as 1-based
  int side      = info[1];      // 0=L, 1=R
  int submodule = info[2] + 1;
  int bar       = info[3] + 1;

  TString sideStr = (side == 0) ? "L" : "R";

  TString canvasTitle = Form("Layer %d | %s | Module %d | Bar %d", layer, sideStr.Data(), submodule, bar);
  // Canvas with 4x4 pads
  TString cname = Form("cTot_%d", pixel_base);
  TCanvas* cTot = new TCanvas(cname, canvasTitle, 1200, 1000);
  cTot->Divide(4, 4, 0.001, 0.001);

  // Histogram array
  TH1D* hTot[nPlots];

  for (int i = 0; i < nPlots; i++) {

    int pixel = pixel_base + i;

    TString hname  = Form("hTot_pix%d", pixel);
    if (raw){
      TString htitle = Form("Pixel %d;TOT (ns);Counts", pixel);
      hTot[i] = new TH1D(hname, htitle, TDCBinNum, totMin, totMax);
    }
    if (!raw){
      TString htitle = Form("Pixel %d w/ TOT > %f;TOT (ns);Counts", pixel, ave_tot[pixel]);
      hTot[i] = new TH1D(hname, htitle, TDCBinNum, totMin, totMax);
    }

    // Fill histogram
    if (raw){
      for (const auto& x : vCDetPaddleRawTot[pixel]) {
        hTot[i]->Fill(x);
      }
    }
    if (!raw){
      for (const auto& x : vCDetPaddleCutTot[pixel]) {
        hTot[i]->Fill(x);
      }
    }

    // Draw
    cTot->cd(i + 1);
    hTot[i]->Draw();

    // If unused pixel: draw a black square in the top-right corner of the pad
    if (kUnusedCDetPixels.count(pixel)) {
      // NDC coordinates: (x1,y1,x2,y2) in [0,1] pad coordinates
      TPaveText* flag = new TPaveText(0.82, 0.82, 0.95, 0.95, "NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetBorderSize(1);
      flag->AddText("");       // empty; just a filled box
      flag->Draw("same");
    }
  }

  cTot->Update();
}

void printOccupancy(int pixel, bool applyTotCut = false){
  if (pixel < 0 || pixel >= 2688) {
    std::cerr << "Invalid CDet pixel ID: " << pixel << std::endl;
    return;
  }

  const double occupancy = applyTotCut
      ? totCutHitOccupancy[pixel]
      : rawHitOccupancy[pixel];

  std::cout << (applyTotCut ? "TOT-cut" : "Raw-hit")
            << " occupancy for pixel " << pixel
            << " = " << occupancy
            << " mean hits/event" << std::endl;
}

TCanvas* plotOccupancyVsID(bool applyTotCut = false, bool savePdf = false){
  if (occupancyEventCount <= 0) {
    std::cerr << "Occupancies are unavailable; run the analysis first."
              << std::endl;
    return nullptr;
  }

  struct OccupancyPanel {
    const char* title;
    int firstChannel;
    int lastChannel;
  };

  const OccupancyPanel panels[4] = {
      {"Layer 1 Left",     0,  671},
      {"Layer 1 Right",  672, 1343},
      {"Layer 2 Left",  1344, 2015},
      {"Layer 2 Right", 2016, 2687}
  };

  const auto& occupancy = applyTotCut
      ? totCutHitOccupancy
      : rawHitOccupancy;
  const auto& hitCount = applyTotCut
      ? totCutHitCount
      : rawHitCount;
  const char* modeName = applyTotCut ? "TotCut" : "Raw";
  const char* modeTitle = applyTotCut ? "TOT-cut" : "Raw-hit";

  const TString canvasName = TString::Format("c%sOccupancyVsID", modeName);
  TCanvas* cOccupancy = new TCanvas(
      canvasName.Data(),
      TString::Format("CDet %s Occupancy vs Channel", modeTitle).Data(),
      1500, 900);
  cOccupancy->Divide(2, 2);

  for (int panel = 0; panel < 4; panel++) {
    cOccupancy->cd(panel + 1);

    const TString frameName = TString::Format(
        "h%sOccupancyPanel%d", modeName, panel);
    TH1D* hPanel = new TH1D(
        frameName.Data(),
        TString::Format("CDet %s %s;Channel ID;Mean accepted hits/event",
                        panels[panel].title, modeTitle).Data(),
        panels[panel].lastChannel - panels[panel].firstChannel + 1,
        panels[panel].firstChannel,
        panels[panel].lastChannel + 1);
    hPanel->SetDirectory(nullptr);
    hPanel->SetStats(0);

    TGraphErrors* gPanel = new TGraphErrors();
    gPanel->SetName(TString::Format(
        "g%sOccupancyPanel%d", modeName, panel));
    gPanel->SetMarkerStyle(20);
    gPanel->SetMarkerSize(0.35);

    double panelMaximum = 0.0;
    for (int channel = panels[panel].firstChannel;
         channel <= panels[panel].lastChannel;
         channel++) {
      if (kUnusedCDetPixels.count(channel)) continue;

      const double occupancyError =
          std::sqrt(static_cast<double>(hitCount[channel])) /
          occupancyEventCount;
      const int point = gPanel->GetN();
      gPanel->SetPoint(point, channel, occupancy[channel]);
      gPanel->SetPointError(point, 0.0, occupancyError);
      panelMaximum = std::max(
          panelMaximum, occupancy[channel] + occupancyError);
    }

    hPanel->SetMinimum(0.0);
    hPanel->SetMaximum(panelMaximum > 0.0 ? 1.10 * panelMaximum : 1.0);
    hPanel->Draw();
    gPanel->Draw("PZ SAME");
  }

  cOccupancy->Update();

  if (savePdf) {
    const TString pdfName = TString::Format(
        "%sOccupancyVsID_run%d.pdf", modeName, gRunNumber);
    cOccupancy->SaveAs(pdfName);
  }

  return cOccupancy;
}

TCanvas* plotRawSinglesRateVsID(bool savePdf = false){
  if (!hRawSinglesRateVsID) {
    std::cerr << "Raw singles rates are unavailable; run the analysis first."
              << std::endl;
    return nullptr;
  }

  struct RatePanel {
    const char* name;
    const char* title;
    int firstChannel;
    int lastChannel;
  };

  const RatePanel panels[4] = {
      {"hRawSinglesRateL1L", "Layer 1 Left",     0,  671},
      {"hRawSinglesRateL1R", "Layer 1 Right",  672, 1343},
      {"hRawSinglesRateL2L", "Layer 2 Left",  1344, 2015},
      {"hRawSinglesRateL2R", "Layer 2 Right", 2016, 2687}
  };

  TCanvas* cRawSinglesRate = new TCanvas(
      "cRawSinglesRateVsID", "CDet Raw Singles Rate vs Channel", 1500, 900);
  cRawSinglesRate->Divide(2, 2);

  for (int panel = 0; panel < 4; panel++) {
    cRawSinglesRate->cd(panel + 1);

    TH1D* hPanel = static_cast<TH1D*>(
        hRawSinglesRateVsID->Clone(panels[panel].name));
    hPanel->Reset("ICES");
    hPanel->SetTitle(TString::Format(
        "CDet %s;Channel ID;Raw singles rate [kHz]",
        panels[panel].title));
    hPanel->GetXaxis()->SetRange(
        panels[panel].firstChannel + 1,
        panels[panel].lastChannel + 1);

    TGraphErrors* gPanel = new TGraphErrors();
    gPanel->SetName(TString::Format("gRawSinglesRatePanel%d", panel));
    gPanel->SetMarkerStyle(20);
    gPanel->SetMarkerSize(0.35);

    double panelMaximumKHz = 0.0;
    for (int channel = panels[panel].firstChannel;
         channel <= panels[panel].lastChannel;
         channel++) {
      if (kUnusedCDetPixels.count(channel)) continue;

      const double rateKHz = rawSinglesRateHz[channel] / 1000.0;
      const double rateErrorKHz = rawSinglesRateErrorHz[channel] / 1000.0;
      const int point = gPanel->GetN();
      gPanel->SetPoint(point, channel, rateKHz);
      gPanel->SetPointError(point, 0.0, rateErrorKHz);
      panelMaximumKHz = std::max(panelMaximumKHz, rateKHz + rateErrorKHz);
    }

    hPanel->SetMinimum(0.0);
    hPanel->SetMaximum(panelMaximumKHz > 0.0 ? 1.10 * panelMaximumKHz : 1.0);
    hPanel->Draw();
    gPanel->Draw("PZ SAME");
  }

  cRawSinglesRate->Update();

  if (savePdf) {
    const TString pdfName =
        TString::Format("RawSinglesRateVsID_run%d.pdf", gRunNumber);
    cRawSinglesRate->SaveAs(pdfName);
  }

  return cRawSinglesRate;
}

TCanvas *plotBarRateHV() {

  TCanvas *daa = new TCanvas("All TDC", "All TDC", 50,50,800,800);

  daa->cd();
  gPad->SetLogx();
  hBarRateHV->Draw();

  return daa;
}

void plotPaddleTOT(
    int paddle_base = 480,
    double width = 1.0,
    double TotMin = 0.0,
    double TotMax = 80.0,
    double fitLow = 0.0,
    double fitHigh = 80.0,
    double binLeLow = 0.0,
    double binLeHigh = 60.0,
    double binTotLow = 0.0,
    double binTotHigh = 80.0,
    double hcalMinTime = 0.0,
    double hcalMaxTime = 250.0,
    const std::vector<double>& externalTotMeans = {},
    const std::vector<double>& externalTotSigmas = {}
)
{
  // Retained for existing callers; this routine uses binTotLow/High and
  // the fitted or external ToT cuts below, rather than these legacy limits.
  (void)TotMin;
  (void)TotMax;
  TH1::AddDirectory(kFALSE);
  int nTotBins = (int)((binTotHigh-binTotLow)/width);
  int nLeBins = (int)((binLeHigh-binLeLow)/width);

  const int nPaddles = 16;

  std::vector<TH1F*> hPaddleTot(nPaddles, nullptr);
  std::vector<TH1F*> hPaddleLe(nPaddles, nullptr);
  std::vector<TF1*>  fPaddleGaussFit(nPaddles, nullptr);

  //identify is external means and sigmas are provided for all paddles
  const bool useExternalTotCuts = externalTotMeans.size() == nPaddles && externalTotSigmas.size() == nPaddles;

  int bar = int(paddle_base/16);
  TCanvas* cPaddlesTOT = new TCanvas("cPaddlesTOT", TString::Format("Bar %d Pixels", bar), 1000, 620);
  cPaddlesTOT->Divide(4, 4, 0.001, 0.001);

  TCanvas* cPaddlesLE = new TCanvas("cPaddlesLE", TString::Format("Bar %d Pixels LE w/ TOT Cut", bar), 1000, 620);
  cPaddlesLE->Divide(4, 4, 0.001, 0.001);

  for (int paddle = 0; paddle < nPaddles ; ++paddle){
    int global_paddle = paddle + paddle_base;
    hPaddleTot[paddle] = new TH1F(TString::Format("hBarTot_Paddle%d", global_paddle),
                               TString::Format("TOT (Paddle %d)", global_paddle),
                               nTotBins, binTotLow, binTotHigh);

    hPaddleLe[paddle] = new TH1F(TString::Format("hBarLe_Paddle%d", global_paddle),
                               TString::Format("LE (Paddle %d)", global_paddle),
                               nLeBins, binLeLow, binLeHigh);
    for (size_t ihit = 0; ihit < vPaddleGoodTot[global_paddle].size(); ++ihit){
      double x = vPaddleGoodTot[global_paddle][ihit];
      double t_HCal = vPaddleMatchHCalTime[global_paddle][ihit];
      if (t_HCal >= hcalMinTime && t_HCal <= hcalMaxTime) hPaddleTot[paddle]->Fill(x);
    }

    double paddleAmpTot = hPaddleTot[paddle]->GetMaximum();
    double paddleMeanTot = hPaddleTot[paddle]->GetMean();
    double paddleRmsTot = hPaddleTot[paddle]->GetRMS();

    //Define & fit Tot spectra
    fPaddleGaussFit[paddle] = new TF1(TString::Format("fPaddleGaussFit_%d", global_paddle), "gaus", fitLow, fitHigh);
    if (hPaddleTot[paddle]->GetEntries() > 20){
      fPaddleGaussFit[paddle]->SetParameters(paddleAmpTot, paddleMeanTot, paddleRmsTot);
      hPaddleTot[paddle]->Fit(fPaddleGaussFit[paddle], "RQ0");
    }
    double fitMean = fPaddleGaussFit[paddle]->GetParameter(1);
    double fitSigma = fPaddleGaussFit[paddle]->GetParameter(2);

    std::cout << "Paddle " << global_paddle << " TOT mean = " << fitMean << " sigma = " << fitSigma << std::endl;
    double cutMean = fitMean;
    double cutSigma = fitSigma;
    bool validTotCut = std::isfinite(cutMean) && std::isfinite(cutSigma) && cutSigma >= 0;
    //examine paddles le spectra using tot cuts - Option of using fit values or external values
    if (useExternalTotCuts) {
      cutMean = externalTotMeans[paddle];
      cutSigma = externalTotSigmas[paddle];
      validTotCut = std::isfinite(cutMean) && std::isfinite(cutSigma) && cutSigma >= 0;
    }

    double totCutLow = 0;
    double totCutHigh = 0;

    if (validTotCut) {
      totCutLow = cutMean - 2.0*cutSigma;
      totCutHigh = cutMean + 2.0*cutSigma;

      for (size_t ihit = 0; ihit < vPaddleGoodTot[global_paddle].size();ihit++){
        double tot = vPaddleGoodTot[global_paddle][ihit];
        double le = vPaddleGoodLe[global_paddle][ihit];
        double t_HCal = vPaddleMatchHCalTime[global_paddle][ihit];
        if (tot > totCutLow && tot < totCutHigh && t_HCal >= hcalMinTime && t_HCal <= hcalMaxTime) hPaddleLe[paddle]->Fill(le);
      }
    }
    // Draw
    cPaddlesTOT->cd(paddle + 1);
    hPaddleTot[paddle]->Draw();
    fPaddleGaussFit[paddle]->Draw("SAME");
    if (validTotCut) {
      const double lineHeight = hPaddleTot[paddle]->GetMaximum();

      TLine* lowLine = new TLine(totCutLow, 0.0, totCutLow, lineHeight);
      TLine* highLine = new TLine(totCutHigh, 0.0, totCutHigh, lineHeight);

      lowLine->SetLineStyle(2);
      highLine->SetLineStyle(2);

      lowLine->Draw("same");
      highLine->Draw("same");
    }
    // If unused pixel: draw a black square in the top-right corner of the pad
    if (kUnusedCDetPixels.count(global_paddle)) {
      // NDC coordinates: (x1,y1,x2,y2) in [0,1] pad coordinates
      TPaveText* flag = new TPaveText(0.82, 0.82, 0.95, 0.95, "NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetBorderSize(1);
      flag->AddText("UNUSED PIXEL");       // empty; just a filled box
      flag->Draw("same");
    }

    //Draw LE spectra with TOT cuts
    cPaddlesLE->cd(paddle + 1);
    hPaddleLe[paddle]->Draw();

    if (kUnusedCDetPixels.count(global_paddle)) {
      // NDC coordinates: (x1,y1,x2,y2) in [0,1] pad coordinates
      TPaveText* flag = new TPaveText(0.82, 0.82, 0.95, 0.95, "NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetTextColor(kWhite);
      flag->SetBorderSize(1);
      flag->AddText("UNUSED PIXEL");       // empty; just a filled box
      flag->Draw("same");
    }
  }

  cPaddlesTOT->Modified();
  cPaddlesTOT->Update();

  cPaddlesLE->Modified();
  cPaddlesLE->Update();
}

static TString CDetPixelCutDirectory(int logicalPixelID)
{
  return TString::Format("pixel_%04d", logicalPixelID);
}

static TCutG *LoadCDetPixelLeTotCut(TFile& input, int logicalPixelID)
{
  if (input.IsZombie()) return nullptr;
  const TString objectPath = CDetPixelCutDirectory(logicalPixelID) + "/cut_le_vs_tot";
  TCutG *stored = dynamic_cast<TCutG*>(input.Get(objectPath));
  if (!stored) return nullptr;
  TCutG *copy = static_cast<TCutG*>(stored->Clone(
      TString::Format("cut_le_vs_tot_pixel_%04d_loaded", logicalPixelID)));
  copy->SetVarX("TOT_ns");
  copy->SetVarY("LE_ns");
  return copy;
}

static TH2D *BuildCDetPixelLeVsTotHistogram(int logicalPixelID,
                                             double totMin, double totMax,
                                             double leMin, double leMax,
                                             double binWidth,
                                             const TString& name)
{
  if (logicalPixelID < 0 || logicalPixelID >= NumCDetPaddles ||
      totMin >= totMax || leMin >= leMax || binWidth <= 0.0) return nullptr;
  const int nTotBins = std::max(1, (int)((totMax - totMin)/binWidth));
  const int nLeBins = std::max(1, (int)((leMax - leMin)/binWidth));
  TH2D *histogram = new TH2D(
      name,
      TString::Format("Good LE vs TOT (Pixel %d);TOT [ns];LE [ns]", logicalPixelID),
      nTotBins, totMin, totMax, nLeBins, leMin, leMax);
  histogram->SetDirectory(nullptr);
  if ((int)vPaddleGoodLe.size() <= logicalPixelID ||
      (int)vPaddleGoodTot.size() <= logicalPixelID) return histogram;
  const size_t nHits = std::min(vPaddleGoodLe[logicalPixelID].size(),
                                vPaddleGoodTot[logicalPixelID].size());
  for (size_t hit = 0; hit < nHits; ++hit)
    histogram->Fill(vPaddleGoodTot[logicalPixelID][hit],
                    vPaddleGoodLe[logicalPixelID][hit]);
  return histogram;
}

void plotCDetPixelLeVsTot(int logicalPixelID, TString cutFile,
                          double totMin, double totMax,
                          double leMin, double leMax, double binWidth)
{
  static unsigned long invocation = 0;
  TH2D *histogram = BuildCDetPixelLeVsTotHistogram(
      logicalPixelID, totMin, totMax, leMin, leMax, binWidth,
      TString::Format("hCDetPixelLeVsTot_%04d_%lu", logicalPixelID, ++invocation));
  if (!histogram) {
    std::cerr << "[CDet pixel LE/TOT] ERROR: invalid pixel ID or histogram range.\n";
    return;
  }
  TCanvas *canvas = new TCanvas(
      TString::Format("cCDetPixelLeVsTot_%04d_%lu", logicalPixelID, invocation),
      TString::Format("Pixel %d LE versus TOT", logicalPixelID), 900, 750);
  canvas->SetLogz();
  histogram->Draw("COLZ");

  TCutG *cut = nullptr;
  if (!cutFile.IsNull() && !gSystem->AccessPathName(cutFile)) {
    TFile input(cutFile, "READ");
    cut = LoadCDetPixelLeTotCut(input, logicalPixelID);
    input.Close();
  }
  if (cut) {
    cut->SetLineColor(kRed + 1);
    cut->SetLineWidth(3);
    cut->SetFillStyle(0);
    cut->Draw("L SAME");
    std::cout << "[CDet pixel LE/TOT] Displaying saved cut for pixel "
              << logicalPixelID << " from " << cutFile << ".\n";
  } else {
    std::cout << "[CDet pixel LE/TOT] No saved cut for pixel "
              << logicalPixelID << " in " << cutFile << ".\n";
  }
  canvas->Modified();
  canvas->Update();
}

void editCDetPixelLeTotCut(int logicalPixelID, TString cutFile,
                           double totMin, double totMax,
                           double leMin, double leMax, double binWidth)
{
  if (gROOT->IsBatch()) {
    std::cerr << "[CDet pixel-cut editor] ERROR: interactive ROOT mode is required.\n";
    return;
  }
  static unsigned long invocation = 0;
  TH2D *histogram = BuildCDetPixelLeVsTotHistogram(
      logicalPixelID, totMin, totMax, leMin, leMax, binWidth,
      TString::Format("hCDetPixelLeVsTotEdit_%04d_%lu", logicalPixelID, ++invocation));
  if (!histogram) {
    std::cerr << "[CDet pixel-cut editor] ERROR: invalid pixel ID or histogram range.\n";
    return;
  }
  if (histogram->GetEntries() <= 0) {
    std::cerr << "[CDet pixel-cut editor] ERROR: no good hits are available for pixel "
              << logicalPixelID << ". Run the main analysis first.\n";
    return;
  }

  TCanvas *canvas = new TCanvas(
      TString::Format("cCDetPixelLeTotCutEdit_%04d_%lu", logicalPixelID, invocation),
      TString::Format("Draw physical-population cut for pixel %d", logicalPixelID),
      950, 800);
  canvas->SetLogz();
  histogram->Draw("COLZ");

  TCutG *oldCut = nullptr;
  if (!cutFile.IsNull() && !gSystem->AccessPathName(cutFile)) {
    TFile oldFile(cutFile, "READ");
    oldCut = LoadCDetPixelLeTotCut(oldFile, logicalPixelID);
    oldFile.Close();
  }
  if (oldCut) {
    oldCut->SetLineColor(kGreen + 2);
    oldCut->SetLineStyle(2);
    oldCut->SetLineWidth(3);
    oldCut->SetFillStyle(0);
    oldCut->Draw("L SAME");
  }
  canvas->Modified();
  canvas->Update();

  std::cout << "\n[CDet pixel-cut editor] Pixel " << logicalPixelID << "\n"
            << "  Draw a polygon around the physical LE-versus-TOT population.\n"
            << "  Left-click to add vertices and double-click to close the polygon.\n"
            << "  The dashed green line, when present, is the previously saved cut.\n"
            << "  Press Escape to cancel without changing the cut file.\n";
  TObject *primitive = gPad->WaitPrimitive("CUTG", "CutG");
  TCutG *drawnCut = dynamic_cast<TCutG*>(primitive);
  if (drawnCut) drawnCut->SetName(
      TString::Format("drawn_cut_pixel_%04d_%lu", logicalPixelID, invocation));
  double twiceArea = 0.0;
  if (drawnCut && drawnCut->GetN() >= 4) {
    for (int point = 0; point < drawnCut->GetN() - 1; ++point)
      twiceArea += drawnCut->GetPointX(point)*drawnCut->GetPointY(point + 1) -
                   drawnCut->GetPointX(point + 1)*drawnCut->GetPointY(point);
  }
  if (!drawnCut || drawnCut->GetN() < 4 || std::fabs(twiceArea) < 1.0e-6) {
    std::cout << "[CDet pixel-cut editor] No valid polygon was drawn. A cut needs at "
              << "least three distinct vertices and nonzero area; nothing saved.\n";
    return;
  }

  size_t selectedHits = 0;
  if ((int)vPaddleGoodLe.size() > logicalPixelID &&
      (int)vPaddleGoodTot.size() > logicalPixelID) {
    const size_t nHits = std::min(vPaddleGoodLe[logicalPixelID].size(),
                                  vPaddleGoodTot[logicalPixelID].size());
    for (size_t hit = 0; hit < nHits; ++hit)
      if (drawnCut->IsInside(vPaddleGoodTot[logicalPixelID][hit],
                             vPaddleGoodLe[logicalPixelID][hit])) ++selectedHits;
  }
  if (selectedHits == 0) {
    std::cout << "[CDet pixel-cut editor] The polygon contains zero pixel hits; "
              << "nothing saved.\n";
    return;
  }

  TCutG *savedCut = static_cast<TCutG*>(drawnCut->Clone("cut_le_vs_tot"));
  savedCut->SetTitle(TString::Format(
      "Pixel %d physical population;TOT [ns];LE [ns]", logicalPixelID));
  savedCut->SetVarX("TOT_ns");
  savedCut->SetVarY("LE_ns");
  savedCut->SetLineColor(kRed + 1);
  savedCut->SetLineStyle(1);
  savedCut->SetLineWidth(3);
  savedCut->SetFillStyle(0);

  TFile output(cutFile, "UPDATE");
  if (output.IsZombie()) {
    std::cerr << "[CDet pixel-cut editor] ERROR: could not update " << cutFile << ".\n";
    return;
  }
  const TString directoryName = CDetPixelCutDirectory(logicalPixelID);
  TDirectory *directory = output.GetDirectory(directoryName);
  if (!directory) directory = output.mkdir(directoryName);
  if (!directory) {
    std::cerr << "[CDet pixel-cut editor] ERROR: could not create directory "
              << directoryName << " in " << cutFile << ".\n";
    output.Close();
    return;
  }
  directory->cd();
  savedCut->Write("cut_le_vs_tot", TObject::kOverwrite);
  histogram->Write("h_le_vs_tot_reference", TObject::kOverwrite);
  TParameter<int>("logical_pixel_id", logicalPixelID).Write("logical_pixel_id", TObject::kOverwrite);
  TParameter<int>("source_run", gRunNumber).Write("source_run", TObject::kOverwrite);
  TParameter<int>("calibration_stage", gCalibrationStage).Write("calibration_stage", TObject::kOverwrite);
  output.Close();

  savedCut->Draw("L SAME");
  canvas->Modified();
  canvas->Update();
  std::cout << "[CDet pixel-cut editor] Saved pixel " << logicalPixelID
            << " cut to " << cutFile << " (" << selectedHits
            << " of " << histogram->GetEntries() << " displayed hits selected).\n";
}

static double ParseCDetFitResultNumber(const std::string& token)
{
  char *end = nullptr;
  const double value = std::strtod(token.c_str(), &end);
  return end && end != token.c_str() ? value : NAN;
}

bool writeCDetCalibrationBaselineFromFitResults(
    TString resultsFile = "CDet_pixel_timing_fit_results.after_2_manual_cuts.dat",
    TString templateCalibration = "CDet_calibration_dt.after_2_manual_cuts.dat",
    TString baselineOutput = "CDet_calibration_dt.before_manual_le_tot.dat")
{
  if (!LoadCalibrationConstants(templateCalibration.Data())) {
    std::cerr << "[CDet baseline recovery] ERROR: could not load calibration template "
              << templateCalibration << ".\n";
    return false;
  }
  std::ifstream input(resultsFile.Data());
  if (!input) {
    std::cerr << "[CDet baseline recovery] ERROR: could not open " << resultsFile << ".\n";
    return false;
  }
  std::vector<double> recovered(NumCDetPaddles, NAN);
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    int pixel = -1, bar = -1, group = -1, entries = 0, manualCut = 0, fitStatus = -1;
    std::string source, broadText, groupText, usedText, usedErrorText,
                sigmaText, correctionText, totalText;
    if (!(row >> pixel >> bar >> group >> entries >> manualCut >> source >>
          broadText >> groupText >> usedText >> usedErrorText >> sigmaText >>
          correctionText >> totalText >> fitStatus)) continue;
    if (pixel < 0 || pixel >= NumCDetPaddles) continue;
    const double correctionValue = ParseCDetFitResultNumber(correctionText);
    const double totalValue = ParseCDetFitResultNumber(totalText);
    if (std::isfinite(correctionValue) && std::isfinite(totalValue))
      recovered[pixel] = totalValue - correctionValue;
  }
  int recoveredCount = 0;
  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
    if (!std::isfinite(recovered[pixel])) continue;
    gPixelToffsetCorr[pixel] = recovered[pixel];
    ++recoveredCount;
  }
  if (recoveredCount != NumCDetPaddles) {
    std::cerr << "[CDet baseline recovery] ERROR: recovered " << recoveredCount
              << " / " << NumCDetPaddles << " offsets; no file written.\n";
    return false;
  }
  if (!WriteCalibrationConstants(baselineOutput.Data())) return false;
  std::cout << "[CDet baseline recovery] Reconstructed the analysis-time baseline in "
            << baselineOutput << " from " << recoveredCount << " pixel rows.\n";
  return true;
}

void buildCDetPixelLeTotReviewQueue(TString resultsFile, TString cutFile,
                                    TString queueOutput, int minEntries,
                                    double centroidDifference,
                                    double sigmaThreshold)
{
  gCDetPixelReviewQueue.clear();
  gCDetPixelReviewIndex = -1;
  gCDetPixelReviewCutFile = cutFile;
  if (minEntries < 1 || centroidDifference <= 0.0 || sigmaThreshold <= 0.0) {
    std::cerr << "[CDet pixel review] ERROR: invalid review threshold.\n";
    return;
  }

  std::ifstream input(resultsFile.Data());
  if (!input) {
    std::cerr << "[CDet pixel review] ERROR: could not open " << resultsFile << ".\n";
    return;
  }

  std::vector<bool> hasSavedCut(NumCDetPaddles, false);
  if (!cutFile.IsNull() && !gSystem->AccessPathName(cutFile)) {
    TFile cutInput(cutFile, "READ");
    if (!cutInput.IsZombie()) {
      for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
        const TString path = CDetPixelCutDirectory(pixel) + "/cut_le_vs_tot";
        hasSavedCut[pixel] = cutInput.Get(path) != nullptr;
      }
    }
    cutInput.Close();
  }

  int rowsRead = 0;
  int alreadyCut = 0;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    int pixel = -1, bar = -1, group = -1, entries = 0, manualCut = 0, fitStatus = -1;
    std::string source, broadText, groupText, usedText, usedErrorText,
                sigmaText, correctionText, totalText;
    if (!(row >> pixel >> bar >> group >> entries >> manualCut >> source >>
          broadText >> groupText >> usedText >> usedErrorText >> sigmaText >>
          correctionText >> totalText >> fitStatus)) continue;
    ++rowsRead;
    if (pixel < 0 || pixel >= NumCDetPaddles || IsUnusedPixel(pixel) ||
        entries < minEntries) continue;
    if (hasSavedCut[pixel] || manualCut != 0) {
      ++alreadyCut;
      continue;
    }

    const double broadMean = ParseCDetFitResultNumber(broadText);
    const double groupMeanValue = ParseCDetFitResultNumber(groupText);
    const double usedMeanValue = ParseCDetFitResultNumber(usedText);
    const double usedSigmaValue = ParseCDetFitResultNumber(sigmaText);
    const double broadGroupDifference =
        std::isfinite(broadMean) && std::isfinite(groupMeanValue)
            ? std::fabs(broadMean - groupMeanValue) : NAN;
    const double usedGroupDifference =
        std::isfinite(usedMeanValue) && std::isfinite(groupMeanValue)
            ? std::fabs(usedMeanValue - groupMeanValue) : NAN;

    double score = 0.0;
    std::vector<std::string> reasons;
    if (source != "individual_fit") {
      score += source == "individual_broad_fallback" ? 100.0 : 80.0;
      reasons.push_back("nonstandard fit source=" + source);
    }
    if (std::isfinite(broadGroupDifference) &&
        broadGroupDifference >= centroidDifference) {
      score += 20.0 + 5.0*std::min(10.0, broadGroupDifference);
      reasons.push_back(TString::Format("broad-group=%.2f ns", broadGroupDifference).Data());
    }
    if (std::isfinite(usedGroupDifference) &&
        usedGroupDifference >= centroidDifference) {
      score += 15.0 + 4.0*std::min(10.0, usedGroupDifference);
      reasons.push_back(TString::Format("used-group=%.2f ns", usedGroupDifference).Data());
    }
    if (std::isfinite(usedSigmaValue) && usedSigmaValue >= sigmaThreshold) {
      score += 10.0 + 4.0*std::min(8.0, usedSigmaValue - sigmaThreshold);
      reasons.push_back(TString::Format("sigma=%.2f ns", usedSigmaValue).Data());
    }
    if (score <= 0.0) continue;

    std::ostringstream reasonText;
    for (size_t reason = 0; reason < reasons.size(); ++reason) {
      if (reason) reasonText << "; ";
      reasonText << reasons[reason];
    }
    gCDetPixelReviewQueue.push_back({pixel, entries, score,
        broadGroupDifference, usedGroupDifference, usedSigmaValue,
        source, reasonText.str()});
  }
  input.close();

  std::sort(gCDetPixelReviewQueue.begin(), gCDetPixelReviewQueue.end(),
      [](const CDetPixelReviewCandidate& left,
         const CDetPixelReviewCandidate& right) {
        if (left.score != right.score) return left.score > right.score;
        return left.pixel < right.pixel;
      });

  std::ofstream output(queueOutput.Data());
  if (output) {
    output << "# Ranked conservative CDet LE-versus-TOT manual-review queue\n"
           << "# source results: " << resultsFile << "\n"
           << "# saved cuts skipped from: " << cutFile << "\n"
           << "# thresholds: entries >= " << minEntries
           << ", centroid difference >= " << centroidDifference
           << " ns, sigma >= " << sigmaThreshold << " ns\n"
           << "# rank pixel bar entries score source broad_group_ns used_group_ns sigma_ns reasons\n";
    for (size_t index = 0; index < gCDetPixelReviewQueue.size(); ++index) {
      const auto& candidate = gCDetPixelReviewQueue[index];
      output << index + 1 << " " << candidate.pixel << " " << candidate.pixel/16
             << " " << candidate.entries << " " << candidate.score << " "
             << candidate.source << " " << candidate.broadGroupDifference << " "
             << candidate.usedGroupDifference << " " << candidate.usedSigma
             << " \"" << candidate.reasons << "\"\n";
    }
  }

  std::cout << "[CDet pixel review]\n"
            << "  fit-result rows read: " << rowsRead << "\n"
            << "  pixels skipped because a manual cut already exists: " << alreadyCut << "\n"
            << "  conservative review candidates: " << gCDetPixelReviewQueue.size() << "\n"
            << "  ranked queue: " << queueOutput << "\n"
            << "Call printCDetPixelLeTotReviewQueue() to inspect the highest-ranked pixels,\n"
            << "then reviewNextCDetPixelLeTotCandidate() to begin.\n";
}

void printCDetPixelLeTotReviewQueue(int maximumRows)
{
  if (gCDetPixelReviewQueue.empty()) {
    std::cout << "[CDet pixel review] The queue is empty. Build it first.\n";
    return;
  }
  const int rows = std::min(std::max(1, maximumRows),
                            (int)gCDetPixelReviewQueue.size());
  std::cout << "[CDet pixel review] Highest-ranked " << rows << " of "
            << gCDetPixelReviewQueue.size() << " candidates:\n";
  for (int index = 0; index < rows; ++index) {
    const auto& candidate = gCDetPixelReviewQueue[index];
    std::cout << "  " << index + 1 << ". pixel " << candidate.pixel
              << " (bar " << candidate.pixel/16 << ", entries "
              << candidate.entries << ", score " << candidate.score << "): "
              << candidate.reasons << "\n";
  }
}

void reviewCurrentCDetPixelLeTotCandidate()
{
  if (gCDetPixelReviewIndex < 0 ||
      gCDetPixelReviewIndex >= (Long64_t)gCDetPixelReviewQueue.size()) {
    std::cout << "[CDet pixel review] No current candidate. Call the Next function.\n";
    return;
  }
  const auto& candidate = gCDetPixelReviewQueue[gCDetPixelReviewIndex];
  std::cout << "[CDet pixel review] Candidate " << gCDetPixelReviewIndex + 1
            << " / " << gCDetPixelReviewQueue.size() << ": pixel "
            << candidate.pixel << " -- " << candidate.reasons << "\n";
  editCDetPixelLeTotCut(candidate.pixel, gCDetPixelReviewCutFile);
}

void reviewNextCDetPixelLeTotCandidate()
{
  if (gCDetPixelReviewQueue.empty()) {
    std::cout << "[CDet pixel review] The queue is empty. Build it first.\n";
    return;
  }
  if (gCDetPixelReviewIndex + 1 >= (Long64_t)gCDetPixelReviewQueue.size()) {
    std::cout << "[CDet pixel review] End of queue.\n";
    return;
  }
  ++gCDetPixelReviewIndex;
  reviewCurrentCDetPixelLeTotCandidate();
}

void reviewPreviousCDetPixelLeTotCandidate()
{
  if (gCDetPixelReviewQueue.empty()) {
    std::cout << "[CDet pixel review] The queue is empty. Build it first.\n";
    return;
  }
  if (gCDetPixelReviewIndex <= 0) {
    gCDetPixelReviewIndex = 0;
    std::cout << "[CDet pixel review] At the beginning of the queue.\n";
  } else {
    --gCDetPixelReviewIndex;
  }
  reviewCurrentCDetPixelLeTotCandidate();
}

void plotPaddles(int bar = 29, double width = 1, double LeMin = 0, double LeMax = 60,
                 double TotCutLow = 0, double TotCutMax = 60,
                 double TotMin = 0, double TotMax = 60, double binLow = 0, double binHigh = 60){
  TH1::AddDirectory(kFALSE);
  if (bar < 0 || bar >= NumCDetPaddles/NumPaddles || width <= 0.0 || TotMin >= TotMax || binLow >= binHigh) {
    std::cerr << "[CDet paddle plots] ERROR: invalid bar, bin width, or histogram range.\n";
    return;
  }
  const int Nbins = (int)((binHigh-binLow)/width);
  const int paddle_base = bar*NumPaddles;
  (void)LeMin;
  (void)LeMax;
  (void)TotCutLow;
  (void)TotCutMax;

  std::vector<TF1*> fPixelGaussFit(NumPaddles, nullptr);
  hPaddleGoodLe.assign(NumPaddles, nullptr);
  hPaddleLEvsTOT.assign(NumPaddles, nullptr);
  TCanvas* cPaddles = new TCanvas(TString::Format("cPaddlesBar%d", bar), TString::Format("Bar %d Pixels", bar), 1200, 1000);
  cPaddles->Divide(4, 4, 0.001, 0.001);
  TCanvas* cPaddles2D = new TCanvas(TString::Format("cPaddles2DBar%d", bar), TString::Format("Bar %d Pixels LE versus TOT", bar), 1200, 1000);
  cPaddles2D->Divide(4, 4, 0.001, 0.001);

  for (int paddle = 0; paddle < NumPaddles; ++paddle){
    const int global_paddle = paddle + paddle_base;
    hPaddleGoodLe[paddle] = new TH1F(TString::Format("hBarGoodLe_Paddle%d", global_paddle), TString::Format("Good LE (Paddle %d)", global_paddle), Nbins, binLow, binHigh);
    hPaddleLEvsTOT[paddle] = new TH2F(TString::Format("hPaddleLEvsTOT_Paddle%d", global_paddle), TString::Format("Good LE vs TOT (Paddle %d);TOT [ns];LE [ns]", global_paddle), Nbins, TotMin, TotMax, Nbins, binLow, binHigh);

    for (double x : vPaddleGoodLe[global_paddle]) hPaddleGoodLe[paddle]->Fill(x);
    const double paddleMeanLe = hPaddleGoodLe[paddle]->GetMean();
    fPixelGaussFit[paddle] = new TF1(TString::Format("fPixelGaussFit_Paddle%d", global_paddle), "gaus", paddleMeanLe - 15.0, paddleMeanLe + 15.0);
    if (hPaddleGoodLe[paddle]->GetEntries() > 20){
      fPixelGaussFit[paddle]->SetParameters(hPaddleGoodLe[paddle]->GetMaximum(), paddleMeanLe, hPaddleGoodLe[paddle]->GetRMS());
      hPaddleGoodLe[paddle]->Fit(fPixelGaussFit[paddle], "RQ0");
    }

    for (size_t ihit = 0; ihit < vPaddleGoodLe[global_paddle].size(); ++ihit) hPaddleLEvsTOT[paddle]->Fill(vPaddleGoodTot[global_paddle][ihit], vPaddleGoodLe[global_paddle][ihit]);

    cPaddles->cd(paddle + 1);
    hPaddleGoodLe[paddle]->Draw();
    if (hPaddleGoodLe[paddle]->GetEntries() > 20) fPixelGaussFit[paddle]->Draw("SAME");
    if (kUnusedCDetPixels.count(global_paddle)) {
      TPaveText* flag = new TPaveText(0.82, 0.82, 0.95, 0.95, "NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetBorderSize(1);
      flag->AddText("UNUSED PIXEL");
      flag->Draw("same");
    }

    cPaddles2D->cd(paddle + 1);
    hPaddleLEvsTOT[paddle]->Draw("COLZ");
    if (kUnusedCDetPixels.count(global_paddle)) {
      TPaveText* flag = new TPaveText(0.82, 0.82, 0.95, 0.95, "NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetBorderSize(1);
      flag->AddText("UNUSED PIXEL");
      flag->Draw("same");
    }
  }
  cPaddles->Update();
  cPaddles->SaveAs(TString::Format("PaddleOffsetBar%d.pdf", bar));
  cPaddles2D->Update();
  cPaddles2D->SaveAs(TString::Format("PaddleOffset2DBar%d.pdf", bar));
}

void plotAllPaddles(double width = 1, double LeMin = 0, double LeMax = 60,
                    double TotMin = 0, double TotMax = 40,
                    double binLow = 0, double binHigh = 60, TString saveTag = "") {

  // Compatibility arguments: the existing routine plots already-selected
  // vectors with binLow/High and does not apply another LE or ToT selection.
  (void)LeMin;
  (void)LeMax;
  (void)TotMin;
  (void)TotMax;
  TH1::AddDirectory(kFALSE);

  const Bool_t previousBatchMode = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  constexpr int nBars          = 168;
  constexpr int pixelsPerBar   = 16;
  constexpr int nTotalPaddles  = nBars * pixelsPerBar;

  const int nBins = static_cast<int>((binHigh - binLow) / width);
  TString saveDir = "CDetPaddleTimes";
  if (!saveTag.IsNull() && saveTag != "") {
    saveDir += "/";
    saveDir += saveTag;
  }
  // Create output directory if it does not already exist.
  gSystem->mkdir(saveDir, kTRUE);

  if (vPaddleGoodLe.size() < nTotalPaddles) {
    std::cerr << "plotPaddles error: vPaddleGoodLe has only "
              << vPaddleGoodLe.size()
              << " entries, but at least "
              << nTotalPaddles
              << " are required."
              << std::endl;
    gROOT->SetBatch(previousBatchMode);
    return;
  }

  for (int bar = 0; bar < nBars; ++bar) {

    const int paddleBase = bar * pixelsPerBar;

    TString canvasName =
        TString::Format("cCDetBar%d_LE", bar);

    TString canvasTitle =
        TString::Format("CDet Bar %d Leading-Edge Times", bar);

    TCanvas* canvas =
        new TCanvas(canvasName, canvasTitle, 1200, 1000);

    canvas->Divide(4, 4, 0.001, 0.001);

    std::vector<TH1F*> histograms(pixelsPerBar, nullptr);
    std::vector<TF1*> fits(pixelsPerBar, nullptr);

    for (int localPaddle = 0;
         localPaddle < pixelsPerBar;
         ++localPaddle) {

      const int globalPaddle = paddleBase + localPaddle;

      TString histName =
          TString::Format("hBar%dGoodLe_Paddle%d",
                          bar,
                          globalPaddle);

      TString histTitle =
          TString::Format(
              "Bar %d, Paddle %d;Leading-edge time [ns];Counts",
              bar,
              globalPaddle);

      histograms[localPaddle] =
          new TH1F(histName,
                   histTitle,
                   nBins,
                   binLow,
                   binHigh);

      for (double le : vPaddleGoodLe[globalPaddle]) {
        histograms[localPaddle]->Fill(le);
      }

      canvas->cd(localPaddle + 1);

      histograms[localPaddle]->SetStats(kTRUE);
      histograms[localPaddle]->Draw();

      if (histograms[localPaddle]->GetEntries() > 20) {

        const double mean = histograms[localPaddle]->GetMean();
        const double rms  = histograms[localPaddle]->GetRMS();

        double fitLow  = mean - 15.0;
        double fitHigh = mean + 15.0;

        // Keep the fitting interval inside the histogram range.
        fitLow  = std::max(fitLow,  binLow);
        fitHigh = std::min(fitHigh, binHigh);

        TString fitName =
            TString::Format("fGauss_Bar%d_Paddle%d",
                            bar,
                            globalPaddle);

        fits[localPaddle] =
            new TF1(fitName, "gaus", fitLow, fitHigh);

        fits[localPaddle]->SetParameters(
            histograms[localPaddle]->GetMaximum(),
            mean,
            rms);

        histograms[localPaddle]->Fit(
            fits[localPaddle],
            "RQ0");

        fits[localPaddle]->Draw("SAME");
      }

      if (kUnusedCDetPixels.count(globalPaddle)) {

        TPaveText* flag =
            new TPaveText(0.72, 0.82, 0.96, 0.94, "NDC");

        flag->SetFillColor(kBlack);
        flag->SetTextColor(kWhite);
        flag->SetLineColor(kBlack);
        flag->SetBorderSize(1);
        flag->SetTextSize(0.035);
        flag->AddText("UNUSED");
        flag->Draw("SAME");
      }
    }

    TString outputName =
        TString::Format(
            "%s/CDetBar%03d_LE.pdf",
            saveDir.Data(),
            bar);

    canvas->SaveAs(outputName);

    // Clean up before moving to the next bar.
    for (int localPaddle = 0;
         localPaddle < pixelsPerBar;
         ++localPaddle) {

      delete fits[localPaddle];
      delete histograms[localPaddle];
    }

    delete canvas;
  }

  std::cout << "Saved CDet leading-edge plots for "
            << nBars
            << " bars in CDetPaddleTimes/"
            << std::endl;

  gROOT->SetBatch(previousBatchMode);
}

TCanvas *plotAllTDC(bool overwrite = false, double width = 1, double binLow = 0, double binHigh = 60, bool savePlots = false, TString saveTag = "", TString saveDir = "tdcPlots"){
  gLastCalibrationFitSucceeded = false;
  TH1::AddDirectory(kFALSE);
  int Nbins = (int)((binHigh-binLow)/width);

  //define histograms
  hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            Nbins, binLow, binHigh);
  hAllRawTe = new TH1F(TString::Format("hRawTe"),
            TString::Format("hRawTe"),
            Nbins, binLow, binHigh+TotBinHigh);
  hAllRawTot = new TH1F(TString::Format("hRawTot"),
            TString::Format("hRawTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllRawPMT = new TH1F(TString::Format("hRawPMT"),
            TString::Format("hRawPMT"),
            nTdc, 0, nTdc);
  hAllRawBar = new TH1F(TString::Format("hRawBar"),
            TString::Format("hRawBar"),
            168, 0, 168);
  hAllGoodLe = new TH1F(TString::Format("hAllGoodLe"),
            TString::Format("hAllGoodLe"),
            Nbins, binLow, binHigh);
  hAllGoodTe = new TH1F(TString::Format("hAllGoodTe"),
            TString::Format("hAllGoodTe"),
            Nbins, binLow, binHigh+TotBinHigh);
  hAllGoodTot = new TH1F(TString::Format("hAllGoodTot"),
            TString::Format("hAllGoodTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllGoodPMT = new TH1F(TString::Format("hAllGoodPMT"),
            TString::Format("hAllGoodPMT"),
            nTdc, 0, nTdc);
  hAllGoodBar = new TH1F(TString::Format("hAllGoodBar"),
            TString::Format("hAllGoodBar"),
            168, 0, 168);

  // Per-pixel good leading-edge histograms (pixelID = GoodElID)
  hPaddleGoodLe.assign(NumCDetPaddles, nullptr);
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    hPaddleGoodLe[pixelID] = new TH1F(TString::Format("hPixelGoodLe_Pixel%d", pixelID), TString::Format("Good LE (logical pixel ID %d)", pixelID), Nbins, binLow, binHigh);
    for (double x : vPaddleGoodLe[pixelID]) hPaddleGoodLe[pixelID]->Fill(x);
  }

  // Per-bar good leading-edge histograms (16 logical pixels per bar).
  // These feed the four 42-bar diagnostic canvases saved below.
  hBarGoodLe.assign(NumPMTs, nullptr);
  for (int bar = 0; bar < NumPMTs; ++bar) {
    hBarGoodLe[bar] = new TH1F(
        TString::Format("hBarGoodLe_Bar%d", bar),
        TString::Format("Good LE (bar %d);Leading-edge time [ns];Counts", bar),
        Nbins, binLow, binHigh);
    if (bar < (int)vBarGoodLe.size()) {
      for (double x : vBarGoodLe[bar]) hBarGoodLe[bar]->Fill(x);
    }
  }

  // fill necessary histograms from vectors
  for (double x : vAllRawLe)   hAllRawLe->Fill(x);
  for (double x : vAllRawTe)   hAllRawTe->Fill(x);
  for (double x : vAllRawTot)  hAllRawTot->Fill(x);
  for (double x : vAllRawPMT)  hAllRawPMT->Fill(x);
  for (double x : vAllRawBar)  hAllRawBar->Fill(x);

  for (double x : vAllGoodLe)  hAllGoodLe->Fill(x);
  for (double x : vAllGoodTe)  hAllGoodTe->Fill(x);
  for (double x : vAllGoodTot) hAllGoodTot->Fill(x);
  for (double x : vAllGoodPMT) hAllGoodPMT->Fill(x);
  for (double x : vAllGoodBar) hAllGoodBar->Fill(x);

  // ------------------------------------------------------------
  // Compute per-pixel residual offsets relative to the global good-LE mean
  // residual[pixel] = mean(all) - mean(pixel)
  // This is the correction to ADD on top of whatever is already applied.
  // Pixels with < 20 entries and known unused pixels get residual = 0.
  // ------------------------------------------------------------
  const double meanAllGoodLe = hAllGoodLe->GetMean();
  if (hAllGoodLe->GetEntries() < 20 || !std::isfinite(meanAllGoodLe)) {
    std::cerr << "[CDet] ERROR: insufficient good-hit statistics for pixel-offset fitting.\n";
    return nullptr;
  }
  TF1 *fGaus_all = new TF1("fGaus_all", "gaus", binLow + 20, binHigh);
  fGaus_all->SetParameters(hAllGoodLe->GetMaximum(),meanAllGoodLe,hAllGoodLe->GetRMS());
  const int globalFitStatus = hAllGoodLe->Fit(fGaus_all, "RQS0");
  double gausFitAllGoodLeMean = fGaus_all->GetParameter(1);
  if (globalFitStatus != 0 || !std::isfinite(gausFitAllGoodLeMean)) {
    std::cerr << "[CDet] ERROR: global pixel-offset fit failed; constants were not written.\n";
    return nullptr;
  }
  std::vector<double> vPixelResidualOffset(NumCDetPaddles, 0.0);
  std::vector<int> vPixelOffsetNhits(NumCDetPaddles, 0);

  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    if (!hPaddleGoodLe[pixelID]) continue;
    const int nent = hPaddleGoodLe[pixelID]->GetEntries();
    vPixelOffsetNhits[pixelID] = nent;
    if (!IsUnusedPixel(pixelID) && nent >= 20) {
      //use gaussian fit to get mean value of histogram
      TF1 *fGaus = new TF1(TString::Format("fGaus_pixel%d", pixelID), "gaus", binLow + 20, binHigh);
      fGaus->SetParameters(hPaddleGoodLe[pixelID]->GetMaximum(), meanAllGoodLe, hPaddleGoodLe[pixelID]->GetRMS());
      const int pixelFitStatus = hPaddleGoodLe[pixelID]->Fit(fGaus,"RQS0");
      double gausFitGoodLeMean = fGaus->GetParameter(1);
      if (pixelFitStatus != 0 || !std::isfinite(gausFitGoodLeMean) ||
          gausFitGoodLeMean > 60 || gausFitGoodLeMean < 0)
        vPixelResidualOffset[pixelID] = meanAllGoodLe - hPaddleGoodLe[pixelID]->GetMean();
      else vPixelResidualOffset[pixelID] = gausFitAllGoodLeMean - gausFitGoodLeMean;
    } else {
      vPixelResidualOffset[pixelID] = 0.0;
    }
  }

  // ------------------------------------------------------------
  // Build total offsets to write out.
  //
  // Case 1: no file was loaded before run
  //   total = residual
  //
  // Case 2: file was loaded and overwrite==true
  //   total = old + residual
  //
  // Case 3: file was loaded and overwrite==false
  //   do not write file
  // ------------------------------------------------------------
  std::vector<double> vPixelTotalOffset(NumCDetPaddles, 0.0);

  if ((int)gPixelToffsetCorr.size() == NumCDetPaddles) {
    for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
      vPixelTotalOffset[pixelID] = gPixelToffsetCorr[pixelID] + vPixelResidualOffset[pixelID];
    }
  } else {
    for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
      vPixelTotalOffset[pixelID] = vPixelResidualOffset[pixelID];
    }
  }

  // Write file if:
  //   - no offsets were loaded originally, OR
  //   - overwrite requested
  const bool shouldWriteOffsets = (!gPixelToffsetLoaded) || overwrite;

  if (shouldWriteOffsets) {
    // Update global vector to the new total values
    gPixelToffsetCorr = vPixelTotalOffset;
    gPixelToffsetNhits = vPixelOffsetNhits;

    if (WriteCalibrationConstants(gCalibrationFile)) {
      if (overwrite && gPixelToffsetLoaded) {
        std::cout << "[CDet] Updated pixel offsets in calibration file using total = old + residual\n";
      } else {
        std::cout << "[CDet] Wrote pixel offsets to calibration file from current pass\n";
      }
    } else return nullptr;
  }

  // Plot offsets vs logical pixel ID
  // Show residual offsets by default, since that is what current spectra imply.
  TCanvas* cOffsets = new TCanvas("cPixelGoodLeOffsets", "Mean LE offsets vs logical pixel ID", 50, 900, 1200, 500);
  std::vector<double> xPixel(NumCDetPaddles), yOff(NumCDetPaddles);
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    xPixel[pixelID] = pixelID;
    yOff[pixelID] = vPixelResidualOffset[pixelID];
  }
  TGraph* grPixelOffsets = new TGraph(NumCDetPaddles, xPixel.data(), yOff.data());
  grPixelOffsets->SetName("grPixelGoodLeOffsets");
  grPixelOffsets->SetTitle(TString::Format("Residual LE offsets vs logical pixel ID;Logical pixel ID;#mu_{all} - #mu_{pixel} (ns)  (global mean = %.3f ns)", gausFitAllGoodLeMean));
  grPixelOffsets->SetMarkerStyle(20);
  grPixelOffsets->SetMarkerSize(0.4);
  grPixelOffsets->Draw("AP");

  TCanvas *caa = new TCanvas("All TDC", "All TDC", 50,50,800,800);
  caa->Divide(2,3,0.01,0.01,0);
  TCanvas *caaa = new TCanvas("All Chan", "All Chan", 850,50,800,800);
  caaa->Divide(2,2,0.01,0.01,0);

  caa->cd(1);
  hAllRawLe->SetFillColor(kRed);
  hAllRawLe->SetMinimum(0.0);
  hAllRawLe->Draw();

  caa->cd(2);
  hAllGoodLe->SetFillColor(kBlue);
  hAllGoodLe->SetMinimum(0.0);
  hAllGoodLe->Draw();

  caa->cd(3);
  hAllRawTe->SetFillColor(kRed);
  hAllRawTe->SetMinimum(0.0);
  hAllRawTe->Draw();

  caa->cd(4);
  hAllGoodTe->SetFillColor(kBlue);
  hAllGoodTe->SetMinimum(0.0);
  hAllGoodTe->Draw();

  caa->cd(5);
  hAllRawTot->SetFillColor(kRed);
  hAllRawTot->Draw();

  caa->cd(6);
  hAllGoodTot->SetFillColor(kBlue);
  hAllGoodTot->Draw();

  caaa->cd(1);
  hAllRawPMT->SetFillColor(kRed);
  hAllRawPMT->Draw();

  caaa->cd(2);
  hAllGoodPMT->SetFillColor(kBlue);
  hAllGoodPMT->Draw();

  caaa->cd(3);
  hAllRawBar->SetFillColor(kRed);
  hAllRawBar->Draw();

  caaa->cd(4);
  hAllGoodBar->SetFillColor(kBlue);
  hAllGoodBar->Draw();

  // ------------------------------------------------------------
  // Draw per-bar Good LE histograms on 4 canvases (42 per canvas)
  // ------------------------------------------------------------
  auto drawBarRange = [&](const char* cname, const char* ctitle, int barStart, int nBars){
    TCanvas* c = new TCanvas(cname, ctitle, 1400, 900);
    c->Divide(7, 6, 0.001, 0.001);
    for(int i = 0; i < nBars; ++i){
      int bar = barStart + i;
      c->cd(i+1);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.02);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.08);

      if(bar < 0 || bar >= (int)hBarGoodLe.size() || !hBarGoodLe[bar]) continue;

      TH1F* h = hBarGoodLe[bar];
      h->SetStats(0);
      h->GetXaxis()->SetTitleSize(0.07);
      h->GetXaxis()->SetLabelSize(0.06);
      h->GetYaxis()->SetTitleSize(0.07);
      h->GetYaxis()->SetLabelSize(0.06);
      h->GetXaxis()->SetTitleOffset(0.85);
      h->GetYaxis()->SetTitleOffset(0.85);
      h->Draw("HIST");
      const int nent = h->GetEntries();

      if (nent >= 20) {
        int maxBin = h->GetMaximumBin();

        TF1 *fGaus = new TF1(TString::Format("fGaus_bar%d_%s", bar, cname),"gaus", 20, 60);
        fGaus->SetParameters(h->GetMaximum(), meanAllGoodLe, h->GetRMS());

        fGaus->SetLineColor(kRed);
        fGaus->SetLineWidth(2);

        int fitStatus = h->Fit(fGaus, "RQ0");

        if (fitStatus == 0) {
          fGaus->Draw("SAME");

          TLatex lat;
          lat.SetNDC();
          lat.SetTextSize(0.055);
          lat.DrawLatex(
              0.18, 0.82,
              TString::Format("#mu = %.2f ns", fGaus->GetParameter(1))
          );
        }
      }
    }

    return c;
  };

  drawBarRange("cBarGoodLe_L1L", "Good LE by Bar: Layer 1 Left (Bars 1-42)",     0,  42);
  drawBarRange("cBarGoodLe_L1R", "Good LE by Bar: Layer 1 Right (Bars 43-84)",   42,  42);
  drawBarRange("cBarGoodLe_L2L", "Good LE by Bar: Layer 2 Left (Bars 85-126)",   84,  42);
  drawBarRange("cBarGoodLe_L2R", "Good LE by Bar: Layer 2 Right (Bars 127-168)",126,  42);
  // 4/15/2026 -- B. Spaude added this to save these plots to folder. Easier to confirm runs okay
  if (savePlots) {
    if (gSystem->AccessPathName(saveDir)) {
      gSystem->mkdir(saveDir, true);
    }

    TString suffix = saveTag;
    if (suffix.Length() > 0) suffix = "_" + suffix;

    cOffsets->SaveAs(TString::Format("%s/cBarGoodLeOffsets%s.pdf",
                                     saveDir.Data(), suffix.Data()));

    caa->SaveAs(TString::Format("%s/AllTDC%s.pdf",
                                saveDir.Data(), suffix.Data()));

    caaa->SaveAs(TString::Format("%s/AllChan%s.pdf",
                                 saveDir.Data(), suffix.Data()));

    TCanvas *c1 = (TCanvas*)gROOT->FindObject("cBarGoodLe_L1L");
    TCanvas *c2 = (TCanvas*)gROOT->FindObject("cBarGoodLe_L1R");
    TCanvas *c3 = (TCanvas*)gROOT->FindObject("cBarGoodLe_L2L");
    TCanvas *c4 = (TCanvas*)gROOT->FindObject("cBarGoodLe_L2R");

    if (c1) c1->SaveAs(TString::Format("%s/cBarGoodLe_L1L%s.pdf",
                                       saveDir.Data(), suffix.Data()));
    if (c2) c2->SaveAs(TString::Format("%s/cBarGoodLe_L1R%s.pdf",
                                       saveDir.Data(), suffix.Data()));
    if (c3) c3->SaveAs(TString::Format("%s/cBarGoodLe_L2L%s.pdf",
                                       saveDir.Data(), suffix.Data()));
    if (c4) c4->SaveAs(TString::Format("%s/cBarGoodLe_L2R%s.pdf",
                                       saveDir.Data(), suffix.Data()));
  }

  gLastCalibrationFitSucceeded = true;
  return caa;
}


void plot2DrefVsLE(double width = 1, double tmin=0, double tmax=60){
  int TDCBinNum = (int)((tmax-tmin)/width);
  TH2D *h2 = new TH2D("h2","t_ref vs t_scint; t_ref; t_scint",TDCBinNum, tmin, tmax, TDCBinNum, tmin, tmax);

  // Use reference times that are aligned with the GOOD-event vectors
  const size_t Nev = std::min(vGoodRefRawLe.size(), vGoodLe.size());
  for (size_t ev = 0; ev < Nev; ev++) { // iterate through good-selected events
    const double t_ref = vGoodRefRawLe[ev];
    if (std::isnan(t_ref)) continue; // if no ref recorded for this entry, skip

    const size_t Nhits = vGoodLe[ev].size();
    for (size_t ihit = 0; ihit < Nhits; ihit++) {
      const double t_det = vGoodLe[ev][ihit];
      h2->Fill(t_ref,t_det);
    }
  }
  TCanvas *crefVSle = new TCanvas("crefVSle", "CDet ref vs LE",900,700);
  h2->Draw("COLZ");
}

TH1* SubtractFitFromHist(const TH1* hIn, TF1* fFit, const char* outName = nullptr,
                         bool clampNegToZero = true, int firstBin = 1, int lastBin = -1) {
  if (!hIn || !fFit) {
    std::cerr << "SubtractFitFromHist ERROR: null input.\n";
    return nullptr;
  }

  TH1* hSub = (TH1*)hIn->Clone(outName ? outName : (TString(hIn->GetName()) + "_sub").Data());
  hSub->SetDirectory(nullptr);

  if (lastBin < 0) lastBin = hSub->GetNbinsX();
  firstBin = std::max(firstBin, 1);
  lastBin  = std::min(lastBin, hSub->GetNbinsX());

  for (int ibin = firstBin; ibin <= lastBin; ++ibin) {
    const double x  = hSub->GetBinCenter(ibin);
    const double bw = hSub->GetBinWidth(ibin);

    const double content = hSub->GetBinContent(ibin);
    const double err     = hSub->GetBinError(ibin);

    // Average fit value over the bin (Integral/binwidth) so subtraction matches histogram binning
    const double fAvg = fFit->Integral(x - bw/2.0, x + bw/2.0) / bw;

    double newContent = content - fAvg;
    if (clampNegToZero && newContent < 0) newContent = 0;

    hSub->SetBinContent(ibin, newContent);
    hSub->SetBinError(ibin, err); // keep original errors (common choice)
  }

  // If you want to subtract everywhere else too, call with firstBin=1,lastBin=nbins.
  return hSub;
}

void plotCDetLayersTimeComp(bool overwrite = false, int pixelBase = 416, double Width = 1,
                            double diffMinCut = -15, double diffMaxCut = 15,
                            double xdiffMinCut = -0.1, double xdiffMaxCut = 0.1,
                            double LeMin = 0.02, double LeMax = 60,
                            double TotMinCut = 0, double TotMaxCut = 70,
                            double DiffMin = -20, double DiffMax = 20,
                            double CDetMin = 0, double CDetMax = 60,
                            double CDetTotMin = 0, double CDetTotMax = 80,
                            double ECalMin = -40, double ECalMax = 40,
                            double tdiffECalCDetMin = -60, double tdiffECalCDetMax = 30,
                            bool allowMultiplePairs = true,
                            double XBinWidth = 0.005, double XMin = -1.5, double XMax = 1.5,
                            double ZBinWidth = 0.01, double ZMin = 0.0, double ZMax = 7.0,
                            double HCalMin = -10.0, double HCalMax = 10.0,
                            bool displayBestHits = false,
                            double bestHitPeakMean = std::numeric_limits<double>::quiet_NaN(),
                            double bestHitPeakSigma = std::numeric_limits<double>::quiet_NaN(),
                            double bestHitNSigma = 3.0,
                            double yCorrectionRefractiveIndex = 0.0){
  gLastCalibrationFitSucceeded = false;
  gLastECalFixedEffectsSlope = std::numeric_limits<double>::quiet_NaN();
  gLastECalFixedEffectsSlopeError = std::numeric_limits<double>::quiet_NaN();
  gLastECalFixedEffectsValid = false;
  gCDetAcceptedPairMeanTimes.clear();
  gCDetProjectedHalfBarPairMeanTimes.clear();

  TH1::AddDirectory(kFALSE);

  if (pixelBase < 0 || pixelBase >= NumCDetPaddles/2) {
    std::cerr << "[plotCDetLayersTimeComp] ERROR: Layer 1 pixel ID must be in [0, " << NumCDetPaddles/2 - 1 << "].\n";
    return;
  }
  if (Width <= 0.0 || CDetMin >= CDetMax || XBinWidth <= 0.0 || XMin >= XMax || ZBinWidth <= 0.0 || ZMin >= ZMax || HCalMin >= HCalMax) {
    std::cerr << "[plotCDetLayersTimeComp] ERROR: invalid timing, x-position, or z-position binning or range.\n";
    return;
  }

  int TDCBinNum = (int)((DiffMax-DiffMin)/Width);
  int NADCBins = (int)((ECalMax-ECalMin)/4); //4ns bins for ECal, since fADC 4ns resolution
  const int NPairedLeBins = std::max(1, (int)((CDetMax-CDetMin)/Width));
  const int NXBins = std::max(1, (int)((XMax-XMin)/XBinWidth));
  const int NZBins = std::max(1, (int)((ZMax-ZMin)/ZBinWidth));
  const bool applyYPropagationDiagnostic = yCorrectionRefractiveIndex > 0.0;
  const double vacuumLightSpeedMPerNs = 0.299792458;
  const double yPropagationSpeedMPerNs = applyYPropagationDiagnostic
      ? vacuumLightSpeedMPerNs/yCorrectionRefractiveIndex : 0.0;
  const double yPropagationSlopeMagnitude = applyYPropagationDiagnostic
      ? 1.0/yPropagationSpeedMPerNs : 0.0;
  const int selectedLayer1BarBase = (pixelBase/NumPaddles)*NumPaddles;
  const int selectedLayer2BarBase = selectedLayer1BarBase + NumCDetPaddles/2;
  const int selectedBarNumber = selectedLayer1BarBase/NumPaddles;
  std::cout << "[plotCDetLayersTimeComp] Layer 1 pixel " << pixelBase
            << " normalized to bar " << selectedBarNumber
            << " (pixels " << selectedLayer1BarBase << "-" << selectedLayer1BarBase + NumPaddles - 1 << ").\n";

  TH1D* hCDetTimeDiff = new TH1D("hCDetTimeDiff", "CDet Layer Time Difference; Time Difference (ns);Counts", TDCBinNum, DiffMin, DiffMax);
  TH1D* hCDetTimeDiffYCorrected = nullptr;
  TH1D* hCDetXDiff = new TH1D("hCDetXDiff", "CDet Layer X Difference; X Diff (m);Counts", 600,-1.5,1.5);
  TH1D* hCDetLe1 = new TH1D("hCDetLe1", "CDet Layer 1 Good Time;Layer 1 LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetLe2 = new TH1D("hCDetLe2", "CDet Layer 2 Good Time;Layer 2 LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetLe1YCorrected = nullptr;
  TH1D* hCDetLe2YCorrected = nullptr;
  TH1D* hCDetProjectedHalfBarLe1 = new TH1D(
      "hCDetProjectedHalfBarLe1",
      "Projected-half-bar hits, Layer 1;Corrected LE time (ns);Hits",
      TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetProjectedHalfBarLe2 = new TH1D(
      "hCDetProjectedHalfBarLe2",
      "Projected-half-bar hits, Layer 2;Corrected LE time (ns);Hits",
      TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetProjectedHalfBarPairMean = new TH1D(
      "hCDetProjectedHalfBarPairMean",
      "Projected-half-bar pair mean time;Corrected (t_{L1}+t_{L2})/2 (ns);Pairs",
      TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDetProjectedHalfBarLe2vs1 = new TH2D(
      "hCDetProjectedHalfBarLe2vs1",
      "Projected-half-bar Layer 2 vs Layer 1 time;Layer 1 corrected LE (ns);Layer 2 corrected LE (ns)",
      TDCBinNum, CDetMin, CDetMax, TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetProjectedHalfBarLe1YCorrected = nullptr;
  TH1D* hCDetProjectedHalfBarLe2YCorrected = nullptr;
  TH1D* hCDetProjectedHalfBarPairMeanYCorrected = nullptr;
  TH2D* hCDetProjectedHalfBarLe2vs1YCorrected = nullptr;
  TH2D* hCDetLe2vs1 = new TH2D("hCDetLe2vs1", "CDet Layer 2 Time vs CDet Layer 1 Time;Layer 1 LE (ns);Layer 2 LE (ns)",TDCBinNum,CDetMin,CDetMax,TDCBinNum,CDetMin,CDetMax);
  TH2D* hCDetLe2vs1YCorrected = nullptr;
  if (applyYPropagationDiagnostic) {
    hCDetTimeDiffYCorrected = new TH1D(
        "hCDetTimeDiffYCorrected",
        "CDet Layer Time Difference after y correction;t_{L2}^{y}-t_{L1}^{y} (ns);Counts",
        TDCBinNum, DiffMin, DiffMax);
    hCDetLe1YCorrected = new TH1D(
        "hCDetLe1YCorrected",
        "CDet Layer 1 Good Time after y correction;Layer 1 corrected LE (ns);Counts",
        TDCBinNum, CDetMin, CDetMax);
    hCDetLe2YCorrected = new TH1D(
        "hCDetLe2YCorrected",
        "CDet Layer 2 Good Time after y correction;Layer 2 corrected LE (ns);Counts",
        TDCBinNum, CDetMin, CDetMax);
    hCDetLe2vs1YCorrected = new TH2D(
        "hCDetLe2vs1YCorrected",
        "CDet Layer 2 vs Layer 1 time after y correction;Layer 1 corrected LE (ns);Layer 2 corrected LE (ns)",
        TDCBinNum, CDetMin, CDetMax, TDCBinNum, CDetMin, CDetMax);
    hCDetProjectedHalfBarLe1YCorrected = new TH1D(
        "hCDetProjectedHalfBarLe1YCorrected",
        "Projected-half-bar hits, Layer 1 after y correction;Corrected LE time (ns);Hits",
        TDCBinNum, CDetMin, CDetMax);
    hCDetProjectedHalfBarLe2YCorrected = new TH1D(
        "hCDetProjectedHalfBarLe2YCorrected",
        "Projected-half-bar hits, Layer 2 after y correction;Corrected LE time (ns);Hits",
        TDCBinNum, CDetMin, CDetMax);
    hCDetProjectedHalfBarPairMeanYCorrected = new TH1D(
        "hCDetProjectedHalfBarPairMeanYCorrected",
        "Projected-half-bar pair mean after y correction;Corrected (t_{L1}^{y}+t_{L2}^{y})/2 (ns);Pairs",
        TDCBinNum, CDetMin, CDetMax);
    hCDetProjectedHalfBarLe2vs1YCorrected = new TH2D(
        "hCDetProjectedHalfBarLe2vs1YCorrected",
        "Projected-half-bar Layer 2 vs Layer 1 after y correction;Layer 1 corrected LE (ns);Layer 2 corrected LE (ns)",
        TDCBinNum, CDetMin, CDetMax, TDCBinNum, CDetMin, CDetMax);
  }
  TH2D* hCDetTot2vs1 = new TH2D("hCDetTot2vs1", "CDet Layer 2 ToT vs CDet Layer 1 ToT;Layer 1 ToT (ns);Layer 2 ToT (ns)",TDCBinNum,CDetTotMin,CDetTotMax,TDCBinNum,CDetTotMin,CDetTotMax);
  TH2D* h2CDetx2VsCDetx1 = new TH2D("h2CDetx2VsCDetx1", "CDet Layer 2 x vs CDet Layer 1 x;CDet Layer 1 x (m);CDet Layer 2 x (m)",600,-1.5,1.5,600,-1.5,1.5);
  TH1I* hNpairPerEvent = new TH1I("hNpairPerEvent", "CDet accepted pairs per event;N_{pairs};Events", 20, 0, 20);
  TH1D* hCDetBarLe1 = new TH1D("hCDetBarLe1", TString::Format("CDet Bar 1L%d Good Time;Layer 1 LE (ns);Counts", selectedBarNumber), TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetBarLe2 = new TH1D("hCDetBarLe2", TString::Format("CDet Bar 2L%d Good Time;Layer 2 LE (ns);Counts", selectedBarNumber), TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet1BarLeVsTot = new TH2D("hCDet1BarLeVsTot", TString::Format("CDet Bar 1L%d LE vs TOT;TOT (ns);LE (ns)", selectedBarNumber), TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet2BarLeVsTot = new TH2D("hCDet2BarLeVsTot", TString::Format("CDet Bar 2L%d LE vs TOT;TOT (ns);LE (ns)", selectedBarNumber), TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet1BarVerticalLeVsTot = new TH2D(
      "hCDet1BarVerticalLeVsTot",
      TString::Format("Bar 1L%d vertical ECal-time component (-3<t_{ECal}<1, 23<t_{CDet}<40 ns);ToT (ns);Corrected LE (ns)", selectedBarNumber),
      TDCBinNum, CDetTotMin, CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet1BarDiagonalLeVsTot = new TH2D(
      "hCDet1BarDiagonalLeVsTot",
      TString::Format("Bar 1L%d diagonal component (1<t_{ECal}<18, |t_{CDet}-27-0.73t_{ECal}|<3 ns);ToT (ns);Corrected LE (ns)", selectedBarNumber),
      TDCBinNum, CDetTotMin, CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDetAllLowerVerticalLeVsTot = new TH2D(
      "hCDetAllLowerVerticalLeVsTot",
      "All CDet bars: lower vertical component (-5<t_{ECal}<1, 23<t_{CDet}<30 ns);ToT (ns);Corrected LE (ns)",
      TDCBinNum, CDetTotMin, CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDetAllUpperVerticalLeVsTot = new TH2D(
      "hCDetAllUpperVerticalLeVsTot",
      "All CDet bars: upper vertical component (-5<t_{ECal}<1, 30<t_{CDet}<40 ns);ToT (ns);Corrected LE (ns)",
      TDCBinNum, CDetTotMin, CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDetAllDiagonalLeVsTot = new TH2D(
      "hCDetAllDiagonalLeVsTot",
      "All CDet bars: diagonal component (t_{ECal}>1, t_{CDet}>30, |t_{CDet}-27-0.73t_{ECal}|<3 ns);ToT (ns);Corrected LE (ns)",
      TDCBinNum, CDetTotMin, CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet1LeVsTot = new TH2D("hCDet1LeVsTot", "CDet Layer 1 Le vs Tot;Tot (ns);LE (ns)", TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet2LeVsTot = new TH2D("hCDet2LeVsTot", "CDet Layer 2 Le vs Tot;Tot (ns);LE (ns)", TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hECalVsCDetDt = new TH2D("hECalVsCDetDt", "ECal Time vs CDet dt;ECal ADC Time (ns);CDet dt_12 (ns)", NADCBins, ECalMin, ECalMax,TDCBinNum,DiffMin,DiffMax);
  TH2D* hECalVsCDetDtSingle = new TH2D("hECalVsCDetDtSingle", "CDet Single dt vs ECal Time;ECal ADC Time (ns);CDet dt_12 (ns)", NADCBins, ECalMin, ECalMax, TDCBinNum,DiffMin,DiffMax);
  TH2D* hECalVsCDetT = new TH2D("hECalVsCDetT", "CDet t vs ECal Time;ECal ADC Time (ns);CDet t (ns)", NADCBins, ECalMin, ECalMax,TDCBinNum,CDetMin,CDetMax);
  hECalVsCDetT->SetDirectory(nullptr);
  gCDetPairedMeanTimeVsECal = hECalVsCDetT;
  TH2D* hECalVsCDetTL1 = new TH2D("hECalVsCDetTL1", "CDet Layer 1 t vs ECal Time;ECal ADC Time (ns);Layer 1 t (ns)", NADCBins, ECalMin, ECalMax,TDCBinNum,CDetMin,CDetMax);
  TH2D* hECalVsCDetTL2 = new TH2D("hECalVsCDetTL2", "CDet Layer 2 t vs ECal Time;ECal ADC Time (ns);Layer 2 t (ns)", NADCBins, ECalMin, ECalMax,TDCBinNum,CDetMin,CDetMax);
  const double beforeP1TimeMin = CDetMin +
      std::min(0.0, std::min(gECalFitP1*ECalMin, gECalFitP1*ECalMax)) - 5.0;
  const double beforeP1TimeMax = CDetMax +
      std::max(0.0, std::max(gECalFitP1*ECalMin, gECalFitP1*ECalMax)) + 5.0;
  const int beforeP1TimeBins = std::max(
      1, (int)std::ceil((beforeP1TimeMax-beforeP1TimeMin)/Width));
  TH2D* hECalVsCDetTBeforeP1 = new TH2D(
      "hECalVsCDetTBeforeP1",
      "CDet pair time before p1 correction;ECal ADC Time (ns);CDet pair mean time (ns)",
      NADCBins, ECalMin, ECalMax, beforeP1TimeBins, beforeP1TimeMin, beforeP1TimeMax);
  TH2D* hECalVsCDetTL1BeforeP1 = new TH2D(
      "hECalVsCDetTL1BeforeP1",
      "CDet Layer 1 time before p1 correction;ECal ADC Time (ns);Layer 1 time (ns)",
      NADCBins, ECalMin, ECalMax, beforeP1TimeBins, beforeP1TimeMin, beforeP1TimeMax);
  TH2D* hECalVsCDetTL2BeforeP1 = new TH2D(
      "hECalVsCDetTL2BeforeP1",
      "CDet Layer 2 time before p1 correction;ECal ADC Time (ns);Layer 2 time (ns)",
      NADCBins, ECalMin, ECalMax, beforeP1TimeBins, beforeP1TimeMin, beforeP1TimeMax);
  const int NHCalBins = std::max(1, (int)((HCalMax-HCalMin)/0.5));
  TH2D* hHCalVsCDetT = new TH2D("hHCalVsCDetT",
      "Corrected CDet pair time vs HCal time;HCal ADC time (ns);Corrected CDet pair mean time (ns)",
      NHCalBins, HCalMin, HCalMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hECalVsHCalAccepted = new TH2D("hECalVsHCalAccepted",
      "ECal vs HCal time for accepted CDet pairs;ECal ADC time (ns);HCal ADC time (ns)",
      NADCBins, ECalMin, ECalMax, NHCalBins, HCalMin, HCalMax);
  const int NECalYBins = 120;
  const double ECalYMin = -0.6;
  const double ECalYMax = 0.6;
  TH2D* hCDetTvsECalY[NumLayers][NumSides][3] = {};
  const char* sideName[NumSides] = {"Left", "Right"};
  const char* sectionGroupName[3] = {"Sections 1+6", "Sections 2+5", "Sections 3+4"};
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int side = 0; side < NumSides; ++side) {
      for (int group = 0; group < 3; ++group) {
        hCDetTvsECalY[layer][side][group] = new TH2D(
            TString::Format("hCDetL%d%sSectionGroup%dTvsECalY", layer + 1, sideName[side], group + 1),
            TString::Format("Corrected CDet Layer %d %s, %s;ECal y position (m);Corrected CDet LE time (ns)",
                            layer + 1, sideName[side], sectionGroupName[group]),
            NECalYBins, ECalYMin, ECalYMax, NPairedLeBins, CDetMin, CDetMax);
      }
    }
  }
  TH2D* hECalVsSelectedBarT = new TH2D("hECalVsSelectedBarT", TString::Format("CDet Layer 1 bar %d t vs ECal Time;ECal ADC Time (ns);Layer 1 t (ns)", selectedBarNumber), NADCBins, ECalMin, ECalMax,TDCBinNum,CDetMin,CDetMax);
  TH2D* hECalVsSelectedBarTBeforeP1 = new TH2D(
      "hECalVsSelectedBarTBeforeP1",
      TString::Format("CDet Layer 1 bar %d time before p1 correction;ECal ADC Time (ns);Layer 1 time (ns)", selectedBarNumber),
      NADCBins, ECalMin, ECalMax, beforeP1TimeBins, beforeP1TimeMin, beforeP1TimeMax);
  TH2D* hECalVsCDetTSingle = new TH2D("hECalVsCDetTSingle", "CDet Single t vs ECal Time;ECal ADC Time (ns);CDet t (ns)", NADCBins, ECalMin, ECalMax, TDCBinNum,CDetMin,CDetMax);

  TH2D* hCDet1IDvs2ID = new TH2D("hCDet1IDvs2ID", "CDet Front Paddle vs Back Paddle ID;Back Paddle;Front Paddle", 1344, 1343.5, 2687.5, 1344, -0.5, 1343.5);
  TH2D* hCDetTimeDiffvsx1 = new TH2D("hCDetTimeDiffvsx1", "CDet Time Diff vs x1 position; x1 pos (m);CDet dt (ns)", 600,-1.5,1.5,TDCBinNum,DiffMin,DiffMax);
  TH2D* hCDetTimeDiffvsx2 = new TH2D("hCDetTimeDiffvsx2", "CDet Time Diff vs x2 position; x2 pos (m);CDet dt (ns)", 600,-1.5,1.5,TDCBinNum,DiffMin,DiffMax);

  TH2D* hCDetTimeDiffvsy1 = new TH2D("hCDetTimeDiffvsy1", "CDet Time Diff vs y1 position; y1 pos (m);CDet dt (ns)", 600,-0.5,0.5,TDCBinNum,DiffMin,DiffMax);
  TH2D* hCDetTimeDiffvsy2 = new TH2D("hCDetTimeDiffvsy2", "CDet Time Diff vs y2 position; y2 pos (m);CDet dt (ns)", 600,-0.5,0.5,TDCBinNum,DiffMin,DiffMax);

  TH1D* hDtCDetECal = new TH1D("hDtCDetECal", "CDet t - ECal t;CDet t - ECal t (ns);Counts", TDCBinNum, tdiffECalCDetMin, tdiffECalCDetMax);
  TH2D* hDtvsDxCDetECal = new TH2D("hDtvsDxCDetECal", "CDet ECal dt vs dx;dx_ECalCDet (m);dt_ECalCDet (ns)", NXDiffBins, XDiffLow, XDiffHigh, TDCBinNum, tdiffECalCDetMin, tdiffECalCDetMax);
  TH2D* hSelectedBarX2VsX1 = new TH2D("hSelectedBarX2VsX1", TString::Format("Matched CDet x positions for Layer 1 bar %d (pixels %d-%d);Layer 1 x (m);Matched Layer 2 x (m)", selectedBarNumber, selectedLayer1BarBase, selectedLayer1BarBase + NumPaddles - 1), NXBins, XMin, XMax, NXBins, XMin, XMax);
  TH2D* hSelectedBarECalXVsX1 = new TH2D("hSelectedBarECalXVsX1", TString::Format("ECal projection at Layer 1 vs CDet x, Layer 1 bar %d;Layer 1 x (m);ECal x projected to Layer 1 (m)", selectedBarNumber), NXBins, XMin, XMax, NXBins, XMin, XMax);
  TH2D* hSelectedBarECalXVsX2 = new TH2D("hSelectedBarECalXVsX2", TString::Format("ECal projection at Layer 2 vs matched CDet x, Layer 1 bar %d;Matched Layer 2 x (m);ECal x projected to Layer 2 (m)", selectedBarNumber), NXBins, XMin, XMax, NXBins, XMin, XMax);
  TH2D* hSelectedBarXByPlane = new TH2D("hSelectedBarXByPlane", TString::Format("Hit x through CDet and ECal, Layer 1 bar %d;Detector plane;x position (m)", selectedBarNumber), 3, 0.5, 3.5, NXBins, XMin, XMax);
  TH2D* hSelectedBarXVsZ = new TH2D("hSelectedBarXVsZ", TString::Format("Hit positions through CDet and ECal, Layer 1 bar %d;z position (m);x position (m)", selectedBarNumber), NZBins, ZMin, ZMax, NXBins, XMin, XMax);
  TH1D* hCDetAllDtCutPixelOffsetLe = new TH1D(
      "hCDetAllDtCutPixelOffsetLe",
      "All CDet hits passing ECal-CDet #Deltat cut;Pixel-offset-only LE (ns);Hits",
      NPairedLeBins, CDetMin, CDetMax);
  std::vector<TH1D*> hSelectedBarDtCutPixelOffsetLe(NumPaddles, nullptr);
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    const int globalPixel = selectedLayer1BarBase + localPixel;
    hSelectedBarDtCutPixelOffsetLe[localPixel] = new TH1D(
        TString::Format("hSelectedBarDtCutPixelOffsetLe_Pixel%d", globalPixel),
        TString::Format("ECal-CDet #Deltat-selected LE, Pixel %d;Pixel-offset-only LE (ns);Hits", globalPixel),
        NPairedLeBins, CDetMin, CDetMax);
  }
  std::vector<TH1D*> hSelectedBarPairedLe(NumPaddles, nullptr);
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    const int globalPixel = selectedLayer1BarBase + localPixel;
    hSelectedBarPairedLe[localPixel] = new TH1D(TString::Format("hSelectedBarPairedLe_Pixel%d", globalPixel), TString::Format("Paired Good LE (Pixel %d);LE (ns);Counts", globalPixel), NPairedLeBins, CDetMin, CDetMax);
  }
  gCDetPairedLeSpectra.assign(NumCDetPaddles, nullptr);
  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
    gCDetPairedLeSpectra[pixel] = new TH1D(
        TString::Format("hPairedGoodLe_Pixel%04d", pixel),
        TString::Format("Paired Good LE (Pixel %d);LE (ns);Counts", pixel),
        NPairedLeBins, CDetMin, CDetMax);
    gCDetPairedLeSpectra[pixel]->SetDirectory(nullptr);
  }
  const int barsPerLayer = NumCDetPaddles/(NumLayers*NumPaddles);
  gCDetBarECalTimingSpectra.assign(NumLayers*barsPerLayer, nullptr);
  // Preserve unbinned accepted-hit coordinates for the fixed-effects timing
  // regression.  The TH2D bank below is intentionally visualization-only:
  // its coarse ECal bins attenuate changes in the applied slope.
  gCDetHalfBarECalTimingSamples.assign(NumLayers*barsPerLayer, {});
  gCDetFrozenECalTimingSamples.clear();
  std::vector<std::vector<std::pair<double,double>>> halfBarHCalTimingSamples(
      NumLayers*barsPerLayer);
  gCDetHalfBarECalYSamples.assign(NumLayers*barsPerLayer, {});
  gCDetYPairedSamples.clear();
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int halfBar = 0; halfBar < barsPerLayer; ++halfBar) {
      const int index = layer*barsPerLayer + halfBar;
      gCDetBarECalTimingSpectra[index] = new TH2D(
          TString::Format("hCDetL%dHalfBar%03dTvsECalT", layer + 1, halfBar),
          TString::Format("CDet Layer %d half-bar %d t vs ECal Time;ECal ADC Time (ns);CDet t (ns)",
                          layer + 1, halfBar),
          NADCBins, ECalMin, ECalMax, TDCBinNum, CDetMin, CDetMax);
      gCDetBarECalTimingSpectra[index]->SetDirectory(nullptr);
    }
  }
  hSelectedBarXByPlane->GetXaxis()->SetBinLabel(1, "CDet L1");
  hSelectedBarXByPlane->GetXaxis()->SetBinLabel(2, "CDet L2");
  hSelectedBarXByPlane->GetXaxis()->SetBinLabel(3, "ECal");

  const size_t Nev = vGoodLe.size();
  pairs_CDet.clear();
  pairs_CDet.resize(Nev);

  // Determine each half-bar's physical footprint from the detector coordinates
  // already stored with the accepted hits.  A half-bar contains 16 adjacent
  // pixels in x and is CDet_y_half_length long on either side of its center.
  // The small x padding extends the observed outer pixel centers to the actual
  // edges of the half-bar.
  struct HalfBarFootprint {
    double xMin = std::numeric_limits<double>::infinity();
    double xMax = -std::numeric_limits<double>::infinity();
    double ySum = 0.0;
    double zSum = 0.0;
    long long samples = 0;
  };
  const int totalHalfBars = NumCDetPaddles/NumPaddles;
  std::vector<std::vector<std::pair<double,double>>> yCorrectionBeforeSamples(totalHalfBars);
  std::vector<std::vector<std::pair<double,double>>> yCorrectionAfterSamples(totalHalfBars);
  std::vector<HalfBarFootprint> halfBarFootprints(totalHalfBars);
  for (size_t ev = 0; ev < vGoodID.size(); ++ev) {
    const size_t nCoordinates = std::min(
        vGoodID[ev].size(),
        std::min(vCDetGoodX[ev].size(), std::min(vCDetGoodY[ev].size(), vCDetGoodZ[ev].size())));
    for (size_t hit = 0; hit < nCoordinates; ++hit) {
      const int pixel = vGoodID[ev][hit];
      if (pixel < 0 || pixel >= NumCDetPaddles) continue;
      HalfBarFootprint &footprint = halfBarFootprints[pixel/NumPaddles];
      footprint.xMin = std::min(footprint.xMin, vCDetGoodX[ev][hit]);
      footprint.xMax = std::max(footprint.xMax, vCDetGoodX[ev][hit]);
      footprint.ySum += vCDetGoodY[ev][hit];
      footprint.zSum += vCDetGoodZ[ev][hit];
      ++footprint.samples;
    }
  }
  const double halfPixelPitch = 0.5*0.00525*0.5*(XCorr1 + XCorr2);
  auto projectionIsInHitHalfBar = [&](int pixel, double ecalX, double ecalY,
                                      double hitZ) {
    if (pixel < 0 || pixel >= NumCDetPaddles) return false;
    const HalfBarFootprint &footprint = halfBarFootprints[pixel/NumPaddles];
    if (footprint.samples <= 0 || !std::isfinite(footprint.xMin) ||
        !std::isfinite(footprint.xMax)) return false;
    const double xProjection = ecalX*hitZ/ECal_dist;
    const double yProjection = ecalY*hitZ/ECal_dist;
    const double yCenter = footprint.ySum/footprint.samples;
    return xProjection >= footprint.xMin - halfPixelPitch &&
           xProjection <= footprint.xMax + halfPixelPitch &&
           std::fabs(yProjection - yCenter) <= CDet_y_half_length;
  };
  auto applyYPropagationCorrection = [&](int pixel, double time, double ecalY,
                                         double hitZ) {
    if (!applyYPropagationDiagnostic || pixel < 0 || pixel >= NumCDetPaddles)
      return time;
    const HalfBarFootprint &footprint = halfBarFootprints[pixel/NumPaddles];
    if (footprint.samples <= 0) return time;
    const double yProjection = ecalY*hitZ/ECal_dist;
    const double yCenter = footprint.ySum/footprint.samples;
    const int pixelsPerLayer = NumCDetPaddles/NumLayers;
    const int withinLayer = pixel%pixelsPerLayer;
    const int side = withinLayer/NumCDetPaddlesPerSide;
    const double signedSlope = side == 0
        ? -yPropagationSlopeMagnitude : yPropagationSlopeMagnitude;
    return time - signedSlope*(yProjection - yCenter);
  };
  long long projectedHalfBarPairs = 0;

  // ------------------------------
  // Tuning parameters
  // ------------------------------
  double dt0   = 0.0;   // ns; start 0, later set to peak of dt histogram
  double dtWin = 0.5 * (diffMaxCut-diffMinCut);  // ns; start wide (10–15)
  double dxWin = 0.5 * (xdiffMaxCut-xdiffMinCut);   // m; ~one bar

  double sigT  = 1.0;   // ns; timing scale for score
  double sigX  = 0.01;   // m; start ~1, later tighten toward 0.5

  // Recover the LE value that existed immediately after the per-pixel offset
  // was applied. Pair finding and the ECal-CDet dt cut continue to use the
  // fully corrected time; this inverse is used only by the two new displays.
  auto pixelOffsetOnlyLe = [&](double correctedTime, int pixel, double tot,
                               double ecalTime) {
    double time = correctedTime;
    if (std::isfinite(gGlobalTimingShift)) time -= gGlobalTimingShift;
    if (gUseTimeWalkCorr)
      time += GetTimeWalkCorrection(GetLayerFromID(pixel), tot);
    if (gUseECalTimeCorr)
      time += gECalFitP0 + gECalFitP1*ecalTime - gECalDeltaShift;
    return time;
  };

  for (size_t ev = 0; ev < Nev; ev++) { //iterate through events
    std::vector<double> vCDet1Time;
    std::vector<double> vCDet2Time;
    std::vector<double> vCDet1Tot;
    std::vector<double> vCDet2Tot;
    std::vector<double> vCDet1x;
    std::vector<double> vCDet2x;
    std::vector<double> vCDet1y;
    std::vector<double> vCDet2y;
    std::vector<double> vCDet1z;
    std::vector<double> vCDet2z;
    std::vector<double> vCDet1ID;
    std::vector<double> vCDet2ID;
    double t_ECal = v_GoodECalAdcTime[ev];
    const double t_HCal = ev < v_GoodHCalAdcTime.size()
        ? v_GoodHCalAdcTime[ev] : std::numeric_limits<double>::quiet_NaN();

    // ------------------------------
    // First pass: split hits into layers
    // ------------------------------
    const size_t Nhits = std::min(vGoodLe[ev].size(), vGoodTot[ev].size());
    for (size_t ihit = 0; ihit < Nhits; ++ihit) {
      const double timingForCuts = vGoodLe[ev][ihit] -
          (gShiftInvariantPairTimingCuts ? gGlobalTimingShift : 0.0);
      if (timingForCuts >= LeMin && timingForCuts <= LeMax && vGoodTot[ev][ihit] >= TotMinCut && vGoodTot[ev][ihit] <= TotMaxCut && t_ECal >= ECalMin && t_ECal <= ECalMax){
        if (vGoodID[ev][ihit] >= 0 && vGoodID[ev][ihit] <= 1343){ //layer 1 hits
          vCDet1Time.push_back(vGoodLe[ev][ihit]);
          vCDet1Tot.push_back(vGoodTot[ev][ihit]);
          vCDet1ID.push_back(vGoodID[ev][ihit]);
          vCDet1x.push_back(vCDetGoodX[ev][ihit]);
          vCDet1y.push_back(vCDetGoodY[ev][ihit]);
          vCDet1z.push_back(vCDetGoodZ[ev][ihit]);
        }
        else if (vGoodID[ev][ihit] >= 1344 && vGoodID[ev][ihit] <= 2687){//layer 2 hits
          vCDet2Time.push_back(vGoodLe[ev][ihit]);
          vCDet2Tot.push_back(vGoodTot[ev][ihit]);
          vCDet2ID.push_back(vGoodID[ev][ihit]);
          vCDet2x.push_back(vCDetGoodX[ev][ihit]);
          vCDet2y.push_back(vCDetGoodY[ev][ihit]);
          vCDet2z.push_back(vCDetGoodZ[ev][ihit]);
        }
      }
    }

    //require at least one hit in each layer to even try pairing:
    if (vCDet1Time.empty() || vCDet2Time.empty()) continue;

    // ------------------------------
    // Second pass: build candidate pairs (all pairs that pass dt/dx cuts)
    // ------------------------------
    std::vector<Cand> cands;
    cands.reserve(vCDet1Time.size() * vCDet2Time.size());

    for (int i1 = 0; i1 < (int)vCDet1Time.size(); i1++) {
      for (int j2 = 0; j2 < (int)vCDet2Time.size(); j2++) {

        const double dt = vCDet2Time[j2] - vCDet1Time[i1];
        const double dx = vCDet2x[j2] - vCDet1x[i1];

        if (dt < diffMinCut || dt > diffMaxCut) continue;
        if (fabs(dx) > dxWin) continue;
        if (fabs(vCDet2y[j2] - vCDet1y[i1]) > 0.08) continue; //require hits to be on same side of CDet (each y-cord is ~7.5cm apart)

        const double score = (((dt - dt0)*(dt - dt0)) / (sigT*sigT)) + ((dx*dx) / (sigX*sigX));

        cands.push_back({i1, j2, dt, dx, score});
      }
    }

    if (cands.empty()) continue;

    // ------------------------------
    // Sort candidates by score (best first)
    // ------------------------------
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) {
                return a.score < b.score;
              });

    // ------------------------------
    // Select non-conflicting pairs (greedy one-to-one matching)
    //   - ensures you don’t reuse the same hit in L1 or L2 twice
    // ------------------------------
    std::vector<char> used1(vCDet1Time.size(), 0);
    std::vector<char> used2(vCDet2Time.size(), 0);

    pairs_CDet[ev].clear();
    pairs_CDet[ev].reserve(std::min(vCDet1Time.size(), vCDet2Time.size()));

    for (const auto &c : cands) {
      if (used1[c.i1] || used2[c.j2]) continue;

      used1[c.i1] = 1;
      used2[c.j2] = 1;
      PairHit p;
      p.t1 = vCDet1Time[c.i1];
      p.tot1 = vCDet1Tot[c.i1];
      p.x1 = vCDet1x[c.i1];
      p.y1 = vCDet1y[c.i1];
      p.z1 = vCDet1z[c.i1];
      p.id1 = vCDet1ID[c.i1];

      p.t2 = vCDet2Time[c.j2];
      p.tot2 = vCDet2Tot[c.j2];
      p.x2 = vCDet2x[c.j2];
      p.y2 = vCDet2y[c.j2];
      p.z2 = vCDet2z[c.j2];
      p.id2 = vCDet2ID[c.j2];

      p.dt = c.dt;
      p.dx = c.dx;
      p.score = c.score;

      pairs_CDet[ev].push_back(p);

      if (!allowMultiplePairs) break;
    }
    hNpairPerEvent->Fill((int)pairs_CDet[ev].size());

    // ------------------------------
    // Third pass: fill histograms using the ACCEPTED pairs
    // ------------------------------
    //if (ev <= 20) {std::cout << "CDet1 Size= " << vCDet1Time.size() << " CDet2 Size= " << vCDet2Time.size() << " Npairs= " << pairs_CDet[ev].size() << "\n";
   // }

    for (size_t ip = 0; ip < pairs_CDet[ev].size(); ip++) {
      const auto &p = pairs_CDet[ev][ip];

      /*if (ev <= 20) {
        if (ip == 0) std::cout << "--------new event ----------\n";
        std::cout << "pair=" << ip
                  << " t1=" << p.t1 << " x1=" << p.x1 << " TOT1= " << p.tot1 << " ID1=" << p.id1 << "\n";
        std::cout << "pair=" << ip
                  << " t2=" << p.t2 << " x2=" << p.x2 << " TOT2= " << p.tot2 << " ID2=" << p.id2 << "\n";
        std::cout << "pair=" << ip
                  << " tdiff=" << p.dt << " xdiff=" << p.dx
                  << " score=" << p.score << "\n";
        std::cout << "===============\n";
      }*/

      // Apply final/tighter cuts (these can be narrower than dtWin/dxWin)
      if (p.dt >= diffMinCut && p.dt <= diffMaxCut && p.dx >= xdiffMinCut && p.dx <= xdiffMaxCut) {
        double t_pair = (p.t1 + p.t2) / 2;
        const double t1YCorrected = applyYPropagationCorrection(
            p.id1, p.t1, v_GoodECalY[ev], p.z1);
        const double t2YCorrected = applyYPropagationCorrection(
            p.id2, p.t2, v_GoodECalY[ev], p.z2);
        const double tPairYCorrected = 0.5*(t1YCorrected + t2YCorrected);
        const double dtYCorrected = t2YCorrected - t1YCorrected;
        const double timingShiftForCuts =
            gShiftInvariantPairTimingCuts ? gGlobalTimingShift : 0.0;
        double dt_EC = t_ECal - (t_pair - timingShiftForCuts);
        if (dt_EC >= tdiffECalCDetMin && dt_EC <= tdiffECalCDetMax){
          const double pixelOffsetLe1 = pixelOffsetOnlyLe(
              p.t1, p.id1, p.tot1, t_ECal);
          const double pixelOffsetLe2 = pixelOffsetOnlyLe(
              p.t2, p.id2, p.tot2, t_ECal);
          hCDetAllDtCutPixelOffsetLe->Fill(pixelOffsetLe1);
          hCDetAllDtCutPixelOffsetLe->Fill(pixelOffsetLe2);
          if (p.id1 >= selectedLayer1BarBase &&
              p.id1 < selectedLayer1BarBase + NumPaddles)
            hSelectedBarDtCutPixelOffsetLe[p.id1 - selectedLayer1BarBase]->Fill(
                pixelOffsetLe1);

          double x_pair = (p.x1 + p.x2) / 2;
          double z_pair = (p.z1 + p.z2) / 2;
          double x_ECal = v_GoodECalX[ev]*z_pair/ECal_dist;
          double dx_EC = x_pair - x_ECal;

          hCDet1IDvs2ID->Fill(p.id2, p.id1);
          hCDetTimeDiff->Fill(p.dt);
          hCDetXDiff->Fill(p.dx);
          hCDetLe2vs1->Fill(p.t1, p.t2);
          hCDetTot2vs1->Fill(p.tot1, p.tot2);
          hCDetLe1->Fill(p.t1);
          hCDetLe2->Fill(p.t2);
          if (applyYPropagationDiagnostic) {
            hCDetTimeDiffYCorrected->Fill(dtYCorrected);
            hCDetLe1YCorrected->Fill(t1YCorrected);
            hCDetLe2YCorrected->Fill(t2YCorrected);
            hCDetLe2vs1YCorrected->Fill(t1YCorrected, t2YCorrected);
          }
          h2CDetx2VsCDetx1->Fill(p.x1,p.x2);
          hCDet1LeVsTot->Fill(p.tot1,p.t1);
          hCDet2LeVsTot->Fill(p.tot2,p.t2);
          auto fillTimingStructure = [&](double time, double tot) {
            if (t_ECal > -5.0 && t_ECal < 1.0) {
              if (time > 23.0 && time < 30.0)
                hCDetAllLowerVerticalLeVsTot->Fill(tot, time);
              else if (time > 30.0 && time < 40.0)
                hCDetAllUpperVerticalLeVsTot->Fill(tot, time);
            }
            if (t_ECal > 1.0 && time > 30.0 &&
                std::fabs(time - (27.0 + 0.73*t_ECal)) < 3.0)
              hCDetAllDiagonalLeVsTot->Fill(tot, time);
          };
          fillTimingStructure(p.t1, p.tot1);
          fillTimingStructure(p.t2, p.tot2);
          hECalVsCDetDt->Fill(t_ECal,p.dt);
          hCDetTimeDiffvsx1->Fill(p.x1, p.dt);
          hCDetTimeDiffvsx2->Fill(p.x2, p.dt);
          hCDetTimeDiffvsy1->Fill(p.y1, p.dt);
          hCDetTimeDiffvsy2->Fill(p.y2, p.dt);

          hECalVsCDetT->Fill(t_ECal,t_pair);
          if (gUseECalTimeCorr) {
            const double p1Term = gECalFitP1*t_ECal;
            hECalVsCDetTBeforeP1->Fill(t_ECal, t_pair + p1Term);
            hECalVsCDetTL1BeforeP1->Fill(t_ECal, p.t1 + p1Term);
            hECalVsCDetTL2BeforeP1->Fill(t_ECal, p.t2 + p1Term);
          }
          if (std::isfinite(t_HCal) && t_HCal >= HCalMin && t_HCal <= HCalMax) {
            hHCalVsCDetT->Fill(t_HCal, t_pair);
            hECalVsHCalAccepted->Fill(t_ECal, t_HCal);
          }
          gCDetAcceptedPairMeanTimes.push_back(t_pair);
          hECalVsCDetTL1->Fill(t_ECal,p.t1);
          hECalVsCDetTL2->Fill(t_ECal,p.t2);
          const bool layer1ProjectionMatch = projectionIsInHitHalfBar(
              p.id1, v_GoodECalX[ev], v_GoodECalY[ev], p.z1);
          const bool layer2ProjectionMatch = projectionIsInHitHalfBar(
              p.id2, v_GoodECalX[ev], v_GoodECalY[ev], p.z2);
          if (layer1ProjectionMatch && layer2ProjectionMatch) {
            hCDetProjectedHalfBarLe1->Fill(p.t1);
            hCDetProjectedHalfBarLe2->Fill(p.t2);
            hCDetProjectedHalfBarPairMean->Fill(t_pair);
            gCDetProjectedHalfBarPairMeanTimes.push_back(t_pair);
            hCDetProjectedHalfBarLe2vs1->Fill(p.t1, p.t2);
            if (applyYPropagationDiagnostic) {
              hCDetProjectedHalfBarLe1YCorrected->Fill(t1YCorrected);
              hCDetProjectedHalfBarLe2YCorrected->Fill(t2YCorrected);
              hCDetProjectedHalfBarPairMeanYCorrected->Fill(tPairYCorrected);
              hCDetProjectedHalfBarLe2vs1YCorrected->Fill(t1YCorrected, t2YCorrected);
            }
            ++projectedHalfBarPairs;
          }
          auto fillCDetTimeVsECalY = [&](int pixel, double time) {
            if (pixel < 0 || pixel >= NumCDetPaddles) return;
            const int pixelsPerLayer = NumCDetPaddles/NumLayers;
            const int layer = pixel/pixelsPerLayer;
            const int withinLayer = pixel%pixelsPerLayer;
            const int side = withinLayer/NumCDetPaddlesPerSide;
            const int withinSide = withinLayer%NumCDetPaddlesPerSide;
            const int module = withinSide/(NumBars*NumPaddles);
            const int pmt = (withinSide%(NumBars*NumPaddles))/NumPaddles;
            const int section = 2*module + pmt/NumHalfBarsPerBank;
            const int sectionGroup = std::min(section, 5 - section);
            hCDetTvsECalY[layer][side][sectionGroup]->Fill(v_GoodECalY[ev], time);
            if (time >= 20.0 && time <= 40.0)
              gCDetHalfBarECalYSamples[pixel/NumPaddles].emplace_back(v_GoodECalY[ev], time);
          };
          fillCDetTimeVsECalY(p.id1, p.t1);
          fillCDetTimeVsECalY(p.id2, p.t2);
          if (applyYPropagationDiagnostic) {
            const int halfBar1 = p.id1/NumPaddles;
            const int halfBar2 = p.id2/NumPaddles;
            const double yProjection1 = v_GoodECalY[ev]*p.z1/ECal_dist;
            const double yProjection2 = v_GoodECalY[ev]*p.z2/ECal_dist;
            yCorrectionBeforeSamples[halfBar1].emplace_back(yProjection1, p.t1);
            yCorrectionAfterSamples[halfBar1].emplace_back(yProjection1, t1YCorrected);
            yCorrectionBeforeSamples[halfBar2].emplace_back(yProjection2, p.t2);
            yCorrectionAfterSamples[halfBar2].emplace_back(yProjection2, t2YCorrected);
          }
          if (p.t1 >= 20.0 && p.t1 <= 40.0 && p.t2 >= 20.0 && p.t2 <= 40.0)
            gCDetYPairedSamples.push_back({p.id1/NumPaddles, p.id2/NumPaddles, p.t1, p.t2});
          if (p.id1 >= 0 && p.id1 < NumCDetPaddles)
            gCDetPairedLeSpectra[(int)p.id1]->Fill(p.t1);
          if (p.id2 >= 0 && p.id2 < NumCDetPaddles)
            gCDetPairedLeSpectra[(int)p.id2]->Fill(p.t2);
          auto fillHalfBarTiming = [&](int pixel, double time, double tot,
                                       int layer) {
            const int pixelsPerLayer = NumCDetPaddles/NumLayers;
            const int halfBar = (pixel%pixelsPerLayer)/NumPaddles;
            const int index = layer*barsPerLayer + halfBar;
            if (index >= 0 && index < (int)gCDetBarECalTimingSpectra.size()) {
              gCDetBarECalTimingSpectra[index]->Fill(t_ECal, time);
              gCDetHalfBarECalTimingSamples[index].push_back(std::make_pair(t_ECal, time));
              gCDetFrozenECalTimingSamples.push_back(
                  {ev, ip, layer, pixel, index, t_ECal, time, tot});
              if (std::isfinite(t_HCal) && t_HCal >= HCalMin && t_HCal <= HCalMax)
                halfBarHCalTimingSamples[index].push_back(std::make_pair(t_HCal, time));
            }
          };
          fillHalfBarTiming(p.id1, p.t1, p.tot1, 0);
          fillHalfBarTiming(p.id2, p.t2, p.tot2, 1);

          hDtCDetECal->Fill(dt_EC);
          hDtvsDxCDetECal->Fill(dx_EC, dt_EC);

          const bool isSelectedLayer1Bar = p.id1 >= selectedLayer1BarBase && p.id1 < selectedLayer1BarBase + NumPaddles;
          const bool isSelectedLayer2Bar = p.id2 >= selectedLayer2BarBase && p.id2 < selectedLayer2BarBase + NumPaddles;
          if (isSelectedLayer1Bar) {
            hECalVsSelectedBarT->Fill(t_ECal,p.t1);
            if (gUseECalTimeCorr)
              hECalVsSelectedBarTBeforeP1->Fill(
                  t_ECal, p.t1 + gECalFitP1*t_ECal);
            hCDetBarLe1->Fill(p.t1);
            hCDet1BarLeVsTot->Fill(p.tot1, p.t1);
            if (t_ECal > -3.0 && t_ECal < 1.0 && p.t1 > 23.0 && p.t1 < 40.0)
              hCDet1BarVerticalLeVsTot->Fill(p.tot1, p.t1);
            const double diagonalResidual = p.t1 - (27.0 + 0.73*t_ECal);
            if (t_ECal > 1.0 && t_ECal < 18.0 && std::fabs(diagonalResidual) < 3.0)
              hCDet1BarDiagonalLeVsTot->Fill(p.tot1, p.t1);
            hSelectedBarPairedLe[(int)p.id1 - selectedLayer1BarBase]->Fill(p.t1);
            const double xECalAtLayer1 = v_GoodECalX[ev]*p.z1/ECal_dist;
            const double xECalAtLayer2 = v_GoodECalX[ev]*p.z2/ECal_dist;
            hSelectedBarX2VsX1->Fill(p.x1, p.x2);
            hSelectedBarECalXVsX1->Fill(p.x1, xECalAtLayer1);
            hSelectedBarECalXVsX2->Fill(p.x2, xECalAtLayer2);
            hSelectedBarXByPlane->Fill(1.0, p.x1);
            hSelectedBarXByPlane->Fill(2.0, p.x2);
            hSelectedBarXByPlane->Fill(3.0, v_GoodECalX[ev]);
            hSelectedBarXVsZ->Fill(p.z1, p.x1);
            hSelectedBarXVsZ->Fill(p.z2, p.x2);
            hSelectedBarXVsZ->Fill(ECal_dist, v_GoodECalX[ev]);
          }
          if (isSelectedLayer2Bar) {
            hCDetBarLe2->Fill(p.t2);
            hCDet2BarLeVsTot->Fill(p.tot2, p.t2);
          }
          if (isSelectedLayer1Bar && isSelectedLayer2Bar) {
            hECalVsCDetDtSingle->Fill(t_ECal,p.dt);
            hECalVsCDetTSingle->Fill(t_ECal,p.t1);
          }
        }//ecal cdet tdiff cut
	    }//fill histogram with cuts
    }// end pair hits loop

  }// end event loop
  // ------------------------------
  // Make Plots
  // ------------------------------
  std::vector<TF1*> fSelectedBarPairedLe(NumPaddles, nullptr);
  std::vector<bool> validSelectedBarPairedLeFit(NumPaddles, false);
  TCanvas *cSelectedBarPairedLe = new TCanvas("cSelectedBarPairedLe", TString::Format("Paired-hit LE spectra for Layer 1 bar %d", selectedBarNumber), 1200,1000);
  cSelectedBarPairedLe->Divide(4,4,0.001,0.001);
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    const int globalPixel = selectedLayer1BarBase + localPixel;
    const double displayFitMin = gTargetMeanLE - 5.0;
    const double displayFitMax = gTargetMeanLE + 5.0;
    fSelectedBarPairedLe[localPixel] = new TF1(TString::Format("fSelectedBarPairedLe_Pixel%d", globalPixel), "gaus", displayFitMin, displayFitMax);
    const int firstFitBin = hSelectedBarPairedLe[localPixel]->GetXaxis()->FindBin(displayFitMin);
    const int lastFitBin = hSelectedBarPairedLe[localPixel]->GetXaxis()->FindBin(displayFitMax);
    const double fitEntries = hSelectedBarPairedLe[localPixel]->Integral(firstFitBin, lastFitBin);
    if (fitEntries > 20) {
      int peakBin = firstFitBin;
      for (int bin = firstFitBin + 1; bin <= lastFitBin; ++bin) {
        if (hSelectedBarPairedLe[localPixel]->GetBinContent(bin) > hSelectedBarPairedLe[localPixel]->GetBinContent(peakBin)) peakBin = bin;
      }
      fSelectedBarPairedLe[localPixel]->SetParameters(hSelectedBarPairedLe[localPixel]->GetBinContent(peakBin), hSelectedBarPairedLe[localPixel]->GetBinCenter(peakBin), 2.0);
      validSelectedBarPairedLeFit[localPixel] =
          hSelectedBarPairedLe[localPixel]->Fit(fSelectedBarPairedLe[localPixel], "RQ0") == 0;
    }

    cSelectedBarPairedLe->cd(localPixel + 1);
    hSelectedBarPairedLe[localPixel]->Draw();
    if (validSelectedBarPairedLeFit[localPixel]) fSelectedBarPairedLe[localPixel]->Draw("SAME");
    if (kUnusedCDetPixels.count(globalPixel)) {
      TPaveText *flag = new TPaveText(0.82,0.82,0.95,0.95,"NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetTextColor(kWhite);
      flag->SetBorderSize(1);
      flag->AddText("UNUSED PIXEL");
      flag->Draw("SAME");
    }
  }
  cSelectedBarPairedLe->Update();

  TCanvas *cCDetAllDtCutPixelOffsetLe = new TCanvas(
      "cCDetAllDtCutPixelOffsetLe",
      "Pixel-offset-only CDet LE after ECal-CDet dt selection", 1000, 700);
  hCDetAllDtCutPixelOffsetLe->Draw();
  cCDetAllDtCutPixelOffsetLe->Update();

  TCanvas *cSelectedBarDtCutPixelOffsetLe = new TCanvas(
      "cSelectedBarDtCutPixelOffsetLe",
      TString::Format("Pixel-offset-only LE for Layer 1 bar %d after ECal-CDet dt selection", selectedBarNumber),
      1200, 1000);
  cSelectedBarDtCutPixelOffsetLe->Divide(4, 4, 0.001, 0.001);
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    cSelectedBarDtCutPixelOffsetLe->cd(localPixel + 1);
    hSelectedBarDtCutPixelOffsetLe[localPixel]->Draw();
  }
  cSelectedBarDtCutPixelOffsetLe->Update();

  TCanvas *cCDetLayerTimes = new TCanvas("cCDetLayerTimes", "CDet Layer 1 and 2 LE",900,700);
  cCDetLayerTimes->Divide(1,2);

  cCDetLayerTimes->cd(1);
  //gPad->SetLogz();
  hCDetLe1->Draw();

  cCDetLayerTimes->cd(2);
  //gPad->SetLogz();
  hCDetLe2->Draw();

  if (applyYPropagationDiagnostic) {
    TCanvas *cCDetLayerTimesYCorrected = new TCanvas(
        "cCDetLayerTimesYCorrected",
        "CDet Layer 1 and 2 LE after y-propagation correction", 900, 700);
    cCDetLayerTimesYCorrected->Divide(1,2);
    cCDetLayerTimesYCorrected->cd(1);
    hCDetLe1YCorrected->Draw();
    cCDetLayerTimesYCorrected->cd(2);
    hCDetLe2YCorrected->Draw();
    cCDetLayerTimesYCorrected->Update();
    cCDetLayerTimesYCorrected->SaveAs("cCDetLayerTimesYCorrected.jpg");
  }

  gCDetProjectedHalfBarTimingCanvas = new TCanvas(
      "cCDetProjectedHalfBarTiming",
      "Corrected timing for pairs in the ECal-projected CDet half-bars", 1200, 900);
  gCDetProjectedHalfBarTimingCanvas->Divide(2,2);
  gCDetProjectedHalfBarTimingCanvas->cd(1);
  hCDetProjectedHalfBarLe1->Draw();
  gCDetProjectedHalfBarTimingCanvas->cd(2);
  hCDetProjectedHalfBarLe2->Draw();
  gCDetProjectedHalfBarTimingCanvas->cd(3);
  hCDetProjectedHalfBarPairMean->Draw();
  gCDetProjectedHalfBarTimingCanvas->cd(4);
  hCDetProjectedHalfBarLe2vs1->Draw("COLZ");
  gCDetProjectedHalfBarTimingCanvas->Update();

  std::cout << "\n[CDet projected-half-bar timing]\n"
            << "  accepted layer pairs in the nominal-origin ECal-projected half-bars: "
            << projectedHalfBarPairs << "\n"
            << "  Layer 1 corrected LE: mean=" << hCDetProjectedHalfBarLe1->GetMean()
            << " ns, RMS=" << hCDetProjectedHalfBarLe1->GetStdDev() << " ns\n"
            << "  Layer 2 corrected LE: mean=" << hCDetProjectedHalfBarLe2->GetMean()
            << " ns, RMS=" << hCDetProjectedHalfBarLe2->GetStdDev() << " ns\n"
            << "  pair mean corrected LE: mean=" << hCDetProjectedHalfBarPairMean->GetMean()
            << " ns, RMS=" << hCDetProjectedHalfBarPairMean->GetStdDev() << " ns\n";
  if (applyYPropagationDiagnostic) {
    TCanvas *cCDetProjectedHalfBarTimingYCorrected = new TCanvas(
        "cCDetProjectedHalfBarTimingYCorrected",
        "Projected-half-bar timing after y-propagation correction", 1200, 900);
    cCDetProjectedHalfBarTimingYCorrected->Divide(2,2);
    cCDetProjectedHalfBarTimingYCorrected->cd(1);
    hCDetProjectedHalfBarLe1YCorrected->Draw();
    cCDetProjectedHalfBarTimingYCorrected->cd(2);
    hCDetProjectedHalfBarLe2YCorrected->Draw();
    cCDetProjectedHalfBarTimingYCorrected->cd(3);
    hCDetProjectedHalfBarPairMeanYCorrected->Draw();
    cCDetProjectedHalfBarTimingYCorrected->cd(4);
    hCDetProjectedHalfBarLe2vs1YCorrected->Draw("COLZ");
    cCDetProjectedHalfBarTimingYCorrected->Update();
    cCDetProjectedHalfBarTimingYCorrected->SaveAs(
        "cCDetProjectedHalfBarTimingYCorrected.jpg");
    std::cout << "[CDet projected-half-bar timing after y correction]\n"
              << "  Layer 1: mean=" << hCDetProjectedHalfBarLe1YCorrected->GetMean()
              << " ns, RMS=" << hCDetProjectedHalfBarLe1YCorrected->GetStdDev() << " ns\n"
              << "  Layer 2: mean=" << hCDetProjectedHalfBarLe2YCorrected->GetMean()
              << " ns, RMS=" << hCDetProjectedHalfBarLe2YCorrected->GetStdDev() << " ns\n"
              << "  pair mean: mean=" << hCDetProjectedHalfBarPairMeanYCorrected->GetMean()
              << " ns, RMS=" << hCDetProjectedHalfBarPairMeanYCorrected->GetStdDev() << " ns\n";
  }

  TCanvas *cCDetLayerTimeVsECalY = new TCanvas(
      "cCDetLayerTimeVsECalY", "Corrected CDet timing versus ECal y by layer, side, and section", 1800, 1200);
  cCDetLayerTimeVsECalY->Divide(3,4);
  const double positionTimingDisplayMin = 20.0;
  const double positionTimingDisplayMax = 40.0;
  const double positionTimingMinimumBinContent = 5.0;
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int side = 0; side < NumSides; ++side) {
      const int row = layer*NumSides + side;
      for (int group = 0; group < 3; ++group) {
        cCDetLayerTimeVsECalY->cd(row*3 + group + 1);
        TH2D *histogram = hCDetTvsECalY[layer][side][group];
        histogram->GetYaxis()->SetRangeUser(positionTimingDisplayMin, positionTimingDisplayMax);
        histogram->SetMinimum(positionTimingMinimumBinContent);
        histogram->Draw("COLZ");
        TProfile *profile = histogram->ProfileX(
            TString::Format("pCDetL%d%sSectionGroup%dTvsECalY", layer + 1, sideName[side], group + 1));
        profile->SetMarkerStyle(20);
        profile->SetMarkerSize(0.55);
        profile->SetMarkerColor(kBlack);
        profile->SetLineColor(kBlack);
        profile->Draw("SAME");

        int firstFitBin = 0;
        int lastFitBin = 0;
        for (int bin = 1; bin <= profile->GetNbinsX(); ++bin) {
          if (profile->GetBinEntries(bin) < positionTimingMinimumBinContent) continue;
          if (firstFitBin == 0) firstFitBin = bin;
          lastFitBin = bin;
        }
        int fitStatus = -1;
        TF1 *fit = nullptr;
        if (firstFitBin > 0 && lastFitBin > firstFitBin) {
          const double fitMin = profile->GetXaxis()->GetBinLowEdge(firstFitBin);
          const double fitMax = profile->GetXaxis()->GetBinUpEdge(lastFitBin);
          fit = new TF1(
              TString::Format("fCDetL%d%sSectionGroup%dTvsECalY", layer + 1, sideName[side], group + 1),
              "pol1", fitMin, fitMax);
          fit->SetLineColor(kRed);
          fitStatus = profile->Fit(fit, "QRS");
          if (fitStatus == 0) fit->Draw("SAME");
        }

        TPaveText *annotation = new TPaveText(0.12, 0.74, 0.50, 0.90, "NDC");
        annotation->SetFillColor(kWhite);
        annotation->SetTextAlign(12);
        annotation->SetTextSize(0.035);
        if (fitStatus == 0) {
          const double slope = fit->GetParameter(1);
          const double slopeError = fit->GetParError(1);
          const double chi2 = fit->GetChisquare();
          const int ndf = fit->GetNDF();
          annotation->AddText(TString::Format("p1 = %.3f #pm %.3f ns/m", slope, slopeError));
          annotation->AddText(TString::Format("#chi^{2}/NDF = %.1f/%d = %.2f", chi2, ndf,
                                              ndf > 0 ? chi2/ndf : 0.0));
          std::cout << "[CDet y-timing] Layer " << layer + 1 << " " << sideName[side]
                    << ", " << sectionGroupName[group]
                    << ": p1=" << slope << " +/- " << slopeError << " ns/m"
                    << ", chi2/NDF=" << chi2 << "/" << ndf
                    << (ndf > 0 ? TString::Format(" = %.3f", chi2/ndf).Data() : "") << "\n";
        } else {
          annotation->AddText("Linear fit unavailable");
        }
        annotation->Draw("SAME");
      }
    }
  }
  cCDetLayerTimeVsECalY->Update();

  struct CDetYFixedEffectsResult {
    double xy = 0.0;
    double xx = 0.0;
    double slope = NAN;
    double slopeError = NAN;
    double residualSquares = 0.0;
    long long entries = 0;
    int halfBars = 0;
    int ndf = 0;
  };
  CDetYFixedEffectsResult yFixedEffects[NumLayers][NumSides][3];
  TH2D *hCDetYWithin[NumLayers][NumSides][3] = {};
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int side = 0; side < NumSides; ++side) {
      for (int group = 0; group < 3; ++group) {
        hCDetYWithin[layer][side][group] = new TH2D(
            TString::Format("hCDetL%d%sSectionGroup%dECalYWithin", layer + 1, sideName[side], group + 1),
            TString::Format("CDet Layer %d %s, %s fixed effects;ECal y - half-bar mean (m);CDet time - half-bar mean (ns)",
                            layer + 1, sideName[side], sectionGroupName[group]),
            120, -0.35, 0.35, 120, -10.0, 10.0);
      }
    }
  }

  const int minCDetYHalfBarEntries = 100;
  for (int halfBar = 0; halfBar < (int)gCDetHalfBarECalYSamples.size(); ++halfBar) {
    const auto &samples = gCDetHalfBarECalYSamples[halfBar];
    if ((int)samples.size() < minCDetYHalfBarEntries) continue;
    const int pixel = halfBar*NumPaddles;
    const int pixelsPerLayer = NumCDetPaddles/NumLayers;
    const int layer = pixel/pixelsPerLayer;
    const int withinLayer = pixel%pixelsPerLayer;
    const int side = withinLayer/NumCDetPaddlesPerSide;
    const int withinSide = withinLayer%NumCDetPaddlesPerSide;
    const int module = withinSide/(NumBars*NumPaddles);
    const int pmt = (withinSide%(NumBars*NumPaddles))/NumPaddles;
    const int section = 2*module + pmt/NumHalfBarsPerBank;
    const int group = std::min(section, 5 - section);
    double sumY = 0.0, sumT = 0.0;
    for (const auto &sample : samples) { sumY += sample.first; sumT += sample.second; }
    const double meanY = sumY/samples.size();
    const double meanT = sumT/samples.size();
    CDetYFixedEffectsResult &result = yFixedEffects[layer][side][group];
    ++result.halfBars;
    result.entries += samples.size();
    for (const auto &sample : samples) {
      const double dy = sample.first - meanY;
      const double dt = sample.second - meanT;
      result.xy += dy*dt;
      result.xx += dy*dy;
      hCDetYWithin[layer][side][group]->Fill(dy, dt);
    }
  }
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int side = 0; side < NumSides; ++side) {
      for (int group = 0; group < 3; ++group) {
        CDetYFixedEffectsResult &result = yFixedEffects[layer][side][group];
        if (result.xx > 0.0) result.slope = result.xy/result.xx;
      }
    }
  }
  for (int halfBar = 0; halfBar < (int)gCDetHalfBarECalYSamples.size(); ++halfBar) {
    const auto &samples = gCDetHalfBarECalYSamples[halfBar];
    if ((int)samples.size() < minCDetYHalfBarEntries) continue;
    const int pixel = halfBar*NumPaddles;
    const int pixelsPerLayer = NumCDetPaddles/NumLayers;
    const int layer = pixel/pixelsPerLayer;
    const int withinLayer = pixel%pixelsPerLayer;
    const int side = withinLayer/NumCDetPaddlesPerSide;
    const int withinSide = withinLayer%NumCDetPaddlesPerSide;
    const int module = withinSide/(NumBars*NumPaddles);
    const int pmt = (withinSide%(NumBars*NumPaddles))/NumPaddles;
    const int section = 2*module + pmt/NumHalfBarsPerBank;
    const int group = std::min(section, 5 - section);
    double sumY = 0.0, sumT = 0.0;
    for (const auto &sample : samples) { sumY += sample.first; sumT += sample.second; }
    const double meanY = sumY/samples.size();
    const double meanT = sumT/samples.size();
    CDetYFixedEffectsResult &result = yFixedEffects[layer][side][group];
    for (const auto &sample : samples) {
      const double residual = (sample.second - meanT) - result.slope*(sample.first - meanY);
      result.residualSquares += residual*residual;
    }
  }

  TCanvas *cCDetECalYFixedEffects = new TCanvas(
      "cCDetECalYFixedEffects", "CDet versus ECal y fixed-effects slopes", 1800, 1200);
  cCDetECalYFixedEffects->Divide(3,4);
  std::cout << "\n[CDet/ECal y fixed-effects fits]\n"
            << "  Half-bar means are removed before estimating each propagation slope.\n";
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int side = 0; side < NumSides; ++side) {
      const int row = layer*NumSides + side;
      for (int group = 0; group < 3; ++group) {
        cCDetECalYFixedEffects->cd(row*3 + group + 1);
        TH2D *histogram = hCDetYWithin[layer][side][group];
        histogram->Draw("COLZ");
        CDetYFixedEffectsResult &result = yFixedEffects[layer][side][group];
        result.ndf = result.entries - result.halfBars - 1;
        if (result.ndf > 0 && result.xx > 0.0)
          result.slopeError = std::sqrt((result.residualSquares/result.ndf)/result.xx);
        TF1 *fit = new TF1(
            TString::Format("fCDetL%d%sSectionGroup%dECalYWithin", layer + 1, sideName[side], group + 1),
            "[0]*x", -0.30, 0.30);
        fit->SetParameter(0, result.slope);
        fit->SetLineColor(kRed);
        fit->Draw("SAME");
        TPaveText *annotation = new TPaveText(0.12, 0.74, 0.54, 0.90, "NDC");
        annotation->SetFillColor(kWhite);
        annotation->SetTextAlign(12);
        annotation->SetTextSize(0.035);
        annotation->AddText(TString::Format("within slope = %.3f #pm %.3f ns/m",
                                            result.slope, result.slopeError));
        annotation->AddText(TString::Format("%lld hits in %d half-bars", result.entries, result.halfBars));
        annotation->Draw("SAME");
        std::cout << "[CDet y fixed effects] Layer " << layer + 1 << " " << sideName[side]
                  << ", " << sectionGroupName[group] << ": slope=" << result.slope
                  << " +/- " << result.slopeError << " ns/m, entries=" << result.entries
                  << ", half-bars=" << result.halfBars << ", NDF=" << result.ndf << "\n";
      }
    }
  }
  cCDetECalYFixedEffects->Update();

  if (applyYPropagationDiagnostic) {
    struct YCorrectionResult {
      double slope = NAN;
      double slopeError = NAN;
      double rms = NAN;
      long long entries = 0;
      int halfBars = 0;
    };
    auto fitCenteredSamples = [&](const std::vector<std::vector<std::pair<double,double>>> &bank,
                                  int wantedLayer, int wantedSide, int wantedGroup) {
      YCorrectionResult result;
      double xx = 0.0, xy = 0.0, centeredTimeSquares = 0.0;
      double residualSquares = 0.0;
      struct CenteredSample { double y, t; };
      std::vector<CenteredSample> centered;
      for (int halfBar = 0; halfBar < (int)bank.size(); ++halfBar) {
        const auto &samples = bank[halfBar];
        if ((int)samples.size() < minCDetYHalfBarEntries) continue;
        const int pixel = halfBar*NumPaddles;
        const int pixelsPerLayer = NumCDetPaddles/NumLayers;
        const int layer = pixel/pixelsPerLayer;
        const int withinLayer = pixel%pixelsPerLayer;
        const int side = withinLayer/NumCDetPaddlesPerSide;
        const int withinSide = withinLayer%NumCDetPaddlesPerSide;
        const int module = withinSide/(NumBars*NumPaddles);
        const int pmt = (withinSide%(NumBars*NumPaddles))/NumPaddles;
        const int section = 2*module + pmt/NumHalfBarsPerBank;
        const int group = std::min(section, 5 - section);
        if (layer != wantedLayer || side != wantedSide || group != wantedGroup) continue;
        double sumY = 0.0, sumT = 0.0;
        for (const auto &sample : samples) { sumY += sample.first; sumT += sample.second; }
        const double meanY = sumY/samples.size();
        const double meanT = sumT/samples.size();
        ++result.halfBars;
        result.entries += samples.size();
        for (const auto &sample : samples) {
          const double dy = sample.first - meanY;
          const double dt = sample.second - meanT;
          centered.push_back({dy, dt});
          xx += dy*dy;
          xy += dy*dt;
          centeredTimeSquares += dt*dt;
        }
      }
      if (xx <= 0.0) return result;
      result.slope = xy/xx;
      for (const auto &sample : centered) {
        const double residual = sample.t - result.slope*sample.y;
        residualSquares += residual*residual;
      }
      const long long ndf = result.entries - result.halfBars - 1;
      if (ndf > 0) {
        result.slopeError = std::sqrt((residualSquares/ndf)/xx);
        result.rms = std::sqrt(centeredTimeSquares/result.entries);
      }
      return result;
    };

    TCanvas *comparison = new TCanvas(
        "cCDetYPropagationCorrectionComparison",
        "Physics-motivated CDet y-propagation correction", 1500, 900);
    comparison->Divide(3,4);
    std::cout << "\n[CDet y-propagation correction comparison]\n"
              << "  n=" << yCorrectionRefractiveIndex
              << ", v=c/n=" << 100.0*yPropagationSpeedMPerNs << " cm/ns"
              << ", |dt/dy|=" << yPropagationSlopeMagnitude << " ns/m\n"
              << "  The accepted event/hit sample is identical before and after correction.\n";
    for (int layer = 0; layer < NumLayers; ++layer) {
      for (int side = 0; side < NumSides; ++side) {
        const int row = layer*NumSides + side;
        for (int group = 0; group < 3; ++group) {
          const YCorrectionResult before = fitCenteredSamples(
              yCorrectionBeforeSamples, layer, side, group);
          const YCorrectionResult after = fitCenteredSamples(
              yCorrectionAfterSamples, layer, side, group);
          comparison->cd(row*3 + group + 1);
          TH1F *frame = gPad->DrawFrame(-0.30, -10.0, 0.30, 10.0,
              TString::Format("Layer %d %s, %s;y_{ECal#rightarrowCDet}-#LTy#GT_{half-bar} (m);t-#LTt#GT_{half-bar} (ns)",
                              layer + 1, sideName[side], sectionGroupName[group]));
          (void)frame;
          TF1 *beforeLine = new TF1(
              TString::Format("fCDetYCorrectionBeforeL%dS%dG%d", layer, side, group),
              "[0]*x", -0.30, 0.30);
          beforeLine->SetParameter(0, before.slope);
          beforeLine->SetLineColor(kBlue + 1);
          beforeLine->SetLineWidth(2);
          beforeLine->Draw("SAME");
          TF1 *afterLine = new TF1(
              TString::Format("fCDetYCorrectionAfterL%dS%dG%d", layer, side, group),
              "[0]*x", -0.30, 0.30);
          afterLine->SetParameter(0, after.slope);
          afterLine->SetLineColor(kRed + 1);
          afterLine->SetLineWidth(2);
          afterLine->Draw("SAME");
          TPaveText *label = new TPaveText(0.12, 0.70, 0.62, 0.90, "NDC");
          label->SetFillColor(kWhite);
          label->SetTextAlign(12);
          label->AddText(TString::Format("before: %.3f #pm %.3f ns/m", before.slope, before.slopeError));
          label->AddText(TString::Format("after:  %.3f #pm %.3f ns/m", after.slope, after.slopeError));
          label->AddText(TString::Format("centered-time RMS: %.3f #rightarrow %.3f ns", before.rms, after.rms));
          label->Draw("SAME");
          std::cout << "  Layer " << layer + 1 << " " << sideName[side]
                    << ", " << sectionGroupName[group]
                    << ": slope " << before.slope << " +/- " << before.slopeError
                    << " -> " << after.slope << " +/- " << after.slopeError
                    << " ns/m; centered-time RMS " << before.rms << " -> " << after.rms
                    << " ns; entries=" << before.entries << "\n";
        }
      }
    }
    comparison->Update();
    comparison->SaveAs("CDet_y_propagation_correction_comparison.pdf");
    comparison->SaveAs("CDet_y_propagation_correction_comparison.png");
  }

  // TCanvas *cCDetLeVsTot = new TCanvas("cCDetLeVsTot", "CDet LE vs Tot",900,700);
  // cCDetLeVsTot->Divide(1,2);

  // cCDetLeVsTot->cd(1);
  // //gPad->SetLogz();
  // hCDet1LeVsTot->Draw();

  // cCDetLeVsTot->cd(2);
  // //gPad->SetLogz();
  // hCDet2LeVsTot->Draw();


  // // ----- plots for selected bars -----
  // TCanvas *cCDetLayerTimes1Bar = new TCanvas("cCDetLayerTimes1Bar", "CDet Selected Bars",900,700);
  // cCDetLayerTimes1Bar->Divide(1,2);

  // cCDetLayerTimes1Bar->cd(1);
  // //gPad->SetLogz();
  // hCDetBarLe1->Draw();

  // cCDetLayerTimes1Bar->cd(2);
  // hCDetBarLe2->Draw();

  TCanvas *cCDetLeVsTotBar = new TCanvas(
      "cCDetLeVsTotBar", "CDet LE vs TOT (Selected Bars)", 1200, 600);
  cCDetLeVsTotBar->Divide(2, 1);
  cCDetLeVsTotBar->cd(1);
  hCDet1BarLeVsTot->Draw("COLZ");
  cCDetLeVsTotBar->cd(2);
  hCDet2BarLeVsTot->Draw("COLZ");
  TCanvas *cCDetSelectedBarTimingStructures = new TCanvas(
      "cCDetSelectedBarTimingStructures",
      "Selected Bar vertical and diagonal timing components", 1400, 650);
  cCDetSelectedBarTimingStructures->Divide(2, 1);
  cCDetSelectedBarTimingStructures->cd(1);
  hCDet1BarVerticalLeVsTot->Draw("COLZ");
  cCDetSelectedBarTimingStructures->cd(2);
  hCDet1BarDiagonalLeVsTot->Draw("COLZ");
  cCDetSelectedBarTimingStructures->Update();
  TCanvas *cCDetAllBarsTimingStructures = new TCanvas(
      "cCDetAllBarsTimingStructures",
      "All-bar CDet timing populations", 1800, 600);
  cCDetAllBarsTimingStructures->Divide(3, 1);
  cCDetAllBarsTimingStructures->cd(1);
  hCDetAllLowerVerticalLeVsTot->Draw("COLZ");
  cCDetAllBarsTimingStructures->cd(2);
  hCDetAllUpperVerticalLeVsTot->Draw("COLZ");
  cCDetAllBarsTimingStructures->cd(3);
  hCDetAllDiagonalLeVsTot->Draw("COLZ");
  cCDetAllBarsTimingStructures->Update();
  auto printTimingStructure = [](const char *label, TH2D *histogram) {
    std::cout << "[CDet timing populations] " << label
              << ": entries=" << histogram->GetEntries()
              << ", mean ToT=" << histogram->GetMean(1)
              << ", ToT RMS=" << histogram->GetStdDev(1)
              << ", mean corrected LE=" << histogram->GetMean(2)
              << ", LE RMS=" << histogram->GetStdDev(2) << " ns\n";
  };
  printTimingStructure("lower vertical", hCDetAllLowerVerticalLeVsTot);
  printTimingStructure("upper vertical", hCDetAllUpperVerticalLeVsTot);
  printTimingStructure("diagonal", hCDetAllDiagonalLeVsTot);
  // ----- End plots for 1 bar -----

  TCanvas *cCDetTDiff = new TCanvas("cCDetTDiff", "CDet Time Diff",900,700);
  // ---- Gaussian fit on the "NoCuts" histogram ----
  // TF1 *fGaus = new TF1("fGaus", "gaus", DiffMin, DiffMax);
  // int maxBin = hCDetTimeDiff->GetMaximumBin();
  // double peakX = hCDetTimeDiff->GetBinCenter(maxBin);
  // fGaus->SetParameters(hCDetTimeDiff->GetMaximum(), peakX, diffMaxCut/2); // amp, mean, sigma guess

  // fGaus->SetRange(diffMinCut, diffMaxCut);
  // hCDetTimeDiff->Fit(fGaus, "R");

  // draw the fit on top of the already-drawn histogram
  // fGaus->SetLineColor(kBlack);

  hCDetTimeDiff->Draw();
  if (applyYPropagationDiagnostic) {
    TCanvas *cCDetTDiffYCorrected = new TCanvas(
        "cCDetTDiffYCorrected", "CDet Time Diff after y-propagation correction", 900, 700);
    hCDetTimeDiffYCorrected->Draw();
    cCDetTDiffYCorrected->Update();
    cCDetTDiffYCorrected->SaveAs("cCDetTDiffYCorrected.jpg");
    std::cout << "[CDet layer timing before/after y correction]\n"
              << "  Layer 1 RMS: " << hCDetLe1->GetStdDev() << " -> "
              << hCDetLe1YCorrected->GetStdDev() << " ns\n"
              << "  Layer 2 RMS: " << hCDetLe2->GetStdDev() << " -> "
              << hCDetLe2YCorrected->GetStdDev() << " ns\n"
              << "  Layer 2 - Layer 1 RMS: " << hCDetTimeDiff->GetStdDev()
              << " -> " << hCDetTimeDiffYCorrected->GetStdDev() << " ns\n";
  }
  // fGaus->Draw("SAME");

  // TCanvas *cCDetXDiff = new TCanvas("cCDetXDiff", "CDet X Diff",900,700);
  // hCDetXDiff->Draw();

  TCanvas *cCDetLayer2v1 = new TCanvas("cCDetLayer2v1", "CDet Layer 2 vs 1",900,700);
  hCDetLe2vs1->SetMinimum(20);
  hCDetLe2vs1->Draw("COLZ");

  TCanvas *cCDetTotLayer2v1 = new TCanvas("cCDetTotLayer2v1", "CDet Tot Layer 2 vs 1",900,700);
  hCDetTot2vs1->SetMinimum(40);
  hCDetTot2vs1->Draw("COLZ");

  TCanvas *cNpair = new TCanvas("cNpair", "CDet accepted pairs per event", 900, 700);
  hNpairPerEvent->Draw();

  // TCanvas *cXplot = new TCanvas("cXplot", "CDet layer 2 vs 1 xposition", 900, 700);
  // h2CDetx2VsCDetx1->Draw("COLZ");

  // TCanvas * cDTvsECalT = new TCanvas("cDTvsECalT", " CDet dt vs ECal t", 900,700);
  // hECalVsCDetDt->Draw("COLZ");

  // TCanvas * cDTvsECalTSingle = new TCanvas("cDTvsECalTSingle", " CDet Single dt vs ECal t", 900,700);
  // hECalVsCDetDtSingle->Draw("COLZ");

  TCanvas * cCDetTvsECalT = new TCanvas("cCDetTvsECalT", " CDet t vs ECal t", 900,700);
  hECalVsCDetT->SetMinimum(40);
  hECalVsCDetT->Draw("COLZ");
  hECalVsCDetT->GetYaxis()->SetRangeUser(20, 40);

  // ---------------------------------------------
  // Linear parametrization: <CDet t> vs ECal ADC t
  // Use a profile along X (ECal time) and fit with a straight line.
  // This yields a robust mean-trend fit for the 2D distribution.
  // ---------------------------------------------
  TProfile* pCDetTvsECalT = hECalVsCDetT->ProfileX("pCDetTvsECalT");
  pCDetTvsECalT->SetMarkerStyle(20);
  pCDetTvsECalT->GetYaxis()->SetRangeUser(20, 40);
  pCDetTvsECalT->SetMarkerSize(0.6);

  TF1* fCDetTvsECalT_lin = new TF1("fCDetTvsECalT_lin", "pol1", ECalMin, ECalMax);
  // Quiet fit, respect range
  const int ecalFitStatus = pCDetTvsECalT->Fit(fCDetTvsECalT_lin, "QRS");

  // Overlay the profile points and fit on top of the COLZ plot
  pCDetTvsECalT->Draw("SAME");
  fCDetTvsECalT_lin->Draw("SAME");

  const double p0 = fCDetTvsECalT_lin->GetParameter(0);
  const double p1 = fCDetTvsECalT_lin->GetParameter(1);
  const double e0 = fCDetTvsECalT_lin->GetParError(0);
  const double e1 = fCDetTvsECalT_lin->GetParError(1);
  const double ecalChi2 = fCDetTvsECalT_lin->GetChisquare();
  const int ecalNdf = fCDetTvsECalT_lin->GetNDF();
  const double ecalChi2Ndf = ecalNdf > 0 ? ecalChi2/ecalNdf : 0.0;
  const double ecalProb = ecalNdf > 0 ? TMath::Prob(ecalChi2, ecalNdf) : 0.0;
  std::cout << "\n[plotCDetLayersTimeComp] Linear fit for <CDet t> vs ECal ADC t:\n"
            << "  <t_CDet> = p0 + p1 * t_ECal\n"
            << "  p0 = " << p0 << " +/- " << e0 << " ns\n"
            << "  p1 = " << p1 << " +/- " << e1 << " (ns/ns)\n"
            << "  chi2/NDF = " << ecalChi2 << "/" << ecalNdf
            << " = " << ecalChi2Ndf << "  probability = " << ecalProb << "\n\n";

  TPaveText* ptFit = new TPaveText(0.14, 0.80, 0.52, 0.92, "NDC");
  ptFit->SetFillColor(0);
  ptFit->SetTextAlign(12);
  ptFit->AddText("<t_{CDet}> = p0 + p1 t_{ECal}");
  ptFit->AddText(Form("p0 = %.3f #pm %.3f ns", p0, e0));
  ptFit->AddText(Form("p1 = %.5f #pm %.5f", p1, e1));
  ptFit->AddText(Form("#chi^{2}/NDF = %.1f/%d = %.2f", ecalChi2, ecalNdf, ecalChi2Ndf));
  ptFit->Draw("SAME");

  // Compare the all-pair trend with the two layers and the selected Layer-1 bar.
  // These pooled fits are diagnostics only.  They mix within-half-bar timing
  // dependence with differences among half-bar means (Simpson's paradox), so
  // the overwrite decision below uses the fixed-effects result instead.
  TCanvas *cCDetTvsECalTDiagnostics = new TCanvas("cCDetTvsECalTDiagnostics", "CDet/ECal timing trend diagnostics", 1200, 900);
  cCDetTvsECalTDiagnostics->Divide(2,2);
  cCDetTvsECalTDiagnostics->cd(1);
  hECalVsCDetT->Draw("COLZ");
  pCDetTvsECalT->Draw("SAME");
  fCDetTvsECalT_lin->Draw("SAME");

  auto drawTimingTrendDiagnostic = [&](int pad, TH2D *hist, const char *profileName,
                                       const char *fitName, const char *label) {
    cCDetTvsECalTDiagnostics->cd(pad);
    hist->SetMinimum(40);
    hist->GetYaxis()->SetRangeUser(20, 40);
    hist->Draw("COLZ");
    TProfile *profile = hist->ProfileX(profileName);
    profile->SetMarkerStyle(20);
    profile->SetMarkerSize(0.6);
    profile->GetYaxis()->SetRangeUser(20, 40);
    TF1 *fit = new TF1(fitName, "pol1", ECalMin, ECalMax);
    const bool enoughData = profile->GetEntries() >= 3;
    int status = -1;
    if (enoughData) status = profile->Fit(fit, "QRS");
    profile->Draw("SAME");
    if (status == 0) fit->Draw("SAME");
    const double chi2 = fit->GetChisquare();
    const int ndf = fit->GetNDF();
    const double chi2Ndf = ndf > 0 ? chi2/ndf : 0.0;
    TPaveText *summary = new TPaveText(0.13, 0.78, 0.58, 0.92, "NDC");
    summary->SetFillColor(0);
    summary->SetTextAlign(12);
    summary->AddText(label);
    if (status == 0) {
      summary->AddText(Form("p1 = %.5f #pm %.5f", fit->GetParameter(1), fit->GetParError(1)));
      summary->AddText(Form("#chi^{2}/NDF = %.1f/%d = %.2f", chi2, ndf, chi2Ndf));
    } else {
      summary->AddText("Insufficient data or fit failed");
    }
    summary->Draw("SAME");
    std::cout << "[plotCDetLayersTimeComp] " << label
              << ": p1=" << fit->GetParameter(1) << " +/- " << fit->GetParError(1)
              << ", chi2/NDF=" << chi2 << "/" << ndf << " = " << chi2Ndf
              << ", fit_status=" << status << "\n";
  };
  drawTimingTrendDiagnostic(2, hECalVsCDetTL1, "pCDetL1TvsECalT", "fCDetL1TvsECalT_lin", "All Layer 1 hits");
  drawTimingTrendDiagnostic(3, hECalVsCDetTL2, "pCDetL2TvsECalT", "fCDetL2TvsECalT_lin", "All Layer 2 hits");
  drawTimingTrendDiagnostic(4, hECalVsSelectedBarT, "pSelectedBarTvsECalT", "fSelectedBarTvsECalT_lin", TString::Format("Layer 1 bar %d", selectedBarNumber).Data());
  cCDetTvsECalTDiagnostics->Update();

  if (gUseECalTimeCorr) {
    TCanvas *cCDetTvsECalTBeforeP1 = new TCanvas(
        "cCDetTvsECalTBeforeP1",
        "CDet/ECal timing before the p1 correction", 1200, 900);
    cCDetTvsECalTBeforeP1->Divide(2,2);
    auto drawBeforeP1 = [&](int pad, TH2D *histogram, const char *profileName,
                            const char *fitName, const char *label) {
      cCDetTvsECalTBeforeP1->cd(pad);
      histogram->SetMinimum(40);
      // Build the profile before restricting the displayed y range.  ROOT's
      // ProfileX honors an active y-axis range, which would otherwise clip
      // tails in an ECal-dependent way and bias this presentation slope.
      histogram->GetYaxis()->SetRange(0, 0);
      TProfile *profile = histogram->ProfileX(profileName);
      histogram->GetYaxis()->SetRangeUser(
          std::max(beforeP1TimeMin, 20.0), std::min(beforeP1TimeMax, 70.0));
      histogram->Draw("COLZ");
      profile->SetMarkerStyle(20);
      profile->SetMarkerSize(0.65);
      profile->SetMarkerColor(kBlack);
      profile->SetLineColor(kBlack);
      TF1 *fit = new TF1(fitName, "pol1", ECalMin, ECalMax);
      int status = -1;
      if (profile->GetEntries() >= 3) status = profile->Fit(fit, "QRS");
      profile->Draw("SAME");
      if (status == 0) fit->Draw("SAME");
      TPaveText *summary = new TPaveText(0.13, 0.75, 0.61, 0.92, "NDC");
      summary->SetFillColor(kWhite);
      summary->SetTextAlign(12);
      summary->AddText(label);
      if (status == 0) {
        const int ndf = fit->GetNDF();
        summary->AddText(Form("p1 = %.5f #pm %.5f", fit->GetParameter(1),
                              fit->GetParError(1)));
        summary->AddText(Form("#chi^{2}/NDF = %.1f/%d = %.2f",
                              fit->GetChisquare(), ndf,
                              ndf > 0 ? fit->GetChisquare()/ndf : 0.0));
        std::cout << "[CDet before-p1 presentation] " << label
                  << ": p1=" << fit->GetParameter(1) << " +/- "
                  << fit->GetParError(1) << " ns/ns\n";
      } else {
        summary->AddText("Insufficient data or fit failed");
      }
      summary->Draw("SAME");
    };
    drawBeforeP1(1, hECalVsCDetTBeforeP1,
                 "pCDetTvsECalTBeforeP1", "fCDetTvsECalTBeforeP1",
                 "All accepted pairs");
    drawBeforeP1(2, hECalVsCDetTL1BeforeP1,
                 "pCDetL1TvsECalTBeforeP1", "fCDetL1TvsECalTBeforeP1",
                 "All Layer 1 hits");
    drawBeforeP1(3, hECalVsCDetTL2BeforeP1,
                 "pCDetL2TvsECalTBeforeP1", "fCDetL2TvsECalTBeforeP1",
                 "All Layer 2 hits");
    drawBeforeP1(4, hECalVsSelectedBarTBeforeP1,
                 "pSelectedBarTvsECalTBeforeP1", "fSelectedBarTvsECalTBeforeP1",
                 TString::Format("Layer 1 bar %d", selectedBarNumber).Data());
    cCDetTvsECalTBeforeP1->cd(1);
    TPaveText *method = new TPaveText(0.13, 0.64, 0.61, 0.73, "NDC");
    method->SetFillColor(kWhite);
    method->SetTextAlign(12);
    method->AddText(Form("p1 correction removed: + %.6f t_{ECal}",
                         gECalFitP1));
    method->Draw("SAME");
    cCDetTvsECalTBeforeP1->Update();
    cCDetTvsECalTBeforeP1->SaveAs("cCDetTvsECalTBeforeP1.jpg");
    cCDetTvsECalTBeforeP1->SaveAs("cCDetTvsECalTBeforeP1.pdf");
  }

  // Fixed-effects regression for the ECal timing slope.  Center x and y
  // independently within every half-bar, then pool those centered samples.
  // This is equivalent to fitting a common slope with a separate intercept
  // for each half-bar and removes the between-half-bar covariance that makes
  // the unstratified profile appear artificially flat.
  struct FixedEffectsResult {
    bool valid = false;
    long long entries = 0;
    int groups = 0;
    double slope = std::numeric_limits<double>::quiet_NaN();
    double slopeError = std::numeric_limits<double>::quiet_NaN();
    double withinCovariance = 0.0;
    double withinVarianceX = 0.0;
    double residualSumSquares = 0.0;
    long long ndf = 0;
  };

  auto calculateFixedEffects = [&](const std::vector<std::vector<std::pair<double,double>>> &sampleBank,
                                   int firstIndex, int lastIndex,
                                   TH2D *centeredHistogram) {
    FixedEffectsResult result;
    struct GroupMean { int index; double weight; double x; double y; };
    std::vector<GroupMean> means;
    for (int index = firstIndex; index < lastIndex; ++index) {
      const std::vector<std::pair<double,double>> &samples =
          sampleBank[index];
      if (samples.size() < 3) continue;
      double weight = 0.0, sumX = 0.0, sumY = 0.0;
      for (const std::pair<double,double> &sample : samples) {
        weight += 1.0;
        sumX += sample.first;
        sumY += sample.second;
      }
      means.push_back({index, weight, sumX/weight, sumY/weight});
      result.entries += std::llround(weight);
    }
    result.groups = means.size();
    for (const GroupMean &group : means) {
      const std::vector<std::pair<double,double>> &samples =
          sampleBank[group.index];
      for (const std::pair<double,double> &sample : samples) {
        const double dx = sample.first - group.x;
        const double dy = sample.second - group.y;
        result.withinCovariance += dx*dy;
        result.withinVarianceX += dx*dx;
        if (centeredHistogram) centeredHistogram->Fill(dx, dy);
      }
    }
    result.ndf = result.entries - result.groups - 1;
    if (result.withinVarianceX <= 0.0 || result.ndf <= 0) return result;
    result.slope = result.withinCovariance/result.withinVarianceX;
    for (const GroupMean &group : means) {
      const std::vector<std::pair<double,double>> &samples =
          sampleBank[group.index];
      for (const std::pair<double,double> &sample : samples) {
        const double dx = sample.first - group.x;
        const double dy = sample.second - group.y;
        const double residual = dy - result.slope*dx;
        result.residualSumSquares += residual*residual;
      }
    }
    result.slopeError = std::sqrt(
        (result.residualSumSquares/result.ndf)/result.withinVarianceX);
    result.valid = std::isfinite(result.slope) &&
                   std::isfinite(result.slopeError) && result.slopeError > 0.0;
    return result;
  };

  TH2D *hCDetL1ECalWithin = new TH2D("hCDetL1ECalWithin",
      "Layer 1 within-half-bar timing;ECal time - half-bar mean (ns);CDet time - half-bar mean (ns)",
      120, -30, 30, 120, -15, 15);
  TH2D *hCDetL2ECalWithin = new TH2D("hCDetL2ECalWithin",
      "Layer 2 within-half-bar timing;ECal time - half-bar mean (ns);CDet time - half-bar mean (ns)",
      120, -30, 30, 120, -15, 15);
  TH2D *hCDetAllECalWithin = new TH2D("hCDetAllECalWithin",
      "Combined within-half-bar timing;ECal time - half-bar mean (ns);CDet time - half-bar mean (ns)",
      120, -30, 30, 120, -15, 15);
  const int halfBarsPerLayer = barsPerLayer;
  const FixedEffectsResult fixedL1 = calculateFixedEffects(
      gCDetHalfBarECalTimingSamples, 0, halfBarsPerLayer, hCDetL1ECalWithin);
  const FixedEffectsResult fixedL2 = calculateFixedEffects(
      gCDetHalfBarECalTimingSamples, halfBarsPerLayer, 2*halfBarsPerLayer, hCDetL2ECalWithin);
  const FixedEffectsResult fixedAll = calculateFixedEffects(
      gCDetHalfBarECalTimingSamples, 0, 2*halfBarsPerLayer, hCDetAllECalWithin);
  gLastECalFixedEffectsSlope = fixedAll.slope;
  gLastECalFixedEffectsSlopeError = fixedAll.slopeError;
  gLastECalFixedEffectsValid = fixedAll.valid;

  auto printFixedEffects = [&](const char *label, const FixedEffectsResult &fit) {
    std::cout << "[plotCDetLayersTimeComp] " << label
              << " fixed-effects slope: p1=" << fit.slope
              << " +/- " << fit.slopeError << " ns/ns"
              << ", entries=" << fit.entries << ", half-bars=" << fit.groups
              << ", NDF=" << fit.ndf << "\n";
  };
  std::cout << "\n[CDet/ECal fixed-effects timing fit]\n"
            << "  Separate intercepts are removed for every half-bar.\n";
  printFixedEffects("Layer 1", fixedL1);
  printFixedEffects("Layer 2", fixedL2);
  printFixedEffects("Combined", fixedAll);
  std::cout << "  The combined fixed-effects slope, not the pooled profile slope, "
            << "is used for calibration updates.\n\n";

  TCanvas *cCDetECalFixedEffects = new TCanvas(
      "cCDetECalFixedEffects", "CDet/ECal within-half-bar fixed-effects timing", 1500, 500);
  cCDetECalFixedEffects->Divide(3, 1);
  auto drawFixedEffects = [&](int pad, TH2D *hist, const FixedEffectsResult &fit,
                              const char *label) {
    cCDetECalFixedEffects->cd(pad);
    hist->Draw("COLZ");
    if (fit.valid) {
      TF1 *line = new TF1(TString::Format("fCDetECalWithin%d", pad), "[0]*x", -30, 30);
      line->SetParameter(0, fit.slope);
      line->SetLineColor(kRed);
      line->Draw("SAME");
    }
    TPaveText *box = new TPaveText(0.13, 0.78, 0.61, 0.92, "NDC");
    box->SetFillColor(0);
    box->SetTextAlign(12);
    box->AddText(label);
    if (fit.valid) {
      box->AddText(Form("within p1 = %.5f #pm %.5f", fit.slope, fit.slopeError));
      box->AddText(Form("%lld hits in %d half-bars", fit.entries, fit.groups));
    } else {
      box->AddText("fixed-effects fit unavailable");
    }
    box->Draw("SAME");
  };
  drawFixedEffects(1, hCDetL1ECalWithin, fixedL1, "Layer 1 fixed effects");
  drawFixedEffects(2, hCDetL2ECalWithin, fixedL2, "Layer 2 fixed effects");
  drawFixedEffects(3, hCDetAllECalWithin, fixedAll, "Combined fixed effects");
  cCDetECalFixedEffects->Update();

  TH2D *hCDetAllHCalWithin = new TH2D("hCDetAllHCalWithin",
      "Combined within-half-bar HCal timing;HCal time - half-bar mean (ns);CDet time - half-bar mean (ns)",
      140, -35, 35, 120, -15, 15);
  const FixedEffectsResult fixedHCalAll = calculateFixedEffects(
      halfBarHCalTimingSamples, 0, 2*halfBarsPerLayer, hCDetAllHCalWithin);
  TProfile *pCDetTvsHCalT = hHCalVsCDetT->ProfileX("pCDetTvsHCalT");
  TF1 *fCDetTvsHCalT = new TF1("fCDetTvsHCalT", "pol1", HCalMin, HCalMax);
  int hcalFitStatus = -1;
  if (pCDetTvsHCalT->GetEntries() >= 3)
    hcalFitStatus = pCDetTvsHCalT->Fit(fCDetTvsHCalT, "QRS");
  TCanvas *cCDetTvsHCalT = new TCanvas(
      "cCDetTvsHCalT", "CDet/HCal timing control", 1500, 900);
  cCDetTvsHCalT->Divide(2, 2);
  cCDetTvsHCalT->cd(1); hHCalVsCDetT->Draw("COLZ");
  pCDetTvsHCalT->SetMarkerStyle(20); pCDetTvsHCalT->Draw("SAME");
  if (hcalFitStatus == 0) fCDetTvsHCalT->Draw("SAME");
  cCDetTvsHCalT->cd(2); hCDetAllHCalWithin->Draw("COLZ");
  if (fixedHCalAll.valid) {
    TF1 *line = new TF1("fCDetHCalWithin", "[0]*x", -35, 35);
    line->SetParameter(0, fixedHCalAll.slope);
    line->SetLineColor(kRed); line->Draw("SAME");
  }
  cCDetTvsHCalT->cd(3); hECalVsHCalAccepted->Draw("COLZ");
  cCDetTvsHCalT->cd(4); hHCalVsCDetT->ProjectionX("hAcceptedHCalTime")->Draw();
  cCDetTvsHCalT->Update();
  std::cout << "\n[CDet/HCal timing control]\n"
            << "  HCal window: " << HCalMin << " to " << HCalMax << " ns\n";
  if (hcalFitStatus == 0)
    std::cout << "  pooled profile slope: " << fCDetTvsHCalT->GetParameter(1)
              << " +/- " << fCDetTvsHCalT->GetParError(1) << " ns/ns\n";
  printFixedEffects("Combined CDet/HCal", fixedHCalAll);

  const bool pooledFitOK = ecalFitStatus == 0 && ecalNdf > 0 &&
                           pCDetTvsECalT->GetEntries() >= 3 &&
                           std::isfinite(p0) && std::isfinite(p1);
  // A residual-slope closure test needs only the fixed-effects result.  The
  // uncorrected absolute-calibration stage also needs the pooled intercept.
  const bool ecalFitOK = fixedAll.valid && (gUseECalTimeCorr || pooledFitOK);
  if (overwrite && ecalFitOK) {
    if (gUseECalTimeCorr) {
      // The corrected distribution is recentered to gTargetMeanLE, so its
      // fitted intercept cannot determine the stored absolute p0.  Only the
      // residual slope remains identifiable in this mode.
      gECalFitP1 += fixedAll.slope;
      std::cout << "[CDet] ECal correction already active: added within-half-bar residual slope "
                << fixedAll.slope << " to p1 and left p0 unchanged (the "
                << gTargetMeanLE << " ns recentering removes p0 sensitivity).\n";
    } else {
      // The fixed-effects regression determines the physical common slope;
      // the pooled fit supplies the absolute intercept in this uncorrected
      // stage because demeaning intentionally removes intercept information.
      gECalFitP0 = p0;
      gECalFitP1 = fixedAll.slope;
      double sumResidual = 0.0;
      long long residualCount = 0;
      const size_t nEvents = std::min(vGoodLe.size(), v_GoodECalAdcTime.size());
      for (size_t event = 0; event < nEvents; ++event) {
        const double ecalTime = v_GoodECalAdcTime[event];
        const size_t nHits = std::min(vGoodLe[event].size(), vGoodID[event].size());
        for (size_t hit = 0; hit < nHits; ++hit) {
          const int bar = vGoodID[event][hit] / NumPaddles;
          if (bar < 0 || bar >= NumPMTs) continue;
          sumResidual += vGoodLe[event][hit] -
                         (gECalFitP0 + gECalFitP1*ecalTime);
          ++residualCount;
        }
      }
      if (residualCount <= 0) {
        std::cerr << "[CDet] ERROR: cannot derive ECal delta from an empty accepted sample.\n";
        gLastCalibrationFitSucceeded = false;
        return;
      }
      gECalDeltaShift = gTargetMeanLE - sumResidual/(double)residualCount;
      gECalDeltaLoaded = true;
      std::cout << "[CDet] ECal correction disabled: stored pooled absolute p0 and "
                << "fixed-effects p1, with detector-wide delta="
                << gECalDeltaShift << " ns.\n";
    }
    gECalParamsLoaded = true;
    gLastCalibrationFitSucceeded = WriteCalibrationConstants(gCalibrationFile);
    std::cout << "[CDet] Updated ECal timing parameters in calibration file: "
              << "p0=" << gECalFitP0 << "  p1=" << gECalFitP1 << "\n";
  } else if (overwrite) {
    std::cerr << "[CDet] ERROR: ECal timing fit failed; constants were not written.\n";
  } else {
    gLastCalibrationFitSucceeded = ecalFitOK;
  }

  TCanvas* cCDetFrontvsBackID = new TCanvas("cCDetFrontvsBackID", "CDet Front Paddle vs Back Paddle ID", 1000, 620);
  hCDet1IDvs2ID->Draw("COLZ");

  // TCanvas * cCDetTvsECalTSingle = new TCanvas("cCDetTvsECalTSingle", " CDet Single t vs ECal t", 900,700);
  // hECalVsCDetTSingle->Draw("COLZ");

  // TCanvas *cDtCDetECal = new TCanvas("cDtCDetECal", "CDet ECal dt",900,700);
  // hDtCDetECal->Draw();

  // TCanvas *cDtvsDxCDetECal = new TCanvas("cDtvsDxCDetECal", "CDet ECal dt vs dx",900,700);
  // hDtvsDxCDetECal->Draw("COLZ");

  TCanvas *cDtCDetLayersvsPos = new TCanvas("cDtCDetLayersvsPos", "CDet dt vs position", 900,700);
  cDtCDetLayersvsPos->Divide(2,2);
  cDtCDetLayersvsPos->cd(1);

  hCDetTimeDiffvsx1->Draw();

  cDtCDetLayersvsPos->cd(2);
  hCDetTimeDiffvsx2->Draw();

  cDtCDetLayersvsPos->cd(3);
  hCDetTimeDiffvsy1->Draw();

  cDtCDetLayersvsPos->cd(4);
  hCDetTimeDiffvsy2->Draw();

  TCanvas *cSelectedBarXMap = new TCanvas("cSelectedBarXMap", TString::Format("CDet and ECal x map for Layer 1 bar %d", selectedBarNumber), 1200,1000);
  cSelectedBarXMap->Divide(2,2);
  cSelectedBarXMap->cd(1);
  hSelectedBarX2VsX1->Draw("COLZ");
  cSelectedBarXMap->cd(2);
  hSelectedBarECalXVsX1->Draw("COLZ");
  cSelectedBarXMap->cd(3);
  hSelectedBarECalXVsX2->Draw("COLZ");
  cSelectedBarXMap->cd(4);
  hSelectedBarXByPlane->Draw("COLZ");

  TCanvas *cSelectedBarXVsZ = new TCanvas("cSelectedBarXVsZ", TString::Format("CDet-to-ECal x-z map for Layer 1 bar %d", selectedBarNumber), 1100,700);
  hSelectedBarXVsZ->Draw("COLZ");

  // The event browser is an interactive diagnostic.  Avoid constructing it in
  // batch calibration jobs, where the extra canvases and retained snapshots
  // would only consume memory.
  if (!gROOT->IsBatch()) {
    BuildCDetEventDisplay(selectedBarNumber, diffMinCut, diffMaxCut,
                          xdiffMinCut, xdiffMaxCut,
                          tdiffECalCDetMin, tdiffECalCDetMax,
                          XMin, XMax, ZMin, ZMax, displayBestHits,
                          bestHitPeakMean, bestHitPeakSigma, bestHitNSigma);
  }
}

void reportCDetPairedTimeResolution(bool draw,
                                    double gaussianFitMin,
                                    double gaussianFitMax)
{
  gLastPairedCoreMean = std::numeric_limits<double>::quiet_NaN();
  gLastPairedCoreMeanError = std::numeric_limits<double>::quiet_NaN();
  gLastPairedCoreSigma = std::numeric_limits<double>::quiet_NaN();
  gLastPairedCoreSigmaError = std::numeric_limits<double>::quiet_NaN();
  gLastPairedTimeEntries = 0.0;
  gLastPairedCoreFitValid = false;
  if (!gCDetPairedMeanTimeVsECal) {
    std::cerr << "[CDet paired timing resolution] ERROR: paired timing histogram "
              << "is not available. Run plotCDetLayersTimeComp first.\n";
    return;
  }
  if (!(gaussianFitMax > gaussianFitMin)) {
    std::cerr << "[CDet paired timing resolution] ERROR: Gaussian fit maximum "
              << "must be greater than its minimum.\n";
    return;
  }

  const double entries = static_cast<double>(gCDetAcceptedPairMeanTimes.size());
  gLastPairedTimeEntries = entries;
  double mean = 0.0;
  for (double value : gCDetAcceptedPairMeanTimes) mean += value;
  if (entries > 0.0) mean /= entries;

  // Anchor the reporting histogram to the unbinned sample mean. A global
  // additive timing shift then moves the samples, bin edges, and fit interval
  // together, so the fitted correction cannot depend on the starting shift's
  // phase relative to a fixed absolute-time bin grid.
  static unsigned long reportNumber = 0;
  ++reportNumber;
  const int projectionBins = gCDetPairedMeanTimeVsECal->GetNbinsY();
  const double projectionWidth =
      gCDetPairedMeanTimeVsECal->GetYaxis()->GetBinWidth(1);
  const double projectionSpan = projectionBins * projectionWidth;
  TH1D *timeProjection = new TH1D(
      TString::Format("hCDetPairedMeanTimeResolution_%lu", reportNumber),
      "Accepted-pair mean CDet time;#LTt_{CDet}#GT (ns);Accepted pairs",
      projectionBins, mean - 0.5*projectionSpan, mean + 0.5*projectionSpan);
  timeProjection->SetDirectory(nullptr);
  for (double value : gCDetAcceptedPairMeanTimes) timeProjection->Fill(value);
  const double rms = timeProjection->GetStdDev();
  const double rmsError = timeProjection->GetStdDevError();

  TF1 *coreFit = new TF1(
      TString::Format("fCDetPairedMeanTimeCore_%lu", reportNumber),
      "gaus", gaussianFitMin, gaussianFitMax);
  int fitStatus = -1;
  if (timeProjection->Integral(
          timeProjection->FindFixBin(gaussianFitMin),
          timeProjection->FindFixBin(gaussianFitMax)) >= 10.0) {
    TFitResultPtr fitResult = timeProjection->Fit(coreFit, "QRS");
    fitStatus = int(fitResult);
  }

  std::cout << "\n[CDet paired timing resolution]\n"
            << "  accepted pairs: " << entries << "\n"
            << "  mean pair time: " << mean << " ns\n"
            << "  detector-wide RMS: " << rms << " +/- " << rmsError
            << " ns\n";
  if (fitStatus == 0) {
    gLastPairedCoreMean = coreFit->GetParameter(1);
    gLastPairedCoreMeanError = coreFit->GetParError(1);
    gLastPairedCoreSigma = coreFit->GetParameter(2);
    gLastPairedCoreSigmaError = coreFit->GetParError(2);
    gLastPairedCoreFitValid =
        std::isfinite(gLastPairedCoreMean) &&
        std::isfinite(gLastPairedCoreMeanError) &&
        std::isfinite(gLastPairedCoreSigma) && gLastPairedCoreSigma > 0.0;
    std::cout << "  Gaussian core mean (" << gaussianFitMin << " to "
              << gaussianFitMax << " ns): " << gLastPairedCoreMean
              << " +/- " << gLastPairedCoreMeanError << " ns\n"
              << "  Gaussian core sigma (" << gaussianFitMin << " to "
              << gaussianFitMax << " ns): " << coreFit->GetParameter(2)
              << " +/- " << coreFit->GetParError(2) << " ns\n"
              << "  Gaussian chi2/NDF: " << coreFit->GetChisquare() << "/"
              << coreFit->GetNDF() << "\n";
  } else {
    std::cout << "  Gaussian core fit unavailable (fit status " << fitStatus
              << ").\n";
  }
  std::cout << "  The RMS is the width of the full accepted-pair mean-time "
            << "distribution; the Gaussian sigma describes only its core.\n";

  if (draw) {
    TCanvas *canvas = new TCanvas(
        TString::Format("cCDetPairedTimeResolution_%lu", reportNumber),
        "CDet paired-layer timing resolution", 900, 700);
    canvas->cd();
    timeProjection->Draw();
    if (fitStatus == 0) coreFit->Draw("SAME");
    canvas->Modified();
    canvas->Update();
  }
}

void writeAllCDetPairedLeSpectra(TString outputFile)
{
  if ((int)gCDetPairedLeSpectra.size() != NumCDetPaddles) {
    std::cerr << "[CDet paired LE writer] ERROR: paired spectra are not available. "
              << "Run plotCDetLayersTimeComp first.\n";
    return;
  }
  if (outputFile.IsNull()) {
    std::cerr << "[CDet paired LE writer] ERROR: output filename is empty.\n";
    return;
  }

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "[CDet paired LE writer] ERROR: cannot create " << outputFile << ".\n";
    return;
  }

  const int barsPerLayer = NumCDetPaddles/(NumLayers*NumPaddles);
  int fittedSpectra = 0;
  for (int bar = 0; bar < barsPerLayer; ++bar) {
    TDirectory *barDirectory = output.mkdir(TString::Format("bar_%03d", bar));
    barDirectory->cd();
    TCanvas *canvas = new TCanvas(
        TString::Format("cPairedLeBar%03d", bar),
        TString::Format("Paired-hit LE spectra for bar %d, Layers 1 and 2", bar),
        1800, 1000);
    canvas->Divide(8, 4, 0.001, 0.001);

    for (int layer = 0; layer < NumLayers; ++layer) {
      TDirectory *layerDirectory = barDirectory->mkdir(TString::Format("layer_%d", layer + 1));
      for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
        const int pixel = layer*(NumCDetPaddles/NumLayers) + bar*NumPaddles + localPixel;
        TH1D *histogram = gCDetPairedLeSpectra[pixel];
        if (!histogram) continue;

        barDirectory->cd();
        canvas->cd(layer*NumPaddles + localPixel + 1);
        histogram->Draw();

        const double fitMin = gTargetMeanLE - 5.0;
        const double fitMax = gTargetMeanLE + 5.0;
        const int firstFitBin = histogram->GetXaxis()->FindBin(fitMin);
        const int lastFitBin = histogram->GetXaxis()->FindBin(fitMax);
        if (histogram->Integral(firstFitBin, lastFitBin) > 20) {
          int peakBin = firstFitBin;
          for (int bin = firstFitBin + 1; bin <= lastFitBin; ++bin) {
            if (histogram->GetBinContent(bin) > histogram->GetBinContent(peakBin)) peakBin = bin;
          }
          TF1 *fit = new TF1(
              TString::Format("fPairedGoodLe_Pixel%04d", pixel), "gaus", fitMin, fitMax);
          fit->SetParameters(histogram->GetBinContent(peakBin), histogram->GetBinCenter(peakBin), 2.0);
          if (histogram->Fit(fit, "RQ0") == 0) {
            fit->Draw("SAME");
            ++fittedSpectra;
          }
        }
        if (kUnusedCDetPixels.count(pixel)) {
          TPaveText *flag = new TPaveText(0.72,0.80,0.95,0.94,"NDC");
          flag->SetFillColor(kBlack);
          flag->SetLineColor(kBlack);
          flag->SetTextColor(kWhite);
          flag->SetBorderSize(1);
          flag->AddText("UNUSED PIXEL");
          flag->Draw("SAME");
        }
        layerDirectory->cd();
        histogram->Write();
        barDirectory->cd();
      }
    }
    barDirectory->cd();
    canvas->Write("paired_le_layers_1_and_2");
    delete canvas;
    output.cd();
  }
  output.Close();
  std::cout << "[CDet paired LE writer] Wrote " << barsPerLayer
            << " two-layer bar canvases and all " << NumCDetPaddles
            << " paired-hit spectra to " << outputFile
            << " (" << fittedSpectra << " successful 25-35 ns fits).\n";
}

void showCDetPairedLeBar(int bar, TString inputFile)
{
  const int barsPerLayer = NumCDetPaddles/(NumLayers*NumPaddles);
  if (bar < 0 || bar >= barsPerLayer) {
    std::cerr << "[CDet paired LE viewer] ERROR: bar must be in [0, "
              << barsPerLayer - 1 << "].\n";
    return;
  }
  TFile input(inputFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "[CDet paired LE viewer] ERROR: cannot open " << inputFile << ".\n";
    return;
  }

  static unsigned long invocation = 0;
  ++invocation;
  std::vector<TH1D*> histograms(NumLayers*NumPaddles, nullptr);
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
      const int pixel = layer*(NumCDetPaddles/NumLayers) + bar*NumPaddles + localPixel;
      const TString path = TString::Format(
          "bar_%03d/layer_%d/hPairedGoodLe_Pixel%04d", bar, layer + 1, pixel);
      TH1D *stored = dynamic_cast<TH1D*>(input.Get(path));
      if (!stored) {
        std::cerr << "[CDet paired LE viewer] ERROR: missing " << path
                  << " in " << inputFile << ".\n";
        return;
      }
      histograms[layer*NumPaddles + localPixel] = static_cast<TH1D*>(stored->Clone(
          TString::Format("hPairedGoodLeView_Pixel%04d_%lu", pixel, invocation)));
      histograms[layer*NumPaddles + localPixel]->SetDirectory(nullptr);
    }
  }
  input.Close();

  TCanvas *canvas = new TCanvas(
      TString::Format("cPairedLeBar%03d_view_%lu", bar, invocation),
      TString::Format("Paired-hit LE spectra for bar %d, Layers 1 and 2", bar),
      1800, 1000);
  canvas->Divide(8, 4, 0.001, 0.001);
  for (int layer = 0; layer < NumLayers; ++layer) {
    for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
      const int index = layer*NumPaddles + localPixel;
      canvas->cd(index + 1);
      histograms[index]->Draw();
    }
  }
  canvas->Modified();
  canvas->Update();
  std::cout << "[CDet paired LE viewer] Displaying bar " << bar
            << " from " << inputFile << ".\n";
}

void writeAllCDetBarECalTimingDiagnostics(TString outputFile, TString summaryFile)
{
  const int barsPerLayer = NumCDetPaddles/(NumLayers*NumPaddles);
  if ((int)gCDetBarECalTimingSpectra.size() != NumLayers*barsPerLayer) {
    std::cerr << "[CDet bar/ECal writer] ERROR: per-bar timing spectra are not available. "
              << "Run plotCDetLayersTimeComp first.\n";
    return;
  }
  if (outputFile.IsNull() || summaryFile.IsNull()) {
    std::cerr << "[CDet bar/ECal writer] ERROR: output filenames must not be empty.\n";
    return;
  }
  TFile output(outputFile, "RECREATE");
  std::ofstream summary(summaryFile.Data());
  if (output.IsZombie() || !summary) {
    std::cerr << "[CDet bar/ECal writer] ERROR: cannot create output files.\n";
    return;
  }
  summary << "# Layer-1 CDet timing versus ECal time by bar\n"
          << "# Same accepted pairs and cuts as plotCDetLayersTimeComp\n"
          << "# bar entries profile_entries fit_status p0_ns p0_err_ns p1_ns_per_ns p1_err chi2 ndf chi2_ndf probability\n";

  int successfulFits = 0;
  for (int bar = 0; bar < barsPerLayer; ++bar) {
    TH2D *histogram = gCDetBarECalTimingSpectra[bar];
    if (!histogram) continue;
    TDirectory *directory = output.mkdir(TString::Format("bar_%03d", bar));
    directory->cd();
    TCanvas *canvas = new TCanvas(
        TString::Format("cCDetBar%03dTvsECalT", bar),
        TString::Format("CDet Layer 1 bar %d t vs ECal Time", bar), 900, 700);
    histogram->SetMinimum(40);
    histogram->GetYaxis()->SetRangeUser(20, 40);
    histogram->Draw("COLZ");
    TProfile *profile = histogram->ProfileX(TString::Format("pCDetBar%03dTvsECalT", bar));
    profile->SetMarkerStyle(20);
    profile->SetMarkerSize(0.7);
    profile->GetYaxis()->SetRangeUser(20, 40);
    TF1 *fit = new TF1(
        TString::Format("fCDetBar%03dTvsECalT", bar), "pol1",
        histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax());
    int fitStatus = -1;
    if (profile->GetEntries() >= 3) fitStatus = profile->Fit(fit, "QRS");
    profile->Draw("SAME");
    if (fitStatus == 0) fit->Draw("SAME");

    const double p0 = fitStatus == 0 ? fit->GetParameter(0) : std::numeric_limits<double>::quiet_NaN();
    const double e0 = fitStatus == 0 ? fit->GetParError(0) : std::numeric_limits<double>::quiet_NaN();
    const double p1 = fitStatus == 0 ? fit->GetParameter(1) : std::numeric_limits<double>::quiet_NaN();
    const double e1 = fitStatus == 0 ? fit->GetParError(1) : std::numeric_limits<double>::quiet_NaN();
    const double chi2 = fitStatus == 0 ? fit->GetChisquare() : std::numeric_limits<double>::quiet_NaN();
    const int ndf = fitStatus == 0 ? fit->GetNDF() : 0;
    const double chi2Ndf = ndf > 0 ? chi2/ndf : std::numeric_limits<double>::quiet_NaN();
    const double probability = ndf > 0 ? TMath::Prob(chi2, ndf) : std::numeric_limits<double>::quiet_NaN();

    TPaveText *annotation = new TPaveText(0.13, 0.77, 0.57, 0.92, "NDC");
    annotation->SetFillColor(0);
    annotation->SetTextAlign(12);
    annotation->AddText(TString::Format("Layer 1 bar %d", bar));
    if (fitStatus == 0) {
      annotation->AddText(Form("p1 = %.5f #pm %.5f", p1, e1));
      annotation->AddText(Form("#chi^{2}/NDF = %.1f/%d = %.2f", chi2, ndf, chi2Ndf));
      ++successfulFits;
    } else {
      annotation->AddText("Insufficient data or fit failed");
    }
    annotation->Draw("SAME");
    canvas->Modified();
    canvas->Update();

    histogram->Write(TString::Format("hCDetL1HalfBar%03dTvsECalT", bar));
    profile->Write();
    fit->Write();
    canvas->Write("timing_vs_ecal");
    summary << bar << " " << histogram->GetEntries() << " " << profile->GetEntries()
            << " " << fitStatus << " " << p0 << " " << e0 << " " << p1 << " " << e1
            << " " << chi2 << " " << ndf << " " << chi2Ndf << " " << probability << "\n";
    delete canvas;
    output.cd();
  }
  output.Close();
  summary.close();
  std::cout << "[CDet bar/ECal writer] Wrote " << barsPerLayer
            << " Layer-1 bar timing canvases to " << outputFile
            << " and fit results to " << summaryFile
            << " (" << successfulFits << " successful fits).\n";
}

void showCDetBarECalTiming(int bar, TString inputFile)
{
  const int barsPerLayer = NumCDetPaddles/(NumLayers*NumPaddles);
  if (bar < 0 || bar >= barsPerLayer) {
    std::cerr << "[CDet bar/ECal viewer] ERROR: bar must be in [0, "
              << barsPerLayer - 1 << "].\n";
    return;
  }
  TFile input(inputFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "[CDet bar/ECal viewer] ERROR: cannot open " << inputFile << ".\n";
    return;
  }
  const TString path = TString::Format("bar_%03d/hCDetL1HalfBar%03dTvsECalT", bar, bar);
  TH2D *stored = dynamic_cast<TH2D*>(input.Get(path));
  if (!stored) {
    std::cerr << "[CDet bar/ECal viewer] ERROR: missing " << path
              << " in " << inputFile << ".\n";
    return;
  }
  static unsigned long invocation = 0;
  ++invocation;
  TH2D *histogram = static_cast<TH2D*>(stored->Clone(
      TString::Format("hCDetBar%03dTvsECalT_view_%lu", bar, invocation)));
  histogram->SetDirectory(nullptr);
  input.Close();

  TCanvas *canvas = new TCanvas(
      TString::Format("cCDetBar%03dTvsECalT_view_%lu", bar, invocation),
      TString::Format("CDet Layer 1 bar %d t vs ECal Time", bar), 900, 700);
  histogram->SetMinimum(40);
  histogram->GetYaxis()->SetRangeUser(20, 40);
  histogram->Draw("COLZ");
  TProfile *profile = histogram->ProfileX(
      TString::Format("pCDetBar%03dTvsECalT_view_%lu", bar, invocation));
  profile->SetMarkerStyle(20);
  profile->SetMarkerSize(0.7);
  profile->GetYaxis()->SetRangeUser(20, 40);
  TF1 *fit = new TF1(
      TString::Format("fCDetBar%03dTvsECalT_view_%lu", bar, invocation), "pol1",
      histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax());
  int fitStatus = -1;
  if (profile->GetEntries() >= 3) fitStatus = profile->Fit(fit, "QRS");
  profile->Draw("SAME");
  if (fitStatus == 0) fit->Draw("SAME");
  TPaveText *annotation = new TPaveText(0.13, 0.77, 0.57, 0.92, "NDC");
  annotation->SetFillColor(0);
  annotation->SetTextAlign(12);
  annotation->AddText(TString::Format("Layer 1 bar %d", bar));
  if (fitStatus == 0) {
    const double chi2 = fit->GetChisquare();
    const int ndf = fit->GetNDF();
    annotation->AddText(Form("p1 = %.5f #pm %.5f", fit->GetParameter(1), fit->GetParError(1)));
    annotation->AddText(Form("#chi^{2}/NDF = %.1f/%d = %.2f", chi2, ndf, ndf > 0 ? chi2/ndf : 0.0));
  } else {
    annotation->AddText("Insufficient data or fit failed");
  }
  annotation->Draw("SAME");
  canvas->Modified();
  canvas->Update();
  std::cout << "[CDet bar/ECal viewer] Displaying Layer-1 bar " << bar
            << " from " << inputFile << ".\n";
}

void plotAllCDetBarECalTiming(TString outputFile, TString summaryFile)
{
  const int halfModulesPerLayer = NumModules*NumSides;
  const int halfBarsPerLayer = halfModulesPerLayer*NumBars;
  if ((int)gCDetBarECalTimingSpectra.size() != NumLayers*halfBarsPerLayer) {
    std::cerr << "[CDet half-bar/ECal overview] ERROR: half-bar spectra are not available. "
              << "Run plotCDetLayersTimeComp first.\n";
    return;
  }
  TFile output(outputFile, "RECREATE");
  std::ofstream summary(summaryFile.Data());
  if (output.IsZombie() || !summary) {
    std::cerr << "[CDet half-bar/ECal overview] ERROR: cannot create output files.\n";
    return;
  }
  summary << "# CDet timing versus ECal time by half-bar (PMT group)\n"
          << "# Four overview canvases: two seven-half-bar banks for each layer\n"
          << "# layer half_module bank half_bar_in_half_module global_half_bar entries fit_status p0_ns p0_err_ns p1_ns_per_ns p1_err chi2 ndf chi2_ndf probability\n";

  int successfulFits = 0;
  for (int layer = 0; layer < NumLayers; ++layer) {
    TDirectory *layerDirectory = output.mkdir(TString::Format("layer_%d", layer + 1));
    for (int bank = 0; bank < 2; ++bank) {
      TCanvas *canvas = new TCanvas(
          TString::Format("cCDetLayer%dHalfBarBank%dECalTiming", layer + 1, bank),
          TString::Format("CDet Layer %d half-bars %d-%d versus ECal time", layer + 1,
                          bank*NumHalfBarsPerBank, bank*NumHalfBarsPerBank + NumHalfBarsPerBank - 1),
          2100, 1500);
      canvas->Divide(NumHalfBarsPerBank, halfModulesPerLayer, 0.001, 0.001);
      for (int halfModule = 0; halfModule < halfModulesPerLayer; ++halfModule) {
        TDirectory *halfModuleDirectory = layerDirectory->GetDirectory(
            TString::Format("half_module_%d", halfModule));
        if (!halfModuleDirectory)
          halfModuleDirectory = layerDirectory->mkdir(TString::Format("half_module_%d", halfModule));
        for (int column = 0; column < NumHalfBarsPerBank; ++column) {
          const int halfBarInHalfModule = bank*NumHalfBarsPerBank + column;
          const int globalHalfBar = halfModule*NumBars + halfBarInHalfModule;
          const int index = layer*halfBarsPerLayer + globalHalfBar;
          TH2D *histogram = gCDetBarECalTimingSpectra[index];
          canvas->cd(halfModule*NumHalfBarsPerBank + column + 1);
        histogram->SetMinimum(1);
        histogram->GetYaxis()->SetRangeUser(20, 40);
        histogram->Draw("COLZ");
        TProfile *profile = histogram->ProfileX(TString::Format(
              "pCDetL%dHalfBar%03dTvsECalT", layer + 1, globalHalfBar));
        profile->SetMarkerStyle(20);
        profile->SetMarkerSize(0.45);
        profile->GetYaxis()->SetRangeUser(20, 40);
        TF1 *fit = new TF1(TString::Format(
              "fCDetL%dHalfBar%03dTvsECalT", layer + 1, globalHalfBar),
            "pol1", histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax());
        int fitStatus = -1;
        if (profile->GetEntries() >= 3) fitStatus = profile->Fit(fit, "QRS");
        profile->Draw("SAME");
        if (fitStatus == 0) fit->Draw("SAME");

        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double p0 = fitStatus == 0 ? fit->GetParameter(0) : nan;
        const double e0 = fitStatus == 0 ? fit->GetParError(0) : nan;
        const double p1 = fitStatus == 0 ? fit->GetParameter(1) : nan;
        const double e1 = fitStatus == 0 ? fit->GetParError(1) : nan;
        const double chi2 = fitStatus == 0 ? fit->GetChisquare() : nan;
        const int ndf = fitStatus == 0 ? fit->GetNDF() : 0;
        const double chi2Ndf = ndf > 0 ? chi2/ndf : nan;
        const double probability = ndf > 0 ? TMath::Prob(chi2, ndf) : nan;

        TPaveText *annotation = new TPaveText(0.10, 0.72, 0.74, 0.91, "NDC");
        annotation->SetFillColor(0);
        annotation->SetTextAlign(12);
        annotation->SetTextSize(0.050);
          annotation->AddText(Form("HM %d, half-bar %d", halfModule, halfBarInHalfModule));
        if (fitStatus == 0) {
          annotation->AddText(Form("p1=%.3f#pm%.3f", p1, e1));
          annotation->AddText(Form("#chi^{2}/NDF=%.1f", chi2Ndf));
          ++successfulFits;
        } else {
          annotation->AddText("fit unavailable");
        }
        annotation->Draw("SAME");

        halfModuleDirectory->cd();
        histogram->Write();
        profile->Write();
        fit->Write();
        layerDirectory->cd();
          summary << layer + 1 << " " << halfModule << " " << bank << " "
                << halfBarInHalfModule << " " << globalHalfBar
                << " " << histogram->GetEntries() << " " << fitStatus
                << " " << p0 << " " << e0 << " " << p1 << " " << e1
                << " " << chi2 << " " << ndf << " " << chi2Ndf
                << " " << probability << "\n";
        }
      }
      canvas->Modified();
      canvas->Update();
      output.cd();
      canvas->Write(TString::Format("layer_%d_bank_%d_overview", layer + 1, bank));
    }
  }
  output.Close();
  summary.close();
  std::cout << "[CDet half-bar/ECal overview] Displayed and wrote four 7x6 canvases to "
            << outputFile << ", with fit results in " << summaryFile
            << " (" << successfulFits << " successful fits of 168).\n";
}

void calibrateCDetHalfBarIntercepts(bool overwrite, double referenceECalTime,
                                    int minEntries, double maxInterceptError,
                                    TString outputRoot, TString outputSummary)
{
  const int halfBarsPerLayer = NumCDetPaddles/(NumLayers*NumPaddles);
  const int totalHalfBars = NumLayers*halfBarsPerLayer;
  if ((int)gCDetHalfBarECalTimingSamples.size() != totalHalfBars) {
    std::cerr << "[CDet half-bar intercepts] ERROR: unbinned samples are unavailable. "
              << "Run plotCDetLayersTimeComp first.\n";
    return;
  }
  if (minEntries < 3 || maxInterceptError <= 0.0) {
    std::cerr << "[CDet half-bar intercepts] ERROR: invalid quality thresholds.\n";
    return;
  }

  struct GroupResult {
    int entries = 0;
    double meanX = NAN, meanY = NAN, intercept = NAN, error = NAN;
    double correction = 0.0;
    bool valid = false;
  };
  std::vector<GroupResult> groups(totalHalfBars);

  // Estimate the common within-half-bar slope directly from unbinned hits.
  double withinXY = 0.0, withinXX = 0.0;
  for (int index = 0; index < totalHalfBars; ++index) {
    const auto &samples = gCDetHalfBarECalTimingSamples[index];
    if ((int)samples.size() < minEntries) continue;
    double sx = 0.0, sy = 0.0;
    for (const auto &sample : samples) { sx += sample.first; sy += sample.second; }
    const double mx = sx/samples.size(), my = sy/samples.size();
    for (const auto &sample : samples) {
      const double dx = sample.first - mx;
      withinXY += dx*(sample.second - my);
      withinXX += dx*dx;
    }
  }
  if (withinXX <= 0.0) {
    std::cerr << "[CDet half-bar intercepts] ERROR: insufficient within-group variation.\n";
    return;
  }
  const double commonSlope = withinXY/withinXX;

  double targetNumerator = 0.0, targetDenominator = 0.0;
  for (int index = 0; index < totalHalfBars; ++index) {
    const auto &samples = gCDetHalfBarECalTimingSamples[index];
    GroupResult &group = groups[index];
    group.entries = samples.size();
    if (samples.empty()) continue;
    double sx = 0.0, sy = 0.0;
    for (const auto &sample : samples) { sx += sample.first; sy += sample.second; }
    group.meanX = sx/samples.size();
    group.meanY = sy/samples.size();
    group.intercept = group.meanY + commonSlope*(referenceECalTime - group.meanX);
    if (samples.size() >= 2) {
      double sse = 0.0;
      for (const auto &sample : samples) {
        const double predicted = group.intercept + commonSlope*(sample.first - referenceECalTime);
        const double residual = sample.second - predicted;
        sse += residual*residual;
      }
      group.error = std::sqrt((sse/(samples.size() - 1))/samples.size());
    }
    group.valid = group.entries >= minEntries && std::isfinite(group.error) &&
                  group.error <= maxInterceptError;
    if (group.valid) {
      targetNumerator += group.entries*group.intercept;
      targetDenominator += group.entries;
    }
  }
  if (targetDenominator <= 0.0) {
    std::cerr << "[CDet half-bar intercepts] ERROR: no half-bars passed validation.\n";
    return;
  }
  const double detectorReference = targetNumerator/targetDenominator;
  double betweenVariance = 0.0, withinResidualSquares = 0.0;
  int validGroups = 0;
  for (int index = 0; index < totalHalfBars; ++index) {
    GroupResult &group = groups[index];
    if (!group.valid) continue;
    group.correction = detectorReference - group.intercept;
    betweenVariance += group.entries*std::pow(group.intercept - detectorReference, 2);
    for (const auto &sample : gCDetHalfBarECalTimingSamples[index]) {
      const double predicted = group.intercept + commonSlope*(sample.first - referenceECalTime);
      withinResidualSquares += std::pow(sample.second - predicted, 2);
    }
    ++validGroups;
  }
  const double betweenRms = std::sqrt(betweenVariance/targetDenominator);
  const double withinRms = std::sqrt(withinResidualSquares/targetDenominator);
  const double totalRmsBefore = std::sqrt(
      (withinResidualSquares + betweenVariance)/targetDenominator);
  const double predictedTotalRmsAfter = withinRms;

  std::ofstream summary(outputSummary.Data());
  TFile rootOutput(outputRoot, "RECREATE");
  if (!summary || rootOutput.IsZombie()) {
    std::cerr << "[CDet half-bar intercepts] ERROR: cannot create diagnostic outputs.\n";
    return;
  }
  summary << "# CDet half-bar intercept alignment at common ECal time\n"
          << "# reference_ecal_time_ns " << referenceECalTime << "\n"
          << "# common_within_slope_ns_per_ns " << commonSlope << "\n"
          << "# detector_reference_ns " << detectorReference << "\n"
          << "# layer half_bar entries mean_ecal_ns mean_cdet_ns projected_cdet_at_reference_ns intercept_error_ns valid proposed_offset_increment_ns\n";

  rootOutput.cd();
  std::vector<double> xValid, yValid, eyValid, correctionValid;
  xValid.reserve(validGroups); yValid.reserve(validGroups);
  eyValid.reserve(validGroups); correctionValid.reserve(validGroups);
  TH1D *hIntercepts = new TH1D("hCDetHalfBarIntercepts",
      "Half-bar times at common ECal reference;CDet time at reference (ns);Half-bars", 120, 15, 45);
  TH1D *hCorrections = new TH1D("hCDetHalfBarCorrections",
      "Proposed half-bar offset increments;Offset increment (ns);Half-bars", 120, -15, 15);
  for (int index = 0; index < totalHalfBars; ++index) {
    const int layer = index/halfBarsPerLayer;
    const int halfBar = index%halfBarsPerLayer;
    const GroupResult &group = groups[index];
    if (group.valid) {
      xValid.push_back(index);
      yValid.push_back(group.intercept);
      eyValid.push_back(group.error);
      correctionValid.push_back(group.correction);
      hIntercepts->Fill(group.intercept);
      hCorrections->Fill(group.correction);
    }
    summary << layer + 1 << " " << halfBar << " " << group.entries << " "
            << group.meanX << " " << group.meanY << " " << group.intercept << " "
            << group.error << " " << group.valid << " " << group.correction << "\n";
  }
  TGraphErrors *gIntercepts = new TGraphErrors(validGroups, xValid.data(), yValid.data(),
                                                nullptr, eyValid.data());
  gIntercepts->SetName("gCDetHalfBarInterceptsAtReference");
  gIntercepts->SetTitle(TString::Format(
      "Half-bar CDet time at ECal reference %.1f ns;Layer/half-bar index;Projected CDet time (ns)",
      referenceECalTime));
  TGraph *gCorrections = new TGraph(validGroups, xValid.data(), correctionValid.data());
  gCorrections->SetName("gCDetHalfBarProposedOffsetIncrements");
  gCorrections->SetTitle("Proposed half-bar offset increments;Layer/half-bar index;Offset increment (ns)");
  TCanvas *canvas = new TCanvas("cCDetHalfBarInterceptAlignment",
      "CDet half-bar intercept alignment", 1200, 900);
  canvas->Divide(2, 2);
  canvas->cd(1); gIntercepts->SetMarkerStyle(20); gIntercepts->SetMarkerSize(0.45); gIntercepts->Draw("AP");
  canvas->cd(2); hIntercepts->Draw();
  canvas->cd(3); gCorrections->SetMarkerStyle(20); gCorrections->SetMarkerSize(0.45); gCorrections->Draw("AP");
  canvas->cd(4); hCorrections->Draw();
  rootOutput.Write();
  rootOutput.Close();
  summary.close();

  std::cout << "\n[CDet half-bar intercept alignment]\n"
            << "  common unbinned within-half-bar slope: " << commonSlope << " ns/ns\n"
            << "  ECal reference time: " << referenceECalTime << " ns\n"
            << "  validated half-bars: " << validGroups << " / " << totalHalfBars << "\n"
            << "  detector reference at common ECal time: " << detectorReference << " ns\n"
            << "  within-half-bar residual RMS: " << withinRms << " ns\n"
            << "  between-half-bar intercept RMS: " << betweenRms << " ns\n"
            << "  total timing RMS before proposed alignment: " << totalRmsBefore << " ns\n"
            << "  predicted total timing RMS after proposed alignment: "
            << predictedTotalRmsAfter << " ns\n"
            << "  diagnostics: " << outputRoot << "\n"
            << "  proposed corrections: " << outputSummary << "\n";

  if (!overwrite) {
    std::cout << "  diagnostic only: calibration constants were not changed.\n";
    return;
  }
  if (!gPixelToffsetLoaded || (int)gPixelToffsetCorr.size() != NumCDetPaddles) {
    std::cerr << "[CDet half-bar intercepts] ERROR: active pixel offsets are unavailable; nothing written.\n";
    return;
  }
  const int pixelsPerLayer = NumCDetPaddles/NumLayers;
  for (int index = 0; index < totalHalfBars; ++index) {
    if (!groups[index].valid) continue;
    const int layer = index/halfBarsPerLayer;
    const int halfBar = index%halfBarsPerLayer;
    const int pixelBase = layer*pixelsPerLayer + halfBar*NumPaddles;
    for (int localPixel = 0; localPixel < NumPaddles; ++localPixel)
      gPixelToffsetCorr[pixelBase + localPixel] += groups[index].correction;
  }
  gLastCalibrationFitSucceeded = WriteCalibrationConstants(gCalibrationFile);
  std::cout << "  applied validated half-bar corrections to " << gCalibrationFile << ".\n";
}

void extractCDetYPositionCalibration(TString outputFile, int minEntries)
{
  const int totalHalfBars = NumCDetPaddles/NumPaddles;
  if ((int)gCDetHalfBarECalYSamples.size() != totalHalfBars || minEntries < 3) {
    std::cerr << "[CDet y calibration] ERROR: samples are unavailable or minEntries is invalid. "
              << "Run plotCDetLayersTimeComp first.\n";
    return;
  }
  struct GroupFit { double xy = 0.0, xx = 0.0, slope = NAN, sse = 0.0, error = NAN; long long n = 0; int h = 0; };
  GroupFit fits[NumLayers][NumSides][3];
  std::vector<double> meanY(totalHalfBars, NAN), meanT(totalHalfBars, NAN);
  std::vector<int> valid(totalHalfBars, 0);
  auto decode = [&](int halfBar, int &layer, int &side, int &group) {
    const int pixel = halfBar*NumPaddles;
    const int pixelsPerLayer = NumCDetPaddles/NumLayers;
    layer = pixel/pixelsPerLayer;
    const int withinLayer = pixel%pixelsPerLayer;
    side = withinLayer/NumCDetPaddlesPerSide;
    const int withinSide = withinLayer%NumCDetPaddlesPerSide;
    const int module = withinSide/(NumBars*NumPaddles);
    const int pmt = (withinSide%(NumBars*NumPaddles))/NumPaddles;
    const int section = 2*module + pmt/NumHalfBarsPerBank;
    group = std::min(section, 5 - section);
  };
  for (int halfBar = 0; halfBar < totalHalfBars; ++halfBar) {
    const auto &samples = gCDetHalfBarECalYSamples[halfBar];
    if ((int)samples.size() < minEntries) continue;
    double sy = 0.0, st = 0.0;
    for (const auto &sample : samples) { sy += sample.first; st += sample.second; }
    meanY[halfBar] = sy/samples.size(); meanT[halfBar] = st/samples.size(); valid[halfBar] = 1;
    int layer, side, group; decode(halfBar, layer, side, group);
    GroupFit &fit = fits[layer][side][group]; ++fit.h; fit.n += samples.size();
    for (const auto &sample : samples) {
      const double dy = sample.first - meanY[halfBar], dt = sample.second - meanT[halfBar];
      fit.xy += dy*dt; fit.xx += dy*dy;
    }
  }
  for (int layer = 0; layer < NumLayers; ++layer)
    for (int side = 0; side < NumSides; ++side)
      for (int group = 0; group < 3; ++group)
        if (fits[layer][side][group].xx > 0.0)
          fits[layer][side][group].slope = fits[layer][side][group].xy/fits[layer][side][group].xx;
  for (int halfBar = 0; halfBar < totalHalfBars; ++halfBar) {
    if (!valid[halfBar]) continue;
    int layer, side, group; decode(halfBar, layer, side, group);
    GroupFit &fit = fits[layer][side][group];
    for (const auto &sample : gCDetHalfBarECalYSamples[halfBar]) {
      const double residual = (sample.second - meanT[halfBar]) - fit.slope*(sample.first - meanY[halfBar]);
      fit.sse += residual*residual;
    }
  }
  for (int layer = 0; layer < NumLayers; ++layer)
    for (int side = 0; side < NumSides; ++side)
      for (int group = 0; group < 3; ++group) {
        GroupFit &fit = fits[layer][side][group];
        const long long ndf = fit.n - fit.h - 1;
        if (ndf > 0 && fit.xx > 0.0) fit.error = std::sqrt((fit.sse/ndf)/fit.xx);
      }
  std::ofstream output(outputFile.Data());
  if (!output) { std::cerr << "[CDet y calibration] ERROR: cannot write " << outputFile << ".\n"; return; }
  output << "# CDet y-position calibration trained against ECal y\n"
         << "# This file is independent of CDet_calibration_dt.dat and does not alter timing constants.\n"
         << "CDetYPosition.MinEntries: " << minEntries << "\n";
  const char *sideKey[NumSides] = {"Left", "Right"};
  for (int layer = 0; layer < NumLayers; ++layer)
    for (int side = 0; side < NumSides; ++side)
      for (int group = 0; group < 3; ++group) {
        const GroupFit &fit = fits[layer][side][group];
        output << "CDetYPosition.Slope.L" << layer + 1 << "." << sideKey[side] << ".G" << group + 1
               << ": " << fit.slope << "\n";
        output << "CDetYPosition.SlopeError.L" << layer + 1 << "." << sideKey[side] << ".G" << group + 1
               << ": " << fit.error << "\n";
      }
  for (int halfBar = 0; halfBar < totalHalfBars; ++halfBar) {
    double intercept = NAN;
    if (valid[halfBar]) {
      int layer, side, group; decode(halfBar, layer, side, group);
      intercept = meanT[halfBar] - fits[layer][side][group].slope*meanY[halfBar];
    }
    output << "CDetYPosition.Valid.HalfBar" << halfBar << ": " << valid[halfBar] << "\n";
    output << "CDetYPosition.Intercept.HalfBar" << halfBar << ": " << intercept << "\n";
  }
  output.close();
  std::cout << "[CDet y calibration] Wrote 12 propagation slopes and "
            << std::count(valid.begin(), valid.end(), 1) << " valid half-bar intercepts to " << outputFile << ".\n"
            << "  CDet timing calibration constants were not changed.\n";
}

void plotCDetYPositionResolution(TString calibrationFile)
{
  const int totalHalfBars = NumCDetPaddles/NumPaddles;
  if ((int)gCDetHalfBarECalYSamples.size() != totalHalfBars) {
    std::cerr << "[CDet y validation] ERROR: run plotCDetLayersTimeComp first.\n"; return;
  }
  TEnv env;
  if (env.ReadFile(calibrationFile, kEnvLocal) < 0) {
    std::cerr << "[CDet y validation] ERROR: cannot read " << calibrationFile << ".\n"; return;
  }
  if (env.GetValue("CDetYPosition.MinEntries", -1) < 0) {
    std::cerr << "[CDet y validation] ERROR: " << calibrationFile
              << " does not contain readable TEnv keys. Regenerate it with "
              << "extractCDetYPositionCalibration().\n";
    return;
  }
  TH2D *correlation[NumLayers] = {
      new TH2D("hCDetL1YRecoVsECalY", "Layer 1 CDet y reconstruction;ECal y (m);CDet reconstructed y (m)", 120,-0.6,0.6,120,-0.6,0.6),
      new TH2D("hCDetL2YRecoVsECalY", "Layer 2 CDet y reconstruction;ECal y (m);CDet reconstructed y (m)", 120,-0.6,0.6,120,-0.6,0.6)};
  TH1D *residual[NumLayers] = {
      new TH1D("hCDetL1YResidual", "Layer 1 CDet y resolution;CDet y - ECal y (m);Hits", 160,-0.4,0.4),
      new TH1D("hCDetL2YResidual", "Layer 2 CDet y resolution;CDet y - ECal y (m);Hits", 160,-0.4,0.4)};
  const char *sideKey[NumSides] = {"Left", "Right"};
  for (int halfBar = 0; halfBar < totalHalfBars; ++halfBar) {
    if (!env.GetValue(TString::Format("CDetYPosition.Valid.HalfBar%d", halfBar), 0)) continue;
    const int pixel = halfBar*NumPaddles, pixelsPerLayer = NumCDetPaddles/NumLayers;
    const int layer = pixel/pixelsPerLayer, withinLayer = pixel%pixelsPerLayer;
    const int side = withinLayer/NumCDetPaddlesPerSide;
    const int withinSide = withinLayer%NumCDetPaddlesPerSide;
    const int module = withinSide/(NumBars*NumPaddles);
    const int pmt = (withinSide%(NumBars*NumPaddles))/NumPaddles;
    const int section = 2*module + pmt/NumHalfBarsPerBank;
    const int group = std::min(section, 5 - section);
    const double slope = env.GetValue(TString::Format("CDetYPosition.Slope.L%d.%s.G%d", layer+1,sideKey[side],group+1), NAN);
    const double intercept = env.GetValue(TString::Format("CDetYPosition.Intercept.HalfBar%d", halfBar), NAN);
    if (!std::isfinite(slope) || std::abs(slope) < 1e-9 || !std::isfinite(intercept)) continue;
    for (const auto &sample : gCDetHalfBarECalYSamples[halfBar]) {
      const double yReco = (sample.second - intercept)/slope;
      correlation[layer]->Fill(sample.first, yReco);
      residual[layer]->Fill(yReco - sample.first);
    }
  }
  TCanvas *canvas = new TCanvas("cCDetYPositionResolution", "CDet y-position validation", 1200, 900);
  canvas->Divide(2,2);
  for (int layer = 0; layer < NumLayers; ++layer) {
    canvas->cd(layer*2 + 1); correlation[layer]->Draw("COLZ");
    TLine *diagonal = new TLine(-0.6,-0.6,0.6,0.6); diagonal->SetLineColor(kRed); diagonal->Draw("SAME");
    canvas->cd(layer*2 + 2); residual[layer]->Draw();
    std::cout << "[CDet y validation] Layer " << layer + 1 << ": entries=" << residual[layer]->GetEntries()
              << ", mean residual=" << residual[layer]->GetMean() << " m, RMS=" << residual[layer]->GetStdDev() << " m\n";
  }
  canvas->Update();

  TH2D *layerCorrelation = new TH2D(
      "hCDetL2YVsL1Y", "CDet layer-to-layer y reconstruction;Layer 1 reconstructed y (m);Layer 2 reconstructed y (m)",
      120, -0.6, 0.6, 120, -0.6, 0.6);
  TH1D *layerDifference = new TH1D(
      "hCDetYLayerDifference", "CDet layer-to-layer y difference;Layer 1 y - Layer 2 y (m);Accepted pairs",
      200, -0.8, 0.8);
  auto reconstructY = [&](int halfBar, double time, double &yReco) {
    if (halfBar < 0 || halfBar >= totalHalfBars ||
        !env.GetValue(TString::Format("CDetYPosition.Valid.HalfBar%d", halfBar), 0)) return false;
    const int pixel = halfBar*NumPaddles, pixelsPerLayer = NumCDetPaddles/NumLayers;
    const int layer = pixel/pixelsPerLayer, withinLayer = pixel%pixelsPerLayer;
    const int side = withinLayer/NumCDetPaddlesPerSide;
    const int withinSide = withinLayer%NumCDetPaddlesPerSide;
    const int module = withinSide/(NumBars*NumPaddles);
    const int pmt = (withinSide%(NumBars*NumPaddles))/NumPaddles;
    const int section = 2*module + pmt/NumHalfBarsPerBank;
    const int group = std::min(section, 5 - section);
    const double slope = env.GetValue(
        TString::Format("CDetYPosition.Slope.L%d.%s.G%d", layer+1,sideKey[side],group+1), NAN);
    const double intercept = env.GetValue(
        TString::Format("CDetYPosition.Intercept.HalfBar%d", halfBar), NAN);
    if (!std::isfinite(slope) || std::abs(slope) < 1e-9 || !std::isfinite(intercept)) return false;
    yReco = (time - intercept)/slope;
    return std::isfinite(yReco);
  };
  for (const auto &sample : gCDetYPairedSamples) {
    double y1 = NAN, y2 = NAN;
    if (!reconstructY(sample.halfBar1, sample.t1, y1) ||
        !reconstructY(sample.halfBar2, sample.t2, y2)) continue;
    layerCorrelation->Fill(y1, y2);
    layerDifference->Fill(y1 - y2);
  }
  TCanvas *layerCanvas = new TCanvas(
      "cCDetYLayerComparison", "CDet layer-to-layer y-position validation", 1200, 600);
  layerCanvas->Divide(2,1);
  layerCanvas->cd(1); layerCorrelation->Draw("COLZ");
  TLine *layerDiagonal = new TLine(-0.6,-0.6,0.6,0.6);
  layerDiagonal->SetLineColor(kRed); layerDiagonal->Draw("SAME");
  layerCanvas->cd(2); layerDifference->Draw();
  layerCanvas->Update();
  const double differenceRms = layerDifference->GetStdDev();
  std::cout << "[CDet y layer comparison] paired entries=" << layerDifference->GetEntries()
            << ", mean(L1-L2)=" << layerDifference->GetMean() << " m"
            << ", RMS(L1-L2)=" << differenceRms << " m"
            << ", equal-independent-layer estimate=" << differenceRms/std::sqrt(2.0) << " m\n";
}

void plotCDetLayersTimeComp(const char *configFile,
                            Int_t overwriteOverride = -1)
{
  TEnv env;
  const char *caller = "CDet layer comparison";
  if (!LoadCDetConfiguration(env, configFile, caller)) return;

  const Bool_t overwrite = overwriteOverride >= 0
                               ? overwriteOverride != 0
                               : env.GetValue("display.overwrite", false);
  const Int_t pixel = env.GetValue("display.pixel", 416);
  const Double_t width = env.GetValue("display.histogram_width", 1.0);
  const Double_t layerDtMin = env.GetValue("display.layer_dt_min", -15.0);
  const Double_t layerDtMax = env.GetValue("display.layer_dt_max", 15.0);
  const Double_t layerDxMin = env.GetValue("display.layer_dx_min", -0.01);
  const Double_t layerDxMax = env.GetValue("display.layer_dx_max", 0.01);
  const Double_t leMin = env.GetValue("display.le_min", 0.02);
  const Double_t leMax = env.GetValue("display.le_max", 60.0);
  const Double_t totMin = env.GetValue("display.tot_min", 0.0);
  const Double_t totMax = env.GetValue("display.tot_max", 70.0);
  const Double_t dtHistMin = env.GetValue("display.dt_hist_min", -20.0);
  const Double_t dtHistMax = env.GetValue("display.dt_hist_max", 20.0);
  const Double_t cdetTimeMin = env.GetValue("display.cdet_time_min", 0.0);
  const Double_t cdetTimeMax = env.GetValue("display.cdet_time_max", 60.0);
  const Double_t cdetTotMin = env.GetValue("display.cdet_tot_min", 0.0);
  const Double_t cdetTotMax = env.GetValue("display.cdet_tot_max", 80.0);
  const Double_t ecalTimeMin = env.GetValue("display.ecal_time_min", -40.0);
  const Double_t ecalTimeMax = env.GetValue("display.ecal_time_max", 40.0);
  const Double_t hcalTimeMin = env.GetValue("display.hcal_time_min", -10.0);
  const Double_t hcalTimeMax = env.GetValue("display.hcal_time_max", 10.0);
  const Double_t ecalCdetDtMin = env.GetValue("display.ecal_cdet_dt_min", -60.0);
  const Double_t ecalCdetDtMax = env.GetValue("display.ecal_cdet_dt_max", 30.0);
  TString hitMode = env.GetValue("display.hit_mode", "all");
  hitMode.ToLower();
  const Bool_t displayBestHits = hitMode == "best";
  const Bool_t hasBestHitPeakMean = env.Defined("display.best_hit_peak_mean");
  const Bool_t hasBestHitPeakSigma = env.Defined("display.best_hit_peak_sigma");
  const Double_t bestHitPeakMean = hasBestHitPeakMean
      ? env.GetValue("display.best_hit_peak_mean", 0.0)
      : std::numeric_limits<double>::quiet_NaN();
  const Double_t bestHitPeakSigma = hasBestHitPeakSigma
      ? env.GetValue("display.best_hit_peak_sigma", 0.0)
      : std::numeric_limits<double>::quiet_NaN();
  const Double_t bestHitNSigma = env.GetValue("display.best_hit_nsigma", 3.0);
  const Bool_t allowMultiplePairs = env.GetValue("display.allow_multiple_pairs", true);
  const Double_t xBinWidth = env.GetValue("display.x_bin_width", 0.005);
  const Double_t xMin = env.GetValue("display.x_min", -1.5);
  const Double_t xMax = env.GetValue("display.x_max", 1.5);
  const Double_t zBinWidth = env.GetValue("display.z_bin_width", 0.01);
  const Double_t zMin = env.GetValue("display.z_min", 0.0);
  const Double_t zMax = env.GetValue("display.z_max", 7.0);
  const Double_t yCorrectionRefractiveIndex =
      env.GetValue("display.y_correction_refractive_index", 0.0);

  if (pixel < 0 || pixel >= NumCDetPaddles/2 || width <= 0.0 ||
      layerDtMin >= layerDtMax || layerDxMin >= layerDxMax ||
      leMin >= leMax || totMin >= totMax || dtHistMin >= dtHistMax ||
      cdetTimeMin >= cdetTimeMax || cdetTotMin >= cdetTotMax ||
      ecalTimeMin >= ecalTimeMax || hcalTimeMin >= hcalTimeMax || ecalCdetDtMin >= ecalCdetDtMax ||
      (hitMode != "all" && hitMode != "best") ||
      hasBestHitPeakMean != hasBestHitPeakSigma ||
      (hasBestHitPeakMean && (!std::isfinite(bestHitPeakMean) ||
                              !std::isfinite(bestHitPeakSigma) || bestHitPeakSigma <= 0.0)) ||
      !std::isfinite(bestHitNSigma) || bestHitNSigma <= 0.0 ||
      xBinWidth <= 0.0 || xMin >= xMax || zBinWidth <= 0.0 || zMin >= zMax ||
      !std::isfinite(yCorrectionRefractiveIndex) || yCorrectionRefractiveIndex < 0.0 ||
      (yCorrectionRefractiveIndex > 0.0 && yCorrectionRefractiveIndex < 1.0)) {
    std::cerr << "[" << caller << "] ERROR: invalid pixel, binning, or min/max range in "
              << configFile << ".\n";
    return;
  }

  plotCDetLayersTimeComp(overwrite, pixel, width, layerDtMin, layerDtMax,
                         layerDxMin, layerDxMax, leMin, leMax, totMin, totMax,
                         dtHistMin, dtHistMax, cdetTimeMin, cdetTimeMax,
                         cdetTotMin, cdetTotMax, ecalTimeMin, ecalTimeMax,
                         ecalCdetDtMin, ecalCdetDtMax, allowMultiplePairs,
                         xBinWidth, xMin, xMax, zBinWidth, zMin, zMax,
                         hcalTimeMin, hcalTimeMax, displayBestHits,
                         bestHitPeakMean, bestHitPeakSigma, bestHitNSigma,
                         yCorrectionRefractiveIndex);
}

void BuildCDetEventDisplay(int selectedBar, double diffMinCut, double diffMaxCut,
                           double xdiffMinCut, double xdiffMaxCut,
                           double tdiffECalCDetMin, double tdiffECalCDetMax,
                           double xMin, double xMax, double zMin, double zMax,
                           bool initialBestHits, double bestHitPeakMean,
                           double bestHitPeakSigma, double bestHitNSigma)
{
  gCDetDisplayEvents.clear();
  gCDetDisplayIndex = -1;
  gCDetDisplayXMin = xMin;
  gCDetDisplayXMax = xMax;
  gCDetDisplayZMin = zMin;
  gCDetDisplayZMax = zMax;
  gCDetDisplayBestHits = false;
  gCDetDisplayBestHitTimingValid = false;
  gCDetDisplayBestHitPeakMean = std::numeric_limits<double>::quiet_NaN();
  gCDetDisplayBestHitPeakSigma = std::numeric_limits<double>::quiet_NaN();
  gCDetDisplayBestHitNSigma = bestHitNSigma;

  if (selectedBar < 0 || selectedBar >= NumCDetPaddles/(2*NumPaddles)) {
    std::cerr << "[CDet event display] ERROR: Layer-1 bar " << selectedBar
              << " is outside the valid range [0, "
              << NumCDetPaddles/(2*NumPaddles) - 1 << "].\n";
    return;
  }

  const int selectedLayer1BarBase = selectedBar * NumPaddles;
  size_t nEvents = pairs_CDet.size();
  nEvents = std::min(nEvents, vGoodID.size());
  nEvents = std::min(nEvents, vCDetGoodX.size());
  nEvents = std::min(nEvents, v_GoodECalX.size());
  nEvents = std::min(nEvents, v_GoodECalY.size());
  nEvents = std::min(nEvents, v_GoodECalE.size());
  nEvents = std::min(nEvents, v_GoodECalAdcTime.size());
  nEvents = std::min(nEvents, vTreeEntry.size());

  for (size_t ev = 0; ev < nEvents; ++ev) {
    std::vector<PairHit> selectedPairs;
    for (const auto& pair : pairs_CDet[ev]) {
      const bool selectedLayer1Bar = pair.id1 >= selectedLayer1BarBase &&
                                     pair.id1 < selectedLayer1BarBase + NumPaddles;
      if (!selectedLayer1Bar) continue;
      if (pair.dt < diffMinCut || pair.dt > diffMaxCut ||
          pair.dx < xdiffMinCut || pair.dx > xdiffMaxCut) continue;

      const double pairTime = 0.5 * (pair.t1 + pair.t2);
      const double ecalCDetDt = v_GoodECalAdcTime[ev] - pairTime;
      if (ecalCDetDt < tdiffECalCDetMin || ecalCDetDt > tdiffECalCDetMax) continue;
      selectedPairs.push_back(pair);
    }
    if (selectedPairs.empty()) continue;

    CDetDisplayEvent displayEvent;
    displayEvent.savedEventIndex = ev;
    displayEvent.treeEntry = vTreeEntry[ev];
    displayEvent.runNumber = gRunNumber;
    displayEvent.selectedBar = selectedBar;
    displayEvent.ecalX = v_GoodECalX[ev];
    displayEvent.ecalY = v_GoodECalY[ev];
    displayEvent.ecalZ = ECal_dist;
    displayEvent.ecalEnergy = v_GoodECalE[ev];
    displayEvent.ecalTime = v_GoodECalAdcTime[ev];
    displayEvent.selectedPairs = selectedPairs;

    const size_t nHits = std::min(
        std::min(vGoodID[ev].size(), vGoodLe[ev].size()),
        std::min(vGoodTot[ev].size(),
                 std::min(vCDetGoodX[ev].size(),
                          std::min(vCDetGoodY[ev].size(), vCDetGoodZ[ev].size()))));
    displayEvent.hits.reserve(nHits);
    for (size_t ihit = 0; ihit < nHits; ++ihit) {
      const int id = vGoodID[ev][ihit];
      if (id < 0 || id >= NumCDetPaddles) continue;
      const int detectorColumn = id / NumCDetPaddlesPerSide;
      const int withinSide = id % NumCDetPaddlesPerSide;
      const int pixelsPerModule = NumBars * NumPaddles;
      CDetDisplayHit hit;
      hit.id = id;
      hit.layer = detectorColumn / NumSides;
      hit.side = detectorColumn % NumSides;
      hit.module = withinSide / pixelsPerModule;
      hit.bar = (withinSide % pixelsPerModule) / NumPaddles;
      hit.pixel = withinSide % NumPaddles;
      // vCDetGoodX already contains the layer-specific corrected coordinate.
      hit.x = vCDetGoodX[ev][ihit];
      hit.y = vCDetGoodY[ev][ihit];
      hit.z = vCDetGoodZ[ev][ihit];
      hit.le = vGoodLe[ev][ihit];
      hit.te = ihit < vGoodTe[ev].size() ? vGoodTe[ev][ihit] : 0.0;
      hit.tot = vGoodTot[ev][ihit];
      displayEvent.hits.push_back(hit);
    }
    gCDetDisplayEvents.push_back(displayEvent);
  }

  std::cout << "[CDet event display] Built " << gCDetDisplayEvents.size()
            << " events containing accepted pairs for Layer-1 bar "
            << selectedBar << ".\n";
  if (gCDetDisplayEvents.empty()) return;

  if (std::isfinite(bestHitPeakMean) && std::isfinite(bestHitPeakSigma) &&
      bestHitPeakSigma > 0.0) {
    gCDetDisplayBestHitPeakMean = bestHitPeakMean;
    gCDetDisplayBestHitPeakSigma = bestHitPeakSigma;
    gCDetDisplayBestHitTimingValid = true;
    std::cout << "[CDet event display] Using configured best-hit ECal-CDet peak: "
              << gCDetDisplayBestHitPeakMean << " +/- "
              << gCDetDisplayBestHitNSigma << "*"
              << gCDetDisplayBestHitPeakSigma << " ns.\n";
  } else {
    std::vector<double> timingSample;
    for (const auto& event : gCDetDisplayEvents) {
      timingSample.reserve(timingSample.size() + 2*event.selectedPairs.size());
      for (const auto& pair : event.selectedPairs) {
        timingSample.push_back(event.ecalTime - pair.t1);
        timingSample.push_back(event.ecalTime - pair.t2);
      }
    }
    if (timingSample.size() >= 20) {
      const auto limits = std::minmax_element(timingSample.begin(), timingSample.end());
      const double sampleMin = std::max(tdiffECalCDetMin, *limits.first - 0.5);
      const double sampleMax = std::min(tdiffECalCDetMax, *limits.second + 0.5);
      const int nBins = std::max(40, std::min(1200,
          static_cast<int>(std::ceil((sampleMax - sampleMin)/0.25))));
      TH1D timingHistogram("hCDetDisplayBestHitTimingFit", "", nBins,
                           sampleMin, sampleMax);
      timingHistogram.SetDirectory(nullptr);
      for (double value : timingSample) timingHistogram.Fill(value);
      const double mode = timingHistogram.GetBinCenter(timingHistogram.GetMaximumBin());
      const double fitMin = std::max(sampleMin, mode - 5.0);
      const double fitMax = std::min(sampleMax, mode + 5.0);
      TF1 peakFit("fCDetDisplayBestHitTimingFit", "gaus", fitMin, fitMax);
      peakFit.SetParameters(timingHistogram.GetMaximum(), mode, 2.0);
      peakFit.SetParLimits(2, 0.1, 10.0);
      const TFitResultPtr fitResult = timingHistogram.Fit(&peakFit, "QRSN");
      const double fittedMean = peakFit.GetParameter(1);
      const double fittedSigma = std::fabs(peakFit.GetParameter(2));
      gCDetDisplayBestHitTimingValid = int(fitResult) == 0 &&
          std::isfinite(fittedMean) && std::isfinite(fittedSigma) &&
          fittedSigma > 0.0 && fittedSigma <= 10.0 &&
          std::fabs(fittedMean - mode) <= 5.0;
      if (gCDetDisplayBestHitTimingValid) {
        gCDetDisplayBestHitPeakMean = fittedMean;
        gCDetDisplayBestHitPeakSigma = fittedSigma;
        std::cout << "[CDet event display] Fitted best-hit ECal-CDet peak from "
                  << timingSample.size() << " accepted-pair hits: mean="
                  << fittedMean << " ns, sigma=" << fittedSigma
                  << " ns; display window is +/- " << gCDetDisplayBestHitNSigma
                  << " sigma.\n";
      }
    }
    if (!gCDetDisplayBestHitTimingValid)
      std::cerr << "[CDet event display] WARNING: best-hit peak fit is unavailable; "
                   "the display will remain in all-hits mode.\n";
  }
  gCDetDisplayBestHits = initialBestHits && gCDetDisplayBestHitTimingValid;

  ShowCDetEvent(0);
  if (!gROOT->IsBatch()) {
    // A GUI close can destroy a TControlBar without clearing our pointer.
    // Build a fresh bar for each newly built event list and never dereference
    // a control pointer retained from an earlier invocation.
    gCDetEventControl = new TControlBar("vertical", "CDet event display");
    gCDetEventControl->AddButton("Previous", "PreviousCDetEvent()", "Show the previous selected event");
    gCDetEventControl->AddButton("Next", "NextCDetEvent()", "Show the next selected event");
    gCDetEventControl->AddButton("All hits", "ShowAllCDetHits()", "Display every retained good CDet hit");
    gCDetEventControl->AddButton("Best CDet hits", "ShowBestCDetHits()", "Display hits within the configured number of sigma of the ECal-CDet peak");
    gCDetEventControl->AddButton("Print", "PrintCDetEvent()", "Print the current event and pair values");
    gCDetEventControl->AddButton("Save PNG", "SaveCDetEvent()", "Save the current four-panel event display");
    gCDetEventControl->Show();
  }
}

static bool CDetDisplayTimeIsBest(double ecalTime, double correctedCDetTime)
{
  if (!gCDetDisplayBestHitTimingValid) return false;
  const double dt = ecalTime - correctedCDetTime;
  return std::fabs(dt - gCDetDisplayBestHitPeakMean) <=
         gCDetDisplayBestHitNSigma*gCDetDisplayBestHitPeakSigma;
}

static bool CDetDisplayHitIsVisible(const CDetDisplayEvent& event,
                                    const CDetDisplayHit& hit)
{
  return !gCDetDisplayBestHits || CDetDisplayTimeIsBest(event.ecalTime, hit.le);
}

static bool CDetDisplayPairIsVisible(const CDetDisplayEvent& event,
                                     const PairHit& pair)
{
  return !gCDetDisplayBestHits ||
         (CDetDisplayTimeIsBest(event.ecalTime, pair.t1) &&
          CDetDisplayTimeIsBest(event.ecalTime, pair.t2));
}

void ShowCDetEvent(Long64_t displayIndex)
{
  if (gCDetDisplayEvents.empty()) {
    std::cerr << "[CDet event display] No display events are available. Run plotCDetLayersTimeComp first.\n";
    return;
  }
  const Long64_t nEvents = static_cast<Long64_t>(gCDetDisplayEvents.size());
  displayIndex %= nEvents;
  if (displayIndex < 0) displayIndex += nEvents;
  gCDetDisplayIndex = displayIndex;
  const auto& event = gCDetDisplayEvents[gCDetDisplayIndex];

  // Closing a ROOT canvas deletes it, but does not update ordinary external
  // pointers. Resolve the canvas through ROOT's live-canvas registry before
  // every redraw so a closed window cannot leave us with a dangling pointer.
  gCDetEventCanvas = dynamic_cast<TCanvas *>(
      gROOT->GetListOfCanvases()->FindObject("cCDetEventDisplay"));
  if (!gCDetEventCanvas) {
    gCDetEventCanvas = new TCanvas("cCDetEventDisplay", "CDet event display", 1500, 950);
    gCDetEventCanvas->Divide(2, 2, 0.005, 0.005);
  }

  std::vector<double> l1z, l1x, l1y, l2z, l2x, l2y;
  for (const auto& hit : event.hits) {
    if (!CDetDisplayHitIsVisible(event, hit)) continue;
    if (hit.layer == 0) {
      l1z.push_back(hit.z); l1x.push_back(hit.x); l1y.push_back(hit.y);
    } else {
      l2z.push_back(hit.z); l2x.push_back(hit.x); l2y.push_back(hit.y);
    }
  }
  double displayYMin = -0.5;
  double displayYMax = 0.5;
  for (const auto& hit : event.hits) {
    if (!CDetDisplayHitIsVisible(event, hit)) continue;
    displayYMin = std::min(displayYMin, hit.y - CDet_y_half_length - 0.05);
    displayYMax = std::max(displayYMax, hit.y + CDet_y_half_length + 0.05);
  }

  const TString title = TString::Format(
      "Run %d, tree entry %lld, selected event %lld/%lld, Layer-1 bar %d [%s]",
      event.runNumber, event.treeEntry, gCDetDisplayIndex + 1, nEvents,
      event.selectedBar, gCDetDisplayBestHits ? "best CDet hits" : "all hits");

  gCDetEventCanvas->cd(1);
  gPad->Clear();
  gPad->DrawFrame(gCDetDisplayZMin, gCDetDisplayXMin,
                  gCDetDisplayZMax, gCDetDisplayXMax,
                  title + ";z position (m);x position (m)");
  TLine *ecalTrajectoryX = new TLine(0.0, 0.0, event.ecalZ, event.ecalX);
  ecalTrajectoryX->SetLineColor(kGray + 2);
  ecalTrajectoryX->SetLineStyle(2);
  ecalTrajectoryX->SetLineWidth(2);
  ecalTrajectoryX->Draw();
  if (!l1z.empty()) {
    TGraph *gL1 = new TGraph(l1z.size(), l1z.data(), l1x.data());
    gL1->SetMarkerStyle(24); gL1->SetMarkerColor(kBlue + 1); gL1->SetMarkerSize(1.2); gL1->Draw("P SAME");
  }
  if (!l2z.empty()) {
    TGraph *gL2 = new TGraph(l2z.size(), l2z.data(), l2x.data());
    gL2->SetMarkerStyle(25); gL2->SetMarkerColor(kGreen + 2); gL2->SetMarkerSize(1.2); gL2->Draw("P SAME");
  }
  for (const auto& pair : event.selectedPairs) {
    if (!CDetDisplayPairIsVisible(event, pair)) continue;
    const double correctedX1 = pair.x1;
    const double correctedX2 = pair.x2;
    TLine *pairLine = new TLine(pair.z1, correctedX1, pair.z2, correctedX2);
    pairLine->SetLineColor(kRed + 1); pairLine->SetLineWidth(3); pairLine->Draw();
    const double projectedX1 = event.ecalX * pair.z1 / event.ecalZ;
    const double projectedX2 = event.ecalX * pair.z2 / event.ecalZ;
    TLine *residual1 = new TLine(pair.z1, projectedX1, pair.z1, correctedX1);
    TLine *residual2 = new TLine(pair.z2, projectedX2, pair.z2, correctedX2);
    residual1->SetLineColor(kMagenta + 1); residual1->SetLineWidth(2); residual1->Draw();
    residual2->SetLineColor(kMagenta + 1); residual2->SetLineWidth(2); residual2->Draw();
  }
  {
    double z = event.ecalZ, x = event.ecalX;
    TGraph *gECal = new TGraph(1, &z, &x);
    gECal->SetMarkerStyle(29); gECal->SetMarkerColor(kBlack); gECal->SetMarkerSize(2.0); gECal->Draw("P SAME");
  }

  gCDetEventCanvas->cd(2);
  gPad->Clear();
  gPad->DrawFrame(gCDetDisplayZMin, displayYMin, gCDetDisplayZMax, displayYMax,
                  "y-z trajectory;z position (m);y position (m)");
  TLine *ecalTrajectoryY = new TLine(0.0, 0.0, event.ecalZ, event.ecalY);
  ecalTrajectoryY->SetLineColor(kGray + 2); ecalTrajectoryY->SetLineStyle(2); ecalTrajectoryY->SetLineWidth(2); ecalTrajectoryY->Draw();
  for (const auto& hit : event.hits) {
    if (!CDetDisplayHitIsVisible(event, hit)) continue;
    TLine *paddleExtent = new TLine(hit.z, hit.y - CDet_y_half_length,
                                    hit.z, hit.y + CDet_y_half_length);
    paddleExtent->SetLineColor(hit.layer == 0 ? kBlue + 1 : kGreen + 2);
    paddleExtent->SetLineWidth(3);
    paddleExtent->Draw();
  }
  for (const auto& pair : event.selectedPairs) {
    if (!CDetDisplayPairIsVisible(event, pair)) continue;
    TLine *selectedL1 = new TLine(pair.z1, pair.y1 - CDet_y_half_length,
                                  pair.z1, pair.y1 + CDet_y_half_length);
    TLine *selectedL2 = new TLine(pair.z2, pair.y2 - CDet_y_half_length,
                                  pair.z2, pair.y2 + CDet_y_half_length);
    selectedL1->SetLineColor(kRed + 1); selectedL1->SetLineWidth(5); selectedL1->Draw();
    selectedL2->SetLineColor(kRed + 1); selectedL2->SetLineWidth(5); selectedL2->Draw();
  }
  {
    double z = event.ecalZ, y = event.ecalY;
    TGraph *gECalY = new TGraph(1, &z, &y);
    gECalY->SetMarkerStyle(29); gECalY->SetMarkerColor(kBlack); gECalY->SetMarkerSize(2.0); gECalY->Draw("P SAME");
  }

  gCDetEventCanvas->cd(3);
  gPad->Clear();
  gPad->DrawFrame(gCDetDisplayXMin, displayYMin, gCDetDisplayXMax, displayYMax,
                  "CDet face view;x position (m);y position (m)");
  for (const auto& hit : event.hits) {
    if (!CDetDisplayHitIsVisible(event, hit)) continue;
    TLine *paddleExtent = new TLine(hit.x, hit.y - CDet_y_half_length,
                                    hit.x, hit.y + CDet_y_half_length);
    paddleExtent->SetLineColor(hit.layer == 0 ? kBlue + 1 : kGreen + 2);
    paddleExtent->SetLineWidth(3);
    paddleExtent->Draw();
  }
  for (const auto& pair : event.selectedPairs) {
    if (!CDetDisplayPairIsVisible(event, pair)) continue;
    double px[2] = {pair.x1, pair.x2};
    double py[2] = {pair.y1, pair.y2};
    TGraph *gPair = new TGraph(2, px, py);
    gPair->SetLineColor(kRed + 1); gPair->SetLineWidth(3); gPair->Draw("L SAME");
    TLine *selectedL1 = new TLine(px[0], pair.y1 - CDet_y_half_length,
                                  px[0], pair.y1 + CDet_y_half_length);
    TLine *selectedL2 = new TLine(px[1], pair.y2 - CDet_y_half_length,
                                  px[1], pair.y2 + CDet_y_half_length);
    selectedL1->SetLineColor(kRed + 1); selectedL1->SetLineWidth(5); selectedL1->Draw();
    selectedL2->SetLineColor(kRed + 1); selectedL2->SetLineWidth(5); selectedL2->Draw();
  }
  const PairHit *firstVisiblePair = nullptr;
  for (const auto& pair : event.selectedPairs) {
    if (CDetDisplayPairIsVisible(event, pair)) {
      firstVisiblePair = &pair;
      break;
    }
  }
  if (firstVisiblePair) {
    // CDet layers have fixed z positions, so the first accepted pair supplies
    // the two layer planes for the ECal-to-target projection.
    const auto& pair = *firstVisiblePair;
    const double projectedX1 = event.ecalX * pair.z1 / event.ecalZ;
    const double projectedY1 = event.ecalY * pair.z1 / event.ecalZ;
    const double projectedX2 = event.ecalX * pair.z2 / event.ecalZ;
    const double projectedY2 = event.ecalY * pair.z2 / event.ecalZ;
    TMarker *projectedECalL1 = new TMarker(projectedX1, projectedY1, 29);
    TMarker *projectedECalL2 = new TMarker(projectedX2, projectedY2, 29);
    projectedECalL1->SetMarkerColor(kBlue + 1); projectedECalL1->SetMarkerSize(2.0); projectedECalL1->Draw();
    projectedECalL2->SetMarkerColor(kGreen + 2); projectedECalL2->SetMarkerSize(2.0); projectedECalL2->Draw();
  }

  gCDetEventCanvas->cd(4);
  gPad->Clear();
  TPaveText *info = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
  info->SetTextAlign(12);
  info->SetTextFont(42);
  info->SetTextSize(0.032);
  info->AddText(title);
  info->AddText(TString::Format("ECal: x=%+.4f m, y=%+.4f m, z=%.4f m", event.ecalX, event.ecalY, event.ecalZ));
  info->AddText(TString::Format("ECal: E=%.4f GeV, ADC time=%.3f ns", event.ecalEnergy, event.ecalTime));
  const size_t visibleHits = std::count_if(
      event.hits.begin(), event.hits.end(),
      [&](const CDetDisplayHit& hit) { return CDetDisplayHitIsVisible(event, hit); });
  const size_t visiblePairs = std::count_if(
      event.selectedPairs.begin(), event.selectedPairs.end(),
      [&](const PairHit& pair) { return CDetDisplayPairIsVisible(event, pair); });
  info->AddText(TString::Format("Displayed CDet hits: %zu/%zu; selected pairs: %zu/%zu",
                                visibleHits, event.hits.size(),
                                visiblePairs, event.selectedPairs.size()));
  if (gCDetDisplayBestHitTimingValid) {
    info->AddText(TString::Format("Best-hit timing: t_{ECal}-t_{CDet} = %.3f +/- %.1f#sigma (sigma=%.3f ns)",
                                  gCDetDisplayBestHitPeakMean,
                                  gCDetDisplayBestHitNSigma,
                                  gCDetDisplayBestHitPeakSigma));
  }
  for (size_t ipair = 0; ipair < event.selectedPairs.size(); ++ipair) {
    const auto& pair = event.selectedPairs[ipair];
    if (!CDetDisplayPairIsVisible(event, pair)) continue;
    const double correctedX1 = pair.x1;
    const double correctedX2 = pair.x2;
    const double pairTime = 0.5 * (pair.t1 + pair.t2);
    const double projectedX1 = event.ecalX * pair.z1 / event.ecalZ;
    const double projectedX2 = event.ecalX * pair.z2 / event.ecalZ;
    info->AddText(TString::Format("Pair %zu: L1 ID %d  L2 ID %d  dt=%+.3f ns  corrected dx=%+.4f m  score=%.3f",
                                  ipair + 1, pair.id1, pair.id2, pair.dt, pair.dx, pair.score));
    info->AddText(TString::Format("  L1: corrected x=%+.4f y=%+.4f z=%.4f LE=%.3f TOT=%.3f residual=%+.4f m",
                                  correctedX1, pair.y1, pair.z1, pair.t1, pair.tot1, correctedX1 - projectedX1));
    info->AddText(TString::Format("  L2: corrected x=%+.4f y=%+.4f z=%.4f LE=%.3f TOT=%.3f residual=%+.4f m",
                                  correctedX2, pair.y2, pair.z2, pair.t2, pair.tot2, correctedX2 - projectedX2));
    info->AddText(TString::Format("  ECal-CDet time=%+.3f ns", event.ecalTime - pairTime));
  }
  info->Draw();

  gCDetEventCanvas->Modified();
  gCDetEventCanvas->Update();
}

void NextCDetEvent()
{
  ShowCDetEvent(gCDetDisplayIndex + 1);
}

void PreviousCDetEvent()
{
  ShowCDetEvent(gCDetDisplayIndex - 1);
}

void ShowAllCDetHits()
{
  gCDetDisplayBestHits = false;
  if (gCDetDisplayIndex >= 0) ShowCDetEvent(gCDetDisplayIndex);
}

void ShowBestCDetHits()
{
  if (!gCDetDisplayBestHitTimingValid) {
    std::cerr << "[CDet event display] Best-hit timing peak is unavailable.\n";
    return;
  }
  gCDetDisplayBestHits = true;
  if (gCDetDisplayIndex >= 0) ShowCDetEvent(gCDetDisplayIndex);
}

void PrintCDetEvent()
{
  if (gCDetDisplayIndex < 0 || gCDetDisplayIndex >= static_cast<Long64_t>(gCDetDisplayEvents.size())) {
    std::cerr << "[CDet event display] No current event.\n";
    return;
  }
  const auto& event = gCDetDisplayEvents[gCDetDisplayIndex];
  const size_t visibleHits = std::count_if(
      event.hits.begin(), event.hits.end(),
      [&](const CDetDisplayHit& hit) { return CDetDisplayHitIsVisible(event, hit); });
  const size_t visiblePairs = std::count_if(
      event.selectedPairs.begin(), event.selectedPairs.end(),
      [&](const PairHit& pair) { return CDetDisplayPairIsVisible(event, pair); });
  std::cout << "[CDet event display] run=" << event.runNumber
            << " tree_entry=" << event.treeEntry
            << " display_index=" << gCDetDisplayIndex
            << " selected_bar=" << event.selectedBar
            << " ecal_x=" << event.ecalX
            << " ecal_y=" << event.ecalY
            << " ecal_E=" << event.ecalEnergy
            << " ecal_t=" << event.ecalTime
            << " hit_mode=" << (gCDetDisplayBestHits ? "best" : "all")
            << " nhits=" << visibleHits << "/" << event.hits.size()
            << " npairs=" << visiblePairs << "/" << event.selectedPairs.size() << "\n";
  for (size_t ipair = 0; ipair < event.selectedPairs.size(); ++ipair) {
    const auto& pair = event.selectedPairs[ipair];
    if (!CDetDisplayPairIsVisible(event, pair)) continue;
    std::cout << "  pair=" << ipair
              << " id1=" << pair.id1 << " id2=" << pair.id2
              << " corrected_x1=" << pair.x1
              << " corrected_x2=" << pair.x2
              << " y1=" << pair.y1 << " y2=" << pair.y2
              << " z1=" << pair.z1 << " z2=" << pair.z2
              << " t1=" << pair.t1 << " t2=" << pair.t2
              << " tot1=" << pair.tot1 << " tot2=" << pair.tot2
              << " dt=" << pair.dt << " corrected_dx=" << pair.dx
              << " score=" << pair.score << "\n";
  }
}

void SaveCDetEvent()
{
  gCDetEventCanvas = dynamic_cast<TCanvas *>(
      gROOT->GetListOfCanvases()->FindObject("cCDetEventDisplay"));
  if (!gCDetEventCanvas || gCDetDisplayIndex < 0 ||
      gCDetDisplayIndex >= static_cast<Long64_t>(gCDetDisplayEvents.size())) {
    std::cerr << "[CDet event display] No current event to save.\n";
    return;
  }
  const auto& event = gCDetDisplayEvents[gCDetDisplayIndex];
  const TString fileName = TString::Format(
      "CDetEventDisplay_run%d_entry%lld_bar%d_%s.png",
      event.runNumber, event.treeEntry, event.selectedBar,
      gCDetDisplayBestHits ? "best_hits" : "all_hits");
  gCDetEventCanvas->SaveAs(fileName);
  std::cout << "[CDet event display] Saved " << fileName << "\n";
}

void plotECalCDetTimeCutStudy(double Width = 1.0, int logicalPixelID = 485,
                              double dtMinCut = -30.0, double dtMaxCut = 0.0,
                              double DiffMin = -60.0, double DiffMax = 30.0,
                              double LeMin = 0.0, double LeMax = 60.0,
                              double TeMin = 0.0, double TeMax = 80.0,
                              double TotMin = 0.0, double TotMax = 80.0,
                              double ECalMin = -40.0, double ECalMax = 40.0,
                              bool drawDtVsTot = false,
                              double ECalEnergyMin = 1.0, double ECalEnergyMax = 12.0,
                              double localFitHalfWidth = 8.0, double NReject = 2.5,
                              double minSigma = 0.5, double maxSigma = 20.0,
                              double maxChi2Ndf = 10.0, double centroidEdgeMargin = 1.0,
                              double PeakSeedMin = -25.0, double PeakSeedMax = -5.0){
  TH1::AddDirectory(kFALSE);
  (void)localFitHalfWidth; // Retained for positional compatibility; timing-study fits now use dtMinCut-dtMaxCut exactly.

  if (logicalPixelID < 0 || logicalPixelID >= NumCDetPaddles) {
    std::cerr << "[CDet timing-cut study] ERROR: logical pixel ID " << logicalPixelID
              << " is outside the physical CDet range [0, " << NumCDetPaddles - 1
              << "].\n";
    return;
  }
  if (dtMinCut >= dtMaxCut) {
    std::cerr << "[CDet timing-cut study] ERROR: invalid delta-t window ["
              << dtMinCut << ", " << dtMaxCut << "] ns.\n";
    return;
  }
  if (Width <= 0.0 || DiffMin >= DiffMax || dtMinCut < DiffMin || dtMaxCut > DiffMax ||
      PeakSeedMin >= PeakSeedMax || PeakSeedMin < dtMinCut || PeakSeedMax > dtMaxCut ||
      minSigma <= 0.0 || minSigma >= maxSigma || maxChi2Ndf <= 0.0 || centroidEdgeMargin < 0.0 ||
      LeMin >= LeMax || TeMin >= TeMax || TotMin >= TotMax || ECalMin >= ECalMax ||
      ECalEnergyMin >= ECalEnergyMax || NReject <= 0.0) {
    std::cerr << "[CDet timing-cut study] ERROR: invalid histogram binning or range.\n";
    return;
  }

  const size_t nEvents = vGoodLe.size();
  if (nEvents == 0) {
    std::cerr << "[CDet timing-cut study] ERROR: no fourth-pass baseline events are available.\n";
    return;
  }
  if (v_GoodECalE.size() != nEvents) {
    std::cerr << "[CDet timing-cut study] ERROR: ECal energy vector is not aligned with the fourth-pass baseline events: events=" << nEvents << ", ECal energy=" << v_GoodECalE.size() << ".\n";
    return;
  }

  int NDiffBins = (int)((DiffMax - DiffMin)/Width);
  int NLeBins = (int)((LeMax - LeMin)/Width);
  int NTeBins = (int)((TeMax - TeMin)/Width);
  int NTotBinsStudy = (int)((TotMax - TotMin)/Width);
  int NADCBins = (int)((ECalMax - ECalMin)/4.0);

  static unsigned long invocation = 0;
  const unsigned long tag = ++invocation;
  auto uniqueName = [tag](const char *base) {
    return TString::Format("%s_%lu", base, tag);
  };

  TH2D *hDtVsPixel = new TH2D(uniqueName("hCDetECalCutStudyDtVsPixel"), "Before #Deltat cut, after ECal energy cut;CDet logical pixel ID;t_{ECal}-t_{CDet,LE} (ns)", NumCDetPaddles, -0.5, NumCDetPaddles - 0.5, NDiffBins, DiffMin, DiffMax);
  TH1D *hDt = new TH1D(uniqueName("hCDetECalCutStudyDt"), "Before #Deltat cut, after ECal energy cut;t_{ECal}-t_{CDet,LE} (ns);Baseline hits", NDiffBins, DiffMin, DiffMax);
  TH1D *hSelectedPixelDt = new TH1D(uniqueName("hCDetECalCutStudySelectedPixelDt"), TString::Format("Before #Deltat cut, after ECal energy cut, logical pixel ID %d;t_{ECal}-t_{CDet,LE} (ns);Baseline hits", logicalPixelID), NDiffBins, DiffMin, DiffMax);
  TH1D *hLePass = new TH1D(uniqueName("hCDetECalCutStudyLePass"), "CDet LE after #Deltat cut;CDet LE time (ns);Passing hits", NLeBins, LeMin, LeMax);
  TH1D *hTePass = new TH1D(uniqueName("hCDetECalCutStudyTePass"), "CDet TE after #Deltat cut;CDet TE time (ns);Passing hits", NTeBins, TeMin, TeMax);
  TH1D *hTotPass = new TH1D(uniqueName("hCDetECalCutStudyTotPass"), "CDet TOT after #Deltat cut;CDet TOT (ns);Passing hits", NTotBinsStudy, TotMin, TotMax);
  TH1D *hSelectedPixelLe = new TH1D(uniqueName("hCDetECalCutStudySelectedPixelLe"), TString::Format("CDet LE after #Deltat cut, logical pixel ID %d;CDet LE time (ns);Passing hits", logicalPixelID), NLeBins, LeMin, LeMax);
  TH2D *hECalVsCDetPass = new TH2D(uniqueName("hCDetECalCutStudyECalVsCDetPass"), "ECal time versus CDet LE after #Deltat cut;CDet LE time (ns);ECal ADC time (ns)", NLeBins, LeMin, LeMax, NADCBins, ECalMin, ECalMax);
  TH2D *hECalVsSelectedPixelPass = new TH2D(uniqueName("hCDetECalCutStudyECalVsSelectedPixelPass"), TString::Format("ECal time versus CDet LE after #Deltat cut, logical pixel ID %d;CDet LE time (ns);ECal ADC time (ns)", logicalPixelID), NLeBins, LeMin, LeMax, NADCBins, ECalMin, ECalMax);
  TH2D *hDtVsTot = nullptr;
  if (drawDtVsTot) {
    hDtVsTot = new TH2D(uniqueName("hCDetECalCutStudyDtVsTot"), TString::Format("Before #Deltat cut, after ECal energy cut, logical pixel ID %d;CDet TOT (ns);t_{ECal}-t_{CDet,LE} (ns)", logicalPixelID), NTotBinsStudy, TotMin, TotMax, NDiffBins, DiffMin, DiffMax);
  }

  // std::vector<unsigned long> baselinePerPixel(NumCDetPaddles, 0);
  // std::vector<unsigned long> passingPerPixel(NumCDetPaddles, 0);
  size_t baselineHitCount = 0;
  size_t passingHitCount = 0;
  size_t selectedPixelPassingCount = 0;
  size_t energySelectedEventCount = 0;

  for (size_t ev = 0; ev < nEvents; ++ev) {
    const double ecalEnergy = v_GoodECalE[ev];
    if (ecalEnergy < ECalEnergyMin || ecalEnergy > ECalEnergyMax) continue;
    ++energySelectedEventCount;
    const double tECal = v_GoodECalAdcTime[ev];
    for (size_t ihit = 0; ihit < vGoodLe[ev].size(); ++ihit) {
      const double le = vGoodLe[ev][ihit];
      const double te = vGoodTe[ev][ihit];
      const double tot = vGoodTot[ev][ihit];
      const int pixelID = vGoodID[ev][ihit];
      const double dt = tECal - le;

      hDtVsPixel->Fill(pixelID, dt);
      hDt->Fill(dt);
      if (pixelID == logicalPixelID) {
        hSelectedPixelDt->Fill(dt);
        if (hDtVsTot) hDtVsTot->Fill(tot, dt);
      }
      ++baselineHitCount;
      // ++baselinePerPixel[pixelID];

      // This is the only additional hit-quality requirement in this routine.
      if (dtMinCut <= dt && dt <= dtMaxCut) {
        ++passingHitCount;
        // ++passingPerPixel[pixelID];
        hLePass->Fill(le);
        hTePass->Fill(te);
        hTotPass->Fill(tot);
        hECalVsCDetPass->Fill(le, tECal);
        if (pixelID == logicalPixelID) {
          ++selectedPixelPassingCount;
          hSelectedPixelLe->Fill(le);
          hECalVsSelectedPixelPass->Fill(le, tECal);
        }
      }
    }
  }

  // TGraphErrors *gEfficiency = new TGraphErrors();
  // gEfficiency->SetName(uniqueName("gCDetECalCutStudyEfficiency"));
  // gEfficiency->SetTitle("ECal-CDet #Deltat-cut efficiency versus CDet logical pixel ID;CDet logical pixel ID;Hit efficiency");
  // for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
  //   const double denominator = baselinePerPixel[pixelID];
  //   if (denominator == 0.0) continue;
  //   const double efficiency = passingPerPixel[pixelID] / denominator;
  //   const double uncertainty = std::sqrt(efficiency * (1.0 - efficiency) / denominator);
  //   const int point = gEfficiency->GetN();
  //   gEfficiency->SetPoint(point, pixelID, efficiency);
  //   gEfficiency->SetPointError(point, 0.0, uncertainty);
  // }
  // gEfficiency->SetMarkerStyle(20);
  // gEfficiency->SetMarkerSize(0.45);
  // gEfficiency->SetMinimum(0.0);
  // gEfficiency->SetMaximum(1.05);

  TCanvas *cDtVsPixel = new TCanvas(uniqueName("cCDetECalCutStudyDtVsPixel"), "Uncut ECal-CDet delta-t versus logical pixel ID", 1100, 700);
  hDtVsPixel->Draw("COLZ");

  TCanvas *cDt = new TCanvas(uniqueName("cCDetECalCutStudyDt"), "Uncut ECal-CDet delta-t", 900, 700);
  hDt->Draw("HIST");
  const double lineTop = std::max(1.0, 1.05 * hDt->GetMaximum());
  TLine *dtMinLine = new TLine(dtMinCut, 0.0, dtMinCut, lineTop);
  TLine *dtMaxLine = new TLine(dtMaxCut, 0.0, dtMaxCut, lineTop);
  dtMinLine->SetLineColor(kRed + 1);
  dtMaxLine->SetLineColor(kRed + 1);
  dtMinLine->SetLineWidth(2);
  dtMaxLine->SetLineWidth(2);
  dtMinLine->SetLineStyle(2);
  dtMaxLine->SetLineStyle(2);
  dtMinLine->Draw("SAME");
  dtMaxLine->Draw("SAME");

  TF1 *fSelectedPixelLocal = nullptr;
  TF1 *fSelectedPixelBackgroundReject = nullptr;
  TF1 *fSelectedPixelBackground = nullptr;
  TF1 *fSelectedPixelClean = nullptr;
  TH1D *hSelectedPixelDtClean = nullptr;
  bool selectedPixelLocalValid = false;
  bool selectedPixelBackgroundValid = false;
  bool selectedPixelCleanValid = false;
  int selectedPixelLocalStatus = -1;
  int selectedPixelBackgroundStatus = -1;
  int selectedPixelCleanStatus = -1;
  double selectedPixelLocalMean = NAN;
  double selectedPixelLocalSigma = NAN;
  double selectedPixelRejectLow = NAN;
  double selectedPixelRejectHigh = NAN;

  const int selectedPixelFitBinMin = hSelectedPixelDt->FindBin(dtMinCut);
  const int selectedPixelFitBinMax = hSelectedPixelDt->FindBin(dtMaxCut);
  if (hSelectedPixelDt->Integral(selectedPixelFitBinMin, selectedPixelFitBinMax) > 0.0) {
    int selectedPixelSearchBinMin = hSelectedPixelDt->FindBin(PeakSeedMin);
    int selectedPixelSearchBinMax = hSelectedPixelDt->FindBin(PeakSeedMax);
    int selectedPixelPeakBin = selectedPixelSearchBinMin;
    for (int bin = selectedPixelSearchBinMin + 1; bin <= selectedPixelSearchBinMax; ++bin) {
      if (hSelectedPixelDt->GetBinContent(bin) > hSelectedPixelDt->GetBinContent(selectedPixelPeakBin)) selectedPixelPeakBin = bin;
    }
    if (hSelectedPixelDt->GetBinContent(selectedPixelPeakBin) <= 0.0) {
      selectedPixelSearchBinMin = selectedPixelFitBinMin;
      selectedPixelSearchBinMax = selectedPixelFitBinMax;
      selectedPixelPeakBin = selectedPixelSearchBinMin;
      for (int bin = selectedPixelSearchBinMin + 1; bin <= selectedPixelSearchBinMax; ++bin) {
        if (hSelectedPixelDt->GetBinContent(bin) > hSelectedPixelDt->GetBinContent(selectedPixelPeakBin)) selectedPixelPeakBin = bin;
      }
    }
    const double selectedPixelPeak = hSelectedPixelDt->GetBinCenter(selectedPixelPeakBin);
    const double localFitMin = dtMinCut;
    const double localFitMax = dtMaxCut;
    const double selectedPixelLocalBackground = 0.5*(hSelectedPixelDt->GetBinContent(hSelectedPixelDt->FindBin(localFitMin)) + hSelectedPixelDt->GetBinContent(hSelectedPixelDt->FindBin(localFitMax)));
    fSelectedPixelLocal = new TF1(uniqueName("fCDetECalCutStudySelectedPixelLocal"), "gaus(0)+pol1(3)", localFitMin, localFitMax);
    fSelectedPixelLocal->SetParameters(std::max(1.0, hSelectedPixelDt->GetBinContent(selectedPixelPeakBin) - selectedPixelLocalBackground), selectedPixelPeak, std::max(Width, (dtMaxCut - dtMinCut)/10.0), selectedPixelLocalBackground, 0.0);
    fSelectedPixelLocal->SetParNames("Local amplitude", "Local mean", "Local sigma", "Local background intercept", "Local background slope");
    selectedPixelLocalStatus = hSelectedPixelDt->Fit(fSelectedPixelLocal, "RQN0");
    selectedPixelLocalMean = fSelectedPixelLocal->GetParameter(1);
    selectedPixelLocalSigma = std::fabs(fSelectedPixelLocal->GetParameter(2));
    selectedPixelLocalValid = selectedPixelLocalStatus == 0 && std::isfinite(fSelectedPixelLocal->GetParameter(0)) && fSelectedPixelLocal->GetParameter(0) > 0.0 && std::isfinite(selectedPixelLocalMean) && std::isfinite(selectedPixelLocalSigma) && std::isfinite(fSelectedPixelLocal->GetParameter(3)) && std::isfinite(fSelectedPixelLocal->GetParameter(4)) && selectedPixelLocalSigma >= minSigma && selectedPixelLocalSigma <= maxSigma && selectedPixelLocalMean > dtMinCut && selectedPixelLocalMean < dtMaxCut;

    /* Background subtraction was removed because the production timing mean
       comes from the local gaus+pol1 fit and the cleaned plots were unused. */
#if 0
    if (selectedPixelLocalValid) {
      selectedPixelRejectLow = std::max(DiffMin, selectedPixelLocalMean - NReject*selectedPixelLocalSigma);
      selectedPixelRejectHigh = std::min(DiffMax, selectedPixelLocalMean + NReject*selectedPixelLocalSigma);
      CDetTimingBackgroundRejectLow = selectedPixelRejectLow;
      CDetTimingBackgroundRejectHigh = selectedPixelRejectHigh;
      double maxSidebandBinContent = 0.0;
      int leftSidebandBins = 0;
      int rightSidebandBins = 0;
      for (int bin = 1; bin <= hSelectedPixelDt->GetNbinsX(); ++bin) {
        const double binCenter = hSelectedPixelDt->GetBinCenter(bin);
        if (binCenter < DiffMin || binCenter > DiffMax || (selectedPixelRejectLow <= binCenter && binCenter <= selectedPixelRejectHigh)) continue;
        maxSidebandBinContent = std::max(maxSidebandBinContent, hSelectedPixelDt->GetBinContent(bin));
        if (binCenter < selectedPixelRejectLow) ++leftSidebandBins;
        if (binCenter > selectedPixelRejectHigh) ++rightSidebandBins;
      }
      if (maxSidebandBinContent > 0.0 && leftSidebandBins >= 3 && rightSidebandBins >= 3) {
        fSelectedPixelBackgroundReject = new TF1(uniqueName("fCDetECalCutStudySelectedPixelBackgroundReject"), CDetTimingBackgroundGaussianReject, DiffMin, DiffMax, 3);
        fSelectedPixelBackgroundReject->SetParameters(0.8*maxSidebandBinContent, 0.5*(DiffMin + DiffMax), std::max(Width, (DiffMax - DiffMin)/3.0));
        fSelectedPixelBackgroundReject->SetParLimits(0, 0.0, maxSidebandBinContent);
        fSelectedPixelBackgroundReject->SetParLimits(1, DiffMin, DiffMax);
        fSelectedPixelBackgroundReject->SetParLimits(2, Width, DiffMax - DiffMin);
        fSelectedPixelBackgroundReject->SetParNames("Background amplitude", "Background mean", "Background sigma");
        selectedPixelBackgroundStatus = hSelectedPixelDt->Fit(fSelectedPixelBackgroundReject, "RQN0");
        const double backgroundSigma = std::fabs(fSelectedPixelBackgroundReject->GetParameter(2));
        selectedPixelBackgroundValid = selectedPixelBackgroundStatus == 0 && std::isfinite(fSelectedPixelBackgroundReject->GetParameter(0)) && std::isfinite(fSelectedPixelBackgroundReject->GetParameter(1)) && std::isfinite(backgroundSigma) && fSelectedPixelBackgroundReject->GetParameter(0) >= 0.0 && backgroundSigma > 0.0;

        if (selectedPixelBackgroundValid) {
          fSelectedPixelBackground = new TF1(uniqueName("fCDetECalCutStudySelectedPixelBackground"), "gaus", DiffMin, DiffMax);
          fSelectedPixelBackground->SetParameters(fSelectedPixelBackgroundReject->GetParameter(0), fSelectedPixelBackgroundReject->GetParameter(1), backgroundSigma);
          hSelectedPixelDtClean = (TH1D*)hSelectedPixelDt->Clone(uniqueName("hCDetECalCutStudySelectedPixelDtClean"));
          hSelectedPixelDtClean->SetTitle(TString::Format("Background-subtracted ECal-CDet #Deltat, logical pixel ID %d;t_{ECal}-t_{CDet,LE} (ns);Signal estimate", logicalPixelID));
          for (int bin = 1; bin <= hSelectedPixelDtClean->GetNbinsX(); ++bin) {
            const double originalError = hSelectedPixelDt->GetBinError(bin);
            hSelectedPixelDtClean->SetBinContent(bin, hSelectedPixelDt->GetBinContent(bin) - fSelectedPixelBackground->Eval(hSelectedPixelDt->GetBinCenter(bin)));
            hSelectedPixelDtClean->SetBinError(bin, originalError);
          }
          const double cleanFitMin = dtMinCut;
          const double cleanFitMax = dtMaxCut;
          fSelectedPixelClean = new TF1(uniqueName("fCDetECalCutStudySelectedPixelClean"), "gaus", cleanFitMin, cleanFitMax);
          fSelectedPixelClean->SetParameters(std::max(1.0, fSelectedPixelLocal->GetParameter(0)), selectedPixelLocalMean, selectedPixelLocalSigma);
          fSelectedPixelClean->SetParNames("Clean amplitude", "Clean mean", "Clean sigma");
          selectedPixelCleanStatus = hSelectedPixelDtClean->Fit(fSelectedPixelClean, "RQN0");
          const double cleanAmplitude = fSelectedPixelClean->GetParameter(0);
          const double cleanAmplitudeError = fSelectedPixelClean->GetParError(0);
          const double cleanMean = fSelectedPixelClean->GetParameter(1);
          const double cleanSigma = std::fabs(fSelectedPixelClean->GetParameter(2));
          const double cleanMeanError = fSelectedPixelClean->GetParError(1);
          const double cleanSigmaError = fSelectedPixelClean->GetParError(2);
          const int cleanNdf = fSelectedPixelClean->GetNDF();
          const double cleanChi2Ndf = cleanNdf > 0 ? fSelectedPixelClean->GetChisquare()/cleanNdf : NAN;
          selectedPixelCleanValid = selectedPixelCleanStatus == 0 && std::isfinite(cleanAmplitude) && std::isfinite(cleanAmplitudeError) && std::isfinite(cleanMean) && std::isfinite(cleanSigma) && std::isfinite(cleanMeanError) && std::isfinite(cleanSigmaError) && cleanMeanError > 0.0 && cleanSigmaError > 0.0 && cleanMean - cleanFitMin > centroidEdgeMargin && cleanFitMax - cleanMean > centroidEdgeMargin && cleanSigma >= minSigma && cleanSigma <= maxSigma && cleanAmplitude > 0.0 && cleanNdf > 0 && std::isfinite(cleanChi2Ndf) && cleanChi2Ndf <= maxChi2Ndf;
        }
      }
    }
#endif
  }

  TCanvas *cSelectedPixelDt = new TCanvas(uniqueName("cCDetECalCutStudySelectedPixelDt"), TString::Format("ECal-CDet delta-t signal fit, logical pixel ID %d", logicalPixelID), 1000, 750);
  hSelectedPixelDt->SetStats(kTRUE);
  hSelectedPixelDt->SetLineColor(kBlack);
  hSelectedPixelDt->SetLineWidth(2);
  double selectedPixelPlotMin = 0.0;
  double selectedPixelPlotMax = std::max(1.0, 1.15*hSelectedPixelDt->GetMaximum());
  hSelectedPixelDt->SetMinimum(selectedPixelPlotMin);
  hSelectedPixelDt->SetMaximum(selectedPixelPlotMax);
  hSelectedPixelDt->Draw("HIST");
  TLegend *selectedPixelLegend = new TLegend(0.53, 0.62, 0.88, 0.88);
  selectedPixelLegend->AddEntry(hSelectedPixelDt, "Original #Deltat", "l");
  if (selectedPixelLocalValid) {
    fSelectedPixelLocal->SetLineColor(kBlue + 1);
    fSelectedPixelLocal->SetLineStyle(2);
    fSelectedPixelLocal->SetLineWidth(2);
    fSelectedPixelLocal->Draw("SAME");
    selectedPixelLegend->AddEntry(fSelectedPixelLocal, "Initial local Gaussian", "l");
  }
  selectedPixelLegend->Draw();
  if (selectedPixelLocalValid && fSelectedPixelLocal->GetNDF() > 0) AddFitResultsToStatsBox(hSelectedPixelDt, selectedPixelLocalMean, fSelectedPixelLocal->GetParError(1), selectedPixelLocalSigma, fSelectedPixelLocal->GetParError(2), fSelectedPixelLocal->GetChisquare()/fSelectedPixelLocal->GetNDF());

  TCanvas *cPassingTiming = new TCanvas(uniqueName("cCDetECalCutStudyPassingTiming"), "CDet timing after ECal-CDet delta-t cut", 1500, 500);
  cPassingTiming->Divide(3, 1);
  cPassingTiming->cd(1);
  hLePass->Draw("HIST");
  cPassingTiming->cd(2);
  hTePass->Draw("HIST");
  cPassingTiming->cd(3);
  hTotPass->Draw("HIST");

  TCanvas *cSelectedPixel = new TCanvas(uniqueName("cCDetECalCutStudySelectedPixel"), TString::Format("CDet logical pixel ID %d after delta-t cut", logicalPixelID), 900, 700);
  hSelectedPixelLe->Draw("HIST");

  TCanvas *cECalVsCDetPass = new TCanvas(uniqueName("cCDetECalCutStudyECalVsCDetPass"), "ECal time versus CDet LE after delta-t cut", 900, 700);
  hECalVsCDetPass->Draw("COLZ");

  TCanvas *cECalVsSelectedPixelPass = new TCanvas(uniqueName("cCDetECalCutStudyECalVsSelectedPixelPass"), TString::Format("ECal time versus CDet logical pixel ID %d after delta-t cut", logicalPixelID), 900, 700);
  hECalVsSelectedPixelPass->Draw("COLZ");

  // TCanvas *cEfficiency = new TCanvas(uniqueName("cCDetECalCutStudyEfficiency"), "ECal-CDet delta-t-cut efficiency", 1100, 700);
  // gEfficiency->Draw("AP");
  // gEfficiency->GetXaxis()->SetLimits(-0.5, NumCDetPaddles - 0.5);
  // if (IsCrossTargetRun(gRunNumber)) {
  //   TPaveText *crossTargetNote = new TPaveText(0.14, 0.82, 0.39, 0.89, "NDC");
  //   crossTargetNote->SetFillColor(0);
  //   crossTargetNote->SetBorderSize(1);
  //   crossTargetNote->SetTextAlign(22);
  //   crossTargetNote->AddText(TString::Format("Cross-target run %d", gRunNumber));
  //   crossTargetNote->Draw();
  // }

  if (hDtVsTot) {
    TCanvas *cDtVsTot = new TCanvas(uniqueName("cCDetECalCutStudyDtVsTot"), TString::Format("Uncut ECal-CDet delta-t versus CDet TOT, logical pixel ID %d", logicalPixelID), 900, 700);
    hDtVsTot->Draw("COLZ");
  }

  // const double overallEfficiency = baselineHitCount > 0 ? static_cast<double>(passingHitCount) / baselineHitCount : 0.0;
  std::cout << "[CDet timing-cut study]\n"
            << "  ECal energy cut: [" << ECalEnergyMin << ", " << ECalEnergyMax << "] GeV\n"
            << "  peak search interval: [" << dtMinCut << ", " << dtMaxCut << "] ns\n"
            << "  broad background fit interval: [" << DiffMin << ", " << DiffMax << "] ns\n"
            << "  local and clean fits use the full peak-search interval; background rejection: " << NReject << " sigma\n"
            << "  events passing ECal energy cut: " << energySelectedEventCount << " / " << nEvents << "\n"
            << "  baseline hits: " << baselineHitCount << "\n"
            << "  passing hits: " << passingHitCount << "\n"
            // << "  overall efficiency: " << overallEfficiency << "\n"
            << "  selected logical pixel ID: " << logicalPixelID << "\n"
            << "  selected-pixel passing hits: " << selectedPixelPassingCount << "\n"
            << "  local peak fit status / mean / sigma: " << selectedPixelLocalStatus << " / " << selectedPixelLocalMean << " / " << selectedPixelLocalSigma << " ns\n"
            << "  background rejection window: [" << selectedPixelRejectLow << ", " << selectedPixelRejectHigh << "] ns (NReject=" << NReject << ")\n"
            << "  broad background fit status: " << selectedPixelBackgroundStatus;
  if (selectedPixelBackgroundValid) std::cout << ", amplitude / mean / sigma: " << fSelectedPixelBackground->GetParameter(0) << " / " << fSelectedPixelBackground->GetParameter(1) << " / " << std::fabs(fSelectedPixelBackground->GetParameter(2));
  std::cout << "\n  final cleaned-peak fit status: " << selectedPixelCleanStatus;
  if (selectedPixelCleanValid) std::cout << ", mean / sigma: " << fSelectedPixelClean->GetParameter(1) << " +/- " << fSelectedPixelClean->GetParError(1) << " / " << std::fabs(fSelectedPixelClean->GetParameter(2)) << " ns";
  std::cout << "\n";
}

void extractCDetBarPixelTimingOffsets(int pixelBase = 480, double Width = 1.0,
                                      double HistMin = -60.0, double HistMax = 30.0,
                                      double FitMin = -30.0, double FitMax = 0.0,
                                      int minEntries = 100, double minSigma = 0.5, double maxSigma = 20.0,
                                      double maxChi2Ndf = 10.0, double centroidEdgeMargin = 1.0,
                                      bool saveFitCanvases = false, TString fitCanvasDir = "CDetPixelTimingFits",
                                      bool saveCandidateTable = false,
                                      TString candidateOutput = "CDet_pixel_timing_offsets_candidate.dat",
                                      double ECalEnergyMin = 1.0, double ECalEnergyMax = 12.0,
                                      double TotMin = 0.0, double TotMax = 80.0,
                                      double localFitHalfWidth = 8.0, double NReject = 2.5,
                                      double PeakSeedMin = -25.0, double PeakSeedMax = -5.0,
                                      bool makeProjectionComparison = true,
                                      double AcceptedTotMin = gCDetDiagnosticAcceptedTotMin,
                                      double AcceptedTotMax = gCDetDiagnosticAcceptedTotMax,
                                      TString pixelCutFile = "CDet_pixel_quality_cuts.root") {
  TH1::AddDirectory(kFALSE);
  (void)localFitHalfWidth; // Retained for positional compatibility; extraction fits now use FitMin-FitMax exactly.

  if (pixelBase < 0 || pixelBase >= NumCDetPaddles) {
    std::cerr << "[CDet pixel timing] ERROR: requested logical pixel ID " << pixelBase << " is outside [0, " << NumCDetPaddles - 1 << "].\n";
    return;
  }
  if (Width <= 0.0 || HistMin >= HistMax || FitMin >= FitMax || FitMin < HistMin || FitMax > HistMax || PeakSeedMin >= PeakSeedMax || PeakSeedMin < FitMin || PeakSeedMax > FitMax || minEntries < 1 || minSigma <= 0.0 || minSigma >= maxSigma || maxChi2Ndf <= 0.0 || centroidEdgeMargin < 0.0 || ECalEnergyMin >= ECalEnergyMax || TotMin >= TotMax || AcceptedTotMin >= AcceptedTotMax || NReject <= 0.0) {
    std::cerr << "[CDet pixel timing] ERROR: invalid histogram, fit, or quality-limit argument.\n";
    return;
  }

  const int requestedPixel = pixelBase;
  pixelBase = (pixelBase/NumPaddles)*NumPaddles;
  const int bar = pixelBase/NumPaddles;
  const int nBins = (int)((HistMax - HistMin)/Width);
  const int nTotBins = (int)((TotMax - TotMin)/Width);
  if (nBins < 1 || nTotBins < 1) {
    std::cerr << "[CDet pixel timing] ERROR: histogram binning produces fewer than one bin.\n";
    return;
  }

  const size_t nEvents = vGoodLe.size();
  if (vGoodID.size() != nEvents || vGoodTot.size() != nEvents || v_GoodECalAdcTime.size() != nEvents || v_GoodECalE.size() != nEvents) {
    std::cerr << "[CDet pixel timing] ERROR: event vectors are not aligned: LE=" << vGoodLe.size() << ", TOT=" << vGoodTot.size() << ", ID=" << vGoodID.size() << ", ECal time=" << v_GoodECalAdcTime.size() << ", ECal energy=" << v_GoodECalE.size() << ".\n";
    return;
  }
  for (size_t ev = 0; ev < nEvents; ++ev) {
    if (vGoodID[ev].size() != vGoodLe[ev].size() || vGoodTot[ev].size() != vGoodLe[ev].size()) {
      std::cerr << "[CDet pixel timing] ERROR: LE, TOT, and ID vectors differ in event " << ev << ": LE=" << vGoodLe[ev].size() << ", TOT=" << vGoodTot[ev].size() << ", ID=" << vGoodID[ev].size() << ".\n";
      return;
    }
  }

  static unsigned long invocation = 0;
  const unsigned long tag = ++invocation;
  auto uniqueName = [tag](const char *base) { return TString::Format("%s_%lu", base, tag); };

  std::vector<TH1D*> hPixelDt(NumPaddles, nullptr);
  std::vector<TH1D*> hPixelDtClean(NumPaddles, nullptr);
  std::vector<TH2D*> hPixelDtVsTot(NumPaddles, nullptr);
  std::vector<TF1*> fPixelLocal(NumPaddles, nullptr);
  std::vector<TF1*> fPixelBackgroundReject(NumPaddles, nullptr);
  std::vector<TF1*> fPixelBackground(NumPaddles, nullptr);
  std::vector<TF1*> fPixelClean(NumPaddles, nullptr);
  std::vector<int> fitEntries(NumPaddles, 0);
  std::vector<int> fitStatus(NumPaddles, -1);
  std::vector<int> validityCode(NumPaddles, 0);
  std::vector<std::string> failureReason(NumPaddles, "not processed");
  std::vector<double> amplitude(NumPaddles, NAN), amplitudeErr(NumPaddles, NAN);
  std::vector<double> centroid(NumPaddles, NAN), centroidErr(NumPaddles, NAN);
  std::vector<double> sigma(NumPaddles, NAN), sigmaErr(NumPaddles, NAN);
  std::vector<double> backgroundAmplitude(NumPaddles, NAN), backgroundMean(NumPaddles, NAN), backgroundSigma(NumPaddles, NAN);
  std::vector<double> backgroundAmplitudeLimit(NumPaddles, NAN);
  std::vector<double> rejectLow(NumPaddles, NAN), rejectHigh(NumPaddles, NAN);
  std::vector<double> chi2(NumPaddles, NAN), chi2Ndf(NumPaddles, NAN);
  std::vector<int> ndf(NumPaddles, 0);
  std::vector<bool> validLocalFit(NumPaddles, false);
  std::vector<bool> validCleanFit(NumPaddles, false);
  std::vector<bool> validFit(NumPaddles, false);
  std::vector<double> correction(NumPaddles, NAN), correctionErr(NumPaddles, NAN);

  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    const int pixelID = pixelBase + localPixel;
    hPixelDt[localPixel] = new TH1D(uniqueName(TString::Format("hCDetPixelTimingDt_%d", pixelID)), TString::Format("After ECal energy cut, logical pixel ID %d;t_{ECal}-t_{CDet,LE} (ns);Baseline hits", pixelID), nBins, HistMin, HistMax);
    hPixelDtVsTot[localPixel] = new TH2D(uniqueName(TString::Format("hCDetPixelTimingDtVsTot_%d", pixelID)), TString::Format("After ECal energy cut, logical pixel ID %d;CDet TOT (ns);t_{ECal}-t_{CDet,LE} (ns)", pixelID), nTotBins, TotMin, TotMax, nBins, HistMin, HistMax);
  }
  TH1D *hBarDt = new TH1D(uniqueName("hCDetBarTimingDt"),
                          TString::Format("After ECal energy cut, CDet bar %d (all instrumented pixels);t_{ECal}-t_{CDet,LE} (ns);Baseline hits", bar),
                          nBins, HistMin, HistMax);
  TH2D *hBarDtVsTot = new TH2D(uniqueName("hCDetBarTimingDtVsTot"),
                               TString::Format("After ECal energy cut, CDet bar %d (all instrumented pixels);CDet TOT (ns);t_{ECal}-t_{CDet,LE} (ns)", bar),
                               nTotBins, TotMin, TotMax, nBins, HistMin, HistMax);
  TH1D *hBarProjectedDt = new TH1D(uniqueName("hCDetBarProjectedTimingDt"),
                                   TString::Format("ECal projection in CDet bar %d;t_{ECal}-t_{CDet,LE} (ns);Trajectory-matched hits", bar),
                                   nBins, HistMin, HistMax);
  TH2D *hBarProjectedDtVsTot = new TH2D(uniqueName("hCDetBarProjectedTimingDtVsTot"),
                                        TString::Format("ECal projection in CDet bar %d;CDet TOT (ns);t_{ECal}-t_{CDet,LE} (ns)", bar),
                                        nTotBins, TotMin, TotMax, nBins, HistMin, HistMax);
  TH1D *hBarProjectedQualityDt = new TH1D(uniqueName("hCDetBarProjectedQualityTimingDt"),
                                          TString::Format("ECal projection + TOT/polygon selection, CDet bar %d;t_{ECal}-t_{CDet,LE} (ns);Selected hits", bar),
                                          nBins, HistMin, HistMax);
  TH2D *hBarProjectedQualityDtVsTot = new TH2D(uniqueName("hCDetBarProjectedQualityTimingDtVsTot"),
                                               TString::Format("ECal projection + TOT/polygon selection, CDet bar %d;CDet TOT (ns);t_{ECal}-t_{CDet,LE} (ns)", bar),
                                               nTotBins, TotMin, TotMax, nBins, HistMin, HistMax);
  TH1D *hValidity = new TH1D(uniqueName("hCDetBarPixelValidity"), "Pixel fit status;CDet logical pixel ID;Status code", NumPaddles, pixelBase - 0.5, pixelBase + NumPaddles - 0.5);
  TH1D *hCentroidDistribution = new TH1D(uniqueName("hCDetBarPixelCentroidDistribution"), "Reliable signal-fit centroids;#mu_{i} (ns);Pixels", nBins, HistMin, HistMax);
  TH1D *hPredictedCentroidDistribution = new TH1D(uniqueName("hCDetBarPixelPredictedCentroidDistribution"), "Predicted centroids after candidate corrections;#mu_{i}-c_{i} (ns);Pixels", nBins, HistMin, HistMax);

  double barXMin = std::numeric_limits<double>::infinity();
  double barXMax = -std::numeric_limits<double>::infinity();
  double barYSum = 0.0;
  size_t barCoordinateCount = 0;
  const bool projectionVectorsAligned =
      vCDetGoodX.size() == nEvents && vCDetGoodY.size() == nEvents && vCDetGoodZ.size() == nEvents &&
      v_GoodECalX.size() == nEvents && v_GoodECalY.size() == nEvents;
  if (makeProjectionComparison && projectionVectorsAligned) {
    for (size_t ev = 0; ev < nEvents; ++ev) {
      const size_t nCoordinates = std::min(vGoodID[ev].size(),
          std::min(vCDetGoodX[ev].size(), std::min(vCDetGoodY[ev].size(), vCDetGoodZ[ev].size())));
      for (size_t ihit = 0; ihit < nCoordinates; ++ihit) {
        const int pixelID = vGoodID[ev][ihit];
        if (pixelID < pixelBase || pixelID >= pixelBase + NumPaddles) continue;
        barXMin = std::min(barXMin, vCDetGoodX[ev][ihit]);
        barXMax = std::max(barXMax, vCDetGoodX[ev][ihit]);
        barYSum += vCDetGoodY[ev][ihit];
        ++barCoordinateCount;
      }
    }
  }
  const bool projectionAvailable = makeProjectionComparison && projectionVectorsAligned && barCoordinateCount > 0;
  const double barYCenter = projectionAvailable ? barYSum/barCoordinateCount : NAN;
  const double halfPixelPitch = 0.5*0.00525*0.5*(XCorr1 + XCorr2);

  std::vector<TCutG*> pixelLeTotCut(NumPaddles, nullptr);
  int loadedBarPixelCuts = 0;
  if (!pixelCutFile.IsNull() && !gSystem->AccessPathName(pixelCutFile)) {
    TFile cutInput(pixelCutFile, "READ");
    if (!cutInput.IsZombie()) {
      for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
        pixelLeTotCut[localPixel] = LoadCDetPixelLeTotCut(cutInput, pixelBase + localPixel);
        if (pixelLeTotCut[localPixel]) ++loadedBarPixelCuts;
      }
    }
    cutInput.Close();
  }

  size_t energySelectedEventCount = 0;
  size_t projectedBarHitCount = 0;
  size_t projectedQualityBarHitCount = 0;
  for (size_t ev = 0; ev < nEvents; ++ev) {
    const double ecalEnergy = v_GoodECalE[ev];
    if (ecalEnergy < ECalEnergyMin || ecalEnergy > ECalEnergyMax) continue;
    ++energySelectedEventCount;
    const double tECal = v_GoodECalAdcTime[ev];
    for (size_t ihit = 0; ihit < vGoodLe[ev].size(); ++ihit) {
      const int pixelID = vGoodID[ev][ihit];
      if (pixelID < pixelBase || pixelID >= pixelBase + NumPaddles) continue;
      const int localPixel = pixelID - pixelBase;
      const double dt = tECal - vGoodLe[ev][ihit];
      hPixelDt[localPixel]->Fill(dt);
      hPixelDtVsTot[localPixel]->Fill(vGoodTot[ev][ihit], dt);
      if (!IsUnusedPixel(pixelID)) {
        hBarDt->Fill(dt);
        hBarDtVsTot->Fill(vGoodTot[ev][ihit], dt);
        if (projectionAvailable && ihit < vCDetGoodX[ev].size() && ihit < vCDetGoodY[ev].size() && ihit < vCDetGoodZ[ev].size()) {
          const double hitZ = vCDetGoodZ[ev][ihit];
          const double xProjection = v_GoodECalX[ev]*hitZ/ECal_dist;
          const double yProjection = v_GoodECalY[ev]*hitZ/ECal_dist;
          const bool projectionInBar =
              xProjection >= barXMin - halfPixelPitch && xProjection <= barXMax + halfPixelPitch &&
              std::fabs(yProjection - barYCenter) <= CDet_y_half_length;
          if (projectionInBar) {
            hBarProjectedDt->Fill(dt);
            hBarProjectedDtVsTot->Fill(vGoodTot[ev][ihit], dt);
            ++projectedBarHitCount;
            const double tot = vGoodTot[ev][ihit];
            const double le = vGoodLe[ev][ihit];
            const bool passesQualitySelection = pixelLeTotCut[localPixel]
                ? pixelLeTotCut[localPixel]->IsInside(tot, le)
                : (tot > AcceptedTotMin && tot < AcceptedTotMax);
            if (passesQualitySelection) {
              hBarProjectedQualityDt->Fill(dt);
              hBarProjectedQualityDtVsTot->Fill(tot, dt);
              ++projectedQualityBarHitCount;
            }
          }
        }
      }
    }
  }

  TF1 *fBarLocal = nullptr;
  bool validBarLocalFit = false;
  int barLocalFitStatus = -1;
  const int barFitBinMin = hBarDt->FindBin(FitMin);
  const int barFitBinMax = hBarDt->FindBin(FitMax);
  const int barFitEntries = (int)hBarDt->Integral(barFitBinMin, barFitBinMax);
  if (barFitEntries >= minEntries) {
    int peakBinMin = hBarDt->FindBin(PeakSeedMin);
    int peakBinMax = hBarDt->FindBin(PeakSeedMax);
    int peakBin = peakBinMin;
    for (int bin = peakBinMin + 1; bin <= peakBinMax; ++bin) {
      if (hBarDt->GetBinContent(bin) > hBarDt->GetBinContent(peakBin)) peakBin = bin;
    }
    if (hBarDt->GetBinContent(peakBin) <= 0.0) {
      peakBin = barFitBinMin;
      for (int bin = barFitBinMin + 1; bin <= barFitBinMax; ++bin) {
        if (hBarDt->GetBinContent(bin) > hBarDt->GetBinContent(peakBin)) peakBin = bin;
      }
    }
    const double peak = hBarDt->GetBinCenter(peakBin);
    const double edgeBackground = 0.5*(hBarDt->GetBinContent(hBarDt->FindBin(FitMin)) + hBarDt->GetBinContent(hBarDt->FindBin(FitMax)));
    fBarLocal = new TF1(uniqueName("fCDetBarTimingLocal"), "gaus(0)+pol1(3)", FitMin, FitMax);
    fBarLocal->SetParameters(std::max(1.0, hBarDt->GetBinContent(peakBin) - edgeBackground), peak,
                             std::max(Width, (FitMax - FitMin)/10.0), edgeBackground, 0.0);
    fBarLocal->SetParNames("Bar amplitude", "Bar mean", "Bar sigma", "Background intercept", "Background slope");
    barLocalFitStatus = hBarDt->Fit(fBarLocal, "RQN0");
    const double barMean = fBarLocal->GetParameter(1);
    const double barSigma = std::fabs(fBarLocal->GetParameter(2));
    validBarLocalFit = barLocalFitStatus == 0 &&
                       std::isfinite(fBarLocal->GetParameter(0)) && fBarLocal->GetParameter(0) > 0.0 &&
                       std::isfinite(barMean) && barMean > FitMin && barMean < FitMax &&
                       std::isfinite(barSigma) && barSigma >= minSigma && barSigma <= maxSigma;
  }

  TF1 *fBarProjectedLocal = nullptr;
  bool validBarProjectedLocalFit = false;
  int barProjectedLocalFitStatus = -1;
  const int barProjectedFitEntries = (int)hBarProjectedDt->Integral(
      hBarProjectedDt->FindBin(FitMin), hBarProjectedDt->FindBin(FitMax));
  if (barProjectedFitEntries >= minEntries) {
    int peakBinMin = hBarProjectedDt->FindBin(PeakSeedMin);
    int peakBinMax = hBarProjectedDt->FindBin(PeakSeedMax);
    int peakBin = peakBinMin;
    for (int bin = peakBinMin + 1; bin <= peakBinMax; ++bin) {
      if (hBarProjectedDt->GetBinContent(bin) > hBarProjectedDt->GetBinContent(peakBin)) peakBin = bin;
    }
    if (hBarProjectedDt->GetBinContent(peakBin) <= 0.0) {
      peakBin = hBarProjectedDt->FindBin(FitMin);
      for (int bin = peakBin + 1; bin <= hBarProjectedDt->FindBin(FitMax); ++bin) {
        if (hBarProjectedDt->GetBinContent(bin) > hBarProjectedDt->GetBinContent(peakBin)) peakBin = bin;
      }
    }
    const double peak = hBarProjectedDt->GetBinCenter(peakBin);
    const double edgeBackground = 0.5*(hBarProjectedDt->GetBinContent(hBarProjectedDt->FindBin(FitMin)) +
                                       hBarProjectedDt->GetBinContent(hBarProjectedDt->FindBin(FitMax)));
    fBarProjectedLocal = new TF1(uniqueName("fCDetBarProjectedTimingLocal"), "gaus(0)+pol1(3)", FitMin, FitMax);
    fBarProjectedLocal->SetParameters(
        std::max(1.0, hBarProjectedDt->GetBinContent(peakBin) - edgeBackground), peak,
        std::max(Width, (FitMax - FitMin)/10.0), edgeBackground, 0.0);
    fBarProjectedLocal->SetParNames("Projected amplitude", "Projected mean", "Projected sigma", "Background intercept", "Background slope");
    barProjectedLocalFitStatus = hBarProjectedDt->Fit(fBarProjectedLocal, "RQN0");
    const double projectedMean = fBarProjectedLocal->GetParameter(1);
    const double projectedSigma = std::fabs(fBarProjectedLocal->GetParameter(2));
    validBarProjectedLocalFit = barProjectedLocalFitStatus == 0 &&
        std::isfinite(fBarProjectedLocal->GetParameter(0)) && fBarProjectedLocal->GetParameter(0) > 0.0 &&
        std::isfinite(projectedMean) && projectedMean > FitMin && projectedMean < FitMax &&
        std::isfinite(projectedSigma) && projectedSigma >= minSigma && projectedSigma <= maxSigma;
  }

  TF1 *fBarProjectedQualityLocal = nullptr;
  bool validBarProjectedQualityLocalFit = false;
  int barProjectedQualityLocalFitStatus = -1;
  const int barProjectedQualityFitEntries = (int)hBarProjectedQualityDt->Integral(
      hBarProjectedQualityDt->FindBin(FitMin), hBarProjectedQualityDt->FindBin(FitMax));
  if (barProjectedQualityFitEntries >= minEntries) {
    int peakBinMin = hBarProjectedQualityDt->FindBin(PeakSeedMin);
    int peakBinMax = hBarProjectedQualityDt->FindBin(PeakSeedMax);
    int peakBin = peakBinMin;
    for (int bin = peakBinMin + 1; bin <= peakBinMax; ++bin) {
      if (hBarProjectedQualityDt->GetBinContent(bin) > hBarProjectedQualityDt->GetBinContent(peakBin)) peakBin = bin;
    }
    if (hBarProjectedQualityDt->GetBinContent(peakBin) <= 0.0) {
      peakBin = hBarProjectedQualityDt->FindBin(FitMin);
      for (int bin = peakBin + 1; bin <= hBarProjectedQualityDt->FindBin(FitMax); ++bin) {
        if (hBarProjectedQualityDt->GetBinContent(bin) > hBarProjectedQualityDt->GetBinContent(peakBin)) peakBin = bin;
      }
    }
    const double peak = hBarProjectedQualityDt->GetBinCenter(peakBin);
    const double edgeBackground = 0.5*(
        hBarProjectedQualityDt->GetBinContent(hBarProjectedQualityDt->FindBin(FitMin)) +
        hBarProjectedQualityDt->GetBinContent(hBarProjectedQualityDt->FindBin(FitMax)));
    fBarProjectedQualityLocal = new TF1(uniqueName("fCDetBarProjectedQualityTimingLocal"),
                                        "gaus(0)+pol1(3)", FitMin, FitMax);
    fBarProjectedQualityLocal->SetParameters(
        std::max(1.0, hBarProjectedQualityDt->GetBinContent(peakBin) - edgeBackground), peak,
        std::max(Width, (FitMax - FitMin)/10.0), edgeBackground, 0.0);
    fBarProjectedQualityLocal->SetParNames("Selected amplitude", "Selected mean", "Selected sigma",
                                           "Background intercept", "Background slope");
    barProjectedQualityLocalFitStatus = hBarProjectedQualityDt->Fit(fBarProjectedQualityLocal, "RQN0");
    const double selectedMean = fBarProjectedQualityLocal->GetParameter(1);
    const double selectedSigma = std::fabs(fBarProjectedQualityLocal->GetParameter(2));
    validBarProjectedQualityLocalFit = barProjectedQualityLocalFitStatus == 0 &&
        std::isfinite(fBarProjectedQualityLocal->GetParameter(0)) && fBarProjectedQualityLocal->GetParameter(0) > 0.0 &&
        std::isfinite(selectedMean) && selectedMean > FitMin && selectedMean < FitMax &&
        std::isfinite(selectedSigma) && selectedSigma >= minSigma && selectedSigma <= maxSigma;
  }

  int instrumentedPixels = 0;
  int sufficientStatistics = 0;
  int successfulFits = 0;
  int lowStatisticsCount = 0;
  int rootFitFailureCount = 0;
  int invalidParameterCount = 0;
  int boundaryCount = 0;
  int sigmaCount = 0;
  int amplitudeCount = 0;
  int ndfCount = 0;
  int chi2Count = 0;
  int peakSeedFallbackCount = 0;

  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    const int pixelID = pixelBase + localPixel;
    if (IsUnusedPixel(pixelID)) {
      validityCode[localPixel] = 1;
      failureReason[localPixel] = "unused pixel";
      continue;
    }
    ++instrumentedPixels;

    const int fitBinMin = hPixelDt[localPixel]->FindBin(FitMin);
    const int fitBinMax = hPixelDt[localPixel]->FindBin(FitMax);
    fitEntries[localPixel] = (int)hPixelDt[localPixel]->Integral(fitBinMin, fitBinMax);
    if (fitEntries[localPixel] < minEntries) {
      validityCode[localPixel] = 2;
      failureReason[localPixel] = "insufficient entries";
      ++lowStatisticsCount;
      continue;
    }
    ++sufficientStatistics;

    int peakBinMin = hPixelDt[localPixel]->FindBin(PeakSeedMin);
    int peakBinMax = hPixelDt[localPixel]->FindBin(PeakSeedMax);
    int peakBin = peakBinMin;
    for (int binIndex = peakBinMin + 1; binIndex <= peakBinMax; ++binIndex) {
      if (hPixelDt[localPixel]->GetBinContent(binIndex) > hPixelDt[localPixel]->GetBinContent(peakBin)) peakBin = binIndex;
    }
    if (hPixelDt[localPixel]->GetBinContent(peakBin) <= 0.0) {
      peakBinMin = fitBinMin;
      peakBinMax = fitBinMax;
      peakBin = peakBinMin;
      for (int binIndex = peakBinMin + 1; binIndex <= peakBinMax; ++binIndex) {
        if (hPixelDt[localPixel]->GetBinContent(binIndex) > hPixelDt[localPixel]->GetBinContent(peakBin)) peakBin = binIndex;
      }
      ++peakSeedFallbackCount;
    }
    const double peak = hPixelDt[localPixel]->GetBinCenter(peakBin);
    const double localFitMin = FitMin;
    const double localFitMax = FitMax;
    const double localBackground = 0.5*(hPixelDt[localPixel]->GetBinContent(hPixelDt[localPixel]->FindBin(localFitMin)) + hPixelDt[localPixel]->GetBinContent(hPixelDt[localPixel]->FindBin(localFitMax)));
    fPixelLocal[localPixel] = new TF1(uniqueName(TString::Format("fCDetPixelTimingLocal_%d", pixelID)), "gaus(0)+pol1(3)", localFitMin, localFitMax);
    fPixelLocal[localPixel]->SetParameters(std::max(1.0, hPixelDt[localPixel]->GetBinContent(peakBin) - localBackground), peak, std::max(Width, (FitMax - FitMin)/10.0), localBackground, 0.0);
    fPixelLocal[localPixel]->SetParNames("Local amplitude", "Local mean", "Local sigma", "Local background intercept", "Local background slope");
    fitStatus[localPixel] = hPixelDt[localPixel]->Fit(fPixelLocal[localPixel], "RQN0");
    const double localMean = fPixelLocal[localPixel]->GetParameter(1);
    const double localSigma = std::fabs(fPixelLocal[localPixel]->GetParameter(2));
    if (fitStatus[localPixel] != 0 || !std::isfinite(fPixelLocal[localPixel]->GetParameter(0)) || fPixelLocal[localPixel]->GetParameter(0) <= 0.0 || !std::isfinite(localMean) || !std::isfinite(localSigma) || !std::isfinite(fPixelLocal[localPixel]->GetParameter(3)) || !std::isfinite(fPixelLocal[localPixel]->GetParameter(4)) || localSigma < minSigma || localSigma > maxSigma || localMean <= FitMin || localMean >= FitMax) {
      validityCode[localPixel] = 3;
      failureReason[localPixel] = "invalid local signal fit";
      ++rootFitFailureCount;
      continue;
    }
    validLocalFit[localPixel] = true;

    amplitude[localPixel] = fPixelLocal[localPixel]->GetParameter(0);
    amplitudeErr[localPixel] = fPixelLocal[localPixel]->GetParError(0);
    centroid[localPixel] = localMean;
    centroidErr[localPixel] = fPixelLocal[localPixel]->GetParError(1);
    sigma[localPixel] = localSigma;
    sigmaErr[localPixel] = fPixelLocal[localPixel]->GetParError(2);
    chi2[localPixel] = fPixelLocal[localPixel]->GetChisquare();
    ndf[localPixel] = fPixelLocal[localPixel]->GetNDF();
    chi2Ndf[localPixel] = ndf[localPixel] > 0 ? chi2[localPixel]/ndf[localPixel] : NAN;

    ++successfulFits;
    if (!std::isfinite(amplitude[localPixel]) || !std::isfinite(amplitudeErr[localPixel]) || !std::isfinite(centroid[localPixel]) || !std::isfinite(centroidErr[localPixel]) || !std::isfinite(sigma[localPixel]) || !std::isfinite(sigmaErr[localPixel]) || !std::isfinite(chi2[localPixel]) || centroidErr[localPixel] <= 0.0 || sigmaErr[localPixel] <= 0.0) {
      validityCode[localPixel] = 4;
      failureReason[localPixel] = "invalid parameter or uncertainty";
      ++invalidParameterCount;
      continue;
    }
    if (centroid[localPixel] - FitMin <= centroidEdgeMargin || FitMax - centroid[localPixel] <= centroidEdgeMargin) {
      validityCode[localPixel] = 5;
      failureReason[localPixel] = "centroid near fit boundary";
      ++boundaryCount;
      continue;
    }
    if (sigma[localPixel] < minSigma || sigma[localPixel] > maxSigma) {
      validityCode[localPixel] = 6;
      failureReason[localPixel] = "sigma outside limits";
      ++sigmaCount;
      continue;
    }
    if (amplitude[localPixel] <= 0.0) {
      validityCode[localPixel] = 7;
      failureReason[localPixel] = "nonpositive Gaussian amplitude";
      ++amplitudeCount;
      continue;
    }
    if (ndf[localPixel] <= 0) {
      validityCode[localPixel] = 8;
      failureReason[localPixel] = "invalid NDF";
      ++ndfCount;
      continue;
    }
    if (!std::isfinite(chi2Ndf[localPixel]) || chi2Ndf[localPixel] > maxChi2Ndf) {
      validityCode[localPixel] = 9;
      failureReason[localPixel] = "chi2/NDF outside limit";
      ++chi2Count;
      continue;
    }

    validityCode[localPixel] = 10;
    failureReason[localPixel] = "valid";
    validFit[localPixel] = true;

    // Background-subtracted diagnostics removed; calibration uses fPixelLocal.
#if 0
    // Calibration centroids and the weighted bar reference come from fPixelLocal.
    rejectLow[localPixel] = std::max(HistMin, localMean - NReject*localSigma);
    rejectHigh[localPixel] = std::min(HistMax, localMean + NReject*localSigma);
    CDetTimingBackgroundRejectLow = rejectLow[localPixel];
    CDetTimingBackgroundRejectHigh = rejectHigh[localPixel];
    double maxSidebandBinContent = 0.0;
    int leftSidebandBins = 0;
    int rightSidebandBins = 0;
    for (int bin = 1; bin <= hPixelDt[localPixel]->GetNbinsX(); ++bin) {
      const double binCenter = hPixelDt[localPixel]->GetBinCenter(bin);
      if (binCenter < HistMin || binCenter > HistMax || (rejectLow[localPixel] <= binCenter && binCenter <= rejectHigh[localPixel])) continue;
      maxSidebandBinContent = std::max(maxSidebandBinContent, hPixelDt[localPixel]->GetBinContent(bin));
      if (binCenter < rejectLow[localPixel]) ++leftSidebandBins;
      if (binCenter > rejectHigh[localPixel]) ++rightSidebandBins;
    }
    backgroundAmplitudeLimit[localPixel] = maxSidebandBinContent;
    if (maxSidebandBinContent <= 0.0 || leftSidebandBins < 3 || rightSidebandBins < 3) {
      continue;
    }
    fPixelBackgroundReject[localPixel] = new TF1(uniqueName(TString::Format("fCDetPixelTimingBackgroundReject_%d", pixelID)), CDetTimingBackgroundGaussianReject, HistMin, HistMax, 3);
    fPixelBackgroundReject[localPixel]->SetParameters(0.8*maxSidebandBinContent, 0.5*(HistMin + HistMax), std::max(Width, (HistMax - HistMin)/3.0));
    fPixelBackgroundReject[localPixel]->SetParLimits(0, 0.0, maxSidebandBinContent);
    fPixelBackgroundReject[localPixel]->SetParLimits(1, HistMin, HistMax);
    fPixelBackgroundReject[localPixel]->SetParLimits(2, Width, HistMax - HistMin);
    const int backgroundFitStatus = hPixelDt[localPixel]->Fit(fPixelBackgroundReject[localPixel], "RQN0");
    backgroundAmplitude[localPixel] = fPixelBackgroundReject[localPixel]->GetParameter(0);
    backgroundMean[localPixel] = fPixelBackgroundReject[localPixel]->GetParameter(1);
    backgroundSigma[localPixel] = std::fabs(fPixelBackgroundReject[localPixel]->GetParameter(2));
    if (backgroundFitStatus != 0 || !std::isfinite(backgroundAmplitude[localPixel]) || !std::isfinite(backgroundMean[localPixel]) || !std::isfinite(backgroundSigma[localPixel]) || backgroundAmplitude[localPixel] < 0.0 || backgroundSigma[localPixel] <= 0.0) {
      continue;
    }

    fPixelBackground[localPixel] = new TF1(uniqueName(TString::Format("fCDetPixelTimingBackground_%d", pixelID)), "gaus", HistMin, HistMax);
    fPixelBackground[localPixel]->SetParameters(backgroundAmplitude[localPixel], backgroundMean[localPixel], backgroundSigma[localPixel]);
    hPixelDtClean[localPixel] = (TH1D*)hPixelDt[localPixel]->Clone(uniqueName(TString::Format("hCDetPixelTimingDtClean_%d", pixelID)));
    hPixelDtClean[localPixel]->SetTitle(TString::Format("Background-subtracted, logical pixel ID %d;t_{ECal}-t_{CDet,LE} (ns);Signal estimate", pixelID));
    for (int bin = 1; bin <= hPixelDtClean[localPixel]->GetNbinsX(); ++bin) {
      const double originalError = hPixelDt[localPixel]->GetBinError(bin);
      hPixelDtClean[localPixel]->SetBinContent(bin, hPixelDt[localPixel]->GetBinContent(bin) - fPixelBackground[localPixel]->Eval(hPixelDt[localPixel]->GetBinCenter(bin)));
      hPixelDtClean[localPixel]->SetBinError(bin, originalError);
    }

    const double cleanFitMin = FitMin;
    const double cleanFitMax = FitMax;
    fPixelClean[localPixel] = new TF1(uniqueName(TString::Format("fCDetPixelTimingClean_%d", pixelID)), "gaus", cleanFitMin, cleanFitMax);
    fPixelClean[localPixel]->SetParameters(std::max(1.0, fPixelLocal[localPixel]->GetParameter(0)), localMean, localSigma);
    const int cleanFitStatus = hPixelDtClean[localPixel]->Fit(fPixelClean[localPixel], "RQN0");
    validCleanFit[localPixel] = cleanFitStatus == 0 && std::isfinite(fPixelClean[localPixel]->GetParameter(0)) && std::isfinite(fPixelClean[localPixel]->GetParameter(1)) && std::isfinite(fPixelClean[localPixel]->GetParameter(2));
#endif
  }

  double weightSum = 0.0;
  double weightedCentroidSum = 0.0;
  int referencePixels = 0;
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    if (!validFit[localPixel]) continue;
    const double weight = 1.0/(centroidErr[localPixel]*centroidErr[localPixel]);
    weightSum += weight;
    weightedCentroidSum += centroid[localPixel]*weight;
    ++referencePixels;
  }

  const bool referenceValid = referencePixels >= 2 && weightSum > 0.0;
  const double referenceCentroid = referenceValid ? weightedCentroidSum/weightSum : NAN;
  const double referenceCentroidErr = referenceValid ? std::sqrt(1.0/weightSum) : NAN;
  if (referenceValid) {
    for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
      if (!validFit[localPixel]) continue;
      correction[localPixel] = centroid[localPixel] - referenceCentroid;
      const double correlatedVariance = centroidErr[localPixel]*centroidErr[localPixel] - referenceCentroidErr*referenceCentroidErr;
      correctionErr[localPixel] = correlatedVariance > 0.0 ? std::sqrt(correlatedVariance) : std::sqrt(centroidErr[localPixel]*centroidErr[localPixel] + referenceCentroidErr*referenceCentroidErr);
    }
  }

  TGraphErrors *gCentroid = new TGraphErrors();
  TGraphErrors *gCorrection = new TGraphErrors();
  TGraphErrors *gSigma = new TGraphErrors();
  TGraph *gEntries = new TGraph();
  TGraph *gChi2Ndf = new TGraph();
  gCentroid->SetName(uniqueName("gCDetBarPixelCentroid"));
  gCorrection->SetName(uniqueName("gCDetBarPixelCorrection"));
  gSigma->SetName(uniqueName("gCDetBarPixelSigma"));
  gEntries->SetName(uniqueName("gCDetBarPixelEntries"));
  gChi2Ndf->SetName(uniqueName("gCDetBarPixelChi2Ndf"));
  gCentroid->SetTitle("Reliable signal+background fitted centroids;CDet logical pixel ID;#mu_{i} (ns)");
  gCorrection->SetTitle("Candidate additive pixel corrections;CDet logical pixel ID;c_{i}=#mu_{i}-#mu_{0} (ns)");
  gSigma->SetTitle("Reliable fitted Gaussian widths;CDet logical pixel ID;Gaussian #sigma_{i} (ns)");
  gEntries->SetTitle("Fit-region entries;CDet logical pixel ID;Entries");
  gChi2Ndf->SetTitle("Fit quality;CDet logical pixel ID;#chi^{2}/NDF");
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    const int pixelID = pixelBase + localPixel;
    hValidity->SetBinContent(localPixel + 1, validityCode[localPixel]);
    int point = gEntries->GetN();
    gEntries->SetPoint(point, pixelID, fitEntries[localPixel]);
    if (std::isfinite(chi2Ndf[localPixel])) gChi2Ndf->SetPoint(gChi2Ndf->GetN(), pixelID, chi2Ndf[localPixel]);
    if (!validFit[localPixel]) continue;
    point = gCentroid->GetN();
    gCentroid->SetPoint(point, pixelID, centroid[localPixel]);
    gCentroid->SetPointError(point, 0.0, centroidErr[localPixel]);
    point = gSigma->GetN();
    gSigma->SetPoint(point, pixelID, sigma[localPixel]);
    gSigma->SetPointError(point, 0.0, sigmaErr[localPixel]);
    hCentroidDistribution->Fill(centroid[localPixel]);
    if (referenceValid) {
      point = gCorrection->GetN();
      gCorrection->SetPoint(point, pixelID, correction[localPixel]);
      gCorrection->SetPointError(point, 0.0, correctionErr[localPixel]);
      hPredictedCentroidDistribution->Fill(centroid[localPixel] - correction[localPixel]);
    }
  }

  for (TGraph *graph : {static_cast<TGraph*>(gCentroid), static_cast<TGraph*>(gCorrection), static_cast<TGraph*>(gSigma), gEntries, gChi2Ndf}) {
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.8);
  }

  // Tile the five diagnostic windows on a 16-inch MacBook display.  The dense
  // canvases deliberately open compactly and can be enlarged for inspection.
  TCanvas *cBarFits = new TCanvas(uniqueName("cCDetBarPixelTimingFits"), TString::Format("CDet bar %d pixel timing fits", bar),
                                  10, 35, 540, 420);
  cBarFits->Divide(4, 4, 0.001, 0.001);
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    cBarFits->cd(localPixel + 1);
    hPixelDt[localPixel]->SetStats(kTRUE);
    hPixelDt[localPixel]->SetLineColor(kBlack);
    if (hPixelDtClean[localPixel]) {
      hPixelDt[localPixel]->SetMinimum(std::min(0.0, 1.15*hPixelDtClean[localPixel]->GetMinimum()));
      hPixelDtClean[localPixel]->SetLineColor(kGreen + 2);
      hPixelDtClean[localPixel]->SetLineWidth(2);
    }
    hPixelDt[localPixel]->Draw("HIST");
    if (fPixelLocal[localPixel]) {
      fPixelLocal[localPixel]->SetLineColor(kBlue + 1);
      fPixelLocal[localPixel]->SetLineStyle(2);
      fPixelLocal[localPixel]->Draw("SAME");
    }
    if (fPixelBackground[localPixel]) {
      fPixelBackground[localPixel]->SetLineColor(kMagenta + 2);
      fPixelBackground[localPixel]->SetLineStyle(3);
      fPixelBackground[localPixel]->Draw("SAME");
    }
    if (hPixelDtClean[localPixel]) hPixelDtClean[localPixel]->Draw("HIST SAME");
    if (validCleanFit[localPixel]) {
      fPixelClean[localPixel]->SetLineColor(kRed + 1);
      fPixelClean[localPixel]->SetLineWidth(2);
      fPixelClean[localPixel]->Draw("SAME");
    }
    if (validLocalFit[localPixel] && fPixelLocal[localPixel]->GetNDF() > 0) AddFitResultsToStatsBox(hPixelDt[localPixel], fPixelLocal[localPixel]->GetParameter(1), fPixelLocal[localPixel]->GetParError(1), std::fabs(fPixelLocal[localPixel]->GetParameter(2)), fPixelLocal[localPixel]->GetParError(2), fPixelLocal[localPixel]->GetChisquare()/fPixelLocal[localPixel]->GetNDF());
  }

  #if 0
  TCanvas *cBarCleanFits = new TCanvas(uniqueName("cCDetBarPixelTimingCleanFits"), TString::Format("CDet bar %d background-subtracted pixel timing fits", bar),
                                       570, 35, 540, 420);
  cBarCleanFits->Divide(4, 4, 0.001, 0.001);
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    cBarCleanFits->cd(localPixel + 1);
    if (!hPixelDtClean[localPixel]) {
      TPaveText *cleanFitUnavailable = new TPaveText(0.12, 0.42, 0.88, 0.58, "NDC");
      cleanFitUnavailable->SetFillColor(0);
      cleanFitUnavailable->SetBorderSize(1);
      cleanFitUnavailable->AddText(TString::Format("Logical pixel ID %d", pixelBase + localPixel));
      cleanFitUnavailable->AddText("Background-subtracted histogram unavailable");
      cleanFitUnavailable->Draw();
      continue;
    }
    hPixelDtClean[localPixel]->SetStats(kTRUE);
    hPixelDtClean[localPixel]->SetLineColor(kGreen + 2);
    hPixelDtClean[localPixel]->SetLineWidth(2);
    hPixelDtClean[localPixel]->SetMinimum(std::min(0.0, 1.15*hPixelDtClean[localPixel]->GetMinimum()));
    hPixelDtClean[localPixel]->Draw("HIST");
    if (validCleanFit[localPixel]) {
      fPixelClean[localPixel]->SetLineColor(kRed + 1);
      fPixelClean[localPixel]->SetLineWidth(2);
      fPixelClean[localPixel]->Draw("SAME");
      if (fPixelClean[localPixel]->GetNDF() > 0) AddFitResultsToStatsBox(hPixelDtClean[localPixel], fPixelClean[localPixel]->GetParameter(1), fPixelClean[localPixel]->GetParError(1), std::fabs(fPixelClean[localPixel]->GetParameter(2)), fPixelClean[localPixel]->GetParError(2), fPixelClean[localPixel]->GetChisquare()/fPixelClean[localPixel]->GetNDF());
    }
  }

  #endif
  TCanvas *cBarDtVsTot = new TCanvas(uniqueName("cCDetBarPixelTimingDtVsTot"), TString::Format("CDet bar %d pixel delta-t versus TOT", bar),
                                     1130, 35, 540, 420);
  cBarDtVsTot->Divide(4, 4, 0.001, 0.001);
  for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
    cBarDtVsTot->cd(localPixel + 1);
    hPixelDtVsTot[localPixel]->SetStats(kFALSE);
    hPixelDtVsTot[localPixel]->Draw("COLZ");
  }

  TCanvas *cBarAmalgamated = new TCanvas(uniqueName("cCDetBarTimingAmalgamated"), TString::Format("CDet bar %d amalgamated timing", bar),
                                         10, 480, 820, 500);
  cBarAmalgamated->Divide(2, 2);
  cBarAmalgamated->cd(1);
  hBarDt->SetStats(kTRUE);
  hBarDt->SetLineColor(kBlack);
  hBarDt->SetLineWidth(2);
  hBarDt->Draw("HIST");
  if (validBarLocalFit) {
    fBarLocal->SetLineColor(kRed + 1);
    fBarLocal->SetLineWidth(3);
    fBarLocal->Draw("SAME");
    if (fBarLocal->GetNDF() > 0) {
      AddFitResultsToStatsBox(hBarDt, fBarLocal->GetParameter(1), fBarLocal->GetParError(1),
                              std::fabs(fBarLocal->GetParameter(2)), fBarLocal->GetParError(2),
                              fBarLocal->GetChisquare()/fBarLocal->GetNDF());
    }
  }
  cBarAmalgamated->cd(2);
  hBarProjectedDt->SetStats(kTRUE);
  hBarProjectedDt->SetLineColor(kBlack);
  hBarProjectedDt->SetLineWidth(2);
  hBarProjectedDt->Draw("HIST");
  if (validBarProjectedLocalFit) {
    fBarProjectedLocal->SetLineColor(kRed + 1);
    fBarProjectedLocal->SetLineWidth(3);
    fBarProjectedLocal->Draw("SAME");
    if (fBarProjectedLocal->GetNDF() > 0) {
      AddFitResultsToStatsBox(hBarProjectedDt, fBarProjectedLocal->GetParameter(1), fBarProjectedLocal->GetParError(1),
                              std::fabs(fBarProjectedLocal->GetParameter(2)), fBarProjectedLocal->GetParError(2),
                              fBarProjectedLocal->GetChisquare()/fBarProjectedLocal->GetNDF());
    }
  }
  if (!projectionAvailable) {
    TPaveText *projectionNote = new TPaveText(0.16, 0.42, 0.84, 0.58, "NDC");
    projectionNote->SetFillColor(0);
    projectionNote->SetBorderSize(1);
    projectionNote->AddText(makeProjectionComparison ? "ECal projection unavailable: coordinate vectors are missing" : "ECal projection comparison disabled");
    projectionNote->Draw();
  }
  cBarAmalgamated->cd(3);
  hBarProjectedQualityDt->SetStats(kTRUE);
  hBarProjectedQualityDt->SetLineColor(kBlack);
  hBarProjectedQualityDt->SetLineWidth(2);
  hBarProjectedQualityDt->Draw("HIST");
  if (validBarProjectedQualityLocalFit) {
    fBarProjectedQualityLocal->SetLineColor(kRed + 1);
    fBarProjectedQualityLocal->SetLineWidth(3);
    fBarProjectedQualityLocal->Draw("SAME");
    if (fBarProjectedQualityLocal->GetNDF() > 0) {
      AddFitResultsToStatsBox(hBarProjectedQualityDt,
                              fBarProjectedQualityLocal->GetParameter(1), fBarProjectedQualityLocal->GetParError(1),
                              std::fabs(fBarProjectedQualityLocal->GetParameter(2)), fBarProjectedQualityLocal->GetParError(2),
                              fBarProjectedQualityLocal->GetChisquare()/fBarProjectedQualityLocal->GetNDF());
    }
  }
  cBarAmalgamated->cd(4);
  hBarProjectedQualityDtVsTot->SetStats(kTRUE);
  hBarProjectedQualityDtVsTot->Draw("COLZ");

  TCanvas *cBarSummary = new TCanvas(uniqueName("cCDetBarPixelTimingSummary"), TString::Format("CDet bar %d timing summary", bar),
                                     850, 480, 820, 500);
  cBarSummary->Divide(4, 2);
  cBarSummary->cd(1); gCentroid->Draw("AP"); gCentroid->GetXaxis()->SetLimits(pixelBase - 0.5, pixelBase + NumPaddles - 0.5);
  cBarSummary->cd(2); gCorrection->Draw("AP"); gCorrection->GetXaxis()->SetLimits(pixelBase - 0.5, pixelBase + NumPaddles - 0.5);
  cBarSummary->cd(3); gSigma->Draw("AP"); gSigma->GetXaxis()->SetLimits(pixelBase - 0.5, pixelBase + NumPaddles - 0.5);
  cBarSummary->cd(4); gEntries->Draw("AP"); gEntries->GetXaxis()->SetLimits(pixelBase - 0.5, pixelBase + NumPaddles - 0.5);
  cBarSummary->cd(5); gChi2Ndf->Draw("AP"); gChi2Ndf->GetXaxis()->SetLimits(pixelBase - 0.5, pixelBase + NumPaddles - 0.5);
  cBarSummary->cd(6); hValidity->SetMinimum(0.0); hValidity->SetMaximum(10.5); hValidity->Draw("HIST");
  cBarSummary->cd(7); hCentroidDistribution->Draw("HIST");
  cBarSummary->cd(8); hPredictedCentroidDistribution->Draw("HIST");
  if (referenceValid) {
    cBarSummary->cd(1);
    TPaveText *referenceNote = new TPaveText(0.14, 0.78, 0.58, 0.89, "NDC");
    referenceNote->SetFillColor(0);
    referenceNote->SetBorderSize(1);
    referenceNote->AddText(TString::Format("Bar-local #mu_{0} = %.4f #pm %.4f ns", referenceCentroid, referenceCentroidErr));
    referenceNote->Draw();
  }

  if (saveFitCanvases) {
    TString saveDir = TString::Format("%s/run_%d_stage_%d", fitCanvasDir.Data(), gRunNumber, gCalibrationStage);
    gSystem->mkdir(saveDir, kTRUE);
    cBarFits->SaveAs(TString::Format("%s/bar_%03d_pixels_%d-%d.pdf", saveDir.Data(), bar, pixelBase, pixelBase + NumPaddles - 1));
    cBarFits->SaveAs(TString::Format("%s/bar_%03d_pixels_%d-%d.png", saveDir.Data(), bar, pixelBase, pixelBase + NumPaddles - 1));
    cBarDtVsTot->SaveAs(TString::Format("%s/bar_%03d_dt_vs_tot.pdf", saveDir.Data(), bar));
    cBarAmalgamated->SaveAs(TString::Format("%s/bar_%03d_amalgamated.pdf", saveDir.Data(), bar));
    cBarAmalgamated->SaveAs(TString::Format("%s/bar_%03d_amalgamated.png", saveDir.Data(), bar));
    cBarSummary->SaveAs(TString::Format("%s/bar_%03d_summary.pdf", saveDir.Data(), bar));
  }

  if (saveCandidateTable) {
    const TString requestedBaseName = gSystem->BaseName(candidateOutput.Data());
    const TString activeBaseName = gSystem->BaseName(gCalibrationFile.c_str());
    if (requestedBaseName == "CDet_calibration.dat" || requestedBaseName == activeBaseName) {
      std::cerr << "[CDet pixel timing] ERROR: refusing to overwrite active calibration file with candidate output '" << candidateOutput << "'.\n";
    } else {
      std::ofstream output(candidateOutput.Data());
      if (!output.is_open()) {
        std::cerr << "[CDet pixel timing] ERROR: could not open candidate output '" << candidateOutput << "'.\n";
      } else {
        output << "# CANDIDATE DIAGNOSTIC PRODUCT -- NOT AN ACTIVE CALIBRATION\n"
               << "# run " << gRunNumber << "\n"
               << "# calibration_stage " << gCalibrationStage << "\n"
               << "# timing_units ns\n"
               << "# selection fourth-pass baseline good-hit vectors; instrumented physical logical IDs only\n"
               << "# ecal_energy_cut_gev " << ECalEnergyMin << " " << ECalEnergyMax << "\n"
               << "# dt = tECal - tCDetLE\n"
               << "# fit_interval_ns " << FitMin << " " << FitMax << "\n"
               << "# signal_fit_interval_ns " << FitMin << " " << FitMax << "\n"
               << "# calibration_centroid_source gaus(0)+pol1(3)_signal_mean\n"
               << "# peak_seed_interval_ns " << PeakSeedMin << " " << PeakSeedMax << "\n"
               << "# reference_scope bar-local\n"
               << "# mu0_ns " << referenceCentroid << " mu0_err_ns " << referenceCentroidErr << "\n"
               << "# correction c_i = mu_i - mu0; proposed convention tCDet_i' = tCDet_i + c_i\n"
               << "# pixel_id entries fit_status valid mu_ns mu_err_ns sigma_ns sigma_err_ns correction_ns correction_err_ns chi2 ndf chi2_ndf amplitude amplitude_err validity_code failure_reason\n";
        for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
          output << pixelBase + localPixel << " " << fitEntries[localPixel] << " " << fitStatus[localPixel] << " " << validFit[localPixel] << " " << centroid[localPixel] << " " << centroidErr[localPixel] << " " << sigma[localPixel] << " " << sigmaErr[localPixel] << " " << correction[localPixel] << " " << correctionErr[localPixel] << " " << chi2[localPixel] << " " << ndf[localPixel] << " " << chi2Ndf[localPixel] << " " << amplitude[localPixel] << " " << amplitudeErr[localPixel] << " " << validityCode[localPixel] << " \"" << failureReason[localPixel] << "\"\n";
        }
      }
    }
  }

  double correctionSum = 0.0;
  double correctionSquareSum = 0.0;
  double correctionMin = 0.0;
  double correctionMax = 0.0;
  int correctionCount = 0;
  if (referenceValid) {
    for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
      if (!validFit[localPixel]) continue;
      if (correctionCount == 0) correctionMin = correctionMax = correction[localPixel];
      correctionMin = std::min(correctionMin, correction[localPixel]);
      correctionMax = std::max(correctionMax, correction[localPixel]);
      correctionSum += correction[localPixel];
      correctionSquareSum += correction[localPixel]*correction[localPixel];
      ++correctionCount;
    }
  }
  const double correctionMean = correctionCount > 0 ? correctionSum/correctionCount : NAN;
  const double correctionRms = correctionCount > 0 ? std::sqrt(std::max(0.0, correctionSquareSum/correctionCount - correctionMean*correctionMean)) : NAN;

  std::cout << "[CDet pixel timing]\n"
            << "  requested logical pixel ID: " << requestedPixel << "\n"
            << "  normalized bar pixel base: " << pixelBase << "\n"
            << "  processed logical pixel IDs: " << pixelBase << "-" << pixelBase + NumPaddles - 1 << "\n"
            << "  ECal energy cut: [" << ECalEnergyMin << ", " << ECalEnergyMax << "] GeV\n"
            << "  peak search interval: [" << FitMin << ", " << FitMax << "] ns\n"
            << "  preferred peak-seed interval: [" << PeakSeedMin << ", " << PeakSeedMax << "] ns; fallbacks: " << peakSeedFallbackCount << "\n"
            << "  broad background fit interval: [" << HistMin << ", " << HistMax << "] ns\n"
            << "  calibration centroids use gaus(0)+pol1(3) over the full peak-search interval\n"
            << "  events passing ECal energy cut: " << energySelectedEventCount << " / " << nEvents << "\n"
            << "  amalgamated bar fit entries/status/valid: " << barFitEntries << " / " << barLocalFitStatus << " / " << validBarLocalFit << "\n"
            << "  amalgamated bar mean/sigma: " << (validBarLocalFit ? fBarLocal->GetParameter(1) : NAN) << " / " << (validBarLocalFit ? std::fabs(fBarLocal->GetParameter(2)) : NAN) << " ns\n"
            << "  ECal projection comparison available: " << projectionAvailable << "\n"
            << "  trajectory-matched bar hits: " << projectedBarHitCount << "\n"
            << "  trajectory-matched fit entries/status/valid: " << barProjectedFitEntries << " / " << barProjectedLocalFitStatus << " / " << validBarProjectedLocalFit << "\n"
            << "  trajectory-matched mean/sigma: " << (validBarProjectedLocalFit ? fBarProjectedLocal->GetParameter(1) : NAN) << " / " << (validBarProjectedLocalFit ? std::fabs(fBarProjectedLocal->GetParameter(2)) : NAN) << " ns\n"
            << "  default accepted TOT interval: (" << AcceptedTotMin << ", " << AcceptedTotMax << ") ns\n"
            << "  manual LE-versus-TOT polygons loaded for this bar: " << loadedBarPixelCuts << " from '" << pixelCutFile << "'\n"
            << "  quality hierarchy: a saved polygon replaces the common TOT interval for that pixel\n"
            << "  trajectory + TOT/polygon selected hits: " << projectedQualityBarHitCount << "\n"
            << "  trajectory + quality fit entries/status/valid: " << barProjectedQualityFitEntries << " / " << barProjectedQualityLocalFitStatus << " / " << validBarProjectedQualityLocalFit << "\n"
            << "  trajectory + quality mean/sigma: " << (validBarProjectedQualityLocalFit ? fBarProjectedQualityLocal->GetParameter(1) : NAN) << " / " << (validBarProjectedQualityLocalFit ? std::fabs(fBarProjectedQualityLocal->GetParameter(2)) : NAN) << " ns\n"
            << "  physical instrumented pixels considered: " << instrumentedPixels << "\n"
            << "  pixels with sufficient statistics: " << sufficientStatistics << "\n"
            << "  successful ROOT fits: " << successfulFits << "\n"
            << "  pixels included in bar-local reference: " << referencePixels << "\n"
            << "  bar-local mu0: " << referenceCentroid << " +/- " << referenceCentroidErr << " ns\n"
            << "  candidate correction mean/RMS: " << correctionMean << " / " << correctionRms << " ns\n"
            << "  candidate correction range: [" << correctionMin << ", " << correctionMax << "] ns\n"
            << "  failures -- low statistics: " << lowStatisticsCount << ", ROOT fit: " << rootFitFailureCount << ", invalid parameters: " << invalidParameterCount << ", boundary: " << boundaryCount << ", sigma: " << sigmaCount << ", amplitude: " << amplitudeCount << ", NDF: " << ndfCount << ", chi2/NDF: " << chi2Count << "\n"
            << "  proposed sign: tCDet_i' = tCDet_i + c_i; no corrections were applied\n";
}

// Survey the CDet-ECal timing peak in every physical bar using the same
// trajectory projection and per-pixel TOT/polygon hierarchy as the focused
// extractCDetBarPixelTimingOffsets diagnostic.  This is diagnostic only.
void surveyCDetBarTimingPeaks(double ECalEnergyMin = 3.0, double ECalEnergyMax = 4.5,
                              bool saveDiagnostics = false,
                              TString diagnosticOutput = "CDet_bar_timing_survey.root",
                              double Width = 1.0,
                              double HistMin = -60.0, double HistMax = 30.0,
                              double FitMin = -30.0, double FitMax = 0.0,
                              int minEntries = 100,
                              double minSigma = 0.5, double maxSigma = 20.0,
                              double maxChi2Ndf = 10.0,
                              double PeakSeedMin = -25.0, double PeakSeedMax = -5.0,
                              double AcceptedTotMin = gCDetDiagnosticAcceptedTotMin,
                              double AcceptedTotMax = gCDetDiagnosticAcceptedTotMax,
                              TString pixelCutFile = "CDet_pixel_quality_cuts.root") {
  TH1::AddDirectory(kFALSE);
  const int barsPerLayer = NumCDetPaddles/(NumLayers*NumPaddles);
  const int totalBars = NumCDetPaddles/NumPaddles;
  const size_t nEvents = vGoodLe.size();

  if (Width <= 0.0 || HistMax <= HistMin || FitMax <= FitMin ||
      ECalEnergyMax <= ECalEnergyMin || AcceptedTotMax <= AcceptedTotMin) {
    std::cerr << "[CDet bar timing survey] ERROR: invalid histogram, fit, energy, or TOT limits.\n";
    return;
  }
  if (nEvents == 0 || vGoodTot.size() != nEvents || vGoodID.size() != nEvents ||
      vCDetGoodX.size() != nEvents || vCDetGoodY.size() != nEvents ||
      vCDetGoodZ.size() != nEvents || v_GoodECalX.size() != nEvents ||
      v_GoodECalY.size() != nEvents || v_GoodECalE.size() != nEvents ||
      v_GoodECalAdcTime.size() != nEvents) {
    std::cerr << "[CDet bar timing survey] ERROR: required good-hit/ECal vectors are empty or misaligned.\n";
    return;
  }

  static unsigned long surveyInvocation = 0;
  const TString suffix = TString::Format("_%lu", ++surveyInvocation);
  const int nBins = std::max(1, int(std::ceil((HistMax-HistMin)/Width)));

  std::vector<double> xMin(totalBars, std::numeric_limits<double>::infinity());
  std::vector<double> xMax(totalBars, -std::numeric_limits<double>::infinity());
  std::vector<double> ySum(totalBars, 0.0);
  std::vector<long long> geometryHits(totalBars, 0);
  for (size_t ev = 0; ev < nEvents; ++ev) {
    const size_t nh = std::min(std::min(vGoodID[ev].size(), vCDetGoodX[ev].size()),
                               std::min(vCDetGoodY[ev].size(), vCDetGoodZ[ev].size()));
    for (size_t ih = 0; ih < nh; ++ih) {
      const int pixel = vGoodID[ev][ih];
      if (pixel < 0 || pixel >= NumCDetPaddles || IsUnusedPixel(pixel)) continue;
      const int bar = pixel/NumPaddles;
      xMin[bar] = std::min(xMin[bar], vCDetGoodX[ev][ih]);
      xMax[bar] = std::max(xMax[bar], vCDetGoodX[ev][ih]);
      ySum[bar] += vCDetGoodY[ev][ih];
      ++geometryHits[bar];
    }
  }

  std::vector<TCutG*> pixelCuts(NumCDetPaddles, nullptr);
  int loadedCuts = 0;
  TFile cutInput(pixelCutFile, "READ");
  if (!cutInput.IsZombie()) {
    for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
      if (IsUnusedPixel(pixel)) continue;
      pixelCuts[pixel] = LoadCDetPixelLeTotCut(cutInput, pixel);
      if (pixelCuts[pixel]) ++loadedCuts;
    }
    cutInput.Close();
  } else {
    std::cout << "[CDet bar timing survey] No readable pixel-cut file '" << pixelCutFile
              << "'; using the common TOT interval for every pixel.\n";
  }

  std::vector<TH1D*> spectra(totalBars, nullptr);
  for (int bar = 0; bar < totalBars; ++bar) {
    const int layer = bar/barsPerLayer + 1;
    const int layerBar = bar%barsPerLayer;
    spectra[bar] = new TH1D(TString::Format("hCDetBarTimingSurveyL%dB%02d%s", layer, layerBar, suffix.Data()),
                            TString::Format("Layer %d bar %d; t_{ECal}-t_{CDet,LE} (ns); Selected hits", layer, layerBar),
                            nBins, HistMin, HistMax);
  }

  const double halfPixelPitch = 0.5*0.00525*0.5*(XCorr1+XCorr2);
  long long energySelectedEvents = 0;
  long long selectedHits = 0;
  for (size_t ev = 0; ev < nEvents; ++ev) {
    if (!std::isfinite(v_GoodECalE[ev]) || v_GoodECalE[ev] < ECalEnergyMin || v_GoodECalE[ev] > ECalEnergyMax) continue;
    ++energySelectedEvents;
    const size_t nh = std::min(std::min(vGoodLe[ev].size(), vGoodTot[ev].size()),
                      std::min(std::min(vGoodID[ev].size(), vCDetGoodX[ev].size()),
                               std::min(vCDetGoodY[ev].size(), vCDetGoodZ[ev].size())));
    for (size_t ih = 0; ih < nh; ++ih) {
      const int pixel = vGoodID[ev][ih];
      if (pixel < 0 || pixel >= NumCDetPaddles || IsUnusedPixel(pixel)) continue;
      const int bar = pixel/NumPaddles;
      if (geometryHits[bar] == 0 || !std::isfinite(vCDetGoodZ[ev][ih])) continue;
      const double projectedX = v_GoodECalX[ev]*vCDetGoodZ[ev][ih]/ECal_dist;
      const double projectedY = v_GoodECalY[ev]*vCDetGoodZ[ev][ih]/ECal_dist;
      const double barY = ySum[bar]/geometryHits[bar];
      if (projectedX < xMin[bar]-halfPixelPitch || projectedX > xMax[bar]+halfPixelPitch ||
          std::fabs(projectedY-barY) > CDet_y_half_length) continue;
      const double le = vGoodLe[ev][ih];
      const double tot = vGoodTot[ev][ih];
      const bool passesQuality = pixelCuts[pixel] ? pixelCuts[pixel]->IsInside(tot, le)
                                                  : (tot > AcceptedTotMin && tot < AcceptedTotMax);
      if (!passesQuality) continue;
      spectra[bar]->Fill(v_GoodECalAdcTime[ev]-le);
      ++selectedHits;
    }
  }

  std::vector<TF1*> fits(totalBars, nullptr);
  std::vector<bool> valid(totalBars, false);
  std::vector<double> mean(totalBars, NAN), meanErr(totalBars, NAN);
  std::vector<double> sigma(totalBars, NAN), sigmaErr(totalBars, NAN);
  std::vector<double> significance(totalBars, NAN), yield(totalBars, NAN), yieldErr(totalBars, NAN);
  // Preserve the deliberately broad fit acceptance above, but use a stricter,
  // separately reported quality classification for detector-wide summaries.
  enum BarSummaryStatus {
    kInsufficientStatistics = 0,
    kRejectedFit = 1,
    kLowSignificance = 2,
    kExcessiveWidth = 3,
    kUncertainCentroid = 4,
    kUnstableYield = 5,
    kRecommendedFit = 6
  };
  const double recommendedMinSignificance = 3.0;
  const double recommendedMaxSigma = 10.0;
  const double recommendedMaxCentroidError = 3.0;
  const double recommendedMaxRelativeYieldError = 1.0;
  int sufficientBars = 0, acceptedBars = 0;
  for (int bar = 0; bar < totalBars; ++bar) {
    TH1D *h = spectra[bar];
    const double fitEntries = h->Integral(h->FindBin(FitMin+1e-9), h->FindBin(FitMax-1e-9));
    if (fitEntries < minEntries) continue;
    ++sufficientBars;
    int seedLo = h->FindBin(std::max(FitMin, PeakSeedMin)+1e-9);
    int seedHi = h->FindBin(std::min(FitMax, PeakSeedMax)-1e-9);
    if (seedHi < seedLo) { seedLo = h->FindBin(FitMin+1e-9); seedHi = h->FindBin(FitMax-1e-9); }
    int seedBin = seedLo;
    for (int bin = seedLo+1; bin <= seedHi; ++bin)
      if (h->GetBinContent(bin) > h->GetBinContent(seedBin)) seedBin = bin;
    const double seedMean = h->GetBinCenter(seedBin);
    const double seedAmplitude = std::max(1.0, h->GetBinContent(seedBin));
    fits[bar] = new TF1(TString::Format("fCDetBarTimingSurvey%03d%s", bar, suffix.Data()), "gaus(0)+pol1(3)", FitMin, FitMax);
    fits[bar]->SetParameters(seedAmplitude, seedMean, 3.0, 0.0, 0.0);
    fits[bar]->SetParLimits(0, 0.0, std::max(10.0, 10.0*seedAmplitude));
    fits[bar]->SetParLimits(1, FitMin, FitMax);
    fits[bar]->SetParLimits(2, minSigma, maxSigma);
    const int status = h->Fit(fits[bar], "RQN0");
    const double amp = fits[bar]->GetParameter(0);
    const double ampErr = fits[bar]->GetParError(0);
    mean[bar] = fits[bar]->GetParameter(1);
    meanErr[bar] = fits[bar]->GetParError(1);
    sigma[bar] = std::fabs(fits[bar]->GetParameter(2));
    sigmaErr[bar] = fits[bar]->GetParError(2);
    const int ndf = fits[bar]->GetNDF();
    const double chi2Ndf = ndf > 0 ? fits[bar]->GetChisquare()/ndf : NAN;
    valid[bar] = status == 0 && std::isfinite(amp) && amp > 0.0 && std::isfinite(ampErr) && ampErr > 0.0 &&
                 std::isfinite(mean[bar]) && mean[bar] > FitMin && mean[bar] < FitMax &&
                 std::isfinite(sigma[bar]) && sigma[bar] >= minSigma && sigma[bar] <= maxSigma &&
                 ndf > 0 && std::isfinite(chi2Ndf) && chi2Ndf <= maxChi2Ndf;
    if (!valid[bar]) continue;
    significance[bar] = amp/ampErr;
    yield[bar] = amp*std::sqrt(2.0*TMath::Pi())*sigma[bar]/Width;
    if (sigmaErr[bar] >= 0.0)
      yieldErr[bar] = yield[bar]*std::sqrt((ampErr/amp)*(ampErr/amp) + (sigmaErr[bar]/sigma[bar])*(sigmaErr[bar]/sigma[bar]));
    ++acceptedBars;
  }

  std::vector<int> summaryStatus(totalBars, kInsufficientStatistics);
  int statusCounts[7] = {0, 0, 0, 0, 0, 0, 0};
  for (int bar = 0; bar < totalBars; ++bar) {
    TH1D *h = spectra[bar];
    const double fitEntries = h->Integral(h->FindBin(FitMin+1e-9), h->FindBin(FitMax-1e-9));
    int status = kInsufficientStatistics;
    if (fitEntries >= minEntries) {
      if (!valid[bar]) {
        status = kRejectedFit;
      } else if (!std::isfinite(significance[bar]) || significance[bar] <= recommendedMinSignificance) {
        status = kLowSignificance;
      } else if (!std::isfinite(sigma[bar]) || sigma[bar] >= recommendedMaxSigma) {
        status = kExcessiveWidth;
      } else if (!std::isfinite(meanErr[bar]) || meanErr[bar] <= 0.0 ||
                 meanErr[bar] >= recommendedMaxCentroidError) {
        status = kUncertainCentroid;
      } else if (!std::isfinite(yield[bar]) || yield[bar] <= 0.0 ||
                 !std::isfinite(yieldErr[bar]) || yieldErr[bar] < 0.0 ||
                 yieldErr[bar]/yield[bar] >= recommendedMaxRelativeYieldError) {
        status = kUnstableYield;
      } else {
        status = kRecommendedFit;
      }
    }
    summaryStatus[bar] = status;
    ++statusCounts[status];
  }

  std::vector<TGraphErrors*> gMean(NumLayers), gSigma(NumLayers), gYield(NumLayers);
  std::vector<TGraph*> gSignificance(NumLayers), gStatus(NumLayers);
  std::vector<TCanvas*> canvases(NumLayers, nullptr);
  for (int layerIndex = 0; layerIndex < NumLayers; ++layerIndex) {
    gMean[layerIndex] = new TGraphErrors();
    gSigma[layerIndex] = new TGraphErrors();
    gYield[layerIndex] = new TGraphErrors();
    gSignificance[layerIndex] = new TGraph();
    gStatus[layerIndex] = new TGraph();
    int point = 0;
    for (int layerBar = 0; layerBar < barsPerLayer; ++layerBar) {
      const int bar = layerIndex*barsPerLayer + layerBar;
      gStatus[layerIndex]->SetPoint(layerBar, layerBar, summaryStatus[bar]);
      if (summaryStatus[bar] != kRecommendedFit) continue;
      gMean[layerIndex]->SetPoint(point, layerBar, mean[bar]);
      gMean[layerIndex]->SetPointError(point, 0.0, meanErr[bar]);
      gSigma[layerIndex]->SetPoint(point, layerBar, sigma[bar]);
      gSigma[layerIndex]->SetPointError(point, 0.0, sigmaErr[bar]);
      gSignificance[layerIndex]->SetPoint(point, layerBar, significance[bar]);
      gYield[layerIndex]->SetPoint(point, layerBar, yield[bar]);
      gYield[layerIndex]->SetPointError(point, 0.0, yieldErr[bar]);
      ++point;
    }
    const int layer = layerIndex+1;
    canvases[layerIndex] = new TCanvas(TString::Format("cCDetBarTimingSurveyLayer%d%s", layer, suffix.Data()),
                                       TString::Format("CDet layer %d bar timing survey", layer), 1250, 1050);
    canvases[layerIndex]->Divide(2,3);
    canvases[layerIndex]->cd(1);
    gMean[layerIndex]->SetTitle(TString::Format("Layer %d recommended timing centroids, %.1f<E_{ECal}<%.1f GeV;Bar in layer;Centroid (ns)", layer, ECalEnergyMin, ECalEnergyMax));
    gMean[layerIndex]->SetMarkerStyle(20); gMean[layerIndex]->Draw("AP");
    canvases[layerIndex]->cd(2);
    gSigma[layerIndex]->SetTitle(TString::Format("Layer %d recommended Gaussian widths;Bar in layer;#sigma (ns)", layer));
    gSigma[layerIndex]->SetMarkerStyle(20); gSigma[layerIndex]->Draw("AP");
    canvases[layerIndex]->cd(3);
    gSignificance[layerIndex]->SetTitle(TString::Format("Layer %d recommended peak significances;Bar in layer;Amplitude / amplitude error", layer));
    gSignificance[layerIndex]->SetMarkerStyle(20); gSignificance[layerIndex]->Draw("AP");
    canvases[layerIndex]->cd(4);
    gYield[layerIndex]->SetTitle(TString::Format("Layer %d recommended Gaussian yields;Bar in layer;Peak yield (hits)", layer));
    gYield[layerIndex]->SetMarkerStyle(20); gYield[layerIndex]->Draw("AP");
    canvases[layerIndex]->cd(5);
    gStatus[layerIndex]->SetTitle(TString::Format("Layer %d bar fit-status classification;Bar in layer;Status code", layer));
    gStatus[layerIndex]->SetMarkerStyle(20);
    gStatus[layerIndex]->SetMarkerSize(0.8);
    gStatus[layerIndex]->Draw("AP");
    gStatus[layerIndex]->GetYaxis()->SetRangeUser(-0.5, 6.5);
    gStatus[layerIndex]->GetYaxis()->SetNdivisions(7, false);
    canvases[layerIndex]->cd(6);
    TPaveText *qualityKey = new TPaveText(0.06, 0.08, 0.94, 0.92, "NDC");
    qualityKey->SetFillStyle(0);
    qualityKey->SetBorderSize(1);
    qualityKey->SetTextAlign(12);
    qualityKey->AddText("Recommended summary requirements:");
    qualityKey->AddText(TString::Format("significance > %.1f", recommendedMinSignificance));
    qualityKey->AddText(TString::Format("#sigma < %.1f ns", recommendedMaxSigma));
    qualityKey->AddText(TString::Format("centroid uncertainty < %.1f ns", recommendedMaxCentroidError));
    qualityKey->AddText(TString::Format("relative yield uncertainty < %.0f%%", 100.0*recommendedMaxRelativeYieldError));
    qualityKey->AddText("Status: 0 low N, 1 rejected fit, 2 low significance,");
    qualityKey->AddText("3 excessive width, 4 uncertain centroid,");
    qualityKey->AddText("5 unstable yield, 6 recommended");
    qualityKey->Draw();
    canvases[layerIndex]->Update();
    canvas_vector.push_back(canvases[layerIndex]);
  }

  if (saveDiagnostics) {
    TDirectory *previousDirectory = gDirectory;
    TFile output(diagnosticOutput, "RECREATE");
    if (output.IsZombie()) {
      std::cerr << "[CDet bar timing survey] ERROR: could not create '" << diagnosticOutput << "'.\n";
    } else {
      for (int layerIndex = 0; layerIndex < NumLayers; ++layerIndex) {
        TDirectory *layerDirectory = output.mkdir(TString::Format("Layer%d", layerIndex+1));
        for (int layerBar = 0; layerBar < barsPerLayer; ++layerBar) {
          const int bar = layerIndex*barsPerLayer + layerBar;
          TDirectory *barDirectory = layerDirectory->mkdir(TString::Format("bar%02d", layerBar));
          barDirectory->cd();
          spectra[bar]->Write();
          if (fits[bar]) fits[bar]->Write();
        }
      }
      output.cd();
      for (int layerIndex = 0; layerIndex < NumLayers; ++layerIndex) canvases[layerIndex]->Write();
      output.Close();
    }
    if (previousDirectory) previousDirectory->cd();
  }

  for (TCutG *cut : pixelCuts) delete cut;
  int geometryBars = 0;
  for (long long count : geometryHits) if (count > 0) ++geometryBars;
  std::cout << "[CDet bar timing survey]\n"
            << "  ECal energy cut: [" << ECalEnergyMin << ", " << ECalEnergyMax << "] GeV\n"
            << "  trajectory selection: projected ECal (x,y) lies in the same CDet half-bar\n"
            << "  default TOT interval: (" << AcceptedTotMin << ", " << AcceptedTotMax << ") ns\n"
            << "  manual pixel polygons loaded: " << loadedCuts << " from '" << pixelCutFile << "'\n"
            << "  events passing energy cut: " << energySelectedEvents << " / " << nEvents << "\n"
            << "  selected hits: " << selectedHits << "\n"
            << "  bars with geometry / sufficient statistics / accepted fits: "
            << geometryBars << " / " << sufficientBars << " / " << acceptedBars << "\n"
            << "  summary classification (low N / rejected / low significance / excessive width /"
               " uncertain centroid / unstable yield / recommended): "
            << statusCounts[kInsufficientStatistics] << " / " << statusCounts[kRejectedFit] << " / "
            << statusCounts[kLowSignificance] << " / " << statusCounts[kExcessiveWidth] << " / "
            << statusCounts[kUncertainCentroid] << " / " << statusCounts[kUnstableYield] << " / "
            << statusCounts[kRecommendedFit] << "\n"
            << "  canvases: " << canvases[0]->GetName() << ", " << canvases[1]->GetName() << "\n";
  if (saveDiagnostics) std::cout << "  ROOT diagnostics: " << diagnosticOutput << "\n";
}

void extractAllCDetPixelTimingOffsets(bool generateOffsets = false, double Width = 1.0,
                                      double HistMin = -60.0, double HistMax = 30.0,
                                      double FitMin = -30.0, double FitMax = 0.0,
                                      int minEntries = 100, double minSigma = 0.5, double maxSigma = 20.0,
                                      double maxChi2Ndf = 10.0, double centroidEdgeMargin = 1.0,
                                      TString calibrationOutput = "", TString fitResultsOutput = "",
                                      double ECalEnergyMin = 1.0, double ECalEnergyMax = 12.0,
                                      double NReject = 2.5,
                                      double PeakSeedMin = -25.0, double PeakSeedMax = -5.0) {
  TH1::AddDirectory(kFALSE);
  gLastCalibrationFitSucceeded = false;

  if (!generateOffsets) {
    std::cout << "[CDet all-pixel timing] Offset generation disabled. Call extractAllCDetPixelTimingOffsets(true) to fit pixels and write updated calibration constants.\n";
    return;
  }

  if (calibrationOutput.IsNull()) calibrationOutput = gCalibrationFile.c_str();
  if (fitResultsOutput.IsNull()) fitResultsOutput = gRunNumber == 0 ? TString::Format("CDet_pixel_timing_fit_results_group%d.dat", ecalClusterGroupIndex) : "CDet_pixel_timing_fit_results.dat";
  const bool hadLoadedPixelOffsets = gPixelToffsetLoaded;

  if (Width <= 0.0 || HistMin >= HistMax || FitMin >= FitMax || FitMin < HistMin || FitMax > HistMax || PeakSeedMin >= PeakSeedMax || PeakSeedMin < FitMin || PeakSeedMax > FitMax || minEntries < 1 || minSigma <= 0.0 || minSigma >= maxSigma || maxChi2Ndf <= 0.0 || centroidEdgeMargin < 0.0 || ECalEnergyMin >= ECalEnergyMax || NReject <= 0.0) {
    std::cerr << "[CDet all-pixel timing] ERROR: invalid histogram, fit, or quality-limit argument.\n";
    return;
  }

  const TString calibrationBaseName = gSystem->BaseName(calibrationOutput.Data());
  const TString fitResultsBaseName = gSystem->BaseName(fitResultsOutput.Data());
  if (calibrationBaseName == "CDet_calibration.dat") {
    std::cerr << "[CDet all-pixel timing] ERROR: refusing to overwrite the legacy-method calibration file '" << calibrationOutput << "'.\n";
    return;
  }
  if (calibrationOutput == fitResultsOutput || calibrationBaseName == fitResultsBaseName) {
    std::cerr << "[CDet all-pixel timing] ERROR: calibration and fit-results outputs must be different files.\n";
    return;
  }

  const int nBins = (int)((HistMax - HistMin)/Width);
  if (nBins < 1) {
    std::cerr << "[CDet all-pixel timing] ERROR: histogram binning produces fewer than one bin.\n";
    return;
  }

  const size_t nEvents = vGoodLe.size();
  if (vGoodID.size() != nEvents || v_GoodECalAdcTime.size() != nEvents || v_GoodECalE.size() != nEvents) {
    std::cerr << "[CDet all-pixel timing] ERROR: event vectors are not aligned: LE=" << vGoodLe.size() << ", ID=" << vGoodID.size() << ", ECal time=" << v_GoodECalAdcTime.size() << ", ECal energy=" << v_GoodECalE.size() << ".\n";
    return;
  }
  for (size_t ev = 0; ev < nEvents; ++ev) {
    if (vGoodID[ev].size() != vGoodLe[ev].size()) {
      std::cerr << "[CDet all-pixel timing] ERROR: LE and ID vectors differ in event " << ev << ".\n";
      return;
    }
  }

  static unsigned long invocation = 0;
  const unsigned long tag = ++invocation;
  auto uniqueName = [tag](const char *base) { return TString::Format("%s_%lu", base, tag); };

  std::vector<TH1D*> hPixelDt(NumCDetPaddles, nullptr);
  std::vector<TH1D*> hPixelDtClean(NumCDetPaddles, nullptr);
  std::vector<TF1*> fPixelLocal(NumCDetPaddles, nullptr);
  std::vector<TF1*> fPixelBackgroundReject(NumCDetPaddles, nullptr);
  std::vector<TF1*> fPixelBackground(NumCDetPaddles, nullptr);
  std::vector<TF1*> fPixelClean(NumCDetPaddles, nullptr);
  std::vector<int> fitEntries(NumCDetPaddles, 0), fitStatus(NumCDetPaddles, -1), validityCode(NumCDetPaddles, 0), ndf(NumCDetPaddles, 0);
  std::vector<std::string> failureReason(NumCDetPaddles, "not processed");
  std::vector<double> amplitude(NumCDetPaddles, NAN), amplitudeErr(NumCDetPaddles, NAN), centroid(NumCDetPaddles, NAN), centroidErr(NumCDetPaddles, NAN), sigma(NumCDetPaddles, NAN), sigmaErr(NumCDetPaddles, NAN);
  std::vector<double> backgroundAmplitude(NumCDetPaddles, NAN), backgroundAmplitudeLimit(NumCDetPaddles, NAN), backgroundMean(NumCDetPaddles, NAN), backgroundSigma(NumCDetPaddles, NAN), rejectLow(NumCDetPaddles, NAN), rejectHigh(NumCDetPaddles, NAN);
  std::vector<double> chi2(NumCDetPaddles, NAN), chi2Ndf(NumCDetPaddles, NAN), correction(NumCDetPaddles, 0.0), correctionErr(NumCDetPaddles, NAN), totalOffset(NumCDetPaddles, 0.0);
  std::vector<bool> validFit(NumCDetPaddles, false);

  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) hPixelDt[pixelID] = new TH1D(uniqueName(TString::Format("hCDetAllPixelTimingDt_%d", pixelID)), TString::Format("After ECal energy cut, logical pixel ID %d;t_{ECal}-t_{CDet,LE} (ns);Baseline hits", pixelID), nBins, HistMin, HistMax);

  size_t energySelectedEventCount = 0;
  for (size_t ev = 0; ev < nEvents; ++ev) {
    if (v_GoodECalE[ev] < ECalEnergyMin || v_GoodECalE[ev] > ECalEnergyMax) continue;
    ++energySelectedEventCount;
    const double tECal = v_GoodECalAdcTime[ev];
    for (size_t ihit = 0; ihit < vGoodLe[ev].size(); ++ihit) {
      const int pixelID = vGoodID[ev][ihit];
      if (pixelID < 0 || pixelID >= NumCDetPaddles) continue;
      hPixelDt[pixelID]->Fill(tECal - vGoodLe[ev][ihit]);
    }
  }

  int validPixelCount = 0;
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    if (IsUnusedPixel(pixelID)) {
      validityCode[pixelID] = 1;
      failureReason[pixelID] = "unused pixel";
      continue;
    }
    const int fitBinMin = hPixelDt[pixelID]->FindBin(FitMin);
    const int fitBinMax = hPixelDt[pixelID]->FindBin(FitMax);
    fitEntries[pixelID] = (int)hPixelDt[pixelID]->Integral(fitBinMin, fitBinMax);
    if (fitEntries[pixelID] < minEntries) {
      validityCode[pixelID] = 2;
      failureReason[pixelID] = "insufficient entries";
      continue;
    }

    int peakBinMin = hPixelDt[pixelID]->FindBin(PeakSeedMin);
    int peakBinMax = hPixelDt[pixelID]->FindBin(PeakSeedMax);
    int peakBin = peakBinMin;
    for (int bin = peakBinMin + 1; bin <= peakBinMax; ++bin) if (hPixelDt[pixelID]->GetBinContent(bin) > hPixelDt[pixelID]->GetBinContent(peakBin)) peakBin = bin;
    if (hPixelDt[pixelID]->GetBinContent(peakBin) <= 0.0) {
      peakBin = fitBinMin;
      for (int bin = fitBinMin + 1; bin <= fitBinMax; ++bin) if (hPixelDt[pixelID]->GetBinContent(bin) > hPixelDt[pixelID]->GetBinContent(peakBin)) peakBin = bin;
    }
    const double peak = hPixelDt[pixelID]->GetBinCenter(peakBin);
    const double localBackground = 0.5*(hPixelDt[pixelID]->GetBinContent(fitBinMin) + hPixelDt[pixelID]->GetBinContent(fitBinMax));
    fPixelLocal[pixelID] = new TF1(uniqueName(TString::Format("fCDetAllPixelTimingLocal_%d", pixelID)), "gaus(0)+pol1(3)", FitMin, FitMax);
    fPixelLocal[pixelID]->SetParameters(std::max(1.0, hPixelDt[pixelID]->GetBinContent(peakBin) - localBackground), peak, std::max(Width, (FitMax - FitMin)/10.0), localBackground, 0.0);
    fitStatus[pixelID] = hPixelDt[pixelID]->Fit(fPixelLocal[pixelID], "RQN0");
    const double localMean = fPixelLocal[pixelID]->GetParameter(1);
    const double localSigma = std::fabs(fPixelLocal[pixelID]->GetParameter(2));
    if (fitStatus[pixelID] != 0 || !std::isfinite(fPixelLocal[pixelID]->GetParameter(0)) || fPixelLocal[pixelID]->GetParameter(0) <= 0.0 || !std::isfinite(localMean) || !std::isfinite(localSigma) || localSigma < minSigma || localSigma > maxSigma || localMean <= FitMin || localMean >= FitMax) {
      validityCode[pixelID] = 3;
      failureReason[pixelID] = "invalid local signal fit";
      continue;
    }

    amplitude[pixelID] = fPixelLocal[pixelID]->GetParameter(0);
    amplitudeErr[pixelID] = fPixelLocal[pixelID]->GetParError(0);
    centroid[pixelID] = localMean;
    centroidErr[pixelID] = fPixelLocal[pixelID]->GetParError(1);
    sigma[pixelID] = localSigma;
    sigmaErr[pixelID] = fPixelLocal[pixelID]->GetParError(2);
    chi2[pixelID] = fPixelLocal[pixelID]->GetChisquare();
    ndf[pixelID] = fPixelLocal[pixelID]->GetNDF();
    chi2Ndf[pixelID] = ndf[pixelID] > 0 ? chi2[pixelID]/ndf[pixelID] : NAN;

    if (!std::isfinite(amplitude[pixelID]) || !std::isfinite(amplitudeErr[pixelID]) || !std::isfinite(centroid[pixelID]) || !std::isfinite(centroidErr[pixelID]) || !std::isfinite(sigma[pixelID]) || !std::isfinite(sigmaErr[pixelID]) || !std::isfinite(chi2[pixelID]) || centroidErr[pixelID] <= 0.0 || sigmaErr[pixelID] <= 0.0) { validityCode[pixelID] = 4; failureReason[pixelID] = "invalid parameter or uncertainty"; continue; }
    if (centroid[pixelID] - FitMin <= centroidEdgeMargin || FitMax - centroid[pixelID] <= centroidEdgeMargin) { validityCode[pixelID] = 5; failureReason[pixelID] = "centroid near fit boundary"; continue; }
    if (sigma[pixelID] < minSigma || sigma[pixelID] > maxSigma) { validityCode[pixelID] = 6; failureReason[pixelID] = "sigma outside limits"; continue; }
    if (amplitude[pixelID] <= 0.0) { validityCode[pixelID] = 7; failureReason[pixelID] = "nonpositive Gaussian amplitude"; continue; }
    if (ndf[pixelID] <= 0) { validityCode[pixelID] = 8; failureReason[pixelID] = "invalid NDF"; continue; }
    if (!std::isfinite(chi2Ndf[pixelID]) || chi2Ndf[pixelID] > maxChi2Ndf) { validityCode[pixelID] = 9; failureReason[pixelID] = "chi2/NDF outside limit"; continue; }
    validityCode[pixelID] = 10;
    failureReason[pixelID] = "valid";
    validFit[pixelID] = true;
    ++validPixelCount;

    // Background-subtracted diagnostics removed; production uses fPixelLocal.
#if 0
    // Production centroids and the detector-wide reference come from fPixelLocal.
    rejectLow[pixelID] = std::max(HistMin, localMean - NReject*localSigma);
    rejectHigh[pixelID] = std::min(HistMax, localMean + NReject*localSigma);
    CDetTimingBackgroundRejectLow = rejectLow[pixelID];
    CDetTimingBackgroundRejectHigh = rejectHigh[pixelID];
    double maxSidebandBinContent = 0.0;
    int leftSidebandBins = 0;
    int rightSidebandBins = 0;
    for (int bin = 1; bin <= hPixelDt[pixelID]->GetNbinsX(); ++bin) {
      const double binCenter = hPixelDt[pixelID]->GetBinCenter(bin);
      if (rejectLow[pixelID] <= binCenter && binCenter <= rejectHigh[pixelID]) continue;
      maxSidebandBinContent = std::max(maxSidebandBinContent, hPixelDt[pixelID]->GetBinContent(bin));
      if (binCenter < rejectLow[pixelID]) ++leftSidebandBins;
      if (binCenter > rejectHigh[pixelID]) ++rightSidebandBins;
    }
    backgroundAmplitudeLimit[pixelID] = maxSidebandBinContent;
    if (maxSidebandBinContent <= 0.0 || leftSidebandBins < 3 || rightSidebandBins < 3) {
      continue;
    }

    fPixelBackgroundReject[pixelID] = new TF1(uniqueName(TString::Format("fCDetAllPixelTimingBackgroundReject_%d", pixelID)), CDetTimingBackgroundGaussianReject, HistMin, HistMax, 3);
    fPixelBackgroundReject[pixelID]->SetParameters(0.8*maxSidebandBinContent, 0.5*(HistMin + HistMax), std::max(Width, (HistMax - HistMin)/3.0));
    fPixelBackgroundReject[pixelID]->SetParLimits(0, 0.0, maxSidebandBinContent);
    fPixelBackgroundReject[pixelID]->SetParLimits(1, HistMin, HistMax);
    fPixelBackgroundReject[pixelID]->SetParLimits(2, Width, HistMax - HistMin);
    const int backgroundStatus = hPixelDt[pixelID]->Fit(fPixelBackgroundReject[pixelID], "RQN0");
    backgroundAmplitude[pixelID] = fPixelBackgroundReject[pixelID]->GetParameter(0);
    backgroundMean[pixelID] = fPixelBackgroundReject[pixelID]->GetParameter(1);
    backgroundSigma[pixelID] = std::fabs(fPixelBackgroundReject[pixelID]->GetParameter(2));
    if (backgroundStatus != 0 || !std::isfinite(backgroundAmplitude[pixelID]) || !std::isfinite(backgroundMean[pixelID]) || !std::isfinite(backgroundSigma[pixelID]) || backgroundAmplitude[pixelID] < 0.0 || backgroundSigma[pixelID] <= 0.0) {
      continue;
    }

    fPixelBackground[pixelID] = new TF1(uniqueName(TString::Format("fCDetAllPixelTimingBackground_%d", pixelID)), "gaus", HistMin, HistMax);
    fPixelBackground[pixelID]->SetParameters(backgroundAmplitude[pixelID], backgroundMean[pixelID], backgroundSigma[pixelID]);
    hPixelDtClean[pixelID] = (TH1D*)hPixelDt[pixelID]->Clone(uniqueName(TString::Format("hCDetAllPixelTimingDtClean_%d", pixelID)));
    for (int bin = 1; bin <= hPixelDtClean[pixelID]->GetNbinsX(); ++bin) {
      hPixelDtClean[pixelID]->SetBinContent(bin, hPixelDt[pixelID]->GetBinContent(bin) - fPixelBackground[pixelID]->Eval(hPixelDt[pixelID]->GetBinCenter(bin)));
      hPixelDtClean[pixelID]->SetBinError(bin, hPixelDt[pixelID]->GetBinError(bin));
    }

    fPixelClean[pixelID] = new TF1(uniqueName(TString::Format("fCDetAllPixelTimingClean_%d", pixelID)), "gaus", FitMin, FitMax);
    fPixelClean[pixelID]->SetParameters(std::max(1.0, fPixelLocal[pixelID]->GetParameter(0)), localMean, localSigma);
    hPixelDtClean[pixelID]->Fit(fPixelClean[pixelID], "RQN0");
#endif
  }

  double referenceWeightSum = 0.0;
  double weightedCentroidSum = 0.0;
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    if (!validFit[pixelID]) continue;
    const double weight = 1.0/(centroidErr[pixelID]*centroidErr[pixelID]);
    referenceWeightSum += weight;
    weightedCentroidSum += centroid[pixelID]*weight;
  }
  const bool detectorReferenceValid = validPixelCount >= 2 && referenceWeightSum > 0.0;
  const double detectorReference = detectorReferenceValid ? weightedCentroidSum/referenceWeightSum : NAN;
  const double detectorReferenceErr = detectorReferenceValid ? std::sqrt(1.0/referenceWeightSum) : NAN;
  if (detectorReferenceValid) {
    for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
      if (!validFit[pixelID]) continue;
      correction[pixelID] = centroid[pixelID] - detectorReference;
      const double correlatedVariance = centroidErr[pixelID]*centroidErr[pixelID] - detectorReferenceErr*detectorReferenceErr;
      correctionErr[pixelID] = correlatedVariance > 0.0 ? std::sqrt(correlatedVariance) : std::sqrt(centroidErr[pixelID]*centroidErr[pixelID] + detectorReferenceErr*detectorReferenceErr);
    }
  }

  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    const double existingOffset = gPixelToffsetLoaded && (int)gPixelToffsetCorr.size() == NumCDetPaddles ? gPixelToffsetCorr[pixelID] : 0.0;
    totalOffset[pixelID] = existingOffset + (validFit[pixelID] && detectorReferenceValid ? correction[pixelID] : 0.0);
  }

  std::ofstream calibration(calibrationOutput.Data());
  if (!calibration.is_open()) {
    std::cerr << "[CDet all-pixel timing] ERROR: could not open calibration candidate '" << calibrationOutput << "'.\n";
    return;
  }
  calibration << "# CDet master calibration constants\n# run " << gRunNumber << "\n# calibration_stage " << gCalibrationStage << "\n# events_processed " << gNumEventsInRun << "\n\n[PixelOffsets]\n";
  calibration.setf(std::ios::fixed);
  calibration.precision(6);
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) calibration << pixelID << " " << totalOffset[pixelID] << " " << fitEntries[pixelID] << "\n";
  calibration << "\n[ECalTiming]\np0 " << gECalFitP0 << "\np1 " << gECalFitP1
              << "\ndelta " << gECalDeltaShift
              << "\n\n[TimeWalk]\np1_L1 " << gTimeWalkP1_L1 << "\np1_L2 " << gTimeWalkP1_L2 << "\ntotref_L1 " << gTimeWalkTotRef_L1 << "\ntotref_L2 " << gTimeWalkTotRef_L2 << "\ntotmin " << gTimeWalkTotMin << "\ntotmax " << gTimeWalkTotMax << "\n";
  calibration.close();

  std::ofstream fitResults(fitResultsOutput.Data());
  if (!fitResults.is_open()) {
    std::cerr << "[CDet all-pixel timing] ERROR: calibration candidate was written, but fit-results file '" << fitResultsOutput << "' could not be opened.\n";
    return;
  }
  fitResults << "# CDet all-pixel ECal-CDet timing-offset extraction\n# run " << gRunNumber << "\n# calibration_stage " << gCalibrationStage << "\n# timing_units ns\n# ecal_energy_cut_gev " << ECalEnergyMin << " " << ECalEnergyMax << "\n# dt = tECal - tCDetLE\n# fit_interval_ns " << FitMin << " " << FitMax << "\n# calibration_centroid_source gaus(0)+pol1(3)_signal_mean\n# peak_seed_interval_ns " << PeakSeedMin << " " << PeakSeedMax << "\n# reference_scope detector-wide\n# detector_mu0_ns " << detectorReference << " detector_mu0_err_ns " << detectorReferenceErr << "\n# residual_correction = mu_i - detector_mu0; total_offset = existing_offset + residual_correction\n# pixel_id bar entries fit_status valid_fit detector_reference_valid mu_ns mu_err_ns sigma_ns sigma_err_ns residual_correction_ns correction_err_ns total_offset_ns detector_mu0_ns detector_mu0_err_ns chi2 ndf chi2_ndf amplitude amplitude_err validity_code failure_reason\n";
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    const int bar = pixelID/NumPaddles;
    fitResults << pixelID << " " << bar << " " << fitEntries[pixelID] << " " << fitStatus[pixelID] << " " << validFit[pixelID] << " " << detectorReferenceValid << " " << centroid[pixelID] << " " << centroidErr[pixelID] << " " << sigma[pixelID] << " " << sigmaErr[pixelID] << " " << correction[pixelID] << " " << correctionErr[pixelID] << " " << totalOffset[pixelID] << " " << detectorReference << " " << detectorReferenceErr << " " << chi2[pixelID] << " " << ndf[pixelID] << " " << chi2Ndf[pixelID] << " " << amplitude[pixelID] << " " << amplitudeErr[pixelID] << " " << validityCode[pixelID] << " \"" << failureReason[pixelID] << "\"\n";
  }
  fitResults.close();
  gPixelToffsetCorr = totalOffset;
  gPixelToffsetNhits = fitEntries;
  gPixelToffsetLoaded = true;
  gCalibrationLoaded = true;
  gLastCalibrationFitSucceeded = true;

  std::cout << "[CDet all-pixel timing]\n"
            << "  events passing ECal energy cut: " << energySelectedEventCount << " / " << nEvents << "\n"
            << "  valid pixel fits: " << validPixelCount << " / " << NumCDetPaddles << "\n"
            << "  calibration centroids: gaus(0)+pol1(3) signal means\n"
            << "  detector-wide weighted centroid: " << detectorReference << " +/- " << detectorReferenceErr << " ns\n"
            << "  calibration file: " << calibrationOutput << "\n"
            << "  detailed fit results: " << fitResultsOutput << "\n"
            << "  existing loaded offsets were " << (hadLoadedPixelOffsets ? "retained and incremented" : "not present; constants start from zero") << "\n";
}

// Diagnostic alternative to extractAllCDetPixelTimingOffsets().  It uses the
// common timing of each contiguous eight-pixel MAPMT group to keep individual
// fits from locking onto a larger cross-talk peak.  It deliberately writes to
// separate candidate/result files and does not modify the active in-memory or
// production calibration.
void extractHierarchicalCDetPixelTimingOffsetsDiagnostic(
    double Width = 1.0,
    double HistMin = -60.0, double HistMax = 30.0,
    double BroadFitMin = -30.0, double BroadFitMax = 10.0,
    int minPixelEntries = 35, int minGroupEntries = 100,
    double minSigma = 0.5, double maxSigma = 8.0,
    double maxChi2Ndf = 15.0, double maxCentroidError = 2.0,
    double minSignalSignificance = 1.5,
    double minimumHalfWindow = 4.0, double maximumHalfWindow = 8.0,
    double ECalEnergyMin = 1.0, double ECalEnergyMax = 12.0,
    TString candidateOutput = "CDet_calibration_dt_hierarchical_candidate.dat",
    TString resultsOutput = "CDet_pixel_timing_fit_results_hierarchical.dat",
    TString rootOutput = "CDet_pixel_timing_hierarchical_diagnostics.root",
    TString plotDirectory = "tdcPlots/hierarchical",
    TString diagnosticBars = "29,79,104,118,139,148",
    bool activateCalibration = false,
    TString pixelCutFile = "") {
  TH1::AddDirectory(kFALSE);
  gLastCalibrationFitSucceeded = false;

  if (Width <= 0.0 || HistMin >= HistMax || BroadFitMin >= BroadFitMax ||
      BroadFitMin < HistMin || BroadFitMax > HistMax || minPixelEntries < 1 ||
      minGroupEntries < 1 || minSigma <= 0.0 ||
      minSigma >= maxSigma || maxChi2Ndf <= 0.0 || maxCentroidError <= 0.0 ||
      minSignalSignificance <= 0.0 ||
      minimumHalfWindow <= 0.0 || minimumHalfWindow > maximumHalfWindow ||
      ECalEnergyMin >= ECalEnergyMax) {
    std::cerr << "[CDet hierarchical timing] ERROR: invalid diagnostic argument.\n";
    return;
  }
  if (!activateCalibration && candidateOutput == gCalibrationFile.c_str()) {
    std::cerr << "[CDet hierarchical timing] ERROR: refusing to overwrite active calibration '"
              << gCalibrationFile << "'.\n";
    return;
  }

  const size_t nEvents = vGoodLe.size();
  if (vGoodID.size() != nEvents || vGoodTot.size() != nEvents ||
      v_GoodECalAdcTime.size() != nEvents ||
      v_GoodECalE.size() != nEvents) {
    std::cerr << "[CDet hierarchical timing] ERROR: event vectors are not aligned.\n";
    return;
  }
  for (size_t ev = 0; ev < nEvents; ++ev) {
    if (vGoodID[ev].size() != vGoodLe[ev].size() ||
        vGoodTot[ev].size() != vGoodLe[ev].size()) {
      std::cerr << "[CDet hierarchical timing] ERROR: LE, TOT, and ID vectors differ in event "
                << ev << ".\n";
      return;
    }
  }

  const int nBins = (int)((HistMax - HistMin)/Width);
  if (nBins < 1) {
    std::cerr << "[CDet hierarchical timing] ERROR: histogram binning is empty.\n";
    return;
  }

  static unsigned long invocation = 0;
  const unsigned long tag = ++invocation;
  auto uname = [tag](const char *base) {
    return TString::Format("%s_hier_%lu", base, tag);
  };
  auto medianOf = [](std::vector<double> values) -> double {
    if (values.empty()) return NAN;
    const size_t middle = values.size()/2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size()%2) return upper;
    std::nth_element(values.begin(), values.begin() + middle - 1, values.end());
    return 0.5*(upper + values[middle - 1]);
  };

  std::vector<TH1D*> pixelHist(NumCDetPaddles, nullptr);
  std::vector<TCutG*> pixelLeTotCut(NumCDetPaddles, nullptr);
  int loadedPixelCuts = 0;
  if (!pixelCutFile.IsNull() && !gSystem->AccessPathName(pixelCutFile)) {
    TFile cutInput(pixelCutFile, "READ");
    if (!cutInput.IsZombie()) {
      for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
        pixelLeTotCut[pixel] = LoadCDetPixelLeTotCut(cutInput, pixel);
        if (pixelLeTotCut[pixel]) ++loadedPixelCuts;
      }
    }
    cutInput.Close();
  }
  if (loadedPixelCuts > 0) {
    std::cout << "[CDet hierarchical timing] Loaded " << loadedPixelCuts
              << " manual LE-versus-TOT pixel cut(s) from " << pixelCutFile << ".\n";
  } else {
    std::cout << "[CDet hierarchical timing] No manual LE-versus-TOT pixel cuts loaded; "
              << "using the automatic path for every pixel.\n";
  }
  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
    pixelHist[pixel] = new TH1D(
        uname(TString::Format("hHierPixelDt_%d", pixel)),
        TString::Format("Pixel %d;t_{ECal}-t_{CDet,LE} (ns);Hits", pixel),
        nBins, HistMin, HistMax);
    pixelHist[pixel]->Sumw2();
  }

  size_t selectedEvents = 0;
  for (size_t ev = 0; ev < nEvents; ++ev) {
    if (v_GoodECalE[ev] < ECalEnergyMin || v_GoodECalE[ev] > ECalEnergyMax) continue;
    ++selectedEvents;
    const double tECal = v_GoodECalAdcTime[ev];
    for (size_t hit = 0; hit < vGoodLe[ev].size(); ++hit) {
      const int pixel = vGoodID[ev][hit];
      if (pixel >= 0 && pixel < NumCDetPaddles &&
          (!pixelLeTotCut[pixel] ||
           pixelLeTotCut[pixel]->IsInside(vGoodTot[ev][hit], vGoodLe[ev][hit])))
        pixelHist[pixel]->Fill(tECal - vGoodLe[ev][hit]);
    }
  }

  const int groupsPerBar = 2;
  const int nBars = NumCDetPaddles/NumPaddles;
  const int nGroups = nBars*groupsPerBar;
  std::vector<TH1D*> groupHist(nGroups, nullptr);
  std::vector<TF1*> groupFit(nGroups, nullptr);
  std::vector<double> groupMean(nGroups, NAN), groupSigma(nGroups, NAN),
                      groupMeanError(nGroups, NAN);
  std::vector<int> groupContributors(nGroups, 0);
  std::vector<bool> groupValid(nGroups, false);

  for (int group = 0; group < nGroups; ++group) {
    const int firstPixel = (group/2)*16 + (group%2)*8;
    groupHist[group] = new TH1D(
        uname(TString::Format("hHierGroupDt_%d", group)),
        TString::Format("Bar %d pixels %d-%d (raw-count sum);t_{ECal}-t_{CDet,LE} (ns);Hits",
                        group/2, firstPixel, firstPixel + 7),
        nBins, HistMin, HistMax);
    groupHist[group]->Sumw2();
    for (int local = 0; local < 8; ++local) {
      const int pixel = firstPixel + local;
      if (IsUnusedPixel(pixel)) continue;
      const int lo = pixelHist[pixel]->FindBin(BroadFitMin);
      const int hi = pixelHist[pixel]->FindBin(BroadFitMax);
      const double integral = pixelHist[pixel]->Integral(lo, hi);
      if (integral <= 0.0) continue;
      groupHist[group]->Add(pixelHist[pixel]);
      ++groupContributors[group];
    }
    const int groupFitBinMin = groupHist[group]->FindBin(BroadFitMin);
    const int groupFitBinMax = groupHist[group]->FindBin(BroadFitMax);
    if (groupHist[group]->Integral(groupFitBinMin, groupFitBinMax) < minGroupEntries) continue;

    const int lo = groupHist[group]->FindBin(BroadFitMin);
    const int hi = groupHist[group]->FindBin(BroadFitMax);
    int peakBin = lo;
    for (int bin = lo + 1; bin <= hi; ++bin)
      if (groupHist[group]->GetBinContent(bin) > groupHist[group]->GetBinContent(peakBin))
        peakBin = bin;
    const double seed = groupHist[group]->GetBinCenter(peakBin);
    const double background = 0.5*(groupHist[group]->GetBinContent(lo) +
                                   groupHist[group]->GetBinContent(hi));
    groupFit[group] = new TF1(uname(TString::Format("fHierGroup_%d", group)),
                              "gaus(0)+pol1(3)", BroadFitMin, BroadFitMax);
    groupFit[group]->SetParameters(
        std::max(0.01, groupHist[group]->GetBinContent(peakBin) - background),
        seed, 3.0, background, 0.0);
    groupFit[group]->SetParLimits(0, 0.0,
        1.5*std::max(1.0, groupHist[group]->GetMaximum()));
    groupFit[group]->SetParLimits(1, BroadFitMin, BroadFitMax);
    groupFit[group]->SetParLimits(2, minSigma, maxSigma);
    const int status = groupHist[group]->Fit(groupFit[group], "RQN0");
    const double mean = groupFit[group]->GetParameter(1);
    const double sigma = std::fabs(groupFit[group]->GetParameter(2));
    const double meanError = groupFit[group]->GetParError(1);
    groupValid[group] = status == 0 && std::isfinite(mean) &&
        std::isfinite(meanError) && meanError > 0.0 &&
        std::isfinite(sigma) && sigma >= minSigma && sigma <= maxSigma &&
        mean > BroadFitMin + Width && mean < BroadFitMax - Width;
    if (groupValid[group]) {
      groupMean[group] = mean;
      groupSigma[group] = sigma;
      groupMeanError[group] = meanError;
    }
  }

  std::vector<TF1*> broadPixelFit(NumCDetPaddles, nullptr),
                    narrowPixelFit(NumCDetPaddles, nullptr);
  std::vector<double> broadMean(NumCDetPaddles, NAN), usedMean(NumCDetPaddles, NAN),
                      usedMeanError(NumCDetPaddles, NAN), usedSigma(NumCDetPaddles, NAN),
                      correction(NumCDetPaddles, 0.0), totalOffset(NumCDetPaddles, 0.0);
  std::vector<int> fitEntries(NumCDetPaddles, 0), fitStatus(NumCDetPaddles, -1),
                   broadFitStatus(NumCDetPaddles, -1);
  std::vector<std::string> source(NumCDetPaddles, "unavailable"), reason(NumCDetPaddles, "not fitted");
  std::vector<double> individualReferenceCentroids;

  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
    const int group = (pixel/16)*2 + (pixel%16)/8;
    const int lo = pixelHist[pixel]->FindBin(BroadFitMin);
    const int hi = pixelHist[pixel]->FindBin(BroadFitMax);
    fitEntries[pixel] = (int)pixelHist[pixel]->Integral(lo, hi);
    if (IsUnusedPixel(pixel)) {
      source[pixel] = "unused";
      reason[pixel] = "known unused pixel";
      continue;
    }
    // Always attempt the original broad individual fit when statistics permit.
    // It is diagnostic for valid groups and becomes the calibration fallback
    // when the eight-pixel group fit is invalid.
    if (fitEntries[pixel] >= minPixelEntries) {
      int peakBin = lo;
      for (int bin = lo + 1; bin <= hi; ++bin)
        if (pixelHist[pixel]->GetBinContent(bin) > pixelHist[pixel]->GetBinContent(peakBin))
          peakBin = bin;
      broadPixelFit[pixel] = new TF1(uname(TString::Format("fHierBroadPixel_%d", pixel)),
                                     "gaus(0)+pol1(3)", BroadFitMin, BroadFitMax);
      broadPixelFit[pixel]->SetParameters(
          std::max(1.0, pixelHist[pixel]->GetBinContent(peakBin)),
          pixelHist[pixel]->GetBinCenter(peakBin), 3.0, 0.0, 0.0);
      broadPixelFit[pixel]->SetParLimits(0, 0.0,
          1.5*std::max(1.0, pixelHist[pixel]->GetMaximum()));
      broadPixelFit[pixel]->SetParLimits(1, BroadFitMin, BroadFitMax);
      broadPixelFit[pixel]->SetParLimits(2, minSigma, maxSigma);
      broadFitStatus[pixel] = pixelHist[pixel]->Fit(broadPixelFit[pixel], "RQN0");
      broadMean[pixel] = broadPixelFit[pixel]->GetParameter(1);
    }

    // A manual two-dimensional cut has already resolved the population
    // ambiguity. Seed its timing fit from the cleaned histogram itself instead
    // of forcing it toward the eight-pixel group centroid.
    if (pixelLeTotCut[pixel]) {
      const double manualSearchLow = BroadFitMin;
      const double manualSearchHigh = std::min(HistMax, std::max(20.0, BroadFitMax));
      const int manualLo = pixelHist[pixel]->FindBin(manualSearchLow);
      const int manualHi = pixelHist[pixel]->FindBin(manualSearchHigh);
      fitEntries[pixel] = (int)pixelHist[pixel]->Integral(manualLo, manualHi);
      if (fitEntries[pixel] >= minPixelEntries) {
        int peakBin = manualLo;
        for (int bin = manualLo + 1; bin <= manualHi; ++bin)
          if (pixelHist[pixel]->GetBinContent(bin) > pixelHist[pixel]->GetBinContent(peakBin))
            peakBin = bin;
        const double seed = pixelHist[pixel]->GetBinCenter(peakBin);
        const double fitLow = std::max(manualSearchLow, seed - maximumHalfWindow);
        const double fitHigh = std::min(manualSearchHigh, seed + maximumHalfWindow);
        const double background = 0.5*(pixelHist[pixel]->GetBinContent(
                                            pixelHist[pixel]->FindBin(fitLow)) +
                                        pixelHist[pixel]->GetBinContent(
                                            pixelHist[pixel]->FindBin(fitHigh)));
        narrowPixelFit[pixel] = new TF1(
            uname(TString::Format("fHierManualPixel_%d", pixel)),
            "gaus(0)+pol1(3)", fitLow, fitHigh);
        narrowPixelFit[pixel]->SetParameters(
            std::max(1.0, pixelHist[pixel]->GetBinContent(peakBin) - background),
            seed, 3.0, background, 0.0);
        narrowPixelFit[pixel]->SetParLimits(
            0, 0.0, 1.5*std::max(1.0, pixelHist[pixel]->GetMaximum()));
        narrowPixelFit[pixel]->SetParLimits(1, fitLow, fitHigh);
        narrowPixelFit[pixel]->SetParLimits(2, minSigma, maxSigma);
        fitStatus[pixel] = pixelHist[pixel]->Fit(narrowPixelFit[pixel], "RQN0");

        const double fitMean = narrowPixelFit[pixel]->GetParameter(1);
        const double fitMeanError = narrowPixelFit[pixel]->GetParError(1);
        const double fitSigma = std::fabs(narrowPixelFit[pixel]->GetParameter(2));
        const double amplitude = narrowPixelFit[pixel]->GetParameter(0);
        const double amplitudeError = narrowPixelFit[pixel]->GetParError(0);
        const int ndf = narrowPixelFit[pixel]->GetNDF();
        const double chi2Ndf = ndf > 0
            ? narrowPixelFit[pixel]->GetChisquare()/ndf : NAN;
        const double significance = amplitudeError > 0.0
            ? amplitude/amplitudeError : NAN;
        const bool nearBoundary = fitMean - fitLow <= 0.5*Width ||
                                  fitHigh - fitMean <= 0.5*Width;
        const bool manualValid = fitStatus[pixel] == 0 && std::isfinite(fitMean) &&
            std::isfinite(fitMeanError) && fitMeanError > 0.0 &&
            fitMeanError <= maxCentroidError && std::isfinite(fitSigma) &&
            fitSigma >= minSigma && fitSigma <= maxSigma && !nearBoundary &&
            amplitude > 0.0 && std::isfinite(significance) &&
            significance >= minSignalSignificance && ndf > 0 &&
            std::isfinite(chi2Ndf) && chi2Ndf <= maxChi2Ndf;
        if (manualValid) {
          source[pixel] = "manual_2d_fit";
          reason[pixel] = "valid fit seeded from manually selected two-dimensional population";
          usedMean[pixel] = fitMean;
          usedMeanError[pixel] = fitMeanError;
          usedSigma[pixel] = fitSigma;
          individualReferenceCentroids.push_back(fitMean);
          continue;
        }
      }

      if (groupValid[group]) {
        source[pixel] = "group_fallback";
        reason[pixel] = fitEntries[pixel] < minPixelEntries
            ? "manual selection has insufficient individual entries"
            : "manual-population fit failed quality criteria";
        usedMean[pixel] = groupMean[group];
        usedMeanError[pixel] = groupMeanError[group];
        usedSigma[pixel] = groupSigma[group];
      } else {
        reason[pixel] = fitEntries[pixel] < minPixelEntries
            ? "manual selection and invalid group have insufficient entries"
            : "manual-population fit failed and group fit is invalid";
      }
      continue;
    }

    if (!groupValid[group]) {
      if (fitEntries[pixel] < minPixelEntries) {
        reason[pixel] = "invalid group fit and insufficient individual entries";
        continue;
      }
      const double fitMean = broadPixelFit[pixel]->GetParameter(1);
      const double fitMeanError = broadPixelFit[pixel]->GetParError(1);
      const double fitSigma = std::fabs(broadPixelFit[pixel]->GetParameter(2));
      const double amplitude = broadPixelFit[pixel]->GetParameter(0);
      const double amplitudeError = broadPixelFit[pixel]->GetParError(0);
      const int ndf = broadPixelFit[pixel]->GetNDF();
      const double chi2Ndf = ndf > 0 ? broadPixelFit[pixel]->GetChisquare()/ndf : NAN;
      const double significance = amplitudeError > 0.0 ? amplitude/amplitudeError : NAN;
      const bool nearBoundary = fitMean - BroadFitMin <= 0.5*Width ||
                                BroadFitMax - fitMean <= 0.5*Width;
      const bool broadValid = broadFitStatus[pixel] == 0 && std::isfinite(fitMean) &&
          std::isfinite(fitMeanError) && fitMeanError > 0.0 && fitMeanError <= maxCentroidError &&
          std::isfinite(fitSigma) && fitSigma >= minSigma && fitSigma <= maxSigma &&
          !nearBoundary && amplitude > 0.0 && std::isfinite(significance) &&
          significance >= minSignalSignificance && ndf > 0 && std::isfinite(chi2Ndf) &&
          chi2Ndf <= maxChi2Ndf;
      fitStatus[pixel] = broadFitStatus[pixel];
      if (broadValid) {
        source[pixel] = "individual_broad_fallback";
        reason[pixel] = "valid broad individual fit after invalid group fit";
        usedMean[pixel] = fitMean;
        usedMeanError[pixel] = fitMeanError;
        usedSigma[pixel] = fitSigma;
        individualReferenceCentroids.push_back(fitMean);
      } else {
        reason[pixel] = "invalid group fit and broad individual fit failed quality criteria";
      }
      continue;
    }

    if (fitEntries[pixel] < minPixelEntries) {
      source[pixel] = "group_fallback";
      reason[pixel] = "insufficient individual entries";
      usedMean[pixel] = groupMean[group];
      usedMeanError[pixel] = groupMeanError[group];
      usedSigma[pixel] = groupSigma[group];
      continue;
    }

    const double halfWindow = std::min(maximumHalfWindow,
        std::max(minimumHalfWindow, 2.0*groupSigma[group]));
    const double fitLow = std::max(BroadFitMin, groupMean[group] - halfWindow);
    const double fitHigh = std::min(BroadFitMax, groupMean[group] + halfWindow);
    narrowPixelFit[pixel] = new TF1(uname(TString::Format("fHierNarrowPixel_%d", pixel)),
                                    "gaus(0)+pol1(3)", fitLow, fitHigh);
    const int groupBin = pixelHist[pixel]->FindBin(groupMean[group]);
    narrowPixelFit[pixel]->SetParameters(
        std::max(1.0, pixelHist[pixel]->GetBinContent(groupBin)),
        groupMean[group], std::min(3.0, groupSigma[group]), 0.0, 0.0);
    narrowPixelFit[pixel]->SetParLimits(0, 0.0, 1.5*std::max(1.0, pixelHist[pixel]->GetMaximum()));
    narrowPixelFit[pixel]->SetParLimits(1, fitLow, fitHigh);
    narrowPixelFit[pixel]->SetParLimits(2, minSigma, maxSigma);
    fitStatus[pixel] = pixelHist[pixel]->Fit(narrowPixelFit[pixel], "RQN0");

    const double fitMean = narrowPixelFit[pixel]->GetParameter(1);
    const double fitMeanError = narrowPixelFit[pixel]->GetParError(1);
    const double fitSigma = std::fabs(narrowPixelFit[pixel]->GetParameter(2));
    const double amplitude = narrowPixelFit[pixel]->GetParameter(0);
    const double amplitudeError = narrowPixelFit[pixel]->GetParError(0);
    const int ndf = narrowPixelFit[pixel]->GetNDF();
    const double chi2Ndf = ndf > 0 ? narrowPixelFit[pixel]->GetChisquare()/ndf : NAN;
    const double significance = amplitudeError > 0.0 ? amplitude/amplitudeError : NAN;
    const bool nearBoundary = fitMean - fitLow <= 0.5*Width || fitHigh - fitMean <= 0.5*Width;
    const bool individualValid = fitStatus[pixel] == 0 && std::isfinite(fitMean) &&
        std::isfinite(fitMeanError) && fitMeanError > 0.0 && fitMeanError <= maxCentroidError &&
        std::isfinite(fitSigma) && fitSigma >= minSigma && fitSigma <= maxSigma &&
        !nearBoundary &&
        amplitude > 0.0 && std::isfinite(significance) && significance >= minSignalSignificance &&
        ndf > 0 && std::isfinite(chi2Ndf) && chi2Ndf <= maxChi2Ndf;
    if (individualValid) {
      source[pixel] = "individual_fit";
      reason[pixel] = "valid constrained fit";
      usedMean[pixel] = fitMean;
      usedMeanError[pixel] = fitMeanError;
      usedSigma[pixel] = fitSigma;
      individualReferenceCentroids.push_back(fitMean);
    } else {
      source[pixel] = "group_fallback";
      reason[pixel] = "individual constrained fit failed quality criteria";
      usedMean[pixel] = groupMean[group];
      usedMeanError[pixel] = groupMeanError[group];
      usedSigma[pixel] = groupSigma[group];
    }
  }

  double detectorReference = medianOf(individualReferenceCentroids);
  if (!std::isfinite(detectorReference)) {
    std::vector<double> validGroups;
    for (int group = 0; group < nGroups; ++group)
      if (groupValid[group]) validGroups.push_back(groupMean[group]);
    detectorReference = medianOf(validGroups);
  }
  if (!std::isfinite(detectorReference)) {
    std::cerr << "[CDet hierarchical timing] ERROR: no valid detector reference.\n";
    return;
  }

  int manualFitCount = 0, individualCount = 0, broadFallbackCount = 0, fallbackCount = 0,
      retainedCount = 0, unavailableCount = 0;
  if (!gAnalysisPixelToffsetSnapshotValid ||
      (int)gAnalysisPixelToffsetCorr.size() != NumCDetPaddles) {
    std::cerr << "[CDet hierarchical timing] ERROR: no valid analysis-time pixel-offset "
              << "snapshot is available. Rerun the main analysis with the desired "
              << "baseline calibration before extracting offsets.\n";
    return;
  }
  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
    const double existing = gAnalysisPixelToffsetCorr[pixel];
    if (source[pixel] == "manual_2d_fit" || source[pixel] == "individual_fit" ||
        source[pixel] == "individual_broad_fallback" ||
        source[pixel] == "group_fallback") {
      correction[pixel] = usedMean[pixel] - detectorReference;
      totalOffset[pixel] = existing + correction[pixel];
      if (source[pixel] == "manual_2d_fit") ++manualFitCount;
      else if (source[pixel] == "individual_fit") ++individualCount;
      else if (source[pixel] == "individual_broad_fallback") ++broadFallbackCount;
      else ++fallbackCount;
    } else {
      totalOffset[pixel] = existing;
      if (gPixelToffsetLoaded || existing != 0.0) {
        if (source[pixel] != "unused") source[pixel] = "retained_existing";
        ++retainedCount;
      } else {
        ++unavailableCount;
      }
    }
  }

  std::ofstream candidate(candidateOutput.Data());
  std::ofstream results(resultsOutput.Data());
  if (!candidate || !results) {
    std::cerr << "[CDet hierarchical timing] ERROR: could not open diagnostic output files.\n";
    return;
  }
  candidate << (activateCalibration
                    ? "# CDet hierarchical eight-pixel-group timing calibration\n"
                    : "# DIAGNOSTIC hierarchical eight-pixel-group CDet calibration candidate\n")
            << "# run " << gRunNumber << "\n# calibration_stage " << gCalibrationStage
            << "\n# detector_reference_median_ns " << detectorReference << "\n\n[PixelOffsets]\n";
  candidate.setf(std::ios::fixed); candidate.precision(6);
  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel)
    candidate << pixel << " " << totalOffset[pixel] << " " << fitEntries[pixel] << "\n";
  candidate << "\n[ECalTiming]\np0 " << gECalFitP0 << "\np1 " << gECalFitP1
            << "\ndelta " << gECalDeltaShift
            << "\n\n[TimeWalk]\np1_L1 " << gTimeWalkP1_L1 << "\np1_L2 " << gTimeWalkP1_L2
            << "\ntotref_L1 " << gTimeWalkTotRef_L1 << "\ntotref_L2 " << gTimeWalkTotRef_L2
            << "\ntotmin " << gTimeWalkTotMin << "\ntotmax " << gTimeWalkTotMax << "\n";

  results << "# CDet hierarchical eight-pixel-group timing diagnostic\n"
          << "# reference median of valid individual centroids\n"
          << "# manual LE-versus-TOT cuts: " << pixelCutFile << " (loaded " << loadedPixelCuts << ")\n"
          << "# pixel bar group entries manual_2d_cut source broad_mu_ns group_mu_ns used_mu_ns used_mu_err_ns used_sigma_ns correction_ns total_offset_ns fit_status reason\n";
  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
    const int group = (pixel/16)*2 + (pixel%16)/8;
    results << pixel << " " << pixel/16 << " " << group << " " << fitEntries[pixel]
            << " " << (pixelLeTotCut[pixel] ? 1 : 0)
            << " " << source[pixel] << " " << broadMean[pixel] << " " << groupMean[group]
            << " " << usedMean[pixel] << " " << usedMeanError[pixel] << " " << usedSigma[pixel]
            << " " << correction[pixel] << " " << totalOffset[pixel] << " " << fitStatus[pixel]
            << " \"" << reason[pixel] << "\"\n";
  }
  candidate.close();
  results.close();

  if (activateCalibration) {
    gPixelToffsetCorr = totalOffset;
    gPixelToffsetNhits = fitEntries;
    gPixelToffsetLoaded = true;
    gCalibrationLoaded = true;
  }

  TFile diagnosticFile(rootOutput, "RECREATE");
  if (!diagnosticFile.IsZombie()) {
    diagnosticFile.mkdir("groups");
    diagnosticFile.cd("groups");
    for (int group = 0; group < nGroups; ++group) {
      if (!groupHist[group]) continue;
      if (groupFit[group]) {
        groupFit[group]->SetLineColor(groupValid[group] ? kGreen + 2 : kMagenta + 1);
        groupFit[group]->SetLineWidth(2);
      }
      groupHist[group]->Write();
    }

    diagnosticFile.mkdir("pixels");
    diagnosticFile.cd("pixels");
    for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
      if (!pixelHist[pixel] || pixelHist[pixel]->GetEntries() <= 0) continue;
      if (broadPixelFit[pixel]) {
        broadPixelFit[pixel]->SetLineColor(kBlue + 1);
        broadPixelFit[pixel]->SetLineStyle(2);
        broadPixelFit[pixel]->SetLineWidth(2);
      }
      if (narrowPixelFit[pixel]) {
        narrowPixelFit[pixel]->SetLineColor(source[pixel] == "individual_fit" ? kRed + 1 : kOrange + 7);
        narrowPixelFit[pixel]->SetLineWidth(2);
      }
      pixelHist[pixel]->Write();
    }
    diagnosticFile.Close();
  } else {
    std::cerr << "[CDet hierarchical timing] WARNING: could not create ROOT diagnostic file '"
              << rootOutput << "'.\n";
  }

  gSystem->mkdir(plotDirectory, kTRUE);
  std::vector<int> plottedBars;
  std::stringstream barStream(diagnosticBars.Data());
  std::string token;
  while (std::getline(barStream, token, ',')) {
    std::stringstream parser(token);
    int bar = -1;
    if (!(parser >> bar) || bar < 0 || bar >= nBars) continue;
    TCanvas *canvas = new TCanvas(
        uname(TString::Format("cHierBar%d", bar)),
        TString::Format("Hierarchical timing diagnostic, bar %d", bar), 1600, 1000);
    canvas->Divide(4, 4);
    for (int local = 0; local < 16; ++local) {
      const int pixel = bar*16 + local;
      const int group = bar*2 + local/8;
      canvas->cd(local + 1);
      pixelHist[pixel]->SetLineColor(kBlack);
      pixelHist[pixel]->Draw("HIST");
      if (broadPixelFit[pixel]) {
        broadPixelFit[pixel]->SetLineColor(kBlue + 1);
        broadPixelFit[pixel]->SetLineStyle(2);
        broadPixelFit[pixel]->Draw("SAME");
      }
      if (narrowPixelFit[pixel]) {
        narrowPixelFit[pixel]->SetLineColor(kRed + 1);
        narrowPixelFit[pixel]->Draw("SAME");
      }
      if (groupValid[group]) {
        TLine *line = new TLine(groupMean[group], 0.0, groupMean[group],
                                std::max(1.0, pixelHist[pixel]->GetMaximum()));
        line->SetLineColor(kGreen + 2);
        line->SetLineWidth(2);
        line->Draw("SAME");
      }
      TLatex *label = new TLatex(0.12, 0.84,
          TString::Format("%s", source[pixel].c_str()));
      label->SetNDC(); label->SetTextSize(0.06); label->Draw();
    }
    canvas->SaveAs(TString::Format("%s/bar%03d_hierarchical.pdf", plotDirectory.Data(), bar));
    canvas->SaveAs(TString::Format("%s/bar%03d_hierarchical.png", plotDirectory.Data(), bar));
    plottedBars.push_back(bar);
  }

  // Store self-contained histogram copies with their fit functions attached.
  // Persisting the interactive canvases themselves is fragile: their pads retain
  // pointers to process-owned objects and can crash when Draw() rebuilds a canvas
  // in a later ROOT session.  A histogram owns its attached functions when read
  // back from a TFile, so these overlays are safe to browse and draw directly.
  TFile overlayFile(rootOutput, "UPDATE");
  if (!overlayFile.IsZombie()) {
    if (!overlayFile.GetDirectory("fit_overlays")) overlayFile.mkdir("fit_overlays");
    for (int bar : plottedBars) {
      TDirectory *overlayRoot = overlayFile.GetDirectory("fit_overlays");
      const TString barDirectory = TString::Format("bar%03d", bar);
      if (!overlayRoot->GetDirectory(barDirectory)) overlayRoot->mkdir(barDirectory);
      TDirectory *barDir = overlayRoot->GetDirectory(barDirectory);
      barDir->cd();

      for (int localGroup = 0; localGroup < 2; ++localGroup) {
        const int group = bar*2 + localGroup;
        if (!groupHist[group]) continue;
        TH1D *groupOverlay = static_cast<TH1D*>(groupHist[group]->Clone(
            TString::Format("group_%d_with_fit", group)));
        groupOverlay->SetDirectory(nullptr);
        if (groupFit[group]) {
          TF1 *fitCopy = static_cast<TF1*>(groupFit[group]->Clone(
              TString::Format("group_%d_fit", group)));
          fitCopy->SetLineColor(groupValid[group] ? kGreen + 2 : kMagenta + 1);
          fitCopy->SetLineWidth(2);
          groupOverlay->GetListOfFunctions()->Add(fitCopy);
        }
        groupOverlay->Write();
      }

      for (int local = 0; local < 16; ++local) {
        const int pixel = bar*16 + local;
        if (!pixelHist[pixel]) continue;
        TH1D *pixelOverlay = static_cast<TH1D*>(pixelHist[pixel]->Clone(
            TString::Format("pixel_%04d_with_fits", pixel)));
        pixelOverlay->SetDirectory(nullptr);
        if (broadPixelFit[pixel]) {
          TF1 *fitCopy = static_cast<TF1*>(broadPixelFit[pixel]->Clone(
              TString::Format("pixel_%04d_broad_fit", pixel)));
          fitCopy->SetLineColor(kBlue + 1);
          fitCopy->SetLineStyle(2);
          fitCopy->SetLineWidth(2);
          pixelOverlay->GetListOfFunctions()->Add(fitCopy);
        }
        if (narrowPixelFit[pixel]) {
          TF1 *fitCopy = static_cast<TF1*>(narrowPixelFit[pixel]->Clone(
              TString::Format("pixel_%04d_constrained_fit", pixel)));
          fitCopy->SetLineColor(kRed + 1);
          fitCopy->SetLineWidth(2);
          pixelOverlay->GetListOfFunctions()->Add(fitCopy);
        }
        pixelOverlay->Write();
      }
    }
    overlayFile.Close();
  }

  gLastCalibrationFitSucceeded = true;
  std::cout << (activateCalibration
                    ? "[CDet hierarchical timing calibration]\n"
                    : "[CDet hierarchical timing diagnostic]\n")
            << "  ECal-selected events: " << selectedEvents << " / " << nEvents << "\n"
            << "  manual LE-versus-TOT cuts applied: " << loadedPixelCuts << "\n"
            << "  robust detector reference: " << detectorReference << " ns\n"
            << "  manual-population individual fits: " << manualFitCount << "\n"
            << "  individual constrained fits: " << individualCount << "\n"
            << "  broad individual fits after invalid groups: " << broadFallbackCount << "\n"
            << "  eight-pixel group fallbacks: " << fallbackCount << "\n"
            << "  retained existing/unavailable: " << retainedCount << " / " << unavailableCount << "\n"
            << (activateCalibration ? "  active calibration: " :
                                      "  candidate calibration (not activated): ")
            << candidateOutput << "\n"
            << "  detailed results: " << resultsOutput << "\n"
            << "  ROOT diagnostics: " << rootOutput << "\n";
}

// Production entry point for the tuned hierarchical pixel-offset method.  The
// output tag keeps all three passes of the full workflow independently
// inspectable while the calibration constants themselves are written to the
// active gCalibrationFile and loaded into the current ROOT session.
void extractHierarchicalCDetPixelTimingOffsets(
    bool generateOffsets = false,
    TString outputTag = "hierarchical",
    TString pixelCutFile = "") {
  gLastCalibrationFitSucceeded = false;
  if (!generateOffsets) {
    std::cout << "[CDet hierarchical timing] Offset generation disabled. Call "
              << "extractHierarchicalCDetPixelTimingOffsets(true) to fit pixels "
              << "and write updated calibration constants.\n";
    return;
  }

  outputTag.ReplaceAll(" ", "_");
  outputTag.ReplaceAll("/", "_");
  if (outputTag.IsNull()) outputTag = "hierarchical";
  const TString resultsOutput = TString::Format(
      "CDet_pixel_timing_fit_results_%s.dat", outputTag.Data());
  const TString rootOutput = TString::Format(
      "CDet_pixel_timing_%s.root", outputTag.Data());
  const TString plotDirectory = TString::Format(
      "tdcPlots/%s", outputTag.Data());

  extractHierarchicalCDetPixelTimingOffsetsDiagnostic(
      1.0, -60.0, 30.0, -30.0, 10.0,
      35, 100, 0.5, 8.0, 15.0, 2.0, 1.5, 4.0, 8.0,
      1.0, 12.0, gCalibrationFile.c_str(), resultsOutput, rootOutput,
      plotDirectory, "29,79,104,118,139,148", true, pixelCutFile);
}

bool ReadCDetPixelOffsetsForComparison(const TString& calibrationFile, std::vector<double>& offsets, std::vector<bool>& found) {
  offsets.assign(NumCDetPaddles, 0.0);
  found.assign(NumCDetPaddles, false);
  std::ifstream input(calibrationFile.Data());
  if (!input.is_open()) {
    std::cerr << "[CDet offset comparison] ERROR: could not open '" << calibrationFile << "'.\n";
    return false;
  }

  std::string line;
  std::string section;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (line[0] == '[') {
      section = line;
      continue;
    }
    if (section != "[PixelOffsets]") continue;
    std::istringstream row(line);
    int pixelID = -1;
    double offset = 0.0;
    if (row >> pixelID >> offset && pixelID >= 0 && pixelID < NumCDetPaddles) {
      offsets[pixelID] = offset;
      found[pixelID] = true;
    }
  }
  const int offsetsRead = std::count(found.begin(), found.end(), true);
  if (offsetsRead != NumCDetPaddles) {
    std::cerr << "[CDet offset comparison] ERROR: incomplete [PixelOffsets] section in '" << calibrationFile << "': " << offsetsRead << " / " << NumCDetPaddles << " rows found.\n";
    return false;
  }
  return true;
}

void plotCDetPixelOffsetMethodDifference(TString regularCalibrationFile = "CDet_calibration.dat", TString dtCalibrationFile = "CDet_calibration_dt.dat") {
  TGraph *gOffsetDifference = new TGraph();
  gOffsetDifference->SetName("gCDetPixelOffsetMethodDifference");
  gOffsetDifference->SetTitle("CDet pixel-offset method comparison;CDet logical pixel ID;#Deltat-method correction - regular correction (ns)");

  std::vector<double> regularOffsets;
  std::vector<double> dtOffsets;
  std::vector<bool> regularFound;
  std::vector<bool> dtFound;
  if (!ReadCDetPixelOffsetsForComparison(regularCalibrationFile, regularOffsets, regularFound) || !ReadCDetPixelOffsetsForComparison(dtCalibrationFile, dtOffsets, dtFound)) return;

  int skippedZero = 0;
  for (int pixelID = 0; pixelID < NumCDetPaddles; ++pixelID) {
    if (!regularFound[pixelID] || !dtFound[pixelID]) continue;
    if (regularOffsets[pixelID] == 0.0 || dtOffsets[pixelID] == 0.0) {
      ++skippedZero;
      continue;
    }
    gOffsetDifference->SetPoint(gOffsetDifference->GetN(), pixelID, dtOffsets[pixelID] - regularOffsets[pixelID]);
  }
  gOffsetDifference->SetMarkerStyle(20);
  gOffsetDifference->SetMarkerSize(0.45);

  TCanvas *cOffsetDifference = new TCanvas("cCDetPixelOffsetMethodDifference", "CDet pixel-offset method comparison", 1200, 700);
  gOffsetDifference->Draw("AP");
  gOffsetDifference->GetXaxis()->SetLimits(-0.5, NumCDetPaddles - 0.5);

  std::cout << "[CDet offset comparison]\n"
            << "  regular calibration: " << regularCalibrationFile << "\n"
            << "  delta-t calibration: " << dtCalibrationFile << "\n"
            << "  plotted pixels: " << gOffsetDifference->GetN() << "\n"
            << "  skipped because either correction was exactly zero: " << skippedZero << "\n";
}

void plotCDetPixelLeAndDtSpectra(int logicalPixelID = 485, double Width = 1.0,
                                 double HistMin = -60.0, double HistMax = 30.0,
                                 double FitMin = -30.0, double FitMax = 0.0,
                                 double ECalEnergyMin = 1.0, double ECalEnergyMax = 12.0,
                                 double LeMin = 0.0, double LeMax = 60.0,
                                 TString regularCalibrationFile = "CDet_calibration.dat",
                                 TString dtCalibrationFile = "CDet_calibration_dt.dat") {
  TH1::AddDirectory(kFALSE);
  if (logicalPixelID < 0 || logicalPixelID >= NumCDetPaddles || Width <= 0.0 || HistMin >= HistMax || FitMin >= FitMax || FitMin < HistMin || FitMax > HistMax || ECalEnergyMin >= ECalEnergyMax || LeMin >= LeMax) {
    std::cerr << "[CDet selected-pixel spectra] ERROR: invalid pixel ID, bin width, or range.\n";
    return;
  }

  const int NLeBins = (int)((LeMax - LeMin)/Width);
  const int NDtBins = (int)((HistMax - HistMin)/Width);
  static unsigned long invocation = 0;
  const unsigned long tag = ++invocation;

  TH1D *hPixelLeSpectrum = new TH1D(TString::Format("hCDetPixelLeSpectrum_%d_%lu", logicalPixelID, tag), TString::Format("CDet good-hit LE, logical pixel ID %d;CDet LE time (ns);Good hits", logicalPixelID), NLeBins, LeMin, LeMax);
  TH1D *hPixelDtSpectrum = new TH1D(TString::Format("hCDetPixelDtSpectrum_%d_%lu", logicalPixelID, tag), TString::Format("ECal-CDet #Deltat, logical pixel ID %d;t_{ECal}-t_{CDet,LE} (ns);Good hits after ECal energy cut", logicalPixelID), NDtBins, HistMin, HistMax);

  for (size_t ev = 0; ev < vGoodLe.size(); ++ev) {
    const bool passesECalEnergyCut = ECalEnergyMin <= v_GoodECalE[ev] && v_GoodECalE[ev] <= ECalEnergyMax;
    for (size_t ihit = 0; ihit < vGoodLe[ev].size(); ++ihit) {
      if (vGoodID[ev][ihit] != logicalPixelID) continue;
      hPixelLeSpectrum->Fill(vGoodLe[ev][ihit]);
      if (passesECalEnergyCut) hPixelDtSpectrum->Fill(v_GoodECalAdcTime[ev] - vGoodLe[ev][ihit]);
    }
  }

  TCanvas *cPixelLeSpectrum = new TCanvas(TString::Format("cCDetPixelLeSpectrum_%d_%lu", logicalPixelID, tag), TString::Format("CDet LE spectrum, logical pixel ID %d", logicalPixelID), 900, 700);
  hPixelLeSpectrum->Draw("HIST");

  TCanvas *cPixelDtSpectrum = new TCanvas(TString::Format("cCDetPixelDtSpectrum_%d_%lu", logicalPixelID, tag), TString::Format("ECal-CDet delta-t spectrum, logical pixel ID %d", logicalPixelID), 900, 700);
  hPixelDtSpectrum->Draw("HIST");
  const double dtLineTop = std::max(1.0, 1.05*hPixelDtSpectrum->GetMaximum());
  TLine *fitMinLine = new TLine(FitMin, 0.0, FitMin, dtLineTop);
  TLine *fitMaxLine = new TLine(FitMax, 0.0, FitMax, dtLineTop);
  fitMinLine->SetLineColor(kRed + 1);
  fitMaxLine->SetLineColor(kRed + 1);
  fitMinLine->SetLineStyle(2);
  fitMaxLine->SetLineStyle(2);
  fitMinLine->SetLineWidth(2);
  fitMaxLine->SetLineWidth(2);
  fitMinLine->Draw("SAME");
  fitMaxLine->Draw("SAME");

  std::vector<double> regularOffsets;
  std::vector<double> dtOffsets;
  std::vector<bool> regularFound;
  std::vector<bool> dtFound;
  if (!ReadCDetPixelOffsetsForComparison(regularCalibrationFile, regularOffsets, regularFound) || !ReadCDetPixelOffsetsForComparison(dtCalibrationFile, dtOffsets, dtFound)) return;

  const double offsetDifference = dtOffsets[logicalPixelID] - regularOffsets[logicalPixelID];
  std::cout << "[CDet selected-pixel spectra]\n"
            << "  logical pixel ID: " << logicalPixelID << "\n"
            << "  LE-spectrum entries: " << hPixelLeSpectrum->GetEntries() << "\n"
            << "  delta-t-spectrum entries after ECal energy cut [" << ECalEnergyMin << ", " << ECalEnergyMax << "] GeV: " << hPixelDtSpectrum->GetEntries() << "\n"
            << "  delta-t fit interval shown: [" << FitMin << ", " << FitMax << "] ns\n"
            << "  regular-method correction from " << regularCalibrationFile << ": " << regularOffsets[logicalPixelID] << " ns\n"
            << "  delta-t-method correction from " << dtCalibrationFile << ": " << dtOffsets[logicalPixelID] << " ns\n"
            << "  delta-t minus regular correction: " << offsetDifference << " ns\n";
}

void plotECalCDetTimeComp(double Width = 1, double diffMinCut = -30, double diffMaxCut = 0,
                          double LeMin = 0.02, double LeMax = 60,
                          double TotMin = 0, double TotMax = 150,
                          double DiffMin = -60, double DiffMax = 30,
                          double CDetTotMin = 0, double CDetTotMax = 80,
                          double CDetMin = 0, double CDetMax = 60,
                          double ECalMin = -40, double ECalMax = 40){
  TH1::AddDirectory(kFALSE);
  int NADCBins = (int)((ECalMax-ECalMin)/4); //4ns bins for ECal, since fADC 4ns resolution
  int TDCBinNum = (int)((DiffMax-DiffMin)/Width);

  TH1D* hECalMinusCDetTime = new TH1D("hECalMinusCDetTime", "ECal-CDet Time;Time Diff (ns);Counts", TDCBinNum, DiffMin, DiffMax);
  TH1D* hECalMinusCDetTimeNoCuts = new TH1D("hECalMinusCDetTimeNoCuts", "ECal-CDet Time;Time Diff (ns);Counts", TDCBinNum, DiffMin, DiffMax);
  TH2D* h2ECalMinusCDetTime = new TH2D("h2ECalMinusCDetTime", "ECal-CDet Time vs CDet Time;CDet 'Good' Time (ns);ECal-CDet Time (ns)", TDCBinNum, CDetMin, CDetMax, TDCBinNum, DiffMin, DiffMax);
  TH2D* h2ECalMinusCDetTot = new TH2D("h2ECalMinusCDetTot", "ECal-CDet Time vs CDet Tot;CDet 'Good' Tot (ns);ECal-CDet Time (ns)", TDCBinNum, CDetTotMin, CDetTotMax, TDCBinNum, DiffMin, DiffMax);
  TH2D* hECalVsCDet = new TH2D("hECalVsCDet", "ECal Time vs CDet Time;CDet LE Time (ns);ECal ADC Time (ns)",TDCBinNum,CDetMin,CDetMax, NADCBins, ECalMin, ECalMax);
  TH2D* h2ECalxVsCDetx = new TH2D("h2ECalxVsCDetx", "ECal Good x vs CDet Good x;CDet Good x (m);ECal Good x (m)",600,-1.5,1.5,200,-1.5,1.5);
  TH2D* h2ECalxVsCDetxNoProject = new TH2D("h2ECalxVsCDetxNoProject", "ECal Actual x vs CDet Good x;CDet Good x (m);ECal Actual x (m)",600,-1.5,1.5,200,-1.5,1.5);
  // TH1I* hGoodHitsPerEvent = new TH1I("hGoodHitsPerEvent", "Good hits per event;N_{good hits};Events", 100, 0, 100);

  const size_t Nev = std::min(vGoodLe.size(),v_GoodECalAdcTime.size());
  int sumNhits = 0;
  int sumGoodHits1 = 0;
  int sumGoodHits2 = 0;
  //vectors for plotting good hits in layers 1 and 2
  // std::vector<int> vGoodHits1PerEvent(Nev, 0);
  // std::vector<int> vGoodHits2PerEvent(Nev, 0);
  // std::vector<int> vGoodHitsPerEvent(Nev, 0);

  for (size_t ev = 0; ev < Nev; ev++) { //iterate through events
    sumNhits += vnhits1[ev]+vnhits2[ev];
    int countGoodHits1 = 0;
    int countGoodHits2 = 0;

    double t_ECal = v_GoodECalAdcTime[ev];
    double x_ECal_actal = v_GoodECalX[ev];

    const size_t Nhits = std::min(vGoodLe[ev].size(), vGoodTot[ev].size());
    for (size_t ihit = 0; ihit < Nhits; ++ihit) {
      double t_CDet = vGoodLe[ev][ihit];
      double tot = vGoodTot[ev][ihit];
      double t_diff = t_ECal-t_CDet;
      hECalMinusCDetTimeNoCuts->Fill(t_diff);

      if (t_CDet >= LeMin && t_CDet <= LeMax && tot >= TotMin && tot <= TotMax && t_diff >= diffMinCut && t_diff <= diffMaxCut){
        double x_CDet = vCDetGoodX[ev][ihit];
        double x_ECal = v_GoodECalX[ev]*vCDetGoodZ[ev][ihit]/ECal_dist;
        hECalMinusCDetTime->Fill(t_diff);
        h2ECalMinusCDetTime->Fill(t_CDet, t_diff);
        h2ECalMinusCDetTot->Fill(tot, t_diff);
        hECalVsCDet->Fill(t_CDet, t_ECal);
        h2ECalxVsCDetx->Fill(x_CDet,x_ECal);
        h2ECalxVsCDetxNoProject->Fill(x_CDet,x_ECal_actal);

        int sbselemid = (Int_t)vGoodID[ev][ihit];
        int sbsrown = sbselemid%672;
        int sbscoln = sbselemid/672;
        int mylayern = sbscoln/2;
        int mypaddlen = sbscoln*672 + sbsrown;

        if (mylayern == 0) countGoodHits1++;
        else countGoodHits2++;

	    }//fill histogram with cuts
    }//finished looking at hits
    // vGoodHitsPerEvent[ev] = countGoodHits1 + countGoodHits2;
    // vGoodHits1PerEvent[ev] = countGoodHits1;
    // vGoodHits2PerEvent[ev] = countGoodHits2;

    sumGoodHits1 += countGoodHits1;
    sumGoodHits2 += countGoodHits2;

  } //finished looking at all events
  double aveNhits = (sumNhits / Nev);
  double aveGoodHits = ((sumGoodHits1+sumGoodHits2) / Nev);
  double aveGoodHits1 = (sumGoodHits1 / Nev);
  double aveGoodHits2 = (sumGoodHits2 / Nev);
  std::cout << "ave nhits/event w/o cuts= " << aveNhits << endl;
  std::cout << "ave nGoodHits/event= " << aveGoodHits << endl;
  std::cout << "ave nGoodHits Layer 1/event= " << aveGoodHits1 << endl;
  std::cout << "ave nGoodHits Layer 2/event= " << aveGoodHits2 << endl;

  //make profile to fit time diff vs CDet time with linear fit
  // TProfile *prof = h2ECalMinusCDetTime->ProfileX("prof");
  // prof->Fit("pol1");

  // TF1 *fit = prof->GetFunction("pol1");  // retrieve automatic fit

  //make canvas and draw hist
  TCanvas *cTimeDiff = new TCanvas("cTimeDiff", "ECal ADCtime Minus CDet Good LE",900,700);

  hECalMinusCDetTime->SetLineColor(kRed);
  hECalMinusCDetTimeNoCuts->SetLineColor(kBlue);
  hECalMinusCDetTimeNoCuts->Draw("HIST");
  hECalMinusCDetTime->Draw("SAME");

  // ---- Gaussian fit on the "NoCuts" histogram ----
  TF1 *fGausNoCuts = new TF1("fGausNoCuts", "gaus", DiffMin, DiffMax);
  int maxBin = hECalMinusCDetTimeNoCuts->GetMaximumBin();
  double peakX = hECalMinusCDetTimeNoCuts->GetBinCenter(maxBin);
  fGausNoCuts->SetParameters(hECalMinusCDetTimeNoCuts->GetMaximum(), peakX, 20.0); // amp, mean, sigma guess

  // (optional) restrict fit range around the peak so you fit the main bump, not the tails/background
  double fitLo = peakX - 40.0;
  double fitHi = peakX + 40.0;
  fGausNoCuts->SetRange(fitLo, fitHi);

  // do the fit ("R" uses the TF1 range, "0" suppresses ROOT fit printout if you want)
  hECalMinusCDetTimeNoCuts->Fit(fGausNoCuts, "R");

  // draw the fit on top of the already-drawn histogram
  fGausNoCuts->SetLineColor(kBlack);
  fGausNoCuts->Draw("SAME");

  std::cout << "NoCuts Gaussian: mean = " << fGausNoCuts->GetParameter(1) << "  sigma = " << fGausNoCuts->GetParameter(2) << std::endl;

  auto leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(hECalMinusCDetTime,"Cuts","l");
  leg->AddEntry(hECalMinusCDetTimeNoCuts,"NoCuts","l");
  leg->Draw();

  //gaussian subtracted hist
  TH1* hNoCutsSub = SubtractFitFromHist(hECalMinusCDetTimeNoCuts, fGausNoCuts);
  hNoCutsSub->SetLineColor(kMagenta);
  hNoCutsSub->Draw("HIST SAME");

  TCanvas *cXComp = new TCanvas("cXComp", "ECal x vs CDet x",900,700);
  cXComp->SetLogz();
  h2ECalxVsCDetx->Draw("COLZ");

  TCanvas *c2DtimeComps = new TCanvas("c2DtimeComps", "ECal-CDet Time Comparisons",900,700);
  c2DtimeComps->Divide(1,3);

  c2DtimeComps->cd(1);
  //gPad->SetLogz();
  hECalVsCDet->Draw("COLZ");

  c2DtimeComps->cd(2);
  //gPad->SetLogz();
  h2ECalMinusCDetTime->Draw("COLZ"); //heatmap

  c2DtimeComps->cd(3);
  //gPad->SetLogz();
  h2ECalMinusCDetTot->Draw("COLZ"); //heatmap

  TCanvas *cXCompActual = new TCanvas("cXCompActual", "ECal x vs CDet x w/o projection",900,700); //uses non projected ecal x
  cXCompActual->SetLogz();
  h2ECalxVsCDetx->Draw("COLZ");

} //end routine

void plotRawXCorrelation(double tDiffMin = 80, double tDiffMax = 100){

  TH2D* h2RawECalxVsCDetx = new TH2D("h2ECalxVsCDetx", "ECal Good x vs CDet Good x;CDet Good x (m);ECal Good x (m)",600,-1.5,1.5,200,-1.5,1.5);

  if (vCDetX.size() != v_ECalX.size()){
    std::cout << "vCDetX and vECalX not the same size" << std::endl;
  }
  const size_t Nev = vCDetX.size();
  for (size_t ev = 0; ev < Nev; ev++){
    double t_e = v_ECalAdcTime[ev];
    //std::cout << "event = " << ev << " " << "size of cdetX = " << vCDetX[ev].size() << std::endl;
    //std::cout << " " <<std::endl;
    //std::cout << "event = " << ev << " " << "size of cdetLE = " << vRawLe[ev].size() << std::endl;
    //std::cout << "event = " << ev << " " << "size of cdetLE = " << vRawTe[ev].size() << std::endl;
    //std::cout << "event = " << ev << " " << "size of cdetLE = " << vRawTot[ev].size() << std::endl;
    const size_t Nhits = vCDetX[ev].size();
    for (size_t ihit = 0; ihit < Nhits; ihit++){
      double t_c = vRawLe[ev][ihit];
      double t_diff = t_e - t_c;
      if (t_diff >= tDiffMin && t_diff <= tDiffMax){
        double x_c = vCDetX[ev][ihit];
        double z_c = vCDetZ[ev][ihit];
        double x_e = v_ECalX[ev]*(vCDetZ[ev][ihit]/ECal_dist);
        //std::cout << "CDetX = " << x_c << " ECalX = " << x_e << "CDetZ = " << z_c << std::endl;
        h2RawECalxVsCDetx->Fill(x_c,x_e);
      }//if statement for cuts
    }//end hit loop
  }//end event loop
  TCanvas *cRawXComp = new TCanvas("cRawXComp", "ECal x vs Raw CDet x",900,700);
  h2RawECalxVsCDetx->Draw("COLZ");
}

void plotTimeECalVsCDet(double Width = 0.0160167/2,
                        double LeMin = 0.02, double LeMax = 60,
                        double TotMin = 0, double TotMax = 150,
                        double CDetMin = 0, double CDetMax = 60,
                        double ECalMin = -40, double ECalMax = 40){
  int NADCBins = (int)((ECalMax-ECalMin)/4); //4ns bins for ECal, since fADC 4ns resolution
  int TDCBinNum = (int)((CDetMax-CDetMin)/Width);
  TH2D* hECalVsCDet = new TH2D("hECalVsCDet", "ECal Time vs CDet Time;CDet LE Time (ns);ECal ADC Time (ns)",TDCBinNum,CDetMin,CDetMax, NADCBins, ECalMin, ECalMax);
  const size_t Nev = std::min(vGoodLe.size(),v_GoodECalAdcTime.size());

  for (size_t ev = 0; ev < Nev; ev++) { //iterate through events
    double t_ECal = v_GoodECalAdcTime[ev];

    const size_t Nhits = std::min(vGoodLe[ev].size(), vGoodTot[ev].size());
    for (size_t ihit = 0; ihit < Nhits; ++ihit) {
      double t_CDet = vGoodLe[ev][ihit];
      double tot = vGoodTot[ev][ihit];
      if (t_CDet >= LeMin && t_CDet <= LeMax && tot >= TotMin && tot <= TotMax){
        hECalVsCDet->Fill(t_CDet, t_ECal);
      }
    }//finished looking at hits
  } //finished looking at all events

  //make canvas and draw hist
  TCanvas *cTimeComp = new TCanvas("cTimeComp", "ECal ADCtime vs CDet LE",900,700);
  hECalVsCDet->Draw("COLZ");
} //end routine

void plotXDiffSections(double le_min = 0, double le_max = 60, double tDiffMin = 80, double tDiffMax = 100){
  //define nonchanging histograms
  TH1D* h1 = new TH1D("h1", "xDiffLayer1;xdiff (m);Counts", NXDiffBins,XDiffLow,XDiffHigh);
  TH1D* h2 = new TH1D("h2", "xDiffLayer2;xdiff (m);Counts", NXDiffBins,XDiffLow,XDiffHigh);
  TH1D* h3 = new TH1D("h3", "good le;le (ns);Counts", NTDCBins, le_min, le_max);

  const size_t Nev = vCDetGoodX.size();
  for (size_t ev = 0; ev < Nev; ++ev) {
    double t_ECal = v_GoodECalAdcTime[ev];
    for (std::size_t n = 0; n < vCDetGoodX[ev].size(); ++n){
      double t_CDet = vGoodLe[ev][n];
      double t_diff = t_ECal-t_CDet;
      if (vGoodLe[ev][n] <= le_max && vGoodLe[ev][n] >= le_min && t_diff >= tDiffMin && t_diff <= tDiffMax){
        if (vGoodLayer[ev][n]==0){ //layer 1
          h1->Fill(vCDetGoodX[ev][n]-(v_GoodECalX[ev]*vCDetGoodZ[ev][n]/ECal_dist));
        }
        else if (vGoodLayer[ev][n]==1){ //layer 2
          h2->Fill(vCDetGoodX[ev][n]-(v_GoodECalX[ev]*vCDetGoodZ[ev][n]/ECal_dist));
        }
        h3->Fill(vGoodLe[ev][n]);
      }
    }
  } // individual time step histograms filled

  TCanvas *cTimeComp = new TCanvas("cTimeComp", "ECal ADCtime vs CDet LE",900,700);
  TCanvas* c1 = new TCanvas("c1", "xDiff Layer 1",900,700);
  h1->Draw("HIST");

  TCanvas* c2 = new TCanvas("c2", "xDiff Layer 2",900,700);
  h2->Draw("HIST");

  TCanvas* c3 = new TCanvas("c3", "le",900,700);
  h3->Draw("HIST");

} //end plotXDiffSections

void plotPMTRates(Int_t mymodule=1, Int_t mylayer=1, Int_t choice = 1){

  Int_t offsetl = (mylayer-1)*1344 + (mymodule-1)*224;
  Int_t offsetr = (mylayer-1)*1344 + 672 + (mymodule-1)*224;
  Int_t xcpos = 50 + (choice-1)*400;

  TCanvas *caPMT = new TCanvas(TString::Format("PMT Rates Left %d",choice), TString::Format("PMT Rates Left %d",choice), xcpos,50,400,550);
  caPMT->Divide(2,7,0.01,0.01,0);
  TCanvas *caPMTT = new TCanvas(TString::Format("PMT Rates Right %d",choice), TString::Format("PMT Rates Right %d",choice), xcpos,650,400,550);
  caPMTT->Divide(2,7,0.01,0.01,0);

  Double_t histymax = 12500;
  Double_t rate_levelu = 1800; // Expect 600 for 50000 events @ 5uA
  Double_t rate_levell = 600; // Expect 600 for 50000 events @ 5uA
  TLine *lineu = new TLine(0, rate_levelu, 2688, rate_levelu);
  TLine *linel = new TLine(0, rate_levell, 2688, rate_levell);

  for (Int_t ii=0; ii<14; ii++) {
    caPMT->cd(ii+1);
    gPad->SetLogy();
    gPad->DrawFrame(offsetl+ii*16,1.,offsetl+ii*16+16,histymax);
    hAllRawPMT->Draw("sames");

    //draw unused pixels on each PMT plot
    double xmin = offsetl + ii * 16;
    double xmax = xmin + 16;
    for (double x : missingPixelBins) {
      if (x < xmin || x >= xmax) continue;

      int bin = hAllRawPMT->FindBin(x);
      double xlow = hAllRawPMT->GetBinLowEdge(bin);
      double xup = xlow + hAllRawPMT->GetBinWidth(bin);
      double yup = hAllRawPMT->GetBinContent(bin);

      TBox *box = new TBox(xlow, 0, xup, yup);
      box->SetFillColor(kBlack);
      box->SetFillStyle(1001);
      box->Draw("sames");
    }

    lineu->Draw("sames");
    linel->Draw("sames");
  }
  caPMT->Update();

  for (Int_t ii=0; ii<14; ii++) {
    caPMTT->cd(ii+1);
    gPad->SetLogy();
    gPad->DrawFrame(offsetr+ii*16,1.,offsetr+ii*16+16.,histymax);
    hAllRawPMT->Draw("sames");

    double xmin = offsetr + ii * 16;
    double xmax = xmin + 16;
    for (double x : missingPixelBins) {
      if (x < xmin || x >= xmax) continue;

      int bin = hAllRawPMT->FindBin(x);
      double xlow = hAllRawPMT->GetBinLowEdge(bin);
      double xup = xlow + hAllRawPMT->GetBinWidth(bin);
      double yup = hAllRawPMT->GetBinContent(bin);

      TBox *boxx = new TBox(xlow, 0, xup, yup);
      boxx->SetFillColor(kBlack);
      boxx->SetFillStyle(1001);
      boxx->Draw("sames");
    }

    lineu->Draw("sames");
    linel->Draw("sames");
  }
  caPMTT->Update();


  return;


}
TCanvas *plotBarRates(){

  TCanvas *cabar = new TCanvas("BarRates", "BarRates", 50,50,800,800);
  cabar->Divide(3,4,0.01,0.01,0);

  Double_t histymax = 200000.0;
  Double_t rate_level = 600; // Expect 600 for 50000 events @ 5uA
  TLine *line = new TLine(0, rate_level, 168, rate_level);

  cabar->cd(1);
  gPad->SetLogy();
  gPad->DrawFrame(0.,1.,14.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(2);
  gPad->SetLogy();
  gPad->DrawFrame(14.,1.,28.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(3);
  gPad->SetLogy();
  gPad->DrawFrame(28.,1.,42.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(4);
  gPad->SetLogy();
  gPad->DrawFrame(42.,1.,56.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(5);
  gPad->SetLogy();
  gPad->DrawFrame(56.,1.,70.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(6);
  gPad->SetLogy();
  gPad->DrawFrame(70.,1.,84.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(7);
  gPad->SetLogy();
  gPad->DrawFrame(84.,1.,98.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(8);
  gPad->SetLogy();
  gPad->DrawFrame(98.,1.,112.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(9);
  gPad->SetLogy();
  gPad->DrawFrame(112.,1.,126.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(10);
  gPad->SetLogy();
  gPad->DrawFrame(126.,1.,140.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(11);
  gPad->SetLogy();
  gPad->DrawFrame(140.,1.,154.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(12);
  gPad->SetLogy();
  gPad->DrawFrame(154.,1.,168.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  return cabar;


}

TCanvas *plotGoodTDC2D(){

  TCanvas *cac = new TCanvas("all2d", "all2d", 50,50,800,800);
  cac->Divide(2,2,0.01,0.01,0);

  cac->cd(1);
  gPad->SetLogz();
  h2AllGoodLe->Draw("colz");
  cac->cd(2);
  gPad->SetLogz();
  h2AllGoodTe->Draw("colz");
  cac->cd(3);
  gPad->SetLogz();
  h2AllGoodTot->Draw("colz");

  return cac;
}

TCanvas *plotRefTDC() {
  hRefRawLe = new TH1F(TString::Format("hRefRawLe"),
            TString::Format("hRefRawLe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefRawTe = new TH1F(TString::Format("hRefRawTe"),
            TString::Format("hRefRawTe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefRawTot = new TH1F(TString::Format("hRefRawTot"),
            TString::Format("hRefRawTot"),
            RefNTotBins, RefTotBinLow, RefTotBinHigh);
  hRefRawPMT = new TH1F(TString::Format("hRefRawPMT"),
            TString::Format("hRefRawPMT"),
            32, 2688, 2720);

  hRefGoodLe = new TH1F(TString::Format("hRefGoodLe"),
            TString::Format("hRefGoodLe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTe = new TH1F(TString::Format("hRefGoodTe"),
            TString::Format("hRefGoodTe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTot = new TH1F(TString::Format("hRefGoodTot"),
            TString::Format("hRefGoodTot"),
            RefNTotBins, RefTotBinLow, RefTotBinHigh);
  hRefGoodPMT = new TH1F(TString::Format("hRefGoodPMT"),
            TString::Format("hRefGoodPMT"),
            32, 2688, 2720);


  //fill histograms
  for (double val : vRefRawLe)   hRefRawLe->Fill(val);
  for (double val : vRefRawTe)   hRefRawTe->Fill(val);
  for (double val : vRefRawTot)  hRefRawTot->Fill(val);
  for (int    val : vRefRawPMT)  hRefRawPMT->Fill(val);

  for (double val : vRefGoodLe)   hRefGoodLe->Fill(val);
  for (double val : vRefGoodTe)   hRefGoodTe->Fill(val);
  for (double val : vRefGoodTot)  hRefGoodTot->Fill(val);
  for (int    val : vRefGoodPMT)  hRefGoodPMT->Fill(val);

  //make canvas
  TCanvas *cbb = new TCanvas("ref", "ref", 850,50, 1200,800);
  cbb->Divide(2,2,0.01,0.01,0);

  cbb->cd(1);
  gPad->SetLogy();
  hRefRawLe->Draw();

  cbb->cd(2);
  gPad->SetLogy();
  hRefRawTe->Draw();

  cbb->cd(3);
  gPad->SetLogy();
  hRefRawTot->Draw();

  cbb->cd(4);
  gPad->SetLogy();
  hRefRawPMT->Draw();

  TCanvas *cbbGood = new TCanvas("refGood", "refGood", 850, 50, 1200, 800);
  cbbGood->Divide(2,2,0.01,0.01,0);

  cbbGood->cd(1);
  gPad->SetLogy();
  hRefGoodLe->Draw();

  cbbGood->cd(2);
  gPad->SetLogy();
  hRefGoodTe->Draw();

  cbbGood->cd(3);
  gPad->SetLogy();
  hRefGoodTot->Draw();

  cbbGood->cd(4);
  gPad->SetLogy();
  hRefGoodPMT->Draw();

  return cbb;
}


TCanvas *plotCDetTDC(){

  TCanvas *canvas[NumHalfModules];
  for (int cmodule=1;cmodule<=NumModules;cmodule++) {
   for (int layer=1;layer<=NumLayers;layer++) {
    for (int side=1;side<=NumSides;side++){
     int xposition1 = 100 + (layer-1)*1000 + (side-1)*500;
     int yposition1 = 100 + (NumModules-cmodule)*400;
     int xposition2 = 150 + (layer-1)*1000 + (side-1)*500;
     int yposition2 = 150 + (NumModules-cmodule)*400;

     int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
     int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;

     TString cname;
     cname.Form("c1_%d_%d_%d",cmodule,layer,side);
     TCanvas *c1 = new TCanvas(cname, cname, xposition1,yposition1,380,280);
     c1->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c1;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c1->cd(ii+1);
        if (ii == 1) {
		cout << "side_group_plot = " << side_group_plot << endl;
	}

	hRawTe[elemID_start + ii ]->Draw();

     }

     cname.Form("c2_%d_%d_%d",cmodule,layer,side);
     TCanvas *c2 = new TCanvas(cname, cname, xposition2,yposition2,380,280);
     c2->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c2;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c2->cd(ii+1);
        if (ii == 1) {
                cout << "side_group_plot = " << side_group_plot << endl;
        }

        hRawLe[elemID_start + ii + 1]->Draw();

     }


    }
   }
  }

  return canvas[0];

}

TCanvas *plotHalfModule(int cmodule = 1, int side = 1, int layer = 1){


  TCanvas *c4 = new TCanvas("c4", "c4", 50,50,1000,1200);
  c4->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
  TCanvas *c4b = new TCanvas("c4b", "c4b", 1050,50,1000,1200);
  c4b->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);

  int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;


  for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c4->cd(ii+1);
  	hRawTe[elemID_start + ii ]->Draw();
	c4b->cd(ii+1);
  	hRawLe[elemID_start + ii ]->Draw();
  }

  return c4;

}

TCanvas *plotBarTDC(int bar = 39, int side = 1, int layer = 1){

  int mymodule = (bar-1)/NumBars+1; //mymodule in this case represnts top, middle, or bottom as opposed todetector labels whe
  int paddle_start = (bar - (mymodule-1)*NumBars - 1)*NumPaddles;
  std::cout << "mymodule = " << mymodule << "  paddle_start = " << paddle_start << std::endl;

  int side_group = 6*(layer-1) + 2*(mymodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (mymodule-1)*NumBars*NumPaddles + paddle_start;


  std::cout << "side_group = " << side_group << std::endl;

  TCanvas *c3 = new TCanvas("c3", "c3", 150,150,600,450);
  c3->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333 = new TCanvas("c333", "c333", 750,150,600,450);
  c333->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c3a = new TCanvas("c3a", "c3a", 150,650,600,450);
  c3a->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333a = new TCanvas("c333a", "c333a", 750,650,600,450);
  c333a->Divide(4,4, 0.01, 0.01, 0);

  for (int ii = 0; ii < NumPaddles; ii++) {

        c3->cd(ii+1);
  	hRawLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c333->cd(ii+1);
  	hGoodLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c3a->cd(ii+1);
  	hRawTe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c333a->cd(ii+1);
  	hGoodTe[elemID_start + ii ]->Draw();
  }

  return c3;

}

TCanvas *plotTOTvsLE(){

  TCanvas *c123 = new TCanvas("c123", "c123", 50,50,1000,1000);

  c123->cd();
  h2TDCTOTvsLE->Draw("colz");

  return c123;

}

TCanvas *plotTOTvsXDiff(){

  TCanvas *c1234 = new TCanvas("c1234", "c1234", 50,50,1000,1000);
  c1234->Divide(2,2,0.01,0.01,0);

  c1234->cd(1);
  h2TOTvsXDiff1->Draw("colz");
  c1234->cd(2);
  h2TOTvsXDiff2->Draw("colz");
  c1234->cd(3);
  h2LEvsXDiff1->Draw("colz");
  c1234->cd(4);
  h2LEvsXDiff2->Draw("colz");

  return c1234;

}


TCanvas *plotRowColLayer(){


  TCanvas *c5 = new TCanvas("c5", "c5", 50,50,800,800);
  c5->Divide(4,2, 0.01, 0.01, 0);

  c5->cd(1);
  gPad->SetLogy();
  hRowLayer1Side1->Draw();
  c5->cd(2);
  gPad->SetLogy();
  hRowLayer1Side2->Draw();
  c5->cd(3);
  gPad->SetLogy();
  hRowLayer2Side1->Draw();
  c5->cd(4);
  gPad->SetLogy();
  hRowLayer2Side2->Draw();
  c5->cd(5);
  gPad->SetLogy();
  hRow->Draw();
  c5->cd(6);
  gPad->SetLogy();
  hCol->Draw();
  c5->cd(7);
  gPad->SetLogy();
  hLayer->Draw();
  c5->cd(8);
  gPad->SetLogy();
  hHitPMT->Draw();

  return c5;

}

TCanvas *plotNhits(){
  TCanvas *c55 = new TCanvas("c55", "c5", 50,50,800,800);
  c55->Divide(3,4, 0.01, 0.01, 0);

  c55->cd(1);
  hnhits1->Draw();
  c55->cd(2);
  hngoodhits1->Draw();
  c55->cd(3);
  hngoodTDChits1->Draw();
  c55->cd(4);
  hnhits2->Draw();
  c55->cd(5);
  hngoodhits2->Draw();
  c55->cd(6);
  hngoodTDChits2->Draw();

  c55->cd(7);
  hnpaddles->Draw();
  c55->cd(8);
  hngoodpaddles->Draw();
  c55->cd(9);
  hngoodTDCpaddles->Draw();

  c55->cd(10);
  hnhits_ev->Draw();
  c55->cd(11);
  hngoodhits_ev->Draw();
  c55->cd(12);
  hngoodTDChits_ev->Draw();

  return c55;

}

TCanvas *plotTDC2d(){

  h2d_RawLE  = new TH2F("h2d_RawLE","Raw LE vs PMT", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTE  = new TH2F("h2d_RawTE","Raw TE vs PMT", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTot = new TH2F("h2d_RawTot","Raw TOT vs PMT", NTotBins,TotBinLow,TotBinHigh,nTdc+1,0,nTdc+1);

  // TH2D *h2d_RawLE = new TH2D("h2d_RawLE", "Raw LE vs PMT", 400, 0, 200, 2700, 0, 2700);
  // TH2D *h2d_RawTE = new TH2D("h2d_RawTE", "Raw TE vs PMT", 400, 0, 200, 2700, 0, 2700);
  // TH2D *h2d_RawTot = new TH2D("h2d_RawTot", "Raw Tot vs PMT", 400, 0, 200, 2700, 0, 2700);

  for (size_t evt = 0; evt < vRawLe.size(); evt++) {
    for (size_t hit = 0; hit < vRawLe[evt].size(); hit++) {
      h2d_RawLE->Fill(vRawLe[evt][hit], vRawID[evt][hit]);
      h2d_RawTE->Fill(vRawTe[evt][hit], vRawID[evt][hit]);
      h2d_RawTot->Fill(vRawTot[evt][hit], vRawID[evt][hit]);
    }
  }



  TCanvas *c6 = new TCanvas("c6", "c6", 50,50,800,800);
  c6->Divide(2,2, 0.01, 0.01, 0);

  c6->cd(1);
  h2d_RawLE->Draw();
  c6->cd(2);
  h2d_RawTE->Draw();
  c6->cd(3);
  h2d_RawTot->Draw();
  c6->cd(4);
  h2d_Mult->Draw();

  return c6;

}

auto *plotEECalCDet() {

  TCanvas *c1717 = new TCanvas("c1717", "c1717", 50,50,800,800);
  c1717->Divide(1,2, 0.01, 0.01, 0);

  c1717->cd(1);
  hEECalCDet1->Draw();
  c1717->cd(2);
  hEECalCDet2->Draw();

  return c1717;
}

auto plotXYZ(){

  hHitX = new TH1F("HitXposition","HitXPosition",1000,-2.0,2.0);
  hHitY = new TH1F("HitYposition","HitYPosition",200,-0.5,0.5);
  hHitZ = new TH1F("HitZposition","HitZPosition",200,5.5,6.0);
   size_t N = vCDetGoodX.size();
    for (size_t ev = 0; ev < N; ev++){
      for (size_t hit = 0; hit < vCDetGoodX[ev].size(); hit++) {
      //if (vRefRawTot[i] >= cutTotMin && vRefRawTot[i] <= cutTotMax) {
        hHitX->Fill(vCDetGoodX[ev][hit]);
        hHitY->Fill(vCDetGoodY[ev][hit]);
        hHitZ->Fill(vCDetGoodZ[ev][hit]);
      //}
      }
    }

   TCanvas *c7 = new TCanvas("c7", "c7", 800,800);
   c7->Draw();
   TPad *p1 = new TPad("p1","p1",0.05,0.0,0.45,1.0);
   p1->Draw();
   p1->Divide(1,3);

   p1->cd(1);
   gPad->SetLogy();
   hHitX->Draw();

   p1->cd(2);
   gPad->SetLogy();
   hHitY->Draw();

   p1->cd(3);
   gPad->SetLogy();
   hHitZ->Draw();

   c7->cd(0);
   TPad *p2 = new TPad("p1","p1",0.55,0.0,0.95,1.0);
   p2->Draw();
   p2->Divide(2,1);

   p2->cd(1);
   gPad->SetLogz();
   hHitXY1->Draw("colz");
   p2->cd(2);
   gPad->SetLogz();
   hHitXY2->Draw("colz");


  return c7;

}

/*auto plotXYECalCDet(){

   TCanvas *c8 = new TCanvas("c8", "c7", 1200,1200);
   c8->Divide(3,3);

   c8->cd(1);
   gPad->SetLogz();
   hXECalCDet1->Draw("colz");

   c8->cd(2);
   gPad->SetLogz();
   hXECalCDet2->Draw("colz");

   c8->cd(3);
   hXECal->Draw();

   c8->cd(4);
   hYECal->Draw();

   c8->cd(5);
    // Define the Gaussian + constant background function
    TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

    // Set initial parameters:
    // [0] amplitude, [1] mean, [2] sigma, [3] constant background
    fitFunc->SetParameters(hXDiffECalCDet1->GetMaximum(), 0.02, 0.01, hXDiffECalCDet1->GetMinimum());
    fitFunc->SetParNames("Amplitude", "Mean", "Sigma", "Background");

    // Optional: set limits on the parameters if needed
    fitFunc->SetParLimits(1, -0.05, 0.05);  // constrain the mean near 0.02
    fitFunc->SetParLimits(2, 0.001, 0.05);  // positive sigma

    // Fit the histogram
    hXDiffECalCDet1->Fit(fitFunc, "R");  // "R" = use function range only

    // Draw the result
    hXDiffECalCDet1->Draw();
    fitFunc->Draw("same");
    // Extract Gaussian parameters
    double A = fitFunc->GetParameter(0); // Amplitude
    double mu = fitFunc->GetParameter(1); // Mean
    double sigma = fitFunc->GetParameter(2); // Sigma
    double bg = fitFunc->GetParameter(3); // Background level

    // Define integration limits
    double x_min = mu - 3*sigma;
    double x_max = mu + 3*sigma;

    // Signal: integral of the Gaussian part only over ±3σ
    TF1* gausOnly = new TF1("gausOnly", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
    gausOnly->SetParameters(A, mu, sigma);
    double signal = gausOnly->Integral(x_min, x_max);

    // Noise: integral of background over same range
    double noise = bg * (x_max - x_min);

    // Compute signal-to-noise ratio
    double snr = (noise > 0) ? signal / noise : 0;

    std::cout << "Signal (Gaussian, ±3σ): " << signal << std::endl;
    std::cout << "Noise (Background, ±3σ): " << noise << std::endl;
    std::cout << "Signal-to-Noise Ratio: " << snr << std::endl;


   c8->cd(6);
    // Define the Gaussian + constant background function
    TF1* fitFunc2 = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

    // Set initial parameters:
    // [0] amplitude, [1] mean, [2] sigma, [3] constant background
    fitFunc2->SetParameters(hXDiffECalCDet2->GetMaximum(), 0.02, 0.01, hXDiffECalCDet2->GetMinimum());
    fitFunc2->SetParNames("Amplitude", "Mean", "Sigma", "Background");

    // Optional: set limits on the parameters if needed
    fitFunc2->SetParLimits(1, -0.05, 0.05);  // constrain the mean near 0.02
    fitFunc2->SetParLimits(2, 0.001, 0.05);  // positive sigma

    // Fit the histogram
    hXDiffECalCDet2->Fit(fitFunc2, "R");  // "R" = use function range only

    // Draw the result
    hXDiffECalCDet2->Draw();
    fitFunc2->Draw("same");

    // Extract Gaussian parameters
    double A2 = fitFunc2->GetParameter(0); // Amplitude
    double mu2 = fitFunc2->GetParameter(1); // Mean
    double sigma2 = fitFunc2->GetParameter(2); // Sigma
    double bg2 = fitFunc2->GetParameter(3); // Background level

    // Define integration limits
    double x_min2 = mu2 - 3*sigma2;
    double x_max2 = mu2 + 3*sigma2;

    // Signal: integral of the Gaussian part only over ±3σ
    TF1* gausOnly2 = new TF1("gausOnly2", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
    gausOnly2->SetParameters(A2, mu2, sigma2);
    double signal2 = gausOnly2->Integral(x_min2, x_max2);

    // Noise: integral of background over same range
    double noise2 = bg2 * (x_max2 - x_min2);

    // Compute signal-to-noise ratio
    double snr2 = (noise2 > 0) ? signal2 / noise2 : 0;

    std::cout << "Signal (Gaussian, ±3σ): " << signal2 << std::endl;
    std::cout << "Noise (Background, ±3σ): " << noise2 << std::endl;
    std::cout << "Signal-to-Noise Ratio: " << snr2 << std::endl;

   c8->cd(7);
   hYECalCDet1->Draw();

   c8->cd(8);
   hYECalCDet2->Draw();


   c8->cd(9);
   hXYECal->Draw("colz");


  return c8;
}*/

auto plotXYECalCDet(){

   // 4x3 layout so we can add the two new "min-hit" histograms cleanly
   TCanvas *c8 = new TCanvas("c8", "plotXYECalCDet", 1600,1200);
   c8->Divide(4,3);

   // ---------------- Row 1: XECal vs XCDet (cut-based and min-hit-based) ----------------

   c8->cd(1);
   gPad->SetLogz();
   hXECalCDet1->SetMinimum(10);
   hXECalCDet1->Draw("colz");

   c8->cd(2);
   gPad->SetLogz();
   hXECalCDet2->SetMinimum(10);
   hXECalCDet2->Draw("colz");

   c8->cd(3);
   gPad->SetLogz();
   hXECalCDet1_min->SetMinimum(10);
   hXECalCDet1_min->Draw("colz");

   c8->cd(4);
   gPad->SetLogz();
   hXECalCDet2_min->SetMinimum(10);
   hXECalCDet2_min->Draw("colz");

   // ---------------- Row 2: ECal X/Y + Xdiff fits ----------------

   c8->cd(5);
   hXECal->Draw();

   c8->cd(6);
   hYECal->Draw();

   // ---- Fit XDiff layer 1 (same code as you had, moved to pad 7)
   c8->cd(7);
   {
     TF1* fitFunc = new TF1("fitFunc1", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

     fitFunc->SetParameters(hXDiffECalCDet1->GetMaximum(), 0.02, 0.01, hXDiffECalCDet1->GetMinimum());
     fitFunc->SetParNames("Amplitude", "Mean", "Sigma", "Background");

     fitFunc->SetParLimits(1, -0.05, 0.05);
     fitFunc->SetParLimits(2, 0.001, 0.05);

     hXDiffECalCDet1->Fit(fitFunc, "R");
     hXDiffECalCDet1->Draw();
     fitFunc->Draw("same");

     double A     = fitFunc->GetParameter(0);
     double mu    = fitFunc->GetParameter(1);
     double sigma = fitFunc->GetParameter(2);
     double bg    = fitFunc->GetParameter(3);

     double x_min = mu - 3*sigma;
     double x_max = mu + 3*sigma;

     TF1* gausOnly = new TF1("gausOnly1", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
     gausOnly->SetParameters(A, mu, sigma);

     double signal = gausOnly->Integral(x_min, x_max);
     double noise  = bg * (x_max - x_min);
     double snr    = (noise > 0) ? signal / noise : 0;

     std::cout << "Layer1: Signal (Gaussian, ±3σ): " << signal << std::endl;
     std::cout << "Layer1: Noise  (Background, ±3σ): " << noise  << std::endl;
     std::cout << "Layer1: Signal-to-Noise Ratio: " << snr << std::endl;
   }

   // ---- Fit XDiff layer 2 (same code as you had, moved to pad 8)
   c8->cd(8);
   {
     TF1* fitFunc2 = new TF1("fitFunc2", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

     fitFunc2->SetParameters(hXDiffECalCDet2->GetMaximum(), 0.02, 0.01, hXDiffECalCDet2->GetMinimum());
     fitFunc2->SetParNames("Amplitude", "Mean", "Sigma", "Background");

     fitFunc2->SetParLimits(1, -0.05, 0.05);
     fitFunc2->SetParLimits(2, 0.001, 0.05);

     hXDiffECalCDet2->Fit(fitFunc2, "R");
     hXDiffECalCDet2->Draw();
     fitFunc2->Draw("same");

     double A2     = fitFunc2->GetParameter(0);
     double mu2    = fitFunc2->GetParameter(1);
     double sigma2 = fitFunc2->GetParameter(2);
     double bg2    = fitFunc2->GetParameter(3);

     double x_min2 = mu2 - 3*sigma2;
     double x_max2 = mu2 + 3*sigma2;

     // IMPORTANT FIX: use (x_min2, x_max2) here (your current file uses x_min/x_max by accident)
     TF1* gausOnly2 = new TF1("gausOnly2", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min2, x_max2);
     gausOnly2->SetParameters(A2, mu2, sigma2);

     double signal2 = gausOnly2->Integral(x_min2, x_max2);
     double noise2  = bg2 * (x_max2 - x_min2);
     double snr2    = (noise2 > 0) ? signal2 / noise2 : 0;

     std::cout << "Layer2: Signal (Gaussian, ±3σ): " << signal2 << std::endl;
     std::cout << "Layer2: Noise  (Background, ±3σ): " << noise2  << std::endl;
     std::cout << "Layer2: Signal-to-Noise Ratio: " << snr2 << std::endl;
   }

   // ---------------- Row 3: Y correlations and x residual vs corrected CDet x ----------------

   c8->cd(9);
   hYECalCDet1->Draw();

   c8->cd(10);
   hYECalCDet2->Draw();

   c8->cd(11);
   pXDiffECalCDet1VsXCDet1->SetMarkerStyle(20);
   pXDiffECalCDet1VsXCDet1->SetMarkerSize(0.7);
   pXDiffECalCDet1VsXCDet1->SetMarkerColor(kBlue + 1);
   pXDiffECalCDet1VsXCDet1->SetLineColor(kBlue + 1);
   pXDiffECalCDet1VsXCDet1->Draw("E1");

   c8->cd(12);
   pXDiffECalCDet2VsXCDet2->SetMarkerStyle(20);
   pXDiffECalCDet2VsXCDet2->SetMarkerSize(0.7);
   pXDiffECalCDet2VsXCDet2->SetMarkerColor(kGreen + 2);
   pXDiffECalCDet2VsXCDet2->SetLineColor(kGreen + 2);
   pXDiffECalCDet2VsXCDet2->Draw("E1");

   return c8;
}


//------------Some routines when Ben getting familiar with branches
auto plotDpp(int nbins = 100, double xmin = -0.1, double xmax =  0.1)
{
    TH1D* h = new TH1D("hHeep_dpp", "heep_dpp;#delta p/p;Counts", nbins, xmin, xmax);
    for (double x : vheep_dpp) if (std::isfinite(x)) h->Fill(x);
    auto* c9 = new TCanvas(Form("c_%s","heep_dpp"), "heep_dpp", 900, 650); h->Draw();
    return c9;
}

//auto find hist range
static std::pair<double,double> MinMaxFlat(const std::vector<std::vector<double>>& vv){
  double mn =  std::numeric_limits<double>::infinity();
  double mx = -std::numeric_limits<double>::infinity();
  for (const auto& v : vv){
    for (double x : v){
      if (std::isfinite(x)) {
        if (x < mn) mn = x;
        if (x > mx) mx = x;
      }
    }
  }
  if (!std::isfinite(mn) || !std::isfinite(mx)) { // empty or all non-finite
    mn = 0.0; mx = 1.0;
  }
  if (mn == mx) { // collapse -> pad a bit
    mn -= 0.5; mx += 0.5;
  }
  return {mn, mx};
}

auto plotECalClusX(double xmin = -1.5, double xmax = 1.5, bool log=true)
{
  int nbin = std::ceil((xmax - xmin)/0.0425);
  TH1D* h_ECal_clus_x = new TH1D("h_ECal_clus_x", "hECalClusX;Clus x (m);Counts", nbin,-1.5,1.5);

  //fill hist from vectors
  for (const auto& vec : v_ECal_clus_x) {
    for (double val : vec){
      h_ECal_clus_x->Fill(val);
    }
  }
  TCanvas* c_ECal_clus_x = new TCanvas("c_ECal_clus_x", "ECal Cluster x",800,600);
  if (log==true){
    c_ECal_clus_x->SetLogy();
  }
  h_ECal_clus_x->Draw("HIST");
}
auto plotECalNclus()
{
  int nclus_max = *std::max_element(v_ECal_nclus.begin(), v_ECal_nclus.end());
  int nclus_min = *std::min_element(v_ECal_nclus.begin(),v_ECal_nclus.end());
  TH1D* h_ECal_nclus = new TH1D("h_ECal_nclus", "hECalNclus;Nclus;Counts", nclus_max+1, nclus_min-0.5, nclus_max + 0.5);
  for (auto val : v_ECal_nclus){
    h_ECal_nclus->Fill(val);
  }
  TCanvas* c_ECal_nclus = new TCanvas("c_ECal_nclus", "ECal Cluster Count",800,600);
  h_ECal_nclus->Draw("HIST");
}

auto plotECalClusE(int binlow = 0, int binhigh=12)
{
  TH1D* h_ECal_clus_e = new TH1D("h_ECal_clus_e", "hECalCluse;Clus E (GeV);Counts", (binhigh+binlow)*100,binlow-0.5,binhigh+0.5);

  //fill hist from vectors
  for (const auto& vec : v_ECal_clus_e) {
    for (double val : vec){
      h_ECal_clus_e->Fill(val);
    }
  }
  TCanvas* c_ECal_clus_e = new TCanvas("c_ECal_clus_e", "ECal Cluster e",800,600);
  h_ECal_clus_e->Draw("HIST");
}
auto plotECalClusAdcTime()
{
  double adctime_min = std::numeric_limits<double>::max();
  double adctime_max = std::numeric_limits<double>::lowest();

  for (const auto& subvec : v_ECal_clus_adctime) {
    if (subvec.empty()) continue; // skip events with no clusters

    // find min and max within this event
    double local_min = *std::min_element(subvec.begin(), subvec.end());
    double local_max = *std::max_element(subvec.begin(), subvec.end());

    // update global range
    if (local_min < adctime_min) adctime_min = local_min;
    if (local_max > adctime_max) adctime_max = local_max;
  }

  TH1D* h_ECal_clus_adctime = new TH1D("h_ECal_clus_adctime", "hECalClusAdctime;Clus adctime;Counts", 1000,adctime_min-10,adctime_max+10);

  //fill hist from vectors
  for (const auto& vec : v_ECal_clus_adctime) {
    for (double val : vec){
      h_ECal_clus_adctime->Fill(val);
    }
  }
  TCanvas* c_ECal_clus_adctime = new TCanvas("c_ECal_clus_adctime", "ECal Cluster adctime",800,600);
  h_ECal_clus_adctime->Draw("HIST");
}
