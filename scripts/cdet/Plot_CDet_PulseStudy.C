// Read the run once; keep vectors available for later plotting, like PlotElastic.
// From gocdetscripts, with OUT_DIR set:
//   .L Plot_CDet_PulseStudy.C+
//   Plot_CDet_PulseStudy(6077, 1000);
//   plotPulseStudy("CDet_run6077_pulse_pixel_timing.conf", "baseline");
// Or override display settings and plot only bar 19:
//   plotPulseStudy("CDet_run6077_pulse_pixel_timing.conf", "bar19", 1., -150., 150., 19, 19);
// Saves four PDFs per bar under PusleStudy/<tag>/ without opening plot windows.
// Interactive study/control plots:
//   plotPulseStudyBars("CDet_run6077_pulse_pixel_timing.conf", "first_pass");
// Position/rx overlays (one visible canvas per bar):
//   plotPulseStudyPositions("CDet_run6077_pulse_pixel_timing.conf", "positions", 19, 30);
//   plotPulseStudyPositionCorrelations("CDet_run6077_pulse_pixel_timing.conf", "correlations", 30);
//   plotPulseStudyBars("CDet_run6077_pulse_pixel_timing.conf", "bars_y", true);
//   plotPulseStudyEventMaps("CDet_run6077_pulse_pixel_timing.conf", "event_maps", -1);
// Background plots, including matched rx canvases for multiple/single pulses in bar 30:
//   plotPulseStudyBackground("CDet_run6077_pulse_pixel_timing.conf", "multiplicity", 19, 30);
// Shares CDet_run<run>_pulse_pixel_timing.conf with the pulse-timing macro.
// Cut arguments override that file for this call only; NAN means use the file.
// Loading always applies TotMin; TotMax and ECalMaxE require their analysis.apply_*_max switches.
// Secondary upper switches are pulse_study.apply_tot_max/apply_ecal_energy_max (default off).
// Secondary ECal-energy and ToT lower cuts default on; their switches can disable secondary cuts.
// The secondary pulse_study.apply_tot_cut switch does not disable the loading lower cut.
// Plot routines save full *_cuts.txt reports and compact energy/timing/position footers.
// Optional y correction uses pulse_study.y_propagation_ns_per_m (default 7 ns/m).
// The final numeric y-correction argument now also has units ns/m, not refractive index.
// Add your pulse readers, vectors and hit cuts in the marked sections below.

#include <TChain.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TList.h>
#include <TLine.h>
#include <TROOT.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
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
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "CDetRunDataset.h"

// Display-only target label. Run numbers and output filenames stay numeric.
TString PulseStudyRunLabel(TEnv& config, int runNumber)
{
  TString target = config.GetValue("analysis.target", "");
  target = target.Strip(TString::kBoth);
  TString normalized = target;
  normalized.ToLower();
  normalized.ReplaceAll("_", " ");
  normalized.ReplaceAll("-", " ");
  if (normalized == "lh2") target = "LH2";
  else if (normalized == "cross target") target = "Cross target";
  TString label = TString::Format("Run %d", runNumber);
  if (!target.IsNull()) label += " (" + target + ")";
  return label;
}

// Keep cut annotations consistent with the optional upper bounds.
TString PulseStudyEnergyCutLabel(double minimum, double maximum, bool applyMaximum)
{
  return applyMaximum ? TString::Format("%.12g < E_{ECal} < %.12g GeV", minimum, maximum)
      : TString::Format("E_{ECal} > %.12g GeV (upper cut off)", minimum);
}

TString PulseStudyToTCutLabel(double minimum, double maximum, bool applyMaximum, bool inclusive = true)
{
  const char* relation = inclusive ? "#leq" : "<";
  return applyMaximum ? TString::Format("%.12g %s ToT %s %.12g ns", minimum, relation, relation, maximum)
      : TString::Format("%.12g %s ToT (ns; upper cut off)", minimum, relation);
}

// One element per ECal-selected event with at least one selected CDet pulse.
// All pulse vectors share [event][hit].

//Global constants
const int pixelsPerBar = 16;
const int numBars = 168;
constexpr double kECalZFromTargetM = 6.144;

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

// Histograms from the most recent plotPulseStudy() call, available in ROOT.
TH1D* hStudyBarDt[numBars] = {};
TH1D* hStudyPixelDt[numBars][pixelsPerBar] = {};
TH2D* hStudyBarLeVsToT[numBars] = {};
TH2D* hStudyPixelLeVsToT[numBars][pixelsPerBar] = {};

// Separate summary histograms: index 0 = study bar, 1 = control bar.
TH1D* hPulseSummaryBarDt[2] = {};
TH1D* hPulseSummaryBarLe[2] = {};
TH1D* hPulseSummaryBarLeCorr[2] = {};
TH1D* hPulseSummaryPixelDt[2][pixelsPerBar] = {};
TH1D* hPulseSummaryPixelLe[2][pixelsPerBar] = {};
TH1D* hPulseSummaryPixelLeCorr[2][pixelsPerBar] = {};
TH2D* hPulseSummaryBarLeVsToT[2] = {};
TH2D* hPulseSummaryBarLeCorrVsToT[2] = {};
TH2D* hPulseSummaryPixelLeVsToT[2][pixelsPerBar] = {};
TH2D* hPulseSummaryPixelLeCorrVsToT[2][pixelsPerBar] = {};
TH2D* hPulseBarDtVsRx[2] = {nullptr, nullptr};
// Stored dt BEFORE optional y correction versus projected ECal y.
TH2D* hPulseBarDtVsECalYProjection[2] = {};
TH2D* hPulsePixelDtVsRx[2][pixelsPerBar] = {};
// Fixed bar-30 ECal-time correlations: 0 = raw LE, 1 = corrected LE (optional y correction).
TH2D* hPulseBar30LeVsECalTime[2] = {};
// Separate position-overlay histograms: index 0 = study bar, 1 = control bar.
TH1D* hPulsePositionBarX[2] = {};
TH1D* hPulsePositionBarRx[2] = {};
// Position correlations: [control/detector/L1/L2][ECal native/projected][raw x/corrected x/DB y].
TH2D* hPulsePositionCorrelation[4][2][3] = {};
// Direct coordinate differences: [L1/L2][raw x/corrected x].
TH1D* hPulseLayerXDifference[2][2] = {};
// Fixed bar 30, independent of the configured control bar: channel/aligned x.
TH1D* hPulseBar30XDifference[2] = {};
// Same-event combinations: bar-30 signal pulse in L1 versus every saved L2 pulse.
TH2D* hPulseLayerX2VsX1[2] = {}; // channel/aligned x
TH1D* hPulseLayerX2MinusX1[2] = {};
TH1D* hPulseBar30Multiplicity = nullptr; // all saved bar-30 pulses; one entry per nonempty event
// Full-detector 1D positions: channel CDet x per pulse, native ECal x per selected event.
TH1D* hPulsePositionCDetX = nullptr;
TH1D* hPulsePositionECalX = nullptr;
// Detector maps and one saved event, split by layer.
TH2Poly* hPulseDetectorMap[2] = {};
TH2Poly* hPulseEventMap[2] = {};
TH1D* hPulseEventBarHits = nullptr;

// Multiplicity is counted per event in the background routine's selected bar.
TH1D* hPulseBackgroundRxMultiple = nullptr;
TH1D* hPulseBackgroundRxSingle = nullptr;
TH2D* hPulseBackgroundDtVsRxMultiple = nullptr;
TH2D* hPulseBackgroundDtVsRxSingle = nullptr;
TH1D* hECalAdcTime = nullptr;
TH2D* hPulseBarDtVsTot[2] = {};
TH2D* hPulseBarDtVsTotNoTotCut[2] = {};
TH1D* hPulseFullBarDtAllCuts[2] = {};
TH1D* hPulseFullBarDtNoLeCut[2] = {};
// Remember the actual loading cuts, including main-macro argument overrides.
int pulseStudyRunNumber = 0;
TString pulseStudyLoadingECalCuts;
TString pulseStudyLoadingPulseCuts;
// Denominator sample: every event passing the loading ECal cuts, INCLUDING
// events with no accepted CDet pulse. These do not share the pulse-vector index.
std::vector<double> vPulseStudyDenominatorECalE, vPulseStudyDenominatorECalTime;
bool pulseStudyDenominatorReady = false;

//-------HCal Vectors ---------
std::vector<double> vHCalAdcTime;

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

// Compact plot labels retain only energy, timing and position selections.
// Loading limits are captured when reading data, so edited configs cannot relabel old cuts.
struct PulseStudyPlotCuts {
  bool known = false;
  double energyMin = NAN, energyMax = NAN, timeMin = NAN, timeMax = NAN;
  double leMin = NAN, leMax = NAN, totMin = NAN, totMax = NAN, rx = NAN, ry = NAN;
};
PulseStudyPlotCuts pulseStudyLoadedPlotCuts;
TString pulseStudyPlotRoutine;
TString pulseStudyPlotEnergy, pulseStudyPlotECalTiming;
TString pulseStudyPlotLoadingTiming, pulseStudyPlotLoadingPosition;
TString pulseStudyPlotSecondaryTiming, pulseStudyPlotSecondaryPosition;
TString pulseStudyPlotDtWindow;
double pulseStudyPlotPeakNSigma = NAN;
bool pulseStudyPlotPositionDt = false;

TString PulseStudyWindow(const char* quantity, double low, double high, const char* unit, bool inclusive = true)
{
  if (!std::isfinite(low)) return TString::Format("%s: unknown", quantity);
  return std::isfinite(high) ? TString::Format("%s: %s%.3g, %.3g%s %s", quantity, inclusive ? "[" : "(", low, high, inclusive ? "]" : ")", unit)
      : TString::Format("%s %s %.3g %s", quantity, inclusive ? "#geq" : ">", low, unit);
}

void PulseStudyPreparePlotCuts(const char* routine, TEnv& config)
{
  pulseStudyPlotRoutine = routine;
  const bool main = pulseStudyPlotRoutine == "Plot_CDet_PulseStudy";
  const bool offsets = pulseStudyPlotRoutine == "extractBarOffsets";
  const bool background = pulseStudyPlotRoutine == "plotPulseStudyBackground";
  const bool maps = pulseStudyPlotRoutine == "plotPulseStudyEventMaps";
  const auto& load = pulseStudyLoadedPlotCuts;
  TString energy = load.known ? PulseStudyWindow("E", load.energyMin, load.energyMax, "GeV", false) : "loading unknown";
  pulseStudyPlotLoadingTiming = load.known ? PulseStudyWindow("t_{E}", load.timeMin, load.timeMax, "ns", false) + "; "
      + PulseStudyWindow("LE", load.leMin, load.leMax, "ns", false) + "; "
      + PulseStudyWindow("ToT", load.totMin, load.totMax, "ns", false) : "loading unknown";
  pulseStudyPlotLoadingPosition = load.known ? TString::Format("|rx| #leq %.3g, |ry| #leq %.3g m", load.rx, load.ry) : "loading unknown";
  if (!main && !offsets && (background || config.GetValue("pulse_study.apply_ecal_energy_cut", 1))) {
    TString secondary = PulseStudyWindow("E", config.GetValue("pulse_study.ecal_energy_min", NAN),
        config.GetValue("pulse_study.apply_ecal_energy_max", 0) ? config.GetValue("pulse_study.ecal_energy_max", NAN) : NAN, "GeV", background);
    if (secondary != energy) energy += "; secondary " + secondary;
  }
  pulseStudyPlotEnergy = "Energy: " + energy;
  TString timing;
  auto append = [&](TString text) { if (!timing.IsNull()) timing += "; "; timing += text; };
  if (!main && (offsets || background || config.GetValue("pulse_study.apply_ecal_time_cut", 0)))
    append(PulseStudyWindow("t_{E}", config.GetValue("pulse_study.ecal_time_min", NAN), config.GetValue("pulse_study.ecal_time_max", NAN), "ns", offsets || background));
  if (!main && !offsets) {
    append(PulseStudyWindow(background ? "LE_{corr}" : "LE", config.GetValue("analysis.le_min", .02), config.GetValue("analysis.le_max", 60.), "ns"));
    if (background || config.GetValue("pulse_study.apply_tot_cut", 1))
      append(PulseStudyWindow("ToT", config.GetValue("analysis.tot_min", 4.),
          config.GetValue("pulse_study.apply_tot_max", 0) ? config.GetValue("analysis.tot_max", 30.) : NAN, "ns"));
  }
  pulseStudyPlotECalTiming = (!main && (offsets || background || config.GetValue("pulse_study.apply_ecal_time_cut", 0)))
      ? PulseStudyWindow("t_{E}", config.GetValue("pulse_study.ecal_time_min", NAN), config.GetValue("pulse_study.ecal_time_max", NAN), "ns", offsets || background) : TString("no secondary ECal time window");
  pulseStudyPlotSecondaryTiming = timing;
  TString position = "none";
  if (!main && !background && !maps) {
    const double radius = config.GetValue("good_pulse.single_cut_radius", -1.);
    if (!offsets && radius > 0.)
      position = TString::Format("((rx-%.3g)/%.3g)^{2}+((dt-%.3g)/%.3g)^{2} #leq %.3g^{2} (m,ns)",
          config.GetValue("good_pulse.single_residual_center_m", 0.), config.GetValue("good_pulse.single_residual_scale_m", .04),
          config.GetValue("good_pulse.single_timing_center_ns", -26.), config.GetValue("good_pulse.single_timing_scale_ns", 5.), radius);
    else if (pulseStudyPlotRoutine != "plotPulseStudy") position = TString::Format("|rx| #leq %.3g m", config.GetValue("analysis.x_difference_max", .08));
  }
  pulseStudyPlotSecondaryPosition = position;
  pulseStudyPlotDtWindow = PulseStudyWindow("dt", config.GetValue("pulse_study.position_dt_min_ns", config.GetValue("pulse_study.fit_signal_min", -36.)),
      config.GetValue("pulse_study.position_dt_max_ns", config.GetValue("pulse_study.fit_signal_max", -20.)), "ns");
  pulseStudyPlotPositionDt = config.GetValue("pulse_study.position_apply_dt_cut", 0) != 0;
}

// A small footer outside the plotting area, shared by single and divided canvases.
// Existing numbered pads stay intact (offset fits revisit them after the first save).
void PulseStudySaveCanvas(TCanvas* canvas, const TString& filename)
{
  canvas->cd();
  const TString name = canvas->GetName();
  const bool position = pulseStudyPlotRoutine == "plotPulseStudyPositionCorrelations";
  const bool layers = position && name.BeginsWith("cPulseLayers");
  const bool multiplicity = position && name == "cPulseBar30Multiplicity";
  const bool eventMap = pulseStudyPlotRoutine == "plotPulseStudyEventMaps";
  const bool singleEvent = eventMap && (name.BeginsWith("cPulseEventMap") || name == "cPulseEventBarHits");
  const bool ecalOnly = name == "cPulseBackgroundECalTime";
  TString secondaryTiming = multiplicity || ecalOnly ? pulseStudyPlotECalTiming : pulseStudyPlotSecondaryTiming;
  if (layers) secondaryTiming += "; bar 30 " + pulseStudyPlotDtWindow + "; L2: loading only";
  else if ((position && !multiplicity && pulseStudyPlotPositionDt) || (eventMap && (singleEvent || pulseStudyPlotPositionDt)))
    secondaryTiming += "; " + pulseStudyPlotDtWindow;
  else if (position || eventMap) secondaryTiming += "; no secondary dt window";
  if (name == "cPulsePeakYieldVsBar") secondaryTiming += TString::Format("; peak: #mu #pm %.3g#sigma", pulseStudyPlotPeakNSigma);
  TString secondaryPosition = multiplicity || ecalOnly ? "none" : pulseStudyPlotSecondaryPosition;
  if (layers) secondaryPosition = "bar 30: " + secondaryPosition + "; L2: none";
  std::vector<TString> lines = {
      pulseStudyPlotEnergy,
      "Timing (load): " + pulseStudyPlotLoadingTiming,
      "Timing (secondary): " + (secondaryTiming.IsNull() ? TString("none") : secondaryTiming),
      "Position (load): " + pulseStudyPlotLoadingPosition + "; secondary: " + secondaryPosition
  };
  if (name.Contains("CutComparison"))
    lines.push_back("Before-LE / before-ToT panels omit that secondary window; loading cuts remain.");
  TPaveText* label = dynamic_cast<TPaveText*>(canvas->GetListOfPrimitives()->FindObject("pulseCompactCuts"));
  if (!label) {
    const double footer = .15;
    bool divided = false;
    for (auto* object : *canvas->GetListOfPrimitives()) {
      if (auto* pad = dynamic_cast<TPad*>(object)) {
        const double x = pad->GetXlowNDC(), y = pad->GetYlowNDC();
        const double w = pad->GetWNDC(), h = pad->GetHNDC();
        pad->SetPad(x, footer + (1.-footer)*y, x+w, footer + (1.-footer)*(y+h));
        pad->Modified(); divided = true;
      }
    }
    if (!divided) canvas->SetBottomMargin(canvas->GetBottomMargin() + footer);
    label = new TPaveText(.015, .005, .985, .14, "NDC");
    label->SetName("pulseCompactCuts"); label->SetBit(kCanDelete);
    label->SetFillColor(kWhite); label->SetBorderSize(0);
    label->SetTextFont(42); label->SetTextAlign(12);
    label->Draw();
  }
  label->Clear();
  label->SetTextSize(0); // Fit the concise lines to the footer width.
  for (const auto& line : lines) label->AddText(line);
  canvas->Modified(); canvas->Update();
  canvas->SaveAs(filename);
}

// Reports use the resolved loading metadata, not a reread of old loading cuts.
// Each routine has its own filename, so several routines may share an output tag.
bool PulseStudyWriteCutReport(const TString& path, const char* routine,
                              const char* configFile, TEnv& config,
                              const std::vector<TString>& selection,
                              const std::vector<TString>& details = {})
{
  PulseStudyPreparePlotCuts(routine, config);
  std::ofstream report(path.Data());
  if (!report) {
    std::cerr << "[CDet pulse study] Cannot write cut report: " << path << '\n';
    return false;
  }
  auto line = [&](TString text) {
    text.ReplaceAll("#leq", "<="); text.ReplaceAll("#geq", ">=");
    text.ReplaceAll("#mu", "mu"); text.ReplaceAll("#sigma", "sigma");
    text.ReplaceAll("#pm", "+/-");
    text.ReplaceAll("E_{ECal}", "E_ECal"); text.ReplaceAll("t_{ECal}", "t_ECal");
    report << text << '\n';
  };
  report << "Routine: " << routine << "\nConfiguration: " << configFile << '\n';
  line(PulseStudyRunLabel(config, pulseStudyRunNumber > 0 ? pulseStudyRunNumber : config.GetValue("analysis.run_number", 0)));
  report << "\nLOADING CUTS (recorded when the saved sample was loaded)\n";
  line(pulseStudyLoadingECalCuts.IsNull() ? "Loading ECal cuts: metadata unavailable; rerun the main macro to record them." : pulseStudyLoadingECalCuts);
  line(pulseStudyLoadingPulseCuts.IsNull() ? "Loading pulse cuts: metadata unavailable; rerun the main macro to record them." : pulseStudyLoadingPulseCuts);
  report << "Loading requires pulse.ecal_eligible and pulse.calib_valid > 0.5; spatial bounds are recorded above.\n"
         << "Saved events have at least one pulse passing loading cuts.\n"
         << "Secondary routines cannot recover pulses/events rejected during loading.\n"
         << "\nSELECTION IN THIS ROUTINE\n";
  for (const auto& text : selection) line(text);
  report << "\nPLOT / FIT SETTINGS AND EXCEPTIONS\n";
  for (const auto& text : details) line(text);
  report << "Display bounds are not cuts; out-of-range values remain in flow bins.\n";
  report.flush();
  if (!report) {
    std::cerr << "[CDet pulse study] Failed writing cut report: " << path << '\n';
    return false;
  }
  std::cout << "[CDet pulse study] Cut report: " << path << '\n';
  return true;
}

// Standard secondary routines share these config cuts (no cut-argument overrides).
// Background and offset extraction have their own explicit reports below.
std::vector<TString> PulseStudySecondaryCuts(TEnv& config, bool rectangularRx = true)
{
  std::vector<TString> cuts;
  const bool energy = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool time = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const bool tot = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  cuts.push_back("Finite ECal energy/time; integer pulse ID in [0,2687]; finite stored dt.");
  cuts.push_back(energy ? "Secondary: " + PulseStudyEnergyCutLabel(config.GetValue("pulse_study.ecal_energy_min", NAN), config.GetValue("pulse_study.ecal_energy_max", NAN), config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0) : "Secondary ECal energy window: off");
  cuts.push_back(time ? TString::Format("Secondary ECal time: %.12g < t_ECal < %.12g ns", config.GetValue("pulse_study.ecal_time_min", NAN), config.GetValue("pulse_study.ecal_time_max", NAN)) : TString("Secondary ECal time window: off"));
  cuts.push_back(TString::Format("Finite raw LE: %.12g <= LE <= %.12g ns", config.GetValue("analysis.le_min", 0.02), config.GetValue("analysis.le_max", 60.0)));
  cuts.push_back("Finite ToT > 0 ns required.");
  cuts.push_back(tot ? "Secondary: " + PulseStudyToTCutLabel(config.GetValue("analysis.tot_min", 4.0), config.GetValue("analysis.tot_max", 30.0), config.GetValue("pulse_study.apply_tot_max", 0) != 0) : "Secondary ToT window: off");
  const double radius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  if (radius > 0.0)
    cuts.push_back(TString::Format("Finite rx; ellipse: ((rx-%.12g)/%.12g)^2 + ((dt-%.12g)/%.12g)^2 <= %.12g^2; rx in m, dt in ns. No additional rectangular rx cut.", config.GetValue("good_pulse.single_residual_center_m", 0.0), config.GetValue("good_pulse.single_residual_scale_m", 0.040), config.GetValue("good_pulse.single_timing_center_ns", -26.0), config.GetValue("good_pulse.single_timing_scale_ns", 5.0), radius));
  else cuts.push_back(rectangularRx ? TString::Format("Ellipse off; finite rx and abs(rx) <= %.12g m", config.GetValue("analysis.x_difference_max", 0.08)) : TString("Ellipse off; no secondary rectangular rx cut in this routine."));
  cuts.push_back("Selections use original stored dt and raw LE, before any optional y correction.");
  return cuts;
}

TString PulseStudyYCorrectionDescription(bool enabled, double delay, const char* scope)
{
  return enabled ? TString::Format("Y correction ON (%s): propagation delay %.12g ns/m; d=abs(y_DB +/- 0.25 m - y_proj), + for LHS / - for RHS; dt += delay*(d-0.25 m), corrected LE -= same. y_proj uses the pulse z. Saved vectors unchanged.", scope, delay)
      : TString("Y correction OFF.");
}

void ClearPulseStudyVectors()
{
  vPulseStudyDenominatorECalE.clear();
  vPulseStudyDenominatorECalTime.clear();
  pulseStudyDenominatorReady = false;
  vHCalAdcTime.clear();
  delete gROOT->GetListOfCanvases()->FindObject("cPulseStudyECalAdcTime");
  delete hECalAdcTime;
  hECalAdcTime = nullptr;
  pulseStudyRunNumber = 0;
  pulseStudyLoadedPlotCuts = PulseStudyPlotCuts{};
  pulseStudyLoadingECalCuts = "";
  pulseStudyLoadingPulseCuts = "";
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
  if (ApplyToTCut == -1) ApplyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  if (std::isnan(ECalMinE)) ECalMinE = config.GetValue("analysis.ecal_energy_min", NAN);
  if (std::isnan(ECalMaxE)) ECalMaxE = config.GetValue("analysis.ecal_energy_max", NAN);
  const bool applyECalEnergyMax = config.GetValue("analysis.apply_ecal_energy_max", 0) != 0;
  const bool applyToTMax = config.GetValue("analysis.apply_tot_max", 0) != 0;
  if (std::isnan(ECalMinT)) ECalMinT = config.GetValue("analysis.ecal_time_min", NAN);
  if (std::isnan(ECalMaxT)) ECalMaxT = config.GetValue("analysis.ecal_time_max", NAN);
  if (std::isnan(EllipseRadius)) EllipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const double RxMax = config.GetValue("analysis.x_difference_max", 0.08);
  const double RyMax = 0.36; // m; retain the previous y window while allowing wider rx.
  const bool applyToTCut = ApplyToTCut != 0;
  const bool useEllipse = EllipseRadius > 0.0;
  // Reject events outside the configured ECal ADC-time window before the hit loop.
  // The per-pulse replay eligibility flag still supplies its cluster/position checks.
  // Optional tighter ADC-time cuts belong to the secondary plotting routines.
  // Display ranges and fit bounds in the config are NOT selection cuts.

  if (RunNumber <= 0 || nevents < -1 || minSeg < -1 || maxSeg < -1 ||
      !std::isfinite(ECalMinE) || (applyECalEnergyMax && (!std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)) ||
      !std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT) {
    std::cerr << "[CDet pulse study] Invalid run, event/segment limit or ECal energy/time bounds.\n";
    return;
  }
  if ((ApplyToTCut != 0 && ApplyToTCut != 1) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (!std::isfinite(TotMin) || TotMin < 0.0 ||
       (applyToTMax && (!std::isfinite(TotMax) || TotMin >= TotMax)))) {
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

  if (!std::isfinite(RxMax) || RxMax <= 0.0) {
    std::cerr << "[CDet pulse study] analysis.x_difference_max must be finite and positive (m).\n";
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

  const double ECalDisplayMin = config.GetValue("display.ecal_time_min", -15.0);
  const double ECalDisplayMax = config.GetValue("display.ecal_time_max", 15.0);
  const double ECalBinWidth = config.GetValue("pulse_study.ecal_time_bin_width_ns", 1.0);
  if (!std::isfinite(ECalDisplayMin) || !std::isfinite(ECalDisplayMax) ||
      ECalDisplayMin >= ECalDisplayMax ||
      !std::isfinite(ECalBinWidth) || ECalBinWidth <= 0.0 ||
      !std::isfinite((ECalDisplayMax - ECalDisplayMin) / ECalBinWidth) ||
      std::ceil((ECalDisplayMax - ECalDisplayMin) / ECalBinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse study] Invalid ECal time display bounds or bin width.\n";
    return;
  }
  pulseStudyLoadedPlotCuts = {true, ECalMinE, applyECalEnergyMax ? ECalMaxE : NAN, ECalMinT, ECalMaxT, LeMin, LeMax, TotMin, applyToTMax ? TotMax : NAN, RxMax, RyMax};
  pulseStudyRunNumber = RunNumber;
  pulseStudyLoadingECalCuts = "Loading: " + PulseStudyEnergyCutLabel(ECalMinE, ECalMaxE, applyECalEnergyMax)
      + TString::Format(", %.12g < t_{ECal} < %.12g ns", ECalMinT, ECalMaxT);
  pulseStudyLoadingPulseCuts = "Loading: " + PulseStudyToTCutLabel(TotMin, TotMax, applyToTMax, false)
      + TString::Format("; %.12g < LE < %.12g ns; finite |rx| <= %.12g m, |ry| <= %.12g m; ecal_eligible and calib_valid > 0.5; replay spatial_pass not required", LeMin, LeMax, RxMax, RyMax);
  hECalAdcTime = new TH1D("hECalAdcTime", TString::Format("%s: ECal ADC time, saved events;ECal ADC time (ns);Events", PulseStudyRunLabel(config, RunNumber).Data()), int(std::ceil((ECalDisplayMax - ECalDisplayMin) / ECalBinWidth)), ECalDisplayMin, ECalDisplayMax);
  hECalAdcTime->SetDirectory(nullptr);

  // 3. SET YOUR TTREEREADERS HERE, before the event loop.
  TTreeReader reader(&T);
  //-------ECal readers--------
  TTreeReaderValue<double> ECalE(reader, "earm.ecal.e");
  TTreeReaderValue<double> ECalAdcTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<double> ECalX(reader, "earm.ecal.x");
  TTreeReaderValue<double> ECalY(reader, "earm.ecal.y");
  TTreeReaderArray<double> ECalXProjection(reader, "earm.cdet.pulse.ecal_x_proj");
  TTreeReaderArray<double> ECalYProjection(reader, "earm.cdet.pulse.ecal_y_proj");
  //------HCal readers--------
  TTreeReaderValue<double> HCalAdcTime(reader, "sbs.hcal.adctime");
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

  vECalE.reserve(limit);
  vPulseStudyDenominatorECalE.reserve(limit);
  vPulseStudyDenominatorECalTime.reserve(limit);
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
            << "  " << PulseStudyEnergyCutLabel(ECalMinE, ECalMaxE, applyECalEnergyMax) << "\n"
            << "  ECal ADC time: " << ECalMinT << " < t < " << ECalMaxT << " ns (before hit loop).\n"
            << "  LE: " << LeMin << " < LE < " << LeMax << " ns\n"
            << "  Loading: " << PulseStudyToTCutLabel(TotMin, TotMax, applyToTMax, false) << "\n"
            << "  Ellipse radius: " << EllipseRadius << " (-1 = off).\n"
            << "  Loading requires ECal energy/time, replay eligibility/calibration flags, explicit rx/ry cuts, ToT, and LE.\n";

  // 4. ONE EVENT LOOP. Save values here; do histogram filling later.
  Long64_t EventCounter = 0;
  while (EventCounter < limit && reader.Next()) {
    ++EventCounter;
    if (EventCounter % 1000 == 0)
      std::cout << "[CDet pulse study] " << EventCounter << '/' << limit << " events\n";

    // ADD YOUR EVENT-LEVEL CUT BOOLEANS HERE.
    //ECal bools
    bool good_ECal_energy = std::isfinite(*ECalE) && *ECalE > ECalMinE && (!applyECalEnergyMax || *ECalE < ECalMaxE);
    bool good_ECal_time = std::isfinite(*ECalAdcTime) && *ECalAdcTime > ECalMinT && *ECalAdcTime < ECalMaxT;
    bool good_ECal_event = good_ECal_energy && good_ECal_time;
    if (!good_ECal_event) continue;

    // Count the ECal opportunity before requiring any CDet hit. Keep energy/time
    // so secondary cuts can be applied to this denominator as well as the pulses.
    vPulseStudyDenominatorECalE.push_back(*ECalE);
    vPulseStudyDenominatorECalTime.push_back(*ECalAdcTime);

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
      // CDet Hit Cuts
      bool good_ecal_eligible = PulseECalEligible[hit] > 0.5;
      const double rx = PulseXResiduals[hit];
      const double ry = PulseYResiduals[hit];
      bool good_spatial_pass = std::isfinite(rx) && std::isfinite(ry) &&
          std::abs(rx) <= RxMax && std::abs(ry) <= RyMax;
      bool good_calibration = PulseGoodCalib[hit] > 0.5;
      bool good_pulse_tot = std::isfinite(PulseToT[hit]) && PulseToT[hit] > TotMin && (!applyToTMax || PulseToT[hit] < TotMax);
      const double leNs = 0.01 * PulseLe[hit];
      bool good_le = std::isfinite(leNs) && leNs > LeMin && leNs < LeMax;
      const bool good_hit = good_ecal_eligible && good_spatial_pass && good_calibration && good_pulse_tot && good_le;
      if (!good_hit) continue;

      // Copy all fields for the SAME passing pulse; convert raw LE to ns.
      thisEventECalXProjection.push_back(ECalXProjection[hit]);
      thisEventECalYProjection.push_back(ECalYProjection[hit]);
      thisEventPulseID.push_back(PulseID[hit]);
      thisEventPulseLe.push_back(leNs); // Same calibrated units used for the loading LE cut.
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
    vHCalAdcTime.push_back(*HCalAdcTime);
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
    hECalAdcTime->Fill(*ECalAdcTime); // Scalar: once per saved event, not per pulse.
  }
  if (EventCounter != limit) {
    std::cerr << "[CDet pulse study] Reader stopped after " << EventCounter << '/' << limit
              << " events; status " << reader.GetEntryStatus()
              << ". Saved vectors contain only ECal-selected events processed so far.\n";
    return;
  }
  std::cout << "[CDet pulse study] Read " << EventCounter << " events; saved "
            << vECalE.size() << " ECal-selected events.\n";
  pulseStudyDenominatorReady = true;
  std::cout << "  ECal denominator (including zero-hit events): "
            << vPulseStudyDenominatorECalE.size() << '\n';
  TCanvas* cECalTime = new TCanvas("cPulseStudyECalAdcTime", TString::Format("%s: ECal ADC time", PulseStudyRunLabel(config, RunNumber).Data()), 1000, 750);
  cECalTime->SetTopMargin(0.10);
  hECalAdcTime->SetStats(false);
  hECalAdcTime->Draw("HIST");
  cECalTime->Modified();
  cECalTime->Update();
  if (!gSystem->AccessPathName("PusleStudy") || gSystem->mkdir("PusleStudy", true) == 0) {
    if (!PulseStudyWriteCutReport(TString::Format("PusleStudy/run%d_loading_cuts.txt", RunNumber), "Plot_CDet_PulseStudy", configFile, config,
        {"Explicit loading rx/ry bounds above replace replay spatial_pass; no ellipse in the main loop.",
         "ToT lower bound is always applied, regardless of ApplyToTCut. ECal energy/time and raw LE windows use strict inequalities.",
         "The ECal histogram counts events with at least one accepted CDet pulse."},
        {TString::Format("Segments %d..%d; requested events %lld; processed events %lld; saved events %zu", minSeg, maxSeg, nevents, EventCounter, vPulseID.size()),
         TString::Format("Input directory: %s", directory)})) return;
    PulseStudySaveCanvas(cECalTime, TString::Format("PusleStudy/run%d_ecal_adc_time.pdf", RunNumber));
  }
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
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const bool applyECalEnergyMax = config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0;
  const bool applyToTMax = config.GetValue("pulse_study.apply_tot_max", 0) != 0;
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if ((applyECalEnergyCut && (!std::isfinite(ECalMinE) || (applyECalEnergyMax && (!std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)))) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || TotMin < 0 || (applyToTMax && (!std::isfinite(TotMax) || TotMin >= TotMax))))) {
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
  if (!PulseStudyWriteCutReport(outputDirectory + "/plotPulseStudy_cuts.txt", "plotPulseStudy", configFile, config,
      PulseStudySecondaryCuts(config, false), {
        TString::Format("Bars %d..%d; dt display [%.12g,%.12g] ns, requested bin width %.12g ns", firstBar, lastBar, DtMin, DtMax, BinWidth),
        TString::Format("Raw LE display [%.12g,%.12g], width %.12g ns; ToT display [%.12g,%.12g], width %.12g ns", LeDisplayMin, LeDisplayMax, LeBinWidth, TotDisplayMin, TotDisplayMax, TotBinWidth),
        "Y correction OFF; dt uses pulse.ecal_residual; raw LE versus ToT."})) return;

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
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && (!applyECalEnergyMax || vECalE[event] < ECalMaxE)));
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
          (!applyToTCut || (totNs >= TotMin && (!applyToTMax || totNs <= TotMax)));
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
    PulseStudySaveCanvas(cBar, TString::Format("%s/bar%03d_dt.pdf", outputDirectory.Data(), bar));
    delete cBar;

    TCanvas* cPixels = new TCanvas(TString::Format("cStudyPixelDt_bar%d", bar), TString::Format("CDet bar %d: pixel dt", bar), 1400, 1000);
    cPixels->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixels->cd(pixel + 1);
      hStudyPixelDt[bar][pixel]->Draw("HIST");
    }
    cPixels->Modified();
    cPixels->Update();
    PulseStudySaveCanvas(cPixels, TString::Format("%s/bar%03d_pixels_dt.pdf", outputDirectory.Data(), bar));
    delete cPixels;

    TCanvas* cLeToTBar = new TCanvas(TString::Format("cStudyBarLeVsToT_bar%d", bar), TString::Format("CDet bar %d: summed raw LE versus ToT", bar), 900, 650);
    hStudyBarLeVsToT[bar]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hStudyBarLeVsToT[bar]->Draw("COLZ");
    PulseStudySaveCanvas(cLeToTBar, TString::Format("%s/bar%03d_le_vs_tot.pdf", outputDirectory.Data(), bar));
    delete cLeToTBar;

    TCanvas* cLeToTPixels = new TCanvas(TString::Format("cStudyPixelLeVsToT_bar%d", bar), TString::Format("CDet bar %d: pixel raw LE versus ToT", bar), 1400, 1000);
    cLeToTPixels->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cLeToTPixels->cd(pixel + 1);
      hStudyPixelLeVsToT[bar][pixel]->SetStats(false);
      gPad->SetRightMargin(0.15);
      hStudyPixelLeVsToT[bar][pixel]->Draw("COLZ");
    }
    PulseStudySaveCanvas(cLeToTPixels, TString::Format("%s/bar%03d_pixels_le_vs_tot.pdf", outputDirectory.Data(), bar));
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

// Additive correction to dt = ECal time - corrected LE, referenced to the bar center.
// Global +y points toward the LHS readout. Logical IDs repeat LHS/RHS in each layer.
// Caller supplies a valid pixel ID, finite positions (m), and positive propagation delay (ns/m).
double PulseStudyYTimingCorrection(int pixelID, double yDB, double yProj, double propagationNsPerMeter)
{
  const double barLength = 0.50; // m
  const bool isLHS = pixelID % 1344 < 672;
  const double yReadout = yDB + (isLHS ? +0.5 : -0.5) * barLength;
  const double distance = std::abs(yReadout - yProj);
  return propagationNsPerMeter * (distance - barLength / 2.0);
}

// Interactive Gaussian + linear-background signal widths; run the main analysis first.
// getBarSigma("CDet_run6077_pulse_pixel_timing.conf", "sigma_study");
// getBarSigma("this.conf", "bar30_y", true); // optional bar-30 y correction
// firstBar/lastBar limit only the bar graphs; the detector spectrum includes all bars.
// Each bar yield counts raw pulses within its fitted mean +/- PeakNSigma*sigma,
// divided by all ECal-selected events (including zero-hit events); it is not Hz.
void getBarSigma(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                 const char* outputTag = nullptr, int firstBar = -1, int lastBar = -1,
                 double BinWidth = NAN, double DtMin = NAN, double DtMax = NAN,
                 double PeakNSigma = NAN, bool ApplyYCorrection = false,
                 double YPropagationNsPerMeter = NAN)
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
  if (firstBar == -1) firstBar = config.GetValue("pulse_study.first_bar", 0);
  if (lastBar == -1) lastBar = config.GetValue("pulse_study.last_bar", numBars - 1);
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const bool applyECalEnergyMax = config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0;
  const bool applyToTMax = config.GetValue("pulse_study.apply_tot_max", 0) != 0;
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double FitSignalMin = config.GetValue("pulse_study.fit_signal_min", -36.0);
  const double FitSignalMax = config.GetValue("pulse_study.fit_signal_max", -20.0);
  const double FitMinEntries = config.GetValue("pulse_study.min_fit_entries", 100.0);
  const double FitSigma = config.GetValue("pulse_study.fit_sigma", 5.0);
  if (std::isnan(PeakNSigma)) PeakNSigma = config.GetValue("pulse_study.peak_nsigma", 2.0);
  // Effective light propagation delay includes indirect photon paths in the bar.
  // Use endpoint distance for either side; this diagnostic still applies only to bar 30.
  const int yCorrectionBar = 30;
  if (std::isnan(YPropagationNsPerMeter))
    YPropagationNsPerMeter = config.GetValue("pulse_study.y_propagation_ns_per_m", 7.0);
  const double propagationNsPerMeter = YPropagationNsPerMeter;
  if (ApplyYCorrection && (!std::isfinite(propagationNsPerMeter) || YPropagationNsPerMeter <= 0.0)) {
    std::cerr << "[CDet bar sigma] Y correction requires a finite, positive propagation delay (ns/m).\n";
    return;
  }
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if (firstBar < 0 || lastBar >= numBars || firstBar > lastBar ||
      controlBar < 0 || controlBar >= numBars ||
      !std::isfinite(BinWidth) || BinWidth <= 0.0 ||
      !std::isfinite(DtMin) || !std::isfinite(DtMax) || DtMin >= DtMax ||
      !std::isfinite((DtMax - DtMin) / BinWidth) ||
      std::ceil((DtMax - DtMin) / BinWidth) > std::numeric_limits<int>::max() - 2 ||
      !std::isfinite(FitSignalMin) || !std::isfinite(FitSignalMax) ||
      FitSignalMin >= FitSignalMax || FitSignalMin < DtMin || FitSignalMax > DtMax ||
      !std::isfinite(FitMinEntries) || FitMinEntries <= 0.0 ||
      !std::isfinite(PeakNSigma) || PeakNSigma <= 0.0 ||
      !std::isfinite(FitSigma) || FitSigma <= 0.0) {
    std::cerr << "[CDet bar sigma] Invalid bar, histogram, or fit settings. Fit windows must lie inside the dt display range.\n";
    return;
  }
  if (!std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      !std::isfinite(xRes) || xRes <= 0.0 ||
      (applyToTCut && (!std::isfinite(TotMin) || TotMin < 0.0 || (applyToTMax && (!std::isfinite(TotMax) || TotMin >= TotMax)))) ||
      (applyECalEnergyCut && (!std::isfinite(ECalMinE) || (applyECalEnergyMax && (!std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)))) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) || ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet bar sigma] Invalid secondary cuts.\n";
    return;
  }
  if (vPulseID.empty()) {
    std::cerr << "[CDet bar sigma] No saved pulse events; run the main macro first.\n";
    return;
  }
  // Do not partially correct a sample if the required saved geometry is missing.
  // With the switch off, these vectors are not required or accessed.
  if (ApplyYCorrection) {
    for (size_t event = 0; event < vPulseID.size(); ++event) {
      for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
        const double id = vPulseID[event][hit];
        if (!std::isfinite(id) || id < yCorrectionBar * pixelsPerBar ||
            id >= (yCorrectionBar + 1) * pixelsPerBar || std::floor(id) != id) continue;
        if (event >= vPulseY.size() || event >= vECalYProjection.size() ||
            hit >= vPulseY[event].size() || hit >= vECalYProjection[event].size() ||
            !std::isfinite(vPulseY[event][hit]) || !std::isfinite(vECalYProjection[event][hit])) {
          std::cerr << "[CDet bar sigma] Missing/nonfinite bar-30 DB y or projected ECal y; rerun the main analysis.\n";
          return;
        }
      }
    }
  }
  const int nDtBins = int(std::ceil((DtMax - DtMin) / BinWidth));
  const double fitSpan = FitSignalMax - FitSignalMin;
  const double sigmaMin = 0.5 * (DtMax - DtMin) / nDtBins;
  const double sigmaMax = fitSpan;
  if (sigmaMin >= sigmaMax || nDtBins < 6) {
    std::cerr << "[CDet bar sigma] Binning is too coarse for a five-parameter fit.\n";
    return;
  }
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
  if (!PulseStudyWriteCutReport(outputDirectory + "/getBarSigma_cuts.txt", "getBarSigma", configFile, config,
      PulseStudySecondaryCuts(config), {
        TString::Format("Control bar %d; graph bars %d..%d; detector histogram includes all bars", controlBar, firstBar, lastBar),
        TString::Format("dt display [%.12g,%.12g] ns, requested bin width %.12g ns", DtMin, DtMax, BinWidth),
        TString::Format("Gaussian + linear fit [%.12g,%.12g] ns; minimum fit-window entries %.12g; initial sigma %.12g ns", FitSignalMin, FitSignalMax, FitMinEntries, FitSigma),
        TString::Format("Yield numerator: raw pulses in each bar's fitted mu +/- %.12g*sigma; includes background; not Hz", PeakNSigma),
        "Yield denominator: all loading ECal-selected events, including zero-hit events, after this routine's secondary ECal cuts; no secondary CDet cuts on denominator.",
        "Mean sigma: inverse-variance weighting of valid per-bar fit sigmas; independent statistical errors only.",
        PulseStudyYCorrectionDescription(ApplyYCorrection, YPropagationNsPerMeter, "bar 30 only, including its detector-spectrum contribution")
      })) return;


  const int RunNumber = pulseStudyRunNumber > 0 ? pulseStudyRunNumber : config.GetValue("analysis.run_number", 0);

  pulseStudyPlotPeakNSigma = PeakNSigma;

  // Clear the previous call before replacing the histograms/functions it displayed.
  const char* canvasNames[] = {"cPulseSigmaControl", "cPulseSigmaDetector", "cPulseSigmaVsBar", "cPulsePeakYieldVsBar"};
  for (const char* name : canvasNames) delete gROOT->GetListOfCanvases()->FindObject(name);
  static std::vector<TH1D*> hBarDt(numBars + 1, nullptr);
  static std::vector<TF1*> fBarDt(numBars + 1, nullptr);
  static std::vector<TF1*> fBarBackgroundDt(numBars + 1, nullptr);
  static TGraphErrors* gSigmaLayer1 = nullptr;
  static TGraphErrors* gSigmaLayer2 = nullptr;
  static TGraphErrors* gSigmaUnavailable = nullptr;
  delete gSigmaLayer1;
  delete gSigmaLayer2;
  delete gSigmaUnavailable;
  gSigmaLayer1 = new TGraphErrors();
  gSigmaLayer2 = new TGraphErrors();
  gSigmaUnavailable = new TGraphErrors();
  gSigmaLayer1->SetName("gPulseSigmaLayer1");
  gSigmaLayer2->SetName("gPulseSigmaLayer2");
  gSigmaUnavailable->SetName("gPulseSigmaUnavailable");

  static TGraphErrors *gYieldLayer1 = nullptr, *gYieldLayer2 = nullptr, *gYieldUnavailable = nullptr;
  delete gYieldLayer1; delete gYieldLayer2; delete gYieldUnavailable;
  gYieldLayer1 = new TGraphErrors();
  gYieldLayer2 = new TGraphErrors();
  gYieldUnavailable = new TGraphErrors();
  gYieldLayer1->SetName("gPulsePeakYieldLayer1");
  gYieldLayer2->SetName("gPulsePeakYieldLayer2");
  gYieldUnavailable->SetName("gPulsePeakYieldUnavailable");

  // The last slot is the full detector; all earlier slots are individual bars.
  // Fit parameters: signal amplitude/mean/sigma, then background intercept/slope.
  for (int bar = 0; bar <= numBars; ++bar) {
    delete hBarDt[bar];
    delete fBarDt[bar];
    delete fBarBackgroundDt[bar];
    TString title = bar == numBars
        ? TString::Format("%s: full-detector dt", PulseStudyRunLabel(config, RunNumber).Data())
        : TString::Format("%s: bar %d dt", PulseStudyRunLabel(config, RunNumber).Data(), bar);
    if (ApplyYCorrection && (bar == yCorrectionBar || bar == numBars))
      title += " (bar 30 y correction on)";
    hBarDt[bar] = new TH1D(TString::Format("hPulseSigmaDt_%d", bar),
        title + (ApplyYCorrection && (bar == yCorrectionBar || bar == numBars)
            ? ";dt (ns), bar 30 y-corrected;Pulses" : ";dt = pulse.ecal_residual (ns);Pulses"), nDtBins, DtMin, DtMax);
    hBarDt[bar]->SetDirectory(nullptr);
    hBarDt[bar]->SetStats(false);
    fBarDt[bar] = new TF1(TString::Format("fPulseSigmaDt_%d", bar),
        "gaus(0) + pol1(3)", FitSignalMin, FitSignalMax);
    fBarDt[bar]->AddToGlobalList(false);
    fBarDt[bar]->SetParNames("Signal amplitude", "Signal mean", "Signal sigma", "Linear background intercept", "Linear background slope");
    fBarDt[bar]->SetNpx(1000);
    fBarDt[bar]->SetLineColor(kRed + 1);
  }

  // Apply identical secondary ECal cuts to the full denominator sample.
  Long64_t ecalEventCount = 0;
  const bool denominatorAvailable = pulseStudyDenominatorReady &&
      vPulseStudyDenominatorECalE.size() == vPulseStudyDenominatorECalTime.size();
  if (denominatorAvailable) {
    for (size_t event = 0; event < vPulseStudyDenominatorECalE.size(); ++event) {
      const double e = vPulseStudyDenominatorECalE[event], t = vPulseStudyDenominatorECalTime[event];
      if (std::isfinite(e) && std::isfinite(t) &&
          (!applyECalEnergyCut || (e > ECalMinE && (!applyECalEnergyMax || e < ECalMaxE))) &&
          (!applyECalTimeCut || (t > ECalMinT && t < ECalMaxT))) ++ecalEventCount;
    }
  }
  // Retain exact per-bar dt values so peak counts are not rounded to histogram bins.
  std::vector<std::vector<double>> barDt(numBars);
  // One pass through the aligned saved vectors. dt is already in ns.
  Long64_t filledPulses = 0;
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    const bool good_ECal_energy = std::isfinite(vECalE[event]) &&
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && (!applyECalEnergyMax || vECalE[event] < ECalMaxE)));
    const bool good_ECal_time = std::isfinite(vECalAdcTime[event]) &&
        (!applyECalTimeCut || (vECalAdcTime[event] > ECalMinT && vECalAdcTime[event] < ECalMaxT));
    if (!good_ECal_energy || !good_ECal_time) continue;
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit];
      const double dt = vPulseDt[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars * pixelsPerBar ||
          std::floor(id) != id || !std::isfinite(dt)) continue;
      const int bar = int(id) / pixelsPerBar;
      // All saved timing vectors used here are in ns.
      const double le = vPulseLe[event][hit];
      const double totNs = vPulseToT[event][hit];
      const double rx = vPulseXResiduals[event][hit];
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(totNs) && totNs > 0.0 &&
          (!applyToTCut || (totNs >= TotMin && (!applyToTMax || totNs <= TotMax)));
      const bool good_x_res = std::isfinite(rx) && std::abs(rx) <= xRes;

      bool good_ellipse = true;
      if (useEllipse) {
        // Use the same branch dt for selection and fitting.
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
      // Selection above always uses the original branch dt and raw LE.
      // Remove propagation delay relative to the bar center: dt_y = dt + propagationNsPerMeter*(d-L/2).
      // The saved projection already uses each pulse's CDet z: y_ECal*z_CDet/z_ECal.
      const double plotDt = ApplyYCorrection && bar == yCorrectionBar
          ? dt + PulseStudyYTimingCorrection(int(id), vPulseY[event][hit],
              vECalYProjection[event][hit], YPropagationNsPerMeter) : dt;
      hBarDt[bar]->Fill(plotDt);
      hBarDt[numBars]->Fill(plotDt);
      barDt[bar].push_back(plotDt);
      ++filledPulses;
    }//hit loop for event
  }//event loop

  // Fit each bar and the detector sum. Never interpret a failed fit as zero width.
  std::vector<double> sigma(numBars + 1, NAN), sigmaError(numBars + 1, NAN);
  std::vector<TString> fitReason(numBars + 1);
  for (int bar = 0; bar <= numBars; ++bar) {
    TH1D* hist = hBarDt[bar];
    TF1* fit = fBarDt[bar];
    double fitEntries = 0.0;
    double height = 0.0;
    int occupiedBins = 0;
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
      if (hist->GetBinCenter(bin) < FitSignalMin || hist->GetBinCenter(bin) > FitSignalMax) continue;
      fitEntries += hist->GetBinContent(bin);
      height = std::max(height, hist->GetBinContent(bin));
      if (hist->GetBinContent(bin) > 0.0) ++occupiedBins;
    }
    if (fitEntries < FitMinEntries || occupiedBins <= 5) {
      fitReason[bar] = TString::Format("Insufficient statistics: %.0f fit-window entries, %d occupied bins",
                                      fitEntries, occupiedBins);
    } else {
      const double peakGuess = std::max(FitSignalMin, std::min(FitSignalMax, ellipseTimeCenter));
      fit->SetParameters(height, peakGuess,
          std::max(sigmaMin, std::min(sigmaMax, FitSigma)), height / 3.0, 0.0);
      fit->SetParLimits(0, 0.0, 10.0 * height);
      fit->SetParLimits(1, FitSignalMin, FitSignalMax);
      fit->SetParLimits(2, sigmaMin, sigmaMax);
      // The linear intercept and slope are free; neither is a timing parameter.
      // S returns diagnostics; N keeps ownership here and prevents automatic overlays.
      TFitResultPtr result = hist->Fit(fit, "QRS0N");
      const double width = fit->GetParameter(2);
      const double error = fit->GetParError(2);
      const double mean = fit->GetParameter(1);
      const bool atLimit = width <= sigmaMin * 1.001 || width >= sigmaMax * 0.999 ||
          mean <= FitSignalMin + 0.001 * (FitSignalMax - FitSignalMin) ||
          mean >= FitSignalMax - 0.001 * (FitSignalMax - FitSignalMin);
      if (int(result) != 0 || !result.Get() || !result->IsValid() || result->CovMatrixStatus() < 2 ||
          !std::isfinite(width) || !std::isfinite(error) || error <= 0.0 || atLimit ||
          !std::isfinite(mean) || fit->GetParameter(0) <= 0.0) {
        fitReason[bar] = "No usable signal fit (convergence, covariance, limits, or signal amplitude)";
      } else {
        sigma[bar] = width;
        sigmaError[bar] = error;
      }
    }
    const TString label = bar == numBars ? TString("Detector") : TString::Format("Bar %d", bar);
    std::cout << "[CDet bar sigma] " << label << ": ";
    if (std::isfinite(sigma[bar])) std::cout << sigma[bar] << " +/- " << sigmaError[bar] << " ns\n";
    else std::cout << fitReason[bar] << '\n';
  }

  // Use one common ECal opportunity count, including events with zero CDet hits.
  // These are raw peak-window yields: no background subtraction or conversion to Hz.
  std::vector<double> peakYield(numBars, NAN);
  std::vector<TString> yieldLabels(numBars);
  for (int bar = 0; bar < numBars; ++bar) {
    if (!denominatorAvailable) {
      yieldLabels[bar] = "Peak yield unavailable: rerun main analysis to save all ECal events";
    } else if (ecalEventCount == 0) {
      yieldLabels[bar] = "Peak yield unavailable: no ECal events pass these cuts";
    } else if (!std::isfinite(sigma[bar])) {
      yieldLabels[bar] = TString::Format("Peak yield unavailable: no usable bar fit; %lld ECal events", ecalEventCount);
    } else {
      const double mu = fBarDt[bar]->GetParameter(1);
      const double peakMin = mu - PeakNSigma * sigma[bar];
      const double peakMax = mu + PeakNSigma * sigma[bar];
      Long64_t peakPulses = 0;
      for (double dt : barDt[bar]) if (dt >= peakMin && dt <= peakMax) ++peakPulses;
      peakYield[bar] = double(peakPulses) / ecalEventCount;
      yieldLabels[bar] = TString::Format("Peak [%.2f, %.2f] ns: %lld / %lld = %.5g pulses/ECal event (raw)",
          peakMin, peakMax, peakPulses, ecalEventCount, peakYield[bar]);
      std::cout << "[CDet bar yield] Bar " << bar << ": " << yieldLabels[bar] << '\n';
    }
  }
  const TString yieldLabel = yieldLabels[controlBar];
  std::cout << "[CDet bar sigma] Control bar " << controlBar << ": " << yieldLabel
            << "; window = fitted mean +/- " << PeakNSigma << " sigma.\n";

  double yieldGraphMax = 0.0;
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    if (!std::isfinite(peakYield[bar])) {
      gYieldUnavailable->SetPoint(gYieldUnavailable->GetN(), bar, 0.0);
      continue;
    }
    TGraphErrors* graph = bar < numBars / 2 ? gYieldLayer1 : gYieldLayer2;
    graph->SetPoint(graph->GetN(), bar, peakYield[bar]);
    yieldGraphMax = std::max(yieldGraphMax, peakYield[bar]);
  }
  if (yieldGraphMax <= 0.0) yieldGraphMax = 1.0;
  gYieldLayer1->SetMarkerColor(kRed); gYieldLayer1->SetLineColor(kRed);
  gYieldLayer2->SetMarkerColor(kBlue); gYieldLayer2->SetLineColor(kBlue);
  gYieldUnavailable->SetMarkerColor(kBlack); gYieldUnavailable->SetLineColor(kBlack);
  for (TGraphErrors* graph : {gYieldLayer1, gYieldLayer2, gYieldUnavailable}) {
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.85);
  }

  double graphMax = 1.0;
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    if (!std::isfinite(sigma[bar])) {
      gSigmaUnavailable->SetPoint(gSigmaUnavailable->GetN(), bar, 0.0);
      continue;
    }
    TGraphErrors* graph = bar < numBars / 2 ? gSigmaLayer1 : gSigmaLayer2;
    const int point = graph->GetN();
    graph->SetPoint(point, bar, sigma[bar]);
    graph->SetPointError(point, 0.0, sigmaError[bar]);
    graphMax = std::max(graphMax, sigma[bar] + sigmaError[bar]);
  }
  gSigmaLayer1->SetMarkerColor(kRed);
  gSigmaLayer1->SetLineColor(kRed);
  gSigmaLayer2->SetMarkerColor(kBlue);
  gSigmaLayer2->SetLineColor(kBlue);
  gSigmaUnavailable->SetMarkerColor(kBlack);
  gSigmaUnavailable->SetLineColor(kBlack);
  for (TGraphErrors* graph : {gSigmaLayer1, gSigmaLayer2, gSigmaUnavailable}) {
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.85);
  }

  // Inverse-variance average of usable individual-bar widths in the plotted range.
  // The detector-sum fit and failed-fit zero placeholders are not measurements here.
  int weightedBars = 0;
  double errorScale = std::numeric_limits<double>::infinity();
  for (int bar = firstBar; bar <= lastBar; ++bar)
    if (std::isfinite(sigma[bar]) && std::isfinite(sigmaError[bar]) && sigmaError[bar] > 0.0)
      errorScale = std::min(errorScale, sigmaError[bar]);
  double weightSum = 0.0, weightedSigmaSum = 0.0;
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    if (!std::isfinite(sigma[bar]) || !std::isfinite(sigmaError[bar]) || sigmaError[bar] <= 0.0) continue;
    // A common rescaling preserves 1/error^2 weights without overflowing them.
    const double relativeError = errorScale / sigmaError[bar];
    const double weight = relativeError * relativeError;
    weightSum += weight;
    weightedSigmaSum += weight * sigma[bar];
    ++weightedBars;
  }
  const double weightedSigma = weightedBars > 0 ? weightedSigmaSum / weightSum : NAN;
  const double weightedSigmaError = weightedBars > 0 ? errorScale / std::sqrt(weightSum) : NAN;
  const TString weightedSigmaLabel = weightedBars > 0
      ? TString::Format("Weighted mean #sigma = %.4f #pm %.4f ns (%d valid bars, range %d-%d)",
          weightedSigma, weightedSigmaError, weightedBars, firstBar, lastBar)
      : TString::Format("Weighted mean #sigma unavailable: no usable fits in bars %d-%d", firstBar, lastBar);
  std::cout << "[CDet bar sigma] " << weightedSigmaLabel
            << "; inverse-variance weights; uncertainty assumes independent bar-fit errors.\n";

  // Separate interactive canvases, each with width/height = the golden ratio.
  const int canvasWidth = 1200;
  const int canvasHeight = int(std::lround(canvasWidth / ((1.0 + std::sqrt(5.0)) / 2.0)));
  for (int plot = 0; plot < 4; ++plot) {
    TCanvas* canvas = new TCanvas(canvasNames[plot], canvasNames[plot], canvasWidth, canvasHeight);
    canvas->SetCanvasSize(canvasWidth, canvasHeight);
    canvas->SetTopMargin(0.26);
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.04);
    canvas->SetBottomMargin(0.12);
    if (plot < 2) {
      const int index = plot == 0 ? controlBar : numBars;
      hBarDt[index]->SetMaximum(std::max(1.0, hBarDt[index]->GetMaximum() * 1.30));
      hBarDt[index]->Draw("HIST");
      if (std::isfinite(sigma[index])) fBarDt[index]->Draw("SAME");
      TPaveText* fitLabel = new TPaveText(0.10, 0.77, 0.96, 0.92, "NDC");
      fitLabel->SetBit(kCanDelete);
      fitLabel->SetFillColor(kWhite);
      fitLabel->SetBorderSize(0);
      fitLabel->SetTextSize(plot == 0 ? 0.020 : 0.025);
      if (std::isfinite(sigma[index])) {
        const int ndf = fBarDt[index]->GetNDF();
        const TString reducedChi2 = ndf > 0
            ? TString::Format("#chi^{2}/ndf = %.3f (ndf = %d)", fBarDt[index]->GetChisquare() / ndf, ndf)
            : TString("#chi^{2}/ndf unavailable (ndf #leq 0)");
        fitLabel->AddText(TString::Format("#sigma_{signal} = %.3f #pm %.3f ns; %s", sigma[index], sigmaError[index], reducedChi2.Data()));
        fitLabel->AddText(TString::Format("Gaussian + linear background; fit [%.1f, %.1f] ns", FitSignalMin, FitSignalMax));
      } else {
        fitLabel->AddText("Signal #sigma unavailable");
        fitLabel->AddText(fitReason[index].BeginsWith("Insufficient") ? "Insufficient fit-window statistics" : "No usable Gaussian + linear fit; see console");
      }
      if (plot == 0) fitLabel->AddText(yieldLabel);
      fitLabel->Draw();
    } else if (plot == 2) {
      TH1F* frame = canvas->DrawFrame(firstBar - 0.5, -0.08 * graphMax,
          lastBar + 0.5, 1.25 * graphMax,
          TString::Format("%s: signal Gaussian width by bar;Bar number (global, zero-based);#sigma_{signal} (ns)", PulseStudyRunLabel(config, RunNumber).Data()));
      frame->SetStats(false);
      gSigmaLayer1->Draw("P SAME");
      gSigmaLayer2->Draw("P SAME");
      gSigmaUnavailable->Draw("P SAME");
      TLine* averageLine = nullptr;
      if (weightedBars > 0) {
        averageLine = new TLine(firstBar - 0.5, weightedSigma, lastBar + 0.5, weightedSigma);
        averageLine->SetBit(kCanDelete);
        averageLine->SetLineColor(kBlack);
        averageLine->SetLineStyle(2);
        averageLine->SetLineWidth(2);
        averageLine->Draw();
      }
      TLegend* legend = new TLegend(0.10, 0.83, 0.96, 0.92);
      legend->SetBit(kCanDelete);
      legend->SetBorderSize(0);
      legend->SetFillColor(kWhite);
      legend->SetTextSize(0.022);
      legend->SetNColumns(2);
      legend->AddEntry(gSigmaLayer1, "Layer 1 (bars 0-83)", "pe");
      legend->AddEntry(gSigmaLayer2, "Layer 2 (bars 84-167)", "pe");
      legend->AddEntry(gSigmaUnavailable, "Unusable fit (marked at 0)", "p");
      if (averageLine) legend->AddEntry(averageLine, "Inverse-variance weighted mean", "l");
      legend->Draw();
    } else {
      TH1F* frame = canvas->DrawFrame(firstBar - 0.5, -0.08 * yieldGraphMax,
          lastBar + 0.5, 1.25 * yieldGraphMax,
          TString::Format("%s: peak pulse yield by bar;Bar number (global, zero-based);Peak pulses / ECal-selected event",
              PulseStudyRunLabel(config, RunNumber).Data()));
      frame->SetStats(false);
      gYieldLayer1->Draw("P SAME");
      gYieldLayer2->Draw("P SAME");
      gYieldUnavailable->Draw("P SAME");
      TLegend* legend = new TLegend(0.10, 0.83, 0.96, 0.92);
      legend->SetBit(kCanDelete); legend->SetBorderSize(0); legend->SetFillColor(kWhite);
      legend->SetTextSize(0.022); legend->SetNColumns(2);
      legend->AddEntry(gYieldLayer1, "Layer 1 (bars 0-83)", "p");
      legend->AddEntry(gYieldLayer2, "Layer 2 (bars 84-167)", "p");
      legend->AddEntry(gYieldUnavailable, "Unavailable fit/denominator (marked at 0, not measured zero)", "p");
      legend->Draw();
    }

    if (plot == 2) {
      TLatex* average = new TLatex(0.12, 0.79, weightedSigmaLabel);
      average->SetNDC(); average->SetTextSize(0.025); average->SetBit(kCanDelete);
      average->Draw();
    }
    canvas->Modified();
    canvas->Update();
    const TString filename = plot == 0 ? TString::Format("bar%03d_control_dt_sigma.pdf", controlBar)
        : TString(plot == 1 ? "detector_dt_sigma.pdf" : plot == 2 ? "sigma_vs_bar.pdf" : "peak_yield_vs_bar.pdf");
    PulseStudySaveCanvas(canvas, outputDirectory + "/" + filename);
  }
  std::cout << "[CDet bar sigma] Selected " << filledPulses << " pulses; saved four PDFs under " << outputDirectory << '\n';
  if (ApplyYCorrection) std::cout << "[CDet bar sigma] Bar 30 y correction on: " << propagationNsPerMeter
      << " ns/m; d=|y_readout-y_proj|, L=0.50 m, center-referenced; original selections and saved vectors preserved.\n";
  if (useEllipse) std::cout << "[CDet bar sigma] The dt-dependent ellipse can narrow the measured spectra.\n";

}//end getBarSigma macro

// Convenient opt-in call; existing numeric firstBar/lastBar calls keep their meaning.
void getBarSigma(const char* configFile, const char* outputTag, bool ApplyYCorrection,
                 double YPropagationNsPerMeter = NAN)
{
  getBarSigma(configFile, outputTag, -1, -1, NAN, NAN, NAN, NAN,
              ApplyYCorrection, YPropagationNsPerMeter);
}


// Interactive plots for the study/control bars: one spectrum or 4x4 grid per canvas.
// The optional ellipse uses branch dt; otherwise the configured x-residual cut applies.
// dt uses vPulseDt. Raw and corrected LE are saved in ns; LE cuts always use raw LE.
// dt starts from pulse.ecal_residual; optional y correction shifts plotted dt and corrected LE.
// Leaves canvases open and does not change ROOT's batch setting. Use an interactive session.
void plotPulseStudyBars(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                        const char* outputTag = nullptr, int studyBar = -1, int controlBar = -1,
                        double BinWidth = NAN, double DtMin = NAN, double DtMax = NAN,
                        double LeBinWidth = NAN, double LeDisplayMin = NAN, double LeDisplayMax = NAN,
                        double TotBinWidth = NAN, double TotDisplayMin = NAN, double TotDisplayMax = NAN,
                        bool ApplyYCorrection = false, double YPropagationNsPerMeter = NAN)
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
  const double ECalDisplayMin = config.GetValue("display.ecal_time_min", -15.0);
  const double ECalDisplayMax = config.GetValue("display.ecal_time_max", 15.0);
  const double ECalBinWidth = config.GetValue("pulse_study.ecal_time_bin_width_ns", 1.0);
  if (!std::isfinite(ECalDisplayMin) || !std::isfinite(ECalDisplayMax) || ECalDisplayMin >= ECalDisplayMax ||
      !std::isfinite(ECalBinWidth) || ECalBinWidth <= 0.0 ||
      !std::isfinite((ECalDisplayMax-ECalDisplayMin)/ECalBinWidth) ||
      std::ceil((ECalDisplayMax-ECalDisplayMin)/ECalBinWidth) > std::numeric_limits<int>::max()-2) {
    std::cerr << "[CDet pulse bar summary] Invalid ECal-time display settings.\n";
    return;
  }
  const int nECalBins = int(std::ceil((ECalDisplayMax-ECalDisplayMin)/ECalBinWidth));
  const double YProjectionMin = config.GetValue("pulse_study.position_min_m", -1.6);
  const double YProjectionMax = config.GetValue("pulse_study.position_max_m", 1.6);
  const double YProjectionBinWidth = config.GetValue("pulse_study.position_bin_width_m", 0.005);
  if (!std::isfinite(YProjectionMin) || !std::isfinite(YProjectionMax) || YProjectionMin >= YProjectionMax ||
      !std::isfinite(YProjectionBinWidth) || YProjectionBinWidth <= 0.0 ||
      !std::isfinite((YProjectionMax-YProjectionMin)/YProjectionBinWidth) ||
      std::ceil((YProjectionMax-YProjectionMin)/YProjectionBinWidth) > std::numeric_limits<int>::max()-2) {
    std::cerr << "[CDet pulse bar summary] Invalid projected-y display settings.\n";
    return;
  }
  const int nYProjectionBins = int(std::ceil((YProjectionMax-YProjectionMin)/YProjectionBinWidth));
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const bool applyECalEnergyMax = config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0;
  const bool applyToTMax = config.GetValue("pulse_study.apply_tot_max", 0) != 0;
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if ((applyECalEnergyCut && (!std::isfinite(ECalMinE) || (applyECalEnergyMax && (!std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)))) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || TotMin < 0 || (applyToTMax && (!std::isfinite(TotMax) || TotMin >= TotMax))))) {
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
  // Projected y is now needed even with the optional correction switched off.
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    if (event >= vECalYProjection.size() || vECalYProjection[event].size() != vPulseID[event].size()) {
      std::cerr << "[CDet pulse bar summary] Missing/misaligned projected ECal y; rerun main analysis.\n";
      return;
    }
  }
  if ((double(nYProjectionBins)+2.0)*(std::ceil((DtMax-DtMin)/BinWidth)+2.0) > std::numeric_limits<int>::max()) {
    std::cerr << "[CDet pulse bar summary] Too many dt-versus-projected-y histogram cells.\n";
    return;
  }
  if (std::isnan(YPropagationNsPerMeter))
    YPropagationNsPerMeter = config.GetValue("pulse_study.y_propagation_ns_per_m", 7.0);
  if (ApplyYCorrection) {
    if (!std::isfinite(YPropagationNsPerMeter) || YPropagationNsPerMeter <= 0.0) {
      std::cerr << "[CDet pulse bar summary] Invalid y propagation delay (ns/m).\n";
      return;
    }
    // Validate before replacing plots; include the fixed bar-30 ECal-time diagnostic.
    for (size_t event = 0; event < vPulseID.size(); ++event) {
      for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
        const double id = vPulseID[event][hit];
        if (!std::isfinite(id) || id < 0 || id >= numBars*pixelsPerBar || std::floor(id) != id) continue;
        const int bar = int(id) / pixelsPerBar;
        if (bar != studyBar && bar != controlBar && bar != 30) continue;
        if (event >= vPulseY.size() || event >= vECalYProjection.size() ||
            hit >= vPulseY[event].size() || hit >= vECalYProjection[event].size() ||
            !std::isfinite(vPulseY[event][hit]) || !std::isfinite(vECalYProjection[event][hit])) {
          std::cerr << "[CDet pulse bar summary] Missing/nonfinite DB y or projected ECal y; rerun main analysis.\n";
          return;
        }
      }
    }
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
  if (!PulseStudyWriteCutReport(outputDirectory + "/plotPulseStudyBars_cuts.txt", "plotPulseStudyBars", configFile, config,
      PulseStudySecondaryCuts(config), {
        TString::Format("Study bar %d; control bar %d; LE versus ECal time plots always use bar 30", studyBar, controlBar),
        "dt-versus-projected-ECal-y diagnostics use original vPulseDt (pulse.ecal_residual), BEFORE optional y correction, with the full usual selection. Finite projected y required only for these new plots.",
        TString::Format("Projected-y display [%.12g,%.12g] m, requested bin width %.12g m; x axis zooms to occupied bins plus two bins of padding", YProjectionMin, YProjectionMax, YProjectionBinWidth),
        "Corrected-LE plots additionally require finite corrected LE; raw LE and dt plots retain pulses with nonfinite corrected LE.",
        "LE comparison: red omits ONLY secondary raw-LE window; blue applies all secondary cuts. Loading LE cut remains in both.",
        "ToT comparison: before omits ONLY secondary ToT window; finite ToT > 0 and secondary raw-LE cut remain. After applies all secondary cuts. Loading ToT cut remains in both.",
        TString::Format("dt display [%.12g,%.12g] ns, requested bin width %.12g ns", DtMin, DtMax, BinWidth),
        TString::Format("LE display [%.12g,%.12g], width %.12g ns; ToT display [%.12g,%.12g], width %.12g ns", LeDisplayMin, LeDisplayMax, LeBinWidth, TotDisplayMin, TotDisplayMax, TotBinWidth),
        PulseStudyYCorrectionDescription(ApplyYCorrection, YPropagationNsPerMeter, "both selected bars and fixed bar 30")
      })) return;


  const char* leECalCanvasNames[] = {"cPulseBar30RawLeVsECalTime", "cPulseBar30LeCorrVsECalTime"};
  for (int kind = 0; kind < 2; ++kind) {
    delete gROOT->GetListOfCanvases()->FindObject(leECalCanvasNames[kind]);
    delete hPulseBar30LeVsECalTime[kind];
    hPulseBar30LeVsECalTime[kind] = nullptr;
  }

  // Replace only this routine's previous canvases, including when the bars change.
  delete gROOT->GetListOfCanvases()->FindObject("cControlBar");
  for (int selected = 0; selected < 2; ++selected) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sCutComparison", roles[selected]));
    delete hPulseBarDtVsTot[selected];
    delete hPulseBarDtVsTotNoTotCut[selected];
    delete hPulseFullBarDtAllCuts[selected];
    delete hPulseFullBarDtNoLeCut[selected];
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%s", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarDt", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarLe", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarLeCorr", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelDt", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelLe", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelLeCorr", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarLeVsToT", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarLeCorrVsToT", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelLeVsToT", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelLeCorrVsToT", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sBarDtVsRx", roles[selected]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sPixelDtVsRx", roles[selected]));
    delete hPulseSummaryBarDt[selected];
    delete hPulseSummaryBarLe[selected];
    delete hPulseSummaryBarLeCorr[selected];
    delete hPulseSummaryBarLeVsToT[selected];
    delete hPulseSummaryBarLeCorrVsToT[selected];
    delete hPulseBarDtVsRx[selected];
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseSummary%sDtVsECalYProjection", roles[selected]));
    delete hPulseBarDtVsECalYProjection[selected];
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      delete hPulseSummaryPixelDt[selected][pixel];
      delete hPulseSummaryPixelLe[selected][pixel];
      delete hPulseSummaryPixelLeCorr[selected][pixel];
      delete hPulseSummaryPixelLeVsToT[selected][pixel];
      delete hPulseSummaryPixelLeCorrVsToT[selected][pixel];
      delete hPulsePixelDtVsRx[selected][pixel];
    }
  }

  const int RunNumber = pulseStudyRunNumber > 0 ? pulseStudyRunNumber : config.GetValue("analysis.run_number", 0);
  for (int kind = 0; kind < 2; ++kind) {
    hPulseBar30LeVsECalTime[kind] = new TH2D(TString::Format("hPulseBar30%sVsECalTime", kind == 0 ? "RawLe" : "LeCorr"),
        TString::Format("%s, bar 30: %s LE versus ECal ADC time;ECal ADC time (ns);%s CDet LE (ns);Pulses",
            PulseStudyRunLabel(config, RunNumber).Data(), kind == 0 ? "raw" : "corrected", kind == 0 ? "Raw" : "Corrected"),
        nECalBins, ECalDisplayMin, ECalDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
    hPulseBar30LeVsECalTime[kind]->SetDirectory(nullptr);
    hPulseBar30LeVsECalTime[kind]->SetStats(false);
  }
  // Define all histograms before filling. Display limits are not additional cuts.
  for (int selected = 0; selected < 2; ++selected) {
    const int bar = bars[selected];
    hPulseBarDtVsTot[selected] = new TH2D(TString::Format("hPulseSummary%sBarDtVsTot", roles[selected]), TString::Format("%s, bar %d: dt versus ToT, secondary ToT cut %s;ToT (ns);Stored pulse.ecal_residual (ns);Pulses", PulseStudyRunLabel(config, RunNumber).Data(), bar, applyToTCut ? "applied" : "off"), nTotBins, TotDisplayMin, TotDisplayMax, nDtBins, DtMin, DtMax);
    hPulseBarDtVsTotNoTotCut[selected] = new TH2D(TString::Format("hPulseSummary%sBarDtVsTotNoTotCut", roles[selected]), TString::Format("%s, bar %d: dt versus ToT, before secondary ToT cut;ToT (ns);Stored pulse.ecal_residual (ns);Pulses", PulseStudyRunLabel(config, RunNumber).Data(), bar), nTotBins, TotDisplayMin, TotDisplayMax, nDtBins, DtMin, DtMax);
    hPulseFullBarDtAllCuts[selected] = new TH1D(TString::Format("hPulseSummary%sFullBarDtAllCuts", roles[selected]), TString::Format("%s, bar %d: dt with LE cut;Stored pulse.ecal_residual (ns);Pulses", PulseStudyRunLabel(config, RunNumber).Data(), bar), nDtBins, DtMin, DtMax);
    hPulseFullBarDtNoLeCut[selected] = new TH1D(TString::Format("hPulseSummary%sFullBarDtNoLeCut", roles[selected]), TString::Format("%s, bar %d: dt without LE cut;Stored pulse.ecal_residual (ns);Pulses", PulseStudyRunLabel(config, RunNumber).Data(), bar), nDtBins, DtMin, DtMax);
    hPulseBarDtVsTot[selected]->SetDirectory(nullptr);
    hPulseBarDtVsTotNoTotCut[selected]->SetDirectory(nullptr);
    hPulseFullBarDtAllCuts[selected]->SetDirectory(nullptr);
    hPulseFullBarDtNoLeCut[selected]->SetDirectory(nullptr);
    hPulseSummaryBarDt[selected] = new TH1D(TString::Format("hPulseSummary%sBarDt", roles[selected]), TString::Format("%s bar %d: branch dt, all pixels;Stored pulse.ecal_residual (ns);Pulses", roles[selected], bar), nDtBins, DtMin, DtMax);
    hPulseSummaryBarLe[selected] = new TH1D(TString::Format("hPulseSummary%sBarLe", roles[selected]), TString::Format("%s bar %d: raw LE, all pixels;Raw CDet LE (ns);Pulses", roles[selected], bar), nLeBins, LeDisplayMin, LeDisplayMax);
    hPulseSummaryBarLeCorr[selected] = new TH1D(TString::Format("hPulseSummary%sBarLeCorr", roles[selected]), TString::Format("%s bar %d: corrected LE, all pixels;Corrected CDet LE (ns);Pulses", roles[selected], bar), nLeBins, LeDisplayMin, LeDisplayMax);
    hPulseSummaryBarDt[selected]->SetDirectory(nullptr);
    hPulseSummaryBarLe[selected]->SetDirectory(nullptr);
    hPulseSummaryBarLeCorr[selected]->SetDirectory(nullptr);
    hPulseSummaryBarLeVsToT[selected] = new TH2D(TString::Format("hPulseSummary%sBarLeVsToT", roles[selected]), TString::Format("%s bar %d: raw LE versus ToT, all pixels;ToT (ns);Raw CDet LE (ns);Pulses", roles[selected], bar), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
    hPulseSummaryBarLeCorrVsToT[selected] = new TH2D(TString::Format("hPulseSummary%sBarLeCorrVsToT", roles[selected]), TString::Format("%s bar %d: corrected LE versus ToT, all pixels;ToT (ns);Corrected CDet LE (ns);Pulses", roles[selected], bar), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
    hPulseSummaryBarLeVsToT[selected]->SetDirectory(nullptr);
    hPulseSummaryBarLeCorrVsToT[selected]->SetDirectory(nullptr);
    hPulseBarDtVsRx[selected] = new TH2D(TString::Format("hPulseSummary%sBarDtVsRx", roles[selected]), TString::Format("%s bar %d: branch dt versus rx, all pixels;rx (m);Stored pulse.ecal_residual (ns);Pulses", roles[selected], bar), 100, -xRes, xRes, nDtBins, DtMin, DtMax);
    hPulseBarDtVsRx[selected]->SetDirectory(nullptr);
    hPulseBarDtVsECalYProjection[selected] = new TH2D(TString::Format("hPulseSummary%sDtVsECalYProjection", roles[selected]),
        TString::Format("%s bar %d: dt before y correction versus ECal y projection;ECal y projected to pixel z (m);dt before y correction (ns);Pulses", roles[selected], bar),
        nYProjectionBins, YProjectionMin, YProjectionMax, nDtBins, DtMin, DtMax);
    hPulseBarDtVsECalYProjection[selected]->SetDirectory(nullptr);
    hPulseBarDtVsECalYProjection[selected]->SetStats(false);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      const int id = bar * pixelsPerBar + pixel;
      hPulseSummaryPixelDt[selected][pixel] = new TH1D(TString::Format("hPulseSummary%sPixel%dDt", roles[selected], pixel), TString::Format("Pixel %d (ID %d): branch dt;Stored pulse.ecal_residual (ns);Pulses", pixel, id), nDtBins, DtMin, DtMax);
      hPulseSummaryPixelLe[selected][pixel] = new TH1D(TString::Format("hPulseSummary%sPixel%dLe", roles[selected], pixel), TString::Format("Pixel %d (ID %d): raw LE;Raw CDet LE (ns);Pulses", pixel, id), nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelLeCorr[selected][pixel] = new TH1D(TString::Format("hPulseSummary%sPixel%dLeCorr", roles[selected], pixel), TString::Format("Pixel %d (ID %d): corrected LE;Corrected CDet LE (ns);Pulses", pixel, id), nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelDt[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLe[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLeCorr[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLeVsToT[selected][pixel] = new TH2D(TString::Format("hPulseSummary%sPixel%dLeVsToT", roles[selected], pixel), TString::Format("Pixel %d (ID %d): raw LE versus ToT;ToT (ns);Raw CDet LE (ns);Pulses", pixel, id), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelLeCorrVsToT[selected][pixel] = new TH2D(TString::Format("hPulseSummary%sPixel%dLeCorrVsToT", roles[selected], pixel), TString::Format("Pixel %d (ID %d): corrected LE versus ToT;ToT (ns);Corrected CDet LE (ns);Pulses", pixel, id), nTotBins, TotDisplayMin, TotDisplayMax, nLeBins, LeDisplayMin, LeDisplayMax);
      hPulseSummaryPixelLeVsToT[selected][pixel]->SetDirectory(nullptr);
      hPulseSummaryPixelLeCorrVsToT[selected][pixel]->SetDirectory(nullptr);
      hPulsePixelDtVsRx[selected][pixel] = new TH2D(TString::Format("hPulseSummary%sPixel%dDtVsRx", roles[selected], pixel), TString::Format("Pixel %d (ID %d): branch dt versus rx;rx (m);Stored pulse.ecal_residual (ns);Pulses", pixel, id), 100, -xRes, xRes, nDtBins, DtMin, DtMax);
      hPulsePixelDtVsRx[selected][pixel]->SetDirectory(nullptr);
    }
  }

  if (ApplyYCorrection) {
    for (int selected = 0; selected < 2; ++selected) {
      for (TH1* hist : {hPulseSummaryBarDt[selected], hPulseFullBarDtNoLeCut[selected], hPulseFullBarDtAllCuts[selected]})
        hist->GetXaxis()->SetTitle("Y-corrected dt (ns)");
      for (TH2* hist : {hPulseBarDtVsTot[selected], hPulseBarDtVsTotNoTotCut[selected], hPulseBarDtVsRx[selected]})
        hist->GetYaxis()->SetTitle("Y-corrected dt (ns)");
      hPulseSummaryBarLeCorr[selected]->GetXaxis()->SetTitle("Replay + y-corrected CDet LE (ns)");
      hPulseSummaryBarLeCorrVsToT[selected]->GetYaxis()->SetTitle("Replay + y-corrected CDet LE (ns)");
      for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
        hPulseSummaryPixelDt[selected][pixel]->GetXaxis()->SetTitle("Y-corrected dt (ns)");
        hPulsePixelDtVsRx[selected][pixel]->GetYaxis()->SetTitle("Y-corrected dt (ns)");
        hPulseSummaryPixelLeCorr[selected][pixel]->GetXaxis()->SetTitle("Replay + y-corrected CDet LE (ns)");
        hPulseSummaryPixelLeCorrVsToT[selected][pixel]->GetYaxis()->SetTitle("Replay + y-corrected CDet LE (ns)");
      }
    }
    hPulseBar30LeVsECalTime[1]->GetYaxis()->SetTitle("Replay + y-corrected CDet LE (ns)");
  }

  // One pass through the aligned saved vectors. dt is already in ns.
  Long64_t filledPulses = 0;
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    const bool good_ECal_energy = std::isfinite(vECalE[event]) &&
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && (!applyECalEnergyMax || vECalE[event] < ECalMaxE)));
    const bool good_ECal_time = std::isfinite(vECalAdcTime[event]) &&
        (!applyECalTimeCut || (vECalAdcTime[event] > ECalMinT && vECalAdcTime[event] < ECalMaxT));
    if (!good_ECal_energy || !good_ECal_time) continue;
    double ecal_y = vECalY[event];
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit];
      const double dt = vPulseDt[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars * pixelsPerBar ||
          std::floor(id) != id || !std::isfinite(dt)) continue;
      const int bar = int(id) / pixelsPerBar;
      const int pixel = int(id) % pixelsPerBar;

      if (bar != studyBar && bar != controlBar && bar != 30) continue;
      // All saved timing vectors used here are in ns.
      const double le = vPulseLe[event][hit];
      const double yDtCorrection = ApplyYCorrection
          ? PulseStudyYTimingCorrection(int(id), vPulseY[event][hit],
              vECalYProjection[event][hit], YPropagationNsPerMeter) : 0.0;
      const double plotDt = dt + yDtCorrection;
      const double leCorr = vPulseLeCorr[event][hit] - yDtCorrection; // ns; saved vectors unchanged.
      const double totNs = vPulseToT[event][hit];
      const double rx = vPulseXResiduals[event][hit];
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(totNs) && totNs > 0.0 &&
          (!applyToTCut || (totNs >= TotMin && (!applyToTMax || totNs <= TotMax)));
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
      if (!std::isfinite(totNs) || totNs <= 0.0 || !good_ellipse) continue;
      // Fixed bar-30 diagnostics, even if neither selected summary bar is 30.
      // Each pulse gets its event's ECal time; corrected LE is already in ns.
      if (bar == 30 && good_le && good_tot) {
        hPulseBar30LeVsECalTime[0]->Fill(vECalAdcTime[event], le);
        if (std::isfinite(leCorr)) hPulseBar30LeVsECalTime[1]->Fill(vECalAdcTime[event], leCorr);
      }
      for (int selected = 0; selected < 2; ++selected) {
        if (bar != bars[selected]) continue;
        // Compare the secondary ToT window with all other cuts, including LE, fixed.
        // Pulses rejected by the main loading ToT window are not in these vectors.
        if (good_le) hPulseBarDtVsTotNoTotCut[selected]->Fill(totNs, plotDt);
        if (!good_tot) continue;
        // Same event, ToT, and spatial/ellipse cuts; omit only the raw-LE cut.
        hPulseFullBarDtNoLeCut[selected]->Fill(plotDt);
        if (!good_le) continue;
        // Deliberately use dt, not plotDt: preserve the propagation slope in this diagnostic.
        if (std::isfinite(vECalYProjection[event][hit]))
          hPulseBarDtVsECalYProjection[selected]->Fill(vECalYProjection[event][hit], dt);
        hPulseSummaryBarDt[selected]->Fill(plotDt);
        hPulseSummaryBarLe[selected]->Fill(le);
        hPulseSummaryPixelDt[selected][pixel]->Fill(plotDt);
        hPulseSummaryPixelLe[selected][pixel]->Fill(le);
        hPulseSummaryBarLeVsToT[selected]->Fill(totNs, le);
        hPulseSummaryPixelLeVsToT[selected][pixel]->Fill(totNs, le);
        // Select with RAW LE above; corrected LE is only another plotted value.
        if (std::isfinite(leCorr)) {
          hPulseSummaryBarLeCorr[selected]->Fill(leCorr);
          hPulseSummaryPixelLeCorr[selected][pixel]->Fill(leCorr);
          hPulseSummaryBarLeCorrVsToT[selected]->Fill(totNs, leCorr);
          hPulseSummaryPixelLeCorrVsToT[selected][pixel]->Fill(totNs, leCorr);
        }
        hPulseBarDtVsRx[selected]->Fill(rx, plotDt);
        hPulsePixelDtVsRx[selected][pixel]->Fill(rx, plotDt);
        hPulseBarDtVsTot[selected]->Fill(totNs, plotDt);
        hPulseFullBarDtAllCuts[selected]->Fill(plotDt);
      }
      if (good_le && good_tot && (bar == studyBar || bar == controlBar)) ++filledPulses;
    }
  }
  // Fourteen canvases per bar; overlay LE cuts on top and compare secondary ToT cuts below.
  for (int selected = 0; selected < 2; ++selected) {
    TCanvas* cComparison = new TCanvas(TString::Format("cPulseSummary%sCutComparison", roles[selected]), TString::Format("%s, %s bar %d: LE and ToT cut comparison", PulseStudyRunLabel(config, RunNumber).Data(), roles[selected], bars[selected]), 1600, 1400);
    // One wide dt overlay above the two ToT comparison panels.
    for (int panel = 0; panel < 3; ++panel) {
      const double x1 = panel == 2 ? 0.50 : 0.0;
      const double x2 = panel == 1 ? 0.50 : 1.0;
      TPad* pad = new TPad(TString::Format("pulseCutComparison%s%d", roles[selected], panel), "",
          x1, panel == 0 ? 0.50 : 0.0, x2, panel == 0 ? 1.0 : 0.50);
      pad->SetNumber(panel + 1);
      pad->Draw();
    }
    TH1* comparisonHistograms[] = {hPulseFullBarDtNoLeCut[selected], hPulseBarDtVsTotNoTotCut[selected], hPulseBarDtVsTot[selected]};
    const double dtMaximum = 1.15 * std::max(1.0, std::max(hPulseFullBarDtAllCuts[selected]->GetMaximum(), hPulseFullBarDtNoLeCut[selected]->GetMaximum()));
    const double dtVsTotMaximum = std::max(1.0, std::max(hPulseBarDtVsTotNoTotCut[selected]->GetMaximum(), hPulseBarDtVsTot[selected]->GetMaximum()));
    for (int panel = 0; panel < 3; ++panel) {
      cComparison->cd(panel + 1);
      gPad->SetTopMargin(0.10);
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(panel > 0 ? 0.16 : 0.05);
      TH1* histogram = comparisonHistograms[panel];
      histogram->SetStats(false);
      histogram->SetMinimum(0.0);
      if (panel == 0) {
        histogram->SetTitle(TString::Format("%s, bar %d: dt before/after secondary raw-LE cut;Stored pulse.ecal_residual (ns);Pulses",
            PulseStudyRunLabel(config, RunNumber).Data(), bars[selected]));
        if (ApplyYCorrection) histogram->GetXaxis()->SetTitle("Y-corrected dt (ns)");
        histogram->SetMaximum(dtMaximum);
        histogram->SetLineWidth(2);
        histogram->SetLineColor(kRed + 1);
        histogram->Draw("HIST");
        TH1D* after = hPulseFullBarDtAllCuts[selected];
        after->SetStats(false);
        after->SetLineWidth(2);
        after->SetLineColor(kBlue + 1);
        after->SetLineStyle(2);
        after->Draw("HIST SAME");
        TLegend* legend = new TLegend(0.66, 0.75, 0.94, 0.90);
        legend->SetBit(kCanDelete);
        legend->SetFillColor(kWhite);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.028);
        legend->AddEntry(histogram, "Before secondary raw-LE cut", "l");
        legend->AddEntry(after, "After secondary raw-LE cut", "l");
        legend->Draw();
      } else {
        histogram->SetMaximum(dtVsTotMaximum);
        histogram->Draw("COLZ");
      }
    }
    cComparison->Modified();
    cComparison->Update();
    PulseStudySaveCanvas(cComparison, TString::Format("%s/bar%03d_%s_le_cut_comparison.pdf", outputDirectory.Data(), bars[selected], roles[selected]));
    TCanvas* cBarDt = new TCanvas(TString::Format("cPulseSummary%sBarDt", roles[selected]), TString::Format("%s bar %d: summed branch dt", roles[selected], bars[selected]), 1000, 850);
    cBarDt->SetTopMargin(0.10);
    hPulseSummaryBarDt[selected]->Draw("HIST");
    cBarDt->Modified();
    cBarDt->Update();

    TCanvas* cBarLe = new TCanvas(TString::Format("cPulseSummary%sBarLe", roles[selected]), TString::Format("%s bar %d: summed raw LE", roles[selected], bars[selected]), 1000, 850);
    cBarLe->SetTopMargin(0.10);
    hPulseSummaryBarLe[selected]->Draw("HIST");
    cBarLe->Modified();
    cBarLe->Update();

    TCanvas* cBarLeCorr = new TCanvas(TString::Format("cPulseSummary%sBarLeCorr", roles[selected]), TString::Format("%s bar %d: summed corrected LE", roles[selected], bars[selected]), 1000, 850);
    cBarLeCorr->SetTopMargin(0.10);
    hPulseSummaryBarLeCorr[selected]->Draw("HIST");
    cBarLeCorr->Modified();
    cBarLeCorr->Update();

    TCanvas* cPixelDt = new TCanvas(TString::Format("cPulseSummary%sPixelDt", roles[selected]), TString::Format("%s bar %d: pixel branch dt", roles[selected], bars[selected]), 1400, 1000);
    cPixelDt->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelDt->cd(pixel + 1);
      hPulseSummaryPixelDt[selected][pixel]->Draw("HIST");
    }
    cPixelDt->Modified();
    cPixelDt->Update();

    TCanvas* cPixelLe = new TCanvas(TString::Format("cPulseSummary%sPixelLe", roles[selected]), TString::Format("%s bar %d: pixel raw LE", roles[selected], bars[selected]), 1400, 1000);
    cPixelLe->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelLe->cd(pixel + 1);
      hPulseSummaryPixelLe[selected][pixel]->Draw("HIST");
    }
    cPixelLe->Modified();
    cPixelLe->Update();

    TCanvas* cPixelLeCorr = new TCanvas(TString::Format("cPulseSummary%sPixelLeCorr", roles[selected]), TString::Format("%s bar %d: pixel corrected LE", roles[selected], bars[selected]), 1400, 1000);
    cPixelLeCorr->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelLeCorr->cd(pixel + 1);
      hPulseSummaryPixelLeCorr[selected][pixel]->Draw("HIST");
    }
    cPixelLeCorr->Modified();
    cPixelLeCorr->Update();

    TCanvas* cBarLeVsToT = new TCanvas(TString::Format("cPulseSummary%sBarLeVsToT", roles[selected]), TString::Format("%s bar %d: summed raw LE versus ToT", roles[selected], bars[selected]), 1000, 850);
    cBarLeVsToT->SetTopMargin(0.10);
    hPulseSummaryBarLeVsToT[selected]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hPulseSummaryBarLeVsToT[selected]->Draw("COLZ");
    cBarLeVsToT->Modified();
    cBarLeVsToT->Update();

    TCanvas* cBarLeCorrVsToT = new TCanvas(TString::Format("cPulseSummary%sBarLeCorrVsToT", roles[selected]), TString::Format("%s bar %d: summed corrected LE versus ToT", roles[selected], bars[selected]), 1000, 850);
    cBarLeCorrVsToT->SetTopMargin(0.10);
    hPulseSummaryBarLeCorrVsToT[selected]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hPulseSummaryBarLeCorrVsToT[selected]->Draw("COLZ");
    cBarLeCorrVsToT->Modified();
    cBarLeCorrVsToT->Update();

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

    TCanvas* cPixelLeCorrVsToT = new TCanvas(TString::Format("cPulseSummary%sPixelLeCorrVsToT", roles[selected]), TString::Format("%s bar %d: pixel corrected LE versus ToT", roles[selected], bars[selected]), 1400, 1000);
    cPixelLeCorrVsToT->Divide(4, 4);
    for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
      cPixelLeCorrVsToT->cd(pixel + 1);
      hPulseSummaryPixelLeCorrVsToT[selected][pixel]->SetStats(false);
      gPad->SetRightMargin(0.15);
      hPulseSummaryPixelLeCorrVsToT[selected][pixel]->Draw("COLZ");
    }
    cPixelLeCorrVsToT->Modified();
    cPixelLeCorrVsToT->Update();

    TCanvas* cBarDtVsRx = new TCanvas(TString::Format("cPulseSummary%sBarDtVsRx", roles[selected]), TString::Format("%s bar %d: summed branch dt versus rx", roles[selected], bars[selected]), 1000, 850);
    cBarDtVsRx->SetTopMargin(0.10);
    hPulseBarDtVsRx[selected]->SetStats(false);
    gPad->SetRightMargin(0.15);
    hPulseBarDtVsRx[selected]->Draw("COLZ");
    cBarDtVsRx->Modified();
    cBarDtVsRx->Update();

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

    TCanvas* cDtVsECalYProjection = new TCanvas(TString::Format("cPulseSummary%sDtVsECalYProjection", roles[selected]),
        hPulseBarDtVsECalYProjection[selected]->GetTitle(), 1200, 742);
    cDtVsECalYProjection->SetRightMargin(0.15);
    TH2D* yHistogram = hPulseBarDtVsECalYProjection[selected];
    int firstYBin = nYProjectionBins+1, lastYBin = 0;
    for (int xbin = 1; xbin <= nYProjectionBins; ++xbin) {
      // Include dt flow bins when finding the occupied projected-y range.
      if (yHistogram->Integral(xbin, xbin, 0, nDtBins+1) <= 0.0) continue;
      firstYBin = std::min(firstYBin, xbin); lastYBin = std::max(lastYBin, xbin);
    }
    if (lastYBin >= firstYBin) yHistogram->GetXaxis()->SetRange(std::max(1, firstYBin-2), std::min(nYProjectionBins, lastYBin+2));
    yHistogram->Draw("COLZ");

    // Full-selection plots keep only a run/bar header; cuts are saved alongside PDFs.
    TCanvas* labeledCanvases[] = {cBarDt, cBarLe, cPixelDt, cPixelLe,
        cBarLeVsToT, cPixelLeVsToT, cBarDtVsRx, cPixelDtVsRx,
        cBarLeCorr, cPixelLeCorr, cBarLeCorrVsToT, cPixelLeCorrVsToT, cDtVsECalYProjection};
    const char* plotNames[] = {"dt", "raw_le", "pixels_dt", "pixels_raw_le",
        "le_vs_tot", "pixels_le_vs_tot", "dt_vs_rx", "pixels_dt_vs_rx",
        "le_corr", "pixels_le_corr", "le_corr_vs_tot", "pixels_le_corr_vs_tot", "dt_before_y_vs_ecal_y_projected"};
    for (int plot = 0; plot < 13; ++plot) {
      TCanvas* canvas = labeledCanvases[plot];
      const bool pixelGrid = plot == 2 || plot == 3 || plot == 5 || plot == 7 || plot == 9 || plot == 11;
      canvas->cd();
      if (pixelGrid) {
        canvas->SetCanvasSize(1400, 1000);
        for (int pixel = 0; pixel < pixelsPerBar; ++pixel) {
          TPad* pad = static_cast<TPad*>(canvas->GetPad(pixel + 1));
          const double x = pad->GetXlowNDC(), y = pad->GetYlowNDC();
          const double w = pad->GetWNDC(), h = pad->GetHNDC();
          pad->SetPad(x, 0.95 * y, x + w, 0.95 * (y + h));
          pad->Modified();
        }
        TLatex* title = new TLatex(0.5, 0.975, TString::Format("%s: %s bar %d", PulseStudyRunLabel(config, RunNumber).Data(), roles[selected], bars[selected]));
        title->SetNDC(); title->SetTextAlign(22); title->SetTextSize(0.025);
        title->SetBit(kCanDelete); title->Draw();
      } else {
        canvas->SetCanvasSize(1200, 742);
        canvas->SetTopMargin(0.10);
        TIter next(canvas->GetListOfPrimitives());
        while (TObject* object = next()) {
          if (auto* histogram = dynamic_cast<TH1*>(object)) {
            histogram->SetTitle(PulseStudyRunLabel(config, RunNumber) + ": " + histogram->GetTitle());
            histogram->SetStats(false);
          }
        }
      }
      canvas->Modified();
      canvas->Update();
      PulseStudySaveCanvas(canvas, TString::Format("%s/bar%03d_%s_%s.pdf", outputDirectory.Data(), bars[selected], roles[selected], plotNames[plot]));
    }
    if (selected == 1) {
      const double colorMax = std::max(1.0, std::max(hPulseBar30LeVsECalTime[0]->GetMaximum(), hPulseBar30LeVsECalTime[1]->GetMaximum()));
      const char* files[] = {"bar030_raw_le_vs_ecal_time.pdf", "bar030_le_corr_vs_ecal_time.pdf"};
      for (int kind = 0; kind < 2; ++kind) {
        TCanvas* canvas = new TCanvas(leECalCanvasNames[kind], hPulseBar30LeVsECalTime[kind]->GetTitle(), 1000, 850);
        canvas->SetCanvasSize(1000, 850);
        canvas->SetTopMargin(0.10);
        canvas->SetLeftMargin(0.12);
        canvas->SetRightMargin(0.15);
        hPulseBar30LeVsECalTime[kind]->SetMinimum(0.0);
        hPulseBar30LeVsECalTime[kind]->SetMaximum(colorMax);
        hPulseBar30LeVsECalTime[kind]->Draw("COLZ");
        canvas->Modified();
        canvas->Update();
        PulseStudySaveCanvas(canvas, outputDirectory + "/" + files[kind]);
      }
    }
  }
  if (ApplyYCorrection) std::cout << "[CDet pulse bar summary] Y correction ON for both selected bars and fixed bar 30; raw LE and stored vectors unchanged.\n";
  std::cout << "[CDet pulse bar summary] Drew study bar " << studyBar << " and control bar " << controlBar
            << "; selected " << filledPulses << " pulses. Thirty canvases remain open (including two pre-y-correction dt/y plots and two bar-30 ECal-time plots).\n"
            << "  dt starts from vPulseDt; corrected LE starts from vPulseLeCorr; optional y correction applies to both.\n"
            << "  LE cuts still use raw vPulseLe; both LE displays share the configured range/binning.\n"
            << "  dt versus rx uses the plotted dt against stored x residual.\n"
            << "  PDF output: " << outputDirectory << '\n';
}

// Convenient opt-in call, with optional study/control bars; numeric bar calls remain compatible.
void plotPulseStudyBars(const char* configFile, const char* outputTag, bool ApplyYCorrection,
                        int studyBar = -1, int controlBar = -1)
{
  plotPulseStudyBars(configFile, outputTag, studyBar, controlBar,
      NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, ApplyYCorrection);
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
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const bool applyECalEnergyMax = config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0;
  const bool applyToTMax = config.GetValue("pulse_study.apply_tot_max", 0) != 0;
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if ((applyECalEnergyCut && (!std::isfinite(ECalMinE) || (applyECalEnergyMax && (!std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)))) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || TotMin < 0 || (applyToTMax && (!std::isfinite(TotMax) || TotMin >= TotMax))))) {
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
  if (!PulseStudyWriteCutReport(outputDirectory + "/plotPulseStudyPositions_cuts.txt", "plotPulseStudyPositions", configFile, config,
      PulseStudySecondaryCuts(config), {
        TString::Format("Study bar %d; control bar %d; x source: %s; finite plotted x required", studyBar, controlBar, UseCorrectedX ? "aligned x" : "channel x"),
        TString::Format("Position display [%.12g,%.12g] m, requested bin width %.12g m", PositionMin, PositionMax, PositionBinWidth),
        "Both x and rx histograms count the same selected pulses. Y timing correction OFF."
      })) return;


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
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && (!applyECalEnergyMax || vECalE[event] < ECalMaxE)));
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
          (!applyToTCut || (totNs >= TotMin && (!applyToTMax || totNs <= TotMax)));
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
    PulseStudySaveCanvas(canvas, TString::Format("%s/bar%03d_%s_x_and_rx.pdf", outputDirectory.Data(), bars[selected], roles[selected]));
  }
  std::cout << "[CDet pulse positions] Drew study bar " << studyBar << " and control bar " << controlBar
            << "; selected " << filledPulses << " pulses. Two position canvases remain open.\n"
            << "  Blue: " << xLabel << "; red: stored rx. Positions are in metres.\n"
            << "  Values outside the display bounds remain in underflow/overflow bins.\n"
            << "  PDF output: " << outputDirectory << '\n';
}

// Native and projected ECal coordinates versus raw/aligned CDet x and DB y.
// Same event/vector/cut architecture as plotPulseStudyBars; additional same-event layer combinations use the bar-30 signal gate.
void plotPulseStudyPositionCorrelations(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                            const char* outputTag = nullptr, int controlBar = -1,
                            double PositionBinWidth = NAN, double PositionMin = NAN, double PositionMax = NAN)
{
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet position correlations] Cannot read configuration with config.version: 1.\n";
    return;
  }
  if (controlBar == -1) controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  if (std::isnan(PositionBinWidth)) PositionBinWidth = config.GetValue("pulse_study.position_bin_width_m", 0.005);
  if (std::isnan(PositionMin)) PositionMin = config.GetValue("pulse_study.position_min_m", -1.6);
  if (std::isnan(PositionMax)) PositionMax = config.GetValue("pulse_study.position_max_m", 1.6);
  const double ControlPositionPadding = config.GetValue("pulse_study.control_position_padding_m", 0.02);
  if (!std::isfinite(ControlPositionPadding) || ControlPositionPadding < 0.0) {
    std::cerr << "[CDet position correlations] Control-position padding must be finite and nonnegative.\n";
    return;
  }
  const bool applyDtCut = config.GetValue("pulse_study.position_apply_dt_cut", 0) != 0;
  const double DtCutMin = config.GetValue("pulse_study.position_dt_min_ns", config.GetValue("pulse_study.fit_signal_min", -36.0));
  const double DtCutMax = config.GetValue("pulse_study.position_dt_max_ns", config.GetValue("pulse_study.fit_signal_max", -20.0));
  if (!std::isfinite(DtCutMin) || !std::isfinite(DtCutMax) || DtCutMin >= DtCutMax) {
    std::cerr << "[CDet position correlations] The bar-30 signal gate requires finite position_dt_min_ns < position_dt_max_ns.\n";
    return;
  }
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  // Loading applied energy/time and replay ECal eligibility; these cuts only narrow saved events.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const bool applyECalEnergyMax = config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0;
  const bool applyToTMax = config.GetValue("pulse_study.apply_tot_max", 0) != 0;
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  const bool useEllipse = ellipseRadius > 0.0;

  if ((applyECalEnergyCut && (!std::isfinite(ECalMinE) || (applyECalEnergyMax && (!std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE)))) ||
      (applyECalTimeCut && (!std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT)) ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToTCut && (!std::isfinite(TotMin) || TotMin < 0 || (applyToTMax && (!std::isfinite(TotMax) || TotMin >= TotMax))))) {
    std::cerr << "[CDet position correlations] Invalid configured ECal, LE, or ToT cuts.\n";
    return;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet position correlations] Invalid ellipse settings; radius -1 disables it.\n";
    return;
  }
  if (controlBar < 0 || controlBar >= numBars ||
      !std::isfinite(PositionBinWidth) || PositionBinWidth <= 0 ||
      !std::isfinite(PositionMin) || !std::isfinite(PositionMax) || PositionMin >= PositionMax ||
      !std::isfinite((PositionMax - PositionMin) / PositionBinWidth) ||
      std::ceil((PositionMax - PositionMin) / PositionBinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet position correlations] Invalid bar IDs or position display settings.\n";
    return;
  }
  if (vPulseID.empty()) {
    std::cerr << "[CDet position correlations] No saved pulse events; run the main macro first.\n";
    return;
  }
  if (!std::isfinite(xRes) || xRes <= 0.0 ||
      std::pow(std::ceil((PositionMax-PositionMin)/PositionBinWidth)+2.0, 2) > std::numeric_limits<int>::max()) {
    std::cerr << "[CDet position correlations] Invalid rx cut or too many 2D histogram cells.\n";
    return;
  }
  // Both coordinate versions use the same saved events and pulse selection.
  if (vECalX.size() != vPulseID.size() || vECalY.size() != vPulseID.size()) {
    std::cerr << "[CDet position correlations] Missing event ECal positions; rerun main analysis.\n";
    return;
  }
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    for (const auto* values : {&vPulseX, &vPulseXCorr, &vPulseY, &vECalXProjection, &vECalYProjection}) {
      if (event >= values->size() || (*values)[event].size() != vPulseID[event].size()) {
        std::cerr << "[CDet position correlations] Position vectors are not event/hit aligned; rerun main analysis.\n";
        return;
      }
    }
  }
  const double DifferenceMin = config.GetValue("pulse_study.position_difference_min_m", -0.20);
  const double DifferenceMax = config.GetValue("pulse_study.position_difference_max_m", 0.20);
  if (!std::isfinite(DifferenceMin) || !std::isfinite(DifferenceMax) || DifferenceMin >= DifferenceMax ||
      !std::isfinite((DifferenceMax-DifferenceMin)/PositionBinWidth) ||
      std::ceil((DifferenceMax-DifferenceMin)/PositionBinWidth) > std::numeric_limits<int>::max()-2) {
    std::cerr << "[CDet position correlations] Invalid position-difference display range.\n";
    return;
  }
  const double LayerDifferenceMin = config.GetValue("pulse_study.layer_difference_min_m", -1.6);
  const double LayerDifferenceMax = config.GetValue("pulse_study.layer_difference_max_m", 1.6);
  if (!std::isfinite(LayerDifferenceMin) || !std::isfinite(LayerDifferenceMax) || LayerDifferenceMin >= LayerDifferenceMax ||
      !std::isfinite((LayerDifferenceMax-LayerDifferenceMin)/PositionBinWidth) ||
      std::ceil((LayerDifferenceMax-LayerDifferenceMin)/PositionBinWidth) > std::numeric_limits<int>::max()-2) {
    std::cerr << "[CDet position correlations] Invalid layer x2-x1 display range.\n";
    return;
  }
  const int nLayerDifferenceBins = int(std::ceil((LayerDifferenceMax-LayerDifferenceMin)/PositionBinWidth));
  const int nDifferenceBins = int(std::ceil((DifferenceMax-DifferenceMin)/PositionBinWidth));
  const int nPositionBins = int(std::ceil((PositionMax - PositionMin) / PositionBinWidth));
  const char* roles[4] = {"Control", "Detector", "Layer1", "Layer2"};
  const TString scopes[4] = {TString::Format("Control bar %d", controlBar),
      "Full detector (both layers)", "CDet layer 1 (bars 0-83)", "CDet layer 2 (bars 84-167)"};
  const TString filePrefixes[4] = {TString::Format("bar%03d_Control", controlBar),
      "detector", "layer1", "layer2"};
  TString tag = outputTag ? outputTag : config.GetValue("pulse_study.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d", config.GetValue("analysis.run_number", 0));
  const TString allowedTagCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowedTagCharacters.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet position correlations] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/PusleStudy/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet position correlations] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }
  if (!PulseStudyWriteCutReport(outputDirectory + "/plotPulseStudyPositionCorrelations_cuts.txt", "plotPulseStudyPositionCorrelations", configFile, config,
      PulseStudySecondaryCuts(config), {
        applyDtCut ? TString::Format("Additional dt cut: %.12g <= pulse.ecal_residual <= %.12g ns, BEFORE any y correction; applies to the single-pulse position plots (layer combinations and multiplicity have the exceptions below)", DtCutMin, DtCutMax) : TString("Additional position-correlation dt cut: OFF"),
        "Fixed bar-30 channel/aligned x-minus-ECal-projection spectra use the same cuts and difference binning as layer spectra, independently of controlBar.",
        TString::Format("Control bar %d; full detector; layer 1 (bars 0..83); layer 2 (84..167)", controlBar),
        TString::Format("Position display [%.12g,%.12g] m, requested bin width %.12g m; control zoom padding %.12g m (at least one bin)", PositionMin, PositionMax, PositionBinWidth, ControlPositionPadding),
        TString::Format("x-minus-projection display [%.12g,%.12g] m", DifferenceMin, DifferenceMax),
        "Finite plotted coordinates required per histogram; each accepted pulse counted once per applicable scope.",
        "ECal projected coordinate = native ECal coordinate * pixel z / ECal z. Native ECal coordinates are unprojected.",
        "Raw/corrected CDet x vs projected ECal x overlays: ProfileX mean CDet x in each ECal-x bin, with errors on the mean, computed from regular TH2 bins (bin-center approximation; excludes under/overflow).",
        "CDet y is the DB channel center; x-minus-projection spectra use direct subtraction, without an extra residual offset.",
        "1D channel CDet x: one entry per accepted pulse. 1D native ECal x: one entry per saved event with an accepted pulse, not weighted by pulse multiplicity.",
        "The additional channel-CDet-x / native-ECal-x overlay retains these counts without normalization; blue counts pulses, red counts events.",
        "Y timing correction OFF. New layer plots are same-event pulse combinations, not replay pair-branch matches.",
        TString::Format("Layer-combination gate: at least one L1 bar-30 pulse passing secondary ECal, raw-LE, ToT and rx/ellipse cuts AND %.12g <= stored dt <= %.12g ns. This gate is ALWAYS applied, regardless of position_apply_dt_cut.", DtCutMin, DtCutMax),
        "Combine EVERY qualifying bar-30 pulse with EVERY saved L2 pulse in the same event (IDs 1344..2687), including L2 pulses outside the dt peak or failing secondary pulse cuts. L2 has loading cuts only; event ECal cuts still apply.",
        "Plot x1 horizontally and x2 vertically; 1D difference is x2-x1. Channel and aligned versions use their respective saved x vectors; finite coordinates required independently. Counts are pulse combinations, not events or matched tracks.",
        TString::Format("Layer x2-x1 display [%.12g,%.12g] m, requested bin width %.12g m; no ECal projection or residual offset subtraction", LayerDifferenceMin, LayerDifferenceMax, PositionBinWidth),
        "Bar-30 multiplicity: all saved bar-30 pulses in each event passing secondary ECal cuts, before secondary pulse/dt cuts. Fill only N>=1; one entry per event. Axis shows total N (additional hits = N-1)."

      })) return;


  const int RunNumber = pulseStudyRunNumber > 0 ? pulseStudyRunNumber : config.GetValue("analysis.run_number", 0);
  const char* coordinateNames[] = {"RawX", "CorrectedX", "Y"};
  const char* coordinateFiles[] = {"raw_x", "corrected_x", "y"};
  const char* cdetLabels[] = {"Channel CDet x (m)", "Aligned CDet x (m)", "DB CDet y (m)"};
  const char* ecalNames[] = {"Native", "Projected"};
  // Replace this routine's plots on each call, then define, fill, and draw.
  for (int selected = 0; selected < 4; ++selected) {
    for (int projected = 0; projected < 2; ++projected) {
      for (int coordinate = 0; coordinate < 3; ++coordinate) {
        delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseCorrelation%s%s%s", roles[selected], ecalNames[projected], coordinateNames[coordinate]));
        delete hPulsePositionCorrelation[selected][projected][coordinate];
        const TString ecalLabel = TString::Format("ECal %s%s (m)", coordinate == 2 ? "y" : "x", projected ? " projected to pixel z" : " at ECal");
        auto*& histogram = hPulsePositionCorrelation[selected][projected][coordinate];
        histogram = new TH2D(TString::Format("hPulseCorrelation%s%s%s", roles[selected], ecalNames[projected], coordinateNames[coordinate]),
            TString::Format("%s, %s: %s vs ECal;%s;%s;Pulses", PulseStudyRunLabel(config, RunNumber).Data(), scopes[selected].Data(), coordinateNames[coordinate], ecalLabel.Data(), cdetLabels[coordinate]),
            nPositionBins, PositionMin, PositionMax, nPositionBins, PositionMin, PositionMax);
        histogram->SetDirectory(nullptr);
        histogram->SetStats(false);
      }
    }
  }

  for (int layer = 0; layer < 2; ++layer) {
    for (int coordinate = 0; coordinate < 2; ++coordinate) {
      delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseLayer%d%sDifference", layer+1, coordinateNames[coordinate]));
      delete hPulseLayerXDifference[layer][coordinate];
      auto*& histogram = hPulseLayerXDifference[layer][coordinate];
      histogram = new TH1D(TString::Format("hPulseLayer%d%sDifference", layer+1, coordinateNames[coordinate]),
          TString::Format("%s, CDet layer %d: %s x - ECal projection;%s - x_{ECal,proj} (m);Pulses",
              PulseStudyRunLabel(config, RunNumber).Data(), layer+1, coordinate ? "aligned" : "channel", coordinate ? "x_{corr}" : "x"),
          nDifferenceBins, DifferenceMin, DifferenceMax);
      histogram->SetDirectory(nullptr);
      histogram->SetStats(false);
    }
  }

  for (int coordinate = 0; coordinate < 2; ++coordinate) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseBar30%sDifference", coordinateNames[coordinate]));
    delete hPulseBar30XDifference[coordinate];
    hPulseBar30XDifference[coordinate] = new TH1D(TString::Format("hPulseBar30%sDifference", coordinateNames[coordinate]),
        TString::Format("%s, bar 30: %s x - ECal projection;%s - x_{ECal,proj} (m);Pulses",
            PulseStudyRunLabel(config, RunNumber).Data(), coordinate ? "aligned" : "channel", coordinate ? "x_{corr}" : "x"),
        nDifferenceBins, DifferenceMin, DifferenceMax);
    hPulseBar30XDifference[coordinate]->SetDirectory(nullptr);
    hPulseBar30XDifference[coordinate]->SetStats(false);
  }

  for (int coordinate = 0; coordinate < 2; ++coordinate) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseLayers%sX2VsX1", coordinateNames[coordinate]));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseLayers%sX2MinusX1", coordinateNames[coordinate]));
    delete hPulseLayerX2VsX1[coordinate];
    delete hPulseLayerX2MinusX1[coordinate];
    hPulseLayerX2VsX1[coordinate] = new TH2D(TString::Format("hPulseLayers%sX2VsX1", coordinateNames[coordinate]),
        TString::Format("%s: %s x_{2} vs x_{1}, bar-30 signal gate;Bar-30 x_{1} (m);Layer-2 x_{2} (m);Pulse combinations", PulseStudyRunLabel(config, RunNumber).Data(), coordinate ? "aligned" : "channel"),
        nPositionBins, PositionMin, PositionMax, nPositionBins, PositionMin, PositionMax);
    hPulseLayerX2MinusX1[coordinate] = new TH1D(TString::Format("hPulseLayers%sX2MinusX1", coordinateNames[coordinate]),
        TString::Format("%s: %s x_{2}-x_{1}, bar-30 signal gate;x_{2}-x_{1} (m);Pulse combinations", PulseStudyRunLabel(config, RunNumber).Data(), coordinate ? "aligned" : "channel"),
        nLayerDifferenceBins, LayerDifferenceMin, LayerDifferenceMax);
    hPulseLayerX2VsX1[coordinate]->SetDirectory(nullptr);
    hPulseLayerX2MinusX1[coordinate]->SetDirectory(nullptr);
    hPulseLayerX2VsX1[coordinate]->SetStats(false);
    hPulseLayerX2MinusX1[coordinate]->SetStats(false);
  }
  delete gROOT->GetListOfCanvases()->FindObject("cPulseBar30Multiplicity");
  delete hPulseBar30Multiplicity;
  hPulseBar30Multiplicity = nullptr;
  // Keep one multiplicity per event, then choose integer bins covering the observed maximum.
  std::vector<int> bar30Multiplicities;
  Long64_t bar30SignalEvents = 0;

  delete gROOT->GetListOfCanvases()->FindObject("cPulsePositionCDetX");
  delete gROOT->GetListOfCanvases()->FindObject("cPulsePositionECalX");
  delete gROOT->GetListOfCanvases()->FindObject("cPulsePositionNativeXOverlay");
  delete hPulsePositionCDetX;
  delete hPulsePositionECalX;
  hPulsePositionCDetX = new TH1D("hPulsePositionCDetX",
      TString::Format("%s, full detector: channel CDet x;Channel CDet x (m);Pulses", PulseStudyRunLabel(config, RunNumber).Data()),
      nPositionBins, PositionMin, PositionMax);
  hPulsePositionECalX = new TH1D("hPulsePositionECalX",
      TString::Format("%s: native ECal x, events with accepted CDet pulses;ECal x at ECal (m);Events", PulseStudyRunLabel(config, RunNumber).Data()),
      nPositionBins, PositionMin, PositionMax);
  hPulsePositionCDetX->SetDirectory(nullptr);
  hPulsePositionECalX->SetDirectory(nullptr);
  hPulsePositionCDetX->SetStats(false);
  hPulsePositionECalX->SetStats(false);

  // One pass through the aligned saved vectors. dt is already in ns.
  Long64_t filledPulses = 0;
  for (size_t event = 0; event < vPulseID.size(); ++event) {
    const bool good_ECal_energy = std::isfinite(vECalE[event]) &&
        (!applyECalEnergyCut || (vECalE[event] > ECalMinE && (!applyECalEnergyMax || vECalE[event] < ECalMaxE)));
    const bool good_ECal_time = std::isfinite(vECalAdcTime[event]) &&
        (!applyECalTimeCut || (vECalAdcTime[event] > ECalMinT && vECalAdcTime[event] < ECalMaxT));
    if (!good_ECal_energy || !good_ECal_time) continue;
    bool eventHasSelectedPulse = false;
    int bar30HitCount = 0;
    std::vector<size_t> bar30SignalHits, layer2Hits;
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit];
      const double dt = vPulseDt[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars * pixelsPerBar || std::floor(id) != id) continue;
      const int bar = int(id) / pixelsPerBar;
      // Collect the event-level populations BEFORE secondary pulse cuts.
      if (bar == 30) ++bar30HitCount;
      if (int(id) >= 1344) layer2Hits.push_back(hit);
      if (!std::isfinite(dt)) continue;
      if (applyDtCut && (dt < DtCutMin || dt > DtCutMax)) continue;

      // All saved timing vectors used here are in ns.
      const double le = vPulseLe[event][hit];
      const double totNs = vPulseToT[event][hit];
      const double rx = vPulseXResiduals[event][hit];
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(totNs) && totNs > 0.0 &&
          (!applyToTCut || (totNs >= TotMin && (!applyToTMax || totNs <= TotMax)));
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
      if (bar == 30 && dt >= DtCutMin && dt <= DtCutMax) bar30SignalHits.push_back(hit);
      eventHasSelectedPulse = true;
      if (std::isfinite(vPulseX[event][hit])) hPulsePositionCDetX->Fill(vPulseX[event][hit]);
      // Fill once per accepted pulse, independent of the control/detector display groups.
      // These are coordinate differences, with no additional residual-offset subtraction.
      const int layer = int(id) / 1344;
      const double xProj = vECalXProjection[event][hit];
      const double xValues[] = {vPulseX[event][hit], vPulseXCorr[event][hit]};
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        const double difference = xValues[coordinate] - xProj;
        if (std::isfinite(xProj) && std::isfinite(xValues[coordinate]) && std::isfinite(difference)) {
          hPulseLayerXDifference[layer][coordinate]->Fill(difference);
          if (bar == 30) hPulseBar30XDifference[coordinate]->Fill(difference);
        }
      }
      for (int selected = 0; selected < 4; ++selected) {
        if (selected == 0 && bar != controlBar) continue;
        // Logical IDs 0..1343 are L1; 1344..2687 are L2. No paired hits required.
        if (selected >= 2 && int(id) / 1344 != selected - 2) continue;
        const double cdet[] = {vPulseX[event][hit], vPulseXCorr[event][hit], vPulseY[event][hit]};
        for (int projected = 0; projected < 2; ++projected) {
          const double ecalX = projected ? vECalXProjection[event][hit] : vECalX[event];
          const double ecalY = projected ? vECalYProjection[event][hit] : vECalY[event];
          for (int coordinate = 0; coordinate < 3; ++coordinate) {
            const double ecal = coordinate == 2 ? ecalY : ecalX;
            if (std::isfinite(ecal) && std::isfinite(cdet[coordinate]))
              hPulsePositionCorrelation[selected][projected][coordinate]->Fill(ecal, cdet[coordinate]);
          }
        }
      }
      ++filledPulses;
    }
    if (bar30HitCount > 0) bar30Multiplicities.push_back(bar30HitCount);
    if (!bar30SignalHits.empty()) ++bar30SignalEvents;
    // All combinations, confined to this event. No L2 secondary dt/rx/LE/ToT cut.
    for (size_t hit1 : bar30SignalHits) {
      for (size_t hit2 : layer2Hits) {
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          const auto& x = coordinate == 0 ? vPulseX[event] : vPulseXCorr[event];
          const double x1 = x[hit1], x2 = x[hit2];
          if (!std::isfinite(x1) || !std::isfinite(x2)) continue;
          hPulseLayerX2VsX1[coordinate]->Fill(x1, x2);
          if (std::isfinite(x2-x1)) hPulseLayerX2MinusX1[coordinate]->Fill(x2-x1);
        }
      }
    }
    // Count each accepted event once, irrespective of pulse multiplicity or display bounds.
    if (eventHasSelectedPulse && std::isfinite(vECalX[event])) hPulsePositionECalX->Fill(vECalX[event]);
  }
  const int maxBar30Hits = bar30Multiplicities.empty() ? 1 : *std::max_element(bar30Multiplicities.begin(), bar30Multiplicities.end());
  hPulseBar30Multiplicity = new TH1D("hPulseBar30Multiplicity",
      TString::Format("%s: saved bar-30 hits per event (N #geq 1);Number of saved bar-30 pulses;Events", PulseStudyRunLabel(config, RunNumber).Data()),
      maxBar30Hits, 0.5, maxBar30Hits+0.5);
  hPulseBar30Multiplicity->SetDirectory(nullptr);
  hPulseBar30Multiplicity->SetStats(false);
  for (int count : bar30Multiplicities) hPulseBar30Multiplicity->Fill(count);

  // Zoom only the control-bar display, preserving the original fine bins and flow counts.
  // Native/projected ECal axes are independent. Raw/aligned x share occupied ranges.
  for (int projected = 0; projected < 2; ++projected) {
    int firstX[3], lastX[3] = {}, firstY[3], lastY[3] = {};
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
      TH2D* hist = hPulsePositionCorrelation[0][projected][coordinate];
      firstX[coordinate] = hist->GetNbinsX() + 1;
      firstY[coordinate] = hist->GetNbinsY() + 1;
      for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
          if (hist->GetBinContent(ix, iy) <= 0.0) continue;
          firstX[coordinate] = std::min(firstX[coordinate], ix);
          lastX[coordinate] = std::max(lastX[coordinate], ix);
          firstY[coordinate] = std::min(firstY[coordinate], iy);
          lastY[coordinate] = std::max(lastY[coordinate], iy);
        }
      }
    }
    firstX[0] = firstX[1] = std::min(firstX[0], firstX[1]);
    lastX[0] = lastX[1] = std::max(lastX[0], lastX[1]);
    firstY[0] = firstY[1] = std::min(firstY[0], firstY[1]);
    lastY[0] = lastY[1] = std::max(lastY[0], lastY[1]);
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
      TH2D* hist = hPulsePositionCorrelation[0][projected][coordinate];
      auto zoomAxis = [&](TAxis* axis, int first, int last) {
        if (last < first) return; // Empty plots retain the configured detector range.
        const int paddingBins = int(std::min(double(axis->GetNbins()),
            std::max(1.0, std::ceil(ControlPositionPadding / axis->GetBinWidth(1)))));
        axis->SetRange(std::max(1, first-paddingBins), std::min(axis->GetNbins(), last+paddingBins));
        axis->SetNdivisions(505, false); // Keep closely spaced decimal tick labels readable.
      };
      zoomAxis(hist->GetXaxis(), firstX[coordinate], lastX[coordinate]);
      zoomAxis(hist->GetYaxis(), firstY[coordinate], lastY[coordinate]);
    }
  }

  // Six interactive canvases for the control bar, detector, and each layer.
  // The detector population includes each accepted pulse once, independent of bar choices.
  for (int selected = 0; selected < 4; ++selected) {
    for (int projected = 0; projected < 2; ++projected) {
      double xColorMax = std::max(1.0, std::max(hPulsePositionCorrelation[selected][projected][0]->GetMaximum(), hPulsePositionCorrelation[selected][projected][1]->GetMaximum()));
      double yColorMax = std::max(1.0, hPulsePositionCorrelation[selected][projected][2]->GetMaximum());
      // Match L1/L2 color scales as well as axes, for direct count comparisons.
      if (selected >= 2) {
        for (int layer = 2; layer < 4; ++layer) {
          for (int coordinate = 0; coordinate < 2; ++coordinate)
            xColorMax = std::max(xColorMax, hPulsePositionCorrelation[layer][projected][coordinate]->GetMaximum());
          yColorMax = std::max(yColorMax, hPulsePositionCorrelation[layer][projected][2]->GetMaximum());
        }
      }
      for (int coordinate = 0; coordinate < 3; ++coordinate) {
        TH2D* histogram = hPulsePositionCorrelation[selected][projected][coordinate];
        TCanvas* canvas = new TCanvas(TString::Format("cPulseCorrelation%s%s%s", roles[selected], ecalNames[projected], coordinateNames[coordinate]), histogram->GetTitle(), 1200, 742);
        canvas->SetCanvasSize(1200, 742);
        canvas->SetTopMargin(0.10);
        canvas->SetLeftMargin(0.12);
        canvas->SetRightMargin(0.15);
        histogram->SetMinimum(0.0);
        if (coordinate < 2) histogram->SetMaximum(xColorMax);
        else if (selected >= 2) histogram->SetMaximum(yColorMax);
        histogram->Draw("COLZ");
        if (projected && coordinate < 2) {
          // Mean CDet x per projected-ECal-x bin; errors are errors on the mean.
          // Use regular y bins only, excluding positions outside the histogram.
          TProfile* profile = histogram->ProfileX(TString::Format("%s_profile", histogram->GetName()),
                                                  1, histogram->GetNbinsY());
          profile->SetDirectory(nullptr);
          profile->SetStats(false);
          profile->SetMarkerStyle(20);
          profile->SetMarkerSize(0.65);
          profile->SetMarkerColor(kBlack);
          profile->SetLineColor(kBlack);
          profile->DrawCopy("E1 SAME"); // Canvas owns the copy, including on repeated calls.
          delete profile;
        }
        canvas->Modified();
        canvas->Update();
        PulseStudySaveCanvas(canvas, TString::Format("%s/%s_%s_vs_ecal_%s.pdf", outputDirectory.Data(), filePrefixes[selected].Data(), coordinateFiles[coordinate], projected ? "projected" : "native"));
      }
    }
  }
  double differenceMaximum = 1.0;
  for (int layer = 0; layer < 2; ++layer)
    for (int coordinate = 0; coordinate < 2; ++coordinate)
      differenceMaximum = std::max(differenceMaximum, hPulseLayerXDifference[layer][coordinate]->GetMaximum());
  for (int layer = 0; layer < 2; ++layer) {
    for (int coordinate = 0; coordinate < 2; ++coordinate) {
      TH1D* histogram = hPulseLayerXDifference[layer][coordinate];
      TCanvas* canvas = new TCanvas(TString::Format("cPulseLayer%d%sDifference", layer+1, coordinateNames[coordinate]), histogram->GetTitle(), 1200, 742);
      canvas->SetCanvasSize(1200, 742);
      canvas->SetTopMargin(0.10);
      canvas->SetLeftMargin(0.12);
      canvas->SetRightMargin(0.05);
      histogram->SetLineColor(layer == 0 ? kRed+1 : kBlue+1);
      histogram->SetLineWidth(2);
      histogram->SetMinimum(0.0);
      histogram->SetMaximum(1.15*differenceMaximum);
      histogram->Draw("HIST");
      canvas->Modified();
      canvas->Update();
      PulseStudySaveCanvas(canvas, TString::Format("%s/layer%d_%s_minus_ecal_projected.pdf", outputDirectory.Data(), layer+1, coordinateFiles[coordinate]));
    }
  }
  const double bar30Maximum = 1.15*std::max(1.0, std::max(hPulseBar30XDifference[0]->GetMaximum(), hPulseBar30XDifference[1]->GetMaximum()));
  for (int coordinate = 0; coordinate < 2; ++coordinate) {
    TH1D* histogram = hPulseBar30XDifference[coordinate];
    TCanvas* canvas = new TCanvas(TString::Format("cPulseBar30%sDifference", coordinateNames[coordinate]), histogram->GetTitle(), 1200, 742);
    canvas->SetCanvasSize(1200, 742);
    canvas->SetTopMargin(0.10); canvas->SetLeftMargin(0.12); canvas->SetRightMargin(0.05);
    histogram->SetLineColor(kRed+1); histogram->SetLineWidth(2);
    histogram->SetMinimum(0.0); histogram->SetMaximum(bar30Maximum);
    histogram->Draw("HIST");
    canvas->Modified(); canvas->Update();
    PulseStudySaveCanvas(canvas, TString::Format("%s/bar030_%s_minus_ecal_projected.pdf", outputDirectory.Data(), coordinateFiles[coordinate]));
  }
  const double layerColorMax = std::max(1.0, std::max(hPulseLayerX2VsX1[0]->GetMaximum(), hPulseLayerX2VsX1[1]->GetMaximum()));
  const double layerDifferenceMaximum = 1.15*std::max(1.0, std::max(hPulseLayerX2MinusX1[0]->GetMaximum(), hPulseLayerX2MinusX1[1]->GetMaximum()));
  for (int coordinate = 0; coordinate < 2; ++coordinate) {
    TCanvas* correlation = new TCanvas(TString::Format("cPulseLayers%sX2VsX1", coordinateNames[coordinate]), hPulseLayerX2VsX1[coordinate]->GetTitle(), 1200, 742);
    correlation->SetCanvasSize(1200, 742);
    correlation->SetTopMargin(0.10); correlation->SetLeftMargin(0.12); correlation->SetRightMargin(0.15);
    hPulseLayerX2VsX1[coordinate]->SetMinimum(0.0);
    hPulseLayerX2VsX1[coordinate]->SetMaximum(layerColorMax);
    hPulseLayerX2VsX1[coordinate]->Draw("COLZ");
    correlation->Modified(); correlation->Update();
    PulseStudySaveCanvas(correlation, TString::Format("%s/bar030_signal_%s_x2_vs_x1.pdf", outputDirectory.Data(), coordinate ? "corrected" : "raw"));
    TCanvas* difference = new TCanvas(TString::Format("cPulseLayers%sX2MinusX1", coordinateNames[coordinate]), hPulseLayerX2MinusX1[coordinate]->GetTitle(), 1200, 742);
    difference->SetCanvasSize(1200, 742);
    difference->SetTopMargin(0.10); difference->SetLeftMargin(0.12); difference->SetRightMargin(0.05);
    hPulseLayerX2MinusX1[coordinate]->SetMinimum(0.0);
    hPulseLayerX2MinusX1[coordinate]->SetMaximum(layerDifferenceMaximum);
    hPulseLayerX2MinusX1[coordinate]->SetLineWidth(2);
    hPulseLayerX2MinusX1[coordinate]->SetLineColor(kBlue+1);
    hPulseLayerX2MinusX1[coordinate]->Draw("HIST");
    difference->Modified(); difference->Update();
    PulseStudySaveCanvas(difference, TString::Format("%s/bar030_signal_%s_x2_minus_x1.pdf", outputDirectory.Data(), coordinate ? "corrected" : "raw"));
  }
  TCanvas* multiplicity = new TCanvas("cPulseBar30Multiplicity", hPulseBar30Multiplicity->GetTitle(), 1200, 742);
  multiplicity->SetCanvasSize(1200, 742);
  multiplicity->SetTopMargin(0.10); multiplicity->SetLeftMargin(0.12); multiplicity->SetRightMargin(0.05);
  if (maxBar30Hits <= 20) {
    for (int count = 1; count <= maxBar30Hits; ++count)
      hPulseBar30Multiplicity->GetXaxis()->SetBinLabel(count, TString::Format("%d", count));
  }
  hPulseBar30Multiplicity->SetMinimum(0.0);
  hPulseBar30Multiplicity->SetLineWidth(2);
  hPulseBar30Multiplicity->SetLineColor(kBlue+1);
  hPulseBar30Multiplicity->Draw("HIST");
  multiplicity->Modified(); multiplicity->Update();
  PulseStudySaveCanvas(multiplicity, outputDirectory + "/bar030_hits_per_event.pdf");
  TH1D* positionHistograms[] = {hPulsePositionCDetX, hPulsePositionECalX};
  const char* positionCanvasNames[] = {"cPulsePositionCDetX", "cPulsePositionECalX"};
  const char* positionFiles[] = {"detector_cdet_x.pdf", "detector_ecal_x_native.pdf"};
  for (int position = 0; position < 2; ++position) {
    TH1D* histogram = positionHistograms[position];
    TCanvas* canvas = new TCanvas(positionCanvasNames[position], histogram->GetTitle(), 1200, 742);
    canvas->SetCanvasSize(1200, 742);
    canvas->SetTopMargin(0.10);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.05);
    histogram->SetLineColor(kBlue+1);
    histogram->SetLineWidth(2);
    histogram->SetMinimum(0.0);
    histogram->SetMaximum(1.15*std::max(1.0, histogram->GetMaximum()));
    histogram->Draw("HIST");
    canvas->Modified();
    canvas->Update();
    PulseStudySaveCanvas(canvas, outputDirectory + "/" + positionFiles[position]);
  }
  TCanvas* xOverlay = new TCanvas("cPulsePositionNativeXOverlay", "Channel CDet x and native ECal x", 1200, 742);
  xOverlay->SetCanvasSize(1200, 742); xOverlay->SetTopMargin(.10); xOverlay->SetLeftMargin(.12); xOverlay->SetRightMargin(.05);
  TH1* cdetX = hPulsePositionCDetX->DrawCopy("HIST", "_overlay");
  TH1* ecalX = hPulsePositionECalX->DrawCopy("HIST SAME", "_overlay");
  cdetX->SetTitle(TString::Format("%s: channel CDet x and native ECal x;Position x (m);Counts (CDet pulses / ECal events)", PulseStudyRunLabel(config, RunNumber).Data()));
  cdetX->SetMaximum(1.45*std::max(1.,std::max(hPulsePositionCDetX->GetBinContent(hPulsePositionCDetX->GetMaximumBin()), hPulsePositionECalX->GetBinContent(hPulsePositionECalX->GetMaximumBin()))));
  cdetX->SetLineColor(kBlue+1); ecalX->SetLineColor(kRed+1);
  TLegend* xLegend = new TLegend(.55,.76,.94,.88);
  xLegend->SetBit(kCanDelete); xLegend->SetBorderSize(0); xLegend->SetFillColor(kWhite);
  xLegend->AddEntry(cdetX,"Channel CDet x (per pulse)","l");
  xLegend->AddEntry(ecalX,"Native ECal x (per event)","l"); xLegend->Draw();
  PulseStudySaveCanvas(xOverlay, outputDirectory + "/detector_cdet_x_and_ecal_x_native.pdf");

  std::cout << "[CDet position correlations] Selected " << filledPulses
            << " detector pulses; thirty-eight canvases remain open (including the native-x overlay, 4 layer-combination plots and bar-30 multiplicity). ECal correlations use ECal horizontally; layer plots use x1 horizontally.\n"
            << "  Bar-30 signal events: " << bar30SignalEvents << "; channel/aligned layer combinations: "
            << hPulseLayerX2VsX1[0]->GetEntries() << "/" << hPulseLayerX2VsX1[1]->GetEntries() << '\n'
            << "  Multiplicity includes " << bar30Multiplicities.size() << " nonempty bar-30 events before secondary pulse cuts. PDF output: " << outputDirectory << '\n';
}

// Detector occupancy and one ECal-selected saved event. No secondary rx/ellipse cut.
// eventIndex is the saved-vector index, NOT the replay event number; -1 chooses the first.
// Single-event maps/counts ALWAYS require the configured position dt coincidence window.
void plotPulseStudyEventMaps(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                            const char* outputTag = nullptr, Long64_t eventIndex = -1)
{
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 || config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet event maps] Cannot read configuration with config.version: 1.\n"; return;
  }
  const bool energyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool energyMaxCut = config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0;
  const bool timeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const bool totCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  const bool totMaxCut = config.GetValue("pulse_study.apply_tot_max", 0) != 0;
  const bool dtCut = config.GetValue("pulse_study.position_apply_dt_cut", 0) != 0;
  const bool correctedX = config.GetValue("pulse_study.position_use_corrected_x", 1) != 0;
  const double emin = config.GetValue("pulse_study.ecal_energy_min", NAN), emax = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double tmin = config.GetValue("pulse_study.ecal_time_min", NAN), tmax = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double lemin = config.GetValue("analysis.le_min", .02), lemax = config.GetValue("analysis.le_max", 60.);
  const double totmin = config.GetValue("analysis.tot_min", 4.), totmax = config.GetValue("analysis.tot_max", 30.);
  const double dtmin = config.GetValue("pulse_study.position_dt_min_ns", config.GetValue("pulse_study.fit_signal_min", -36.));
  const double dtmax = config.GetValue("pulse_study.position_dt_max_ns", config.GetValue("pulse_study.fit_signal_max", -20.));
  const double minimum = config.GetValue("pulse_study.position_min_m", -1.6), maximum = config.GetValue("pulse_study.position_max_m", 1.6);
  const double width = config.GetValue("pulse_study.position_bin_width_m", .005);
  auto validWindow = [](double a, double b) { return std::isfinite(a) && std::isfinite(b) && a < b; };
  if (!validWindow(dtmin, dtmax) || !validWindow(lemin, lemax) ||
      (energyCut && (!std::isfinite(emin) || (energyMaxCut && !validWindow(emin, emax)))) ||
      (timeCut && !validWindow(tmin, tmax)) ||
      (totCut && (!std::isfinite(totmin) || totmin < 0 || (totMaxCut && !validWindow(totmin, totmax)))) ||
      !validWindow(minimum, maximum) || !std::isfinite(width) || width <= 0 ||
      !std::isfinite((maximum-minimum)/width) || std::ceil((maximum-minimum)/width) > std::numeric_limits<int>::max()-2) {
    std::cerr << "[CDet event maps] Invalid cuts or position binning.\n"; return;
  }
  const size_t nEvents = vPulseID.size();
  if (!nEvents || eventIndex < -1 || (eventIndex >= 0 && size_t(eventIndex) >= nEvents)) {
    std::cerr << "[CDet event maps] Run the main macro first; choose saved event index 0..N-1, or -1 for automatic selection.\n"; return;
  }
  if (vECalE.size() != nEvents || vECalAdcTime.size() != nEvents) {
    std::cerr << "[CDet event maps] ECal vectors are not aligned; reload the main macro.\n"; return;
  }
  auto& positions = correctedX ? vPulseXCorr : vPulseX;
  for (const auto* values : {&vPulseDt, &vPulseLe, &vPulseToT, &positions, &vPulseY}) {
    if (values->size() != nEvents) { std::cerr << "[CDet event maps] Missing pulse vectors; reload.\n"; return; }
    for (size_t event = 0; event < nEvents; ++event)
      if ((*values)[event].size() != vPulseID[event].size()) {
        std::cerr << "[CDet event maps] Pulse vectors are not event/hit aligned; reload.\n"; return;
      }
  }
  auto goodECal = [&](size_t event) {
    const double energy = vECalE[event], time = vECalAdcTime[event];
    return std::isfinite(energy) && std::isfinite(time) &&
        (!energyCut || (energy > emin && (!energyMaxCut || energy < emax))) &&
        (!timeCut || (time > tmin && time < tmax));
  };
  if (eventIndex == -1) {
    for (size_t event = 0; event < nEvents; ++event) if (goodECal(event)) { eventIndex = event; break; }
  }
  if (eventIndex < 0 || !goodECal(size_t(eventIndex))) {
    std::cerr << "[CDet event maps] No ECal-selected event at the requested index (or in the sample).\n"; return;
  }
  // Full DB geometry also supplies bars with no saved hits. Staggered y centers
  // require polygon cells: every cell spans its own bar center +/- 25 cm.
  TString geometryFile = config.GetValue("pulse_study.map_database_file", "../../DB/db_earm.cdet.dat");
  gSystem->ExpandPathName(geometryFile);
  if (!geometryFile.BeginsWith("/")) geometryFile = TString(gSystem->DirName(configFile)) + "/" + geometryFile;
  std::ifstream geometry(geometryFile.Data());
  if (!geometry) { std::cerr << "[CDet event maps] Cannot read geometry DB: " << geometryFile << '\n'; return; }
  std::map<std::string, std::vector<double>> db = {{"xpos", {}}, {"ypos", {}},
      {"selection.x_scale_l1", {}}, {"selection.x_scale_l2", {}},
      {"selection.x_alignment_l1", {}}, {"selection.x_alignment_l2", {}}};
  std::map<std::string, int> blocks;
  std::string line, reading;
  while (std::getline(geometry, line)) {
    line = line.substr(0, line.find('#'));
    const auto equals = line.find('=');
    if (equals != std::string::npos) {
      std::istringstream keyInput(line.substr(0, equals));
      std::string key; keyInput >> key;
      reading.clear();
      const std::string prefix = "earm.cdet.";
      if (key.compare(0, prefix.size(), prefix) == 0 && db.count(key.substr(prefix.size()))) {
        reading = key.substr(prefix.size()); ++blocks[reading];
      }
      line = line.substr(equals+1);
    }
    if (line.find('[') != std::string::npos) reading.clear();
    if (reading.empty()) continue;
    std::istringstream numbers(line);
    double value; while (numbers >> value) db[reading].push_back(value);
  }
  for (const auto& item : db) {
    const bool positions = item.first == "xpos" || item.first == "ypos";
    if (!positions && !correctedX) continue;
    if (blocks[item.first] != 1 || (positions ? item.second.size() < size_t(numBars*pixelsPerBar) : item.second.size() != 1)) {
      std::cerr << "[CDet event maps] Missing/ambiguous geometry DB key: " << item.first << '\n'; return;
    }
  }
  struct MapBar { double minX, maxX, lowX, highX, y; };
  std::vector<MapBar> mapBars(numBars);
  double yMinimum = std::numeric_limits<double>::infinity(), yMaximum = -yMinimum;
  for (int bar = 0; bar < numBars; ++bar) {
    MapBar& b = mapBars[bar];
    b.minX = std::numeric_limits<double>::infinity(); b.maxX = -b.minX; b.y = NAN;
    const int layer = bar / 84;
    const double scale = correctedX ? db[layer ? "selection.x_scale_l2" : "selection.x_scale_l1"][0] : 1.;
    const double alignment = correctedX ? db[layer ? "selection.x_alignment_l2" : "selection.x_alignment_l1"][0] : 0.;
    if (!std::isfinite(scale) || scale <= 0 || !std::isfinite(alignment)) {
      std::cerr << "[CDet event maps] Invalid x alignment in geometry DB.\n"; return;
    }
    int active = 0;
    for (int id = bar*pixelsPerBar; id < (bar+1)*pixelsPerBar; ++id) {
      const double rawX = db["xpos"][id], y = db["ypos"][id];
      if (!std::isfinite(rawX) || !std::isfinite(y) || std::abs(rawX) >= 900. || std::abs(y) >= 900.) continue;
      if (std::isfinite(b.y) && std::abs(b.y-y) > 1.e-9) {
        std::cerr << "[CDet event maps] Inconsistent DB y centers within bar " << bar << '\n'; return;
      }
      const double x = scale*rawX-alignment;
      b.minX = std::min(b.minX,x); b.maxX = std::max(b.maxX,x); b.y = y; ++active;
    }
    if (active < 2 || b.maxX <= b.minX) { std::cerr << "[CDet event maps] Missing geometry for bar " << bar << '\n'; return; }
    const double halfPitch = .5*(b.maxX-b.minX)/(active-1);
    b.lowX = b.minX-halfPitch; b.highX = b.maxX+halfPitch;
    yMinimum = std::min(yMinimum,b.y-.25); yMaximum = std::max(yMaximum,b.y+.25);
  }
  // Neighboring bars on the same side meet halfway between outer pixel centers.
  // This avoids overlapping polygon cells at a rounded DB coordinate boundary.
  for (int side = 0; side < 4; ++side) {
    std::vector<int> order;
    for (int bar = side*42; bar < (side+1)*42; ++bar) order.push_back(bar);
    std::sort(order.begin(),order.end(),[&](int a,int b){return mapBars[a].minX < mapBars[b].minX;});
    for (size_t i = 1; i < order.size(); ++i) {
      MapBar& low = mapBars[order[i-1]]; MapBar& high = mapBars[order[i]];
      const double boundary = .5*(low.maxX+high.minX);
      low.highX = boundary; high.lowX = boundary;
    }
  }
  TString tag = outputTag ? outputTag : config.GetValue("pulse_study.output_tag", "");
  const int run = pulseStudyRunNumber > 0 ? pulseStudyRunNumber : config.GetValue("analysis.run_number", 0);
  if (tag.IsNull()) tag = TString::Format("run%d", run);
  const TString allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) if (allowed.First(tag[i]) == kNPOS) {
    std::cerr << "[CDet event maps] Output tag may contain only letters, numbers, underscores and hyphens.\n"; return;
  }
  const TString directory = TString::Format("%s/PusleStudy/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(directory)) gSystem->mkdir(directory, true);
  FileStat_t info;
  if (gSystem->GetPathInfo(directory, info) != 0 || !R_ISDIR(info.fMode) || gSystem->AccessPathName(directory, kWritePermission)) {
    std::cerr << "[CDet event maps] Cannot write output directory.\n"; return;
  }
  std::vector<TString> cuts = {
      energyCut ? PulseStudyEnergyCutLabel(emin, emax, energyMaxCut) : TString("No secondary ECal energy window."),
      timeCut ? TString::Format("ECal time: %.12g < t_ECal < %.12g ns", tmin, tmax) : TString("No secondary ECal time window."),
      TString::Format("Finite raw LE: %.12g <= LE <= %.12g ns", lemin, lemax),
      totCut ? PulseStudyToTCutLabel(totmin, totmax, totMaxCut) : TString("No secondary ToT window; finite ToT > 0 required."),
      "NO secondary rx, ry or ellipse cuts; original loading cuts still apply.",
      dtCut ? TString::Format("Aggregate maps: %.12g <= stored dt <= %.12g ns", dtmin, dtmax) : TString("Aggregate maps: finite dt, no secondary dt window."),
      TString::Format("Single-event maps AND bar counts ALWAYS require %.12g <= stored dt <= %.12g ns, before y correction.", dtmin, dtmax)
  };
  if (!PulseStudyWriteCutReport(directory + TString::Format("/plotPulseStudyEventMaps_event%lld_cuts.txt", eventIndex), "plotPulseStudyEventMaps", configFile, config, cuts, {
      TString::Format("Saved-vector index %lld (NOT replay event number); selected using ECal cuts only, without requiring bar 30 or any coincident CDet pulse.", eventIndex),
      TString::Format("ECal energy %.12g GeV; ADC time %.12g ns", vECalE[eventIndex], vECalAdcTime[eventIndex]),
      TString::Format("Horizontal CDet y = DB bar center (not reconstructed hit position along bar); vertical CDet x = %s.", correctedX ? "pulse.x_corr" : "pulse.x"),
      TString::Format("Map display: CDet x [%.12g,%.12g] m, requested bin width %.12g m; CDet y [%.12g,%.12g] m. Polygon cells span each bar DB y center +/- 0.25 m; physical gaps are unfilled.", minimum, maximum, width, yMinimum, yMaximum),
      "Geometry DB: " + geometryFile,
      "Dashed guides lie at each bar midpoint in x (midpoint of its outer valid pixel centers), extending only over that bar y center +/- 0.25 m. Aligned-x mode also aligns the geometry guides.",
      "One count per accepted pulse, not unique pixel. Bar counts include every bar, including zero-hit bars; coordinates do not affect counts.",
      "Both layers are displayed separately. Position maps require finite coordinates; overflow remains outside displayed axes. Counts use no display-bound cut."
  })) return;

  for (int layer = 0; layer < 2; ++layer) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseDetectorMapL%d", layer+1));
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPulseEventMapL%d", layer+1));
    delete hPulseDetectorMap[layer]; delete hPulseEventMap[layer];
  }
  delete gROOT->GetListOfCanvases()->FindObject("cPulseEventBarHits"); delete hPulseEventBarHits;
  const int nCDetXBins = int(std::ceil((maximum-minimum)/width));
  const double xStep = (maximum-minimum)/nCDetXBins;
  const TString runLabel = PulseStudyRunLabel(config, run);
  const char* xLabel = correctedX ? "Aligned CDet x (m)" : "Channel CDet x (m)";
  for (int layer = 0; layer < 2; ++layer) {
    hPulseDetectorMap[layer] = new TH2Poly(TString::Format("hPulseDetectorMapL%d", layer+1),
        TString::Format("%s: L%d pulse occupancy;CDet y, DB bar center (m);%s;Pulses", runLabel.Data(), layer+1, xLabel), yMinimum, yMaximum, minimum, maximum);
    hPulseEventMap[layer] = new TH2Poly(TString::Format("hPulseEventMapL%d", layer+1),
        TString::Format("%s, saved event %lld: L%d coincident pulses;CDet y, DB bar center (m);%s;Pulses", runLabel.Data(), eventIndex, layer+1, xLabel), yMinimum, yMaximum, minimum, maximum);
    hPulseDetectorMap[layer]->SetDirectory(nullptr); hPulseEventMap[layer]->SetDirectory(nullptr);
    for (int bar = layer*84; bar < (layer+1)*84; ++bar) {
      const auto& b = mapBars[bar];
      const double low = std::max(minimum,b.lowX), high = std::min(maximum,b.highX);
      if (low >= high) continue;
      const int first = std::max(0,int(std::floor((low-minimum)/xStep)));
      const int last = std::min(nCDetXBins-1,int(std::ceil((high-minimum)/xStep))-1);
      for (int ix = first; ix <= last; ++ix) {
        const double x1 = std::max(low,minimum+ix*xStep), x2 = std::min(high,minimum+(ix+1)*xStep);
        if (x1 >= x2) continue;
        hPulseDetectorMap[layer]->AddBin(b.y-.25,x1,b.y+.25,x2);
        hPulseEventMap[layer]->AddBin(b.y-.25,x1,b.y+.25,x2);
      }
    }
  }
  hPulseEventBarHits = new TH1D("hPulseEventBarHits", TString::Format("%s, saved event %lld: coincident pulses;CDet bar number;Pulses in event", runLabel.Data(), eventIndex), numBars, -.5, numBars-.5);
  hPulseEventBarHits->SetDirectory(nullptr);
  for (size_t event = 0; event < nEvents; ++event) {
    if (!goodECal(event)) continue;
    for (size_t hit = 0; hit < vPulseID[event].size(); ++hit) {
      const double id = vPulseID[event][hit], dt = vPulseDt[event][hit];
      const double le = vPulseLe[event][hit], tot = vPulseToT[event][hit];
      if (!std::isfinite(id) || id < 0 || id >= numBars*pixelsPerBar || std::floor(id) != id || !std::isfinite(dt) ||
          !std::isfinite(le) || le < lemin || le > lemax || !std::isfinite(tot) || tot <= 0 ||
          (totCut && (tot < totmin || (totMaxCut && tot > totmax)))) continue;
      const int layer = int(id)/1344, bar = int(id)/pixelsPerBar;
      const bool coincidence = dt >= dtmin && dt <= dtmax;
      const double x = positions[event][hit], y = vPulseY[event][hit];
      if (std::isfinite(x) && std::isfinite(y)) {
        if (!dtCut || coincidence) hPulseDetectorMap[layer]->Fill(y, x);
        if (Long64_t(event) == eventIndex && coincidence) hPulseEventMap[layer]->Fill(y, x);
      }
      if (Long64_t(event) == eventIndex && coincidence) hPulseEventBarHits->Fill(bar);
    }
  }
  for (int single = 0; single < 2; ++single) {
    TH2Poly** histograms = single ? hPulseEventMap : hPulseDetectorMap;
    const double colorMax = std::max(1., std::max(histograms[0]->GetMaximum(), histograms[1]->GetMaximum()));
    for (int layer = 0; layer < 2; ++layer) {
      TH2Poly* hist = histograms[layer];
      TCanvas* canvas = new TCanvas(TString::Format(single ? "cPulseEventMapL%d" : "cPulseDetectorMapL%d", layer+1), hist->GetTitle(), 1200, 742);
      canvas->SetCanvasSize(1200, 742); canvas->SetLeftMargin(.12); canvas->SetRightMargin(.15); canvas->SetTopMargin(.10);
      hist->SetStats(false); hist->SetMinimum(0); hist->SetMaximum(colorMax); hist->Draw("COLZ0"); // Leave empty geometry cells white.
      for (int bar = layer*84; bar < (layer+1)*84; ++bar) {
        const auto& b = mapBars[bar];
        const double x = .5*(b.minX+b.maxX);
        if (x < minimum || x > maximum) continue;
        TLine* guide = new TLine(b.y-.25,x,b.y+.25,x);
        guide->SetUniqueID(bar+1); guide->SetBit(kCanDelete);
        guide->SetLineStyle(2); guide->SetLineColor(kGray+2); guide->SetLineWidth(1); guide->Draw();
      }
      if (single) {
        // Outline occupied bins so individual event hits remain visible at detector scale.
        TGraph* locations = new TGraph();
        locations->SetName(TString::Format("gPulseEventLocationsL%d", layer+1));
        locations->SetBit(kCanDelete); locations->SetMarkerStyle(24);
        locations->SetMarkerSize(.9); locations->SetMarkerColor(kBlack);
        for (auto* object : *hist->GetBins()) {
          auto* bin = static_cast<TH2PolyBin*>(object);
          if (hist->GetBinContent(bin->GetBinNumber()) > 0)
            locations->SetPoint(locations->GetN(), .5*(bin->GetXMin()+bin->GetXMax()), .5*(bin->GetYMin()+bin->GetYMax()));
        }
        if (locations->GetN()) locations->Draw("P SAME");
        else delete locations;
      }
      PulseStudySaveCanvas(canvas, directory + (single ? TString::Format("/event%lld_layer%d_xy.pdf", eventIndex, layer+1) : TString::Format("/detector_layer%d_xy.pdf", layer+1)));
    }
  }
  TCanvas* counts = new TCanvas("cPulseEventBarHits", hPulseEventBarHits->GetTitle(), 1200, 742);
  counts->SetCanvasSize(1200, 742); counts->SetLeftMargin(.10); counts->SetTopMargin(.10);
  hPulseEventBarHits->SetStats(false); hPulseEventBarHits->SetMinimum(0);
  hPulseEventBarHits->SetMaximum(std::max(1., 1.2*hPulseEventBarHits->GetMaximum()));
  hPulseEventBarHits->SetLineColor(kBlack); hPulseEventBarHits->SetLineWidth(2); hPulseEventBarHits->Draw("HIST");
  TLine* layerBoundary = new TLine(83.5, 0., 83.5, hPulseEventBarHits->GetMaximum());
  layerBoundary->SetLineStyle(2); layerBoundary->SetBit(kCanDelete); layerBoundary->Draw();
  TLatex label; label.SetNDC(); label.SetTextSize(.028);
  label.DrawLatex(.25, .84, "Layer 1 (0-83)"); label.DrawLatex(.65, .84, "Layer 2 (84-167)");
  PulseStudySaveCanvas(counts, directory + TString::Format("/event%lld_hits_vs_bar.pdf", eventIndex));
  std::cout << "[CDet event maps] Saved event index " << eventIndex << " (not replay event number): " << hPulseEventBarHits->GetEntries()
            << " coincident pulses. Five canvases remain open. No secondary rx/ellipse cuts; loading cuts remain.\n  PDF output: " << directory << '\n';
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
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double controlFitMin = config.GetValue("pulse_study.control_fit_min", -40.0);
  const double controlFitMax = config.GetValue("pulse_study.control_fit_max", 0.0);
  const double FitMin = config.GetValue("pulse_study.fit_min", -50.0);
  const double FitMax = config.GetValue("pulse_study.fit_max", 0.0);
  const double DtMin = config.GetValue("display.ecal_cdet_dt_min", -50.0);
  const double DtMax = config.GetValue("display.ecal_cdet_dt_max", 50.0);
  const double BinWidth = config.GetValue("pulse_study.dt_bin_width_ns", 1.0);
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", -10);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", 10);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double xRes = config.GetValue("analysis.x_difference_max", 0.08);

  if (!std::isfinite(DtMin) || !std::isfinite(DtMax) || DtMin >= DtMax ||
      !std::isfinite(BinWidth) || BinWidth <= 0.0 ||
      !std::isfinite((DtMax - DtMin) / BinWidth) ||
      std::ceil((DtMax - DtMin) / BinWidth) > std::numeric_limits<int>::max() - 2) {
    std::cerr << "[CDet pulse offsets] Invalid dt display range or bin width.\n";
    return;
  }
  const int nDtBins = int(std::ceil((DtMax - DtMin) / BinWidth));
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
    hBarDt[pixel] = new TH1D(TString::Format("hBar%dPixel%dDt", myBar, pixel), TString::Format("Bar %d pixel %d;dt [ns];Pulses", myBar, pixel), nDtBins, DtMin, DtMax);
    hBarDt[pixel]->SetDirectory(nullptr);
    hControlBarDt[pixel] = new TH1D(TString::Format("hControlBar%dPixel%dDt", controlBar, pixel), TString::Format("Control bar %d pixel %d;dt [ns];Pulses", controlBar, pixel), nDtBins, DtMin, DtMax);
    hControlBarDt[pixel]->SetDirectory(nullptr);
  }

  Int_t nevents = vECalAdcTime.size();
  for (Int_t ev = 0; ev < nevents; ev++){
    Int_t nhits = vPulseID[ev].size();
    double ECalT = vECalAdcTime[ev];
    if (!std::isfinite(ECalT) || ((ECalT < ECalMinT || ECalT > ECalMaxT))) continue; //apply ECal time cut
    for (Int_t hit = 0; hit < nhits; hit++){
      Int_t bar = Int_t(vPulseID[ev][hit])/pixelsPerBar; //determine which bar
      Int_t pixel = Int_t(vPulseID[ev][hit])%pixelsPerBar; //convert pixelID into pixel within bar
      double rx = vPulseXResiduals[ev][hit];
      if (!std::isfinite(rx) || std::abs(rx) > xRes || !std::isfinite(vPulseDt[ev][hit])) continue; //apply x residual cut
      if (bar == myBar) hBarDt[pixel]->Fill(vPulseDt[ev][hit]);
      if (bar == controlBar) hControlBarDt[pixel]->Fill(vPulseDt[ev][hit]);
    } //end hit loop
  } //end event loop

  if (gSystem->AccessPathName("PusleStudy") && gSystem->mkdir("PusleStudy", true) != 0) {
    std::cerr << "[CDet pulse offsets] Cannot create PusleStudy output directory.\n";
    return;
  }
  const int run = config.GetValue("analysis.run_number", 0);
  const char* plotMode = viewOnly ? "view" : "fit";
  TString selectedPixels = goodPixels.empty() ? "all pixels" : "";
  for (int pixel : goodPixels) selectedPixels += TString::Format("%d ", pixel);
  if (!PulseStudyWriteCutReport(TString::Format("PusleStudy/run%d_bar%03d_extractBarOffsets_%s_cuts.txt", run, myBar, plotMode), "extractBarOffsets", configFile, config, {
        TString::Format("Finite ECal time: %.12g <= t_ECal <= %.12g ns; always applied here regardless of apply_ecal_time_cut", ECalMinT, ECalMaxT),
        TString::Format("Finite rx and dt; abs(rx) <= %.12g m; study bar %d or control bar %d", xRes, myBar, controlBar),
        "No additional ECal energy, LE, ToT or ellipse cuts in this routine (loading cuts remain)."
      }, {
        TString::Format("Mode: %s; fit pixels: %s; grouping: %s", plotMode, selectedPixels.Data(), fullNINO ? "full bar" : "eight-pixel halves"),
        TString::Format("dt display [%.12g,%.12g] ns, requested bin width %.12g ns; study fit [%.12g,%.12g]; control fit [%.12g,%.12g] ns", DtMin, DtMax, BinWidth, FitMin, FitMax, controlFitMin, controlFitMax),
        "All 16 pixels displayed regardless of fit selection; unused-channel masks retained. Y correction OFF."
      })) return;
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
      if (hist->GetMaximum() > 0.0) hist->SetMaximum(1.5 * hist->GetMaximum()); // Keep peaks below the stats box.
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
  PulseStudySaveCanvas(cBarDt, studyPDF);
  PulseStudySaveCanvas(cControlDt, controlPDF);
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
    PulseStudySaveCanvas(cBarDt, studyPDF);
    PulseStudySaveCanvas(cControlDt, controlPDF);

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
                             double ECalBinWidth = NAN, double ECalDisplayMin = NAN, double ECalDisplayMax = NAN,
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

  if (std::isnan(ECalBinWidth)) ECalBinWidth = config.GetValue("pulse_study.ecal_time_bin_width_ns", 1.0);
  if (std::isnan(ECalDisplayMin)) ECalDisplayMin = config.GetValue("display.ecal_time_min", -15.0);
  if (std::isnan(ECalDisplayMax)) ECalDisplayMax = config.GetValue("display.ecal_time_max", 15.0);

  // Existing secondary-cut settings.
  // Saved vectors already have the main macro's loading cuts applied.
  const bool applyECalEnergyCut = config.GetValue("pulse_study.apply_ecal_energy_cut", 1) != 0;
  const bool applyECalTimeCut = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double ECalMinE = config.GetValue("pulse_study.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const bool applyECalEnergyMax = config.GetValue("pulse_study.apply_ecal_energy_max", 0) != 0;
  const bool applyToTMax = config.GetValue("pulse_study.apply_tot_max", 0) != 0;
  const double ECalMinT = config.GetValue("pulse_study.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("pulse_study.ecal_time_max", NAN);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_study.apply_tot_cut", 1) != 0;
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
      !std::isfinite(ECalMinE) || (applyECalEnergyMax && (!std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE))) {
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
  const double TotDisplayMin = config.GetValue("display.tot_min", 0.0);
  const double TotDisplayMax = config.GetValue("display.tot_max", 60.0);
  const double TotBinWidth = config.GetValue("pulse_study.tot_bin_width_ns", 1.0);
  if (!std::isfinite(TotDisplayMin) || !std::isfinite(TotDisplayMax) || TotDisplayMin >= TotDisplayMax ||
      !std::isfinite(TotBinWidth) || TotBinWidth <= 0.0 ||
      !std::isfinite((TotDisplayMax-TotDisplayMin)/TotBinWidth) ||
      std::ceil((TotDisplayMax-TotDisplayMin)/TotBinWidth) > std::numeric_limits<int>::max()-2 ||
      !std::isfinite(TotMin) || TotMin < 0.0 ||
      (applyToTMax && (!std::isfinite(TotMax) || TotMin >= TotMax))) {
    std::cerr << "[CDet pulse background] Invalid ToT display or cut settings.\n";
    return;
  }
  const int nCDetTotBins = int(std::ceil((TotDisplayMax - TotDisplayMin) / TotBinWidth));
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
  if (!PulseStudyWriteCutReport(outputDirectory + "/plotPulseStudyBackground_cuts.txt", "plotPulseStudyBackground", configFile, config, {
        TString::Format("ECal time: %.12g <= t_ECal <= %.12g ns (always applied here)", ECalMinT, ECalMaxT),
        TString::Format("Finite ECal energy >= %.12g GeV; upper bound %s", ECalMinE, applyECalEnergyMax ? TString::Format("<= %.12g GeV", ECalMaxE).Data() : "off"),
        TString::Format("Corrected LE (not raw): %.12g <= LE_corr <= %.12g ns", LeMin, LeMax),
        TString::Format("Finite ToT >= %.12g ns; upper bound %s", TotMin, applyToTMax ? TString::Format("<= %.12g ns", TotMax).Data() : "off"),
        "ECal energy/time and ToT windows above are always applied here, regardless of their secondary apply_*_cut toggles; apply_*_max toggles still control upper energy/ToT bounds.",
        "No secondary rx or ellipse cut. Loading rx/ry cuts still apply.",
        TString::Format("Valid integer pulse ID; bar == %d. Multiplicity counted over ALL selected pixels in this bar before the pixel filter.", controlBar),
        TString::Format("Single/multiple rx and dt-vs-rx plus ToT-vs-rx plots retain local pixel == %d; finite rx (and finite dt for dt-vs-rx). N=1 versus N>=2.", controlPixel)
      }, {
        "ECal time histogram counts saved events passing secondary ECal cuts, before this routine's CDet cuts.",
        "Other bar histograms include all selected pixels; studyBar argument is unused. Y correction OFF.",
        TString::Format("dt display [%.12g,%.12g] ns, width %.12g ns; rx display [%.12g,%.12g] m, width %.12g m", DtMin, DtMax, BinWidth, RxMin, RxMax, RxBinWidth)
      })) return;


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
  hCDetTotVsRx = new TH2D("hPulseBackgroundCDetTotVsRx", TString::Format("Bar %d Pixel %d: CDet Tot versus rx;rx (m);CDet Tot (ns);Pulses", controlBar, controlPixel), nRxBins, RxMin, RxMax, nCDetTotBins, TotDisplayMin, TotDisplayMax);
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
    bool good_ecal_energy = std::isfinite(e_ECal) && e_ECal >= ECalMinE && (!applyECalEnergyMax || e_ECal <= ECalMaxE);
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
      bool good_cdet_tot = std::isfinite(tot) && tot >= TotMin && (!applyToTMax || tot <= TotMax);
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
  PulseStudySaveCanvas(cCDetTotVsRx, TString::Format("%s/bar%03d_pixel%02d_background_cdet_tot_vs_rx.pdf", outputDirectory.Data(), controlBar, controlPixel));

  TCanvas* cECalT = new TCanvas("cPulseBackgroundECalTime", TString::Format("ECal ADC time: %.1f <= t <= %.1f ns", ECalMinT, ECalMaxT), 900, 650);
  hECalT->Draw("HIST");
  cECalT->Modified();
  cECalT->Update();
  PulseStudySaveCanvas(cECalT, TString::Format("%s/background_ecal_adc_time.pdf", outputDirectory.Data()));

  TCanvas* cDtVsRx = new TCanvas("cPulseBackgroundDtVsRx", TString::Format("Bar %d: dt versus rx", controlBar), 900, 650);
  hDtVsRx->SetStats(false);
  gPad->SetRightMargin(0.15);
  hDtVsRx->Draw("COLZ");
  cDtVsRx->Modified();
  cDtVsRx->Update();
  PulseStudySaveCanvas(cDtVsRx, TString::Format("%s/bar%03d_background_dt_vs_rx.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cDt = new TCanvas("cPulseBackgroundDt", TString::Format("Bar %d: dt", controlBar), 900, 650);
  hDt->Draw("HIST");
  cDt->Modified();
  cDt->Update();
  PulseStudySaveCanvas(cDt, TString::Format("%s/bar%03d_background_dt.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cCDetX = new TCanvas("cPulseBackgroundCDetX", TString::Format("Bar %d: channel x", controlBar), 900, 650);
  hCDetX->Draw("HIST");
  cCDetX->Modified();
  cCDetX->Update();
  PulseStudySaveCanvas(cCDetX, TString::Format("%s/bar%03d_background_x.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cCDetXCorr = new TCanvas("cPulseBackgroundCDetXCorr", TString::Format("Bar %d: aligned x", controlBar), 900, 650);
  hCDetXCorr->Draw("HIST");
  cCDetXCorr->Modified();
  cCDetXCorr->Update();
  PulseStudySaveCanvas(cCDetXCorr, TString::Format("%s/bar%03d_background_x_corr.pdf", outputDirectory.Data(), controlBar));

  TCanvas* cCDetXCorrVsX = new TCanvas("cPulseBackgroundCDetXCorrVsX", TString::Format("Bar %d: aligned x versus channel x", controlBar), 900, 650);
  hCDetXCorrVsX->SetStats(false);
  gPad->SetRightMargin(0.15);
  hCDetXCorrVsX->Draw("COLZ");
  cCDetXCorrVsX->Modified();
  cCDetXCorrVsX->Update();
  PulseStudySaveCanvas(cCDetXCorrVsX, TString::Format("%s/bar%03d_background_x_corr_vs_x.pdf", outputDirectory.Data(), controlBar));

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
    PulseStudySaveCanvas(canvas, TString::Format("%s/bar%03d_background_rx_%s.pdf", outputDirectory.Data(), controlBar, categories[category]));
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
    PulseStudySaveCanvas(canvas, TString::Format("%s/bar%03d_background_dt_vs_rx_%s.pdf", outputDirectory.Data(), controlBar, categories[category]));
  }

  std::cout << "[CDet pulse background] Filled ECal time for " << hECalT->GetEntries()
            << " selected events and CDet histograms for " << hDt->GetEntries()
            << " selected pulses in bar " << controlBar << ". Ten canvases remain open.\n"
            << "  rx split after corrected-LE/ToT cuts: " << multiplePulseEvents
            << " events with >=2 pulses, " << singlePulseEvents << " events with exactly 1 pulse.\n"
            << "  rx histograms count pulses, not events; repeated pulses on one pixel count separately.\n"
            << "  PDF output: " << outputDirectory << '\n';
}
