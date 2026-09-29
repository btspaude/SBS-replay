// First-pass, read-only comparison of pixel timing offsets across cross-target
// runs.  This macro deliberately writes diagnostic/candidate files only; it
// never activates a calibration and never overwrites CDet_calibration_dt.dat.
//
// Usage from the cdet script directory:
//   root -l
//   .L CDet_CrossTargetPixelOffsetScan.C+
//   CDet_CrossTargetPixelOffsetScan();
//
// The run list contains one run number per line.  A second, optional column is
// an authoritative configuration path.  If omitted, the macro uses
// CDet_run<run>_projection.conf in the current directory.

#include "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C"

#include <TSystem.h>
#include <TString.h>

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

namespace CDetCrossTargetPixelScan {

struct RunSpec {
  int run = -1;
  std::string config;
};

struct PixelResult {
  double offset = std::numeric_limits<double>::quiet_NaN();
  bool valid = false;
  std::string source;
};

using PixelResults = std::vector<PixelResult>;

// The inventory in CDet_CROSS_TARGET_RUN_INVENTORY.md, excluding Run 5710.
const int kDefaultRuns[] = {
  3573, 3575, 3602, 3603, 3604, 3648, 3682, 3685, 3686, 3694,
  3757, 3758, 3786, 3788, 3844, 3845, 3846, 4003, 4006,
  4344, 4345, 4400, 4722, 4735, 4926, 4929, 4985, 5066, 5068,
  5290, 5292, 5294, 5423, 5722, 5723, 5724, 5726, 5794, 5795, 5992
};

bool ReadRunList(const char *fileName, std::vector<RunSpec>& runs) {
  std::ifstream input(fileName);
  if (!input) {
    std::cerr << "[CDet offset scan] ERROR: cannot open run list '"
              << fileName << "'.\n";
    return false;
  }

  std::string line;
  while (std::getline(input, line)) {
    const std::size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;
    std::istringstream parser(line.substr(first));
    RunSpec spec;
    if (!(parser >> spec.run) || spec.run <= 0 || spec.run == 5710) continue;
    parser >> spec.config;
    if (spec.config.empty())
      spec.config = TString::Format("CDet_run%d_projection.conf", spec.run).Data();
    runs.push_back(spec);
  }
  return !runs.empty();
}

void DefaultRunList(std::vector<RunSpec>& runs) {
  runs.clear();
  for (int run : kDefaultRuns)
    runs.push_back({run, TString::Format("CDet_run%d_projection.conf", run).Data()});
}

bool ReadRun5710Results(const char *fileName, PixelResults& results) {
  results.assign(NumCDetPaddles, PixelResult());
  std::ifstream input(fileName);
  if (!input) {
    std::cerr << "[CDet offset scan] ERROR: cannot open Run 5710 result file '"
              << fileName << "'.\n";
    return false;
  }

  // Result columns are documented by the diagnostic writer:
  // pixel bar group entries manual_cut source broad_mu group_mu used_mu
  // used_mu_err used_sigma correction total_offset fit_status reason.
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    int pixel = -1, bar = 0, group = 0, entries = 0, manualCut = 0;
    std::string source;
    double broad = 0.0, groupMu = 0.0, usedMu = 0.0, usedErr = 0.0;
    double usedSigma = 0.0, correction = 0.0, totalOffset = 0.0;
    int fitStatus = -1;
    if (!(row >> pixel >> bar >> group >> entries >> manualCut >> source
              >> broad >> groupMu >> usedMu >> usedErr >> usedSigma
              >> correction >> totalOffset >> fitStatus)) continue;
    if (pixel < 0 || pixel >= NumCDetPaddles || !std::isfinite(totalOffset)) continue;

    // A nonzero retained value can come from the later Run 5710 half-bar
    // alignment.  Zero retained rows have no direct Run 5710 offset and are
    // intentionally marked unavailable for the requested comparison.
    const bool hasOffset = source != "unused" &&
                           (source != "retained_existing" || totalOffset != 0.0);
    results[pixel].offset = totalOffset;
    results[pixel].valid = hasOffset;
    results[pixel].source = source;
  }
  return true;
}

bool ReadDiagnosticResults(const char *fileName, PixelResults& results) {
  results.assign(NumCDetPaddles, PixelResult());
  std::ifstream input(fileName);
  if (!input) return false;

  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    int pixel = -1, bar = 0, group = 0, entries = 0, manualCut = 0;
    std::string source;
    double broad = 0.0, groupMu = 0.0, usedMu = 0.0, usedErr = 0.0;
    double usedSigma = 0.0, correction = 0.0, totalOffset = 0.0;
    int fitStatus = -1;
    if (!(row >> pixel >> bar >> group >> entries >> manualCut >> source
              >> broad >> groupMu >> usedMu >> usedErr >> usedSigma
              >> correction >> totalOffset >> fitStatus)) continue;
    if (pixel < 0 || pixel >= NumCDetPaddles || !std::isfinite(totalOffset)) continue;
    // The diagnostic extractor marks known unused pixels explicitly.  Other
    // rows are valid only when the fit source produced a finite candidate.
    results[pixel].offset = totalOffset;
    results[pixel].valid = source != "unused" &&
                           source != "retained_existing" &&
                           (fitStatus == 0 || source == "group_fallback" ||
                            source == "individual_broad_fallback" ||
                            source == "individual_fit" || source == "manual_2d_fit");
    results[pixel].source = source;
  }
  return true;
}

bool WriteComparisonTable(const char *fileName,
                          const PixelResults& run5710,
                          const std::vector<RunSpec>& runs,
                          const std::map<int, PixelResults>& perRun) {
  std::ofstream output(fileName);
  if (!output) {
    std::cerr << "[CDet offset scan] ERROR: cannot create output table '"
              << fileName << "'.\n";
    return false;
  }

  output << "# Read-only first-pass cross-target pixel-offset comparison\n"
         << "# Run 5710 reference file is not used to modify calibration constants\n"
         << "# Difference convention: run_offset_ns - run5710_offset_ns\n"
         << "pixel_id\trun5710_offset_ns";
  for (const RunSpec& spec : runs)
    output << "\trun_" << spec.run << "_offset_ns\trun_" << spec.run
           << "_minus_5710_ns";
  output << "\n";

  output << std::setprecision(9);
  for (int pixel = 0; pixel < NumCDetPaddles; ++pixel) {
    output << pixel << "\t";
    if (run5710[pixel].valid) output << run5710[pixel].offset;
    else output << "n/a";
    for (const RunSpec& spec : runs) {
      const PixelResults& values = perRun.at(spec.run);
      output << "\t";
      if (!values[pixel].valid) {
        output << "n/a\tn/a";
      } else {
        output << values[pixel].offset << "\t";
        if (!run5710[pixel].valid) output << "n/a";
        else output << values[pixel].offset - run5710[pixel].offset;
      }
    }
    output << "\n";
  }
  return true;
}

} // namespace CDetCrossTargetPixelScan

void CDet_CrossTargetPixelOffsetScan(
    const char *runList = "",
    const char *run5710Results =
      "CDet_pixel_timing_fit_results_run5710_reviewed_102_cut_pass.dat",
    const char *outputTable = "CDet_cross_target_pixel_offset_scan.tsv",
    const char *workDirectory = "CDet_cross_target_pixel_offset_scan_work",
    Int_t events = -1) {
  using namespace CDetCrossTargetPixelScan;

  std::vector<RunSpec> runs;
  if (runList && runList[0] != '\0') {
    if (!ReadRunList(runList, runs)) return;
  } else {
    DefaultRunList(runs);
  }

  PixelResults run5710;
  if (!ReadRun5710Results(run5710Results, run5710)) return;

  gSystem->mkdir(workDirectory, kTRUE);
  std::map<int, PixelResults> perRun;
  std::vector<RunSpec> successfulRuns;

  for (const RunSpec& spec : runs) {
    const TString tag = TString::Format("run%d", spec.run);
    const TString resultFile = TString::Format("%s/%s_results.dat",
                                                workDirectory, tag.Data());
    const TString candidateFile = TString::Format("%s/%s_candidate.dat",
                                                   workDirectory, tag.Data());
    const TString rootFile = TString::Format("%s/%s_diagnostics.root",
                                              workDirectory, tag.Data());
    const TString plotDirectory = TString::Format("%s/%s_plots",
                                                    workDirectory, tag.Data());

    std::cout << "[CDet offset scan] Processing Run " << spec.run
              << " using " << spec.config << "\n";

    // Stage 0 disables active pixel offsets.  The diagnostic extractor below
    // writes a candidate only, so no production calibration is changed.
    PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
        spec.config.c_str(), 0, events);
    if (!gLastCalibrationStageSucceeded) {
      std::cerr << "[CDet offset scan] WARNING: Run " << spec.run
                << " did not produce a valid analysis sample; skipping.\n";
      continue;
    }

    extractHierarchicalCDetPixelTimingOffsetsDiagnostic(
        1.0, -60.0, 30.0, -30.0, 10.0,
        35, 100, 0.5, 8.0, 15.0, 2.0, 1.5, 4.0, 8.0,
        1.0, 12.0, candidateFile, resultFile, rootFile,
        plotDirectory, "29,79,104,118,139,148", false, "");

    PixelResults values;
    if (!ReadDiagnosticResults(resultFile, values)) {
      std::cerr << "[CDet offset scan] WARNING: no readable result file for Run "
                << spec.run << "; skipping.\n";
      continue;
    }
    perRun[spec.run] = values;
    successfulRuns.push_back(spec);
  }

  if (successfulRuns.empty()) {
    std::cerr << "[CDet offset scan] ERROR: no runs completed successfully.\n";
    return;
  }
  if (!WriteComparisonTable(outputTable, run5710, successfulRuns, perRun)) return;

  std::cout << "[CDet offset scan] Wrote " << outputTable << "\n"
            << "[CDet offset scan] Work files and diagnostic plots are in "
            << workDirectory << "\n"
            << "[CDet offset scan] No production calibration file was activated or overwritten.\n";
}
