#include <TDatime.h>
#include <TString.h>
#include <TSystem.h>

#include "THaRun.h"
#include "THaRunBase.h"

#include <fstream>
#include <iomanip>
#include <iostream>

// Read only the CODA prestart metadata from stream 0, segment 0 for each run.
// No detector initialization or event replay is performed.
//
// Example from an analyzer session on the JLab farm:
//   .L Run_CDet_Extract_CODATimestamps.C+
//   Run_CDet_Extract_CODATimestamps(
//       3573, 6088, "/cache/halla/sbs/GEp/raw",
//       "CDet_5pass_CODA_timestamps.csv")
void Run_CDet_Extract_CODATimestamps(
    int firstRun = 3573, int lastRun = 6088,
    const char *dataDirectory = "/cache/halla/sbs/GEp/raw",
    const char *outputCSV = "CDet_5pass_CODA_timestamps.csv",
    const char *filePrefix = "gep5") {
  if (firstRun <= 0 || lastRun < firstRun || !dataDirectory ||
      !dataDirectory[0] || !outputCSV || !outputCSV[0] || !filePrefix ||
      !filePrefix[0]) {
    std::cerr << "[CDet CODA timestamps] Invalid run range, directory, "
                 "output filename, or CODA prefix.\n";
    return;
  }

  std::ofstream output(outputCSV);
  if (!output) {
    std::cerr << "[CDet CODA timestamps] Cannot create " << outputCSV << '\n';
    return;
  }
  output << "requested_run,coda_run,timestamp,unix_time,db_validity_header,"
            "filename,status\n";

  int found = 0;
  int initialized = 0;
  int missing = 0;
  int failed = 0;
  int mismatched = 0;
  for (int requestedRun = firstRun; requestedRun <= lastRun; ++requestedRun) {
    const TString filename = TString::Format(
        "%s/%s_%d.evio.0.0", dataDirectory, filePrefix, requestedRun);
    if (gSystem->AccessPathName(filename, kReadPermission)) {
      output << requestedRun << ",,,,,\"" << filename << "\",MISSING\n";
      ++missing;
      continue;
    }
    ++found;

    THaRun run(filename.Data(), "CDet CODA timestamp extraction");
    run.SetDataRequired(THaRunBase::kDate | THaRunBase::kRunNumber);
    const int status = run.Init();
    if (status != THaRunBase::READ_OK) {
      output << requestedRun << ",,,,,\"" << filename << "\",INIT_FAILED\n";
      ++failed;
      run.Close();
      continue;
    }

    const int codaRun = static_cast<int>(run.GetNumber());
    const TDatime &date = run.GetDate();
    const TString timestamp = date.AsSQLString();
    const TString dbHeader = TString::Format("-------[ %s ]", timestamp.Data());
    const bool numberMatches = codaRun == requestedRun;
    output << requestedRun << ',' << codaRun << ",\"" << timestamp << "\","
           << date.Convert() << ",\"" << dbHeader << "\",\"" << filename
           << "\"," << (numberMatches ? "OK" : "RUN_MISMATCH") << '\n';
    ++initialized;
    if (!numberMatches)
      ++mismatched;
    run.Close();

    if ((requestedRun - firstRun + 1) % 100 == 0)
      std::cout << "[CDet CODA timestamps] Scanned through run "
                << requestedRun << "; initialized " << initialized
                << " files.\n";
  }

  output.close();
  std::cout << "[CDet CODA timestamps] Requested runs: "
            << lastRun - firstRun + 1 << '\n'
            << "[CDet CODA timestamps] Files found/initialized: " << found
            << '/' << initialized << '\n'
            << "[CDet CODA timestamps] Missing/init failures/run mismatches: "
            << missing << '/' << failed << '/' << mismatched << '\n'
            << "[CDet CODA timestamps] Output: " << outputCSV << std::endl;
}
