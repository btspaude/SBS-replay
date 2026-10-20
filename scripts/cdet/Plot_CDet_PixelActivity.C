// Count all stored decoder slots, without ECal, timing, ToT, or spatial cuts.
// From scripts/cdet, with OUT_DIR set:
//   .L Plot_CDet_PixelActivity.C+
//   Plot_CDet_PixelActivity("CDet_run6077_pixel_activity.conf");
// Optional: override events (-1 = all), input directory, and output tag:
//   Plot_CDet_PixelActivity("CDet_run6077_pixel_activity.conf", 10000, nullptr, "test");
// Produces two visible layer maps, a CSV for all pixels, and a zero-activity list.
// Zero activity in the selected data is not proof of a dead pixel.

#include <TCanvas.h>
#include <TChain.h>
#include <TEnv.h>
#include <TH2D.h>
#include <TList.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "CDetRunDataset.h"

namespace CDetPixelActivity {
constexpr int nPixels = 2688;
constexpr int pixelsPerBar = 16;
constexpr int barsPerLayer = 84;
std::array<Long64_t, nPixels> slotCounts = {};
std::array<Long64_t, nPixels> eventCounts = {};
Long64_t eventsRead = 0;
TH2D* hSlots[2] = {};
}

void Plot_CDet_PixelActivity(
    const char* configFile = "CDet_run6077_pixel_activity.conf",
    Long64_t eventsOverride = -2, const char* inputDirectory = nullptr,
    const char* outputTag = nullptr)
{
  using namespace CDetPixelActivity;
  slotCounts.fill(0);
  eventCounts.fill(0);
  eventsRead = 0;
  for (int layer = 0; layer < 2; ++layer) {
    delete gROOT->GetListOfCanvases()->FindObject(TString::Format("cPixelActivityLayer%d", layer + 1));
    delete hSlots[layer];
    hSlots[layer] = nullptr;
  }

  // Only run-loading and output settings are read. No analysis cuts are used.
  TEnv config;
  if (!configFile || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet pixel activity] Cannot read configuration with config.version: 1.\n";
    return;
  }
  const int run = config.GetValue("analysis.run_number", 0);
  const int minSegment = config.GetValue("analysis.min_segment", -1);
  const int maxSegment = config.GetValue("analysis.max_segment", -1);
  const Long64_t requested = eventsOverride == -2
      ? config.GetValue("analysis.events", -1) : eventsOverride;
  if (run <= 0 || requested < -1 || minSegment < -1 || maxSegment < -1) {
    std::cerr << "[CDet pixel activity] Invalid run or event/segment limit.\n";
    return;
  }
  const char* directory = inputDirectory && inputDirectory[0]
      ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0]) {
    std::cerr << "[CDet pixel activity] Set OUT_DIR or pass inputDirectory.\n";
    return;
  }
  TChain tree("T");
  if (CDetRunDataset::AddToChain(&tree, run, directory, minSegment, maxSegment) <= 0) return;
  const Long64_t entries = tree.GetEntries();
  const Long64_t limit = requested < 0 ? entries : std::min(requested, entries);
  if (limit <= 0) {
    std::cerr << "[CDet pixel activity] No events to inspect; no zero-activity classification made.\n";
    return;
  }
  tree.LoadTree(0);
  if (!tree.GetBranch("earm.cdet.hits.TDCelemID")) {
    std::cerr << "[CDet pixel activity] Missing earm.cdet.hits.TDCelemID. Replay must store raw hits.\n";
    return;
  }

  TString tag = outputTag ? outputTag : config.GetValue("pixel_activity.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d", run);
  const TString allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) {
    if (allowed.First(tag[i]) == kNPOS) {
      std::cerr << "[CDet pixel activity] Output tag may contain only letters, numbers, underscores, and hyphens.\n";
      return;
    }
  }
  const TString outputDirectory = TString::Format("%s/CDetPixelActivity/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  FileStat_t outputInfo;
  if (gSystem->GetPathInfo(outputDirectory, outputInfo) != 0 || !R_ISDIR(outputInfo.fMode) ||
      gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet pixel activity] Cannot write output directory: " << outputDirectory << '\n';
    return;
  }

  // Fixed detector-index bins: these are bar/pixel indices, not a geometry map.
  for (int layer = 0; layer < 2; ++layer) {
    const int firstBar = layer * barsPerLayer;
    hSlots[layer] = new TH2D(TString::Format("hPixelActivityLayer%d", layer + 1), TString::Format("Run %d, layer %d: all stored TDC slots;Global bar ID;Pixel within bar;Decoder slots", run, layer + 1), barsPerLayer, firstBar - 0.5, firstBar + barsPerLayer - 0.5, pixelsPerBar, -0.5, pixelsPerBar - 0.5);
    hSlots[layer]->SetDirectory(nullptr);
  }

  // One decoder slot can contain LE, TE, or both. Count it even if incomplete.
  // Reading IDs alone avoids any implicit time window or LE/TE-pair requirement.
  TTreeReader reader(&tree);
  TTreeReaderArray<double> rawID(reader, "earm.cdet.hits.TDCelemID");
  std::array<Long64_t, nPixels> lastSeen;
  lastSeen.fill(-1);
  Long64_t outsidePixelIDs = 0;
  Long64_t malformedIDs = 0;
  Long64_t totalSlots = 0;
  while (eventsRead < limit && reader.Next()) {
    if (rawID.GetSetupStatus() < 0) break;
    for (size_t hit = 0; hit < rawID.GetSize(); ++hit) {
      const double id = rawID[hit];
      if (!std::isfinite(id) || std::floor(id) != id) {
        ++malformedIDs;
        continue;
      }
      if (id < 0 || id >= nPixels) {
        ++outsidePixelIDs; // Auxiliary/reference channels are not CDet pixels.
        continue;
      }
      const int pixelID = int(id);
      const int bar = pixelID / pixelsPerBar;
      const int pixel = pixelID % pixelsPerBar;
      ++slotCounts[pixelID];
      ++totalSlots;
      hSlots[bar / barsPerLayer]->Fill(bar, pixel);
      if (lastSeen[pixelID] != eventsRead) {
        ++eventCounts[pixelID];
        lastSeen[pixelID] = eventsRead;
      }
    }
    ++eventsRead;
    if (eventsRead % 1000 == 0 || eventsRead == limit)
      std::cout << "[CDet pixel activity] Event " << eventsRead << " / " << limit << '\n';
  }
  if (eventsRead != limit) {
    std::cerr << "[CDet pixel activity] Read stopped early at " << eventsRead
              << " / " << limit << "; no report written. Check input branches/files.\n";
    return;
  }

  std::ofstream csv(TString::Format("%s/pixel_activity.csv", outputDirectory.Data()).Data());
  std::ofstream zeroList(TString::Format("%s/zero_activity_pixels.txt", outputDirectory.Data()).Data());
  if (!csv || !zeroList) {
    std::cerr << "[CDet pixel activity] Cannot open output reports.\n";
    return;
  }
  csv << "run,events_read,pixel_id,layer,bar,pixel_in_bar,decoder_slots,events_with_activity,event_fraction,zero_activity\n";
  csv << std::setprecision(12);
  zeroList << "Run " << run << ": " << eventsRead << " / " << entries << " selected-chain entries inspected.\n"
           << "No event, LE, ToT, ECal, calibration, or spatial cuts.\n"
           << "Source: earm.cdet.hits.TDCelemID (all stored decoder slots, including incomplete pulses).\n"
           << "Zero activity is sample-specific; it does not establish hardware failure.\n"
           << "Non-pixel slots: " << outsidePixelIDs << "; malformed IDs: " << malformedIDs << ".\n"
           << "Input files:\n";
  for (int file = 0; file < tree.GetListOfFiles()->GetEntries(); ++file)
    zeroList << "  " << tree.GetListOfFiles()->At(file)->GetTitle() << '\n';
  zeroList << "\npixel_id layer bar pixel_in_bar\n";
  int nZero = 0;
  for (int id = 0; id < nPixels; ++id) {
    const int bar = id / pixelsPerBar;
    const int layer = bar / barsPerLayer + 1;
    csv << run << ',' << eventsRead << ',' << id << ',' << layer << ',' << bar << ','
        << id % pixelsPerBar << ',' << slotCounts[id] << ',' << eventCounts[id] << ','
        << double(eventCounts[id]) / eventsRead << ',' << (slotCounts[id] == 0) << '\n';
    if (slotCounts[id] == 0) {
      ++nZero;
      zeroList << id << ' ' << layer << ' ' << bar << ' ' << id % pixelsPerBar << '\n';
    }
  }
  zeroList << "\nZero-activity pixels: " << nZero << " / " << nPixels << '\n';
  csv.close();
  zeroList.close();
  if (!csv || !zeroList) {
    std::cerr << "[CDet pixel activity] Error writing reports; outputs may be incomplete.\n";
    return;
  }

  // Draw after filling. Preserve ROOT's mode and leave the two canvases open.
  for (int layer = 0; layer < 2; ++layer) {
    TCanvas* canvas = new TCanvas(TString::Format("cPixelActivityLayer%d", layer + 1), TString::Format("CDet pixel activity: layer %d", layer + 1), 1600, 650);
    gPad->SetRightMargin(0.13);
    hSlots[layer]->SetStats(false);
    hSlots[layer]->Draw("COLZ");
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(TString::Format("%s/layer%d_activity.pdf", outputDirectory.Data(), layer + 1));
  }
  std::cout << "[CDet pixel activity] " << totalSlots << " stored pixel slots in " << eventsRead
            << " events; " << nZero << " / " << nPixels << " pixels with zero recorded activity.\n"
            << "  Non-pixel slots: " << outsidePixelIDs << "; malformed IDs: " << malformedIDs << '\n'
            << "  Output: " << outputDirectory << '\n';
  if (totalSlots == 0)
    std::cout << "  No pixel activity anywhere: verify raw-hit storage and input data before interpreting zeros.\n";
}
