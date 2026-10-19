// Load-once copy of Plot_CDet_PulsePixelTiming.C (88fea7e6 working file).
// See Plot_CDet_PulsePixelTiming_Vectors.md for mappings, units and examples.
// .L Plot_CDet_PulsePixelTiming_Vectors.C+
// Load_CDet_PulsePixelTiming_Vectors("CDet_run6077_pulse_pixel_timing_vectors.conf");
// Plot_CDet_PulseTiming_Vectors("CDet_run6077_pulse_pixel_timing_vectors.conf");
// Cuts are evaluated on every analysis call; stored branches are never changed.
// Dataset changes require an explicit Load call. No calibration files are written.
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

#include <TFile.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <functional>
#include <memory>
#include <limits>
#include <stdexcept>

namespace CDetPulseVectors {

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

// Every branch value is saved unchanged. Named views below refer to Double_t
// leaves in the inspected Podd trees (including IDs, flags, and counters).
// Additional pulse-prefixed leaves are discovered, never silently skipped.
// Typed maps also retain additional numeric leaves without reducing precision.
template<class T> std::map<std::string, std::vector<std::vector<T>>> pulseArrays;
template<class T> std::map<std::string, std::vector<T>> pulseScalars;
std::vector<Long64_t> vEntry;
std::vector<double> vECalX, vECalY, vECalTime, vECalNclus, vECalEnergy;
bool cacheReady = false;
unsigned long long loadGeneration = 0;
std::map<std::string, std::string> loadedSchema;
std::vector<std::string> loadedFiles;
std::string loadedDatasetSettings, loadedDirectory;
Long64_t loadedEventLimit = -1;
std::map<std::string, std::vector<std::string>> canvasNames;

const auto& vPulseID = pulseArrays<Double_t>["earm.cdet.pulse.pmtnum"];
const auto& vPulseIndex = pulseArrays<Double_t>["earm.cdet.pulse.index"];
const auto& vPulseLEIndex = pulseArrays<Double_t>["earm.cdet.pulse.le_index"];
const auto& vPulseTEIndex = pulseArrays<Double_t>["earm.cdet.pulse.te_index"];
const auto& vPulseRow = pulseArrays<Double_t>["earm.cdet.pulse.row"];
const auto& vPulseCol = pulseArrays<Double_t>["earm.cdet.pulse.col"];
const auto& vPulseLayer = pulseArrays<Double_t>["earm.cdet.pulse.layer"];
const auto& vPulseX = pulseArrays<Double_t>["earm.cdet.pulse.x"];
const auto& vPulseY = pulseArrays<Double_t>["earm.cdet.pulse.y"];
const auto& vPulseZ = pulseArrays<Double_t>["earm.cdet.pulse.z"];
const auto& vPulseLe = pulseArrays<Double_t>["earm.cdet.pulse.tdc_le"];
const auto& vPulseTe = pulseArrays<Double_t>["earm.cdet.pulse.tdc_te"];
const auto& vPulseToT = pulseArrays<Double_t>["earm.cdet.pulse.tdc_tot"];
const auto& vPulseLeRaw = pulseArrays<Double_t>["earm.cdet.pulse.tdc_le_raw"];
const auto& vPulseTeRaw = pulseArrays<Double_t>["earm.cdet.pulse.tdc_te_raw"];
const auto& vPulseToTRaw = pulseArrays<Double_t>["earm.cdet.pulse.tdc_tot_raw"];
const auto& vPulseLeCorr = pulseArrays<Double_t>["earm.cdet.pulse.tdc_le_corr"];
const auto& vPulseTeCorr = pulseArrays<Double_t>["earm.cdet.pulse.tdc_te_corr"];
const auto& vPulseToTNs = pulseArrays<Double_t>["earm.cdet.pulse.tdc_tot_ns"];
const auto& vPulseDt = pulseArrays<Double_t>["earm.cdet.pulse.ecal_residual"];
const auto& vPulseCalibValid = pulseArrays<Double_t>["earm.cdet.pulse.calib_valid"];
const auto& vPulseXCorr = pulseArrays<Double_t>["earm.cdet.pulse.x_corr"];
const auto& vPulseECalXProj = pulseArrays<Double_t>["earm.cdet.pulse.ecal_x_proj"];
const auto& vPulseECalYProj = pulseArrays<Double_t>["earm.cdet.pulse.ecal_y_proj"];
const auto& vPulseDx = pulseArrays<Double_t>["earm.cdet.pulse.ecal_x_residual"];
const auto& vPulseDy = pulseArrays<Double_t>["earm.cdet.pulse.ecal_y_residual"];
const auto& vPulseBroadQualityPass = pulseArrays<Double_t>["earm.cdet.pulse.broad_quality_pass"];
const auto& vPulseECalEligible = pulseArrays<Double_t>["earm.cdet.pulse.ecal_eligible"];
const auto& vPulseSpatialPass = pulseArrays<Double_t>["earm.cdet.pulse.spatial_pass"];
const auto& vPulseN = pulseScalars<Double_t>["earm.cdet.pulse.n"];
const auto& vPulseNLEOnly = pulseScalars<Double_t>["earm.cdet.pulse.n_le_only"];
const auto& vPulseNTEOnly = pulseScalars<Double_t>["earm.cdet.pulse.n_te_only"];
const auto& vPulseNInvalidPair = pulseScalars<Double_t>["earm.cdet.pulse.n_invalid_pair"];

template<class T> void ClearTypedVectors()
{
  // Keep map nodes alive: named const references continue to refer to them.
  for (auto& field : pulseArrays<T>) decltype(field.second){}.swap(field.second);
  for (auto& field : pulseScalars<T>) decltype(field.second){}.swap(field.second);
}

void ClearCache()
{
  cacheReady = false;
  ClearTypedVectors<Double_t>(); ClearTypedVectors<Float_t>();
  ClearTypedVectors<Int_t>(); ClearTypedVectors<UInt_t>();
  ClearTypedVectors<Long64_t>(); ClearTypedVectors<ULong64_t>();
  ClearTypedVectors<Short_t>(); ClearTypedVectors<UShort_t>();
  ClearTypedVectors<Char_t>(); ClearTypedVectors<UChar_t>();
  ClearTypedVectors<Bool_t>(); ClearTypedVectors<Long_t>(); ClearTypedVectors<ULong_t>();
  std::vector<Long64_t>().swap(vEntry);
  for (auto* v : {&vECalX, &vECalY, &vECalTime, &vECalNclus, &vECalEnergy})
    std::vector<double>().swap(*v);
  loadedSchema.clear(); loadedFiles.clear(); loadedDatasetSettings.clear();
  loadedDirectory.clear(); loadedEventLimit = -1;
}

bool ReadConfig(const char* path, TEnv& config)
{
  if (!path || !path[0] || config.ReadFile(path, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet vectors] Need a readable config.version: 1 configuration.\n";
    return false;
  }
  return true;
}

TString ConfigDirectory(TEnv& config)
{
  TString directory = config.GetValue("analysis.input_directory", "");
  if (directory.IsNull()) directory = gSystem->Getenv("OUT_DIR");
  if (!directory.IsNull()) {
    gSystem->ExpandPathName(directory);
    if (!gSystem->IsAbsoluteFileName(directory))
      directory = TString(gSystem->WorkingDirectory()) + "/" + directory;
  }
  return directory;
}

std::string DatasetSettings(TEnv& config)
{
  std::ostringstream settings;
  settings << config.GetValue("analysis.run_number", 0) << '\n'
           << config.GetValue("analysis.events", -1) << '\n'
           << config.GetValue("analysis.min_segment", -1) << '\n'
           << config.GetValue("analysis.max_segment", -1) << '\n'
           << ConfigDirectory(config);
  return settings.str();
}

bool ReadCachedConfig(const char* path, TEnv& config)
{
  if (!ReadConfig(path, config)) return false;
  if (!cacheReady || DatasetSettings(config) != loadedDatasetSettings) {
    std::cerr << "[CDet vectors] No complete cache for these dataset settings. Call "
              << "Load_CDet_PulsePixelTiming_Vectors(config, directory, eventLimit) explicitly.\n";
    return false;
  }
  return true;
}

// Inspect all selected files, including rollover parts, before allocating data.
// Schema encodes primitive type plus count-leaf name, or 'scalar'. Count leaves
// themselves are consumed by TTreeReaderArray and retained by ROOT, not plotted.
std::map<std::string, std::string> InspectSchema(TTree& tree)
{
  std::map<std::string, std::string> schema;
  for (auto* object : *tree.GetListOfBranches()) {
    auto* branch = static_cast<TBranch*>(object);
    const std::string name = branch->GetName();
    const bool isPulse = name.compare(0, 15, "earm.cdet.pulse") == 0;
    const bool isECal = name == "earm.ecal.x" || name == "earm.ecal.y" ||
        name == "earm.ecal.adctime" || name == "earm.ecal.nclus" || name == "earm.ecal.e";
    if (!isPulse && !isECal) continue;
    if (branch->GetListOfLeaves()->GetEntries() != 1)
      throw std::runtime_error("Unsupported multi-leaf branch: " + name);
    auto* leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
    if (!leaf->GetLeafCount() && leaf->GetLenStatic() != 1)
      throw std::runtime_error("Unsupported fixed array: " + name);
    schema[name] = std::string(leaf->GetTypeName()) + ":" +
        (leaf->GetLeafCount() ? leaf->GetLeafCount()->GetName() : "scalar");
  }
  return schema;
}

template<class T> void BindField(TTreeReader& reader, const std::string& name,
                                 bool array, Long64_t limit,
                                 std::vector<std::function<void(size_t)>>& copyFields)
{
  if (array) {
    auto input = std::make_shared<TTreeReaderArray<T>>(reader, name.c_str());
    auto& output = pulseArrays<T>[name];
    output.reserve(limit);
    copyFields.emplace_back([input, &output, name](size_t pulses) {
      const size_t size = input->GetSize();
      if (input->GetSetupStatus() < 0 || size != pulses)
        throw std::runtime_error("Reader/array alignment failure: " + name);
      output.emplace_back();
      auto& event = output.back();
      event.reserve(size);
      for (size_t i = 0; i < size; ++i) event.push_back((*input)[i]);
      // TLeafReader updates its value-reader status, not its collection status.
      if (input->ROOT::Internal::TTreeReaderValueBase::GetReadStatus() ==
          ROOT::Internal::TTreeReaderValueBase::kReadError)
        throw std::runtime_error("Array read failure: " + name);
    });
  } else {
    auto input = std::make_shared<TTreeReaderValue<T>>(reader, name.c_str());
    auto& output = pulseScalars<T>[name];
    output.reserve(limit);
    copyFields.emplace_back([input, &output, name](size_t) {
      if (input->GetSetupStatus() < 0 || !input->Get())
        throw std::runtime_error("Scalar reader failure: " + name);
      output.push_back(**input);
    });
  }
}

template<class T> size_t TypedStorageBytes()
{
  size_t bytes = 0;
  for (const auto& field : pulseArrays<T>) {
    bytes += field.second.capacity() * sizeof(std::vector<T>);
    for (const auto& event : field.second) bytes += event.capacity() * sizeof(T);
  }
  for (const auto& field : pulseScalars<T>) bytes += field.second.capacity() * sizeof(T);
  return bytes;
}

size_t StorageBytes()
{
  size_t bytes = vEntry.capacity() * sizeof(Long64_t);
  for (auto* v : {&vECalX, &vECalY, &vECalTime, &vECalNclus, &vECalEnergy})
    bytes += v->capacity() * sizeof(double);
  return bytes + TypedStorageBytes<Double_t>() + TypedStorageBytes<Float_t>() +
      TypedStorageBytes<Int_t>() + TypedStorageBytes<UInt_t>() +
      TypedStorageBytes<Long64_t>() + TypedStorageBytes<ULong64_t>() +
      TypedStorageBytes<Short_t>() + TypedStorageBytes<UShort_t>() +
      TypedStorageBytes<Char_t>() + TypedStorageBytes<UChar_t>() +
      TypedStorageBytes<Bool_t>() + TypedStorageBytes<Long_t>() + TypedStorageBytes<ULong_t>();
}

bool Load(const char* configFile, const char* inputDirectory, Long64_t eventsOverride)
{
  TStopwatch clock;
  ClearCache(); // An unsuccessful reload leaves no usable partial/old cache.
  try {
    TEnv config;
    if (!ReadConfig(configFile, config)) return false;
    const int run = config.GetValue("analysis.run_number", 0);
    const Long64_t events = eventsOverride == -2 ? config.GetValue("analysis.events", -1) : eventsOverride;
    TString directory = inputDirectory && inputDirectory[0] ? TString(inputDirectory) : ConfigDirectory(config);
    gSystem->ExpandPathName(directory);
    if (run <= 0 || events < -1 || directory.IsNull())
      throw std::runtime_error("Invalid run/event limit or missing input directory (set OUT_DIR).");
    if (!gSystem->IsAbsoluteFileName(directory)) directory = TString(gSystem->WorkingDirectory()) + "/" + directory;
    TChain chain("T");
    const auto files = CDetRunDataset::Select(run, directory,
        config.GetValue("analysis.min_segment", -1), config.GetValue("analysis.max_segment", -1));
    if (files.empty()) throw std::runtime_error("No selected replay files.");
    for (const auto& path : files) {
      std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
      TTree* tree = file && !file->IsZombie() ? file->Get<TTree>("T") : nullptr;
      if (!tree) throw std::runtime_error("Cannot read T from " + std::string(path.Data()));
      const auto schema = InspectSchema(*tree);
      if (loadedSchema.empty()) loadedSchema = schema;
      else if (schema != loadedSchema) throw std::runtime_error("Pulse/ECal schema mismatch in " + std::string(path.Data()));
      loadedFiles.push_back(path.Data());
      if (!chain.Add(path)) throw std::runtime_error("Cannot chain " + std::string(path.Data()));
    }
    // The documented 33 fields are a required baseline, not an allow-list.
    const char* requiredArrays[] = {"pmtnum", "index", "le_index", "te_index", "row", "col", "layer", "x", "y", "z", "tdc_le", "tdc_te", "tdc_tot", "tdc_le_raw", "tdc_te_raw", "tdc_tot_raw", "tdc_le_corr", "tdc_te_corr", "tdc_tot_ns", "ecal_residual", "calib_valid", "x_corr", "ecal_x_proj", "ecal_y_proj", "ecal_x_residual", "ecal_y_residual", "broad_quality_pass", "ecal_eligible", "spatial_pass"};
    for (const char* suffix : requiredArrays) {
      const std::string name = std::string("earm.cdet.pulse.") + suffix;
      const auto found = loadedSchema.find(name);
      if (found == loadedSchema.end() || found->second.find("Double_t:") != 0 || found->second == "Double_t:scalar")
        throw std::runtime_error("Missing or changed baseline Double_t pulse array: " + name);
    }
    for (const char* suffix : {"n", "n_le_only", "n_te_only", "n_invalid_pair"}) {
      const std::string name = std::string("earm.cdet.pulse.") + suffix;
      if (!loadedSchema.count(name) || loadedSchema.at(name) != "Double_t:scalar")
        throw std::runtime_error("Missing or changed baseline Double_t pulse counter: " + name);
    }
    for (const char* suffix : {"x", "y", "adctime", "nclus", "e"}) {
      const std::string name = std::string("earm.ecal.") + suffix;
      if (!loadedSchema.count(name) || loadedSchema.at(name) != "Double_t:scalar")
        throw std::runtime_error("Missing or changed Double_t ECal scalar: " + name);
    }
    const Long64_t total = chain.GetEntries();
    if (total <= 0) throw std::runtime_error("Selected dataset has no entries.");
    const Long64_t limit = events >= 0 ? std::min(events, total) : total;
    TTreeReader reader(&chain);
    std::vector<std::function<void(size_t)>> copyFields;
    for (const auto& field : loadedSchema) {
      if (field.first.compare(0, 15, "earm.cdet.pulse") != 0) continue;
      const auto separator = field.second.find(':');
      const std::string type = field.second.substr(0, separator);
      const bool array = field.second.substr(separator + 1) != "scalar";
      std::cout << "[CDet vectors] " << field.first << " -> " << field.second << '\n';
      if (type == "Double_t") BindField<Double_t>(reader, field.first, array, limit, copyFields);
      else if (type == "Float_t") BindField<Float_t>(reader, field.first, array, limit, copyFields);
      else if (type == "Int_t") BindField<Int_t>(reader, field.first, array, limit, copyFields);
      else if (type == "UInt_t") BindField<UInt_t>(reader, field.first, array, limit, copyFields);
      else if (type == "Long64_t") BindField<Long64_t>(reader, field.first, array, limit, copyFields);
      else if (type == "ULong64_t") BindField<ULong64_t>(reader, field.first, array, limit, copyFields);
      else if (type == "Short_t") BindField<Short_t>(reader, field.first, array, limit, copyFields);
      else if (type == "UShort_t") BindField<UShort_t>(reader, field.first, array, limit, copyFields);
      else if (type == "Char_t") BindField<Char_t>(reader, field.first, array, limit, copyFields);
      else if (type == "UChar_t") BindField<UChar_t>(reader, field.first, array, limit, copyFields);
      else if (type == "Bool_t") BindField<Bool_t>(reader, field.first, array, limit, copyFields);
      else if (type == "Long_t") BindField<Long_t>(reader, field.first, array, limit, copyFields);
      else if (type == "ULong_t") BindField<ULong_t>(reader, field.first, array, limit, copyFields);
      else throw std::runtime_error("Unsupported pulse type (nothing omitted): " + field.first + " " + type);
    }
    TTreeReaderValue<Double_t> pulseN(reader, "earm.cdet.pulse.n");
    TTreeReaderValue<Double_t> x(reader, "earm.ecal.x"), y(reader, "earm.ecal.y");
    TTreeReaderValue<Double_t> time(reader, "earm.ecal.adctime"), nclus(reader, "earm.ecal.nclus"), energy(reader, "earm.ecal.e");
    vEntry.reserve(limit);
    for (auto* v : {&vECalX, &vECalY, &vECalTime, &vECalNclus, &vECalEnergy}) v->reserve(limit);
    size_t totalPulses = 0, maxPulses = 0, zeroEvents = 0;
    for (Long64_t entry = 0; entry < limit; ++entry) {
      if (!reader.Next()) throw std::runtime_error("Reader stopped at entry " + std::to_string(entry) + ", status " + std::to_string(reader.GetEntryStatus()));
      if (!pulseN.Get() || !std::isfinite(*pulseN) || *pulseN < 0 || std::floor(*pulseN) != *pulseN)
        throw std::runtime_error("Invalid pulse.n at entry " + std::to_string(entry));
      const size_t pulses = size_t(*pulseN);
      for (auto& copy : copyFields) copy(pulses);
      for (auto* value : {&x, &y, &time, &nclus, &energy})
        if (value->GetSetupStatus() < 0 || !value->Get()) throw std::runtime_error("ECal scalar read failure.");
      vEntry.push_back(reader.GetCurrentEntry());
      vECalX.push_back(*x); vECalY.push_back(*y); vECalTime.push_back(*time);
      vECalNclus.push_back(*nclus); vECalEnergy.push_back(*energy);
      totalPulses += pulses; maxPulses = std::max(maxPulses, pulses);
      if (!pulses) ++zeroEvents;
    }
    if (Long64_t(vEntry.size()) != limit) throw std::runtime_error("Incomplete cache.");
    if (limit == total && (reader.Next() || reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd))
      throw std::runtime_error("Unexpected reader completion status.");
    // Every field appended once per event; every array checked against pulse.n.
    loadedDatasetSettings = DatasetSettings(config);
    loadedDirectory = directory.Data(); loadedEventLimit = events;
    cacheReady = true; ++loadGeneration;
    clock.Stop();
    std::cout << "[CDet vectors] Complete load: " << vEntry.size() << " events, " << totalPulses
              << " pulses; mean/max multiplicity " << (limit ? double(totalPulses)/limit : 0) << '/' << maxPulses
              << "; zero-pulse events " << zeroEvents << "; " << loadedSchema.size()-5 << " pulse fields.\n"
              << "  Load wall time " << clock.RealTime() << " s; estimated vector allocation "
              << StorageBytes()/1048576.0 << " MiB (capacity + nested-vector objects; excludes allocator/map/ROOT overhead; bool capacity is conservative).\n"
              << "  Cache generation " << loadGeneration << "; " << files.size() << " selected files; entry limit " << events << ".\n";
    return true;
  } catch (const std::bad_alloc&) {
    ClearCache(); std::cerr << "[CDet vectors] Insufficient memory for all fields; cache cleared. No precision/field reduction attempted.\n";
  } catch (const std::exception& error) {
    ClearCache(); std::cerr << "[CDet vectors] Load failed: " << error.what() << "; cache cleared.\n";
  }
  return false;
}

bool PrepareOutput(TEnv& config, int& bar, const char* outputTag,
                   const char* family, TString& outputDirectory)
{
  if (bar == -1) bar = config.GetValue("pulse_pixel_timing.study_bar", 19);
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const int runNumber = config.GetValue("analysis.run_number", 0);
  if (bar < 0 || bar >= 168 || controlBar < 0 || controlBar >= 168) {
    std::cerr << "[CDet vectors] Bars must be in 0..167.\n"; return false;
  }
  TString tag = outputTag ? outputTag : config.GetValue("pulse_pixel_timing.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d_bar%d_control%d", runNumber, bar, controlBar);
  const TString allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i=0; i<tag.Length(); ++i) if (allowed.First(tag[i]) == kNPOS) {
    std::cerr << "[CDet vectors] Invalid output tag.\n"; return false;
  }
  outputDirectory = TString::Format("%s/PulsePixelTimingVectors/%s", gSystem->WorkingDirectory(), tag.Data());
  if (gSystem->AccessPathName(outputDirectory)) gSystem->mkdir(outputDirectory, true);
  if (gSystem->AccessPathName(outputDirectory, kWritePermission)) {
    std::cerr << "[CDet vectors] Cannot write output directory.\n"; return false;
  }
  // Close only this family's old canvases before replacing its histograms.
  for (const auto& name : canvasNames[family]) delete gROOT->GetListOfCanvases()->FindObject(name.c_str());
  canvasNames[family].clear();
  return true;
}

void SaveCanvases(const char* family, const TString& directory, std::initializer_list<TCanvas*> canvases)
{
  for (auto* canvas : canvases) {
    canvasNames[family].push_back(canvas->GetName());
    canvas->SaveAs(TString::Format("%s/%s.pdf", directory.Data(), canvas->GetName()));
  }
}

bool GoodEnergy(size_t event, double minimum, double maximum)
{
  return std::isfinite(vECalEnergy[event]) && vECalEnergy[event] > minimum && vECalEnergy[event] < maximum;
}

bool GoodBeforeEllipse(size_t event, size_t pulse, double energyMin, double energyMax,
                       double leMin, double leMax, bool applyToT, double totMin, double totMax)
{
  const double leNs = 0.01*vPulseLe[event][pulse]; // Legacy units only; no correction applied.
  const double totNs = vPulseToTNs[event][pulse];
  const bool good_values = std::isfinite(vPulseID[event][pulse]) && std::isfinite(vPulseDx[event][pulse]) &&
      std::isfinite(vPulseDt[event][pulse]) && std::isfinite(leNs) && std::isfinite(totNs) && totNs > 0.0;
  const bool good_replay = vPulseCalibValid[event][pulse] > 0.5 &&
      vPulseECalEligible[event][pulse] > 0.5 && vPulseSpatialPass[event][pulse] > 0.5;
  const bool good_le = leNs >= leMin && leNs <= leMax;
  const bool good_tot = !applyToT || (totNs >= totMin && totNs <= totMax);
  // Deliberately no broad_quality_pass or new ECal time gate, matching original.
  return GoodEnergy(event, energyMin, energyMax) && good_values && good_replay && good_le && good_tot;
}

bool GoodEllipse(size_t event, size_t pulse, double radius, double xCenter,
                 double timeCenter, double xScale, double timeScale)
{
  if (radius == -1.0) return true;
  const double x = (vPulseDx[event][pulse] - xCenter)/xScale;
  const double t = (vPulseDt[event][pulse] - timeCenter)/timeScale;
  return x*x + t*t <= radius*radius;
}

bool ValidRange(double minimum, double maximum)
{
  return std::isfinite(minimum) && std::isfinite(maximum) && minimum < maximum;
}

bool ValidBins(double width, double minimum, double maximum)
{
  return ValidRange(minimum, maximum) && std::isfinite(width) && width > 0 &&
      std::ceil((maximum-minimum)/width) <= std::numeric_limits<int>::max();
}

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
  const TString frameName = TString::Format("pv_projectionFrame_%d_%d", firstPixel, int(fractions));
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


bool Load_CDet_PulsePixelTiming_Vectors(
    const char* configFile = "CDet_run6077_pulse_pixel_timing_vectors.conf",
    const char* inputDirectory = nullptr, Long64_t eventsOverride = -2)
{
  return CDetPulseVectors::Load(configFile, inputDirectory, eventsOverride);
}

bool Plot_CDet_PulseTiming_Vectors(
    const char* configFile = "CDet_run6077_pulse_pixel_timing_vectors.conf", int bar = -1,
    double ecalDTBinWidthNs = NAN, double ecalDTMinNs = NAN, double ecalDTMaxNs = NAN,
    double ecalDXBinWidthM = NAN, double ecalDXMinM = NAN, double ecalDXMaxM = NAN,
    double dtVsDxMinNs = NAN, double dtVsDxMaxNs = NAN, const char* outputTag = nullptr)
{
  TStopwatch analysisClock;
  TEnv config;
  if (!ReadCachedConfig(configFile, config)) return false;
  const int runNumber = config.GetValue("analysis.run_number", 0);
  const int pixelsPerBar = 16, nPixels = 2688;
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double energyMin = config.GetValue("analysis.ecal_energy_min", NAN);
  const double energyMax = config.GetValue("analysis.ecal_energy_max", NAN);
  if (!ValidRange(energyMin, energyMax)) { std::cerr << "[CDet vectors] Invalid ECal energy interval.\n"; return false; }
  const double leCutMinNs = config.GetValue("analysis.le_min", 0.02);
  const double leCutMaxNs = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_pixel_timing.apply_tot_cut", 0) != 0;
  const double totCutMinNs = config.GetValue("analysis.tot_min", 4.0);
  const double totCutMaxNs = config.GetValue("analysis.tot_max", 30.0);
  const double fitMin = config.GetValue("pulse_pixel_timing.fit_min_ns", -55.0);
  const double fitMax = config.GetValue("pulse_pixel_timing.fit_max_ns", -10.0);
  const int minFitEntries = config.GetValue("pulse_pixel_timing.min_fit_entries", 35);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const bool useEllipse = ellipseRadius > 0.0;
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  if (!std::isfinite(leCutMinNs) || !std::isfinite(leCutMaxNs) || leCutMinNs >= leCutMaxNs ||
      (applyToTCut && (!std::isfinite(totCutMinNs) || !std::isfinite(totCutMaxNs) ||
       totCutMinNs < 0.0 || totCutMinNs >= totCutMaxNs))) {
    std::cerr << "[CDet pulse pixel timing] Invalid analysis.le_min/max or enabled analysis.tot_min/max (ns).\n";
    return false;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse pixel timing] Ellipse radius must be -1 (off) or positive; enabled ellipses require finite centers and positive finite scales.\n";
    return false;
  }
  if (std::isnan(ecalDTBinWidthNs)) ecalDTBinWidthNs = config.GetValue("vectors.display.dt_bin_width_ns", 1.0);
  if (std::isnan(ecalDXBinWidthM)) ecalDXBinWidthM = config.GetValue("vectors.display.dx_bin_width_m", 0.005);
  if (std::isnan(ecalDXMinM)) ecalDXMinM = config.GetValue("vectors.display.dx_min_m", -0.10);
  if (std::isnan(ecalDXMaxM)) ecalDXMaxM = config.GetValue("vectors.display.dx_max_m", 0.10);
  if (std::isnan(ecalDTMinNs)) ecalDTMinNs = config.GetValue("display.ecal_cdet_dt_min", -150.0);
  if (std::isnan(ecalDTMaxNs)) ecalDTMaxNs = config.GetValue("display.ecal_cdet_dt_max", 150.0);
  if (std::isnan(dtVsDxMinNs)) dtVsDxMinNs = config.GetValue("vectors.display.dt_vs_dx_min_ns", ecalDTMinNs);
  if (std::isnan(dtVsDxMaxNs)) dtVsDxMaxNs = config.GetValue("vectors.display.dt_vs_dx_max_ns", ecalDTMaxNs);
  if (!ValidBins(ecalDTBinWidthNs, ecalDTMinNs, ecalDTMaxNs) ||
      !ValidBins(ecalDXBinWidthM, ecalDXMinM, ecalDXMaxM) ||
      !ValidBins(ecalDTBinWidthNs, dtVsDxMinNs, dtVsDxMaxNs) ||
      !ValidRange(fitMin, fitMax) || fitMin < ecalDTMinNs || fitMax > ecalDTMaxNs ||
      fitMax-fitMin <= 2*ecalDTBinWidthNs || minFitEntries < 6) {
    std::cerr << "[CDet vectors] Invalid timing bins/fit settings.\n"; return false;
  }
  TString outputDirectory;
  if (!PrepareOutput(config, bar, outputTag, "timing", outputDirectory)) return false;
  const int firstPixel = pixelsPerBar*bar, firstControlPixel = pixelsPerBar*controlBar;
  const int nECalDT = int(std::ceil((ecalDTMaxNs-ecalDTMinNs)/ecalDTBinWidthNs));
  const int nECalDX = int(std::ceil((ecalDXMaxM-ecalDXMinM)/ecalDXBinWidthM));
  const int nDtVsDx = int(std::ceil((dtVsDxMaxNs-dtVsDxMinNs)/ecalDTBinWidthNs));
  const char* selectionLabel = useEllipse ? "ellipse selected" : "ellipse off";
  for (int i = 0; i < pixelsPerBar; ++i) {
    const int pixelID = firstPixel + i;
    const TString name = TString::Format("pv_hPixelDt_run%d_pixel%d", runNumber, pixelID);
    const TString label = TString::Format("t_{ECal} - t_{CDet,%d}", pixelID);
    const TString title = TString::Format("Bar %d: %s (%s);%s (ns);Counts", bar, label.Data(), selectionLabel, label.Data());
    delete hPixelDt[i]; // Replace the previous histogram when rerunning.
    hPixelDt[i] = new TH1D(name, title, nECalDT, ecalDTMinNs, ecalDTMaxNs);
    hPixelDt[i]->SetDirectory(nullptr);
    hPixelDt[i]->SetStats(false);

    const int controlPixelID = firstControlPixel + i;
    const TString controlName = TString::Format("pv_hControlPixelDt_run%d_pixel%d", runNumber, controlPixelID);
    const TString controlLabel = TString::Format("t_{ECal} - t_{CDet,%d}", controlPixelID);
    const TString controlTitle = TString::Format("Control bar %d: %s (%s);%s (ns);Counts", controlBar, controlLabel.Data(), selectionLabel, controlLabel.Data());
    delete hControlPixelDt[i];
    hControlPixelDt[i] = new TH1D(controlName, controlTitle, nECalDT, ecalDTMinNs, ecalDTMaxNs);
    hControlPixelDt[i]->SetDirectory(nullptr);
    hControlPixelDt[i]->SetStats(false);

  }
  delete hBarDt;
  delete hControlBarDt;
  hBarDt = new TH1D(TString::Format("pv_hBarDt_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: all 16 pixels (%s);t_{ECal} - t_{CDet} (ns);Counts", runNumber, bar, selectionLabel), nECalDT, ecalDTMinNs, ecalDTMaxNs);
  hControlBarDt = new TH1D(TString::Format("pv_hControlBarDt_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: all 16 pixels (%s);t_{ECal} - t_{CDet} (ns);Counts", runNumber, controlBar, selectionLabel), nECalDT, ecalDTMinNs, ecalDTMaxNs);
  for (TH1D* hist : {hBarDt, hControlBarDt}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

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
  hBarDtVsDxBefore = new TH2D(TString::Format("pv_hBarDtVsDxBefore_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: before ellipse%s", runNumber, bar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  hControlBarDtVsDxBefore = new TH2D(TString::Format("pv_hControlBarDtVsDxBefore_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: before ellipse%s", runNumber, controlBar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  hDetectorDtVsDxBefore = new TH2D(TString::Format("pv_hDetectorDtVsDxBefore_run%d", runNumber), TString::Format("Run %d, all CDet pixels: before ellipse%s", runNumber, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  if (useEllipse) {
    hBarDtVsDxAfter = new TH2D(TString::Format("pv_hBarDtVsDxAfter_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: after ellipse%s", runNumber, bar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
    hControlBarDtVsDxAfter = new TH2D(TString::Format("pv_hControlBarDtVsDxAfter_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: after ellipse%s", runNumber, controlBar, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
    hDetectorDtVsDxAfter = new TH2D(TString::Format("pv_hDetectorDtVsDxAfter_run%d", runNumber), TString::Format("Run %d, all CDet pixels: after ellipse%s", runNumber, axes), nECalDX, ecalDXMinM, ecalDXMaxM, nDtVsDx, dtVsDxMinNs, dtVsDxMaxNs);
  }
  for (TH2D* hist : {hBarDtVsDxBefore, hControlBarDtVsDxBefore, hBarDtVsDxAfter, hControlBarDtVsDxAfter, hDetectorDtVsDxBefore, hDetectorDtVsDxAfter}) {
    if (!hist) continue;
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

  Long64_t detectorPulses = 0, selectedPulses = 0, controlPulses = 0;
  for (size_t event=0; event<vEntry.size(); ++event) {
    for (size_t i=0; i<vPulseID[event].size(); ++i) {
      const bool good_before_ellipse = GoodBeforeEllipse(event, i, energyMin, energyMax, leCutMinNs, leCutMaxNs, applyToTCut, totCutMinNs, totCutMaxNs);
      if (!good_before_ellipse) continue;
      const int id = int(vPulseID[event][i]);
      const int iPixel = id-firstPixel, iControlPixel = id-firstControlPixel;
      const bool inBar = iPixel>=0 && iPixel<16, inControlBar = iControlPixel>=0 && iControlPixel<16;
      const double dt = vPulseDt[event][i], dx = vPulseDx[event][i];
      hDetectorDtVsDxBefore->Fill(dx, dt);
      if (inBar) hBarDtVsDxBefore->Fill(dx, dt);
      if (inControlBar) hControlBarDtVsDxBefore->Fill(dx, dt);
      const bool good_final = GoodEllipse(event, i, ellipseRadius, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale);
      if (!good_final) continue;
      ++detectorPulses;
      if (inBar) { ++selectedPulses; hPixelDt[iPixel]->Fill(dt); }
      if (inControlBar) { ++controlPulses; hControlPixelDt[iControlPixel]->Fill(dt); }
      if (useEllipse) {
        hDetectorDtVsDxAfter->Fill(dx, dt);
        if (inBar) hBarDtVsDxAfter->Fill(dx, dt);
        if (inControlBar) hControlBarDtVsDxAfter->Fill(dx, dt);
      }
    }
  }
  for (int i=0; i<16; ++i) { hBarDt->Add(hPixelDt[i]); hControlBarDt->Add(hControlPixelDt[i]); }
  // 6. Fit and draw both bars: one pixel canvas and one summed canvas each.
  TCanvas *cPixelDt = new TCanvas(TString::Format("pv_cPixelDt_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: pulse timing", runNumber, bar), 1400, 1000);
  cPixelDt->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cPixelDt->cd(i + 1);
    DrawSpectrum(hPixelDt[i], fitMin, fitMax, minFitEntries);
  }
  cPixelDt->Modified();
  cPixelDt->Update();

  TCanvas *cBarDt = new TCanvas(TString::Format("pv_cBarDt_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: summed pulse timing", runNumber, bar), 900, 650);
  DrawSpectrum(hBarDt, fitMin, fitMax, minFitEntries);
  cBarDt->Modified();
  cBarDt->Update();

  TCanvas *cControlPixelDt = new TCanvas(TString::Format("pv_cControlPixelDt_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: pulse timing", runNumber, controlBar), 1400, 1000);
  cControlPixelDt->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cControlPixelDt->cd(i + 1);
    DrawSpectrum(hControlPixelDt[i], fitMin, fitMax, minFitEntries);
  }
  cControlPixelDt->Modified();
  cControlPixelDt->Update();

  TCanvas *cControlBarDt = new TCanvas(TString::Format("pv_cControlBarDt_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: summed pulse timing", runNumber, controlBar), 900, 650);
  DrawSpectrum(hControlBarDt, fitMin, fitMax, minFitEntries);
  cControlBarDt->Modified();
  cControlBarDt->Update();

  // 7. Before/after panels use identical axes and color scales within each bar.
  // Reusing canvas names replaces any previous after panels when rerunning off.
  TCanvas *cBarDtVsDx = new TCanvas(TString::Format("pv_cBarDtVsDx_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: dt vs dx", runNumber, bar), useEllipse ? 1400 : 800, 650);
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

  TCanvas *cControlBarDtVsDx = new TCanvas(TString::Format("pv_cControlBarDtVsDx_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: dt vs dx", runNumber, controlBar), useEllipse ? 1400 : 800, 650);
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

  TCanvas *cDetectorDtVsDx = new TCanvas(TString::Format("pv_cDetectorDtVsDx_run%d", runNumber), TString::Format("CDet run %d, entire detector: dt vs dx", runNumber), useEllipse ? 1400 : 800, 650);
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

  SaveCanvases("timing", outputDirectory, {cPixelDt,cBarDt,cControlPixelDt,cControlBarDt,cBarDtVsDx,cControlBarDtVsDx,cDetectorDtVsDx});
  std::cout << "[CDet vectors timing] study/control/detector selected pulses " << selectedPulses << '/' << controlPulses << '/' << detectorPulses << '\n';
  analysisClock.Stop();
  std::cout << "[CDet vectors timing] " << vEntry.size() << " cached events; wall time "
            << analysisClock.RealTime() << " s including drawing/export; generation " << loadGeneration << '\n';
  return true;
}

bool Plot_CDet_PulseLEToT_Vectors(
    const char* configFile = "CDet_run6077_pulse_pixel_timing_vectors.conf", int bar = -1,
    double leBinWidthNs = NAN, double leMinNs = NAN, double leMaxNs = NAN,
    double totBinWidthNs = NAN, double totMinNs = NAN, double totMaxNs = NAN,
    const char* outputTag = nullptr)
{
  TStopwatch analysisClock;
  TEnv config;
  if (!ReadCachedConfig(configFile, config)) return false;
  const int runNumber = config.GetValue("analysis.run_number", 0);
  const int pixelsPerBar = 16, nPixels = 2688;
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double energyMin = config.GetValue("analysis.ecal_energy_min", NAN);
  const double energyMax = config.GetValue("analysis.ecal_energy_max", NAN);
  if (!ValidRange(energyMin, energyMax)) { std::cerr << "[CDet vectors] Invalid ECal energy interval.\n"; return false; }
  const double leCutMinNs = config.GetValue("analysis.le_min", 0.02);
  const double leCutMaxNs = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_pixel_timing.apply_tot_cut", 0) != 0;
  const double totCutMinNs = config.GetValue("analysis.tot_min", 4.0);
  const double totCutMaxNs = config.GetValue("analysis.tot_max", 30.0);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const bool useEllipse = ellipseRadius > 0.0;
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  if (!std::isfinite(leCutMinNs) || !std::isfinite(leCutMaxNs) || leCutMinNs >= leCutMaxNs ||
      (applyToTCut && (!std::isfinite(totCutMinNs) || !std::isfinite(totCutMaxNs) ||
       totCutMinNs < 0.0 || totCutMinNs >= totCutMaxNs))) {
    std::cerr << "[CDet pulse pixel timing] Invalid analysis.le_min/max or enabled analysis.tot_min/max (ns).\n";
    return false;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse pixel timing] Ellipse radius must be -1 (off) or positive; enabled ellipses require finite centers and positive finite scales.\n";
    return false;
  }
  if (std::isnan(leBinWidthNs)) leBinWidthNs = config.GetValue("vectors.display.le_bin_width_ns", 1.0);
  if (std::isnan(leMinNs)) leMinNs = config.GetValue("vectors.display.le_min_ns", -20.0);
  if (std::isnan(leMaxNs)) leMaxNs = config.GetValue("vectors.display.le_max_ns", 100.0);
  if (std::isnan(totBinWidthNs)) totBinWidthNs = config.GetValue("vectors.display.tot_bin_width_ns", 1.0);
  if (std::isnan(totMinNs)) totMinNs = config.GetValue("vectors.display.tot_min_ns", 0.0);
  if (std::isnan(totMaxNs)) totMaxNs = config.GetValue("vectors.display.tot_max_ns", 100.0);
  if (!ValidBins(leBinWidthNs, leMinNs, leMaxNs) || !ValidBins(totBinWidthNs, totMinNs, totMaxNs)) {
    std::cerr << "[CDet vectors] Invalid LE/ToT binning.\n"; return false;
  }
  TString outputDirectory;
  if (!PrepareOutput(config, bar, outputTag, "le_tot", outputDirectory)) return false;
  const int firstPixel = pixelsPerBar*bar, firstControlPixel = pixelsPerBar*controlBar;
  const int nLe = int(std::ceil((leMaxNs-leMinNs)/leBinWidthNs));
  const int nToT = int(std::ceil((totMaxNs-totMinNs)/totBinWidthNs));
  const char* selectionLabel = useEllipse ? "ellipse selected" : "ellipse off";
  for (int i=0; i<16; ++i) {
    const int pixelID=firstPixel+i, controlPixelID=firstControlPixel+i;
    delete hPixelLe[i];
    delete hControlPixelLe[i];
    delete hPixelLeVsToT[i];
    delete hControlPixelLeVsToT[i];
    hPixelLe[i] = new TH1D(TString::Format("pv_hPixelLe_run%d_pixel%d", runNumber, pixelID), TString::Format("Bar %d, pixel %d: LE (%s);t_{LE} (ns);Counts", bar, pixelID, selectionLabel), nLe, leMinNs, leMaxNs);
    hControlPixelLe[i] = new TH1D(TString::Format("pv_hControlPixelLe_run%d_pixel%d", runNumber, controlPixelID), TString::Format("Control bar %d, pixel %d: LE (%s);t_{LE} (ns);Counts", controlBar, controlPixelID, selectionLabel), nLe, leMinNs, leMaxNs);
    hPixelLeVsToT[i] = new TH2D(TString::Format("pv_hPixelLeVsToT_run%d_pixel%d", runNumber, pixelID), TString::Format("Bar %d, pixel %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", bar, pixelID, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
    hControlPixelLeVsToT[i] = new TH2D(TString::Format("pv_hControlPixelLeVsToT_run%d_pixel%d", runNumber, controlPixelID), TString::Format("Control bar %d, pixel %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", controlBar, controlPixelID, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
    TH1* pixelPulsePlots[] = {hPixelLe[i], hControlPixelLe[i], hPixelLeVsToT[i], hControlPixelLeVsToT[i]};
    for (TH1* hist : pixelPulsePlots) {
      hist->SetDirectory(nullptr);
      hist->SetStats(false);
    }
  }
  delete hBarLe;
  delete hControlBarLe;
  delete hBarLeVsToT;
  delete hControlBarLeVsToT;
  hBarLe = new TH1D(TString::Format("pv_hBarLe_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: LE (%s);t_{LE} (ns);Counts", runNumber, bar, selectionLabel), nLe, leMinNs, leMaxNs);
  hControlBarLe = new TH1D(TString::Format("pv_hControlBarLe_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: LE (%s);t_{LE} (ns);Counts", runNumber, controlBar, selectionLabel), nLe, leMinNs, leMaxNs);
  hBarLeVsToT = new TH2D(TString::Format("pv_hBarLeVsToT_run%d_bar%d", runNumber, bar), TString::Format("Run %d, bar %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", runNumber, bar, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
  hControlBarLeVsToT = new TH2D(TString::Format("pv_hControlBarLeVsToT_run%d_bar%d", runNumber, controlBar), TString::Format("Run %d, control bar %d: LE vs ToT (%s);ToT (ns);t_{LE} (ns);Pulses", runNumber, controlBar, selectionLabel), nToT, totMinNs, totMaxNs, nLe, leMinNs, leMaxNs);
  TH1* barPulsePlots[] = {hBarLe, hControlBarLe, hBarLeVsToT, hControlBarLeVsToT};
  for (TH1* hist : barPulsePlots) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

  Long64_t detectorPulses = 0;
  for (size_t event=0; event<vEntry.size(); ++event) {
    for (size_t i=0; i<vPulseID[event].size(); ++i) {
      const bool good_final = GoodBeforeEllipse(event, i, energyMin, energyMax, leCutMinNs, leCutMaxNs, applyToTCut, totCutMinNs, totCutMaxNs) && GoodEllipse(event, i, ellipseRadius, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale);
      if (!good_final) continue;
      ++detectorPulses;
      const int id=int(vPulseID[event][i]);
      const double leNs=0.01*vPulseLe[event][i], totNs=vPulseToTNs[event][i];
      if (id>=firstPixel && id<firstPixel+16) {
        hPixelLe[id-firstPixel]->Fill(leNs); hPixelLeVsToT[id-firstPixel]->Fill(totNs,leNs);
      }
      if (id>=firstControlPixel && id<firstControlPixel+16) {
        hControlPixelLe[id-firstControlPixel]->Fill(leNs); hControlPixelLeVsToT[id-firstControlPixel]->Fill(totNs,leNs);
      }
    }
  }
  for (int i=0; i<16; ++i) {
    hBarLe->Add(hPixelLe[i]); hControlBarLe->Add(hControlPixelLe[i]);
    hBarLeVsToT->Add(hPixelLeVsToT[i]); hControlBarLeVsToT->Add(hControlPixelLeVsToT[i]);
  }
  // 8. Pre-offset LE and LE-versus-ToT for the exact pulses in the 1D dt plots.
  TCanvas *cPixelLe = new TCanvas(TString::Format("pv_cPixelLe_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: LE before pixel offsets", runNumber, bar), 1400, 1000);
  cPixelLe->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cPixelLe->cd(i + 1);
    DrawPulseDistribution(hPixelLe[i], false);
  }
  cPixelLe->Modified();
  cPixelLe->Update();

  TCanvas *cControlPixelLe = new TCanvas(TString::Format("pv_cControlPixelLe_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: LE before pixel offsets", runNumber, controlBar), 1400, 1000);
  cControlPixelLe->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cControlPixelLe->cd(i + 1);
    DrawPulseDistribution(hControlPixelLe[i], false);
  }
  cControlPixelLe->Modified();
  cControlPixelLe->Update();

  TCanvas *cPixelLeVsToT = new TCanvas(TString::Format("pv_cPixelLeVsToT_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: LE vs ToT", runNumber, bar), 1400, 1000);
  cPixelLeVsToT->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cPixelLeVsToT->cd(i + 1);
    DrawPulseDistribution(hPixelLeVsToT[i], true);
  }
  cPixelLeVsToT->Modified();
  cPixelLeVsToT->Update();

  TCanvas *cControlPixelLeVsToT = new TCanvas(TString::Format("pv_cControlPixelLeVsToT_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: LE vs ToT", runNumber, controlBar), 1400, 1000);
  cControlPixelLeVsToT->Divide(4, 4);
  for (int i = 0; i < pixelsPerBar; ++i) {
    cControlPixelLeVsToT->cd(i + 1);
    DrawPulseDistribution(hControlPixelLeVsToT[i], true);
  }
  cControlPixelLeVsToT->Modified();
  cControlPixelLeVsToT->Update();

  TCanvas *cBarLeToT = new TCanvas(TString::Format("pv_cBarLeToT_run%d_bar%d", runNumber, bar), TString::Format("CDet run %d, bar %d: summed LE and LE vs ToT", runNumber, bar), 1400, 650);
  cBarLeToT->Divide(2, 1);
  cBarLeToT->cd(1);
  DrawPulseDistribution(hBarLe, false);
  cBarLeToT->cd(2);
  DrawPulseDistribution(hBarLeVsToT, true);
  cBarLeToT->Modified();
  cBarLeToT->Update();

  TCanvas *cControlBarLeToT = new TCanvas(TString::Format("pv_cControlBarLeToT_run%d_bar%d", runNumber, controlBar), TString::Format("CDet run %d, control bar %d: summed LE and LE vs ToT", runNumber, controlBar), 1400, 650);
  cControlBarLeToT->Divide(2, 1);
  cControlBarLeToT->cd(1);
  DrawPulseDistribution(hControlBarLe, false);
  cControlBarLeToT->cd(2);
  DrawPulseDistribution(hControlBarLeVsToT, true);
  cControlBarLeToT->Modified();
  cControlBarLeToT->Update();

  SaveCanvases("le_tot", outputDirectory, {cPixelLe,cControlPixelLe,cPixelLeVsToT,cControlPixelLeVsToT,cBarLeToT,cControlBarLeToT});
  std::cout << "[CDet vectors LE/ToT] selected detector pulses " << detectorPulses << '\n';
  analysisClock.Stop();
  std::cout << "[CDet vectors LE/ToT] " << vEntry.size() << " cached events; wall time "
            << analysisClock.RealTime() << " s including drawing/export; generation " << loadGeneration << '\n';
  return true;
}

bool Plot_CDet_PulseOccupancy_Vectors(
    const char* configFile = "CDet_run6077_pulse_pixel_timing_vectors.conf", int bar = -1,
    const char* outputTag = nullptr)
{
  TStopwatch analysisClock;
  TEnv config;
  if (!ReadCachedConfig(configFile, config)) return false;
  const int runNumber = config.GetValue("analysis.run_number", 0);
  const int pixelsPerBar = 16, nPixels = 2688;
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double energyMin = config.GetValue("analysis.ecal_energy_min", NAN);
  const double energyMax = config.GetValue("analysis.ecal_energy_max", NAN);
  if (!ValidRange(energyMin, energyMax)) { std::cerr << "[CDet vectors] Invalid ECal energy interval.\n"; return false; }
  const bool occupancyUseECalBaseline = config.GetValue("pulse_pixel_timing.occupancy_use_ecal_baseline",1)!=0;
  const double occupancyDxMaxM = config.GetValue("pulse_pixel_timing.occupancy_dx_max_m",0.08);
  if (!std::isfinite(occupancyDxMaxM) || occupancyDxMaxM<=0) {
    std::cerr << "[CDet vectors] Invalid occupancy x window.\n"; return false;
  }
  TString outputDirectory;
  if (!PrepareOutput(config, bar, outputTag, "occupancy", outputDirectory)) return false;
  const int firstPixel = pixelsPerBar*bar, firstControlPixel = pixelsPerBar*controlBar;
  delete hPixelOccupancyRaw;
  delete hPixelOccupancyECal;
  delete hPixelOccupancyXMatched;
  hPixelOccupancyRaw = new TH1D(TString::Format("pv_hPixelOccupancyRaw_run%d", runNumber), "All stored complete pulses;Logical pixel number;Events with #geq1 pulse", nPixels, -0.5, nPixels - 0.5);
  hPixelOccupancyECal = new TH1D(TString::Format("pv_hPixelOccupancyECal_run%d", runNumber), "ECal-selected pulses before x match;Logical pixel number;Events with #geq1 pulse", nPixels, -0.5, nPixels - 0.5);
  hPixelOccupancyXMatched = new TH1D(TString::Format("pv_hPixelOccupancyXMatched_run%d", runNumber), "ECal-selected pulses with x match;Logical pixel number;Events with #geq1 pulse", nPixels, -0.5, nPixels - 0.5);
  for (TH1D* hist : {hPixelOccupancyRaw, hPixelOccupancyECal, hPixelOccupancyXMatched}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue + 1);
  }
  hPixelOccupancyXMatched->SetLineColor(kRed + 1);

  for (size_t event=0; event<vEntry.size(); ++event) {
    const bool good_energy = GoodEnergy(event,energyMin,energyMax);
    std::array<bool,2688> seenRaw{}, seenECal{}, seenXMatched{};
    for (size_t i=0; i<vPulseID[event].size(); ++i) {
      if (!std::isfinite(vPulseID[event][i])) continue;
      const int id=int(vPulseID[event][i]);
      if (id<0 || id>=nPixels) continue;
      if (!seenRaw[id]) { hPixelOccupancyRaw->Fill(id); seenRaw[id]=true; }
      const bool good_ecal = good_energy && vPulseECalEligible[event][i]>0.5 && std::isfinite(vPulseDx[event][i]);
      if (good_ecal && !seenECal[id]) { hPixelOccupancyECal->Fill(id); seenECal[id]=true; }
      const bool good_x = good_ecal && std::fabs(vPulseDx[event][i])<=occupancyDxMaxM;
      if (good_x && !seenXMatched[id]) { hPixelOccupancyXMatched->Fill(id); seenXMatched[id]=true; }
    }
  }
  TCanvas *cPixelOccupancy = new TCanvas(TString::Format("pv_cPixelOccupancy_run%d_bar%d_control%d", runNumber, bar, controlBar), "CDet event occupancy before and after ECal x match", 1600, 1000);
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

  SaveCanvases("occupancy", outputDirectory, {cPixelOccupancy});
  analysisClock.Stop();
  std::cout << "[CDet vectors occupancy] " << vEntry.size() << " cached events; wall time "
            << analysisClock.RealTime() << " s including drawing/export; generation " << loadGeneration << '\n';
  return true;
}

bool Plot_CDet_PulseProjectionResponse_Vectors(
    const char* configFile = "CDet_run6077_pulse_pixel_timing_vectors.conf", int bar = -1,
    const char* outputTag = nullptr)
{
  TStopwatch analysisClock;
  TEnv config;
  if (!ReadCachedConfig(configFile, config)) return false;
  const int runNumber = config.GetValue("analysis.run_number", 0);
  const int pixelsPerBar = 16, nPixels = 2688;
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double energyMin = config.GetValue("analysis.ecal_energy_min", NAN);
  const double energyMax = config.GetValue("analysis.ecal_energy_max", NAN);
  if (!ValidRange(energyMin, energyMax)) { std::cerr << "[CDet vectors] Invalid ECal energy interval.\n"; return false; }
  const double leCutMinNs = config.GetValue("analysis.le_min", 0.02);
  const double leCutMaxNs = config.GetValue("analysis.le_max", 60.0);
  const bool applyToTCut = config.GetValue("pulse_pixel_timing.apply_tot_cut", 0) != 0;
  const double totCutMinNs = config.GetValue("analysis.tot_min", 4.0);
  const double totCutMaxNs = config.GetValue("analysis.tot_max", 30.0);
  const double ellipseRadius = config.GetValue("good_pulse.single_cut_radius", -1.0);
  const bool useEllipse = ellipseRadius > 0.0;
  const double ellipseXCenter = config.GetValue("good_pulse.single_residual_center_m", 0.0);
  const double ellipseTimeCenter = config.GetValue("good_pulse.single_timing_center_ns", -26.0);
  const double ellipseXScale = config.GetValue("good_pulse.single_residual_scale_m", 0.040);
  const double ellipseTimeScale = config.GetValue("good_pulse.single_timing_scale_ns", 5.0);
  if (!std::isfinite(leCutMinNs) || !std::isfinite(leCutMaxNs) || leCutMinNs >= leCutMaxNs ||
      (applyToTCut && (!std::isfinite(totCutMinNs) || !std::isfinite(totCutMaxNs) ||
       totCutMinNs < 0.0 || totCutMinNs >= totCutMaxNs))) {
    std::cerr << "[CDet pulse pixel timing] Invalid analysis.le_min/max or enabled analysis.tot_min/max (ns).\n";
    return false;
  }
  if (!std::isfinite(ellipseRadius) || (!useEllipse && ellipseRadius != -1.0) ||
      (useEllipse && (!std::isfinite(ellipseXCenter) || !std::isfinite(ellipseTimeCenter) ||
       !std::isfinite(ellipseXScale) || !std::isfinite(ellipseTimeScale) ||
       ellipseXScale <= 0.0 || ellipseTimeScale <= 0.0))) {
    std::cerr << "[CDet pulse pixel timing] Ellipse radius must be -1 (off) or positive; enabled ellipses require finite centers and positive finite scales.\n";
    return false;
  }
  const double ecalTimeMinNs=config.GetValue("analysis.ecal_time_min",NAN);
  const double ecalTimeMaxNs=config.GetValue("analysis.ecal_time_max",NAN);
  const double projectionECalZM = config.GetValue("projection.ecal_z_m", NAN);
  const TString geometryTimestamp = config.GetValue("projection.database_timestamp", "");
  TString geometryFile = config.GetValue("projection.database_file", "../../DB/db_earm.cdet.dat");
  if (!gSystem->IsAbsoluteFileName(geometryFile))
    geometryFile = TString::Format("%s/%s", gSystem->DirName(configFile), geometryFile.Data());
  if (!std::isfinite(ecalTimeMinNs) || !std::isfinite(ecalTimeMaxNs) || ecalTimeMinNs >= ecalTimeMaxNs ||
      !std::isfinite(projectionECalZM) || projectionECalZM <= 0.0) {
    std::cerr << "[CDet projection] Set valid analysis.ecal_time_min/max and projection.ecal_z_m.\n";
    return false;
  }
  std::map<std::string, std::vector<double>> geometryDB;
  if (!ReadProjectionDatabase(geometryFile, geometryTimestamp, geometryDB)) return false;
  for (const auto& item : geometryDB) {
    for (double value : item.second) {
      if (!std::isfinite(value)) {
        std::cerr << "[CDet projection] Nonfinite database value: " << item.first << '\n';
        return false;
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
    return false;
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
      return false;
    }
    geometryX[id] = xScale * x - xAlignment;
    geometryY[id] = y;
    geometryZ[id] = z;
    instrumented[id] = true;
    ++instrumentedPerLayer[layer];
  }
  if (instrumentedPerLayer[0] == 0 || instrumentedPerLayer[1] == 0) {
    std::cerr << "[CDet projection] No instrumented geometry in one or both layers.\n";
    return false;
  }
  TString outputDirectory;
  if (!PrepareOutput(config, bar, outputTag, "projection", outputDirectory)) return false;
  const int firstPixel = pixelsPerBar*bar, firstControlPixel = pixelsPerBar*controlBar;
  const double ecalXMinM=-1.5, ecalXMaxM=1.5, ecalXBinWidthM=0.0425;
  const int nECalX=int(std::ceil((ecalXMaxM-ecalXMinM)/ecalXBinWidthM));
  delete hECalX;
  delete hECalXReference;
  hECalX = new TH1D(TString::Format("pv_hECalX_run%d", runNumber), TString::Format("Run %d: event with a final-selected CDet pulse;ECal x (m);Events", runNumber), nECalX, ecalXMinM, ecalXMaxM);
  hECalXReference = new TH1D(TString::Format("pv_hECalXReference_run%d", runNumber), TString::Format("Run %d: ECal energy/time reference, no CDet requirement;ECal x (m);Events", runNumber), nECalX, ecalXMinM, ecalXMaxM);
  for (TH1D* hist : {hECalX, hECalXReference}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
  }

  delete hProjectionOpportunities;
  delete hProjectionComplete;
  delete hProjectionSelected;
  delete hProjectionInstrumented;
  hProjectionOpportunities = new TH1D(TString::Format("pv_hProjectionOpportunities_run%d", runNumber), "ECal projection opportunities;Logical pixel number;Events", nPixels, -0.5, nPixels - 0.5);
  hProjectionComplete = new TH1D(TString::Format("pv_hProjectionComplete_run%d", runNumber), "Projection opportunity with a complete pulse;Logical pixel number;Events", nPixels, -0.5, nPixels - 0.5);
  hProjectionSelected = new TH1D(TString::Format("pv_hProjectionSelected_run%d", runNumber), "Projection opportunity with a final-selected pulse;Logical pixel number;Events", nPixels, -0.5, nPixels - 0.5);
  hProjectionInstrumented = new TH1D(TString::Format("pv_hProjectionInstrumented_run%d", runNumber), "DB instrumented pixel mask;Logical pixel number;Instrumented", nPixels, -0.5, nPixels - 0.5);
  for (TH1D* hist : {hProjectionOpportunities, hProjectionComplete, hProjectionSelected, hProjectionInstrumented}) {
    hist->SetDirectory(nullptr);
    hist->SetStats(false);
    hist->SetLineWidth(2);
  }
  hProjectionOpportunities->SetLineColor(kBlack);
  hProjectionComplete->SetLineColor(kBlue + 1);
  hProjectionSelected->SetLineColor(kRed + 1);
  for (int id = 0; id < nPixels; ++id) hProjectionInstrumented->SetBinContent(id + 1, instrumented[id]);

  Long64_t referenceEvents=0, referenceWithoutPulses=0, referenceWithoutSelectedPulses=0;
  for (size_t event=0; event<vEntry.size(); ++event) {
    const bool good_ECal_reference = GoodEnergy(event,energyMin,energyMax) &&
        std::isfinite(vECalTime[event]) && vECalTime[event]>ecalTimeMinNs && vECalTime[event]<ecalTimeMaxNs &&
        std::isfinite(vECalNclus[event]) && vECalNclus[event]>0 &&
        std::isfinite(vECalX[event]) && std::isfinite(vECalY[event]);
    bool eventHasGoodPulse=false;
    std::array<bool,2688> seenRaw{}, seenSelected{};
    for (size_t i=0; i<vPulseID[event].size(); ++i) {
      const bool good_id=std::isfinite(vPulseID[event][i]);
      const int id=good_id ? int(vPulseID[event][i]) : -1;
      if (id>=0 && id<nPixels) seenRaw[id]=true;
      const bool good_final = GoodBeforeEllipse(event, i, energyMin, energyMax, leCutMinNs, leCutMaxNs, applyToTCut, totCutMinNs, totCutMaxNs) && GoodEllipse(event, i, ellipseRadius, ellipseXCenter, ellipseTimeCenter, ellipseXScale, ellipseTimeScale);
      if (good_final) {
        eventHasGoodPulse=true;
        if (id>=0 && id<nPixels) seenSelected[id]=true;
      }
    }
    if (eventHasGoodPulse) hECalX->Fill(vECalX[event]);
    if (!good_ECal_reference) continue;
    ++referenceEvents;
    hECalXReference->Fill(vECalX[event]);
    if (vPulseID[event].empty()) ++referenceWithoutPulses;
    if (!eventHasGoodPulse) ++referenceWithoutSelectedPulses;
    for (int id=0; id<nPixels; ++id) {
      if (!instrumented[id]) continue;
      const double dx=geometryX[id]-vECalX[event]*geometryZ[id]/projectionECalZM-projectionXOffsetM;
      const double dy=geometryY[id]-vECalY[event]*geometryZ[id]/projectionECalZM-projectionYOffsetM;
      const bool good_opportunity=std::fabs(dx)<=projectionDxMaxM && std::fabs(dy)<=projectionDyMaxM;
      if (good_opportunity) {
        hProjectionOpportunities->Fill(id);
        if (seenRaw[id]) hProjectionComplete->Fill(id);
        if (seenSelected[id]) hProjectionSelected->Fill(id);
      }
    }
  }
  // Binomial intervals for each conditional fraction; bins with no denominator
  // are omitted from the graphs, rather than displayed as zero response.
  delete gProjectionCompleteFraction;
  delete gProjectionSelectedFraction;
  gProjectionCompleteFraction = new TGraphAsymmErrors;
  gProjectionSelectedFraction = new TGraphAsymmErrors;
  gProjectionCompleteFraction->SetName(TString::Format("pv_gProjectionCompleteFraction_run%d", runNumber));
  gProjectionSelectedFraction->SetName(TString::Format("pv_gProjectionSelectedFraction_run%d", runNumber));
  gProjectionCompleteFraction->Divide(hProjectionComplete, hProjectionOpportunities, "cl=0.683 cp");
  gProjectionSelectedFraction->Divide(hProjectionSelected, hProjectionOpportunities, "cl=0.683 cp");
  for (TGraphAsymmErrors* graph : {gProjectionCompleteFraction, gProjectionSelectedFraction}) {
    const int color = graph == gProjectionCompleteFraction ? kBlue + 1 : kRed + 1;
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(graph == gProjectionCompleteFraction ? 20 : 24);
    graph->SetMarkerSize(0.5);
  }
  TCanvas *cProjectionCounts = new TCanvas(TString::Format("pv_cProjectionCounts_run%d_bar%d_control%d", runNumber, bar, controlBar), "ECal projection opportunities and CDet responses", 1600, 1000);
  TCanvas *cProjectionFractions = new TCanvas(TString::Format("pv_cProjectionFractions_run%d_bar%d_control%d", runNumber, bar, controlBar), "Conditional CDet response fractions", 1600, 1000);
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
  TCanvas *cECalXSelection = new TCanvas(TString::Format("pv_cECalXSelection_run%d", runNumber), "ECal x: independent reference and final CDet selection", 1400, 650);
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

  SaveCanvases("projection", outputDirectory, {cProjectionCounts,cProjectionFractions,cECalXSelection});
  std::cout << "[CDet vectors projection] " << referenceEvents << " reference events; "
            << referenceWithoutPulses << " without pulses; " << referenceWithoutSelectedPulses << " without selected pulses.\n";
  analysisClock.Stop();
  std::cout << "[CDet vectors projection] " << vEntry.size() << " cached events; wall time "
            << analysisClock.RealTime() << " s including drawing/export; generation " << loadGeneration << '\n';
  return true;
}

// Full-suite entry point: first call loads; later calls reuse the same dataset.
// Dataset changes are never implicitly reloaded. Use Load_... explicitly.
bool Plot_CDet_PulsePixelTiming_Vectors(
    const char* configFile = "CDet_run6077_pulse_pixel_timing_vectors.conf", int bar = -1,
    const char* inputDirectory = nullptr, Long64_t eventsOverride = -2,
    double ecalDTBinWidthNs = NAN, double ecalDTMinNs = NAN, double ecalDTMaxNs = NAN,
    double ecalDXBinWidthM = NAN, double ecalDXMinM = NAN, double ecalDXMaxM = NAN,
    double dtVsDxMinNs = NAN, double dtVsDxMaxNs = NAN,
    double leBinWidthNs = NAN, double leMinNs = NAN, double leMaxNs = NAN,
    double totBinWidthNs = NAN, double totMinNs = NAN, double totMaxNs = NAN,
    const char* outputTag = nullptr)
{
  if (!cacheReady && !Load_CDet_PulsePixelTiming_Vectors(configFile,inputDirectory,eventsOverride)) return false;
  TEnv config;
  if (!ReadCachedConfig(configFile,config)) return false;
  TString directory=inputDirectory && inputDirectory[0] ? inputDirectory : loadedDirectory.c_str();
  gSystem->ExpandPathName(directory);
  if (!gSystem->IsAbsoluteFileName(directory)) directory=TString(gSystem->WorkingDirectory())+"/"+directory;
  if (directory!=loadedDirectory.c_str() || (eventsOverride!=-2 && eventsOverride!=loadedEventLimit)) {
    std::cerr << "[CDet vectors] Dataset override differs from cache; explicitly reload first.\n"; return false;
  }
  return Plot_CDet_PulseTiming_Vectors(configFile,bar,ecalDTBinWidthNs,ecalDTMinNs,ecalDTMaxNs,
             ecalDXBinWidthM,ecalDXMinM,ecalDXMaxM,dtVsDxMinNs,dtVsDxMaxNs,outputTag) &&
      Plot_CDet_PulseLEToT_Vectors(configFile,bar,leBinWidthNs,leMinNs,leMaxNs,totBinWidthNs,totMinNs,totMaxNs,outputTag) &&
      Plot_CDet_PulseOccupancy_Vectors(configFile,bar,outputTag) &&
      Plot_CDet_PulseProjectionResponse_Vectors(configFile,bar,outputTag);
}

bool Plot_CDet_PulsePixelTiming_Vectors(
    int runNumber, int bar=-1, const char* inputDirectory=nullptr, Long64_t eventsOverride=-2,
    double ecalDTBinWidthNs=NAN, double ecalDTMinNs=NAN, double ecalDTMaxNs=NAN,
    double ecalDXBinWidthM=NAN, double ecalDXMinM=NAN, double ecalDXMaxM=NAN,
    double dtVsDxMinNs=NAN, double dtVsDxMaxNs=NAN,
    double leBinWidthNs=NAN, double leMinNs=NAN, double leMaxNs=NAN,
    double totBinWidthNs=NAN, double totMinNs=NAN, double totMaxNs=NAN,
    const char* outputTag=nullptr)
{
  const TString config=TString::Format("CDet_run%d_pulse_pixel_timing_vectors.conf",runNumber);
  return Plot_CDet_PulsePixelTiming_Vectors(config,bar,inputDirectory,eventsOverride,
      ecalDTBinWidthNs,ecalDTMinNs,ecalDTMaxNs,ecalDXBinWidthM,ecalDXMinM,ecalDXMaxM,
      dtVsDxMinNs,dtVsDxMaxNs,leBinWidthNs,leMinNs,leMaxNs,totBinWidthNs,totMinNs,totMaxNs,outputTag);
}

bool Plot_CDet_PulsePixelTiming_Vectors(const char* configFile, const char* outputTag,
    int bar=-1, const char* inputDirectory=nullptr, Long64_t eventsOverride=-2)
{
  return Plot_CDet_PulsePixelTiming_Vectors(configFile,bar,inputDirectory,eventsOverride,
      NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,outputTag);
}

bool Plot_CDet_PulsePixelTiming_Vectors(int runNumber, const char* outputTag,
    int bar=-1, const char* inputDirectory=nullptr, Long64_t eventsOverride=-2)
{
  const TString config=TString::Format("CDet_run%d_pulse_pixel_timing_vectors.conf",runNumber);
  return Plot_CDet_PulsePixelTiming_Vectors(config,outputTag,bar,inputDirectory,eventsOverride);
}

} // namespace CDetPulseVectors

using CDetPulseVectors::Load_CDet_PulsePixelTiming_Vectors;
using CDetPulseVectors::Plot_CDet_PulseTiming_Vectors;
using CDetPulseVectors::Plot_CDet_PulseLEToT_Vectors;
using CDetPulseVectors::Plot_CDet_PulseOccupancy_Vectors;
using CDetPulseVectors::Plot_CDet_PulseProjectionResponse_Vectors;
using CDetPulseVectors::Plot_CDet_PulsePixelTiming_Vectors;
