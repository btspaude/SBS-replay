// Regression driver: runs the unchanged original and the vector copy together.
// Use absolute paths; all generated plots/configs go in outputDirectory.
#include "Plot_CDet_PulsePixelTiming.C"
#include "Plot_CDet_PulsePixelTiming_Vectors.C"
#include <cstring>

namespace CDetVectorValidation {
void Require(bool pass, const std::string& what)
{
  if (!pass) throw std::runtime_error("VALIDATION FAILED: " + what);
}
void Compare(TH1* a, TH1* b)
{
  Require(bool(a)==bool(b), "histogram availability");
  if (!a) return;
  const std::string name=a->GetName();
  Require(a->GetNcells()==b->GetNcells() && a->GetDimension()==b->GetDimension(), name+" dimensions");
  Require(a->GetEntries()==b->GetEntries(), name+" entries");
  for (auto axes : {std::make_pair(a->GetXaxis(),b->GetXaxis()),std::make_pair(a->GetYaxis(),b->GetYaxis())}) {
    Require(axes.first->GetNbins()==axes.second->GetNbins(),name+" axis bins");
    for (int i=1;i<=axes.first->GetNbins()+1;++i)
      Require(axes.first->GetBinLowEdge(i)==axes.second->GetBinLowEdge(i),name+" axis edges");
  }
  for (int bin=0;bin<a->GetNcells();++bin) {
    Require(a->GetBinContent(bin)==b->GetBinContent(bin),name+" content bin "+std::to_string(bin));
    Require(a->GetBinError(bin)==b->GetBinError(bin),name+" error bin "+std::to_string(bin));
  }
}
std::vector<std::pair<TH1*,TH1*>> Pairs()
{
  std::vector<std::pair<TH1*,TH1*>> pairs;
  for (int i=0;i<16;++i) pairs.emplace_back(::hPixelDt[i], CDetPulseVectors::hPixelDt[i]);
  for (int i=0;i<16;++i) pairs.emplace_back(::hControlPixelDt[i], CDetPulseVectors::hControlPixelDt[i]);
  pairs.emplace_back(::hBarDt, CDetPulseVectors::hBarDt);
  pairs.emplace_back(::hControlBarDt, CDetPulseVectors::hControlBarDt);
  pairs.emplace_back(::hBarDtVsDxBefore, CDetPulseVectors::hBarDtVsDxBefore);
  pairs.emplace_back(::hBarDtVsDxAfter, CDetPulseVectors::hBarDtVsDxAfter);
  pairs.emplace_back(::hControlBarDtVsDxBefore, CDetPulseVectors::hControlBarDtVsDxBefore);
  pairs.emplace_back(::hControlBarDtVsDxAfter, CDetPulseVectors::hControlBarDtVsDxAfter);
  pairs.emplace_back(::hDetectorDtVsDxBefore, CDetPulseVectors::hDetectorDtVsDxBefore);
  pairs.emplace_back(::hDetectorDtVsDxAfter, CDetPulseVectors::hDetectorDtVsDxAfter);
  for (int i=0;i<16;++i) pairs.emplace_back(::hPixelLe[i], CDetPulseVectors::hPixelLe[i]);
  for (int i=0;i<16;++i) pairs.emplace_back(::hControlPixelLe[i], CDetPulseVectors::hControlPixelLe[i]);
  for (int i=0;i<16;++i) pairs.emplace_back(::hPixelLeVsToT[i], CDetPulseVectors::hPixelLeVsToT[i]);
  for (int i=0;i<16;++i) pairs.emplace_back(::hControlPixelLeVsToT[i], CDetPulseVectors::hControlPixelLeVsToT[i]);
  pairs.emplace_back(::hBarLe, CDetPulseVectors::hBarLe);
  pairs.emplace_back(::hControlBarLe, CDetPulseVectors::hControlBarLe);
  pairs.emplace_back(::hBarLeVsToT, CDetPulseVectors::hBarLeVsToT);
  pairs.emplace_back(::hControlBarLeVsToT, CDetPulseVectors::hControlBarLeVsToT);
  pairs.emplace_back(::hPixelOccupancyRaw, CDetPulseVectors::hPixelOccupancyRaw);
  pairs.emplace_back(::hPixelOccupancyECal, CDetPulseVectors::hPixelOccupancyECal);
  pairs.emplace_back(::hPixelOccupancyXMatched, CDetPulseVectors::hPixelOccupancyXMatched);
  pairs.emplace_back(::hECalX, CDetPulseVectors::hECalX);
  pairs.emplace_back(::hECalXReference, CDetPulseVectors::hECalXReference);
  pairs.emplace_back(::hProjectionOpportunities, CDetPulseVectors::hProjectionOpportunities);
  pairs.emplace_back(::hProjectionComplete, CDetPulseVectors::hProjectionComplete);
  pairs.emplace_back(::hProjectionSelected, CDetPulseVectors::hProjectionSelected);
  pairs.emplace_back(::hProjectionInstrumented, CDetPulseVectors::hProjectionInstrumented);

  return pairs;
}
void CompareAll()
{
  int count=0;
  for (const auto& pair:Pairs()) { Compare(pair.first,pair.second); if(pair.first) ++count; }
  for (auto graphs : {std::make_pair(::gProjectionCompleteFraction,CDetPulseVectors::gProjectionCompleteFraction),
                     std::make_pair(::gProjectionSelectedFraction,CDetPulseVectors::gProjectionSelectedFraction)}) {
    Require(graphs.first->GetN()==graphs.second->GetN(),"fraction graph size");
    for (int i=0;i<graphs.first->GetN();++i) {
      Require(graphs.first->GetPointX(i)==graphs.second->GetPointX(i) && graphs.first->GetPointY(i)==graphs.second->GetPointY(i),"fraction graph points");
      Require(graphs.first->GetErrorYlow(i)==graphs.second->GetErrorYlow(i) && graphs.first->GetErrorYhigh(i)==graphs.second->GetErrorYhigh(i),"fraction graph intervals");
    }
  }
  std::cout<<"VALIDATED "<<count<<" histograms (all cells/errors/entries/axes) and two fraction graphs\n";
}
void CheckStoredValues()
{
  using namespace CDetPulseVectors;
  TChain chain("T"); for (const auto& file:loadedFiles) chain.Add(file.c_str());
  TTreeReader reader(&chain);
  std::map<std::string,std::shared_ptr<TTreeReaderArray<Double_t>>> arrays;
  std::map<std::string,std::shared_ptr<TTreeReaderValue<Double_t>>> scalars;
  for (const auto& field:loadedSchema) {
    if(field.second.find("Double_t:")!=0) continue;
    if(field.second=="Double_t:scalar") scalars[field.first]=std::make_shared<TTreeReaderValue<Double_t>>(reader,field.first.c_str());
    else arrays[field.first]=std::make_shared<TTreeReaderArray<Double_t>>(reader,field.first.c_str());
  }
  size_t nanCount=0, zeroCount=0;
  auto same=[&](double a,double b) {if(std::isnan(a)) ++nanCount; return std::memcmp(&a,&b,sizeof(a))==0;};
  for(size_t e=0;e<vEntry.size();++e) {
    Require(reader.Next() && vEntry[e]==reader.GetCurrentEntry(),"chain entry order");
    if(vPulseID[e].empty()) ++zeroCount;
    for(const auto& field:arrays) {
      const auto& values=pulseArrays<Double_t>.at(field.first).at(e);
      Require(values.size()==field.second->GetSize(),field.first+" size");
      for(size_t p=0;p<values.size();++p) Require(same(values[p],(*field.second)[p]),field.first+" exact value");
    }
    for(const auto& field:scalars) {
      double value=0;
      if(field.first=="earm.ecal.x") value=vECalX[e];
      else if(field.first=="earm.ecal.y") value=vECalY[e];
      else if(field.first=="earm.ecal.adctime") value=vECalTime[e];
      else if(field.first=="earm.ecal.nclus") value=vECalNclus[e];
      else if(field.first=="earm.ecal.e") value=vECalEnergy[e];
      else value=pulseScalars<Double_t>.at(field.first).at(e);
      Require(same(value,**field.second),field.first+" exact scalar");
    }
  }
  std::cout<<"VALIDATED exact Double_t cache values, including "<<nanCount<<" NaNs and "<<zeroCount<<" zero-pulse events\n";
}
void WriteConfig(const std::string& source,const std::string& path,const char* geometry,Long64_t events,bool changed=false)
{
  std::istringstream lines(source); std::ofstream out(path); std::string line;
  // Replace keys, rather than duplicate TEnv definitions.
  while(std::getline(lines,line)) {
    if(line.find("analysis.events:")==0) line="analysis.events: "+std::to_string(events);
    if(line.find("analysis.min_segment:")==0) line="analysis.min_segment: 0";
    if(line.find("analysis.max_segment:")==0) line="analysis.max_segment: 0";
    if(line.find("projection.database_file:")==0) line=std::string("projection.database_file: ")+geometry;
    if(changed && line.find("good_pulse.single_cut_radius:")==0) line="good_pulse.single_cut_radius: 2";
    if(changed && line.find("pulse_pixel_timing.apply_tot_cut:")==0) line="pulse_pixel_timing.apply_tot_cut: 1";
    if(changed && line.find("analysis.le_max:")==0) line="analysis.le_max: 45";
    if(changed && line.find("analysis.ecal_energy_min:")==0) line="analysis.ecal_energy_min: 3.2";
    out<<line<<'\n';
  }
}
}

bool Validate_CDet_PulsePixelTiming_Vectors(const char* sourceConfig,const char* inputDirectory,
    const char* geometryFile,const char* outputDirectory,Long64_t events=1000)
{
  using namespace CDetVectorValidation;
  const TString originalDirectory=gSystem->WorkingDirectory();
  try {
    Require(gSystem->IsAbsoluteFileName(sourceConfig) && gSystem->IsAbsoluteFileName(inputDirectory) &&
            gSystem->IsAbsoluteFileName(geometryFile) && gSystem->IsAbsoluteFileName(outputDirectory),"use absolute paths");
    std::ifstream in(sourceConfig); Require(bool(in),"source config readable");
    const std::string source((std::istreambuf_iterator<char>(in)),{});
    gSystem->mkdir(outputDirectory,true); Require(gSystem->ChangeDirectory(outputDirectory),"output directory");
    const std::string config=std::string(outputDirectory)+"/validation.conf";
    WriteConfig(source,config,geometryFile,events);
    TStopwatch watch;
    Plot_CDet_PulsePixelTiming(config.c_str(),-1,inputDirectory);
    watch.Stop(); const double originalTime=watch.RealTime();
    Require(Load_CDet_PulsePixelTiming_Vectors(config.c_str(),inputDirectory),"initial cache load");
    CheckStoredValues();
    const auto generation=CDetPulseVectors::loadGeneration;
    watch.Start(); Require(Plot_CDet_PulsePixelTiming_Vectors(config.c_str()),"full vector suite"); watch.Stop();
    const double firstTime=watch.RealTime(); CompareAll();
    watch.Start(); Require(Plot_CDet_PulsePixelTiming_Vectors(config.c_str()),"repeat suite"); watch.Stop();
    const double repeatTime=watch.RealTime(); CompareAll();
    Require(CDetPulseVectors::loadGeneration==generation,"no reload on repeat");
    auto* projection=CDetPulseVectors::hProjectionOpportunities;
    auto* occupancy=CDetPulseVectors::hPixelOccupancyRaw;
    const int canvases=gROOT->GetListOfCanvases()->GetSize();
    Require(Plot_CDet_PulseTiming_Vectors(config.c_str()),"independent timing call");
    Require(projection==CDetPulseVectors::hProjectionOpportunities && occupancy==CDetPulseVectors::hPixelOccupancyRaw &&
            canvases==gROOT->GetListOfCanvases()->GetSize(),"timing does not rerun other families or accumulate canvases");
    CompareAll();
    WriteConfig(source,config,geometryFile,events,true);
    Plot_CDet_PulsePixelTiming(config.c_str(),-1,inputDirectory);
    Require(Plot_CDet_PulsePixelTiming_Vectors(config.c_str()),"changed cuts suite"); CompareAll();
    Require(CDetPulseVectors::loadGeneration==generation,"cut edit must reuse cache");
    WriteConfig(source,config,geometryFile,std::max<Long64_t>(1,events/2));
    Require(!Plot_CDet_PulseTiming_Vectors(config.c_str()),"dataset edit rejected");
    Require(Load_CDet_PulsePixelTiming_Vectors(config.c_str(),inputDirectory),"explicit reload");
    Require(CDetPulseVectors::vEntry.size()==size_t(std::max<Long64_t>(1,events/2)),"reload replaces cache");
    Require(CDetPulseVectors::loadGeneration==generation+1,"reload generation");
    Require(!Load_CDet_PulsePixelTiming_Vectors(config.c_str(),"/nonexistent/cdet_validation"),"missing dataset fails");
    Require(!CDetPulseVectors::cacheReady && CDetPulseVectors::vEntry.empty(),"failure clears cache");
    std::cout<<"VALIDATION PASSED; wall seconds original full="<<originalTime<<", cached full="<<firstTime
             <<", cached repeat="<<repeatTime<<" (all include fits/drawing/export; not a controlled benchmark)\n";
    gSystem->ChangeDirectory(originalDirectory); return true;
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n'; gSystem->ChangeDirectory(originalDirectory); return false;
  }
}
