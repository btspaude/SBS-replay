// Synthetic zero-pulse, duplicate-pixel, NaN and extra-type fixtures.
// Load this driver alone in a fresh ROOT session (it includes the parity driver).
#include "Validate_CDet_PulsePixelTiming_Vectors.C"

bool Validate_CDet_PulseVectors_Fixture(const char* sourceConfig, const char* geometryFile,
    const char* directory)
{
  using namespace CDetVectorValidation;
  try {
    const std::string root=directory;
    gSystem->mkdir((root+"/input").c_str(),true);
    const std::string file=root+"/input/cdet_6077_stream0_2_seg0_0.root";
    {
      TFile output((root+"/input/fixture_source.root").c_str(),"RECREATE"); TTree tree("T","Synthetic pulse schema");
      Int_t count=0;
      tree.Branch("pulseCount",&count,"pulseCount/I");
      std::map<std::string,std::array<Double_t,8>> arrays;
      for(const char* suffix : {"pmtnum","index","le_index","te_index","row","col","layer","x","y","z","tdc_le","tdc_te","tdc_tot","tdc_le_raw","tdc_te_raw","tdc_tot_raw","tdc_le_corr","tdc_te_corr","tdc_tot_ns","ecal_residual","calib_valid","x_corr","ecal_x_proj","ecal_y_proj","ecal_x_residual","ecal_y_residual","broad_quality_pass","ecal_eligible","spatial_pass"}) {
        const std::string name=std::string("earm.cdet.pulse.")+suffix;
        tree.Branch(name.c_str(),arrays[suffix].data(),(name+"[pulseCount]/D").c_str());
      }
      std::map<std::string,Double_t> scalars;
      for(const char* suffix : {"n","n_le_only","n_te_only","n_invalid_pair"}) {
        const std::string name=std::string("earm.cdet.pulse.")+suffix;
        tree.Branch(name.c_str(),&scalars[suffix],(name+"/D").c_str());
      }
      for(const char* suffix : {"x","y","adctime","nclus","e"}) {
        const std::string name=std::string("earm.ecal.")+suffix;
        tree.Branch(name.c_str(),&scalars[suffix],(name+"/D").c_str());
      }
      ULong64_t wide=9007199254740993ULL;
      std::array<Float_t,8> floats{};
      tree.Branch("earm.cdet.pulse.extra_wide",&wide,"earm.cdet.pulse.extra_wide/l");
      tree.Branch("earm.cdet.pulse.extra_float",floats.data(),"earm.cdet.pulse.extra_float[pulseCount]/F");
      for(int e=0;e<10;++e) {
        count=e==0 ? 0 : 2;
        for(auto& field:arrays) field.second.fill(0);
        scalars["n"]=count; scalars["n_le_only"]=e; scalars["n_te_only"]=e+1; scalars["n_invalid_pair"]=e+2;
        scalars["x"]=.1; scalars["y"]=.2; scalars["adctime"]=0; scalars["e"]=3.5; scalars["nclus"]=1;
        if(e==3) scalars["e"]=3.0; // Strict boundary: rejected.
        if(e==7) scalars["adctime"]=10; // Reference rejects; stored eligible remains set.
        if(e==8) scalars["nclus"]=0;
        for(int p=0;p<count;++p) {
          arrays["pmtnum"][p]=e==2 ? 480 : 304;
          arrays["index"][p]=p;
          arrays["tdc_le"][p]=p==0 ? 2 : 6000; // Inclusive LE endpoints.
          arrays["tdc_tot_ns"][p]=p==0 ? 4 : 40;
          arrays["ecal_residual"][p]=e==4 ? (p==0 ? -150 : 50) : -26;
          arrays["ecal_x_residual"][p]=e==4 ? .2 : 0;
          arrays["calib_valid"][p]=1;
          arrays["ecal_eligible"][p]=1;
          arrays["spatial_pass"][p]=1;
          if(e==5) arrays["ecal_residual"][p]=NAN;
          floats[p]=p==0 ? .25f : NAN;
        }
        ++wide; tree.Fill();
      }
      tree.Write();
    }
    // Split the ten events across a base file and its rollover part, so the
    // same comparisons also exercise TTreeReader's chain transition.
    {
      TFile input((root+"/input/fixture_source.root").c_str(),"READ");
      auto* source=input.Get<TTree>("T");
      for(int part=0;part<2;++part) {
        const std::string path=part==0 ? file : root+"/input/cdet_6077_stream0_2_seg0_0_1.root";
        TFile output(path.c_str(),"RECREATE");
        auto* copy=source->CopyTree("", "", 5, 5*part); copy->Write();
      }
    }
    Require(Validate_CDet_PulsePixelTiming_Vectors(sourceConfig,(root+"/input").c_str(),geometryFile,(root+"/plots").c_str(),10),"fixture parity suite");
    const std::string config=root+"/plots/validation.conf";
    Require(Load_CDet_PulsePixelTiming_Vectors(config.c_str(),(root+"/input").c_str(),10),"fixture reload");
    Require(CDetPulseVectors::vPulseID[0].empty() && CDetPulseVectors::vECalEnergy[0]==3.5,"zero-pulse event retained");
    Require(CDetPulseVectors::pulseScalars<ULong64_t>.at("earm.cdet.pulse.extra_wide")[0]==9007199254740994ULL &&
            CDetPulseVectors::pulseScalars<ULong64_t>.at("earm.cdet.pulse.extra_wide")[1]==9007199254740995ULL,"extra integer precision");
    Require(CDetPulseVectors::pulseArrays<Float_t>.at("earm.cdet.pulse.extra_float")[1][0]==.25f &&
            std::isnan(CDetPulseVectors::pulseArrays<Float_t>.at("earm.cdet.pulse.extra_float")[1][1]),"extra float NaN");
    const TString cwd=gSystem->WorkingDirectory(); gSystem->ChangeDirectory((root+"/plots").c_str());
    Require(Plot_CDet_PulseProjectionResponse_Vectors(config.c_str()),"zero-pulse reference analysis");
    Require(CDetPulseVectors::hECalXReference->GetEntries()==7,"reference includes zero-pulse event");
    Require(Plot_CDet_PulseOccupancy_Vectors(config.c_str()),"duplicate occupancy analysis");
    Require(CDetPulseVectors::hPixelOccupancyRaw->GetBinContent(305)==8,"duplicate pulses count once per pixel/event");
    gSystem->ChangeDirectory(cwd);
    // Omit one additional branch in a second file; reject inconsistent schemas.
    {
      TFile out((root+"/input/cdet_6077_stream0_2_seg1_1.root").c_str(),"RECREATE");
      TFile input(file.c_str(),"READ"); auto* source=input.Get<TTree>("T");
      source->SetBranchStatus("earm.cdet.pulse.extra_float",0);
      out.cd(); auto* copy=source->CloneTree(); copy->Write();
    }
    std::ifstream source(sourceConfig); const std::string configText((std::istreambuf_iterator<char>(source)),{});
    WriteConfig(configText,config,geometryFile,10);
    std::ifstream configIn(config); std::string changed((std::istreambuf_iterator<char>(configIn)),{}); configIn.close();
    const auto pos=changed.find("analysis.max_segment: 0"); changed.replace(pos,std::strlen("analysis.max_segment: 0"),"analysis.max_segment: 1");
    std::ofstream(config)<<changed;
    Require(!Load_CDet_PulsePixelTiming_Vectors(config.c_str(),(root+"/input").c_str()),"inconsistent file schemas rejected");
    Require(!CDetPulseVectors::cacheReady && CDetPulseVectors::vEntry.empty(),"schema failure clears cache");
    std::cout<<"SYNTHETIC VALIDATION PASSED (zero pulses, duplicates, cut boundaries, flow bins, NaNs, extra types, rollover transition, schema mismatch)\n";
    return true;
  } catch(const std::exception& error) {std::cerr<<error.what()<<'\n'; return false;}
}
