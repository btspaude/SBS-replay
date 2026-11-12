// Load once, then plot saved event-indexed vectors, like Plot_CDet_PulseStudy.
// .L Plot_CDet_PulseSideband.C+
// Plot_CDet_PulseSideband("CDet_run6077_pulse_pixel_timing.conf");
// getSigmaBars("CDet_run6077_pulse_pixel_timing.conf", "background");
// No config changes required. Default sidebands: xCut < |rx| <= 2*xCut.
// Optional loader arguments: events, input directory, outer |rx|, maximum |ry|.
// Pass INFINITY for outer |rx| to accept all pulses outside the central window.
// The replay spatial flag and dt ellipse are NOT required: keep ry explicitly,
// replace the central rx window with the sidebands, and preserve the other cuts.
// Fits use pulse_study.fit_min/max, not fit_signal_min/max. PDFs go to
// PusleStudy/sidebands/<tag>/; canvases remain interactive. Gaussian sigma here
// describes the BACKGROUND shape, not signal timing resolution.

#include <TCanvas.h>
#include <TChain.h>
#include <TEnv.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TList.h>
#include <TPaveText.h>
#include <TROOT.h>
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

// Distinct names let the normal study and sideband study coexist in ROOT.
const int sidebandPixelsPerBar = 16;
const int sidebandNumBars = 168;
int sidebandRunNumber = 0;
TString sidebandLoadingECalCuts, sidebandLoadingPulseCuts, sidebandLoadingSpatialCuts;
std::vector<double> vSidebandECalE, vSidebandECalAdcTime;
std::vector<std::vector<double>> vSidebandPulseID, vSidebandPulseLe,
    vSidebandPulseToT, vSidebandPulseRx, vSidebandPulseDt;

void ClearPulseSidebandVectors()
{
  sidebandRunNumber = 0;
  sidebandLoadingECalCuts = sidebandLoadingPulseCuts = sidebandLoadingSpatialCuts = "";
  vSidebandECalE.clear();
  vSidebandECalAdcTime.clear();
  vSidebandPulseID.clear();
  vSidebandPulseLe.clear();
  vSidebandPulseToT.clear();
  vSidebandPulseRx.clear();
  vSidebandPulseDt.clear();
}

void Plot_CDet_PulseSideband(
    const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
    Long64_t eventsOverride = -2, const char* inputDirectory = nullptr,
    double RxOuter = NAN, double RyMax = 0.36)
{
  ClearPulseSidebandVectors();
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet sidebands] Cannot read config.version: 1 configuration.\n";
    return;
  }
  const int RunNumber = config.GetValue("analysis.run_number", 0);
  const int minSeg = config.GetValue("analysis.min_segment", -1);
  const int maxSeg = config.GetValue("analysis.max_segment", -1);
  const Long64_t nevents = eventsOverride == -2 ? config.GetValue("analysis.events", -1) : eventsOverride;
  const double ECalMinE = config.GetValue("analysis.ecal_energy_min", NAN);
  const double ECalMaxE = config.GetValue("analysis.ecal_energy_max", NAN);
  const double ECalMinT = config.GetValue("analysis.ecal_time_min", NAN);
  const double ECalMaxT = config.GetValue("analysis.ecal_time_max", NAN);
  const double LeMin = config.GetValue("analysis.le_min", 0.02);
  const double LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0);
  const double TotMax = config.GetValue("analysis.tot_max", 30.0);
  const double RxInner = config.GetValue("analysis.x_difference_max", 0.08);
  if (std::isnan(RxOuter)) RxOuter = 2.0 * RxInner;
  // RyMax defaults to the current replay DB selection.y_residual_max (metres).
  if (RunNumber <= 0 || minSeg < -1 || maxSeg < -1 || nevents < -1 ||
      !std::isfinite(ECalMinE) || !std::isfinite(ECalMaxE) || ECalMinE >= ECalMaxE ||
      !std::isfinite(ECalMinT) || !std::isfinite(ECalMaxT) || ECalMinT >= ECalMaxT ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      !std::isfinite(TotMin) || !std::isfinite(TotMax) || TotMin < 0.0 || TotMin >= TotMax ||
      !std::isfinite(RxInner) || RxInner <= 0.0 || !(RxOuter > RxInner) ||
      !std::isfinite(RyMax) || RyMax <= 0.0) {
    std::cerr << "[CDet sidebands] Invalid run, event/segment, or cut settings.\n";
    return;
  }

  const char* directory = inputDirectory && inputDirectory[0] ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0]) {
    std::cerr << "[CDet sidebands] Set OUT_DIR or pass inputDirectory.\n";
    return;
  }
  TChain T("T");
  if (CDetRunDataset::AddToChain(&T, RunNumber, directory, minSeg, maxSeg) <= 0) return;
  const Long64_t entries = T.GetEntries();
  if (entries <= 0) { std::cerr << "[CDet sidebands] No input events.\n"; return; }
  const Long64_t limit = nevents >= 0 ? std::min(nevents, entries) : entries;

  // Readers before the single event loop; only fields needed for this study.
  TTreeReader reader(&T);
  TTreeReaderValue<double> ECalE(reader, "earm.ecal.e");
  TTreeReaderValue<double> ECalAdcTime(reader, "earm.ecal.adctime");
  TTreeReaderArray<double> PulseID(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<double> PulseLe(reader, "earm.cdet.pulse.tdc_le");
  TTreeReaderArray<double> PulseToT(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<double> PulseRx(reader, "earm.cdet.pulse.ecal_x_residual");
  TTreeReaderArray<double> PulseRy(reader, "earm.cdet.pulse.ecal_y_residual");
  TTreeReaderArray<double> PulseDt(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<double> PulseECalEligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<double> PulseGoodCalib(reader, "earm.cdet.pulse.calib_valid");

  sidebandRunNumber = RunNumber;
  sidebandLoadingECalCuts = TString::Format("Loading: %.3g < E_{ECal} < %.3g GeV; %.3g < t_{ECal} < %.3g ns", ECalMinE, ECalMaxE, ECalMinT, ECalMaxT);
  sidebandLoadingPulseCuts = TString::Format("Loading: %.3g < LE < %.3g ns; %.3g < ToT < %.3g ns; eligible + calib flags", LeMin, LeMax, TotMin, TotMax);
  sidebandLoadingSpatialCuts = std::isfinite(RxOuter)
      ? TString::Format("Sidebands: %.3g < |rx| #leq %.3g m; |ry| #leq %.3g m; no dt ellipse", RxInner, RxOuter, RyMax)
      : TString::Format("Sidebands: |rx| > %.3g m (no outer bound); |ry| #leq %.3g m; no dt ellipse", RxInner, RyMax);
  std::cout << "[CDet sidebands] Run " << RunNumber << ": reading " << limit << " events\n"
            << sidebandLoadingECalCuts << '\n' << sidebandLoadingPulseCuts << '\n' << sidebandLoadingSpatialCuts << '\n';

  Long64_t EventCounter = 0, savedPulses = 0;
  while (EventCounter < limit && reader.Next()) {
    ++EventCounter;
    if (EventCounter % 1000 == 0) std::cout << "[CDet sidebands] " << EventCounter << '/' << limit << " events\n";
    const bool good_ecal_event = std::isfinite(*ECalE) && *ECalE > ECalMinE && *ECalE < ECalMaxE &&
        std::isfinite(*ECalAdcTime) && *ECalAdcTime > ECalMinT && *ECalAdcTime < ECalMaxT;
    if (!good_ecal_event) continue;
    const size_t nhits = PulseID.GetSize();
    if (PulseLe.GetSize() != nhits || PulseToT.GetSize() != nhits || PulseRx.GetSize() != nhits ||
        PulseRy.GetSize() != nhits || PulseDt.GetSize() != nhits ||
        PulseECalEligible.GetSize() != nhits || PulseGoodCalib.GetSize() != nhits) {
      std::cerr << "[CDet sidebands] Inconsistent pulse-array lengths; clearing incomplete sample.\n";
      ClearPulseSidebandVectors();
      return;
    }
    std::vector<double> thisEventID, thisEventLe, thisEventToT, thisEventRx, thisEventDt;
    for (size_t hit = 0; hit < nhits; ++hit) {
      const double id = PulseID[hit], leNs = 0.01 * PulseLe[hit];
      const double rx = PulseRx[hit], ry = PulseRy[hit];
      const bool good_id = std::isfinite(id) && id >= 0 && id < sidebandNumBars * sidebandPixelsPerBar && std::floor(id) == id;
      const bool good_ecal_eligible = PulseECalEligible[hit] > 0.5;
      const bool good_calibration = PulseGoodCalib[hit] > 0.5;
      const bool good_le = std::isfinite(leNs) && leNs > LeMin && leNs < LeMax;
      const bool good_tot = std::isfinite(PulseToT[hit]) && PulseToT[hit] > TotMin && PulseToT[hit] < TotMax;
      const bool good_sideband = std::isfinite(rx) && std::abs(rx) > RxInner && std::abs(rx) <= RxOuter;
      const bool good_y = std::isfinite(ry) && std::abs(ry) <= RyMax;
      if (!good_id || !good_ecal_eligible || !good_calibration || !good_le || !good_tot ||
          !good_sideband || !good_y || !std::isfinite(PulseDt[hit])) continue;
      thisEventID.push_back(id);
      thisEventLe.push_back(leNs);
      thisEventToT.push_back(PulseToT[hit]);
      thisEventRx.push_back(rx);
      thisEventDt.push_back(PulseDt[hit]);
      ++savedPulses;
    }
    if (thisEventID.empty()) continue;
    vSidebandECalE.push_back(*ECalE);
    vSidebandECalAdcTime.push_back(*ECalAdcTime);
    vSidebandPulseID.push_back(thisEventID);
    vSidebandPulseLe.push_back(thisEventLe);
    vSidebandPulseToT.push_back(thisEventToT);
    vSidebandPulseRx.push_back(thisEventRx);
    vSidebandPulseDt.push_back(thisEventDt);
  }
  if (EventCounter != limit) {
    std::cerr << "[CDet sidebands] Reader stopped at " << EventCounter << '/' << limit
              << "; clearing incomplete sample. Status " << reader.GetEntryStatus() << '\n';
    ClearPulseSidebandVectors();
    return;
  }
  std::cout << "[CDet sidebands] Saved " << savedPulses << " pulses in " << vSidebandPulseID.size() << " events.\n";
}

void Plot_CDet_PulseSideband(int run, Long64_t events = -2, const char* inputDirectory = nullptr,
                            double RxOuter = NAN, double RyMax = 0.36)
{
  Plot_CDet_PulseSideband(TString::Format("CDet_run%d_pulse_pixel_timing.conf", run), events, inputDirectory, RxOuter, RyMax);
}

void getSigmaBars(const char* configFile = "CDet_run6077_pulse_pixel_timing.conf",
                  const char* outputTag = nullptr, int firstBar = -1, int lastBar = -1,
                  double BinWidth = NAN, double DtMin = NAN, double DtMax = NAN)
{
  TEnv config;
  if (!configFile || !configFile[0] || config.ReadFile(configFile, kEnvLocal) != 0 ||
      config.GetValue("config.version", 0) != 1) {
    std::cerr << "[CDet sideband fits] Cannot read configuration.\n"; return;
  }
  if (vSidebandPulseID.empty()) {
    std::cerr << "[CDet sideband fits] Run Plot_CDet_PulseSideband first; no saved pulses.\n"; return;
  }
  if (config.GetValue("analysis.run_number", 0) != sidebandRunNumber) {
    std::cerr << "[CDet sideband fits] Config run does not match the loaded run.\n"; return;
  }
  if (std::isnan(BinWidth)) BinWidth = config.GetValue("pulse_study.dt_bin_width_ns", 1.0);
  if (std::isnan(DtMin)) DtMin = config.GetValue("display.ecal_cdet_dt_min", -150.0);
  if (std::isnan(DtMax)) DtMax = config.GetValue("display.ecal_cdet_dt_max", 150.0);
  if (firstBar == -1) firstBar = config.GetValue("pulse_study.first_bar", 0);
  if (lastBar == -1) lastBar = config.GetValue("pulse_study.last_bar", sidebandNumBars - 1);
  const int controlBar = config.GetValue("pulse_pixel_timing.control_bar", 30);
  const double FitMin = config.GetValue("pulse_study.fit_min", -50.0);
  const double FitMax = config.GetValue("pulse_study.fit_max", 0.0);
  const double MinEntries = config.GetValue("pulse_study.min_fit_entries", 100.0);
  const double LeMin = config.GetValue("analysis.le_min", 0.02), LeMax = config.GetValue("analysis.le_max", 60.0);
  const double TotMin = config.GetValue("analysis.tot_min", 4.0), TotMax = config.GetValue("analysis.tot_max", 30.0);
  const bool applyToT = config.GetValue("pulse_study.apply_tot_cut", 0) != 0;
  const bool applyEnergy = config.GetValue("pulse_study.apply_ecal_energy_cut", 0) != 0;
  const bool applyTime = config.GetValue("pulse_study.apply_ecal_time_cut", 0) != 0;
  const double Emin = config.GetValue("pulse_study.ecal_energy_min", NAN), Emax = config.GetValue("pulse_study.ecal_energy_max", NAN);
  const double Tmin = config.GetValue("pulse_study.ecal_time_min", NAN), Tmax = config.GetValue("pulse_study.ecal_time_max", NAN);
  if (firstBar < 0 || lastBar >= sidebandNumBars || firstBar > lastBar || controlBar < 0 || controlBar >= sidebandNumBars ||
      !std::isfinite(BinWidth) || BinWidth <= 0 || !std::isfinite(DtMin) || !std::isfinite(DtMax) || DtMin >= DtMax ||
      !std::isfinite((DtMax-DtMin)/BinWidth) || std::ceil((DtMax-DtMin)/BinWidth) > std::numeric_limits<int>::max()-2 ||
      !std::isfinite(FitMin) || !std::isfinite(FitMax) || FitMin >= FitMax || FitMin < DtMin || FitMax > DtMax ||
      !std::isfinite(MinEntries) || MinEntries <= 0 ||
      !std::isfinite(LeMin) || !std::isfinite(LeMax) || LeMin >= LeMax ||
      (applyToT && (!std::isfinite(TotMin) || !std::isfinite(TotMax) || TotMin < 0 || TotMin >= TotMax)) ||
      (applyEnergy && (!std::isfinite(Emin) || !std::isfinite(Emax) || Emin >= Emax)) ||
      (applyTime && (!std::isfinite(Tmin) || !std::isfinite(Tmax) || Tmin >= Tmax))) {
    std::cerr << "[CDet sideband fits] Invalid plot, fit, or cut settings.\n"; return;
  }
  const int nBins = int(std::ceil((DtMax-DtMin)/BinWidth));
  const double sigmaMin = 0.5*(DtMax-DtMin)/nBins, sigmaMax = 2.0*(FitMax-FitMin);
  TString tag = outputTag ? outputTag : config.GetValue("pulse_study.output_tag", "");
  if (tag.IsNull()) tag = TString::Format("run%d",sidebandRunNumber);
  const TString allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";
  for (int i = 0; i < tag.Length(); ++i) if (allowed.First(tag[i]) == kNPOS) {
    std::cerr << "[CDet sideband fits] Invalid output tag.\n"; return;
  }
  const TString out = TString::Format("%s/PusleStudy/sidebands/%s",gSystem->WorkingDirectory(),tag.Data());
  if (gSystem->AccessPathName(out)) gSystem->mkdir(out,true);
  FileStat_t info;
  if (gSystem->GetPathInfo(out,info) || !R_ISDIR(info.fMode) || gSystem->AccessPathName(out,kWritePermission)) {
    std::cerr << "[CDet sideband fits] Cannot write " << out << '\n'; return;
  }

  const char* canvasNames[] = {"cSidebandControlDt","cSidebandDetectorDt","cSidebandSigmaVsBar"};
  for (const char* name : canvasNames) delete gROOT->GetListOfCanvases()->FindObject(name);
  static std::vector<TH1D*> histograms(sidebandNumBars+1,nullptr);
  static std::vector<std::vector<TF1*>> fits(sidebandNumBars+1,std::vector<TF1*>(3,nullptr));
  static TGraphErrors *layer1 = nullptr, *layer2 = nullptr, *unavailable = nullptr;
  delete layer1; delete layer2; delete unavailable;
  layer1 = new TGraphErrors(); layer2 = new TGraphErrors(); unavailable = new TGraphErrors();
  layer1->SetName("gSidebandSigmaL1"); layer2->SetName("gSidebandSigmaL2"); unavailable->SetName("gSidebandSigmaUnavailable");
  const char* formulas[] = {"pol2","pol4","gaus"};
  const char* modelNames[] = {"Quadratic","Fourth-order polynomial","Gaussian"};
  const int colors[] = {kRed+1,kBlue+1,kGreen+2};
  for (int bar = 0; bar <= sidebandNumBars; ++bar) {
    delete histograms[bar];
    const TString label = bar == sidebandNumBars ? TString("Full detector") : TString::Format("Bar %d",bar);
    histograms[bar] = new TH1D(TString::Format("hSidebandDt_%d",bar),
        TString::Format("Run %d: %s sideband dt;dt = pulse.ecal_residual (ns);Pulses",sidebandRunNumber,label.Data()),nBins,DtMin,DtMax);
    histograms[bar]->SetDirectory(nullptr); histograms[bar]->SetStats(false);
    for (int model = 0; model < 3; ++model) {
      delete fits[bar][model];
      fits[bar][model] = new TF1(TString::Format("fSidebandDt_%d_%d",bar,model),formulas[model],FitMin,FitMax);
      fits[bar][model]->AddToGlobalList(false);
      fits[bar][model]->SetLineColor(colors[model]);
      fits[bar][model]->SetLineStyle(model+1);
      fits[bar][model]->SetNpx(1000);
    }
  }

  Long64_t filledPulses = 0;
  for (size_t event = 0; event < vSidebandPulseID.size(); ++event) {
    const double e = vSidebandECalE[event], t = vSidebandECalAdcTime[event];
    const bool good_ecal = std::isfinite(e) && std::isfinite(t) &&
        (!applyEnergy || (e > Emin && e < Emax)) && (!applyTime || (t > Tmin && t < Tmax));
    if (!good_ecal) continue;
    for (size_t hit = 0; hit < vSidebandPulseID[event].size(); ++hit) {
      const double le = vSidebandPulseLe[event][hit], tot = vSidebandPulseToT[event][hit];
      const bool good_le = std::isfinite(le) && le >= LeMin && le <= LeMax;
      const bool good_tot = std::isfinite(tot) && tot > 0 && (!applyToT || (tot >= TotMin && tot <= TotMax));
      if (!good_le || !good_tot) continue;
      const int bar = int(vSidebandPulseID[event][hit])/sidebandPixelsPerBar;
      const double dt = vSidebandPulseDt[event][hit];
      histograms[bar]->Fill(dt); histograms[sidebandNumBars]->Fill(dt); ++filledPulses;
    }
  }

  // Independent alternatives fitted to the SAME bins; no signal exclusion in dt.
  std::vector<std::vector<bool>> valid(sidebandNumBars+1,std::vector<bool>(3,false));
  std::vector<std::vector<TString>> reasons(sidebandNumBars+1,std::vector<TString>(3));
  for (int bar = 0; bar <= sidebandNumBars; ++bar) {
    TH1D* h = histograms[bar];
    double entries = 0, height = 0, sumX = 0, sumX2 = 0;
    int occupied = 0;
    for (int bin = 1; bin <= nBins; ++bin) {
      const double x = h->GetBinCenter(bin), n = h->GetBinContent(bin);
      if (x < FitMin || x > FitMax) continue;
      entries += n; sumX += n*x; sumX2 += n*x*x;
      height = std::max(height,n); if (n > 0) ++occupied;
    }
    for (int model = 0; model < 3; ++model) {
      TF1* f = fits[bar][model];
      if (entries < MinEntries || occupied <= f->GetNpar()) {
        reasons[bar][model] = "Insufficient fit-window statistics"; continue;
      }
      if (model == 2) {
        const double mu = sumX/entries;
        const double rms = std::sqrt(std::max(0.0,sumX2/entries-mu*mu));
        f->SetParameters(height,mu,std::max(sigmaMin*1.01,std::min(sigmaMax*.99,rms)));
        f->SetParLimits(0,0,10*height); f->SetParLimits(1,FitMin,FitMax); f->SetParLimits(2,sigmaMin,sigmaMax);
      }
      // B keeps Gaussian seeds/limits; polynomial coefficients use ROOT's linear fitter.
      TFitResultPtr result = h->Fit(f,model == 2 ? "QRS0NB" : "QRS0N");
      bool good_fit = int(result)==0 && result.Get() && result->IsValid() && result->CovMatrixStatus()>=2 && f->GetNDF()>0;
      for (int p = 0; p < f->GetNpar(); ++p) good_fit = good_fit && std::isfinite(f->GetParameter(p)) && std::isfinite(f->GetParError(p));
      if (model == 2) {
        const double sigma = f->GetParameter(2), mu = f->GetParameter(1), tol = .001*(FitMax-FitMin);
        good_fit = good_fit && f->GetParameter(0)>0 && f->GetParError(2)>0 &&
            sigma>sigmaMin*1.001 && sigma<sigmaMax*.999 && mu>FitMin+tol && mu<FitMax-tol;
      }
      valid[bar][model] = good_fit;
      if (!good_fit) reasons[bar][model] = "Unusable fit (status, covariance, or limits)";
      std::cout << "[CDet sideband fits] " << (bar==sidebandNumBars ? TString("Detector") : TString::Format("Bar %d",bar))
                << " " << modelNames[model] << ": " << (good_fit ? TString::Format("chi2/ndf = %.2f/%d",f->GetChisquare(),f->GetNDF()) : reasons[bar][model]) << '\n';
    }
  }

  double graphMax = 1;
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    if (!valid[bar][2]) { unavailable->SetPoint(unavailable->GetN(),bar,0); continue; }
    const double sigma = fits[bar][2]->GetParameter(2), error = fits[bar][2]->GetParError(2);
    TGraphErrors* graph = bar < 84 ? layer1 : layer2;
    const int point = graph->GetN(); graph->SetPoint(point,bar,sigma); graph->SetPointError(point,0,error);
    graphMax = std::max(graphMax,sigma+error);
  }
  layer1->SetMarkerColor(kRed); layer1->SetLineColor(kRed);
  layer2->SetMarkerColor(kBlue); layer2->SetLineColor(kBlue);
  unavailable->SetMarkerColor(kBlack); unavailable->SetLineColor(kBlack);
  for (TGraphErrors* g : {layer1,layer2,unavailable}) { g->SetMarkerStyle(20); g->SetMarkerSize(.85); }
  for (int plot = 0; plot < 3; ++plot) {
    TCanvas* canvas = new TCanvas(canvasNames[plot],canvasNames[plot],1200,742);
    canvas->SetCanvasSize(1200,742);
    canvas->SetTopMargin(.43); canvas->SetBottomMargin(.12); canvas->SetLeftMargin(.10); canvas->SetRightMargin(.04);
    TLegend* legend = new TLegend(.10,.585,.96,.715);
    legend->SetBit(kCanDelete); legend->SetFillColor(kWhite); legend->SetBorderSize(0); legend->SetTextSize(.022);
    if (plot < 2) {
      const int bar = plot==0 ? controlBar : sidebandNumBars;
      TH1D* h = histograms[bar];
      h->SetMaximum(std::max(1.0,1.15*h->GetMaximum())); h->Draw("HIST");
      for (int model = 0; model < 3; ++model) {
        TF1* f = fits[bar][model];
        TString label = modelNames[model];
        if (valid[bar][model]) {
          f->Draw("SAME");
          label += TString::Format(": #chi^{2}/ndf = %.1f/%d",f->GetChisquare(),f->GetNDF());
          if (model == 2) label += TString::Format("; #sigma_{background} = %.2f #pm %.2f ns",f->GetParameter(2),f->GetParError(2));
        } else label += ": " + reasons[bar][model];
        legend->AddEntry(f,label,"l");
      }
    } else {
      canvas->DrawFrame(firstBar-.5,-.08*graphMax,lastBar+.5,1.2*graphMax,
          TString::Format("Run %d: sideband Gaussian width (not signal resolution);Bar number (global, zero-based);#sigma_{background} (ns)",sidebandRunNumber));
      layer1->Draw("P SAME"); layer2->Draw("P SAME"); unavailable->Draw("P SAME");
      legend->SetNColumns(2);
      legend->AddEntry(layer1,"Layer 1 (0-83)","pe"); legend->AddEntry(layer2,"Layer 2 (84-167)","pe");
      legend->AddEntry(unavailable,"Low statistics / unusable Gaussian fit: 0","p");
    }
    legend->Draw();
    TPaveText* cuts = new TPaveText(.10,.73,.96,.94,"NDC");
    cuts->SetBit(kCanDelete); cuts->SetFillColor(kWhite); cuts->SetBorderSize(0); cuts->SetTextAlign(12); cuts->SetTextSize(.020);
    cuts->AddText(TString::Format("Fit [%.3g, %.3g] ns; bins %.3g ns; both rx sidebands combined",FitMin,FitMax,(DtMax-DtMin)/nBins));
    cuts->AddText(sidebandLoadingECalCuts); cuts->AddText(sidebandLoadingPulseCuts); cuts->AddText(sidebandLoadingSpatialCuts);
    const TString energy = applyEnergy ? TString::Format("%.3g < E < %.3g GeV",Emin,Emax) : TString("energy off");
    const TString time = applyTime ? TString::Format("%.3g < t < %.3g ns",Tmin,Tmax) : TString("ADC time off");
    cuts->AddText("Secondary ECal: " + energy + "; " + time);
    cuts->AddText(TString::Format("Secondary: %.3g #leq LE #leq %.3g ns; %s",LeMin,LeMax,
        applyToT ? TString::Format("%.3g #leq ToT #leq %.3g ns; ToT > 0",TotMin,TotMax).Data() : "finite ToT > 0 (window off)"));
    cuts->Draw(); canvas->Modified(); canvas->Update();
    const TString name = plot==0 ? TString::Format("bar%03d_background_dt.pdf",controlBar)
        : TString(plot==1 ? "detector_background_dt.pdf" : "background_sigma_vs_bar.pdf");
    canvas->SaveAs(out+"/"+name);
  }
  std::cout << "[CDet sideband fits] " << filledPulses << " selected pulses; PDFs in " << out
            << ". Polynomial fits have no sigma; graph shows the Gaussian background width only.\n";
}
