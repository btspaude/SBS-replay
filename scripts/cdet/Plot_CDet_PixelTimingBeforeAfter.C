// Six bar canvases and two full-detector timing-correlation canvases.
// ROOT (from scripts/cdet, in a fresh session):
//   .L Plot_CDet_PixelTimingBeforeAfter.C+
//   Plot_CDet_PixelTimingBeforeAfter("CDet_run5710_projection.conf", 469);
// Replay-pair DT versus DX and the all-pulse Bar-30 timing reproduction:
//   Plot_CDet_PairDTvsDXAndBarTiming("CDet_run6077_projection.conf");
// Any logical pixel 0..2687 selects its containing bar (pixel/16).
// Bar canvases show its middle eight pixels (bar*16+4 through bar*16+11).
// Saved LE/ToT polygons select one common sample in their pixel-aligned
// drawing coordinates, before either timing state is filled.
// Reuse the master's configuration/calibration readers, constants and bad-pixel
// list, but NOT its main routine or calibration-writing routines.
// Plotting only: fill matched before/after histograms in one event pass;
// no event vectors or calibration products are exported for downstream analysis.
#include "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C"
#include <array>
#include <memory>
#include <TPad.h>

namespace CDetPixelTimingBeforeAfter {
// Plot-only form of the automatic path in
// extractHierarchicalCDetPixelTimingOffsetsDiagnostic(): same group seeds,
// constrained/broad fits, quality requirements and detector-reference median.
// The input spectra already include any requested polygon selection.
// No offset updates or calibration-file writes occur here.
double FitDTSpectra(const std::vector<std::unique_ptr<TH1D>>& spectra,
                   double fitMin, double fitMax,
                   std::vector<std::unique_ptr<TF1>>& fits,
                   std::vector<TString>& notes, TString& referenceDescription)
{
  const double width = spectra[0]->GetBinWidth(1);
  const int minPixelEntries = 35, minGroupEntries = 100;
  const double minSigma = 0.5, maxSigma = 8.0;
  fits.resize(NumCDetPaddles); notes.resize(NumCDetPaddles);
  std::vector<double> individualCentroids, groupCentroids;
  auto makeFit = [&](TH1D& h, double lo, double hi, double amplitude,
                     double mean, double sigma, double background) {
    std::unique_ptr<TF1> fit(new TF1(TString(h.GetName())+"_fit", "gaus(0)+pol1(3)", lo, hi));
    fit->SetParameters(amplitude, mean, sigma, background, 0.0);
    fit->SetParLimits(0, 0.0, 1.5*std::max(1.0, h.GetMaximum()));
    fit->SetParLimits(1, lo, hi); fit->SetParLimits(2, minSigma, maxSigma);
    return fit;
  };
  auto peakBin = [](TH1D& h, int lo, int hi) {
    int peak = lo;
    for (int bin = lo+1; bin <= hi; ++bin)
      if (h.GetBinContent(bin) > h.GetBinContent(peak)) peak = bin;
    return peak;
  };
  for (int first = 0; first < NumCDetPaddles; first += 8) {
    std::unique_ptr<TH1D> group(static_cast<TH1D*>(spectra[first]->Clone(TString(spectra[first]->GetName())+"_group")));
    group->SetDirectory(nullptr); group->Reset();
    // ROOT assigns x == histogram maximum to overflow. A fit may end at
    // that edge, but its statistics, seed and background must use visible bins.
    const int lo = std::max(1, group->FindFixBin(fitMin));
    const int hi = std::min(group->GetNbinsX(), group->FindFixBin(fitMax));
    for (int channel = first; channel < first+8; ++channel)
      if (!IsUnusedPixel(channel) && spectra[channel]->Integral(lo, hi) > 0)
        group->Add(spectra[channel].get()); // Raw counts, not normalized pixels.

    bool groupValid = false;
    double groupMean = NAN, groupSigma = NAN;
    if (group->Integral(lo, hi) >= minGroupEntries) {
      const int peak = peakBin(*group, lo, hi);
      const double background = 0.5*(group->GetBinContent(lo)+group->GetBinContent(hi));
      auto fit = makeFit(*group, fitMin, fitMax,
          std::max(0.01, group->GetBinContent(peak)-background),
          group->GetBinCenter(peak), 3.0, background);
      const int status = group->Fit(fit.get(), "RQN0");
      groupMean = fit->GetParameter(1); groupSigma = std::fabs(fit->GetParameter(2));
      const double error = fit->GetParError(1);
      groupValid = status == 0 && std::isfinite(groupMean) &&
          std::isfinite(error) && error > 0 && std::isfinite(groupSigma) &&
          groupSigma >= minSigma && groupSigma <= maxSigma &&
          groupMean > fitMin+width && groupMean < fitMax-width;
      if (groupValid) groupCentroids.push_back(groupMean);
    }

    for (int channel = first; channel < first+8; ++channel) {
      TH1D& h = *spectra[channel];
      if (IsUnusedPixel(channel)) { notes[channel] = "Unused pixel"; continue; }
      const bool enough = h.Integral(lo, hi) >= minPixelEntries;
      notes[channel] = enough ? "Pixel fit rejected" : "Low statistics in fit window";
      if (groupValid)
        notes[channel] = TString::Format("%s; group #mu = %.2f ns", enough ? "Fit rejected" : "Low N", groupMean);
      if (!enough) continue;

      // The broad individual fit is diagnostic when the group is valid; it
      // becomes the pixel fallback only when that group fit is invalid.
      const int peak = peakBin(h, lo, hi);
      auto fit = makeFit(h, fitMin, fitMax, std::max(1.0, h.GetBinContent(peak)),
                        h.GetBinCenter(peak), 3.0, 0.0);
      int status = h.Fit(fit.get(), "RQN0");
      double low = fitMin, high = fitMax;
      if (groupValid) {
        const double halfWindow = std::min(8.0, std::max(4.0, 2.0*groupSigma));
        low = std::max(fitMin, groupMean-halfWindow);
        high = std::min(fitMax, groupMean+halfWindow);
        fit.reset();
        fit = makeFit(h, low, high,
            std::max(1.0, h.GetBinContent(h.FindBin(groupMean))),
            groupMean, std::min(3.0, groupSigma), 0.0);
        status = h.Fit(fit.get(), "RQN0");
      }
      const double mean = fit->GetParameter(1), error = fit->GetParError(1);
      const double sigma = std::fabs(fit->GetParameter(2));
      const double amplitude = fit->GetParameter(0), amplitudeError = fit->GetParError(0);
      const double significance = amplitudeError > 0 ? amplitude/amplitudeError : NAN;
      const int ndf = fit->GetNDF();
      const double chi2Ndf = ndf > 0 ? fit->GetChisquare()/ndf : NAN;
      const bool valid = status == 0 && std::isfinite(mean) &&
          std::isfinite(error) && error > 0 && error <= 2.0 &&
          std::isfinite(sigma) && sigma >= minSigma && sigma <= maxSigma &&
          mean-low > 0.5*width && high-mean > 0.5*width && amplitude > 0 &&
          std::isfinite(significance) && significance >= 1.5 && ndf > 0 &&
          std::isfinite(chi2Ndf) && chi2Ndf <= 15.0;
      if (!valid) continue;
      individualCentroids.push_back(mean);
      notes[channel] = groupValid ? "individual_fit" : "individual_broad_fallback";
      fit->SetLineColor(kRed+1); fit->SetLineWidth(2); fit->SetNpx(500);
      fits[channel] = std::move(fit);
    }
  }
  // A group fallback is never presented as an independently measured mu_i.
  // Only when there are no valid pixel fits do groups define mu_0, as in master.
  const bool useGroups = individualCentroids.empty();
  std::vector<double>& centroids = useGroups ? groupCentroids : individualCentroids;
  referenceDescription = TString::Format("median of %zu valid %s fits", centroids.size(), useGroups ? "group" : "pixel");
  if (centroids.empty()) return NAN;
  std::sort(centroids.begin(), centroids.end());
  const size_t middle = centroids.size()/2;
  return centroids.size()%2 ? centroids[middle] : 0.5*(centroids[middle-1]+centroids[middle]);
}
} // namespace CDetPixelTimingBeforeAfter

void Plot_CDet_PixelTimingBeforeAfter(
    const char *configFile = "CDet_run5710_projection.conf", int pixel = 469,
    const char *outputDirectory = "cdet_timing_before_after",
    const char *inputDirectory = nullptr, Long64_t eventsOverride = -2,
    double leBinWidth = 1, double leMin = 0, double leMax = 60,
    double totBinWidth = 1, double totMin = 0, double totMax = 40,
    bool savePlots = true,
    double dtBinWidth = 1, double dtMin = -40, double dtMax = 10,
    double ecalBinWidth = 1, double ecalMin = 5, double ecalMax = 40,
    double dtBeforeFitMin = -30, double dtBeforeFitMax = 10,
    double dtAfterFitMin = -30, double dtAfterFitMax = 10,
    const char *pixelCutFile = "CDet_run5710_halfbar_aligned_final_archive/"
                              "CDet_pixel_quality_cuts_run5710_halfbar_aligned_final.root")
{
  ApplyCDetPlotStyle();
  const int nLE = leBinWidth > 0 ? int(std::ceil((leMax-leMin)/leBinWidth)) : 0;
  const int nToT = totBinWidth > 0 ? int(std::ceil((totMax-totMin)/totBinWidth)) : 0;
  const int nDT = dtBinWidth > 0 ? int(std::ceil((dtMax-dtMin)/dtBinWidth)) : 0;
  const int nECal = ecalBinWidth > 0 ? int(std::ceil((ecalMax-ecalMin)/ecalBinWidth)) : 0;
  if (pixel < 0 || pixel >= NumCDetPaddles || nLE <= 0 || nToT <= 0 || nDT <= 0 || nECal <= 0 ||
      (savePlots && (!outputDirectory || !outputDirectory[0]))) {
    std::cerr << "[CDet before/after] Invalid pixel, histogram range or output directory.\n";
    return;
  }
  const double fitMin[] = {dtBeforeFitMin, dtAfterFitMin};
  const double fitMax[] = {dtBeforeFitMax, dtAfterFitMax};
  for (int state = 0; state < 2; ++state) {
    if (!std::isfinite(fitMin[state]) || !std::isfinite(fitMax[state]) ||
        fitMin[state] >= fitMax[state] || fitMin[state] < dtMin || fitMax[state] > dtMax) {
      std::cerr << "[CDet before/after] " << (state ? "After" : "Before")
                << " DT fit range [" << fitMin[state] << ", " << fitMax[state]
                << "] ns must be ordered and inside the DT histogram ["
                << dtMin << ", " << dtMax << "] ns (equal endpoints are allowed).\n";
      return;
    }
  }
  TEnv env;
  if (!LoadCDetConfiguration(env, configFile, "CDet before/after")) return;
  const int run = env.GetValue("analysis.run_number", 5710);
  const Long64_t events = eventsOverride == -2
      ? env.GetValue("analysis.events", -1) : eventsOverride;
  if (run <= 0) {
    std::cerr << "[CDet before/after] Select one run so its timing shift is unambiguous.\n";
    return;
  }

  // Load existing constants only. Missing calibration must not silently look
  // like a successfully calibrated 'after' sample. No constants are fitted/written.
  if (!LoadCalibrationConstants("CDet_calibration_dt.dat") ||
      !gPixelToffsetLoaded || !gECalParamsLoaded || !gECalDeltaLoaded ||
      !gTimeWalkParamsLoaded) {
    std::cerr << "[CDet before/after] Need complete PixelOffsets, ECalTiming (including delta), and TimeWalk sections.\n";
    return;
  }
  LoadRunTimingConstants(TString::Format("CDet_run%d.dat", run).Data());
  if (!gGlobalTimingLoaded) {
    std::cerr << "[CDet before/after] Need shift_ns in CDet_run" << run << ".dat.\n";
    return;
  }

  // Use the same reviewed file as the accepted Run-5710 workflow by default.
  // Clone the cuts so they remain valid after the read-only file closes.
  // An explicit empty/null path disables polygons for an uncut comparison.
  std::vector<std::unique_ptr<TCutG>> pixelCuts(NumCDetPaddles);
  int loadedPixelCuts = 0;
  if (pixelCutFile && pixelCutFile[0]) {
    TFile cutInput(pixelCutFile, "READ");
    if (cutInput.IsZombie()) {
      std::cerr << "[CDet before/after] Cannot read requested pixel-cut file: "
                << pixelCutFile << ". No plots produced.\n";
      return;
    }
    for (int channel = 0; channel < NumCDetPaddles; ++channel) {
      pixelCuts[channel].reset(LoadCDetPixelLeTotCut(cutInput, channel));
      if (!pixelCuts[channel]) continue;
      const auto *sourceRun = dynamic_cast<TParameter<int>*>(
          cutInput.Get(CDetPixelCutDirectory(channel)+"/source_run"));
      if (sourceRun && sourceRun->GetVal() != run) {
        std::cerr << "[CDet before/after] Pixel " << channel
                  << " polygon was drawn for run " << sourceRun->GetVal()
                  << ", not requested run " << run
                  << ". Supply matching cuts or explicitly disable them; no plots produced.\n";
        return;
      }
      ++loadedPixelCuts;
    }
    if (loadedPixelCuts == 0) {
      std::cerr << "[CDet before/after] No pixel LE/ToT polygons found in "
                << pixelCutFile << ". No plots produced.\n";
      return;
    }
    std::cout << "[CDet before/after] Loaded " << loadedPixelCuts
              << " pixel LE/ToT polygons from " << pixelCutFile
              << "; applied to BOTH samples at (ToT, pixel-aligned LE),"
                 " before ECal/time-walk/run-shift corrections.\n";
  } else {
    std::cout << "[CDet before/after] Pixel polygons explicitly disabled.\n";
  }

  // The defaults, endpoints and run-file timing-window override match the
  // cross-target master's configuration entry point and event loop.
  const double cutLEMin = env.GetValue("analysis.le_min", 0.02);
  const double cutLEMax = env.GetValue("analysis.le_max", 100.0);
  const double cutToTMin = env.GetValue("analysis.tot_min", 0.02);
  const double cutToTMax = env.GetValue("analysis.tot_max", 150.0);
  const double ecalTimeMin = gRunECalTimeWindowLoaded ? gRunECalTimeMin : env.GetValue("analysis.ecal_time_min", 10.0);
  const double ecalTimeMax = gRunECalTimeWindowLoaded ? gRunECalTimeMax : env.GetValue("analysis.ecal_time_max", 35.0);
  const double ecalEnergyMin = env.GetValue("analysis.ecal_energy_min", -1.0e9);
  const double ecalEnergyMax = env.GetValue("analysis.ecal_energy_max", 1.0e9);
  const int minHits[2] = {env.GetValue("analysis.layer1_hits_min", 1), env.GetValue("analysis.layer2_hits_min", 0)};
  const int maxHits[2] = {env.GetValue("analysis.layer1_hits_max", 100), env.GetValue("analysis.layer2_hits_max", 100)};
  const double xDiffMax = env.GetValue("analysis.x_difference_max", 0.05);
  const double xOffset = env.GetValue("analysis.x_offset", 0.0);
  const double yOffset = env.GetValue("analysis.y_offset", 0.1);
  const int layerChoice = env.GetValue("analysis.layer_choice", 3);
  const bool suppressBad = env.GetValue("analysis.suppress_bad", false);
  const bool useReference = env.GetValue("analysis.use_reference_timing", false);
  if (cutLEMin >= cutLEMax || cutToTMin >= cutToTMax ||
      ecalTimeMin >= ecalTimeMax || ecalEnergyMin >= ecalEnergyMax ||
      minHits[0] < 0 || minHits[0] > maxHits[0] ||
      minHits[1] < 0 || minHits[1] > maxHits[1] || xDiffMax <= 0 ||
      layerChoice < 1 || layerChoice > 3) {
    std::cerr << "[CDet before/after] Invalid analysis cuts.\n";
    return;
  }
  const char *directory = inputDirectory && inputDirectory[0]
      ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0]) {
    std::cerr << "[CDet before/after] Set OUT_DIR or pass inputDirectory.\n";
    return;
  }
  TChain chain("T");
  AddRunFilesToChain(&chain, directory, run,
      env.GetValue("analysis.min_segment", -1), env.GetValue("analysis.max_segment", -1));
  const Long64_t total = chain.GetEntries();
  if (total <= 0) return;
  const Long64_t limit = events > 0 ? std::min(events, total) : total;

  // Same good-hit collection as the calibration master, not replay pulse/pair selections.
  TTreeReader reader(&chain);
  TTreeReaderArray<double> id(reader, "earm.cdet.hit.pmtnum");
  TTreeReaderArray<double> le(reader, "earm.cdet.hit.tdc_le");
  TTreeReaderArray<double> tot(reader, "earm.cdet.hit.tdc_tot");
  TTreeReaderArray<double> x(reader, "earm.cdet.hit.xhit");
  TTreeReaderArray<double> y(reader, "earm.cdet.hit.yhit");
  TTreeReaderArray<double> z(reader, "earm.cdet.hit.zhit");
  TTreeReaderArray<double> mult(reader, "earm.cdet.tdc_mult");
  TTreeReaderArray<double> rawID(reader, "earm.cdet.hits.TDCelemID");
  TTreeReaderArray<double> rawLE(reader, "earm.cdet.hits.t");
  TTreeReaderArray<double> rawToT(reader, "earm.cdet.hits.t_tot");
  TTreeReaderValue<double> ecalX(reader, "earm.ecal.x"), ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<double> ecalE(reader, "earm.ecal.e"), ecalT(reader, "earm.ecal.adctime");

  // Define all histograms together, before the event/hit filling loops.
  const int bar = pixel / 16, base = 16 * bar;
  static int invocation = 0;
  const TString tag = TString::Format("run%d_bar%03d_%d", run, bar, ++invocation);
  std::array<std::array<std::unique_ptr<TH1D>, 16>, 2> hLE;
  // All pixels contribute to the run/sample's reference, regardless of which
  // bar is displayed. These are temporary histograms, not event vectors.
  std::array<std::vector<std::unique_ptr<TH1D>>, 2> hDT;
  std::array<std::array<std::unique_ptr<TH2D>, 16>, 2> hLEToT;
  std::array<std::unique_ptr<TH2D>, 2> hCDetECal;
  const char *stage[] = {"before", "after"};
  for (int state = 0; state < 2; ++state) {
    hCDetECal[state].reset(new TH2D(TString::Format("hCDetECal_%s_%s_detector", tag.Data(), stage[state]), ";ECal ADC time (ns);CDet LE time (ns);Hits", nECal, ecalMin, ecalMax, nLE, leMin, leMax));
    hCDetECal[state]->SetDirectory(nullptr); hCDetECal[state]->SetStats(false);
    hCDetECal[state]->SetStatOverflows(TH1::kConsider);
    hDT[state].resize(NumCDetPaddles);
    for (int channel = 0; channel < NumCDetPaddles; ++channel) {
      hDT[state][channel].reset(new TH1D(TString::Format("hDT_%s_%s_p%d", tag.Data(), stage[state], channel), TString::Format("Pixel %d;t_{ECal} - t_{CDet,LE} (ns);Counts / %.3g ns", channel, (dtMax-dtMin)/nDT), nDT, dtMin, dtMax));
      TH1D *hist = hDT[state][channel].get();
      hist->SetDirectory(nullptr); hist->SetStats(false); hist->Sumw2();
      hist->SetStatOverflows(TH1::kConsider);
      hist->SetLineWidth(2); hist->SetLineColor(state ? kBlue+1 : kGray+2);
    }
    for (int p = 0; p < 16; ++p) {
      hLE[state][p].reset(new TH1D(TString::Format("hLE_%s_%s_p%d", tag.Data(), stage[state], base+p), TString::Format("Pixel %d;LE (ns);Counts / %.3g ns", base+p, (leMax-leMin)/nLE), nLE, leMin, leMax));
      hLEToT[state][p].reset(new TH2D(TString::Format("hLEToT_%s_%s_p%d", tag.Data(), stage[state], base+p), TString::Format("Pixel %d;ToT (ns);LE (ns)", base+p), nToT, totMin, totMax, nLE, leMin, leMax));
      for (TH1 *hist : {static_cast<TH1*>(hLE[state][p].get()), static_cast<TH1*>(hLEToT[state][p].get())}) {
        hist->SetDirectory(nullptr); hist->SetStats(false);
        hist->SetStatOverflows(TH1::kConsider);
      }
      hLE[state][p]->SetLineWidth(2);
      hLE[state][p]->SetLineColor(state ? kBlue+1 : kGray+2);
    }
  }

  Long64_t processed = 0, goodEvents = 0, selectedBarHits = 0;
  Long64_t polygonTestedHits = 0, polygonRejectedHits = 0, polygonRejectedBarHits = 0;
  while (processed < limit && reader.Next()) {
    ++processed;
    if (processed % 1000 == 0)
      std::cout << "[CDet before/after] " << processed << "/" << limit
                << " (file entries " << total << "); good events "
                << goodEvents << "; bar hits " << selectedBarHits << '\n';
    const size_t n = id.GetSize();
    if (le.GetSize() != n || tot.GetSize() != n || x.GetSize() != n ||
        y.GetSize() != n || z.GetSize() != n || mult.GetSize() < n ||
        rawLE.GetSize() != rawID.GetSize() || rawToT.GetSize() != rawID.GetSize()) {
      std::cerr << "[CDet before/after] Mismatched arrays at entry " << reader.GetCurrentEntry() << '\n';
      return;
    }
    // Start with the master's event and hit selections. Saved polygons are
    // an additional common hit mask below; no elastic/pair/ellipse cut is added.
    if (!(*ecalX > -1.5 && *ecalX < 1.5 && *ecalX != 0 &&
          *ecalY > -1.2 && *ecalY < 1.2 && *ecalY != 0 &&
          *ecalT > ecalTimeMin && *ecalT < ecalTimeMax &&
          *ecalE > ecalEnergyMin && *ecalE < ecalEnergyMax)) continue;

    // First count accepted hits in each layer, then apply event occupancies.
    int goodHits[2] = {0, 0};
    std::vector<size_t> candidates;
    for (size_t i = 0; i < n; ++i) {
      const int channel = int(id[i]);
      if (channel < 0 || channel >= NumCDetPaddles) continue;
      const int layer = channel / 1344;
      const double correctedX = CorrectCDetX(x[i], layer);
      const double dx = correctedX - *ecalX*z[i]/ECal_dist - xOffset;
      const double dy = y[i] - *ecalY*z[i]/ECal_dist - yOffset;
      if (!(le[i] >= cutLEMin/TDC_calib_to_ns && le[i] <= cutLEMax/TDC_calib_to_ns &&
            tot[i] >= cutToTMin/TDC_calib_to_ns && tot[i] <= cutToTMax/TDC_calib_to_ns &&
            mult[i] < TDCmult_cut && correctedX < xcut &&
            dx >= -xDiffMax && dx <= xDiffMax &&
            dy >= -1.2*CDet_y_half_length && dy <= 1.2*CDet_y_half_length)) continue;
      candidates.push_back(i);
      if (!check_bad(channel, suppressBad)) ++goodHits[layer];
    }
    if (candidates.empty() || goodHits[0] < minHits[0] || goodHits[0] > maxHits[0] ||
        goodHits[1] < minHits[1] || goodHits[1] > maxHits[1]) continue;
    ++goodEvents;

    double reference = 0;
    if (useReference && !check_bad(2696, suppressBad)) {
      for (size_t i = 0; i < rawID.GetSize(); ++i)
        if (int(rawID[i]) == 2696 && rawLE[i] > 0 && rawLE[i] <= 100/TDC_calib_to_ns &&
            rawToT[i] >= 0 && rawToT[i] <= 200/TDC_calib_to_ns &&
            int(rawLE[i]) > 0 && int(rawToT[i]) > 0)
          reference = rawLE[i]*TDC_calib_to_ns; // Last qualifying reference, as in master.
    }
    for (size_t i : candidates) {
      const int channel = int(id[i]), layer = channel / 1344;
      if (check_bad(channel, suppressBad) ||
          !((layerChoice == 1 && layer == 0) || (layerChoice == 2 && layer == 1) ||
            (layerChoice == 3 && goodHits[0] >= 1 && goodHits[1] >= 1))) continue;
      const double rawTime = le[i]*TDC_calib_to_ns, width = tot[i]*TDC_calib_to_ns;
      const double walkP1 = layer == 0 ? gTimeWalkP1_L1 : gTimeWalkP1_L2;
      const double totRef = layer == 0 ? gTimeWalkTotRef_L1 : gTimeWalkTotRef_L2;
      // The detector correlation's before view aligns pixels before examining
      // the ECal dependence. Other before plots retain completely raw timing.
      const double pixelTime = rawTime - reference + GetPixelToffsetCorr(channel);
      // The cut editor fills vPaddleGoodLe with this pixel-aligned time; that
      // vector is not updated by the later ECal, time-walk or run corrections.
      // Evaluate the polygon once in its drawing coordinates, not separately
      // on the raw and final displayed times. Uncut pixels retain the base cuts.
      if (pixelCuts[channel]) {
        ++polygonTestedHits;
        if (!pixelCuts[channel]->IsInside(width, pixelTime)) {
          ++polygonRejectedHits;
          if (channel >= base && channel < base+16) ++polygonRejectedBarHits;
          continue;
        }
      }
      // Same order and signs as stage 7; no cuts are reapplied to the new time.
      double time = pixelTime - (gECalFitP0 + gECalFitP1*(*ecalT)) + gECalDeltaShift;
      if (std::isfinite(width) && width > 0 && std::isfinite(totRef) && totRef > 0)
        time -= walkP1*(1/std::sqrt(width) - 1/std::sqrt(totRef));
      time += gGlobalTimingShift;
      hDT[0][channel]->Fill(*ecalT-rawTime);
      hDT[1][channel]->Fill(*ecalT-time);
      hCDetECal[0]->Fill(*ecalT, pixelTime);
      hCDetECal[1]->Fill(*ecalT, time);
      // Full-detector spectra are filled before restricting the bar plots.
      const int p = channel - base;
      if (p < 0 || p >= 16) continue;
      ++selectedBarHits;
      const double times[] = {rawTime, time};
      for (int state = 0; state < 2; ++state) {
        hLE[state][p]->Fill(times[state]);
        hLEToT[state][p]->Fill(width, times[state]);
      }
    }
  }
  if (processed != limit) {
    std::cerr << "[CDet before/after] Tree reading stopped early; check branches/input files.\n";
    return;
  }

  std::cout << "[CDet before/after] " << processed << " events; " << goodEvents
            << " good events; " << hCDetECal[0]->GetEntries() << " detector hits; "
            << selectedBarHits << " hits in bar " << bar
            << " (identical before/after sample).\n";
  if (loadedPixelCuts > 0)
    std::cout << "[CDet before/after] Polygons tested " << polygonTestedHits
              << " hits; rejected " << polygonRejectedHits << " detector hits ("
              << polygonRejectedBarHits << " in bar " << bar
              << "). Event occupancy cuts were evaluated before polygons.\n";
  if (savePlots && gSystem->mkdir(outputDirectory, true) != 0 &&
      gSystem->AccessPathName(outputDirectory)) return;

  std::array<std::vector<std::unique_ptr<TF1>>, 2> dtFits;
  std::array<std::vector<TString>, 2> dtFitNotes;
  double mu0[2]; TString referenceDescription[2];
  for (int state = 0; state < 2; ++state) {
    std::cout << "[CDet before/after] Fitting " << stage[state] << " detector DT spectra...\n";
    mu0[state] = CDetPixelTimingBeforeAfter::FitDTSpectra(
        hDT[state], fitMin[state], fitMax[state], dtFits[state],
        dtFitNotes[state], referenceDescription[state]);
    std::cout << "  mu_0 = " << mu0[state] << " ns (" << referenceDescription[state] << ")\n";
  }

  // Each LE/DT panel follows its own peak, including the displayed DT fit.
  // Statistics sit above the frame, so only 5% count-axis headroom is needed.
  // Keep the before/after 2D color scales matched.
  auto histogram = [&](int kind, int state, int p) -> TH1* {
    if (kind == 0) return hLE[state][p].get();
    if (kind == 1) return hLEToT[state][p].get();
    return hDT[state][base+p].get();
  };
  double displayMax[2][3][16];
  for (int kind = 0; kind < 3; ++kind)
    for (int p = 0; p < 16; ++p)
      for (int state = 0; state < 2; ++state) {
        double top = std::max(1.0, histogram(kind, state, p)->GetMaximum());
        if (kind == 1) top = std::max(top, histogram(kind, 1-state, p)->GetMaximum());
        if (kind == 2 && dtFits[state][base+p])
          top = std::max(top, dtFits[state][base+p]->GetMaximum());
        displayMax[state][kind][p] = (kind == 1 ? 1.0 : 1.05)*top;
      }
  const double correlationMax = std::max(1.0, std::max(hCDetECal[0]->GetMaximum(), hCDetECal[1]->GetMaximum()));
  const char *kind[] = {"le", "le_vs_tot", "ecal_cdet_dt"};
  const char *description[] = {"LE", "LE vs ToT", "ECal - CDet #Deltat"};
  const TString polygonLabel = loadedPixelCuts > 0
      ? TString::Format("%d saved polygons at pixel-aligned LE; same hits in both samples", loadedPixelCuts)
      : "Pixel polygons disabled; same hits in both samples";
  for (int state = 0; state < 2; ++state) {
    for (int plot = 0; plot < 3; ++plot) {
      const bool twoD = plot == 1;
      const TString name = TString::Format("CDet_%s_%s_%s", tag.Data(), kind[plot], stage[state]);
      TCanvas *canvas = new TCanvas(name, name, 1600, 820);
      TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(0.033);
      label.DrawLatex(0.035, 0.965, TString::Format("Run %d | Layer %d, bar %d | %s | %s timing calibration", run, pixel/1344+1, bar, description[plot], state ? "After full" : "Before"));
      label.SetTextSize(0.021);
      label.DrawLatex(0.035, 0.937, TString::Format("LE [%.3g, %.3g] ns; ToT [%.3g, %.3g] ns; ECal time (%.3g, %.3g) ns; |#Deltax| #leq %.3g m; |#Deltay| #leq 0.36 m", cutLEMin, cutLEMax, cutToTMin, cutToTMax, ecalTimeMin, ecalTimeMax, xDiffMax));
      label.DrawLatex(0.035, 0.910, polygonLabel);
      if (plot == 0)
        label.DrawLatex(0.035, 0.883, "Independent count-axis maxima; identical before/after x ranges and binning");
      if (plot == 2) {
        const TString reference = std::isfinite(mu0[state])
            ? TString::Format("#mu_{0} = %.2f ns (%s)", mu0[state], referenceDescription[state].Data())
            : "#mu_{0} = n/a (no reliable detector reference)";
        label.DrawLatex(0.035, 0.883, TString::Format("%s | Red: Gaussian + linear background; broad fit [%.3g, %.3g] ns | Independent count axes", reference.Data(), fitMin[state], fitMax[state]));
      }
      TPad *grid = new TPad(name+"_grid", "", 0, 0, 1, plot == 1 ? 0.887 : 0.860);
      grid->Draw(); grid->Divide(4, 2, 0.002, 0.002);
      // Original 4x4 rows two and three: retain their logical pixel labels.
      for (int p = 4; p < 12; ++p) {
        grid->cd(p-3);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(twoD ? 0.15 : 0.04);
        gPad->SetTopMargin(twoD ? 0.11 : plot == 2 ? 0.26 : 0.17);
        TH1 *hist = histogram(plot, state, p);
        hist->SetMinimum(0); hist->SetMaximum(displayMax[state][plot][p]);
        hist->DrawCopy(twoD ? "COLZ" : "HIST");
        if (plot == 2 && !IsUnusedPixel(base+p)) {
          TF1 *fit = dtFits[state][base+p].get();
          if (fit) fit->DrawCopy("SAME");
          // A compact header above the axes keeps fit text clear of the data.
          label.SetTextSize(0.035); label.SetTextAlign(11);
          const TString counts = hist->GetEntries() > 1
              ? TString::Format("N = %.0f   SD = %.2f ns", hist->GetEntries(), hist->GetStdDev())
              : TString::Format("N = %.0f   SD = n/a", hist->GetEntries());
          label.DrawLatex(0.16, 0.89, counts);
          if (fit) {
            label.DrawLatex(0.16, 0.84, TString::Format("#mu_{i} = %.2f #pm %.2f ns", fit->GetParameter(1), fit->GetParError(1)));
            label.DrawLatex(0.16, 0.79, TString::Format("#sigma_{fit} = %.2f #pm %.2f ns", std::fabs(fit->GetParameter(2)), fit->GetParError(2)));
            label.SetTextSize(0.032); label.SetTextAlign(31);
            label.DrawLatex(0.94, 0.89, TString::Format("#chi^{2}/NDF = %.2f", fit->GetChisquare()/fit->GetNDF()));
            if (dtFitNotes[state][base+p] == "individual_broad_fallback")
              label.DrawLatex(0.94, 0.84, "Broad fit");
          } else {
            label.DrawLatex(0.16, 0.84, "#mu_{i}, #sigma_{fit} = n/a");
            label.SetTextSize(0.032);
            label.DrawLatex(0.16, 0.79, dtFitNotes[state][base+p]);
          }
          label.SetTextAlign(11);
        } else if (!twoD) {
          // Keep LE statistics above the frame as well.
          label.SetTextSize(0.039); label.SetTextAlign(31);
          const TString counts = hist->GetEntries() > 1
              ? TString::Format("N = %.0f   SD = %.2f ns", hist->GetEntries(), hist->GetStdDev())
              : TString::Format("N = %.0f   SD = n/a", hist->GetEntries());
          label.DrawLatex(0.94, 0.87, counts);
          label.SetTextAlign(11);
        }
        if (IsUnusedPixel(base+p)) {
          label.SetTextSize(0.065); label.DrawLatex(0.27, 0.53, "Unused pixel");
        }
      }
      canvas->Update();
      if (savePlots) {
        const TString prefix = TString::Format("%s/CDet_run%d_bar%03d_%s_%s", outputDirectory, run, bar, kind[plot], stage[state]);
        canvas->SaveAs(prefix+".pdf"); canvas->SaveAs(prefix+".png");
      }
    }

    // One aggregate correlation per state, including every selected detector
    // hit. The selected bar controls only the six 4x2 canvases above.
    const TString name = TString::Format("CDet_%s_detector_cdet_t_vs_ecal_t_%s", tag.Data(), stage[state]);
    TCanvas *canvas = new TCanvas(name, name, 1600, 1000);
    TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(0.027);
    label.DrawLatex(0.04, 0.963, TString::Format("Run %d | Full detector | CDet t vs ECal t | %s", run, state ? "After: full timing calibration" : "Before: pixel offsets applied"));
    label.SetTextSize(0.017);
    label.DrawLatex(0.04, 0.931, TString::Format("LE [%.3g, %.3g] ns; ToT [%.3g, %.3g] ns; ECal time (%.3g, %.3g) ns; |#Deltax| #leq %.3g m; |#Deltay| #leq 0.36 m", cutLEMin, cutLEMax, cutToTMin, cutToTMax, ecalTimeMin, ecalTimeMax, xDiffMax));
    label.DrawLatex(0.04, 0.904, TString::Format("N = %.0f accepted hits | %s | Black points: mean CDet time in displayed LE bins", hCDetECal[state]->GetEntries(), layerChoice == 3 ? "Layers 1 and 2" : layerChoice == 1 ? "Layer 1" : "Layer 2"));
    label.DrawLatex(0.04, 0.877, polygonLabel);
    TPad *pad = new TPad(name+"_plot", "", 0, 0, 1, 0.853);
    pad->Draw(); pad->cd();
    pad->SetLeftMargin(0.12); pad->SetRightMargin(0.18);
    pad->SetBottomMargin(0.14); pad->SetTopMargin(0.04);
    TH2D *hist = hCDetECal[state].get();
    for (TAxis *axis : {hist->GetXaxis(), hist->GetYaxis(), hist->GetZaxis()}) {
      axis->SetLabelSize(0.036); axis->SetTitleSize(0.042);
    }
    hist->GetZaxis()->SetTitleOffset(1.0);
    hist->SetMinimum(0); hist->SetMaximum(correlationMax);
    hist->DrawCopy("COLZ");
    if (hist->GetEntries() > 0) {
      // Same ProfileX marker convention as cCDetTvsECalT in the master.
      // This pooled trend is a visual summary, not a new calibration fit.
      std::unique_ptr<TProfile> profile(hist->ProfileX(TString::Format("pCDetECal_%s_%s_detector", tag.Data(), stage[state]), 1, nLE));
      profile->SetDirectory(nullptr); profile->SetStats(false);
      profile->SetMarkerStyle(20); profile->SetMarkerSize(0.6);
      profile->SetMarkerColor(kBlack); profile->SetLineColor(kBlack);
      profile->DrawCopy("P E1 SAME");
    }
    canvas->Update();
    if (savePlots) {
      const TString prefix = TString::Format("%s/CDet_run%d_detector_cdet_t_vs_ecal_t_%s", outputDirectory, run, stage[state]);
      canvas->SaveAs(prefix+".pdf"); canvas->SaveAs(prefix+".png");
    }
  }
}

// Two independent samples from the same replay event pass. The 2D plot uses
// stored pair.* membership and ECal-minus-pair-mean timing. Its x coordinate
// is the L1-minus-L2 trajectory residual, i.e. -pair.trajectory_residual.
// The bar spectrum reproduces
// hCDetBar30ECalMinusCDet_ProjectedQuality in the good-pulse diagnostics and
// includes unpaired pulses. All plotted timing is already corrected, in ns.
// This entry point does not load or reapply the macro's calibration/polygons.
// The optional y propagation diagnostic follows plotCDetLayersTimeComp(),
// adding comparisons on the same selected sample. Set the final argument to
// zero to retain only the two replay-timing plots.
void Plot_CDet_PairDTvsDXAndBarTiming(
    const char *configFile = "CDet_run6077_projection.conf", int bar = 30,
    const char *outputDirectory = "cdet_pair_timing",
    const char *inputDirectory = nullptr, Long64_t eventsOverride = -2,
    double dxBinWidth = 0.002, double dxMin = -0.16, double dxMax = 0.16,
    double dtBinWidth = 1.0, double dtMin = -60.0, double dtMax = 30.0,
    double fitMin = -55.0, double fitMax = -10.0, bool savePlots = true,
    double yCorrectionRefractiveIndex = 1.59)
{
  ApplyCDetPlotStyle();
  const int nDX = std::isfinite(dxBinWidth) && dxBinWidth > 0 &&
      std::isfinite(dxMin) && std::isfinite(dxMax) && dxMax > dxMin
      ? int(std::ceil((dxMax-dxMin)/dxBinWidth)) : 0;
  const int nDT = std::isfinite(dtBinWidth) && dtBinWidth > 0 &&
      std::isfinite(dtMin) && std::isfinite(dtMax) && dtMax > dtMin
      ? int(std::ceil((dtMax-dtMin)/dtBinWidth)) : 0;
  if (bar < 0 || bar >= NumCDetPaddles/16 || nDX <= 0 || nDT <= 0 ||
      !std::isfinite(fitMin) || !std::isfinite(fitMax) || fitMin >= fitMax ||
      fitMin < dtMin || fitMax > dtMax || eventsOverride < -2 ||
      !std::isfinite(yCorrectionRefractiveIndex) || yCorrectionRefractiveIndex < 0 ||
      (yCorrectionRefractiveIndex > 0 && yCorrectionRefractiveIndex < 1) ||
      (savePlots && (!outputDirectory || !outputDirectory[0]))) {
    std::cerr << "[CDet pair/bar timing] Invalid bar, event limit, histogram/fit range, refractive index or output directory.\n";
    return;
  }
  const bool applyYCorrection = yCorrectionRefractiveIndex > 0;
  const double ySlopeMagnitude = yCorrectionRefractiveIndex/0.299792458; // ns/m
  TEnv env;
  if (!LoadCDetConfiguration(env, configFile, "CDet pair/bar timing")) return;
  const int run = env.GetValue("analysis.run_number", 6077);
  const Long64_t events = eventsOverride == -2
      ? env.GetValue("analysis.events", -1) : eventsOverride;
  const double energyMin = env.GetValue("analysis.ecal_energy_min", 3.0);
  const double energyMax = env.GetValue("analysis.ecal_energy_max", 4.5);
  if (run <= 0 || !std::isfinite(energyMin) || !std::isfinite(energyMax) ||
      energyMin >= energyMax) {
    std::cerr << "[CDet pair/bar timing] Invalid run or ECal energy interval.\n";
    return;
  }
  const char *directory = inputDirectory && inputDirectory[0]
      ? inputDirectory : gSystem->Getenv("OUT_DIR");
  if (!directory || !directory[0]) {
    std::cerr << "[CDet pair/bar timing] Set OUT_DIR or pass inputDirectory.\n";
    return;
  }
  TChain chain("T");
  AddRunFilesToChain(&chain, directory, run,
      env.GetValue("analysis.min_segment", -1), env.GetValue("analysis.max_segment", -1));
  const Long64_t total = chain.GetEntries();
  if (total <= 0) return;
  const Long64_t limit = events > 0 ? std::min(events, total) : total;
  chain.LoadTree(0);
  for (const char *name : {"earm.ecal.e", "earm.cdet.pulse.pmtnum",
       "earm.cdet.pulse.ecal_residual", "earm.cdet.pulse.tdc_tot_ns",
       "earm.cdet.pulse.calib_valid", "earm.cdet.pulse.ecal_eligible",
       "earm.cdet.pulse.spatial_pass", "earm.cdet.pulse.broad_quality_pass",
       "earm.cdet.pulse.x_corr", "earm.cdet.pulse.ecal_x_proj",
       "earm.cdet.pair.pulse_index_l1", "earm.cdet.pair.pulse_index_l2",
       "earm.cdet.pair.ecal_residual"}) {
    if (!chain.GetBranch(name)) {
      std::cerr << "[CDet pair/bar timing] Missing required branch: " << name << ". No plots produced.\n";
      return;
    }
  }
  if (applyYCorrection)
    for (const char *name : {"earm.cdet.pulse.y", "earm.cdet.pulse.ecal_y_proj"})
      if (!chain.GetBranch(name)) {
        std::cerr << "[CDet pair/bar timing] Y correction requires branch: " << name << ". No plots produced.\n";
        return;
      }
  TTreeReader reader(&chain);
  TTreeReaderValue<double> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderArray<double> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<double> pulseDT(reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<double> tot(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<double> calibValid(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<double> ecalEligible(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<double> spatialPass(reader, "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<double> broadQuality(reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<double> x(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<double> projectedX(reader, "earm.cdet.pulse.ecal_x_proj");
  TTreeReaderArray<double> pairIndexL1(reader, "earm.cdet.pair.pulse_index_l1");
  TTreeReaderArray<double> pairIndexL2(reader, "earm.cdet.pair.pulse_index_l2");
  TTreeReaderArray<double> pairDT(reader, "earm.cdet.pair.ecal_residual");
  std::unique_ptr<TTreeReaderArray<double>> yCenter, projectedY;
  if (applyYCorrection) {
    yCenter.reset(new TTreeReaderArray<double>(reader, "earm.cdet.pulse.y"));
    projectedY.reset(new TTreeReaderArray<double>(reader, "earm.cdet.pulse.ecal_y_proj"));
    std::cout << "[CDet pair/bar timing] Y comparison enabled: n = "
              << yCorrectionRefractiveIndex << ", |dt/dy| = " << ySlopeMagnitude
              << " ns/m; same pulses and stored pairs before/after.\n";
  }
  // pulse.y is the database half-bar center (constant across its paddles).
  // Use the stored projection itself, not ecal_y_residual, which additionally
  // subtracts the spatial-selection offset. This matches the master's model:
  // t_y = t - b_side*(y_proj-y_center), hence (ECal-t)_y = (ECal-t) + b_side*dy.
  auto yResidualShift = [&](size_t i) {
    if (!applyYCorrection) return 0.0;
    if (!std::isfinite((*yCenter)[i]) || !std::isfinite((*projectedY)[i]) ||
        std::fabs((*yCenter)[i]) >= 900 || std::fabs((*projectedY)[i]) >= 900)
      return std::numeric_limits<double>::quiet_NaN();
    const int side = (int(pixel[i])%1344)/672;
    const double signedSlope = side == 0 ? -ySlopeMagnitude : ySlopeMagnitude;
    return signedSlope*((*projectedY)[i]-(*yCenter)[i]);
  };

  static unsigned int invocation = 0;
  const TString tag = TString::Format("run%d_bar%03d_%u", run, bar, ++invocation);
  TH2D hPairDTvsDX("hCDetPairDTvsDX_"+tag, ";(x_{1,corr}-x_{2,corr}) - (x_{ECal}/z_{ECal})(z_{1}-z_{2}) (m);t_{ECal} - <t_{CDet,corr}>_{pair} (ns);Pairs", nDX, dxMin, dxMax, nDT, dtMin, dtMax);
  TH1D hBarDT("hCDetBarProjectedQualityDT_"+tag, TString::Format(";t_{ECal} - t_{CDet,corr} (ns);Pulses / %.3g ns", (dtMax-dtMin)/nDT), nDT, dtMin, dtMax);
  TH2D hPairYDTvsDX("hCDetPairYCorrectedDTvsDX_"+tag, ";(x_{1,corr}-x_{2,corr}) - (x_{ECal}/z_{ECal})(z_{1}-z_{2}) (m);t_{ECal} - <t_{CDet,y corr}>_{pair} (ns);Pairs", nDX, dxMin, dxMax, nDT, dtMin, dtMax);
  TH1D hBarYDT("hCDetBarYCorrectedDT_"+tag, TString::Format(";t_{ECal} - t_{CDet} (ns);Pulses / %.3g ns", (dtMax-dtMin)/nDT), nDT, dtMin, dtMax);
  for (TH1 *hist : {static_cast<TH1*>(&hPairDTvsDX), static_cast<TH1*>(&hBarDT), static_cast<TH1*>(&hPairYDTvsDX), static_cast<TH1*>(&hBarYDT)}) {
    hist->SetDirectory(nullptr); hist->SetStats(false);
    hist->SetStatOverflows(TH1::kConsider);
  }
  hBarDT.SetLineColor(kBlack); hBarDT.SetLineWidth(2);
  hBarYDT.SetLineColor(kBlue+1); hBarYDT.SetLineWidth(2);

  Long64_t processed = 0, energyEvents = 0, pairEvents = 0, barEvents = 0;
  Long64_t malformedPairs = 0;
  while (processed < limit && reader.Next()) {
    ++processed;
    if (processed % 1000 == 0)
      std::cout << "[CDet pair/bar timing] Run " << run << ": entry "
                << processed << "/" << limit << " (chain entries "
                << total << ")" << std::endl;
    const size_t n = pixel.GetSize(), np = pairDT.GetSize();
    if (pulseDT.GetSize() != n || tot.GetSize() != n || calibValid.GetSize() != n ||
        ecalEligible.GetSize() != n || spatialPass.GetSize() != n ||
        broadQuality.GetSize() != n || x.GetSize() != n || projectedX.GetSize() != n ||
        pairIndexL1.GetSize() != np || pairIndexL2.GetSize() != np ||
        (applyYCorrection && (yCenter->GetSize() != n || projectedY->GetSize() != n))) {
      std::cerr << "[CDet pair/bar timing] Mismatched arrays at entry " << reader.GetCurrentEntry() << ". No plots produced.\n";
      return;
    }
    if (!std::isfinite(*ecalEnergy) || *ecalEnergy < energyMin || *ecalEnergy > energyMax) continue;
    ++energyEvents;
    bool hasPair = false, hasBarPulse = false;
    for (size_t p = 0; p < np; ++p) {
      const double index1 = pairIndexL1[p], index2 = pairIndexL2[p];
      if (!std::isfinite(index1) || !std::isfinite(index2) ||
          index1 < 0 || index2 < 0 || index1 >= n || index2 >= n ||
          index1 != std::floor(index1) || index2 != std::floor(index2)) {
        ++malformedPairs; continue;
      }
      const size_t i1 = size_t(index1), i2 = size_t(index2);
      // Both differences use L1 minus L2: projectedX1-projectedX2 is
      // (x_ECal/z_ECal)*(z1-z2), with signed dz (not the positive spacing).
      const double dx = (x[i1]-x[i2]) - (projectedX[i1]-projectedX[i2]);
      if (!(pixel[i1] >= 0 && pixel[i1] < 1344 && pixel[i2] >= 1344 && pixel[i2] < NumCDetPaddles) ||
          !std::isfinite(dx) || !std::isfinite(pairDT[p])) {
        ++malformedPairs; continue;
      }
      // Correct the two members separately; opposite-side pairs have different
      // propagation signs. Do not rerun pairing or its timing/ellipse cuts.
      const double pairShift = 0.5*(yResidualShift(i1)+yResidualShift(i2));
      if (!std::isfinite(pairShift)) {
        std::cerr << "[CDet pair/bar timing] Invalid y geometry in selected pair at entry " << reader.GetCurrentEntry() << ". No plots produced.\n";
        return;
      }
      hPairDTvsDX.Fill(dx, pairDT[p]);
      if (applyYCorrection) hPairYDTvsDX.Fill(dx, pairDT[p]+pairShift);
      hasPair = true;
    }
    // Exact bar population from Plot_CDet_GoodPulseCandidates_AllTDC.C:
    // energy interval is inclusive; ECal position/time eligibility, projection
    // and broad pulse quality come from replay flags. No pair or ellipse cut,
    // additional time window, two-layer occupancy, or local pixel polygon.
    for (size_t i = 0; i < n; ++i) {
      if (!(calibValid[i] > 0.5 && ecalEligible[i] > 0.5 &&
            spatialPass[i] > 0.5 && broadQuality[i] > 0.5) ||
          !std::isfinite(pixel[i]) || !std::isfinite(pulseDT[i]) || !std::isfinite(tot[i])) continue;
      const int channel = int(std::lround(pixel[i]));
      if (channel < bar*16 || channel >= (bar+1)*16) continue;
      const double shift = yResidualShift(i);
      if (!std::isfinite(shift)) {
        std::cerr << "[CDet pair/bar timing] Invalid y geometry in selected bar pulse at entry " << reader.GetCurrentEntry() << ". No plots produced.\n";
        return;
      }
      hBarDT.Fill(pulseDT[i]);
      if (applyYCorrection) hBarYDT.Fill(pulseDT[i]+shift);
      hasBarPulse = true;
    }
    if (hasPair) ++pairEvents;
    if (hasBarPulse) ++barEvents;
  }
  if (processed != limit) {
    std::cerr << "[CDet pair/bar timing] Tree reading stopped early; check branches/input files. No plots produced.\n";
    return;
  }
  std::cout << "[CDet pair/bar timing] " << processed << " entries; " << energyEvents
            << " ECal-energy-selected events; " << hPairDTvsDX.GetEntries()
            << " stored pairs in " << pairEvents << " events; " << hBarDT.GetEntries()
            << " projection + quality pulses in bar " << bar << " from " << barEvents
            << " events; " << malformedPairs << " malformed pairs skipped.\n";
  if (savePlots && gSystem->mkdir(outputDirectory, true) != 0 &&
      gSystem->AccessPathName(outputDirectory)) return;

  // Same Gaussian + linear-background diagnostic and default fit/seed windows
  // as the historical Bar-30 panel. A failed fit is not reported as a result.
  auto fitBarSpectrum = [&](TH1D& hist, const TString& name, int color) {
    std::unique_ptr<TF1> fit;
    if (hist.GetEntries() < 20) return fit;
    const double seedMin = std::max(fitMin, -40.0), seedMax = std::min(fitMax, -15.0);
    const int low = std::max(1, hist.FindFixBin(seedMin < seedMax ? seedMin : fitMin));
    const int high = std::min(nDT, hist.FindFixBin(seedMin < seedMax ? seedMax : fitMax));
    int peak = low;
    for (int bin = low+1; bin <= high; ++bin)
      if (hist.GetBinContent(bin) > hist.GetBinContent(peak)) peak = bin;
    const int fitLow = std::max(1, hist.FindFixBin(fitMin));
    const int fitHigh = std::min(nDT, hist.FindFixBin(fitMax));
    const double background = 0.5*(hist.GetBinContent(fitLow)+hist.GetBinContent(fitHigh));
    fit.reset(new TF1(name, "gaus(0)+pol1(3)", fitMin, fitMax));
    fit->SetParameters(std::max(1.0, hist.GetBinContent(peak)-background), hist.GetBinCenter(peak), 3.0, background, 0.0);
    const int status = hist.Fit(fit.get(), "RQN0");
    if (status != 0 || !std::isfinite(fit->GetParameter(1)) ||
        !std::isfinite(fit->GetParameter(2)) || fit->GetParameter(0) <= 0 ||
        std::fabs(fit->GetParameter(2)) <= 0 || fit->GetNDF() <= 0) fit.reset();
    else { fit->SetLineColor(color); fit->SetLineWidth(2); fit->SetNpx(500); }
    return fit;
  };
  auto fit = fitBarSpectrum(hBarDT, "fCDetBarProjectedQualityDT_"+tag, kRed+1);
  std::unique_ptr<TF1> yFit;
  if (applyYCorrection) {
    yFit = fitBarSpectrum(hBarYDT, "fCDetBarYCorrectedDT_"+tag, kMagenta+1);
    std::cout << "[CDet pair/bar timing] Bar " << bar << " y comparison: "
              << hBarDT.GetEntries() << " -> " << hBarYDT.GetEntries()
              << " pulses; mean " << hBarDT.GetMean() << " -> " << hBarYDT.GetMean()
              << " ns; SD (including tails) " << hBarDT.GetStdDev() << " -> "
              << hBarYDT.GetStdDev() << " ns.\n";
    if (fit && yFit)
      std::cout << "[CDet pair/bar timing] Bar peak sigma " << std::fabs(fit->GetParameter(2))
                << " -> " << std::fabs(yFit->GetParameter(2)) << " ns.\n";
  }

  for (int plot = 0; plot < 2; ++plot) {
    const TString suffix = plot == 0 ? "pair_ecal_cdet_dt_vs_dx"
        : TString::Format("bar%03d_ecal_cdet_dt_projected_quality", bar);
    const TString name = "CDet_"+tag+"_"+suffix;
    TCanvas *canvas = new TCanvas(name, name, 1400, 900);
    TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(0.030);
    label.DrawLatex(0.04, 0.964, plot == 0
        ? TString::Format("Run %d | Stored CDet pairs | ECal-CDet #Deltat vs x trajectory residual", run)
        : TString::Format("Run %d | ECal projection + pulse-quality selection | Bar %d", run, bar));
    label.SetTextSize(0.021);
    label.DrawLatex(0.04, 0.928, TString::Format("ECal energy [%.3g, %.3g] GeV | Timing and alignment from replay", energyMin, energyMax));
    label.DrawLatex(0.04, 0.896, plot == 0
        ? TString::Format("N = %.0f stored pairs in %lld events | Full detector | One entry per pair", hPairDTvsDX.GetEntries(), pairEvents)
        : TString::Format("N = %.0f pulses in %lld events | All pixels %d-%d | Pairing not required", hBarDT.GetEntries(), barEvents, 16*bar, 16*bar+15));
    if (plot == 1) {
      label.SetTextSize(0.020);
      label.DrawLatex(0.04, 0.864, fit
          ? TString::Format("Gaussian + linear background [%.3g, %.3g] ns: #mu = %.2f #pm %.2f ns; #sigma = %.2f #pm %.2f ns", fitMin, fitMax, fit->GetParameter(1), fit->GetParError(1), std::fabs(fit->GetParameter(2)), fit->GetParError(2))
          : "Gaussian + linear background: fit unavailable (low statistics or failed fit)");
    }
    TPad *pad = new TPad(name+"_plot", "", 0, 0, 1, plot == 0 ? 0.865 : 0.832);
    pad->Draw(); pad->cd();
    pad->SetLeftMargin(0.13); pad->SetRightMargin(plot == 0 ? 0.16 : 0.04);
    pad->SetBottomMargin(0.15); pad->SetTopMargin(0.04);
    if (plot == 0) {
      hPairDTvsDX.GetXaxis()->SetTitleSize(0.037);
      hPairDTvsDX.GetXaxis()->SetTitleOffset(1.4);
      hPairDTvsDX.GetYaxis()->SetTitleSize(0.040);
      hPairDTvsDX.DrawCopy("COLZ");
    } else {
      hBarDT.SetMinimum(0);
      hBarDT.SetMaximum(1.05*std::max({1.0, hBarDT.GetMaximum(), fit ? fit->GetMaximum() : 0.0}));
      hBarDT.DrawCopy("HIST");
      if (fit) fit->DrawCopy("SAME");
    }
    canvas->Update();
    if (savePlots) {
      const TString prefix = TString::Format("%s/CDet_run%d_%s", outputDirectory, run, suffix.Data());
      canvas->SaveAs(prefix+".pdf"); canvas->SaveAs(prefix+".png");
    }
  }

  if (!applyYCorrection) return;
  const double colorMax = std::max({1.0, hPairDTvsDX.GetMaximum(), hPairYDTvsDX.GetMaximum()});
  const TString pairName = "CDet_"+tag+"_pair_y_comparison";
  TCanvas *pairCanvas = new TCanvas(pairName, pairName, 1900, 950);
  TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(0.029);
  label.DrawLatex(0.035, 0.961, TString::Format("Run %d | Stored-pair ECal-CDet #Deltat vs x trajectory residual | y propagation comparison", run));
  label.SetTextSize(0.021);
  label.DrawLatex(0.035, 0.921, TString::Format("Same %.0f pairs | ECal energy [%.3g, %.3g] GeV | n = %.3g; |dt/dy| = %.3f ns/m | Matched color scales", hPairDTvsDX.GetEntries(), energyMin, energyMax, yCorrectionRefractiveIndex, ySlopeMagnitude));
  TPad *pairGrid = new TPad(pairName+"_grid", "", 0, 0, 1, 0.89);
  pairGrid->Draw(); pairGrid->Divide(2, 1, 0.005, 0.005);
  for (int state = 0; state < 2; ++state) {
    pairGrid->cd(state+1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.18);
    gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.09);
    TH2D& hist = state ? hPairYDTvsDX : hPairDTvsDX;
    hist.GetXaxis()->SetTitleSize(0.031); hist.GetYaxis()->SetTitleSize(0.036);
    hist.GetXaxis()->SetTitleOffset(1.6);
    hist.GetZaxis()->SetTitleSize(0.038); hist.GetZaxis()->SetTitleOffset(1.4);
    for (TAxis *axis : {hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()}) axis->SetLabelSize(0.037);
    hist.SetMinimum(0); hist.SetMaximum(colorMax); hist.DrawCopy("COLZ");
    label.SetTextSize(0.033);
    label.DrawLatex(0.16, 0.946, state ? "After y propagation correction" : "Replay timing before y correction");
  }
  pairCanvas->Update();

  const TString barName = "CDet_"+tag+"_bar_y_comparison";
  TCanvas *barCanvas = new TCanvas(barName, barName, 1400, 950);
  label.SetTextSize(0.030);
  label.DrawLatex(0.04, 0.964, TString::Format("Run %d | Bar %d | ECal projection + pulse quality | y propagation comparison", run, bar));
  label.SetTextSize(0.021);
  label.DrawLatex(0.04, 0.927, TString::Format("Same %.0f pulses in %lld events | ECal energy [%.3g, %.3g] GeV | n = %.3g; |dt/dy| = %.3f ns/m", hBarDT.GetEntries(), barEvents, energyMin, energyMax, yCorrectionRefractiveIndex, ySlopeMagnitude));
  for (int state = 0; state < 2; ++state) {
    TF1 *currentFit = state ? yFit.get() : fit.get();
    label.DrawLatex(0.04, 0.893-0.033*state, currentFit
        ? TString::Format("%s: #mu = %.2f #pm %.2f ns; #sigma = %.2f #pm %.2f ns | Gaussian + linear background [%.3g, %.3g] ns", state ? "After y" : "Before y", currentFit->GetParameter(1), currentFit->GetParError(1), std::fabs(currentFit->GetParameter(2)), currentFit->GetParError(2), fitMin, fitMax)
        : TString::Format("%s: fit unavailable", state ? "After y" : "Before y"));
  }
  TPad *barPad = new TPad(barName+"_plot", "", 0, 0, 1, 0.830);
  barPad->Draw(); barPad->cd();
  barPad->SetLeftMargin(0.13); barPad->SetRightMargin(0.04);
  barPad->SetBottomMargin(0.15); barPad->SetTopMargin(0.04);
  hBarDT.GetXaxis()->SetTitle("t_{ECal} - t_{CDet} (ns)");
  hBarDT.SetMaximum(1.05*std::max({1.0, hBarDT.GetBinContent(hBarDT.GetMaximumBin()), hBarYDT.GetMaximum(), fit ? fit->GetMaximum() : 0.0, yFit ? yFit->GetMaximum() : 0.0}));
  TH1 *beforeDraw = hBarDT.DrawCopy("HIST"), *afterDraw = hBarYDT.DrawCopy("HIST SAME");
  TF1 *beforeFitDraw = fit ? fit->DrawCopy("SAME") : nullptr;
  TF1 *afterFitDraw = yFit ? yFit->DrawCopy("SAME") : nullptr;
  TLegend legend(0.61, 0.73, 0.92, 0.92);
  legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.027);
  legend.AddEntry(beforeDraw, "Before y correction", "l");
  legend.AddEntry(afterDraw, "After y correction", "l");
  if (beforeFitDraw) legend.AddEntry(beforeFitDraw, "Before-y peak fit", "l");
  if (afterFitDraw) legend.AddEntry(afterFitDraw, "After-y peak fit", "l");
  legend.DrawClone();
  barCanvas->Update();
  if (savePlots) {
    const TString pairPrefix = TString::Format("%s/CDet_run%d_pair_ecal_cdet_dt_vs_dx_y_comparison", outputDirectory, run);
    pairCanvas->SaveAs(pairPrefix+".pdf"); pairCanvas->SaveAs(pairPrefix+".png");
    const TString barPrefix = TString::Format("%s/CDet_run%d_bar%03d_ecal_cdet_dt_y_comparison", outputDirectory, run, bar);
    barCanvas->SaveAs(barPrefix+".pdf"); barCanvas->SaveAs(barPrefix+".png");
  }
}
