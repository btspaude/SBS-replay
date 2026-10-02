// Eight bar canvases for a cross-target timing-calibration comparison.
// ROOT (from scripts/cdet, in a fresh session):
//   .L Plot_CDet_PixelTimingBeforeAfter.C+
//   Plot_CDet_PixelTimingBeforeAfter("CDet_run5710_projection.conf", 485);
// Any logical pixel 0..2687 selects its containing 16-pixel bar (pixel/16).
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
// No manual polygons, offset updates, or calibration-file writes occur here.
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
    const char *configFile = "CDet_run5710_projection.conf", int pixel = 485,
    const char *outputDirectory = "cdet_timing_before_after",
    const char *inputDirectory = nullptr, Long64_t eventsOverride = -2,
    double leBinWidth = 1, double leMin = 0, double leMax = 60,
    double totBinWidth = 1, double totMin = 0, double totMax = 40,
    bool savePlots = true,
    double dtBinWidth = 1, double dtMin = -40, double dtMax = 10,
    double ecalBinWidth = 1, double ecalMin = 5, double ecalMax = 40,
    double dtBeforeFitMin = -30, double dtBeforeFitMax = 10,
    double dtAfterFitMin = -30, double dtAfterFitMax = 10)
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
  std::array<std::array<std::unique_ptr<TH2D>, 16>, 2> hLEToT, hCDetECal;
  const char *stage[] = {"before", "after"};
  for (int state = 0; state < 2; ++state) {
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
      hCDetECal[state][p].reset(new TH2D(TString::Format("hCDetECal_%s_%s_p%d", tag.Data(), stage[state], base+p), TString::Format("Pixel %d;ECal ADC time (ns);CDet LE time (ns)", base+p), nECal, ecalMin, ecalMax, nLE, leMin, leMax));
      for (TH1 *hist : {static_cast<TH1*>(hLE[state][p].get()), static_cast<TH1*>(hLEToT[state][p].get()), static_cast<TH1*>(hCDetECal[state][p].get())}) {
        hist->SetDirectory(nullptr); hist->SetStats(false);
        hist->SetStatOverflows(TH1::kConsider);
      }
      hLE[state][p]->SetLineWidth(2);
      hLE[state][p]->SetLineColor(state ? kBlue+1 : kGray+2);
    }
  }

  Long64_t processed = 0, goodEvents = 0, selectedBarHits = 0;
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
    // No additional elastic, pair, ellipse, or saved polygon cut: neither the
    // master's vGoodLe population nor plotPaddles applies those selections.
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
      // Same order and signs as stage 7; no cuts are reapplied to the new time.
      double time = rawTime - reference + GetPixelToffsetCorr(channel);
      time = time - (gECalFitP0 + gECalFitP1*(*ecalT)) + gECalDeltaShift;
      if (std::isfinite(width) && width > 0 && std::isfinite(totRef) && totRef > 0)
        time -= walkP1*(1/std::sqrt(width) - 1/std::sqrt(totRef));
      time += gGlobalTimingShift;
      hDT[0][channel]->Fill(*ecalT-rawTime);
      hDT[1][channel]->Fill(*ecalT-time);
      // Full-detector DT spectra define mu_0; the other plots need one bar only.
      const int p = channel - base;
      if (p < 0 || p >= 16) continue;
      ++selectedBarHits;
      const double times[] = {rawTime, time};
      for (int state = 0; state < 2; ++state) {
        hLE[state][p]->Fill(times[state]);
        hLEToT[state][p]->Fill(width, times[state]);
        hCDetECal[state][p]->Fill(*ecalT, times[state]);
      }
    }
  }
  if (processed != limit) {
    std::cerr << "[CDet before/after] Tree reading stopped early; check branches/input files.\n";
    return;
  }

  std::cout << "[CDet before/after] " << processed << " events; " << goodEvents
            << " good events; " << selectedBarHits << " hits in bar " << bar
            << " (identical before/after sample).\n";
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

  // Matched axes/scales and small labels make before/after directly comparable.
  auto histogram = [&](int kind, int state, int p) -> TH1* {
    if (kind == 0) return hLE[state][p].get();
    if (kind == 1) return hLEToT[state][p].get();
    if (kind == 2) return hDT[state][base+p].get();
    return hCDetECal[state][p].get();
  };
  double displayMax[4][16];
  for (int kind = 0; kind < 4; ++kind)
    for (int p = 0; p < 16; ++p) {
      double top = std::max(1.0,
          std::max(histogram(kind, 0, p)->GetMaximum(), histogram(kind, 1, p)->GetMaximum()));
      if (kind == 2)
        for (int state = 0; state < 2; ++state)
          if (dtFits[state][base+p]) top = std::max(top, dtFits[state][base+p]->GetMaximum());
      displayMax[kind][p] = (kind % 2 ? 1.0 : kind == 2 ? 1.65 : 1.25) * top;
    }
  const char *kind[] = {"le", "le_vs_tot", "ecal_cdet_dt", "cdet_t_vs_ecal_t"};
  const char *description[] = {"LE", "LE vs ToT", "ECal - CDet #Deltat", "CDet t vs ECal t"};
  for (int state = 0; state < 2; ++state) {
    for (int plot = 0; plot < 4; ++plot) {
      const bool twoD = plot == 1 || plot == 3;
      const TString name = TString::Format("CDet_%s_%s_%s", tag.Data(), kind[plot], stage[state]);
      TCanvas *canvas = new TCanvas(name, name, 1600, 1200);
      TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(0.024);
      label.DrawLatex(0.035, 0.965, TString::Format("Run %d | Layer %d, bar %d | %s | %s timing calibration", run, pixel/1344+1, bar, description[plot], state ? "After full" : "Before"));
      label.SetTextSize(0.015);
      label.DrawLatex(0.035, 0.937, TString::Format("Same hits: LE [%.3g, %.3g] ns; ToT [%.3g, %.3g] ns; ECal time (%.3g, %.3g) ns; |#Deltax| #leq %.3g m; |#Deltay| #leq 0.36 m", cutLEMin, cutLEMax, cutToTMin, cutToTMax, ecalTimeMin, ecalTimeMax, xDiffMax));
      if (plot == 2) {
        const TString reference = std::isfinite(mu0[state])
            ? TString::Format("#mu_{0} = %.2f ns (%s)", mu0[state], referenceDescription[state].Data())
            : "#mu_{0} = n/a (no reliable detector reference)";
        label.DrawLatex(0.035, 0.910, TString::Format("%s | %s sample | Red: Gaussian + linear background; broad fit [%.3g, %.3g] ns", reference.Data(), stage[state], fitMin[state], fitMax[state]));
      }
      TPad *grid = new TPad(name+"_grid", "", 0, 0, 1, plot == 2 ? 0.887 : 0.915);
      grid->Draw(); grid->Divide(4, 4, 0.002, 0.002);
      for (int p = 0; p < 16; ++p) {
        grid->cd(p+1);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(twoD ? 0.15 : 0.04); gPad->SetTopMargin(0.11);
        TH1 *hist = histogram(plot, state, p);
        hist->SetMinimum(0); hist->SetMaximum(displayMax[plot][p]);
        hist->DrawCopy(twoD ? "COLZ" : "HIST");
        if (plot == 3 && hist->GetEntries() > 0) {
          // Same ProfileX marker convention as cCDetTvsECalT in the master.
          // This is a visual mean trend, not a calibration fit or update.
          std::unique_ptr<TProfile> profile(hCDetECal[state][p]->ProfileX(TString::Format("pCDetECal_%s_%s_p%d", tag.Data(), stage[state], base+p), 1, nLE));
          profile->SetDirectory(nullptr); profile->SetStats(false);
          profile->SetMarkerStyle(20); profile->SetMarkerSize(0.6);
          profile->SetMarkerColor(kBlack); profile->SetLineColor(kBlack);
          profile->DrawCopy("P E1 SAME");
        }
        if (plot == 2 && !IsUnusedPixel(base+p)) {
          TF1 *fit = dtFits[state][base+p].get();
          if (fit) fit->DrawCopy("SAME");
          label.SetTextSize(0.039); label.SetTextAlign(31);
          const TString counts = hist->GetEntries() > 1
              ? TString::Format("N = %.0f   SD = %.2f ns", hist->GetEntries(), hist->GetStdDev())
              : TString::Format("N = %.0f   SD = n/a", hist->GetEntries());
          label.DrawLatex(0.94, 0.855, counts);
          if (fit) {
            label.DrawLatex(0.94, 0.795, TString::Format("#mu_{i} = %.2f #pm %.2f ns", fit->GetParameter(1), fit->GetParError(1)));
            label.DrawLatex(0.94, 0.735, TString::Format("#sigma_{fit} = %.2f #pm %.2f ns", std::fabs(fit->GetParameter(2)), fit->GetParError(2)));
            label.DrawLatex(0.94, 0.675, TString::Format("#chi^{2}/NDF = %.2f%s", fit->GetChisquare()/fit->GetNDF(), dtFitNotes[state][base+p] == "individual_broad_fallback" ? " (broad fit)" : ""));
          } else {
            label.DrawLatex(0.94, 0.795, "#mu_{i}, #sigma_{fit} = n/a");
            label.SetTextSize(0.035);
            label.DrawLatex(0.94, 0.735, dtFitNotes[state][base+p]);
          }
          label.SetTextAlign(11);
        } else if (!twoD) {
          label.SetTextSize(0.044);
          label.DrawLatex(0.60, 0.82, TString::Format("N = %.0f", hist->GetEntries()));
          if (hist->GetEntries() > 1)
            label.DrawLatex(0.60, 0.755, TString::Format("SD = %.2f ns", hist->GetStdDev()));
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
  }
}
