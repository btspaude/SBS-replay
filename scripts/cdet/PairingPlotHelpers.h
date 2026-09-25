#ifndef CDET_PAIRING_PLOT_HELPERS_H
#define CDET_PAIRING_PLOT_HELPERS_H

#include <TH1D.h>
#include <TPaveText.h>
#include <TString.h>
#include <cmath>
#include <iostream>

namespace CDetPairingPlots {
inline int Bins(double width, double low, double high) {
  if (!std::isfinite(width) || !std::isfinite(low) || !std::isfinite(high) ||
      width <= 0 || high <= low || (high-low)/width > 1000000)
    return 0;
  return static_cast<int>(std::ceil((high-low)/width));
}

inline void Prepare(TH1D &h) {
  h.SetDirectory(nullptr);
  h.SetStats(false);
  // Preserve full selected-sample moments even outside the display bounds.
  h.SetStatOverflows(TH1::kConsider);
  h.SetLineWidth(2);
}

inline void Report(const TH1D &h) {
  std::cout << h.GetName() << ": N=" << h.GetEntries();
  if (h.GetEntries() >= 2)
    std::cout << ", mean=" << h.GetMean() << " ns, std dev=" << h.GetStdDev()
              << " +/- " << h.GetStdDevError() << " ns";
  else
    std::cout << ", insufficient entries for a width";
  std::cout << ", underflow=" << h.GetBinContent(0)
            << ", overflow=" << h.GetBinContent(h.GetNbinsX()+1) << '\n';
}

inline TString WidthLabel(const TH1D &h) {
  if (h.GetEntries() < 2)
    return "insufficient entries";
  return TString::Format("SD = %.3f #pm %.3f ns", h.GetStdDev(), h.GetStdDevError());
}

inline void Annotate(const TH1D &h) {
  auto *box = new TPaveText(0.48, 0.71, 0.89, 0.89, "NDC");
  box->SetFillColor(0);
  box->SetBorderSize(0);
  box->AddText(TString::Format("N = %.0f", h.GetEntries()));
  box->AddText(WidthLabel(h));
  box->AddText(TString::Format("Under/overflow: %.0f / %.0f",
      h.GetBinContent(0), h.GetBinContent(h.GetNbinsX()+1)));
  box->Draw();
}

inline bool ValidCuts(double totMin, double totMax, double dtMax, double radius) {
  return std::isfinite(totMin) && std::isfinite(totMax) &&
      std::isfinite(dtMax) && std::isfinite(radius) && totMin >= 0 &&
      totMax > totMin && (dtMax == -1 || dtMax > 0) &&
      (radius == -1 || radius > 0);
}

inline bool PassCuts(double tot1, double tot2, double dt, double score,
                     double totMin, double totMax, double dtMax, double radius) {
  return std::isfinite(tot1) && std::isfinite(tot2) && std::isfinite(dt) &&
      tot1 >= totMin && tot1 <= totMax && tot2 >= totMin && tot2 <= totMax &&
      (dtMax < 0 || std::fabs(dt) <= dtMax) &&
      (radius < 0 || (std::isfinite(score) && score >= 0 && score <= radius*radius));
}

inline bool Index(double value, size_t size, size_t &index) {
  if (!std::isfinite(value) || value < 0 || value >= static_cast<double>(size) ||
      std::floor(value) != value)
    return false;
  index = static_cast<size_t>(value);
  return true;
}

inline TString CutLabel(double totMin, double totMax, double dtMax, double radius) {
  return TString::Format("Stored pairs; member ToT [%.3g, %.3g] ns; |dt| max %.3g ns; R max %.3g (-1: off)",
                         totMin, totMax, dtMax, radius);
}
} // namespace CDetPairingPlots
#endif
