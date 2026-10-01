#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TFile.h>
#include <TH3D.h>
#include <TMarker3DBox.h>
#include <TMath.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TRotation.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>
#include <TView.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace CDetFTROIEventDisplay {

// Run 5711 was recorded 2025-08-07 and therefore uses the GEP-3 database
// blocks. Distances are metres and angles are degrees. These values come from
// DB/db_run.dat, DB/db_sbs.dat, and DB/db_sbs.hcal.dat.
const double kEArmThetaDeg = 27.0;
const double kPArmThetaDeg = -18.6;
const double kECalZ = 6.144;
const TVector3 kGEMOriginSBS(-0.147446, -0.00258276, 3.61066);
const TVector3 kGEMAnglesDeg(-0.0290593, -0.342306, 0.309045);
const TVector3 kHCalOriginGEM(0.192, -0.0215, 6.696);

struct Axes {
  TVector3 x;
  TVector3 y;
  TVector3 z;
};

struct HypothesisRank {
  int hypothesis;
  int vertex;
  double chi2NDF;
};

Axes ArmAxes(double thetaDeg)
{
  const double theta = thetaDeg * TMath::DegToRad();
  Axes axes;
  axes.z.SetXYZ(std::sin(theta), 0.0, std::cos(theta));
  axes.x.SetXYZ(0.0, -1.0, 0.0); // analyzer transport +x points down
  axes.y = axes.z.Cross(axes.x).Unit();
  return axes;
}

Axes SurveyedGEMAxes()
{
  TRotation rotation;
  rotation.RotateY(kGEMAnglesDeg.Y() * TMath::DegToRad());
  rotation.RotateZ(kGEMAnglesDeg.Z() * TMath::DegToRad());
  rotation.RotateX(kGEMAnglesDeg.X() * TMath::DegToRad());
  Axes axes;
  axes.x.SetXYZ(rotation.XX(), rotation.YX(), rotation.ZX());
  axes.y.SetXYZ(rotation.XY(), rotation.YY(), rotation.ZY());
  axes.z.SetXYZ(rotation.XZ(), rotation.YZ(), rotation.ZZ());
  return axes;
}

TVector3 ToHall(const TVector3 &local, const Axes &axes,
                const TVector3 &origin = TVector3())
{
  return origin + local.X() * axes.x + local.Y() * axes.y +
      local.Z() * axes.z;
}

TVector3 GEMToHall(const TVector3 &gemLocal, const Axes &protonAxes,
                   const Axes &gemAxes)
{
  const TVector3 inSBS = kGEMOriginSBS + gemLocal.X() * gemAxes.x +
      gemLocal.Y() * gemAxes.y + gemLocal.Z() * gemAxes.z;
  return ToHall(inSBS, protonAxes);
}

double SafeValue(const TTreeReaderArray<Double_t> &values, int index)
{
  return index >= 0 && index < static_cast<int>(values.GetSize()) ?
      values[index] : std::numeric_limits<double>::quiet_NaN();
}

void DrawLine(const TVector3 &a, const TVector3 &b, int color, int width = 2,
              int style = 1)
{
  auto *line = new TPolyLine3D(2);
  line->SetPoint(0, a.X(), a.Y(), a.Z());
  line->SetPoint(1, b.X(), b.Y(), b.Z());
  line->SetLineColor(color);
  line->SetLineWidth(width);
  line->SetLineStyle(style);
  line->Draw("same");
}

void DrawPlane(const TVector3 &center, const Axes &axes, double halfX,
               double halfY, int color, int width = 1, int style = 1)
{
  auto *outline = new TPolyLine3D(5);
  const double signs[5][2] = {
      {-1., -1.}, {1., -1.}, {1., 1.}, {-1., 1.}, {-1., -1.}};
  for (int i = 0; i < 5; ++i) {
    const TVector3 p = center + signs[i][0] * halfX * axes.x +
        signs[i][1] * halfY * axes.y;
    outline->SetPoint(i, p.X(), p.Y(), p.Z());
  }
  outline->SetLineColor(color);
  outline->SetLineWidth(width);
  outline->SetLineStyle(style);
  outline->Draw("same");
}

void DrawMarker(const TVector3 &point, int color, int style, double size)
{
  auto *marker = new TPolyMarker3D(1);
  marker->SetPoint(0, point.X(), point.Y(), point.Z());
  marker->SetMarkerColor(color);
  marker->SetMarkerStyle(style);
  marker->SetMarkerSize(size);
  marker->Draw("same");
}

TString Format(const char *format, double value)
{
  return TString::Format(format, value);
}

} // namespace CDetFTROIEventDisplay

// Interactive Hall-coordinate event display for the Run 5711 FTROI/GEM sample.
//
// selectedHypothesis < 0 chooses the hypothesis with the smallest x chi2/NDF.
// maxHypotheses > 1 overlays that many ranked FTROI rays (the selected one is
// always first). The display intentionally uses skeletal detector outlines;
// it is an event-inspection tool, not a GEANT geometry model.
void Display_CDet_FTROI_GEMEvent(
    Long64_t requestedEventNumber,
    int selectedHypothesis = -1,
    int maxHypotheses = 1,
    const char *inputPattern =
        "/Users/brash/CDet_replay/sbs/Rootfiles/FTROI_step5/rootfiles/"
        "gep5_replayed_5711_stream0_2_seg*_firstevent0_nevent100000*.root",
    const char *outputPNG = "")
{
  using namespace CDetFTROIEventDisplay;

  TChain chain("T");
  const int inputFiles = chain.Add(inputPattern);
  if (inputFiles <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[FTROI display] No readable T trees match " << inputPattern
              << std::endl;
    return;
  }

  TTreeReader reader(&chain);
  TTreeReaderValue<Double_t> eventNumber(reader, "g.evnum");
  TTreeReaderValue<Double_t> ecalE(reader, "earm.ecal.e");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<Double_t> hcalE(reader, "sbs.hcal.e");
  TTreeReaderValue<Double_t> hcalX(reader, "sbs.hcal.x");
  TTreeReaderValue<Double_t> hcalY(reader, "sbs.hcal.y");
  TTreeReaderValue<Double_t> gemNTrack(reader, "sbs.gemFT.track.ntrack");
  TTreeReaderValue<Double_t> gemBestTrack(reader, "sbs.gemFT.track.besttrack");
  TTreeReaderArray<Double_t> gemX(reader, "sbs.gemFT.track.x");
  TTreeReaderArray<Double_t> gemY(reader, "sbs.gemFT.track.y");
  TTreeReaderArray<Double_t> gemXP(reader, "sbs.gemFT.track.xp");
  TTreeReaderArray<Double_t> gemYP(reader, "sbs.gemFT.track.yp");
  TTreeReaderArray<Double_t> gemNHits(reader, "sbs.gemFT.track.nhits");
  TTreeReaderArray<Double_t> gemChi2NDF(reader, "sbs.gemFT.track.chi2ndf");
  TTreeReaderValue<Double_t> heepValid(reader, "heep.datavalid");
  TTreeReaderValue<Double_t> heepDPE(reader, "heep.dpe");
  TTreeReaderValue<Double_t> heepDPP(reader, "heep.dpp");
  TTreeReaderValue<Double_t> heepDPhi(reader, "heep.dphi");
  TTreeReaderArray<Double_t> pulseLayer(reader, "earm.cdet.pulse.layer");
  TTreeReaderArray<Double_t> pulseX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> pulseY(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderArray<Double_t> pulseEligible(reader,
                                           "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> hypPulseL1(reader,
                                        "FTROI.cdet.hyp.pulse_index_l1");
  TTreeReaderArray<Double_t> hypPulseL2(reader,
                                        "FTROI.cdet.hyp.pulse_index_l2");
  TTreeReaderArray<Double_t> hypTopology(reader,
                                        "FTROI.cdet.hyp.y_topology");
  TTreeReaderArray<Double_t> hypScore(reader,
                                      "FTROI.cdet.hyp.source_score");
  TTreeReaderArray<Double_t> vertexBin(reader, "FTROI.cdet.vertex.bin");
  TTreeReaderArray<Double_t> vertexHyp(reader,
                                       "FTROI.cdet.vertex.hyp_index");
  TTreeReaderArray<Double_t> vertexZ(reader, "FTROI.cdet.vertex.z");
  TTreeReaderArray<Double_t> vertexXChi2(reader,
                                         "FTROI.cdet.vertex.xchi2");
  TTreeReaderArray<Double_t> vertexXNDF(reader, "FTROI.cdet.vertex.xndf");
  TTreeReaderArray<Double_t> vertexYCompatible(
      reader, "FTROI.cdet.vertex.ycompatible");
  TTreeReaderArray<Double_t> vertexYSeam(
      reader, "FTROI.cdet.vertex.yseam_compatible");

  while (reader.Next()) {
    if (std::lround(*eventNumber) != requestedEventNumber)
      continue;

    const int nHyp = static_cast<int>(hypPulseL1.GetSize());
    if (nHyp == 0) {
      std::cerr << "[FTROI display] Event " << requestedEventNumber
                << " has no FTROI hypotheses" << std::endl;
      return;
    }

    std::vector<HypothesisRank> ranking;
    for (int ih = 0; ih < nHyp; ++ih) {
      HypothesisRank rank = {ih, -1, std::numeric_limits<double>::infinity()};
      for (int iv = 0; iv < static_cast<int>(vertexZ.GetSize()); ++iv) {
        if (std::lround(vertexHyp[iv]) != ih || vertexXNDF[iv] <= 0.0)
          continue;
        const double value = vertexXChi2[iv] / vertexXNDF[iv];
        if (std::isfinite(value) && value < rank.chi2NDF) {
          rank.vertex = iv;
          rank.chi2NDF = value;
        }
      }
      if (rank.vertex >= 0)
        ranking.push_back(rank);
    }
    std::sort(ranking.begin(), ranking.end(),
              [](const HypothesisRank &a, const HypothesisRank &b) {
                return a.chi2NDF < b.chi2NDF;
              });
    if (ranking.empty()) {
      std::cerr << "[FTROI display] Event " << requestedEventNumber
                << " has no valid FTROI vertex fits" << std::endl;
      return;
    }

    if (selectedHypothesis < 0)
      selectedHypothesis = ranking.front().hypothesis;
    auto selectedIt = std::find_if(
        ranking.begin(), ranking.end(), [selectedHypothesis](const HypothesisRank &r) {
          return r.hypothesis == selectedHypothesis;
        });
    if (selectedIt == ranking.end()) {
      std::cerr << "[FTROI display] Hypothesis " << selectedHypothesis
                << " is absent or has no valid fit in event "
                << requestedEventNumber << std::endl;
      return;
    }
    std::rotate(ranking.begin(), selectedIt, selectedIt + 1);
    maxHypotheses = std::max(1, std::min(maxHypotheses,
                                         static_cast<int>(ranking.size())));
    const HypothesisRank selected = ranking.front();

    const Axes eAxes = ArmAxes(kEArmThetaDeg);
    const Axes pAxes = ArmAxes(kPArmThetaDeg);
    const Axes gemAxes = SurveyedGEMAxes();
    const TVector3 ecalPoint = ToHall(TVector3(*ecalX, *ecalY, kECalZ), eAxes);

    gStyle->SetCanvasPreferGL(kTRUE);
    auto *canvas = new TCanvas(
        TString::Format("c_ftroi_event_%lld", requestedEventNumber),
        TString::Format("Run 5711 event %lld", requestedEventNumber),
        1500, 900);
    auto *viewPad = new TPad("ftroi_view", "Hall-coordinate view", 0.0, 0.0,
                             0.78, 1.0);
    auto *textPad = new TPad("ftroi_text", "Event summary", 0.78, 0.0,
                             1.0, 1.0);
    viewPad->SetFillColor(kWhite);
    textPad->SetFillColor(kWhite);
    viewPad->Draw();
    textPad->Draw();
    viewPad->cd();

    auto *frame = new TH3D(
        TString::Format("ftroi_frame_%lld", requestedEventNumber),
        TString::Format("Run 5711 event %lld;Hall x (m);Hall y (m);Hall z (m)",
                        requestedEventNumber),
        // Equal 12 m spans preserve physical angles and perpendicularity in
        // the ROOT 3D projection. The unused Hall-y space is intentional.
        10, -6.0, 6.0, 10, -6.0, 6.0, 10, -0.75, 11.25);
    frame->SetStats(kFALSE);
    frame->Draw();

    // Hall beam line and nominal central axes. The solid black segment is the
    // electron-arm transport z axis, ending at the center of the ECal plane.
    DrawLine(TVector3(0, 0, -0.5), TVector3(0, 0, 10.5), kGray + 2, 2, 2);
    DrawLine(TVector3(), ToHall(TVector3(0, 0, kECalZ), eAxes), kBlack, 4, 1);
    DrawLine(TVector3(), ToHall(TVector3(0, 0, 10.2), pAxes), kOrange + 7,
             1, 3);

    // ECal and the observed ECal cluster.
    DrawPlane(ToHall(TVector3(0, 0, kECalZ), eAxes), eAxes, 1.481, 0.625,
              kBlue + 1, 2);
    DrawMarker(ecalPoint, kBlue + 1, 29, 1.7);

    // CDet layer planes are placed at their observed pulse z values. All
    // pulses are shown faintly; ECal-eligible and selected pulses stand out.
    double layerZ[2] = {0.0, 0.0};
    int layerN[2] = {0, 0};
    for (int ip = 0; ip < static_cast<int>(pulseZ.GetSize()); ++ip) {
      const int layer = static_cast<int>(std::lround(pulseLayer[ip]));
      if (layer >= 0 && layer < 2 && std::isfinite(pulseZ[ip])) {
        layerZ[layer] += pulseZ[ip];
        ++layerN[layer];
      }
      const TVector3 p = ToHall(TVector3(pulseX[ip], pulseY[ip], pulseZ[ip]),
                                eAxes);
      DrawMarker(p, pulseEligible[ip] == 1.0 ? kCyan + 2 : kGray + 1,
                 pulseEligible[ip] == 1.0 ? 20 : 7,
                 pulseEligible[ip] == 1.0 ? 0.45 : 0.25);
    }
    for (int layer = 0; layer < 2; ++layer) {
      if (layerN[layer] > 0) {
        layerZ[layer] /= layerN[layer];
        DrawPlane(ToHall(TVector3(0, 0, layerZ[layer]), eAxes), eAxes,
                  1.481, 0.625, kCyan + 2, 2, layer == 0 ? 1 : 2);
      }
    }

    // Selected pulse pair and one or more ranked FTROI rays.
    const int rayColors[] = {kRed + 1, kMagenta + 1, kViolet + 1,
                             kPink + 7, kOrange + 10};
    for (int ir = 0; ir < maxHypotheses; ++ir) {
      const HypothesisRank &rank = ranking[ir];
      const TVector3 target(0.0, 0.0, vertexZ[rank.vertex]);
      DrawLine(target, ecalPoint, rayColors[ir % 5], ir == 0 ? 4 : 2,
               ir == 0 ? 1 : 2);
    }
    const int selectedPulses[2] = {
        static_cast<int>(std::lround(hypPulseL1[selectedHypothesis])),
        static_cast<int>(std::lround(hypPulseL2[selectedHypothesis]))};
    bool havePairMean = false;
    TVector3 pairMeanLocal;
    for (int i = 0; i < 2; ++i) {
      const int ip = selectedPulses[i];
      if (ip >= 0 && ip < static_cast<int>(pulseZ.GetSize())) {
        DrawMarker(ToHall(TVector3(pulseX[ip], pulseY[ip], pulseZ[ip]), eAxes),
                   i == 0 ? kRed + 1 : kOrange + 7, 29, 2.0);
      }
    }
    if (selectedPulses[0] >= 0 && selectedPulses[1] >= 0 &&
        selectedPulses[0] < static_cast<int>(pulseZ.GetSize()) &&
        selectedPulses[1] < static_cast<int>(pulseZ.GetSize())) {
      pairMeanLocal = 0.5 *
          (TVector3(pulseX[selectedPulses[0]], pulseY[selectedPulses[0]],
                    pulseZ[selectedPulses[0]]) +
           TVector3(pulseX[selectedPulses[1]], pulseY[selectedPulses[1]],
                    pulseZ[selectedPulses[1]]));
      havePairMean = true;
      DrawMarker(ToHall(pairMeanLocal, eAxes), kMagenta + 2, 34, 2.2);
    }

    // Target scan points for the selected hypothesis: green is y-compatible,
    // red is incompatible, and the best-x point is enlarged.
    for (int iv = 0; iv < static_cast<int>(vertexZ.GetSize()); ++iv) {
      if (std::lround(vertexHyp[iv]) != selectedHypothesis)
        continue;
      DrawMarker(TVector3(0, 0, vertexZ[iv]),
                 vertexYCompatible[iv] == 1.0 ? kGreen + 2 : kRed + 1,
                 iv == selected.vertex ? 29 : 20,
                 iv == selected.vertex ? 1.5 : 0.7);
    }

    // GEM FT skeletal planes, best track, HCAL plane, and HCAL cluster.
    const double gemPlaneZ[] = {0.000, 0.075, 0.202, 0.341,
                                0.484, 0.624, 0.765, 0.890};
    for (double z : gemPlaneZ) {
      const TVector3 center = GEMToHall(TVector3(0, 0, z), pAxes, gemAxes);
      Axes hallGEMAxes = {
          ToHall(gemAxes.x, pAxes) - ToHall(TVector3(), pAxes),
          ToHall(gemAxes.y, pAxes) - ToHall(TVector3(), pAxes),
          ToHall(gemAxes.z, pAxes) - ToHall(TVector3(), pAxes)};
      DrawPlane(center, hallGEMAxes, z < 0.7 ? 0.75 : 1.03,
                z < 0.7 ? 0.20 : 0.31, kOrange + 7, 1);
    }
    const int nTrack = static_cast<int>(std::lround(*gemNTrack));
    int bestTrack = static_cast<int>(std::lround(*gemBestTrack));
    if (bestTrack < 0 || bestTrack >= nTrack)
      bestTrack = nTrack > 0 ? 0 : -1;
    if (bestTrack >= 0) {
      const double x = SafeValue(gemX, bestTrack);
      const double y = SafeValue(gemY, bestTrack);
      const double xp = SafeValue(gemXP, bestTrack);
      const double yp = SafeValue(gemYP, bestTrack);
      DrawLine(GEMToHall(TVector3(x - 0.10 * xp, y - 0.10 * yp, -0.10),
                         pAxes, gemAxes),
               GEMToHall(TVector3(x + 1.05 * xp, y + 1.05 * yp, 1.05),
                         pAxes, gemAxes),
               kOrange + 1, 4);
    }
    const TVector3 hcalCenter = GEMToHall(kHCalOriginGEM, pAxes, gemAxes);
    Axes hallGEMAxes = {
        ToHall(gemAxes.x, pAxes) - ToHall(TVector3(), pAxes),
        ToHall(gemAxes.y, pAxes) - ToHall(TVector3(), pAxes),
        ToHall(gemAxes.z, pAxes) - ToHall(TVector3(), pAxes)};
    DrawPlane(hcalCenter, hallGEMAxes, 0.93, 0.47, kGreen + 3, 2);
    if (std::isfinite(*hcalX) && std::isfinite(*hcalY))
      DrawMarker(GEMToHall(kHCalOriginGEM + TVector3(*hcalX, *hcalY, 0),
                           pAxes, gemAxes),
                 kGreen + 3, 33, 1.6);

    viewPad->Modified();
    viewPad->Update();
    if (viewPad->GetView()) {
      int viewError = 0;
      viewPad->GetView()->SetView(-55.0, 72.0, 0.0, viewError);
      viewPad->Modified();
    }

    textPad->cd();
    auto *summary = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
    summary->SetFillColor(kWhite);
    summary->SetBorderSize(0);
    summary->SetTextAlign(12);
    summary->SetTextFont(42);
    summary->SetTextSize(0.031);
    summary->AddText(TString::Format("Run 5711, event %lld", requestedEventNumber));
    summary->AddText(TString::Format("file: %s", chain.GetCurrentFile() ?
        gSystem->BaseName(chain.GetCurrentFile()->GetName()) : "unknown"));
    summary->AddText(" ");
    summary->AddText("Electron arm (blue/cyan/red)");
    summary->AddText(TString::Format("ECal: E=%.3f GeV", *ecalE));
    summary->AddText(TString::Format("       x=%.3f, y=%.3f m", *ecalX, *ecalY));
    summary->AddText(TString::Format("CDet pulses: %d",
        static_cast<int>(pulseZ.GetSize())));
    summary->AddText(TString::Format("selected FTROI hypothesis: %d",
                                     selectedHypothesis));
    summary->AddText(TString::Format("pulse indices: %d / %d",
                                     selectedPulses[0], selectedPulses[1]));
    if (havePairMean)
      summary->AddText(TString::Format("pair mean x/y/z: %.4f / %.4f / %.4f m",
          pairMeanLocal.X(), pairMeanLocal.Y(), pairMeanLocal.Z()));
    summary->AddText(TString::Format("topology: %s",
        std::lround(hypTopology[selectedHypothesis]) == 1 ? "seam" :
        (std::lround(hypTopology[selectedHypothesis]) == 0 ? "same-side" :
                                                            "single")));
    summary->AddText(TString::Format("source score: %.4f",
                                     hypScore[selectedHypothesis]));
    summary->AddText(TString::Format("best z: %.4f m",
                                     vertexZ[selected.vertex]));
    summary->AddText(TString::Format("best x #chi^{2}/NDF: %.6g",
                                     selected.chi2NDF));
    summary->AddText(TString::Format("y-compatible / seam: %d / %d",
        static_cast<int>(std::lround(vertexYCompatible[selected.vertex])),
        static_cast<int>(std::lround(vertexYSeam[selected.vertex]))));
    summary->AddText(" ");
    summary->AddText("Proton arm (orange/green)");
    summary->AddText(TString::Format("GEM tracks: %d, best: %d", nTrack,
                                     bestTrack));
    if (bestTrack >= 0) {
      summary->AddText(TString::Format("GEM hits / #chi^{2}/NDF: %.0f / %.3f",
          SafeValue(gemNHits, bestTrack), SafeValue(gemChi2NDF, bestTrack)));
      summary->AddText(TString::Format("x,y = %.4f, %.4f m",
          SafeValue(gemX, bestTrack), SafeValue(gemY, bestTrack)));
      summary->AddText(TString::Format("x',y' = %.4f, %.4f",
          SafeValue(gemXP, bestTrack), SafeValue(gemYP, bestTrack)));
    }
    summary->AddText(TString::Format("HCal: E=%.3f GeV", *hcalE));
    summary->AddText(" ");
    summary->AddText(TString::Format("HEEP valid: %d",
                                     static_cast<int>(std::lround(*heepValid))));
    summary->AddText(TString::Format("dpe=%.4f  dpp=%.4f", *heepDPE, *heepDPP));
    summary->AddText(TString::Format("dphi=%.4f rad", *heepDPhi));
    summary->AddText(" ");
    summary->AddText("Dashed gray: Hall beam axis");
    summary->AddText("Solid black: ECal/CDet z axis");
    summary->AddText("Dashed orange: proton central ray");
    summary->AddText("Magenta cross: selected CDet pair mean");
    summary->AddText("Target dots: green compatible, red not");
    summary->AddText("Arms share Hall coordinates; slopes do not.");
    summary->Draw();

    canvas->cd();
    canvas->Modified();
    canvas->Update();
    if (outputPNG && TString(outputPNG).Length() > 0)
      canvas->SaveAs(outputPNG);

    std::cout << std::setprecision(10)
              << "[FTROI display] event=" << requestedEventNumber
              << " selected_hypothesis=" << selectedHypothesis
              << " best_vertex_bin=" << std::lround(vertexBin[selected.vertex])
              << " best_z=" << vertexZ[selected.vertex]
              << " xchi2ndf=" << selected.chi2NDF
              << " overlays=" << maxHypotheses << std::endl;
    std::cout << "[FTROI display] Geometry: GEP-3 Hall coordinates, "
                 "skeletal detector outlines. Drag the 3D pad to rotate."
              << std::endl;
    return;
  }

  std::cerr << "[FTROI display] Event " << requestedEventNumber
            << " was not found in the input chain" << std::endl;
}
