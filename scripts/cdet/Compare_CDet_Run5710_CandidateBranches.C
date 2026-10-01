#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "CDetRunDataset.h"

namespace {

struct ReconstructedPair {
  int pulse1;
  int pulse2;
  double timeMean;
  double dt;
  double dx;
  double dy;
  double score;
  double ecalResidual;
  double trajectoryResidual;
};

void StyleOldNew(TH1D &oldHist, TH1D &newHist) {
  oldHist.SetLineColor(kBlue + 1);
  oldHist.SetLineWidth(3);
  newHist.SetLineColor(kRed + 1);
  newHist.SetLineStyle(2);
  newHist.SetLineWidth(3);
  oldHist.SetStats(false);
  newHist.SetStats(false);
}

void DrawOverlay(TH1D &oldHist, TH1D &newHist, TLegend &legend) {
  StyleOldNew(oldHist, newHist);
  const double maximum = std::max(oldHist.GetMaximum(), newHist.GetMaximum());
  oldHist.SetMaximum(maximum > 0.0 ? 1.12 * maximum : 1.0);
  oldHist.Draw("hist");
  newHist.Draw("hist same");
  legend.DrawClone();
}

} // namespace

// Standalone Run 5710 parity study. "Old" reconstructs the pair hypotheses
// ex post facto from pulse.*. "New" reads only SBSCDet pair_candidate.*.
void Compare_CDet_Run5710_CandidateBranches(
    const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_run5710_candidate_branch_comparison") {
  constexpr int runNumber = 5710;
  constexpr double dtCenterNs = 0.0;
  constexpr double dtMaxNs = 15.0;
  constexpr double dxMaxM = 0.15;
  constexpr double dyMaxM = 0.08;
  constexpr double timeScaleNs = 1.0;
  constexpr double xScaleM = 0.01;
  constexpr double tolerance = 1.0e-10;

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[candidate comparison] inputDirectory or OUT_DIR is required.\n";
    return;
  }

  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[candidate comparison] No Run 5710 events found in "
              << input << ".\n";
    return;
  }

  const char *required[] = {
      "earm.cdet.pulse.pmtnum", "earm.cdet.pulse.tdc_le_corr",
      "earm.cdet.pulse.ecal_residual", "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass",
      "earm.cdet.pulse.ecal_eligible", "earm.cdet.pulse.spatial_pass",
      "earm.cdet.pulse.x_corr", "earm.cdet.pulse.y",
      "earm.cdet.pulse.ecal_x_proj",
      "earm.cdet.pair_candidate.pulse_index_l1",
      "earm.cdet.pair_candidate.pulse_index_l2",
      "earm.cdet.pair_candidate.time_mean",
      "earm.cdet.pair_candidate.dt", "earm.cdet.pair_candidate.dx",
      "earm.cdet.pair_candidate.dy", "earm.cdet.pair_candidate.score",
      "earm.cdet.pair_candidate.ecal_residual",
      "earm.cdet.pair_candidate.trajectory_residual"};
  chain.LoadTree(0);
  for (const char *branch : required) {
    if (!chain.GetBranch(branch)) {
      std::cerr << "[candidate comparison] Missing branch " << branch << ".\n";
      return;
    }
  }

  gSystem->mkdir(outputDirectory, true);
  gStyle->SetOptStat(0);

  TH1D oldMultiplicity("hOldPairCandidateMultiplicity",
      "Pair hypotheses per event;hypotheses;events", 31, -0.5, 30.5);
  TH1D newMultiplicity("hNewPairCandidateMultiplicity",
      "Pair hypotheses per event;hypotheses;events", 31, -0.5, 30.5);
  TH1D oldMeanTime("hOldPairCandidateMeanTime",
      "Pair mean corrected LE;mean corrected LE (ns);hypotheses", 120, 0, 60);
  TH1D newMeanTime("hNewPairCandidateMeanTime",
      "Pair mean corrected LE;mean corrected LE (ns);hypotheses", 120, 0, 60);
  TH1D oldDT("hOldPairCandidateDT",
      "Inter-layer timing;#Deltat = t_{L2}-t_{L1} (ns);hypotheses", 120, -15, 15);
  TH1D newDT("hNewPairCandidateDT",
      "Inter-layer timing;#Deltat = t_{L2}-t_{L1} (ns);hypotheses", 120, -15, 15);
  TH1D oldDX("hOldPairCandidateDX",
      "Inter-layer aligned x;#Deltax = x_{L2}-x_{L1} (m);hypotheses", 120, -0.15, 0.15);
  TH1D newDX("hNewPairCandidateDX",
      "Inter-layer aligned x;#Deltax = x_{L2}-x_{L1} (m);hypotheses", 120, -0.15, 0.15);
  TH1D oldDY("hOldPairCandidateDY",
      "Inter-layer y;#Deltay = y_{L2}-y_{L1} (m);hypotheses", 100, -0.10, 0.10);
  TH1D newDY("hNewPairCandidateDY",
      "Inter-layer y;#Deltay = y_{L2}-y_{L1} (m);hypotheses", 100, -0.10, 0.10);
  TH1D oldECalResidual("hOldPairCandidateECalResidual",
      "ECal minus pair mean time;t_{ECal}-<t_{CDet}> (ns);hypotheses", 160, -40, 40);
  TH1D newECalResidual("hNewPairCandidateECalResidual",
      "ECal minus pair mean time;t_{ECal}-<t_{CDet}> (ns);hypotheses", 160, -40, 40);
  TH1D oldTrajectoryResidual("hOldPairCandidateTrajectoryResidual",
      "ECal-trajectory residual;#Deltax_{pair}-#Deltax_{ECal proj.} (m);hypotheses",
      160, -0.20, 0.20);
  TH1D newTrajectoryResidual("hNewPairCandidateTrajectoryResidual",
      "ECal-trajectory residual;#Deltax_{pair}-#Deltax_{ECal proj.} (m);hypotheses",
      160, -0.20, 0.20);
  TH1D oldScore("hOldPairCandidateScore",
      "CDet-only pair score;(#Deltat/1 ns)^{2}+(#Deltax/0.01 m)^{2};hypotheses",
      160, 0, 320);
  TH1D newScore("hNewPairCandidateScore",
      "CDet-only pair score;(#Deltat/1 ns)^{2}+(#Deltax/0.01 m)^{2};hypotheses",
      160, 0, 320);
  TH2D oldTrajectoryVsTiming("hOldPairTrajectoryVsTiming",
      "Ex-post-facto pulse reconstruction;trajectory residual (m);"
      "t_{ECal}-<t_{CDet}> (ns)", 120, -0.15, 0.15, 120, -40, 20);
  TH2D newTrajectoryVsTiming("hNewPairTrajectoryVsTiming",
      "SBSCDet pair_candidate.*;trajectory residual (m);"
      "t_{ECal}-<t_{CDet}> (ns)", 120, -0.15, 0.15, 120, -40, 20);

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pmt(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> pulseECalResidual(
      reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> calibValid(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> ecalEligible(
      reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatialPass(reader, "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<Double_t> x(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> y(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<Double_t> projectedX(reader, "earm.cdet.pulse.ecal_x_proj");
  TTreeReaderArray<Double_t> newPulse1(
      reader, "earm.cdet.pair_candidate.pulse_index_l1");
  TTreeReaderArray<Double_t> newPulse2(
      reader, "earm.cdet.pair_candidate.pulse_index_l2");
  TTreeReaderArray<Double_t> newTimeMean(
      reader, "earm.cdet.pair_candidate.time_mean");
  TTreeReaderArray<Double_t> newPairDT(reader, "earm.cdet.pair_candidate.dt");
  TTreeReaderArray<Double_t> newPairDX(reader, "earm.cdet.pair_candidate.dx");
  TTreeReaderArray<Double_t> newPairDY(reader, "earm.cdet.pair_candidate.dy");
  TTreeReaderArray<Double_t> newPairScore(reader, "earm.cdet.pair_candidate.score");
  TTreeReaderArray<Double_t> newPairECalResidual(
      reader, "earm.cdet.pair_candidate.ecal_residual");
  TTreeReaderArray<Double_t> newPairTrajectoryResidual(
      reader, "earm.cdet.pair_candidate.trajectory_residual");

  Long64_t events = 0;
  Long64_t malformedEvents = 0;
  Long64_t multiplicityMismatchEvents = 0;
  Long64_t identityMismatchEvents = 0;
  Long64_t valueMismatches = 0;
  Long64_t oldPairs = 0;
  Long64_t newPairs = 0;

  while (reader.Next()) {
    ++events;
    const size_t nPulses = std::min(
        {pmt.GetSize(), le.GetSize(), pulseECalResidual.GetSize(),
         calibValid.GetSize(), broadQuality.GetSize(), ecalEligible.GetSize(),
         spatialPass.GetSize(), x.GetSize(), y.GetSize(), projectedX.GetSize()});
    if (nPulses != pmt.GetSize() || nPulses != le.GetSize() ||
        nPulses != pulseECalResidual.GetSize() ||
        nPulses != calibValid.GetSize() || nPulses != broadQuality.GetSize() ||
        nPulses != ecalEligible.GetSize() || nPulses != spatialPass.GetSize() ||
        nPulses != x.GetSize() || nPulses != y.GetSize() ||
        nPulses != projectedX.GetSize())
      ++malformedEvents;

    std::vector<int> layer1;
    std::vector<int> layer2;
    for (size_t i = 0; i < nPulses; ++i) {
      if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
            ecalEligible[i] > 0.5 && spatialPass[i] > 0.5) ||
          !std::isfinite(pmt[i]) || !std::isfinite(le[i]) ||
          !std::isfinite(x[i]) || !std::isfinite(y[i]) ||
          !std::isfinite(projectedX[i]) ||
          !std::isfinite(pulseECalResidual[i]))
        continue;
      const int pixel = int(std::lround(pmt[i]));
      if (pixel >= 0 && pixel < 1344)
        layer1.push_back(int(i));
      else if (pixel >= 1344 && pixel < 2688)
        layer2.push_back(int(i));
    }

    std::vector<ReconstructedPair> reconstructed;
    for (int pulse1 : layer1) {
      for (int pulse2 : layer2) {
        const double dt = le[pulse2] - le[pulse1];
        const double dx = x[pulse2] - x[pulse1];
        const double dy = y[pulse2] - y[pulse1];
        if (std::fabs(dt - dtCenterNs) > dtMaxNs ||
            std::fabs(dx) > dxMaxM || std::fabs(dy) > dyMaxM)
          continue;
        const double timePull = (dt - dtCenterNs) / timeScaleNs;
        const double xPull = dx / xScaleM;
        reconstructed.push_back({
            pulse1, pulse2, 0.5 * (le[pulse1] + le[pulse2]), dt, dx, dy,
            timePull * timePull + xPull * xPull,
            0.5 * (pulseECalResidual[pulse1] + pulseECalResidual[pulse2]),
            dx - (projectedX[pulse2] - projectedX[pulse1])});
      }
    }
    std::stable_sort(reconstructed.begin(), reconstructed.end(),
        [](const ReconstructedPair &left, const ReconstructedPair &right) {
          return left.score < right.score;
        });

    const size_t nNew = std::min(
        {newPulse1.GetSize(), newPulse2.GetSize(), newTimeMean.GetSize(),
         newPairDT.GetSize(), newPairDX.GetSize(), newPairDY.GetSize(),
         newPairScore.GetSize(), newPairECalResidual.GetSize(),
         newPairTrajectoryResidual.GetSize()});
    oldMultiplicity.Fill(reconstructed.size());
    newMultiplicity.Fill(nNew);
    oldPairs += reconstructed.size();
    newPairs += nNew;
    if (reconstructed.size() != nNew)
      ++multiplicityMismatchEvents;

    std::map<std::pair<int, int>, ReconstructedPair> oldByIdentity;
    for (const ReconstructedPair &pair : reconstructed) {
      oldByIdentity[{pair.pulse1, pair.pulse2}] = pair;
      oldMeanTime.Fill(pair.timeMean);
      oldDT.Fill(pair.dt);
      oldDX.Fill(pair.dx);
      oldDY.Fill(pair.dy);
      oldScore.Fill(pair.score);
      oldECalResidual.Fill(pair.ecalResidual);
      oldTrajectoryResidual.Fill(pair.trajectoryResidual);
      oldTrajectoryVsTiming.Fill(pair.trajectoryResidual, pair.ecalResidual);
    }

    bool eventIdentityMismatch = false;
    for (size_t i = 0; i < nNew; ++i) {
      const int pulse1 = int(std::lround(newPulse1[i]));
      const int pulse2 = int(std::lround(newPulse2[i]));
      newMeanTime.Fill(newTimeMean[i]);
      newDT.Fill(newPairDT[i]);
      newDX.Fill(newPairDX[i]);
      newDY.Fill(newPairDY[i]);
      newScore.Fill(newPairScore[i]);
      newECalResidual.Fill(newPairECalResidual[i]);
      newTrajectoryResidual.Fill(newPairTrajectoryResidual[i]);
      newTrajectoryVsTiming.Fill(newPairTrajectoryResidual[i],
                                 newPairECalResidual[i]);
      const auto found = oldByIdentity.find({pulse1, pulse2});
      if (found == oldByIdentity.end()) {
        eventIdentityMismatch = true;
        continue;
      }
      const ReconstructedPair &old = found->second;
      const double values[][2] = {
          {old.timeMean, newTimeMean[i]}, {old.dt, newPairDT[i]},
          {old.dx, newPairDX[i]}, {old.dy, newPairDY[i]},
          {old.score, newPairScore[i]},
          {old.ecalResidual, newPairECalResidual[i]},
          {old.trajectoryResidual, newPairTrajectoryResidual[i]}};
      for (const auto &value : values) {
        if (!std::isfinite(value[0]) || !std::isfinite(value[1]) ||
            std::fabs(value[0] - value[1]) > tolerance)
          ++valueMismatches;
      }
      oldByIdentity.erase(found);
    }
    if (!oldByIdentity.empty())
      eventIdentityMismatch = true;
    if (eventIdentityMismatch)
      ++identityMismatchEvents;
  }

  TLegend legend(0.54, 0.73, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  StyleOldNew(oldDT, newDT);
  legend.AddEntry(&oldDT, "Old: reconstructed from pulse.*", "l");
  legend.AddEntry(&newDT, "New: SBSCDet pair_candidate.*", "l");

  TCanvas overview("cCDetCandidateBranchComparison",
                   "Run 5710 candidate branch comparison", 1800, 1100);
  overview.Divide(3, 2);
  overview.cd(1); DrawOverlay(oldMultiplicity, newMultiplicity, legend);
  overview.cd(2); DrawOverlay(oldDT, newDT, legend);
  overview.cd(3); DrawOverlay(oldDX, newDX, legend);
  overview.cd(4); DrawOverlay(oldDY, newDY, legend);
  overview.cd(5); DrawOverlay(oldECalResidual, newECalResidual, legend);
  overview.cd(6); DrawOverlay(oldTrajectoryResidual, newTrajectoryResidual, legend);
  overview.SaveAs(Form("%s/CDet_Run5710_CandidateBranchComparison.pdf",
                       outputDirectory));
  overview.SaveAs(Form("%s/CDet_Run5710_CandidateBranchComparison.png",
                       outputDirectory));

  TCanvas detail("cCDetCandidateBranchComparisonDetail",
                 "Run 5710 candidate branch comparison details", 1500, 1000);
  detail.Divide(2, 2);
  detail.cd(1); DrawOverlay(oldMeanTime, newMeanTime, legend);
  detail.cd(2); DrawOverlay(oldScore, newScore, legend);
  detail.cd(3); oldTrajectoryVsTiming.Draw("colz");
  detail.cd(4); newTrajectoryVsTiming.Draw("colz");
  detail.SaveAs(Form("%s/CDet_Run5710_CandidateBranchComparison_Detail.pdf",
                     outputDirectory));
  detail.SaveAs(Form("%s/CDet_Run5710_CandidateBranchComparison_Detail.png",
                     outputDirectory));

  TFile output(Form("%s/CDet_Run5710_CandidateBranchComparison.root",
                    outputDirectory), "RECREATE");
  oldMultiplicity.Write(); newMultiplicity.Write();
  oldMeanTime.Write(); newMeanTime.Write();
  oldDT.Write(); newDT.Write(); oldDX.Write(); newDX.Write();
  oldDY.Write(); newDY.Write(); oldScore.Write(); newScore.Write();
  oldECalResidual.Write(); newECalResidual.Write();
  oldTrajectoryResidual.Write(); newTrajectoryResidual.Write();
  oldTrajectoryVsTiming.Write(); newTrajectoryVsTiming.Write();
  output.Close();

  std::ofstream summary(
      Form("%s/CDet_Run5710_CandidateBranchComparison.txt", outputDirectory));
  summary << "Run 5710 SBSCDet candidate-branch parity study\n"
          << "files " << filesAdded << "\n"
          << "events " << events << "\n"
          << "old_hypotheses " << oldPairs << "\n"
          << "new_hypotheses " << newPairs << "\n"
          << "malformed_pulse_array_events " << malformedEvents << "\n"
          << "multiplicity_mismatch_events " << multiplicityMismatchEvents
          << "\nidentity_mismatch_events " << identityMismatchEvents
          << "\nvalue_mismatches " << valueMismatches
          << "\nvalue_tolerance " << tolerance << "\n";
  summary.close();

  std::cout << "\n[candidate comparison] files/events: " << filesAdded << "/"
            << events << "\n"
            << "[candidate comparison] old/new hypotheses: " << oldPairs
            << "/" << newPairs << "\n"
            << "[candidate comparison] malformed pulse-array events: "
            << malformedEvents << "\n"
            << "[candidate comparison] multiplicity-mismatch events: "
            << multiplicityMismatchEvents << "\n"
            << "[candidate comparison] identity-mismatch events: "
            << identityMismatchEvents << "\n"
            << "[candidate comparison] value mismatches (tolerance "
            << tolerance << "): " << valueMismatches << "\n"
            << "[candidate comparison] output: " << outputDirectory
            << std::endl;
}
