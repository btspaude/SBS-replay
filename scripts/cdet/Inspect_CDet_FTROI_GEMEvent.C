#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

double SafeValue(const TTreeReaderArray<Double_t> &values, int index)
{
  return index >= 0 && index < static_cast<int>(values.GetSize()) ?
      values[index] : std::numeric_limits<double>::quiet_NaN();
}

const char *TopologyName(int topology)
{
  if (topology == 0)
    return "same-side pair";
  if (topology == 1)
    return "opposite-side seam pair";
  if (topology == -1)
    return "single layer";
  return "unknown";
}

} // namespace

// Print a detailed, one-event audit for the Run 5711 GEM-track catalogue.
// sbs.gemFT is in the proton arm and FTROI/CDet is in the electron arm, so
// their slopes are not directly subtracted or ranked here.
void Inspect_CDet_FTROI_GEMEvent(
    Long64_t requestedEventNumber,
    const char *inputPattern =
        "/Users/brash/CDet_replay/sbs/Rootfiles/FTROI_step5/rootfiles/"
        "gep5_replayed_5711_stream0_2_seg*_firstevent0_nevent100000*.root")
{
  TChain chain("T");
  const int inputFiles = chain.Add(inputPattern);
  if (inputFiles <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[CDet FTROI event] No readable T trees match "
              << inputPattern << std::endl;
    return;
  }

  TTreeReader reader(&chain);
  TTreeReaderValue<Double_t> eventNumber(reader, "g.evnum");
  TTreeReaderValue<Double_t> gemNTrack(reader, "sbs.gemFT.track.ntrack");
  TTreeReaderValue<Double_t> gemBestTrack(reader, "sbs.gemFT.track.besttrack");
  TTreeReaderArray<Double_t> gemNHits(reader, "sbs.gemFT.track.nhits");
  TTreeReaderArray<Double_t> gemNGoodHits(reader, "sbs.gemFT.track.ngoodhits");
  TTreeReaderArray<Double_t> gemChi2NDF(reader, "sbs.gemFT.track.chi2ndf");
  TTreeReaderArray<Double_t> gemX(reader, "sbs.gemFT.track.x");
  TTreeReaderArray<Double_t> gemY(reader, "sbs.gemFT.track.y");
  TTreeReaderArray<Double_t> gemXP(reader, "sbs.gemFT.track.xp");
  TTreeReaderArray<Double_t> gemYP(reader, "sbs.gemFT.track.yp");
  TTreeReaderArray<Double_t> gemT0(reader, "sbs.gemFT.track.t0");
  TTreeReaderValue<Double_t> heepDataValid(reader, "heep.datavalid");
  TTreeReaderValue<Double_t> heepDPE(reader, "heep.dpe");
  TTreeReaderValue<Double_t> heepDPP(reader, "heep.dpp");
  TTreeReaderValue<Double_t> heepDPhi(reader, "heep.dphi");
  TTreeReaderValue<Double_t> heepAcoplanarity(reader, "heep.acoplanarity");
  TTreeReaderValue<Double_t> heepDXECal(reader, "heep.dxECAL");
  TTreeReaderValue<Double_t> heepDYECal(reader, "heep.dyECAL");
  TTreeReaderValue<Double_t> heepDtADC(reader, "heep.dt_ADC");
  TTreeReaderValue<Double_t> ecalE(reader, "earm.ecal.e");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<Double_t> timingStatus(reader, "FTROI.cdet.timing_status");
  TTreeReaderValue<Double_t> roiStatus(reader, "FTROI.cdet.roi_status");
  TTreeReaderValue<Double_t> nPulse(reader, "FTROI.cdet.npulse");
  TTreeReaderValue<Double_t> nPairCandidate(reader, "FTROI.cdet.npair_candidate");
  TTreeReaderValue<Double_t> nSingleCandidate(reader, "FTROI.cdet.nsingle_candidate");
  TTreeReaderArray<Double_t> hypSourceType(reader, "FTROI.cdet.hyp.source_type");
  TTreeReaderArray<Double_t> hypSourceIndex(reader, "FTROI.cdet.hyp.source_index");
  TTreeReaderArray<Double_t> hypPulseL1(reader, "FTROI.cdet.hyp.pulse_index_l1");
  TTreeReaderArray<Double_t> hypPulseL2(reader, "FTROI.cdet.hyp.pulse_index_l2");
  TTreeReaderArray<Double_t> hypYTopology(reader, "FTROI.cdet.hyp.y_topology");
  TTreeReaderArray<Double_t> hypSourceScore(reader, "FTROI.cdet.hyp.source_score");
  TTreeReaderArray<Double_t> vertexBin(reader, "FTROI.cdet.vertex.bin");
  TTreeReaderArray<Double_t> vertexHyp(reader, "FTROI.cdet.vertex.hyp_index");
  TTreeReaderArray<Double_t> vertexZ(reader, "FTROI.cdet.vertex.z");
  TTreeReaderArray<Double_t> vertexXChi2(reader, "FTROI.cdet.vertex.xchi2");
  TTreeReaderArray<Double_t> vertexXNDF(reader, "FTROI.cdet.vertex.xndf");
  TTreeReaderArray<Double_t> vertexYCompatible(reader, "FTROI.cdet.vertex.ycompatible");
  TTreeReaderArray<Double_t> vertexYSeamCompatible(
      reader, "FTROI.cdet.vertex.yseam_compatible");
  TTreeReaderArray<Double_t> vertexTheta(reader, "FTROI.cdet.vertex.theta_global");
  TTreeReaderArray<Double_t> vertexPhi(reader, "FTROI.cdet.vertex.phi_global");

  while (reader.Next()) {
    if (std::lround(*eventNumber) != requestedEventNumber)
      continue;

    const Long64_t globalEntry = reader.GetCurrentEntry();
    const Long64_t localEntry = chain.GetTree() ?
        chain.GetTree()->GetReadEntry() : -1;
    const TString sourceFile = chain.GetCurrentFile() ?
        chain.GetCurrentFile()->GetName() : "";
    const int nTrack = static_cast<int>(std::lround(*gemNTrack));
    const int nHyp = static_cast<int>(hypSourceType.GetSize());

    std::cout << std::setprecision(10);
    std::cout << "[CDet FTROI event] event=" << requestedEventNumber
              << " global_entry=" << globalEntry
              << " local_entry=" << localEntry
              << " tree=" << chain.GetTreeNumber() << std::endl;
    std::cout << "[CDet FTROI event] file=" << gSystem->BaseName(sourceFile)
              << std::endl;
    std::cout << "[CDet FTROI event] GEM tracks=" << nTrack
              << " best_index=" << std::lround(*gemBestTrack) << std::endl;
    for (int it = 0; it < nTrack; ++it) {
      std::cout << "  GEM[" << it << "] nhits=" << SafeValue(gemNHits, it)
                << " ngoodhits=" << SafeValue(gemNGoodHits, it)
                << " chi2ndf=" << SafeValue(gemChi2NDF, it)
                << " x=" << SafeValue(gemX, it)
                << " y=" << SafeValue(gemY, it)
                << " xp=" << SafeValue(gemXP, it)
                << " yp=" << SafeValue(gemYP, it)
                << " t0=" << SafeValue(gemT0, it) << std::endl;
    }

    std::cout << "[CDet FTROI event] ECal E/x/y=" << *ecalE << "/"
              << *ecalX << "/" << *ecalY << std::endl;
    std::cout << "[CDet FTROI event] HEEP valid=" << std::lround(*heepDataValid)
              << " dpe=" << *heepDPE << " dpp=" << *heepDPP
              << " dphi=" << *heepDPhi
              << " acoplanarity=" << *heepAcoplanarity
              << " dxECAL=" << *heepDXECal << " dyECAL=" << *heepDYECal
              << " dt_ADC=" << *heepDtADC << std::endl;
    std::cout << "[CDet FTROI event] CDet timing/ROI="
              << std::lround(*timingStatus) << "/" << std::lround(*roiStatus)
              << " pulses=" << std::lround(*nPulse)
              << " pair_candidates=" << std::lround(*nPairCandidate)
              << " single_candidates=" << std::lround(*nSingleCandidate)
              << " hypotheses=" << nHyp << std::endl;

    std::vector<double> bestChi2NDF(
        nHyp, std::numeric_limits<double>::infinity());
    std::vector<int> bestVertex(nHyp, -1);
    std::vector<int> compatible(nHyp, 0);
    std::vector<int> associations(nHyp, 0);
    for (int iv = 0; iv < static_cast<int>(vertexZ.GetSize()); ++iv) {
      const int ih = static_cast<int>(std::lround(vertexHyp[iv]));
      if (ih < 0 || ih >= nHyp)
        continue;
      ++associations[ih];
      if (vertexYCompatible[iv] == 1.0)
        ++compatible[ih];
      const double ndf = vertexXNDF[iv];
      const double chi2NDF = ndf > 0.0 ? vertexXChi2[iv] / ndf :
          std::numeric_limits<double>::infinity();
      if (std::isfinite(chi2NDF) && chi2NDF < bestChi2NDF[ih]) {
        bestChi2NDF[ih] = chi2NDF;
        bestVertex[ih] = iv;
      }
    }

    for (int ih = 0; ih < nHyp; ++ih) {
      const int iv = bestVertex[ih];
      const int topology = static_cast<int>(std::lround(hypYTopology[ih]));
      std::cout << "  FTROI[" << ih << "] source_type="
                << std::lround(hypSourceType[ih])
                << " source_index=" << std::lround(hypSourceIndex[ih])
                << " pulse_l1/l2=" << std::lround(hypPulseL1[ih]) << "/"
                << std::lround(hypPulseL2[ih])
                << " topology=" << topology << " (" << TopologyName(topology)
                << ") source_score=" << hypSourceScore[ih]
                << " ycompatible=" << compatible[ih] << "/"
                << associations[ih];
      if (iv >= 0) {
        std::cout << " best_bin=" << std::lround(vertexBin[iv])
                  << " best_z=" << vertexZ[iv]
                  << " best_xchi2ndf=" << bestChi2NDF[ih]
                  << " best_ycompatible=" << std::lround(vertexYCompatible[iv])
                  << " best_yseam=" << std::lround(vertexYSeamCompatible[iv])
                  << " theta_lab=" << vertexTheta[iv]
                  << " phi_lab=" << vertexPhi[iv];
      }
      std::cout << std::endl;
    }
    std::cout << "[CDet FTROI event] NOTE: GEM is proton-arm and FTROI/CDet is "
                 "electron-arm; use HEEP/elastic closure, not direct slope "
                 "subtraction, for cross-arm comparison."
              << std::endl;
    return;
  }

  std::cerr << "[CDet FTROI event] Event " << requestedEventNumber
            << " was not found in the input chain" << std::endl;
}
