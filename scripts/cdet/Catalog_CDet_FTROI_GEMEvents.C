#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace {

int ParseSegment(const TString &fileName)
{
  const Ssiz_t marker = fileName.Index("_seg");
  if (marker == kNPOS)
    return -1;
  int segment = -1;
  std::sscanf(fileName.Data() + marker + 4, "%d", &segment);
  return segment;
}

std::string JoinTag(std::string tags, const char *tag)
{
  if (!tags.empty())
    tags += ";";
  tags += tag;
  return tags;
}

double ArrayValue(const TTreeReaderArray<Double_t> &values, int index)
{
  return index >= 0 && index < static_cast<int>(values.GetSize()) ?
      values[index] : std::numeric_limits<double>::quiet_NaN();
}

void WriteNumber(std::ostream &output, double value)
{
  if (std::isfinite(value))
    output << value;
}

} // namespace

// Create an event-level CSV catalogue for every event with at least one
// reconstructed front-tracker GEM track. GEM quantities describe the proton
// arm; FTROI/CDet quantities describe the electron arm. They are intentionally
// kept side by side without treating their slopes as directly comparable.
void Catalog_CDet_FTROI_GEMEvents(
    const char *inputPattern,
    const char *outputCSV = "CDet_run5711_FTROI_GEM_event_catalog.csv")
{
  TChain chain("T");
  const int inputFiles = chain.Add(inputPattern);
  if (inputFiles <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[CDet FTROI GEM catalogue] No readable T trees match "
              << inputPattern << std::endl;
    return;
  }

  const char *required[] = {
      "g.evnum", "sbs.gemFT.track.ntrack", "sbs.gemFT.track.besttrack",
      "sbs.gemFT.track.nhits", "sbs.gemFT.track.ngoodhits",
      "sbs.gemFT.track.chi2ndf", "sbs.gemFT.track.x",
      "sbs.gemFT.track.y", "sbs.gemFT.track.xp", "sbs.gemFT.track.yp",
      "sbs.gemFT.track.t0", "sbs.tr.n", "heep.datavalid", "heep.dpe",
      "heep.dpp", "heep.dphi", "heep.acoplanarity", "heep.dxECAL",
      "heep.dyECAL", "heep.dt_ADC", "earm.ecal.e", "earm.ecal.x",
      "earm.ecal.y", "FTROI.cdet.timing_status", "FTROI.cdet.roi_status",
      "FTROI.cdet.npulse", "FTROI.cdet.npair_candidate",
      "FTROI.cdet.nsingle_candidate", "FTROI.cdet.hyp.source_type",
      "FTROI.cdet.hyp.y_topology", "FTROI.cdet.vertex.bin",
      "FTROI.cdet.vertex.hyp_index", "FTROI.cdet.vertex.z",
      "FTROI.cdet.vertex.xchi2", "FTROI.cdet.vertex.xndf",
      "FTROI.cdet.vertex.ycompatible"};
  for (const char *name : required) {
    if (!chain.GetBranch(name)) {
      std::cerr << "[CDet FTROI GEM catalogue] Required branch is missing: "
                << name << std::endl;
      return;
    }
  }

  std::ofstream output(outputCSV);
  if (!output) {
    std::cerr << "[CDet FTROI GEM catalogue] Cannot create " << outputCSV
              << std::endl;
    return;
  }
  output << "catalog_index,event_number,global_entry,local_entry,tree_number,segment,source_file,cohort,review_tier,review_tags,"
            "gem_ntrack,gem_best_index,gem_best_nhits,gem_best_ngoodhits,gem_best_chi2ndf,"
            "gem_best_x,gem_best_y,gem_best_xp,gem_best_yp,gem_best_t0,sbs_track_n,"
            "heep_datavalid,heep_dpe,heep_dpp,heep_dphi,heep_acoplanarity,"
            "heep_dxECAL,heep_dyECAL,heep_dt_ADC,ecal_e,ecal_x,ecal_y,"
            "cdet_timing_status,cdet_roi_status,cdet_npulse,cdet_npair_candidate,"
            "cdet_nsingle_candidate,ftroi_nhyp,ftroi_npair,ftroi_nsame_side,"
            "ftroi_nseam,ftroi_nsingle,ftroi_ycompatible_associations,"
            "ftroi_total_associations,best_hypothesis,best_topology,best_target_bin,"
            "best_target_z,best_xchi2ndf,best_ycompatible,best_at_scan_boundary\n";
  output << std::setprecision(10);

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
  TTreeReaderValue<Double_t> sbsTrackN(reader, "sbs.tr.n");
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
  TTreeReaderArray<Double_t> hypYTopology(reader, "FTROI.cdet.hyp.y_topology");
  TTreeReaderArray<Double_t> vertexBin(reader, "FTROI.cdet.vertex.bin");
  TTreeReaderArray<Double_t> vertexHyp(reader, "FTROI.cdet.vertex.hyp_index");
  TTreeReaderArray<Double_t> vertexZ(reader, "FTROI.cdet.vertex.z");
  TTreeReaderArray<Double_t> vertexXChi2(reader, "FTROI.cdet.vertex.xchi2");
  TTreeReaderArray<Double_t> vertexXNDF(reader, "FTROI.cdet.vertex.xndf");
  TTreeReaderArray<Double_t> vertexYCompatible(reader, "FTROI.cdet.vertex.ycompatible");

  Long64_t catalogIndex = 0;
  Long64_t overlapEvents = 0;
  Long64_t gemOnlyEvents = 0;
  Long64_t multiTrackEvents = 0;
  Long64_t overlapWithSeam = 0;
  Long64_t overlapWithSingle = 0;
  Long64_t overlapBoundaryBest = 0;
  Long64_t gemOnlyNoCDetCandidates = 0;
  Long64_t gemOnlyPairStatus = 0;
  Long64_t gemOnlyL1Status = 0;
  Long64_t gemOnlyL2Status = 0;

  while (reader.Next()) {
    const int nTrack = static_cast<int>(std::lround(*gemNTrack));
    if (nTrack <= 0)
      continue;

    const int nHyp = static_cast<int>(hypSourceType.GetSize());
    const bool overlap = nHyp > 0;
    if (overlap)
      ++overlapEvents;
    else
      ++gemOnlyEvents;
    if (nTrack > 1)
      ++multiTrackEvents;

    int nPair = 0, nSameSide = 0, nSeam = 0, nSingle = 0;
    for (int ih = 0; ih < nHyp; ++ih) {
      const int sourceType = static_cast<int>(std::lround(hypSourceType[ih]));
      const int topology = static_cast<int>(std::lround(hypYTopology[ih]));
      if (sourceType == 1) {
        ++nPair;
        if (topology == 0)
          ++nSameSide;
        else if (topology == 1)
          ++nSeam;
      } else if (sourceType == 2 || sourceType == 3) {
        ++nSingle;
      }
    }
    if (overlap && nSeam > 0)
      ++overlapWithSeam;
    if (overlap && nSingle > 0)
      ++overlapWithSingle;

    int compatibleAssociations = 0;
    int totalAssociations = 0;
    int bestHypothesis = -1;
    int bestTopology = -999;
    int bestTargetBin = -1;
    int bestYCompatible = -1;
    double bestTargetZ = std::numeric_limits<double>::quiet_NaN();
    double bestXChi2NDF = std::numeric_limits<double>::infinity();
    for (int iv = 0; iv < static_cast<int>(vertexZ.GetSize()); ++iv) {
      const int ih = static_cast<int>(std::lround(vertexHyp[iv]));
      if (ih < 0 || ih >= nHyp)
        continue;
      ++totalAssociations;
      if (vertexYCompatible[iv] == 1.0)
        ++compatibleAssociations;
      const double ndf = vertexXNDF[iv];
      const double chi2NDF = ndf > 0.0 ? vertexXChi2[iv] / ndf :
          std::numeric_limits<double>::infinity();
      if (std::isfinite(chi2NDF) && chi2NDF < bestXChi2NDF) {
        bestXChi2NDF = chi2NDF;
        bestHypothesis = ih;
        bestTopology = static_cast<int>(std::lround(hypYTopology[ih]));
        bestTargetBin = static_cast<int>(std::lround(vertexBin[iv]));
        bestTargetZ = vertexZ[iv];
        bestYCompatible = static_cast<int>(std::lround(vertexYCompatible[iv]));
      }
    }
    const bool boundaryBest = bestTargetBin == 0 || bestTargetBin == 23;
    if (overlap && boundaryBest)
      ++overlapBoundaryBest;

    const int status = static_cast<int>(std::lround(*roiStatus));
    if (!overlap) {
      if (status == 0)
        ++gemOnlyNoCDetCandidates;
      else if (status == 1)
        ++gemOnlyPairStatus;
      else if (status == 2)
        ++gemOnlyL1Status;
      else if (status == 3)
        ++gemOnlyL2Status;
    }

    int bestTrack = static_cast<int>(std::lround(*gemBestTrack));
    if (bestTrack < 0 || bestTrack >= nTrack)
      bestTrack = 0;

    const Long64_t globalEntry = reader.GetCurrentEntry();
    const Long64_t localEntry = chain.GetTree() ?
        chain.GetTree()->GetReadEntry() : -1;
    const int treeNumber = chain.GetTreeNumber();
    const TString sourceFile = chain.GetCurrentFile() ?
        chain.GetCurrentFile()->GetName() : "";
    const TString baseName = gSystem->BaseName(sourceFile);
    const int segment = ParseSegment(baseName);

    std::string tags;
    tags = JoinTag(tags, overlap ? "overlap" : "gem_only");
    if (nTrack > 1)
      tags = JoinTag(tags, "multi_track");
    if (nSeam > 0)
      tags = JoinTag(tags, "seam");
    if (nSingle > 0)
      tags = JoinTag(tags, "single");
    if (boundaryBest)
      tags = JoinTag(tags, "boundary_best_z");
    if (!overlap && status != 0)
      tags = JoinTag(tags, "candidate_without_hypothesis");

    std::string reviewTier;
    if (overlap && !boundaryBest)
      reviewTier = "1_overlap_interior_best_z";
    else if (overlap && nSingle > 0)
      reviewTier = "2_overlap_single";
    else if (overlap && nSeam > 0)
      reviewTier = "3_overlap_seam";
    else if (overlap && nTrack > 1)
      reviewTier = "4_overlap_multi_track";
    else if (overlap)
      reviewTier = "5_overlap_standard";
    else if (nTrack > 1)
      reviewTier = "6_gem_only_multi_track";
    else
      reviewTier = "7_gem_only_no_candidate";

    output << catalogIndex << ',' << std::lround(*eventNumber) << ','
           << globalEntry << ',' << localEntry << ',' << treeNumber << ','
           << segment << ','
           << baseName << ',' << (overlap ? "overlap" : "gem_only") << ','
           << reviewTier << ',' << tags << ',' << nTrack << ',' << bestTrack
           << ',';
    WriteNumber(output, ArrayValue(gemNHits, bestTrack)); output << ',';
    WriteNumber(output, ArrayValue(gemNGoodHits, bestTrack)); output << ',';
    WriteNumber(output, ArrayValue(gemChi2NDF, bestTrack)); output << ',';
    WriteNumber(output, ArrayValue(gemX, bestTrack)); output << ',';
    WriteNumber(output, ArrayValue(gemY, bestTrack)); output << ',';
    WriteNumber(output, ArrayValue(gemXP, bestTrack)); output << ',';
    WriteNumber(output, ArrayValue(gemYP, bestTrack)); output << ',';
    WriteNumber(output, ArrayValue(gemT0, bestTrack)); output << ',';
    output << std::lround(*sbsTrackN) << ',' << std::lround(*heepDataValid) << ',';
    WriteNumber(output, *heepDPE); output << ',';
    WriteNumber(output, *heepDPP); output << ',';
    WriteNumber(output, *heepDPhi); output << ',';
    WriteNumber(output, *heepAcoplanarity); output << ',';
    WriteNumber(output, *heepDXECal); output << ',';
    WriteNumber(output, *heepDYECal); output << ',';
    WriteNumber(output, *heepDtADC); output << ',';
    WriteNumber(output, *ecalE); output << ',';
    WriteNumber(output, *ecalX); output << ',';
    WriteNumber(output, *ecalY); output << ',';
    output << std::lround(*timingStatus) << ',' << status << ','
           << std::lround(*nPulse) << ',' << std::lround(*nPairCandidate) << ','
           << std::lround(*nSingleCandidate) << ',' << nHyp << ',' << nPair << ','
           << nSameSide << ',' << nSeam << ',' << nSingle << ','
           << compatibleAssociations << ',' << totalAssociations << ','
           << bestHypothesis << ',' << bestTopology << ',' << bestTargetBin << ',';
    WriteNumber(output, bestTargetZ); output << ',';
    WriteNumber(output, bestXChi2NDF); output << ',' << bestYCompatible << ','
           << (boundaryBest ? 1 : 0) << '\n';
    ++catalogIndex;
  }

  output.close();
  std::cout << "[CDet FTROI GEM catalogue] input files: " << inputFiles << std::endl;
  std::cout << "[CDet FTROI GEM catalogue] rows: " << catalogIndex << std::endl;
  std::cout << "[CDet FTROI GEM catalogue] overlap/gem-only: " << overlapEvents
            << "/" << gemOnlyEvents << std::endl;
  std::cout << "[CDet FTROI GEM catalogue] multi-track: " << multiTrackEvents
            << std::endl;
  std::cout << "[CDet FTROI GEM catalogue] overlap seam/single/boundary-best: "
            << overlapWithSeam << "/" << overlapWithSingle << "/"
            << overlapBoundaryBest << std::endl;
  std::cout << "[CDet FTROI GEM catalogue] gem-only ROI status 0/1/2/3: "
            << gemOnlyNoCDetCandidates << "/" << gemOnlyPairStatus << "/"
            << gemOnlyL1Status << "/" << gemOnlyL2Status << std::endl;
  std::cout << "[CDet FTROI GEM catalogue] wrote " << outputCSV << std::endl;
}
