#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>      // for sscanf
#include <algorithm>   // for std::sort
#include <TMath.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TPaveText.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>

std::vector<TCanvas*> canvas_vector;

static const int TDCmult_cut = 100;
static const double xcut = 998.0;
//static const int nhitcutlow = 2;
//static const int nhitcuthigh = 20;
static const double TDC_calib_to_ns = 0.01;
static const double HotChannelRatio = .01;

static const int NumPaddles = 16;
static const int NumBars = 14;
static const int NumLayers = 2;
static const int NumSides = 2;
static const int NumModules = 3;
static const int NumHalfModules = NumModules*NumSides*NumLayers;

static const int NumCDetPaddles = NumHalfModules*NumBars*NumPaddles; //2688
static const int nRef = 4;
static const int NumRefPaddles = 4;
static const int nTdc = NumCDetPaddles+NumRefPaddles; //2704

static const int NumSidesTotal = NumSides*NumLayers;
static const int NumCDetPaddlesPerSide = NumCDetPaddles/NumSidesTotal; //672
static const int NumLogicalPaddlesPerSide = NumCDetPaddlesPerSide+nRef; //676

static const int nBarsADC = 0;
static const double ADCCUT = 150.;   //100.0

static const double ECal_dist = 6.6;
static const double CDet_y_half_length = 0.30;

int NXDiffBins;
double XDiffLow;
double XDiffHigh;

// For generating sorting cdet layer 1 and 2 hits into pairs
struct Cand {
  int i1;       // index into layer-1 per-event vectors
  int j2;       // index into layer-2 per-event vectors
  double dt;    // t2 - t1
  double dx;    // x2 - x1
  double score; // ranking metric (lower = better)
};

// For storing pairs for cdet layer 1 and layer 2 hits, eventually need to add y,z, tot
struct PairHit {
  double t1, tot1, x1, y1, z1;
  int    id1;

  double t2, tot2, x2, y2, z2;
  int    id2;

  double dt, dx;
  double score;
};

std::vector<std::vector<PairHit>> pairs_CDet;

// List of x-positions (or bins) for unused pixels ----- 1/19 verified correct
static std::vector<double> missingPixelBins = {
3, 13, 28, 31, 41, 42, 57, 59, 65, 79, 83, 95, 109, 111, 115, 127,
140, 143, 145, 156, 172, 175, 176, 188, 195, 199, 213, 220, 236, 239, 244, 255,
268, 271, 284, 287, 300, 303, 307, 319, 332, 335, 339, 351, 354, 364, 371, 381,
384, 396, 401, 410, 419, 423, 435, 436, 451, 461, 465, 479, 480, 483, 508, 511,
512, 515, 540, 543, 546, 559, 563, 573, 576, 589, 596, 605, 609, 610, 627, 638,
643, 655, 656, 665, 674, 675, 696, 703, 707, 709, 725, 729, 738, 748, 752, 766,
777, 780, 784, 791, 800, 812, 818, 828, 844, 847, 850, 860, 867, 868, 884, 885,
900, 904, 912, 927, 940, 943, 945, 947, 967, 971, 986, 991, 1005, 1007, 1011, 1023,
1027, 1028, 1043, 1050, 1066, 1068, 1072, 1075, 1088, 1102, 1106, 1119, 1121, 1135, 1148, 1151,
1162, 1166, 1178, 1182, 1184, 1186, 1203, 1215, 1228, 1231, 1235, 1247, 1249, 1252, 1267, 1274,
1282, 1285, 1299, 1310, 1317, 1321, 1340, 1341, 1349, 1359, 1372, 1375, 1376, 1391, 1392, 1405,
1409, 1420, 1428, 1439, 1443, 1455, 1468, 1471, 1486, 1487, 1500, 1503, 1516, 1519, 1520, 1523,
1536, 1551, 1557, 1567, 1568, 1583, 1584, 1597, 1603, 1615, 1617, 1629, 1646, 1647, 1648, 1654,
1676, 1679, 1692, 1695, 1708, 1711, 1715, 1725, 1732, 1743, 1744, 1757, 1761, 1770, 1778, 1786,
1804, 1807, 1820, 1823, 1836, 1839, 1854, 1855, 1856, 1868, 1877, 1887, 1902, 1903, 1916, 1919,
1934, 1935, 1942, 1951, 1964, 1967, 1973, 1983, 1988, 1999, 2000, 2013, 2028, 2031, 2034, 2047,
2048, 2051, 2064, 2067, 2080, 2085, 2099, 2104, 2112, 2122, 2131, 2143, 2144, 2147, 2160, 2163,
2177, 2188, 2202, 2207, 2208, 2221, 2227, 2239, 2243, 2254, 2259, 2271, 2279, 2283, 2300, 2303,
2307, 2316, 2320, 2334, 2339, 2348, 2355, 2367, 2369, 2383, 2384, 2395, 2405, 2409, 2416, 2422,
2432, 2435, 2448, 2451, 2464, 2479, 2483, 2493, 2499, 2508, 2512, 2513, 2531, 2537, 2544, 2547,
2563, 2570, 2576, 2591, 2592, 2607, 2611, 2621, 2633, 2636, 2643, 2650, 2656, 2657, 2675, 2679};

static const std::unordered_set<int> kUnusedCDetPixels = {
3, 13, 28, 31, 41, 42, 57, 59, 65, 79, 83, 95, 109, 111, 115, 127,
140, 143, 145, 156, 172, 175, 176, 188, 195, 199, 213, 220, 236, 239, 244, 255,
268, 271, 284, 287, 300, 303, 307, 319, 332, 335, 339, 351, 354, 364, 371, 381,
384, 396, 401, 410, 419, 423, 435, 436, 451, 461, 465, 479, 480, 483, 508, 511,
512, 515, 540, 543, 546, 559, 563, 573, 576, 589, 596, 605, 609, 610, 627, 638,
643, 655, 656, 665, 674, 675, 696, 703, 707, 709, 725, 729, 738, 748, 752, 766,
777, 780, 784, 791, 800, 812, 818, 828, 844, 847, 850, 860, 867, 868, 884, 885,
900, 904, 912, 927, 940, 943, 945, 947, 967, 971, 986, 991, 1005, 1007, 1011, 1023,
1027, 1028, 1043, 1050, 1066, 1068, 1072, 1075, 1088, 1102, 1106, 1119, 1121, 1135, 1148, 1151,
1162, 1166, 1178, 1182, 1184, 1186, 1203, 1215, 1228, 1231, 1235, 1247, 1249, 1252, 1267, 1274,
1282, 1285, 1299, 1310, 1317, 1321, 1340, 1341, 1349, 1359, 1372, 1375, 1376, 1391, 1392, 1405,
1409, 1420, 1428, 1439, 1443, 1455, 1468, 1471, 1486, 1487, 1500, 1503, 1516, 1519, 1520, 1523,
1536, 1551, 1557, 1567, 1568, 1583, 1584, 1597, 1603, 1615, 1617, 1629, 1646, 1647, 1648, 1654,
1676, 1679, 1692, 1695, 1708, 1711, 1715, 1725, 1732, 1743, 1744, 1757, 1761, 1770, 1778, 1786,
1804, 1807, 1820, 1823, 1836, 1839, 1854, 1855, 1856, 1868, 1877, 1887, 1902, 1903, 1916, 1919,
1934, 1935, 1942, 1951, 1964, 1967, 1973, 1983, 1988, 1999, 2000, 2013, 2028, 2031, 2034, 2047,
2048, 2051, 2064, 2067, 2080, 2085, 2099, 2104, 2112, 2122, 2131, 2143, 2144, 2147, 2160, 2163,
2177, 2188, 2202, 2207, 2208, 2221, 2227, 2239, 2243, 2254, 2259, 2271, 2279, 2283, 2300, 2303,
2307, 2316, 2320, 2334, 2339, 2348, 2355, 2367, 2369, 2383, 2384, 2395, 2405, 2409, 2416, 2422,
2432, 2435, 2448, 2451, 2464, 2479, 2483, 2493, 2499, 2508, 2512, 2513, 2531, 2537, 2544, 2547,
2563, 2570, 2576, 2591, 2592, 2607, 2611, 2621, 2633, 2636, 2643, 2650, 2656, 2657, 2675, 2679};

inline bool IsUnusedPixel(int elID) {
    return kUnusedCDetPixels.count(elID) != 0;
}


//const TString REPLAYED_DIR = TString(gSystem->Getenv("OUT_DIR")) + "/wrongdbRootfiles";
const TString REPLAYED_DIR = TString(gSystem->Getenv("OUT_DIR"));

// const TString ANALYSED_DIR = gSystem->Getenv("ANALYSED_DIR");
//const TString REPLAYED_DIR = "/volatile/halla/sbs/btspaude/cdet/rootfiles";
const TString ANALYSED_DIR = TString(gSystem->Getenv("OUT_DIR"))+"cdetFiles";

// Parse the "segX_Y" part: returns true and fills firstSeg/lastSeg if found.
bool GetSegRange(const TString& fname, int& firstSeg, int& lastSeg) {
  // Find "_seg"
  Ssiz_t pos = fname.Index("_seg");
  if (pos == kNPOS) return false;

  // Tail looks like "9_9.root" or "9_9_1.root"
  TString tail = fname(pos + 4, fname.Length() - (pos + 4));

  // Extract first two ints; ignore any further suffix
  int a = -1, b = -1;
  if (sscanf(tail.Data(), "%d_%d", &a, &b) == 2) {
    firstSeg = a; lastSeg = b;
    return true;
  }
  return false;
}

void AddRunFilesToChain(TChain *chain, const char *dir, int runnum, int segMin = -1, int segMax = -1) {
  TString prefix = dir;
  std::vector<TString> runfiles;

  TSystemDirectory directory(prefix, prefix);
  TList *files = directory.GetListOfFiles();

  if (files) {
    TIter next(files);
    TSystemFile *f;
    while ((f = (TSystemFile*) next())) {
      if (f->IsDirectory()) continue; // skip dirs like "." and ".."

      TString fname = f->GetName();
      if (!fname.BeginsWith(Form("cdet_%d_", runnum))) continue;
      if (!fname.EndsWith(".root")) continue;

      // Range filtering enabled only if segMin/segMax are set
      if (segMin >= 0 || segMax >= 0) {
        if (segMin < 0) segMin = segMax;
        if (segMax < 0) segMax = segMin;
        if (segMin > segMax) std::swap(segMin, segMax);

        int firstSeg = -1, lastSeg = -1;
        if (!GetSegRange(fname, firstSeg, lastSeg)) continue;

        // accept if [firstSeg,lastSeg] overlaps [segMin,segMax]
        if (lastSeg < segMin || firstSeg > segMax) continue;
      }

      runfiles.push_back(prefix + "/" + fname);
    }
  }

  std::sort(runfiles.begin(), runfiles.end());

  std::cout << "Adding " << runfiles.size() << " files for run " << runnum << "...\n";
  for (auto &file : runfiles) {
    std::cout << "  " << file << "\n";
    chain->Add(file);
  }
}

/* Create globals for vectors */
// Scalars (1D vectors)
std::vector<double> vheep_dpp;
std::vector<double> vheep_dt_ADC;
std::vector<double> vheep_ECalo;
std::vector<double> vheep_eprime_eth;
std::vector<double> vheep_dxECAL;
std::vector<double> vearm_ECal_x;
std::vector<double> vsbs_gemFPP_track_ntrack;
std::vector<double> vheep_dyECAL;

// Arrays (2D vectors)
std::vector<std::vector<double>> vsbs_tr_vz;
std::vector<std::vector<double>> vsbs_gemFPP_track_sclose;
std::vector<std::vector<double>> vsbs_gemFT_track_nhits;
std::vector<std::vector<double>> vsbs_gemFT_track_ngoodhits;

/* CDet & ECal Vectors */
//1D vectors
std::vector<double> vRefRawLe;
std::vector<double> vGoodRefRawLe; // ref LE time aligned to GOOD-event vectors (same indexing as vGoodLe)
std::vector<double> vRefRawTe;
std::vector<double> vRefRawTot;
std::vector<int>    vRefRawPMT;

std::vector<double> vRefGoodLe;
std::vector<double> vRefGoodTe;
std::vector<double> vRefGoodTot;
std::vector<int>    vRefGoodPMT;

std::vector<double> vAllRawLe;
std::vector<double> vAllRawTe;
std::vector<double> vAllRawTot;
std::vector<int> vAllRawPMT;
std::vector<int> vAllRawBar;

std::vector<double> vAllGoodLe;
std::vector<double> vAllGoodECalT;                 // per-hit ECal ADC time aligned with vAllGoodLe/Te
std::vector<std::vector<double>> vBarGoodLeECalT;        // per-bar per-hit ECal time aligned with vBarGoodLe

// --- ECal-time linear correction: remove correlation t_CDet vs t_ECal, then apply global shift
double gECalFitP0 = -40.303;    // p0 from fit: <t_CDet> = p0 + p1*t_ECal
double gECalFitP1 =  0.63615;   // p1 from fit
double gTargetMeanLE = 30.0;    // desired mean corrected LE (ns)
double gECalDeltaShift = 0.0;   // computed shift applied after removing correlation
bool   gUseECalTimeCorr = true; // enable/disable ECal-time correction
std::vector<std::vector<double>> vBarGoodLe;
std::vector<TH1F*> hBarGoodLe; // one histogram per bar (PMT group), built in plotAllTDC()
static const int NumPMTs = NumHalfModules*NumBars; // 168 (does not include 4 ref paddles)

std::vector<double> gBarToffsetCorr;      // size NumPMTs, correction to ADD to LE/TE: (mean_all - mean_bar)
bool gBarToffsetLoaded = false;           // true if offsets were read from file
std::string gBarToffsetFile = "CDet_bar_toffsets.dat";

inline double GetBarToffsetCorr(int elID) {
  const int bar = elID / 16;
  if (0 <= bar && bar < NumPMTs && (int)gBarToffsetCorr.size() == NumPMTs) return gBarToffsetCorr[bar];
  return 0.0;
}

// Reads: bar(1..168)  toffset(ns)  [optional entries]
// Lines starting with '#' are ignored.
// Returns true if file opened and at least one bar was read.
bool LoadBarToffsets(const std::string& fname) {
  gBarToffsetCorr.assign(NumPMTs, 0.0);
  gBarToffsetLoaded = false;

  std::ifstream fin(fname.c_str());
  if (!fin) {
    std::cout << "[CDet] No offset file '" << fname << "' found; proceeding with zero bar offsets.\n";
    return false;
  }

  int nread = 0;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    std::istringstream iss(line);
    int bar1 = -1;
    double dt = 0.0;
    double dummyEntries = 0.0;
    if (!(iss >> bar1 >> dt)) continue;
    if (bar1 >= 1 && bar1 <= NumPMTs) {
      gBarToffsetCorr[bar1-1] = dt;
      nread++;
    }
  }

  if (nread > 0) {
    gBarToffsetLoaded = true;
    std::cout << "[CDet] Loaded " << nread << " bar offsets from '" << fname << "'.\n";
    return true;
  }

  std::cout << "[CDet] Offset file '" << fname << "' contained no readable offsets; using zeros.\n";
  return false;
}

std::vector<double> vAllGoodTe;
std::vector<double> vAllGoodTot;
std::vector<int> vAllGoodPMT;
std::vector<int> vAllGoodBar;

std::vector<int> vhitCDetPMT;
std::vector<int> vRow;
std::vector<std::vector<int>> vGoodCol;
std::vector<std::vector<int>> vGoodLayer;

std::vector<std::vector<double>> vCDetX;
std::vector<std::vector<double>> vCDetY;
std::vector<std::vector<double>> vCDetZ;

std::vector<std::vector<double>> vCDetGoodX;
std::vector<std::vector<double>> vCDetGoodY;
std::vector<std::vector<double>> vCDetGoodZ;
std::vector<Long64_t> vTreeEntry; // maps saved passing-event index -> TTree entry

std::vector<int> vRowLayer1Side1;
std::vector<int> vRowLayer2Side1;
std::vector<int> vRowLayer1Side2;
std::vector<int> vRowLayer2Side2;

std::vector<int> vnhits1;
std::vector<int> vnpaddles;
std::vector<int> vngoodpaddles;
std::vector<int> vngoodTDCpaddles;
std::vector<int> vngoodhits1;
std::vector<int> vngoodTDChits1;
std::vector<int> vnhits2;
std::vector<int> vngoodhits2;
std::vector<int> vngoodTDChits2;

std::vector<std::vector<double>> vCDetPaddleRawTot;
std::vector<std::vector<double>> vCDetPaddleCutTot;

//2D vectors
std::vector<std::vector<double>> vRawLe;
std::vector<std::vector<double>> vRawTe;
std::vector<std::vector<double>> vRawTot;
std::vector<std::vector<int>> vRawID;

std::vector<std::vector<double>> vGoodLe;
std::vector<std::vector<double>> vGoodTe;
std::vector<std::vector<double>> vGoodTot;
std::vector<std::vector<int>> vGoodID;

// per-event 1D
std::vector<double> vCDetMultAll;          // like hMultiplicity
std::vector<int>    vCDetMultAllPMT;       // ids belonging to those mult values

// 2D vectors for all ECal cluster arrays
std::vector<std::vector<double>> v_ECal_clus_adctime;
std::vector<std::vector<double>> v_ECal_clus_again;
std::vector<std::vector<double>> v_ECal_clus_atimeblk;
std::vector<std::vector<double>> v_ECal_clus_col;
std::vector<std::vector<double>> v_ECal_clus_col_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_e;
std::vector<std::vector<double>> v_ECal_clus_e_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_eblk;
std::vector<std::vector<double>> v_ECal_clus_eblk_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_id;
std::vector<std::vector<double>> v_ECal_clus_id_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_nblk;
std::vector<std::vector<double>> v_ECal_clus_nblk_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_row;
std::vector<std::vector<double>> v_ECal_clus_row_goodtdc;
std::vector<std::vector<double>> v_ECal_clus_tdctime;
std::vector<std::vector<double>> v_ECal_clus_tdctime_tw;
std::vector<std::vector<double>> v_ECal_clus_tdctimeblk;
std::vector<std::vector<double>> v_ECal_clus_tdctimeblk_tw;
std::vector<std::vector<double>> v_ECal_clus_x;
std::vector<std::vector<double>> v_ECal_clus_y;

// 1D (event-wise) vectors ECal for scalars:
std::vector<double> v_ECal_nclus;
std::vector<double> v_ECalX;
std::vector<double> v_ECalY;
std::vector<double> v_ECalE;
std::vector<double> v_ECalAdcTime;
//1D Good ECal Events
std::vector<double> v_GoodECalX;
std::vector<double> v_GoodECalY;
std::vector<double> v_GoodECalE;
std::vector<double> v_GoodECalAdcTime;

struct CDetHit {
  int    id;     // pixelID
  double le_ns;  // LE
  double tot_ns; // TOT
  double te_ns;  // TE
};

struct AdjPair {
  int event;      // event number/index
  int id1, id2;   // id1 < id2
  double le1, te1, tot1;
  double le2, te2, tot2;
  int i1, i2;     // hit indices within that pixel for this event (for dedupe)
};

static std::vector<AdjPair> vAdjPairs;

// using Hit = std::pair<int,double>; // (pixelID, tot_ns)
std::vector<std::vector<CDetHit>> vEventHits; // [event][hit]
std::vector<std::vector<CDetHit>> vGoodEventHits;

std::vector<int> rawRate(2688, 0); 
int rateEvTrack = 0;
std::vector<double> chanRates(2688,0);
std::vector<int> cutRate(2688, 0); 
int cutRateEvTrack = 0;
std::vector<double> cutChanRates(2688,0);
std::vector<double> ave_tot(2688,0);
std::vector<int> vNumRawAdjacentHits;
std::vector<int> vNumGoodAdjacentHits;

//copy a TTreeReaderArray<double> into a std::vector<double>, makes it easier to fill the 2D vector
inline std::vector<double> copyArray(const TTreeReaderArray<double>& arr) {
  std::vector<double> v;
  v.reserve(arr.GetSize());
  for (size_t i = 0; i < arr.GetSize(); ++i) v.push_back(arr[i]);
  return v;
}


// per-event vectors grouped by PMT index (2D)
std::vector<std::vector<double>> vCDetMultPerPMT(nTdc); 


/*
namespace TCDet {
  Int_t NdataMult;
  Double_t TDCmult[nTdc*2];

  Int_t NdataRawElID;
  Double_t RawElID[nTdc*2];
  Int_t NdataRawElLE;
  Double_t RawElLE[nTdc*2];
  Int_t NdataRawElTE;
  Double_t RawElTE[nTdc*2];
  Int_t NdataRawElTot;
  Double_t RawElTot[nTdc*2];
  
  Int_t NdataGoodRow;
  Double_t GoodRow[nTdc*2];
  Int_t NdataGoodCol;
  Double_t GoodCol[nTdc*2];
  Int_t NdataGoodLayer;
  Double_t GoodLayer[nTdc*2];

  Int_t NdataGoodElID;
  Double_t GoodElID[nTdc*2];
  Int_t NdataGoodElLE;
  Double_t GoodElLE[nTdc*2];
  Int_t NdataGoodElTE;
  Double_t GoodElTE[nTdc*2];
  Int_t NdataGoodElTot;
  Double_t GoodElTot[nTdc*2];

  Int_t NdataGoodX;
  Double_t GoodX[nTdc*2];
  Int_t NdataGoodY;
  Double_t GoodY[nTdc*2];
  Int_t NdataGoodZ;
  Double_t GoodZ[nTdc*2];

  Double_t ECalX;
  Double_t ECalY;
  Double_t ECalE;
  Double_t nhits;
  Double_t ngoodhits;
  Double_t ngoodTDChits;



};*/

//===================================================== Globals for paddle hits
Double_t nhits_paddles[nTdc*2];
Double_t ngoodhits_paddles[nTdc*2];
Double_t ngoodTDChits_paddles[nTdc*2];
Double_t npaddles;
Double_t ngoodpaddles;
Double_t ngoodTDCpaddles;

TChain *T = 0;
  
//===================================================== Histogram Declarations
// number of histo bins
const int NTotBins = 200;
const double TotBinLow = 1.;
const double TotBinHigh = 51.;
const int RefNTotBins = 800;
const double RefTotBinLow = 1.;
const double RefTotBinHigh = 201.;

double TDCBinLow;
double TDCBinHigh;
int NTDCBins;
double RefTDCBinLow;
double RefTDCBinHigh;
int RefNTDCBins;


//const int num_bad = 0;
//

const int num_bad = 5;

const int bad_channels[] = {
	1161, 1472, 1670, 1696, 1896  
};

//const int num_bad = 24;
//
//const int bad_channels[] = {
//	61, 
//	64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
//	408,417,1286,
//	2124, 2404,2406,2414
//};

//const int num_bad = 63;

//const int bad_channels[] = {
//	35,40,42,45, 
//	80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,
//	136,169,183,184,192,193,194,198,199,
//	314,
//	390,391,392,395,397,433,506,507,
//	678,800,902,933,949,
//	1177,1205,1215,1255,
//	1296,1297,1299,1302,1303,1304,1307,1311,
//	1912,2140,
//	2420,2422,2430,
//	2629,2630,2662
//};
	

/*
const int bad_channels[] = {
		45,   224, 226, 234, 241, 284, 240, 604, 690, 698, 702,
		864, 1174,1437,1488,1726,1912,2131,2132,2213,2214,
		2411,2417,2425,2426,2427,2434,2436,2437,2438,2448, 2457,
		2466,2467,2468,2469,2470,2471,2472,2473,2474,2475,
		2476,2477,2478,2479,2502,2517,2518,2519,2520,2522,
		2592,2607,2613,2614,2616,2620,2624,2628,2632,2633,
		2634,2656,2661,2664,2665,2666,2667,2673,2675,2681,
		2683,2687};
*/

// Raw hits ie all hits
TH1F *hRawLe[nTdc];
TH1F *hRawTe[nTdc];
TH1F *hRawTot[nTdc];
TH1F *hGoodLe[nTdc];
TH1F *hGoodTe[nTdc];
TH1F *hGoodTot[nTdc];

TH1F *hAllRawLe;
TH1F *hAllRawTe;
TH1F *hAllRawTot;
TH1F *hAllRawPMT;
TH1F *hAllRawBar;
TH1F *hAllGoodLe;
TH1F *hAllGoodTe;
TH1F *hAllGoodTot;
TH1F *hAllGoodPMT;
TH1F *hAllGoodBar;

TH2F *h2AllGoodLe;
TH2F *h2AllGoodTe;
TH2F *h2AllGoodTot;

TH2F *h2TDCTOTvsLE;
TH2F *h2CDetX1vsX2;

TH2F *h2TOTvsXDiff1;
TH2F *h2TOTvsXDiff2;
TH2F *h2LEvsXDiff1;
TH2F *h2LEvsXDiff2;

TH2F *hBarRateHV;

TH1F *hRefRawLe;
TH1F *hRefRawTe;
TH1F *hRefRawTot;
TH1F *hRefRawPMT;
TH1F *hRefGoodLe;
TH1F *hRefGoodTe;
TH1F *hRefGoodTot;
TH1F *hRefGoodPMT;

TH1F *hMultiplicityL[nTdc];
TH1F *hMultiplicity;


// hit channel id
TH1F *hHitPMT;
TH1F *hRow;
TH1F *hRowLayer1Side1;
TH1F *hRowLayer1Side2;
TH1F *hRowLayer2Side1;
TH1F *hRowLayer2Side2;
TH1F *hLayer;
TH1F *hCol;

TH1F *hnhits1;
TH1F *hngoodhits1;
TH1F *hngoodTDChits1;
TH1F *hnhits2;
TH1F *hngoodhits2;
TH1F *hngoodTDChits2;
TH1F *hnhits_ev;
TH1F *hngoodhits_ev;
TH1F *hngoodTDChits_ev;

TH1F *hnpaddles;
TH1F *hngoodpaddles;
TH1F *hngoodTDCpaddles;

TH1F *hHitX;
TH1F *hHitY;
TH1F *hHitZ;

TH2F *hHitXY1;
TH2F *hHitXY2;

TH1F *hXECal;
TH1F *hYECal;
TH1F *hEECal;

TH2F *hXECalCDet1;
TH2F *hXECalCDet2;
TH2F *hXECalCDet1_min;
TH2F *hXECalCDet2_min;
TH2F *hYECalCDet1;
TH2F *hYECalCDet2;
TH2F *hEECalCDet1;
TH2F *hEECalCDet2;

TH2F *hXYECal;

TH1F *hXDiffECalCDet1;
TH1F *hXPlusECalCDet1;
TH1F *hXDiffECalCDet2;
TH1F *hXPlusECalCDet2;

TH2F *hXCDet1CDet2;
  
// 2D histograms
TH2F* h2d_RawLE;
TH2F* h2d_RawTE;
TH2F* h2d_RawTot;

TH2F* h2d_GoodLE;
TH2F* h2d_GoodTE;
TH2F* h2d_GoodTot;

TH2F* h2d_Mult;

using namespace std;

bool check_bad(int pmt, bool suppress_bad) {
	bool flag = false;
	if (!suppress_bad) return flag;
	for (int i=0;i<num_bad;i++) {
		if (pmt == bad_channels[i])  flag = true;
	}
	return flag;
}

std::vector<double> extractBinContents(const TH1* hist) {
    if (!hist) {
        throw std::invalid_argument("Null histogram pointer passed.");
    }

    int nBins = hist->GetNbinsX();
    std::vector<double> contents;

    // Loop over all *visible* bins (skip underflow bin 0 and overflow bin nBins+1)
    for (int i = 1; i <= nBins; ++i) {
        contents.push_back(hist->GetBinContent(i));
    }

    return contents;
}


vector<vector<double>> readDataFromFiles(const vector<string>& filenames) {
    const int NUM_VALUES = 42;
    vector<vector<double>> allData;

    for (const auto& filename : filenames) {
        ifstream infile(filename);
        if (!infile) {
            cerr << "Error opening file: " << filename << endl;
            continue; // Skip this file and move on
        }

        vector<double> fileData;
        double value;
        for (int i = 0; i < NUM_VALUES; ++i) {
            infile >> value;
            if (!infile) {
                cerr << "Error reading value " << i << " from file: " << filename << endl;
                break; // Stop reading this file
            }
            fileData.push_back(value);
        }

        if (fileData.size() == NUM_VALUES) {
            allData.push_back(fileData);
        } else {
            cerr << "Incomplete data in file: " << filename << endl;
        }

        infile.close();
    }

    return allData;
}

int getPixelID(int layerNum, int sideNum, int submoduleNum, int pmtNum, int pixelNum){
  // Calculate paddle number, note that missing pixels are included here
  // Validate inputs
  if (layerNum < 1 || layerNum > 2 ||
    sideNum < 1 || sideNum > 2 ||
    submoduleNum < 1 || submoduleNum > 3 ||
    pmtNum < 1 || pmtNum > 14 ||
    pixelNum < 1 || pixelNum > 16) {
    std::cerr << "Error: Invalid input values.\n"
              << "  layerNum must be 1 or 2\n"
              << "  sideNum must be 1 or 2\n"
              << "  submoduleNum must be 1 to 3\n"
              << "  pmtNum must be 1 to 14\n"
              << "  pixelNum must be 1 to 16\n";
    return -1;  // Error code
}

  int pixel = (layerNum - 1) * 1344 + //1344 pixels per layer, 0-1343 in layer 1, 1344-2687 in layer 2
  (submoduleNum - 1) * 224 + //224 pixels per side of a module
  (sideNum - 1) * 672 + //672 pixels per side
  (pmtNum-1) * 16 + //16 pixels per pmt
  pixelNum;
  return pixel - 1;
}

std::vector<int> getLocation(int pixelID) {
  // Check valid range
  if (pixelID < 0 || pixelID > 2687) {
      std::cerr << "Error: pixelID must be in the range 0 to 2687.\n";
      return {};  // return empty vector to signal error
  }

  int layerNum      = pixelID / 1344; //1344 pixels per layer
  pixelID               %= 1344;

  int sideNum       = pixelID / 672;  //672 pixels per side
  pixelID               %= 672;

  int submoduleNum  = pixelID / 224; //224 pixels per side of module
  pixelID               %= 224;

  int pmtNum        = pixelID / 16; //16 pixels per bar
  pixelID               %= 16;
  int pixelNum      = pixelID % 16;

  return {layerNum, sideNum, submoduleNum, pmtNum, pixelNum};
}
//Used to fill 2D arrays 
template<typename T>
std::vector<T> fill2D(const TTreeReaderArray<T>& arr) {
  std::vector<T> tmp;
  int n = arr.GetSize();
  tmp.reserve(n);
  for (int i=0; i<n; i++) tmp.push_back(arr[i]);
  return tmp;
}

void PlotElastic_with_barhist_offsets_io_ecalcorr_apply(Int_t RunNumber1=5811, Int_t nevents=50000, Int_t elastic = 0, Int_t minSeg = -1, Int_t maxSeg = -1,
	Double_t LeMin = 0.02, Double_t LeMax = 60.0,
	Double_t TotMin = 1.0, Double_t TotMax = 150.0, 
	Int_t nhitcutlow1 = 1, Int_t nhitcuthigh1 = 100,
	Int_t nhitcutlow2 = 1, Int_t nhitcuthigh2 = 100,
	Double_t XDiffCut = 0.02, Double_t XOffset = 0.02, Double_t YOffset = 0.1,
        Int_t layer_choice=3,	
	bool suppress_bad = false,
	Int_t nruns=30, Int_t maxstream = 2, Int_t firstevent = 1)
{
  Int_t nseg = nruns/(maxstream+1);
	Double_t RefLeMin = 1.0;
	Double_t RefLeMax = 251.0;
	RefNTDCBins = (RefLeMax-RefLeMin)/4;
	Double_t RefTotMin = 1.0;
	Double_t RefTotMax = 251.0;

	NTDCBins = 2*(LeMax-LeMin)/.0160167; // 4 ns is the trigger time, 0.018 ns is the expected time resolution, if we use a reference TDC ? 
					// 4 ns resolution is the best we can hope for, I think, using only the module trigger time.

	NXDiffBins = (int)((2*XDiffCut)/0.0073);
	XDiffLow = XOffset-XDiffCut;
	XDiffHigh = XOffset+XDiffCut;

  // InFile is the input file without absolute path and without .root suffix
  // nevents is how many events to analyse, -1 for all
  
  // To execute
  // root -l
  // .L PlotRawTDC2D.C+
  // PlotRawTDC2D("filename", -1)
  TDCBinLow = LeMin;
  TDCBinHigh = LeMax;
  RefTDCBinLow = RefLeMin;
  RefTDCBinHigh = RefLeMax;


// Load bar timing offsets if available (applied to CDet LE/TE times).
// If the file does not exist, offsets default to 0 and the macro will generate it from the data.
LoadBarToffsets(gBarToffsetFile);

  
  // hit channel id

  hHitPMT = new TH1F("hHitPMT","hHitPMT",nTdc,0,nTdc);

  hnhits1 = new TH1F("hnhits1","hnhits1",150,1,151);
  hngoodhits1 = new TH1F("hngoothits1","hngoodhits1",75,1,76);
  hngoodTDChits1 = new TH1F("hngoodTDChits1","hngoodTDChits1",75,1,76);

  hnhits2 = new TH1F("hnhits2","hnhits2",150,1,151);
  hngoodhits2 = new TH1F("hngoothits2","hngoodhits2",75,1,76);
  hngoodTDChits2 = new TH1F("hngoodTDChits2","hngoodTDChits2",75,1,76);

  hnhits_ev = new TH1F("hnhits_ev","hnhits_ev",500,0,50000);
  hngoodhits_ev = new TH1F("hngoothits_ev","hngoodhits_ev",500,0,50000);
  hngoodTDChits_ev = new TH1F("hngoodTDChits_ev","hngoodTDChits_ev",500,0,50000);

  hnpaddles = new TH1F("hnpaddles","hnpaddles",200,1,201);
  hngoodpaddles = new TH1F("hngoodpaddles","hnTDCpaddles",75,1,71);
  hngoodTDCpaddles = new TH1F("hngoodTDCpaddles","hngoodTDCpaddles",75,1,76);

  hRow = new TH1F("RowNumber","RowNumber",680, 0, 680);
  hRowLayer1Side1 = new TH1F("RowNumberL1S1","RowNumberL1S1",680, 0, 680);
  hRowLayer1Side2 = new TH1F("RowNumberL1S2","RowNumberL1S2",680, 0, 680);
  hRowLayer2Side1 = new TH1F("RowNumberL2S1","RowNumberL2S1",680, 0, 680);
  hRowLayer2Side2 = new TH1F("RowNumberL2S2","RowNumberL2S2",680, 0, 680);
  hLayer = new TH1F("LayerNumber","LayerNumber",3, 0, 3);
  hCol = new TH1F("ColNumber","ColNumber",3, 0, 3);

  hHitX = new TH1F("HitXposition","HitXPosition",1000,-2.0,2.0);
  hHitY = new TH1F("HitYposition","HitYPosition",200,-0.5,0.5);
  hHitZ = new TH1F("HitZposition","HitZPosition",200,7.5,8.0);
  
  hHitXY1 = new TH2F("HitXY1position","HitXY1Position",9,-1.0,1.0,800,-2.0,2.0);
  hHitXY2 = new TH2F("HitXY2position","HitXY2Position",9,-1.0,1.0,800,-2.0,2.0);
  
  hXECal = new TH1F("XECal","XECal",200,-1.5,1.5);
  hYECal = new TH1F("YECal","YECal",200,-1.0,1.0);
  hEECal = new TH1F("EECal","EECal",200,0.0,20.0);
  
  hXECalCDet1 = new TH2F("XECalCDet1","XECalCDet1",100,-2.0,2.0,100,-2.0,2.0);
  hXECalCDet2 = new TH2F("XECalCDet2","XECalCDet2",100,-2.0,2.0,100,-2.0,2.0);
  hXECalCDet1_min = new TH2F("XECalCDet1_min","XECalCDet1_min (min |x_{CDet}-x_{ECal->CDet}| per event)",100,-2.0,2.0,100,-2.0,2.0);
hXECalCDet2_min = new TH2F("XECalCDet2_min","XECalCDet2_min (min |x_{CDet}-x_{ECal->CDet}| per event)",100,-2.0,2.0,100,-2.0,2.0);
  hYECalCDet1 = new TH2F("YECalCDet1","YECalCDet1",100,-1.0,1.0,9,-1.0,1.0);
  hYECalCDet2 = new TH2F("YECalCDet2","YECalCDet2",100,-1.0,1.0,9,-1.0,1.0);
  hEECalCDet1 = new TH2F("EECalCDet1","EECalCDet1",100,0.0,20.0,100,-2.0,2.0);
  hEECalCDet2 = new TH2F("EECalCDet2","EECalCDet2",100,0.0,20.0,100,-2.0,2.0);
  
  hXYECal = new TH2F("XYECal","XYECal",200,-2.0,2.0,200,-2.0,2.0);
  
  hXDiffECalCDet1 = new TH1F("XDiffECalCDet1","XDiffECalCDet1",NXDiffBins,XDiffLow,XDiffHigh);
  hXPlusECalCDet1 = new TH1F("XPlusECalCDet1","XPlusECalCDet1",NXDiffBins,XDiffLow,XDiffHigh);
  hXDiffECalCDet2 = new TH1F("XDiffECalCDet2","XDiffECalCDet2",NXDiffBins,XDiffLow,XDiffHigh);
  hXPlusECalCDet2 = new TH1F("XPlusECalCDet2","XPlusECalCDet2",NXDiffBins,XDiffLow,XDiffHigh);
  
  hXCDet1CDet2 = new TH2F("XCDet1CDet2","XCDet1CDet2",200,-0.5,0.5,200,-0.5,0.5);
  
  // 2D histograms
  h2d_RawLE  = new TH2F("h2d_RawLE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTE  = new TH2F("h2d_RawTE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTot = new TH2F("h2d_RawTot","", NTotBins,TotBinLow,TotBinHigh,nTdc+1,0,nTdc+1);
  h2d_Mult   = new TH2F("h2d_Mult","", 100,0,100,nTdc+1,0,nTdc+1);
  
  hMultiplicity = new TH1F("hMultiplicity","hMultiplicity",20,0,20);
  
  for(Int_t tdc=0; tdc<nTdc; tdc++){
    hMultiplicityL[tdc] =  new TH1F(TString::Format("hMultiplicity_Bar%d",tdc),
		      TString::Format("hMultiplicity_Bar%d",tdc),
		      10, 0, 10);
  }// element loop

  hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllRawTe = new TH1F(TString::Format("hRawTe"),
            TString::Format("hRawTe"),
            NTDCBins, TDCBinLow, TDCBinHigh+TotBinHigh);
  hAllRawTot = new TH1F(TString::Format("hRawTot"),
            TString::Format("hRawTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllRawPMT = new TH1F(TString::Format("hRawPMT"),
            TString::Format("hRawPMT"),
            nTdc, 0, nTdc);
  hAllRawBar = new TH1F(TString::Format("hRawBar"),
            TString::Format("hRawBar"),
            168, 0, 168);
  hAllGoodLe = new TH1F(TString::Format("hAllGoodLe"),
            TString::Format("hAllGoodLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllGoodTe = new TH1F(TString::Format("hAllGoodTe"),
            TString::Format("hAllGoodTe"),
            NTDCBins, TDCBinLow, TDCBinHigh+TotBinHigh);
  hAllGoodTot = new TH1F(TString::Format("hAllGoodTot"),
            TString::Format("hAllGoodTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllGoodPMT = new TH1F(TString::Format("hAllGoodPMT"),
            TString::Format("hAllGoodPMT"),
            nTdc, 0, nTdc);
  hAllGoodBar = new TH1F(TString::Format("hAllGoodBar"),
            TString::Format("hAllGoodBar"),
            168, 0, 168);
  h2AllGoodLe = new TH2F(TString::Format("h2AllGoodLe"),
            TString::Format("h2AllGoodLe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTe = new TH2F(TString::Format("h2AllGoodTe"),
            TString::Format("h2AllGoodTe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTot = new TH2F(TString::Format("h2AllGoodTot"),
            TString::Format("h2AllGoodTot"),nTdc,0,nTdc,
            NTotBins, TotBinLow, TotBinHigh);
  hBarRateHV = new TH2F(TString::Format("hBarRateHV"),
	    TString::Format("hBarRateHV"), 250,1,250000,
	    50,600,800);

  h2TDCTOTvsLE = new TH2F(TString::Format("h2TDCTOTvsLE"),
            TString::Format("h2TDCTOTvsLE"),NTotBins,TotBinLow,TotBinHigh,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2CDetX1vsX2 = new TH2F(TString::Format("h2CDetX1vsX2"),
            TString::Format("h2CDetX1vsX2"),1000, -2.0, 2.0, 
            1000, -2.0, 2.0);
  h2TOTvsXDiff1 = new TH2F(TString::Format("h2TOTvsXDiff1"),
            TString::Format("h2TOTvsXDiff1"),NTotBins,TotBinLow,TotBinHigh,
            1000, -0.3, 0.3);
  h2TOTvsXDiff2 = new TH2F(TString::Format("h2TOTvsXDiff2"),
            TString::Format("h2TOTvsXDiff2"),NTotBins,TotBinLow,TotBinHigh,
            1000, -0.3, 0.3);
  h2LEvsXDiff1 = new TH2F(TString::Format("h2LEvsXDiff1"),
            TString::Format("h2LEvsXDiff1"),NTDCBins,TDCBinLow,TDCBinHigh,
            1000, -0.3, 0.3);
  h2LEvsXDiff2 = new TH2F(TString::Format("h2LEvsXDiff2"),
            TString::Format("h2LEvsXDiff2"),NTDCBins,TDCBinLow,TDCBinHigh,
            1000, -0.3, 0.3);
  // hRefRawLe = new TH1F(TString::Format("hRefRawLe"),
  //           TString::Format("hRefRawLe"),
  //           RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  // hRefRawTe = new TH1F(TString::Format("hRefRawTe"),
  //           TString::Format("hRefRawTe"),
  //           RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  // hRefRawTot = new TH1F(TString::Format("hRefRawTot"),
  //           TString::Format("hRefRawTot"),
  //           RefNTotBins, RefTotBinLow, RefTotBinHigh);
  // hRefRawPMT = new TH1F(TString::Format("hRefRawPMT"),
  //           TString::Format("hRefRawPMT"),
  //           32, 2688, 2720);
  hRefGoodLe = new TH1F(TString::Format("hRefGoodLe"),
            TString::Format("hRefGoodLe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTe = new TH1F(TString::Format("hRefGoodTe"),
            TString::Format("hRefGoodTe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTot = new TH1F(TString::Format("hRefGoodTot"),
            TString::Format("hRefGoodTot"),
            RefNTotBins, RefTotBinLow, RefTotBinHigh);
  hRefGoodPMT = new TH1F(TString::Format("hRefGoodPMT"),
            TString::Format("hRefGoodPMT"),
            2720, 0, 2720);
  
  
  for(Int_t bar=0; bar<(nTdc); bar++){
    // raw hits
    // leading edge
    hRawLe[bar] = new TH1F(TString::Format("hRawLe_Bar%d",bar),
 	    TString::Format("hRawLe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // trailing edge 
    hRawTe[bar] = new TH1F(TString::Format("hRawTe_Bar%d",bar),
	    TString::Format("hRawTe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // tot 
    hRawTot[bar] = new TH1F(TString::Format("hRawTot_Bar%d",bar),
			     TString::Format("hRawTot_Bar%d",bar),
			     NTotBins, TotBinLow, TotBinHigh);
    // good hits
    // leading edge
    hGoodLe[bar] = new TH1F(TString::Format("hGoodLe_Bar%d",bar),
 	    TString::Format("hGoodLe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // trailing edge 
    hGoodTe[bar] = new TH1F(TString::Format("hGoodTe_Bar%d",bar),
	    TString::Format("hGoodTe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // tot 
    hGoodTot[bar] = new TH1F(TString::Format("hGoodTot_Bar%d",bar),
			     TString::Format("hGoodTot_Bar%d",bar),
			     NTotBins, TotBinLow, TotBinHigh);
  }// bar loop



  //========================================================= Get data from tree
  if (!T) {
  T = new TChain("T");

  int runnum = RunNumber1;

  vCDetPaddleRawTot.assign(2688, std::vector<double>{});
  vCDetPaddleCutTot.assign(2688, std::vector<double>{});
  // Per-bar storage for good leading-edge times (bar = GoodElID/16)
  vBarGoodLe.assign(NumPMTs, std::vector<double>{});
  vBarGoodLeECalT.assign(NumPMTs, std::vector<double>{});
  vAllGoodECalT.clear();
  //int onlySegment = -1; // set to >=0 to pick just one

  AddRunFilesToChain(T, REPLAYED_DIR.Data(), runnum, minSeg, maxSeg);
}

  TTreeReader reader(T);
  
  //Set TTreeReaders
  /* ----- Earm ----- */ 

  // ******CDet******
  // ----- CDet arrays -----
  TTreeReaderArray<double> TDCmult(reader, "earm.cdet.tdc_mult");

  TTreeReaderArray<double> RawElID   (reader, "earm.cdet.hits.TDCelemID");
  TTreeReaderArray<double> RawElLE   (reader, "earm.cdet.hits.t");
  TTreeReaderArray<double> RawElTE   (reader, "earm.cdet.hits.t_te");
  TTreeReaderArray<double> RawElTot  (reader, "earm.cdet.hits.t_tot");

  TTreeReaderArray<double> GoodElID  (reader, "earm.cdet.hit.pmtnum");
  TTreeReaderArray<double> GoodElLE  (reader, "earm.cdet.hit.tdc_le");
  TTreeReaderArray<double> GoodElTE  (reader, "earm.cdet.hit.tdc_te");
  TTreeReaderArray<double> GoodElTot (reader, "earm.cdet.hit.tdc_tot");

  TTreeReaderArray<double> GoodX     (reader, "earm.cdet.hit.xhit");
  TTreeReaderArray<double> GoodY     (reader, "earm.cdet.hit.yhit");
  TTreeReaderArray<double> GoodZ     (reader, "earm.cdet.hit.zhit");

  TTreeReaderArray<double> GoodCol   (reader, "earm.cdet.hit.row");
  TTreeReaderArray<double> GoodRow   (reader, "earm.cdet.hit.col");
  TTreeReaderArray<double> GoodLayer (reader, "earm.cdet.hit.layer");

  // ----- cdet scalars ----- 
  TTreeReaderValue<double> nhits        (reader, "earm.cdet.nhits");
  TTreeReaderValue<double> ngoodhits    (reader, "earm.cdet.ngoodhits");
  TTreeReaderValue<double> ngoodTDChits (reader, "earm.cdet.ngoodTDChits");

  //------ECal-------
  // Cluster arrays
  // TTreeReaderArray<double> ECal_clus_adctime      (reader, "earm.ecal.clus.adctime");
  // TTreeReaderArray<double> ECal_clus_again        (reader, "earm.ecal.clus.again");
  // TTreeReaderArray<double> ECal_clus_atimeblk     (reader, "earm.ecal.clus.atimeblk");
  // TTreeReaderArray<double> ECal_clus_col          (reader, "earm.ecal.clus.col");
  // TTreeReaderArray<double> ECal_clus_e            (reader, "earm.ecal.clus.e");
  // TTreeReaderArray<double> ECal_clus_eblk         (reader, "earm.ecal.clus.eblk");
  // TTreeReaderArray<double> ECal_clus_id           (reader, "earm.ecal.clus.id");
  // TTreeReaderArray<double> ECal_clus_nblk         (reader, "earm.ecal.clus.nblk");
  // TTreeReaderArray<double> ECal_clus_row          (reader, "earm.ecal.clus.row");
  // TTreeReaderArray<double> ECal_clus_x            (reader, "earm.ecal.clus.x");
  // TTreeReaderArray<double> ECal_clus_y            (reader, "earm.ecal.clus.y");

  // Cluster count (scalar)
  //TTreeReaderValue<double> ECal_nclus(reader, "earm.ECal.nclus");

  //event-level ECal branches
  TTreeReaderValue<double> ECalX       (reader, "earm.ecal.x");
  TTreeReaderValue<double> ECalY       (reader, "earm.ecal.y");
  TTreeReaderValue<double> ECalE       (reader, "earm.ecal.e");
  TTreeReaderValue<double> ECalAdcTime (reader, "earm.ecal.adctime");

  /* ----- SBS branches ----- 
    ------- comment out for now ---------
  // Scalars
  TTreeReaderValue<double> heep_dpp(reader, "heep.dpp");
  TTreeReaderValue<double> heep_dt_ADC(reader, "heep.dt_ADC");
  TTreeReaderValue<double> heep_ECalo(reader, "heep.ECalo");
  TTreeReaderValue<double> heep_eprime_eth(reader, "heep.eprime_eth");
  TTreeReaderValue<double> heep_dxECAL(reader, "heep.dxECAL");
  TTreeReaderValue<double> earm_ECal_x(reader, "earm.ECal.x");
  TTreeReaderValue<double> sbs_gemFPP_track_ntrack(reader, "sbs.gemFPP.track.ntrack");
  TTreeReaderValue<double> heep_dyECAL(reader, "heep.dyECAL");

  // Arrays
  TTreeReaderArray<double> sbs_tr_vz(reader, "sbs.tr.vz");
  TTreeReaderArray<double> sbs_gemFPP_track_sclose(reader, "sbs.gemFPP.track.sclose");
  TTreeReaderArray<double> sbs_gemFT_track_nhits(reader, "sbs.gemFT.track.nhits");
  TTreeReaderArray<double> sbs_gemFT_track_ngoodhits(reader, "sbs.gemFT.track.ngoodhits");
  */
  //========================================================= Check no of events
  

  //Should likely just have Nev = T->GetEntries();, so it just grabs all events ------- NOT WORKING --------
  Int_t Nev = T->GetEntries();
  cout << "N entries in tree is " << Nev << endl;
  Int_t NEventsAnalysis;// = Nev;
  if(nevents==-1) NEventsAnalysis = Nev;
  else NEventsAnalysis = nevents;
  cout << "Running analysis for " << NEventsAnalysis << " events" << endl;
  

  /*
  //==================================================== Create output root file
  // root file for viewing fits
  TString subfile; 
  subfile = TString::Format("gep5_replayed_nogems_%d_50k_events.root",RunNumber1);
  TString outrootfile = ANALYSED_DIR + "/RawTDC_" + subfile;
  TFile *f = new TFile(outrootfile, "RECREATE");
  */



  //================================================================= Event Loop
  // variables outside event loop
  Int_t EventCounter = 0;

  // DEBUG controls
  const bool DBG = false;        // master on/off
  const long DBG_ENTRY = -1;    // set to a specific tree entry, or -1 for all

  cout << "Starting Event Loop" << endl;

    int eff_denominator = 0;
    int eff_numerator_layer1 = 0;
    int eff_numerator_layer2 = 0;
    int eff_numerator = 0;

  // event loop start
  Int_t event = 0;
  while(reader.Next()){
    event++;
    event = event - 1;
    EventCounter++;
    // Only stop early if nevents > 0
    if (nevents > 0 && EventCounter > nevents) {
        break;
    }
    Int_t nh = *nhits;
    Int_t ngh = *ngoodhits;
    Int_t ngth = *ngoodTDChits;
    
    if (EventCounter % 1000 == 0) {
	cout << EventCounter << "/" << NEventsAnalysis << "/ Nhits = " << (Int_t)nh << endl;
    	for (Int_t nfill=0; nfill<nh; nfill++) {hnhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<ngh; nfill++) {hngoodhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<ngth; nfill++) {hngoodTDChits_ev->Fill(EventCounter);}
    }

    /* Fill ECal cluster vectors */ /////-------------- These need to get moved into the if statement after second pass, rawEventCounter>=1
    // ---- Per-event filling ----
    //v_ECal_clus_adctime.push_back(copyArray(ECal_clus_adctime));
    //v_ECal_clus_again.push_back(copyArray(ECal_clus_again));
    // v_ECal_clus_atimeblk.push_back(copyArray(ECal_clus_atimeblk));

    //v_ECal_clus_col.push_back(copyArray(ECal_clus_col));
    //v_ECal_clus_row.push_back(copyArray(ECal_clus_row));

    //v_ECal_clus_e.push_back(copyArray(ECal_clus_e));
    //v_ECal_clus_eblk.push_back(copyArray(ECal_clus_eblk));

    // v_ECal_clus_id.push_back(copyArray(ECal_clus_id));

    // v_ECal_clus_nblk.push_back(copyArray(ECal_clus_nblk));
    
    // v_ECal_clus_x.push_back(copyArray(ECal_clus_x));
    // v_ECal_clus_y.push_back(copyArray(ECal_clus_y));

    // Event-level scalars
    //v_ECal_nclus.push_back(*ECal_nclus);
    
    

    bool good_elastic;
    if (elastic == 0) good_elastic = true; //abs(*heep_dt_ADC)<10 && abs(sbs_tr_vz[0]+0.1)<0.18 && *heep_ECalo/(*heep_eprime_eth) > 0.7 && abs(*heep_dxECAL - 0.01 + 0.025 * (*earm_ECal_x)) < 0.05 && *sbs_gemFPP_track_ntrack > 0 && abs(*heep_dyECAL - 0.01) < 0.06 && sbs_gemFPP_track_sclose[0] < 0.01 && (sbs_gemFT_track_nhits[0] > 4 || sbs_gemFT_track_ngoodhits[0] > 2);
    //currently not using full replays, so elastic cuts dont work, just assume all elastic
    else if (elastic == 1) good_elastic = true; //incase one does not want to use the elastic cut
    if (good_elastic){

      //fill vectors we wish to make cuts on for selecting elastics
      // Scalars — push_back the dereferenced values
      /* Comment out heep and gems vectors, do not use root files witht them yet
      vheep_dpp.push_back(*heep_dpp);
      vheep_dt_ADC.push_back(*heep_dt_ADC);
      vheep_ECalo.push_back(*heep_ECalo);
      vheep_eprime_eth.push_back(*heep_eprime_eth);
      vheep_dxECAL.push_back(*heep_dxECAL);
      vearm_ECal_x.push_back(*earm_ECal_x);
      vsbs_gemFPP_track_ntrack.push_back(*sbs_gemFPP_track_ntrack);
      vheep_dyECAL.push_back(*heep_dyECAL);

      //2D arrays
      vsbs_tr_vz.push_back(fill2D(sbs_tr_vz));
      vsbs_gemFPP_track_sclose.push_back(fill2D(sbs_gemFPP_track_sclose));
      vsbs_gemFT_track_nhits.push_back(fill2D(sbs_gemFT_track_nhits));
      vsbs_gemFT_track_ngoodhits.push_back(fill2D(sbs_gemFT_track_ngoodhits));
      */
            // Per-entry reference time (ns) for aligning with GOOD-event vectors
      double thisEvent_refRawLe_ns = std::nan("");

// First pass through hits:  purpose is to get reference LE TDC Value for this event
      
      double event_ref_tdc = 0.0;
      double ref_int = 0;
      double ref_corr = 0;
      for(Int_t el=0; el<RawElID.GetSize(); el++) {
        if ((Int_t)RawElID[el] == 2696) {  // only look at ref PMT 
          bool good_ref_le_time = RawElLE[el] > 0.0/TDC_calib_to_ns && RawElLE[el] <= 100.0/TDC_calib_to_ns;
          bool good_ref_tot = RawElTot[el] >= 0.0/TDC_calib_to_ns && RawElTot[el] <= 200.0/TDC_calib_to_ns;
          bool good_ref_event = good_ref_le_time && good_ref_tot;
          if ( good_ref_event ) {

            //if (RawElID[el] > 2687) {
            //	cout << "el = " << el << " Raw ID = " << RawElID[el] << " raw le = " << 
          //	RawElLE[el] << " raw te = " << RawElTE[el] << " raw tot = " << 
          //	RawElTot[el] << " CDet X = " << GoodX[el] << " ECal X = " << ECalX << endl;
            //}
            if ( !check_bad(RawElID[el],suppress_bad) ) {
            //cout << " el = " << el << endl;
            //cout << " tdc = " << RawElLE[el]*TDC_calib_to_ns << endl;
              if ( (Int_t)RawElID[el] == 2696 && (Int_t)RawElLE[el]>0 && (Int_t)RawElTot[el]>0 ) {
                //cout << " Ref  ID = " << (Int_t)RawElID[el] << " el = " << el << "    LE = " << RawElLE[el]*TDC_calib_to_ns 
                //		<< "    TE = " << RawElTE[el]*TDC_calib_to_ns << "    ToT = " << RawElTot[el]*TDC_calib_to_ns << endl;
                
                thisEvent_refRawLe_ns = RawElLE[el] * TDC_calib_to_ns;
                vRefRawLe.push_back(thisEvent_refRawLe_ns);
                vRefRawTe.push_back(RawElTE[el] * TDC_calib_to_ns);
                vRefRawTot.push_back(RawElTot[el] * TDC_calib_to_ns);
                vRefRawPMT.push_back((int)RawElID[el]);
                //hRefRawLe->Fill(RawElLE[el]*TDC_calib_to_ns);
                //hRefRawTe->Fill(RawElTE[el]*TDC_calib_to_ns);
                //hRefRawTot->Fill(RawElTot[el]*TDC_calib_to_ns);
                //hRefRawPMT->Fill(RawElID[el]);

                event_ref_tdc = RawElLE[el]*TDC_calib_to_ns - 48 ;
                ref_int = std::floor(event_ref_tdc);
                ref_corr = event_ref_tdc - ref_int;

              }
            }
          }
        }
      }// end ref TDC loop

      // second pass: fill raw CDet TDC histos
      std::vector<double> thisEvent_LE;
      std::vector<double> thisEvent_TE;
      std::vector<double> thisEvent_TOT;
      std::vector<int> thisEvent_ID;

      std::vector<double> thisEvent_CDetX;
      std::vector<double> thisEvent_CDetY;
      std::vector<double> thisEvent_CDetZ;
      std::vector<CDetHit> eventHits;
      rateEvTrack++;

      // Build lookup from PMT id -> index in Good* arrays for this TTree entry
      std::unordered_map<int,int> goodIdx;
      goodIdx.reserve(GoodElID.GetSize());
      for (int ig = 0; ig < (int)GoodElID.GetSize(); ig++) {
        goodIdx[(int)GoodElID[ig]] = ig;
      }

      int rawEventCounter = 0;
      for(Int_t el=0; el<RawElID.GetSize(); el++){

      const int raw_pmt = (int)RawElID[el];
      auto itGood = goodIdx.find(raw_pmt);
      const bool hasGood = (itGood != goodIdx.end());
      const int ig = hasGood ? itGood->second : -1;
      const double gx = hasGood ? GoodX[ig] : 1.0e9;
      const double gy = hasGood ? GoodY[ig] : 1.0e9;
      const double gz = hasGood ? GoodZ[ig] : 1.0e9;

      bool good_raw_le_time = RawElLE[el] >= LeMin/TDC_calib_to_ns && RawElLE[el] <= LeMax/TDC_calib_to_ns;
      bool good_raw_tot = RawElTot[el] >= TotMin/TDC_calib_to_ns && RawElTot[el] <= TotMax/TDC_calib_to_ns;
      bool good_mult = TDCmult[el] < TDCmult_cut;
      bool good_CDet_X = hasGood && (fabs(gx) < xcut);
      // bool good_ECal_diff_x = (GoodX[el]-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) <= XDiffCut && 
      //     (GoodX[el]-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) >= -1.0*XDiffCut;
      // bool good_ECal_diff_y = (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) <= 1.2*CDet_y_half_length && 
      //     (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) >= -1.2*CDet_y_half_length;

      bool good_raw_event = good_raw_le_time && good_raw_tot && good_mult && good_CDet_X ;//&& good_ECal_diff_x && good_ECal_diff_y;

      //if ((Int_t)RawElID[el] > 1000) cout << "el = " << el << " Hit ID = " << (Int_t)RawElID[el] << "    TDC = " << RawElLE[el]*TDC_calib_to_ns << endl;
      //cout << "Raw ID = " << RawElID[el] << " raw le = " << RawElLE[el] << " raw te = " << RawElTE[el] << " raw tot = " << RawElTot[el] << endl;
      if ( good_raw_event ) {
        rawEventCounter++;
        int idx = RawElID[el];
        if (0 <= idx && idx < 2688) {
          double le_ns = RawElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr(idx);
          double tot_ns = RawElTot[el]*TDC_calib_to_ns;
          double te_ns = RawElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr(idx);
          rawRate[idx]++;
          vCDetPaddleRawTot[idx].push_back(tot_ns);
          eventHits.push_back({idx, le_ns, tot_ns, te_ns});
        } //getting rates and tot for pixels
        if ( !check_bad(RawElID[el],suppress_bad) ) {
        //cout << " el = " << el << endl;
        //cout << " tdc = " << RawElLE[el]*TDC_calib_to_ns << endl;
          if ( (Int_t)RawElID[el] < 2688 ) {
            //if ((Int_t)RawElID[el] > nTdc) cout << " CDet ID = " << (Int_t)RawElID[el] << "    TDC = " << RawElLE[el]*TDC_calib_to_ns << endl;
            
            //fill this events vectors
            thisEvent_LE.push_back(RawElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)RawElID[el])); 
            thisEvent_TE.push_back(RawElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)RawElID[el]));
            thisEvent_TOT.push_back(RawElTot[el]*TDC_calib_to_ns); //- event_ref_tdc);
            thisEvent_ID.push_back((Int_t)RawElID[el]);
          
            //fill all hits vectors
            vAllRawLe.push_back(RawElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)RawElID[el]));
            vAllRawTe.push_back(RawElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)RawElID[el]));
            vAllRawTot.push_back(RawElTot[el]*TDC_calib_to_ns);
            vAllRawPMT.push_back(RawElID[el]);
            vAllRawBar.push_back((Int_t)(RawElID[el]/16));

            thisEvent_CDetX.push_back(gx);
            thisEvent_CDetY.push_back(gy);
            thisEvent_CDetZ.push_back(gz);
            //if (fabs(GoodX[el]) == 999 && GoodZ[el] != 999){
            if (DBG && (DBG_ENTRY < 0 || reader.GetCurrentEntry() == DBG_ENTRY) && rawEventCounter<20) {
            std::cout << "event = " << rawEventCounter << " " << "cdetX = " << gx << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetY = " << gy << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetZ = " << gz << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetID = " << (Int_t)RawElID[el] << std::endl;
            std::cout << " " <<std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetLE = " << RawElLE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTE = " << RawElTE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTot = " << RawElTot[el]*TDC_calib_to_ns << std::endl;
            std::cout << "-------------------- " <<std::endl;
            }
            /*}
            if (fabs(GoodX[el]) == 999 && GoodZ[el] != -999){
            std::cout << "event = " << rawEventCounter << " " << "cdetX = " << gx << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetZ = " << gz << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetID = " << (Int_t)RawElID[el] << std::endl;
            std::cout << " " <<std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetLE = " << RawElLE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTE = " << RawElTE[el]*TDC_calib_to_ns - event_ref_tdc << std::endl;
            std::cout << "event = " << rawEventCounter << " " << "cdetTot = " << RawElTot[el]*TDC_calib_to_ns << std::endl;
            }*/
            /* Comment out histograms for now
            hRawLe[(Int_t)RawElID[el]]->Fill(RawElLE[el]*TDC_calib_to_ns-event_ref_tdc);
            hRawTe[(Int_t)RawElID[el]]->Fill(RawElTE[el]*TDC_calib_to_ns-event_ref_tdc);
            hRawTot[(Int_t)RawElID[el]]->Fill(RawElTot[el]*TDC_calib_to_ns);
            hAllRawLe->Fill(RawElLE[el]*TDC_calib_to_ns-event_ref_tdc);
            hAllRawTe->Fill(RawElTE[el]*TDC_calib_to_ns-event_ref_tdc);
            hAllRawTot->Fill(RawElTot[el]*TDC_calib_to_ns);
            hAllRawPMT->Fill(RawElID[el]);
            hAllRawBar->Fill((Int_t)(RawElID[el]/16));



            h2d_RawLE->Fill(RawElLE[el]*TDC_calib_to_ns-event_ref_tdc, (Int_t)RawElID[el]);
            h2d_RawTE->Fill(RawElTE[el]*TDC_calib_to_ns-event_ref_tdc, (Int_t)RawElID[el]);
            h2d_RawTot->Fill(RawElTot[el]*TDC_calib_to_ns, (Int_t)RawElID[el]);
            */
          }
        }
      }
    }// all raw tdc hit loop
    if (rawEventCounter >= 1){
      vRawLe.push_back(thisEvent_LE);
      vRawTe.push_back(thisEvent_TE);
      vRawTot.push_back(thisEvent_TOT);
      vRawID.push_back(thisEvent_ID);
      vCDetX.push_back(thisEvent_CDetX);
      vCDetY.push_back(thisEvent_CDetY);
      vCDetZ.push_back(thisEvent_CDetZ);
      v_ECalX.push_back(*ECalX);
      v_ECalY.push_back(*ECalY);
      v_ECalE.push_back(*ECalE);
      v_ECalAdcTime.push_back(*ECalAdcTime);
      vEventHits.push_back(eventHits);
    }
    //check nadjacent pairs for each event 
    auto is_unused = [&](int id) -> bool {
      return kUnusedCDetPixels.count(id) != 0;
    };

    auto is_adjacent_with_skip = [&](int a, int b) -> bool {
      if (b == a + 1) {
        // adjacent normally, but only if that neighbor isn't unused
        return !is_unused(b);
      }
      if (b == a + 2) {
        // treat as adjacent if the in-between pixel is unused
        return is_unused(a + 1) && !is_unused(b);
      }
      return false;
    };
    if (rawEventCounter >= 1) {
      std::vector<int> ids;
      const auto& currentEvent = vEventHits.back();
      ids.reserve(currentEvent.size());

      for (const auto& hit : currentEvent) {
        int id = hit.id;
        // if (260 <= id && id <= 270) ids.push_back(id);
        ids.push_back(id);
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

      int nAdjacentHits = 0;
      for (size_t i = 0; i + 1 < ids.size(); i++) {
        if (is_adjacent_with_skip(ids[i], ids[i + 1])) {
          nAdjacentHits++;
        }
      }

      vNumRawAdjacentHits.push_back(nAdjacentHits);
    }

    // Third pass:  Get layer occupancies
    
    int nhitsc1 = 0;
    int nhitsc2 = 0;
    int ngoodhitsc1 = 0;
    int ngoodhitsc2 = 0;
    int ngoodTDChitsc1 = 0;
    int ngoodTDChitsc2 = 0;
    for (int j=0; j<nTdc; j++) {
      nhits_paddles[j]=0;
      ngoodhits_paddles[j]=0;
      ngoodTDChits_paddles[j]=0;
    }
    npaddles=0;
    ngoodpaddles=0;
    ngoodTDCpaddles=0;
    
    for(Int_t el=0; el<GoodElID.GetSize(); el++){
      int sbselemid = (Int_t)GoodElID[el];
      int sbsrown = sbselemid%672;
      int sbscoln = sbselemid/672;
      //int sbsrown = (Int_t)GoodRow[el];
      //int sbscoln = (Int_t)GoodCol[el];
      int mylayern = sbscoln/2;
      int mypaddlen = sbscoln*672 + sbsrown;

      if (mylayern == 0) {
        nhitsc1++;
      } else {
        nhitsc2++;
      }
      nhits_paddles[mypaddlen]++;

      bool good_ECal_reconstruction = *ECalY > -1.2 && *ECalY < 1.2 &&
                                      *ECalX > -1.5 && *ECalX < 1.5 &&
                                      *ECalX != 0.00 && *ECalY != 0.00 ;
      bool good_le_time = GoodElLE[el] >= LeMin/TDC_calib_to_ns && GoodElLE[el] <= LeMax/TDC_calib_to_ns;
      bool good_tot = GoodElTot[el] >= TotMin/TDC_calib_to_ns && GoodElTot[el] <= TotMax/TDC_calib_to_ns;
      bool good_hit_mult = TDCmult[el] < TDCmult_cut;
      bool good_CDet_X = GoodX[el] < xcut;
      bool good_ECal_diff_x = (GoodX[el]-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) <= XDiffCut && 
          (GoodX[el]-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) >= -1.0*XDiffCut;
      bool good_ECal_diff_y = (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) <= 1.2*CDet_y_half_length && 
          (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) >= -1.2*CDet_y_half_length;
      
      
      bool good_CDet_event = good_ECal_reconstruction && good_ECal_diff_x && good_ECal_diff_y && good_le_time && good_tot && good_hit_mult && good_CDet_X;


      if (good_CDet_event) {
        if ( !check_bad(GoodElID[el], suppress_bad) ) {
          if ( (Int_t)GoodElID[el]%NumSidesTotal < NumCDetPaddlesPerSide )  {

            //cout << "Hit number " << el << ":    Paddle = " << mypaddlen << " Row = " << sbsrown  << " Col = " << sbscoln  << " hits = " << ngoodTDChits_paddles[mypaddlen] << endl;
            //cout << "el = " << el << " Good ID = " << GoodElID[el] << " Good le = " << 
        //	GoodElLE[el] << " Good te = " << GoodElTE[el] << " Good tot = " << 
        //	GoodElTot[el] << " CDet X = " << GoodX[el] << " ECal X = " << ECalX << endl;
            if (mylayern == 0) {
                ngoodhitsc1++;
            } else {
                ngoodhitsc2++;
            }
            ngoodhits_paddles[mypaddlen]++;
          }
        }
      }
    }
    for (int j=0; j<nTdc; j++) {
      if (nhits_paddles[j] > 0) {
        npaddles++;
        //cout << "Paddle = " << j <<  "  nhits = " << ngoodTDChits_paddles[j] << endl;
      }
    }
    for (int j=0; j<nTdc; j++) {
      if (ngoodhits_paddles[j] > 0) {
        ngoodpaddles++;
        //cout << "Paddle = " << j <<  "  nhits = " << ngoodTDChits_paddles[j] << endl;
      }
    }
        //cout << "event " << event << endl;
        //cout << "Number of good layer 1 hits: " << ngoodTDChitsc1 << endl;
        //cout << "Number of good layer 2 hits: " << ngoodTDChitsc2 << endl;
        //cout << "Layer 1 Hit Cut " << nhitcutlow1 << " " << nhitcuthigh1 << endl;
        //cout << "Layer 2 Hit Cut " << nhitcutlow2 << " " << nhitcuthigh2 << endl;
    vnpaddles.push_back(npaddles);
    vngoodpaddles.push_back(ngoodpaddles);
    //hnpaddles->Fill(npaddles);
    //hngoodpaddles->Fill(ngoodpaddles);

    // Fourth pass:  use layer occupancies to apply additional cuts
    
    //before loop temp vectors to fill into vGoodLE, etc.
    std::vector<double> thisEvent_GoodLE;
    std::vector<double> thisEvent_GoodTE;
    std::vector<double> thisEvent_GoodTOT;
    std::vector<int> thisEvent_GoodID;

    std::vector<double> thisEvent_GoodX;
    std::vector<double> thisEvent_GoodY;
    std::vector<double> thisEvent_GoodZ;
    std::vector<int> thisEvent_GoodLayer;
    std::vector<int> thisEvent_GoodCol;
    std::vector<CDetHit> goodEventHits;

    // --- NEW: per-event best (smallest |x-diff|) hit in each layer
    double bestAbsXDiff[2] = {1e99, 1e99};
    double bestXCDet[2]    = {0.0, 0.0};
    double bestXECalProj[2]= {0.0, 0.0};
    bool   foundBest[2]    = {false, false};

    int CDetPassedBoolCount = 0;

    for(Int_t el=0; el<GoodElID.GetSize(); el++){
      bool goodhit_ECal_reconstruction = *ECalY > -1.2 && *ECalY < 1.2 &&
                                         *ECalX > -1.5 && *ECalX < 1.5 &&
                                         *ECalX != 0.00 && *ECalY != 0.00;
      bool goodhit_le_time = GoodElLE[el] >= LeMin/TDC_calib_to_ns && GoodElLE[el] <= LeMax/TDC_calib_to_ns;
      bool goodhit_tot = GoodElTot[el] >= TotMin/TDC_calib_to_ns && GoodElTot[el] <= TotMax/TDC_calib_to_ns;
      bool goodhit_hit_mult = TDCmult[el] < TDCmult_cut;
      bool goodhit_CDet_X = GoodX[el] < xcut;
      bool goodhit_low = ngoodhitsc1 >= nhitcutlow1  && ngoodhitsc2 >= nhitcutlow2;
      bool goodhit_high  = ngoodhitsc1 <= nhitcuthigh1 && ngoodhitsc2 <= nhitcuthigh2; 
      bool goodhit_ECal_diff_x = (GoodX[el]-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) <= XDiffCut && 
          (GoodX[el]-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset) >= -1.0*XDiffCut;
      bool goodhit_ECal_diff_y = (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) <= 1.2*CDet_y_half_length && 
          (GoodY[el]-((*ECalY)*(GoodZ[el])/ECal_dist)-YOffset) >= -1.2*CDet_y_half_length;
      bool goodhit_CDet_event = goodhit_ECal_reconstruction && goodhit_ECal_diff_x && goodhit_ECal_diff_y && goodhit_le_time && goodhit_tot 
        && goodhit_hit_mult && goodhit_CDet_X && goodhit_low && goodhit_high;

      if (goodhit_CDet_event) {
        // GoodX[el]-((*ECalX)*(GoodZ[el])/ECal_dist)-XOffset
        // std::cout << " gx = " << GoodX[el] << " & ECalX_Proj = " << (*ECalX)*GoodZ[el]/ECal_dist - XOffset <<std::endl;
        CDetPassedBoolCount++;
        int idx = GoodElID[el];
        if (0 <= idx && idx < 2688) {
          double le_ns = GoodElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr(idx);
          double tot_ns = GoodElTot[el]*TDC_calib_to_ns;
          double te_ns = GoodElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr(idx);
          goodEventHits.push_back({idx, le_ns, tot_ns, te_ns});
        } //getting rates and tot for pixels
        if ( !check_bad(GoodElID[el], suppress_bad) ) {
          if ( (Int_t)GoodElID[el]%NumSidesTotal < NumCDetPaddlesPerSide )  {
                //cout << "event " << event << endl;
            //cout << "el = " << el << " Good ID = " << GoodElID[el] << " Good le = " << 
          //GoodElLE[el] << " Good te = " << GoodElTE[el] << " Good tot = " << 
          //GoodElTot[el] << " CDet X = " << GoodX[el] << " ECal X = " << ECalX << endl;

            //cout << "Filling good timing histos ... " << ngoodTDChitsc1 << " " << endl;
            
            //std::cout << "Layer = " << (Int_t)GoodLayer[el] << " Side = " << (Int_t)GoodCol[el] << std::endl;
            
            int sbselem = (Int_t)GoodElID[el];
            int sbsrow = sbselem%672;
            int sbscol = sbselem/672;
            //int sbscol = (Int_t)GoodCol[el];
            //int sbsrow = (Int_t)GoodRow[el];
            int myside = sbscol%2;
            int mylayer = sbscol/2;
            int mypaddle = sbscol*672 + sbsrow;

            if (mylayer == 0) {
              ngoodTDChitsc1++;
              eff_numerator_layer1++;
            } 
            else {
              ngoodTDChitsc2++;
              eff_numerator_layer2++;
            }

            ngoodTDChits_paddles[mypaddle]++;
            
            if ( (layer_choice == 1 && mylayer == 0) || (layer_choice == 2 && mylayer == 1) || 
            (layer_choice == 3 && ngoodhitsc1>=1 && ngoodhitsc2 >= 1) ) {
              eff_numerator++;

              thisEvent_GoodLE.push_back(GoodElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)GoodElID[el]));
              thisEvent_GoodTE.push_back(GoodElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)GoodElID[el]));
              thisEvent_GoodTOT.push_back(GoodElTot[el]*TDC_calib_to_ns);
              thisEvent_GoodID.push_back((Int_t)GoodElID[el]);

              // hGoodLe[(Int_t)GoodElID[el]]->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hGoodTe[(Int_t)GoodElID[el]]->Fill(GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hGoodTot[(Int_t)GoodElID[el]]->Fill(GoodElTot[el]*TDC_calib_to_ns);

              double t_ECal_event = *ECalAdcTime;
              double tLE_bar = GoodElLE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)GoodElID[el]);
              double tTE_bar = GoodElTE[el]*TDC_calib_to_ns - event_ref_tdc + GetBarToffsetCorr((Int_t)GoodElID[el]);
              vAllGoodLe.push_back(tLE_bar);
              vAllGoodTe.push_back(tTE_bar);
              vAllGoodTot.push_back(GoodElTot[el]*TDC_calib_to_ns);
              vAllGoodPMT.push_back(GoodElID[el]);
              vAllGoodBar.push_back((Int_t)(GoodElID[el]/16));
              vAllGoodECalT.push_back(t_ECal_event);
              {
                const int bar = (int)(GoodElID[el] / 16);
                if (0 <= bar && bar < NumPMTs) {
                  vBarGoodLe[bar].push_back(tLE_bar);
                  vBarGoodLeECalT[bar].push_back(t_ECal_event);
                }
              }

              // hAllGoodLe->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hAllGoodTe->Fill(GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // hAllGoodTot->Fill(GoodElTot[el]*TDC_calib_to_ns);
              // hAllGoodPMT->Fill(GoodElID[el]);
              // hAllGoodBar->Fill((Int_t)(GoodElID[el]/16));

              // h2AllGoodLe->Fill(GoodElID[el],GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // h2AllGoodTe->Fill(GoodElID[el],GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              // h2AllGoodTot->Fill(GoodElID[el],GoodElTot[el]*TDC_calib_to_ns);

              // h2TDCTOTvsLE->Fill(GoodElTot[el]*TDC_calib_to_ns,GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
              
              vhitCDetPMT.push_back((Int_t)GoodElID[el]);
              //hHitPMT->Fill((Int_t)GoodElID[el]);
              vRow.push_back((Int_t)GoodRow[el]);
              //hRow->Fill((Int_t)GoodRow[el]);
            }
        
            if (myside == 0) {
              if (mylayer == 0) {
                vRowLayer1Side1.push_back((Int_t)GoodRow[el]);
                //hRowLayer1Side1->Fill((Int_t)GoodRow[el]);
              } 
              else {
                vRowLayer2Side1.push_back((Int_t)GoodRow[el]);
                //hRowLayer2Side1->Fill((Int_t)GoodRow[el]);
              }
            }
            else {
              if(mylayer == 0) {
                vRowLayer1Side2.push_back((Int_t)GoodRow[el]);
                //hRowLayer1Side2->Fill((Int_t)GoodRow[el]);
              } 
              else {
                vRowLayer2Side2.push_back((Int_t)GoodRow[el]);
                //hRowLayer2Side2->Fill((Int_t)GoodRow[el]);
              }
            }	
            thisEvent_GoodCol.push_back(myside);
            //hCol->Fill(myside);
            thisEvent_GoodLayer.push_back(mylayer);
            //hLayer->Fill(mylayer);

            thisEvent_GoodX.push_back(GoodX[el]);
            thisEvent_GoodY.push_back(GoodY[el]);
            thisEvent_GoodZ.push_back(GoodZ[el]);

            // --- NEW: compute projected ECal X at this hit's Z, and update per-layer best if this is smallest |x-diff|
            if (*ECalX != 0.00) {
                const double xECalProj = (*ECalX) * (GoodZ[el]) / ECal_dist;
                const double xdiff     = GoodX[el] - xECalProj;
                const double axdiff    = fabs(xdiff);

                if (axdiff < bestAbsXDiff[mylayer]) {
                    bestAbsXDiff[mylayer] = axdiff;
                    bestXCDet[mylayer]    = GoodX[el];
                    bestXECalProj[mylayer]= xECalProj;
                    foundBest[mylayer]    = true;
                }
            }

//------------------------------------------------------- replace hist below
             if (mylayer==0) { //layer 1 "good" histograms & higher level
               //i think we can remove these histograms from here, and put them in their own plot routine, they just need vectors for GoodX positions from CDet and ECal
               h2TOTvsXDiff1->Fill(GoodElTot[el]*TDC_calib_to_ns,GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
               h2LEvsXDiff1->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0,GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
               hHitXY1->Fill(GoodY[el],GoodX[el]);
               hXECalCDet1->Fill(GoodX[el],(*ECalX)*(GoodZ[el])/ECal_dist);
               hYECalCDet1->Fill(GoodY[el],(*ECalY)*(GoodZ[el])/ECal_dist);
               hXDiffECalCDet1->Fill(GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
               hXPlusECalCDet1->Fill(GoodX[el]+(*ECalX)*(GoodZ[el])/ECal_dist);
               hEECalCDet1->Fill(*ECalE,GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
             } 
             else { //layer 2
               h2TOTvsXDiff2->Fill(GoodElTot[el]*TDC_calib_to_ns,GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
               h2LEvsXDiff2->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0,GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
               hHitXY2->Fill(GoodY[el],GoodX[el]);
               hXECalCDet2->Fill(GoodX[el],(*ECalX)*(GoodZ[el])/ECal_dist);
               hYECalCDet2->Fill(GoodY[el],(*ECalY)*(GoodZ[el])/ECal_dist);
               hXDiffECalCDet2->Fill(GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
               hXPlusECalCDet2->Fill(GoodX[el]+(*ECalX)*(GoodZ[el])/ECal_dist);
               hEECalCDet2->Fill(*ECalE,GoodX[el]-(*ECalX)*(GoodZ[el])/ECal_dist);
             }


          } 
          else {
            if (GoodElID[el]==2696){
              vRefGoodLe.push_back(GoodElLE[el]*TDC_calib_to_ns);
              vRefGoodTe.push_back(GoodElTE[el]*TDC_calib_to_ns);
              vRefGoodTot.push_back(GoodElTot[el]*TDC_calib_to_ns);
              vRefGoodPMT.push_back(GoodElID[el]);
            }

            // hRefGoodLe->Fill(GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
            // hRefGoodTe->Fill(GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
            // hRefGoodTot->Fill(GoodElTot[el]*TDC_calib_to_ns);
            // hRefGoodPMT->Fill(GoodElID[el]);
          }
        }
      }
    }// all good tdc hit loop

    // --- NEW: fill ONLY the single best-matching hit per layer for this event
    if (foundBest[0]) hXECalCDet1_min->Fill(bestXCDet[0], bestXECalProj[0]);
    if (foundBest[1]) hXECalCDet2_min->Fill(bestXCDet[1], bestXECalProj[1]);

    if (CDetPassedBoolCount >= 1){
      v_GoodECalX.push_back(*ECalX);
      v_GoodECalY.push_back(*ECalY);
      v_GoodECalE.push_back(*ECalE);
      v_GoodECalAdcTime.push_back(*ECalAdcTime);

      vGoodCol.push_back(thisEvent_GoodCol);
      //hCol->Fill(myside);
      vGoodLayer.push_back(thisEvent_GoodLayer);
      //hLayer->Fill(mylayer);

      vCDetGoodX.push_back(thisEvent_GoodX);
      vCDetGoodY.push_back(thisEvent_GoodY);
      vCDetGoodZ.push_back(thisEvent_GoodZ);
      vTreeEntry.push_back(reader.GetCurrentEntry());
      vGoodRefRawLe.push_back(thisEvent_refRawLe_ns);

      vGoodLe.push_back(thisEvent_GoodLE);
      vGoodTe.push_back(thisEvent_GoodTE);
      vGoodTot.push_back(thisEvent_GoodTOT);
      vGoodID.push_back(thisEvent_GoodID);

      vGoodEventHits.push_back(goodEventHits);
    }
    if (CDetPassedBoolCount >= 1) {
      std::vector<int> ids;
      const auto& currentEvent = vGoodEventHits.back();
      ids.reserve(currentEvent.size());

      for (const auto& hit : currentEvent) {
        int id = hit.id;
        // if (260 <= id && id <= 270) ids.push_back(id);
        ids.push_back(id);
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

      int nAdjacentHits = 0;
      for (size_t i = 0; i + 1 < ids.size(); i++) {
        if (is_adjacent_with_skip(ids[i], ids[i + 1])) {
          nAdjacentHits++;
        }
      }

      vNumGoodAdjacentHits.push_back(nAdjacentHits);
    }//end adjacent check

    if (*ECalX != 0.00 && *ECalY != 0.00) {//double check this later, probably want to fill vectors with ECal hit position
      eff_denominator++;
      hXYECal->Fill(*ECalY,*ECalX);
      hXECal->Fill(*ECalX);
      hYECal->Fill(*ECalY);
      hEECal->Fill(*ECalE);
    };
        
    // vnhits1.push_back(nhitsc1);
    // vngoodhits1.push_back(ngoodhitsc1);
    // vngoodTDChits1.push_back(ngoodTDChitsc1);
    // vnhits2.push_back(nhitsc2);
    // vngoodhits2.push_back(ngoodhitsc2);
    // vngoodTDChits2.push_back(ngoodTDChitsc2);
    //
    hnhits1->Fill(nhitsc1);
    hngoodhits1->Fill(ngoodhitsc1);
    hngoodTDChits1->Fill(ngoodTDChitsc1);
    hnhits2->Fill(nhitsc2);
    hngoodhits2->Fill(ngoodhitsc2);
    hngoodTDChits2->Fill(ngoodTDChitsc2);

    for (int j=0; j<nTdc; j++) {
      if (ngoodTDChits_paddles[j] > 0) {
        ngoodTDCpaddles++;
        //cout << "Paddle = " << j <<  "  nhits = " << ngoodTDChits_paddles[j] << endl;
      }
    }
    // vngoodTDCpaddles.push_back(ngoodTDCpaddles);
        hngoodTDCpaddles->Fill(ngoodTDCpaddles);

        //cout << "Element loop: " << NdataMult << endl;
    for(Int_t tdc=0; tdc<TDCmult.GetSize(); tdc++){
      if (!check_bad(RawElID[tdc],suppress_bad)) {
        hMultiplicity->Fill(TDCmult[tdc]);
        hMultiplicityL[(Int_t)RawElID[tdc]]->Fill(TDCmult[tdc]);
        if( TDCmult[tdc] != 0 ){
          h2d_Mult->Fill(TDCmult[tdc], (Int_t)RawElID[tdc] );
        }
      }
    }// element loop
    }//good elastic bool
  }//event loop 
  std::cout << "nevents = " << rateEvTrack << std::endl;
  for (int i = 0; i < 2688; i++){
    chanRates[i] = (double)rawRate[i] / (rateEvTrack);
    //std::cout << "triggered Rate in Pixel " << 417 + i << " = " << chanRates << " & with time window Rate = " << chanRates / winWidth <<std::endl;
  }

  //Second Pass over all events for tot_ave calc

  for (Int_t idx = 0; idx < 2688; idx++){
    Int_t ihits = vCDetPaddleRawTot[idx].size();
    double sumTot = 0;
    for (Int_t i = 0; i < ihits; i++){
      sumTot += vCDetPaddleRawTot[idx][i];
    }
    ave_tot[idx] = sumTot / ihits;
  }

  for (Int_t idx = 0; idx < 2688; idx++){
    Int_t ihits = vCDetPaddleRawTot[idx].size();
    for (Int_t i = 0; i < ihits; i++){
      double hit_tot = vCDetPaddleRawTot[idx][i];
      if (hit_tot > ave_tot[idx]){
        vCDetPaddleCutTot[idx].push_back(hit_tot);
        cutRate[idx]++;
      }
    }
    cutChanRates[idx] = (double)cutRate[idx] / rateEvTrack;
  }

  std::cout << "Candidate Events = " << eff_denominator << std::endl;
  std::cout << "Layer 1 Events = " << eff_numerator_layer1 << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator_layer1/eff_denominator  << std::endl;
  std::cout << "Layer 2 Events = " << eff_numerator_layer2 << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator_layer2/eff_denominator <<  std::endl;
  std::cout << "One Good Layer Events = " << eff_numerator << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator/eff_denominator <<  std::endl;

/*
  for (Int_t b=0; b<NumCDetPaddles; b++) {
	//if (hRawLe[b]->GetEntries() > EventCounter/HotChannelRatio) {
    if (hRawLe[b]->GetEntries() > EventCounter/NumCDetPaddles*2*1000) {
      int myhotlayer = b/1344 + 1;
      int myhotside = (b%1344)/672 + 1;
      int myhotmodule = (b%672)/244 + 1;
      int myhotbar = (b%672)%224/16 + 1;
      int myhotpaddle = ((b%672)%224)%16 + 1;
      int mycable = b/16;
      //std::cout << "Hot PMT!! ID = " << b << "  layer = " << myhotlayer <<
      //"   side = " << myhotside << "   module = " << myhotmodule <<
      //"   bar = " << myhotbar << "   paddle_PMT = " << myhotpaddle << " CDet Cable = " << mycable << "   Entries = " << hRawLe[b]->GetEntries() << std::endl;
    }
  }
  
    /// Get rid of this whole chunk?
    //========================================================== Write histos
  for(Int_t b=0; b<nTdc; b++){
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawLe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawLe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawLe[b]->Write();
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawTe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTe[b]->Write();
    // hRawTe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTot[b]->GetXaxis()->SetTitle("tot (ns)");
    // hRawTot[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTot[b]->Write();
    // hMultiplicityL[b]->GetXaxis()->SetLabelSize(0.06);
    hMultiplicityL[b]->GetXaxis()->SetTitle("tdc ref hit mult");
    hMultiplicityL[b]->Write();
  }
  // hMultiplicity->GetXaxis()->SetLabelSize(0.06);
  hMultiplicity->GetXaxis()->SetTitle("tdc hit multiplicity");
  hMultiplicity->Write();
  // hHitPMT->GetXaxis()->SetLabelSize(0.06);
  hHitPMT->GetXaxis()->SetTitle("Bar ID of Left PMT Hit");
  hHitPMT->SetTitle("");
  hHitPMT->Write();

  // 2D histograms
  h2d_RawLE->GetXaxis()->SetTitle("TDC Leading Edge Time [ns]");
  h2d_RawLE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawLE->SetTitle("");
  h2d_RawLE->Write();
  h2d_RawTE->GetXaxis()->SetTitle("TDC Trailing Edge Time [ns]");
  h2d_RawTE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTE->SetTitle("");
  h2d_RawTE->Write();

  h2d_RawTot->GetXaxis()->SetTitle("TDC Time-over-threshold [ns]");
  h2d_RawTot->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTot->SetTitle("");
  h2d_RawTot->Write();

  h2d_Mult->GetXaxis()->SetTitle("TDC Multiplicity [ns]");
  h2d_Mult->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_Mult->SetTitle("");
  h2d_Mult->Write();
  */
  // Get HV values
  vector<string> HVfilenames = {"l1Left.dat", "l1Right.dat", "l2Left.dat", "l2Right.dat"};

  vector<vector<double>> data = readDataFromFiles(HVfilenames);
  std::vector<double> barRateContents = extractBinContents(hAllRawBar);
  Double_t xval,yval;
  for (int ii=0;ii<4;ii++) {
    for (int jj=0;jj<42;jj++) {
      xval = data[ii][jj];
            yval = barRateContents[ii*42+jj];
      //cout << "Contents:  " << xval << " "  << yval << endl;
      hBarRateHV->Fill(yval,-xval);
    }
  }

  //========================================================== Close output file
  //f->Close();




  //==================================================== ECal-time correction
  // Remove linear correlation between CDet time and ECal ADC time:
  //   <t_CDet> = p0 + p1*t_ECal
  // Then apply a global shift so that mean(corrected LE) = gTargetMeanLE.
  if (gUseECalTimeCorr) {
    double sumResLE = 0.0;
    long long nResLE = 0;

    const size_t NevCorr = std::min(vGoodLe.size(), v_GoodECalAdcTime.size());
    for (size_t ev = 0; ev < NevCorr; ++ev) {
      const double tE = v_GoodECalAdcTime[ev];
      const size_t Nh = std::min(vGoodLe[ev].size(), vGoodID[ev].size());
      for (size_t ih = 0; ih < Nh; ++ih) {
        const int bar = vGoodID[ev][ih] / 16;
        if (bar < 0 || bar >= NumPMTs) continue;
        const double tLE = vGoodLe[ev][ih]; // already includes bar offset correction
        const double res = tLE - (gECalFitP0 + gECalFitP1*tE);
        sumResLE += res;
        nResLE++;
      }
    }

    const double muResLE = (nResLE > 0) ? (sumResLE / (double)nResLE) : 0.0;
    gECalDeltaShift = gTargetMeanLE - muResLE;

    std::cout << "[CDet] ECal-time corr enabled. p0=" << gECalFitP0 << " p1=" << gECalFitP1
              << "  muResLE=" << muResLE << " ns  => delta=" << gECalDeltaShift << " ns\n";

    // --- Apply to ALL stored CDet LE/TE times (good-hit level vectors) ---
    // Convention: corrected time = (t_barcorr - (p0+p1*tE)) + delta
    for (size_t i = 0; i < vAllGoodLe.size() && i < vAllGoodECalT.size() && i < vAllGoodBar.size(); ++i) {
      const int bar = vAllGoodBar[i];
      if (bar < 0 || bar >= NumPMTs) continue;
      const double tE = vAllGoodECalT[i];
      vAllGoodLe[i] = (vAllGoodLe[i] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
      vAllGoodTe[i] = (vAllGoodTe[i] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
    }

    // per-bar LE storage used for bar histograms
    if (vBarGoodLe.size() == (size_t)NumPMTs && vBarGoodLeECalT.size() == (size_t)NumPMTs) {
      for (int bar = 0; bar < NumPMTs; ++bar) {
        const size_t Nh = std::min(vBarGoodLe[bar].size(), vBarGoodLeECalT[bar].size());
        for (size_t j = 0; j < Nh; ++j) {
          const double tE = vBarGoodLeECalT[bar][j];
          vBarGoodLe[bar][j] = (vBarGoodLe[bar][j] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
        }
      }
    }

    // event-level good-hit vectors (used by many timing comparison plots)
    for (size_t ev = 0; ev < NevCorr; ++ev) {
      const double tE = v_GoodECalAdcTime[ev];
      const size_t Nh = std::min(vGoodLe[ev].size(), vGoodID[ev].size());
      for (size_t ih = 0; ih < Nh; ++ih) {
        vGoodLe[ev][ih] = (vGoodLe[ev][ih] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
        vGoodTe[ev][ih] = (vGoodTe[ev][ih] - (gECalFitP0 + gECalFitP1*tE)) + gECalDeltaShift;
      }
    }
  }
  //================================================================== End Macro
}// end main

void plotNumAdjacent(int nbins = 50){
  TH1::AddDirectory(kFALSE);
  TH1D* hNumRawAdjacentHits = new TH1D("hNumRawAdjacentHits", "Number Raw Hits in Adjacent Pixels", nbins, 0, nbins);
  TH1D* hNumGoodAdjacentHits = new TH1D("hNumGoodAdjacentHits", "Number Raw Hits in Adjacent Pixels", nbins, 0, nbins);

  for (const auto& hit : vNumRawAdjacentHits){
    if (hit > 0) hNumRawAdjacentHits->Fill(hit);
  }
  for (const auto& hit : vNumGoodAdjacentHits){
    if (hit > 0) hNumGoodAdjacentHits->Fill(hit);
  }

  TCanvas* cNumAdjacentHits = new TCanvas("cNumAdjacentHits", "Number of Adjacent Hits", 900,700);
  cNumAdjacentHits->Divide(1,2);
  cNumAdjacentHits->cd(1);
  hNumRawAdjacentHits->Draw();
  cNumAdjacentHits->cd(2);
  hNumGoodAdjacentHits->Draw();

}

void plotAveTotPerPixel() {
  const int nPixels = ave_tot.size();

  TH1::AddDirectory(kFALSE);

  TH1D* hAveTot = new TH1D("hAveTot","Average TOT per CDet pixel;Pixel ID;Average TOT (ns)",nPixels, -0.5, nPixels - 0.5);

  for (int p = 0; p < nPixels; ++p) {
    if (ave_tot[p] > 0) {   // optional guard
      hAveTot->SetBinContent(p + 1, ave_tot[p]);
    }
  }

  TCanvas* c = new TCanvas("cAveTot", "Average TOT per Pixel", 1200, 500);
  hAveTot->Draw("HIST");
}

void plotSingleTot(int pixel_base = 0, bool raw = true, double width = 1, double totMin=1, double totMax=80){
  TH1::AddDirectory(kFALSE);
  if (pixel_base % 16 != 0) {
    Error("plotSingleTot", "pixel_base = %d is not a multiple of 16", pixel_base);
    return;
  }

  const int nPlots = 16;
  int TDCBinNum = (int)((totMax-totMin)/width);

  // Decode pixel → {layer, side, submodule, pmt, pixel}
  auto info = getLocation(pixel_base);

  int layer     = info[0] + 1;  // display as 1-based
  int side      = info[1];      // 0=L, 1=R
  int submodule = info[2] + 1;
  int bar       = info[3] + 1;

  TString sideStr = (side == 0) ? "L" : "R";

  TString canvasTitle = Form("Layer %d | %s | Module %d | Bar %d", layer, sideStr.Data(), submodule, bar);
  // Canvas with 4x4 pads
  TString cname = Form("cTot_%d", pixel_base);
  TCanvas* cTot = new TCanvas(cname, canvasTitle, 1200, 1000);
  cTot->Divide(4, 4, 0.001, 0.001);

  // Histogram array
  TH1D* hTot[nPlots];

  for (int i = 0; i < nPlots; i++) {

    int pixel = pixel_base + i;

    TString hname  = Form("hTot_pix%d", pixel);
    if (raw){
      TString htitle = Form("Pixel %d;TOT (ns);Counts", pixel);
      hTot[i] = new TH1D(hname, htitle, TDCBinNum, totMin, totMax);
    }
    if (!raw){
      TString htitle = Form("Pixel %d w/ TOT > %f;TOT (ns);Counts", pixel, ave_tot[pixel]);
      hTot[i] = new TH1D(hname, htitle, TDCBinNum, totMin, totMax);
    }

    // Fill histogram
    if (raw){
      for (const auto& x : vCDetPaddleRawTot[pixel]) {
        hTot[i]->Fill(x);
      }
    }
    if (!raw){
      for (const auto& x : vCDetPaddleCutTot[pixel]) {
        hTot[i]->Fill(x);
      }
    }

    // Draw
    cTot->cd(i + 1);
    hTot[i]->Draw();

    // If unused pixel: draw a black square in the top-right corner of the pad
    if (kUnusedCDetPixels.count(pixel)) {
      // NDC coordinates: (x1,y1,x2,y2) in [0,1] pad coordinates
      TPaveText* flag = new TPaveText(0.82, 0.82, 0.95, 0.95, "NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetBorderSize(1);
      flag->AddText("");       // empty; just a filled box
      flag->Draw("same");
    }
  }

  cTot->Update();
}

void getRate(int pixel, bool cut = false){
  if (!cut) std:: cout << "Rate in Pixel " << pixel << " = " << chanRates[pixel] << std::endl;
  if (cut) std:: cout << "Rate in Pixel " << pixel << " (with Tot Cut) = " << cutChanRates[pixel] << std::endl;
}

void plotRateVsID(bool raw = true){
  TH1::AddDirectory(kFALSE);
    // --- constants ---
  const int NCHAN_TOTAL = 2688;
  const int NCHAN_LAYER = 1344;
  const int NCHAN_SIDE  = 672;
  const int NMOD        = 3;
  const int NCHAN_MOD   = NCHAN_SIDE / NMOD; // 224

  // Helper: create one segment histogram and fill from chanRates
  auto MakeRateHist = [&](const char* hname,
                          const char* htitle,
                          int idStart, int idEnd,
                          double yMax = 1.5) -> TH1D* {
    const int nbins = idEnd - idStart + 1; // inclusive
    TH1D* h = new TH1D(hname, htitle, nbins, idStart, idEnd + 1); // [start, end+1)
    h->SetStats(0);
    h->SetMinimum(0.0);
    h->SetMaximum(yMax);

    for (int id = idStart; id <= idEnd; id++) {
      const int bin = h->FindBin(id);
      if (raw){
        h->SetBinContent(bin, chanRates[id]);
      }
      if (!raw){
        h->SetBinContent(bin, cutChanRates[id]);
      }
    }
    return h;
  };

  struct Seg { int layer; int mod; const char* side; int start; int end; };

  std::vector<Seg> segs;
  segs.reserve(12);

  auto AddLayerSegs_LeftThenRight = [&](int layer, int base) {
    const int L0 = base + 0;
    const int R0 = base + NCHAN_SIDE;

    // Left side: M1, M2, M3
    for (int m = 0; m < NMOD; m++) {
      int s = L0 + m*NCHAN_MOD;
      int e = s + NCHAN_MOD - 1;
      segs.push_back({layer, m+1, "L", s, e});
    }
    // Right side: M1, M2, M3
    for (int m = 0; m < NMOD; m++) {
      int s = R0 + m*NCHAN_MOD;
      int e = s + NCHAN_MOD - 1;
      segs.push_back({layer, m+1, "R", s, e});
    }
  };

  // Layer 1: IDs 0..1343
  AddLayerSegs_LeftThenRight(1, 0);
  // Layer 2: IDs 1344..2687
  AddLayerSegs_LeftThenRight(2, 1344);

  // --- build histograms ---
  TH1D* hRateSeg[12] = {nullptr};

  for (int i = 0; i < 12; i++) {
    const auto& s = segs[i];
    if (raw) {
      TString name  = Form("hRateVsIDL%dM%d%s", s.layer, s.mod, s.side);
      TString title = Form("CDet L%d %s M%d Rate vs Paddle ID;Paddle ID;Rate",
                          s.layer, s.side, s.mod);
      hRateSeg[i] = MakeRateHist(name.Data(), title.Data(), s.start, s.end, 1.5);
    }
    if (!raw){
      TString name  = Form("hRateVsIDL%dM%d%s", s.layer, s.mod, s.side);
      TString title = Form("CDet L%d %s M%d Rate w/Cut vs Paddle ID;Paddle ID;Rate",
                          s.layer, s.side, s.mod);
      hRateSeg[i] = MakeRateHist(name.Data(), title.Data(), s.start, s.end, 1.5);
    }
  }

  // --- Draw: Layer 1 canvas (Left M1-3 then Right M1-3) ---
  TCanvas* cRateL1 = new TCanvas("cRateL1", "CDet Rate vs ID (Layer 1)", 1400, 800);
  cRateL1->Divide(3,2); // top row: left M1-3, bottom row: right M1-3

  int pad = 1;
  for (int i = 0; i < 12; i++) {
    if (segs[i].layer != 1) continue;
    cRateL1->cd(pad++);
    hRateSeg[i]->Draw("HIST");
  }

  // --- Draw: Layer 2 canvas (Left M1-3 then Right M1-3) ---
  TCanvas* cRateL2 = new TCanvas("cRateL2", "CDet Rate vs ID (Layer 2)", 1400, 800);
  cRateL2->Divide(3,2);

  pad = 1;
  for (int i = 0; i < 12; i++) {
    if (segs[i].layer != 2) continue;
    cRateL2->cd(pad++);
    hRateSeg[i]->Draw("HIST");
  }
}

TCanvas *plotBarRateHV() {
  
  TCanvas *daa = new TCanvas("All TDC", "All TDC", 50,50,800,800);

  daa->cd();
  gPad->SetLogx();
  hBarRateHV->Draw();

  return daa;
}	

TCanvas *plotAllTDC(double width = 1, double binLow=0, double binHigh=60){
  double Nbins = ((binHigh-binLow)/width);
  //define histograms
  hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            Nbins, binLow, binHigh);
  hAllRawTe = new TH1F(TString::Format("hRawTe"),
            TString::Format("hRawTe"),
            Nbins, binLow, binHigh+TotBinHigh);
  hAllRawTot = new TH1F(TString::Format("hRawTot"),
            TString::Format("hRawTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllRawPMT = new TH1F(TString::Format("hRawPMT"),
            TString::Format("hRawPMT"),
            nTdc, 0, nTdc);
  hAllRawBar = new TH1F(TString::Format("hRawBar"),
            TString::Format("hRawBar"),
            168, 0, 168);
  hAllGoodLe = new TH1F(TString::Format("hAllGoodLe"),
            TString::Format("hAllGoodLe"),
            Nbins, binLow, binHigh);
  hAllGoodTe = new TH1F(TString::Format("hAllGoodTe"),
            TString::Format("hAllGoodTe"),
            Nbins, binLow, binHigh+TotBinHigh);
  hAllGoodTot = new TH1F(TString::Format("hAllGoodTot"),
            TString::Format("hAllGoodTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllGoodPMT = new TH1F(TString::Format("hAllGoodPMT"),
            TString::Format("hAllGoodPMT"),
            nTdc, 0, nTdc);
  hAllGoodBar = new TH1F(TString::Format("hAllGoodBar"),
            TString::Format("hAllGoodBar"),
            168, 0, 168);

  // Per-bar good leading-edge histograms (bar = GoodElID/16)
  hBarGoodLe.assign(NumPMTs, nullptr);
  for (int bar = 0; bar < NumPMTs; ++bar) {
    hBarGoodLe[bar] = new TH1F(TString::Format("hBarGoodLe_Bar%d", bar),
                               TString::Format("Good LE (Bar %d)", bar),
                               Nbins, binLow, binHigh);
    for (double x : vBarGoodLe[bar]) hBarGoodLe[bar]->Fill(x);
  }


  //fill necessary histograms from vectors
  for (double x : vAllRawLe) hAllRawLe->Fill(x);
  for (double x : vAllRawTe) hAllRawTe->Fill(x);
  for (double x : vAllRawTot) hAllRawTot->Fill(x);
  for (double x : vAllRawPMT) hAllRawPMT->Fill(x);
  for (double x : vAllRawBar) hAllRawBar->Fill(x);

  for (double x : vAllGoodLe) hAllGoodLe->Fill(x);
  for (double x : vAllGoodTe) hAllGoodTe->Fill(x);
  for (double x : vAllGoodTot) hAllGoodTot->Fill(x);
  for (double x : vAllGoodPMT) hAllGoodPMT->Fill(x);
  for (double x : vAllGoodBar) hAllGoodBar->Fill(x);

  // ------------------------------------------------------------
  // Compute per-bar mean offsets relative to the global good-LE mean
  // Correction convention used here:  toffset[bar] = mean(all) - mean(bar)
  // (this is the value that should be ADDED to LE/TE for that bar)
  // Bars with < 20 entries get offset = 0.
  // ------------------------------------------------------------
  const double meanAllGoodLe = hAllGoodLe->GetMean();
  std::vector<double> vBarGoodLeOffset(NumPMTs, 0.0);
  for (int bar = 0; bar < NumPMTs; ++bar) {
    if (!hBarGoodLe[bar]) continue;
    const double nent = hBarGoodLe[bar]->GetEntries();
    if (nent >= 20.0) {
      vBarGoodLeOffset[bar] = meanAllGoodLe - hBarGoodLe[bar]->GetMean();
    } else {
      vBarGoodLeOffset[bar] = 0.0;
    }
  }


// If offsets were not loaded from file, update global vector and write them out
// (so a first run can generate CDet_bar_toffsets.dat for subsequent corrected replays).
if (!gBarToffsetLoaded) {
  gBarToffsetCorr = vBarGoodLeOffset;

  std::ofstream fout(gBarToffsetFile.c_str());
  if (fout) {
    fout << "# CDet bar time offsets (ns) to ADD to LE/TE times: toffset = mean_all - mean_bar\n";
    fout << "# Generated by PlotElastic.C (plotAllTDC)\n";
    fout << "# Columns: bar(1..168)  toffset_ns  entries\n";
    fout.setf(std::ios::fixed); fout.precision(6);
    for (int bar = 0; bar < NumPMTs; ++bar) {
      const double ent = (hBarGoodLe[bar] ? hBarGoodLe[bar]->GetEntries() : 0.0);
      fout << (bar+1) << " " << vBarGoodLeOffset[bar] << " " << ent << "\n";
    }
    std::cout << "[CDet] Wrote bar offsets to '" << gBarToffsetFile << "'\n";
  } else {
    std::cout << "[CDet] ERROR: could not write offset file '" << gBarToffsetFile << "'\n";
  }
}

  // Plot offsets vs bar number (bars are numbered 1..168 on the x-axis)
  TCanvas* cOffsets = new TCanvas("cBarGoodLeOffsets", "Mean LE offsets vs Bar", 50, 900, 1200, 500);
  std::vector<double> xBar(NumPMTs), yOff(NumPMTs);
  for (int bar = 0; bar < NumPMTs; ++bar) {
    xBar[bar] = bar + 1;
    yOff[bar] = vBarGoodLeOffset[bar];
  }
  TGraph* grBarOffsets = new TGraph(NumPMTs, xBar.data(), yOff.data());
  grBarOffsets->SetName("grBarGoodLeOffsets");
  grBarOffsets->SetTitle(TString::Format(
      "Mean LE offsets vs Bar;Bar number;#mu_{all} - #mu_{bar} (ns)  (global mean = %.3f ns)",
      meanAllGoodLe));
  grBarOffsets->SetMarkerStyle(20);
  grBarOffsets->SetMarkerSize(0.8);
  grBarOffsets->Draw("AP");

  TCanvas *caa = new TCanvas("All TDC", "All TDC", 50,50,800,800);
  caa->Divide(2,3,0.01,0.01,0);
  TCanvas *caaa = new TCanvas("All Chan", "All Chan", 850,50,800,800);
  caaa->Divide(2,2,0.01,0.01,0);

  caa->cd(1);

  hAllRawLe->SetFillColor(kRed);
  hAllRawLe->SetMinimum(0.0);
  hAllRawLe->Draw();
  
  caa->cd(2);
  hAllGoodLe->SetFillColor(kBlue);
  hAllGoodLe->SetMinimum(0.0);
  hAllGoodLe->Draw();
  
  caa->cd(3);
  hAllRawTe->SetFillColor(kRed);
  hAllRawTe->SetMinimum(0.0);
  hAllRawTe->Draw();
  
  caa->cd(4);
  hAllGoodTe->SetFillColor(kBlue);
  hAllGoodTe->SetMinimum(0.0);
  hAllGoodTe->Draw();
  

  caa->cd(5);
  //gPad->SetLogy();
  hAllRawTot->SetFillColor(kRed);
  hAllRawTot->Draw();

  caa->cd(6);
  //gPad->SetLogy();
  hAllGoodTot->SetFillColor(kBlue);
  hAllGoodTot->Draw();

  caaa->cd(1);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllRawPMT->SetFillColor(kRed);
  hAllRawPMT->Draw();
  
  caaa->cd(2);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllGoodPMT->SetFillColor(kBlue);
  hAllGoodPMT->Draw();
  
  caaa->cd(3);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllRawBar->SetFillColor(kRed);
  hAllRawBar->Draw();
  
  caaa->cd(4);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllGoodBar->SetFillColor(kBlue);
  hAllGoodBar->Draw();



// ------------------------------------------------------------
// Draw per-bar Good LE histograms on 4 canvases (42 per canvas)
// Geometry:
//   Layer 1: bars 1-84   (index 0-83)
//   Layer 2: bars 85-168 (index 84-167)
//   Left side : bars 1-42,   85-126  (index 0-41,   84-125)
//   Right side: bars 43-84, 127-168  (index 42-83, 126-167)
// Each canvas: 6 rows x 7 cols = 42 pads
// ------------------------------------------------------------
auto drawBarRange = [&](const char* cname, const char* ctitle, int barStart, int nBars){
  TCanvas* c = new TCanvas(cname, ctitle, 1400, 900);
  c->Divide(7, 6, 0.001, 0.001);
  for(int i = 0; i < nBars; ++i){
    int bar = barStart + i;
    c->cd(i+1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);

    if(bar < 0 || bar >= (int)hBarGoodLe.size() || !hBarGoodLe[bar]) continue;

    TH1F* h = hBarGoodLe[bar];
    h->SetStats(0);

    // Make text readable in small pads
    h->GetXaxis()->SetTitleSize(0.07);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleSize(0.07);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.85);
    h->GetYaxis()->SetTitleOffset(0.85);

    // Optional: uncomment if you want log-y for visibility in low-stat pads
    // gPad->SetLogy();

    h->Draw("HIST");
  }
};

drawBarRange("cBarGoodLe_L1L", "Good LE by Bar: Layer 1 Left (Bars 1-42)",    0,   42);
drawBarRange("cBarGoodLe_L1R", "Good LE by Bar: Layer 1 Right (Bars 43-84)",  42,  42);
drawBarRange("cBarGoodLe_L2L", "Good LE by Bar: Layer 2 Left (Bars 85-126)",  84,  42);
drawBarRange("cBarGoodLe_L2R", "Good LE by Bar: Layer 2 Right (Bars 127-168)",126, 42);

  return caa;
}
void plot2DrefVsLE(double width = 1, double tmin=0, double tmax=60){
  int TDCBinNum = (int)((tmax-tmin)/width);
  TH2D *h2 = new TH2D("h2","t_ref vs t_scint; t_ref; t_scint",TDCBinNum, tmin, tmax, TDCBinNum, tmin, tmax);

  // Use reference times that are aligned with the GOOD-event vectors
  const size_t Nev = std::min(vGoodRefRawLe.size(), vGoodLe.size());
  for (size_t ev = 0; ev < Nev; ev++) { // iterate through good-selected events
    const double t_ref = vGoodRefRawLe[ev];
    if (std::isnan(t_ref)) continue; // if no ref recorded for this entry, skip

    const size_t Nhits = vGoodLe[ev].size();
    for (size_t ihit = 0; ihit < Nhits; ihit++) {
      const double t_det = vGoodLe[ev][ihit];
      h2->Fill(t_ref,t_det);
    }
  }
  TCanvas *crefVSle = new TCanvas("crefVSle", "CDet ref vs LE",900,700);
  h2->Draw("COLZ");
}

TH1* SubtractFitFromHist(const TH1* hIn, TF1* fFit, const char* outName = nullptr, bool clampNegToZero = true, int firstBin = 1, int lastBin = -1) {
  if (!hIn || !fFit) {
    std::cerr << "SubtractFitFromHist ERROR: null input.\n";
    return nullptr;
  }

  TH1* hSub = (TH1*)hIn->Clone(outName ? outName : (TString(hIn->GetName()) + "_sub").Data());
  hSub->SetDirectory(nullptr);

  if (lastBin < 0) lastBin = hSub->GetNbinsX();
  firstBin = std::max(firstBin, 1);
  lastBin  = std::min(lastBin, hSub->GetNbinsX());

  for (int ibin = firstBin; ibin <= lastBin; ++ibin) {
    const double x  = hSub->GetBinCenter(ibin);
    const double bw = hSub->GetBinWidth(ibin);

    const double content = hSub->GetBinContent(ibin);
    const double err     = hSub->GetBinError(ibin);

    // Average fit value over the bin (Integral/binwidth) so subtraction matches histogram binning
    const double fAvg = fFit->Integral(x - bw/2.0, x + bw/2.0) / bw;

    double newContent = content - fAvg;
    if (clampNegToZero && newContent < 0) newContent = 0;

    hSub->SetBinContent(ibin, newContent);
    hSub->SetBinError(ibin, err); // keep original errors (common choice)
  }

  // If you want to subtract everywhere else too, call with firstBin=1,lastBin=nbins.
  return hSub;
}

void plotCDetLayersTimeComp(double Width = 1, double diffMinCut = -15, double diffMaxCut = 15, double xdiffMinCut = -0.01, double xdiffMaxCut = 0.01, double LeMin = 0.02, double LeMax = 60, double TotMinCut = 0, double TotMaxCut = 70, double DiffMin = -20, double DiffMax = 20, double CDetMin = 0, double CDetMax = 60, double CDetTotMin = 0, double CDetTotMax = 80, double ECalMin = 62, double ECalMax = 130, double tdiffECalCDetMin = -100, double tdiffECalCDetMax = 100,bool allowMultiplePairs = true){
  
  TH1::AddDirectory(kFALSE);
  
  int TDCBinNum = (int)((DiffMax-DiffMin)/Width);
  int NADCBins = (int)((ECalMax-ECalMin)/4); //4ns bins for ECal, since fADC 4ns resolution

  TH1D* hCDetTimeDiff = new TH1D("hCDetTimeDiff", "CDet Layer Time Difference; Time Difference (ns);Counts", TDCBinNum, DiffMin, DiffMax);
  TH1D* hCDetXDiff = new TH1D("hCDetXDiff", "CDet Layer X Difference; X Diff (m);Counts", 600,-1.5,1.5);
  TH1D* hCDetLe1 = new TH1D("hCDetLe1", "CDet Layer 1 Good Time;Layer 1 LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetLe2 = new TH1D("hCDetLe2", "CDet Layer 2 Good Time;Layer 2 LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDetLe2vs1 = new TH2D("hCDetLe2vs1", "CDet Layer 2 Time vs CDet Layer 1 Time;Layer 1 LE (ns);Layer 2 LE (ns)",TDCBinNum,CDetMin,CDetMax,TDCBinNum,CDetMin,CDetMax);
  TH2D* hCDetTot2vs1 = new TH2D("hCDetTot2vs1", "CDet Layer 2 Tot vs CDet Layer 1 Tot;Layer 1 LE (ns);Layer 2 LE (ns)",TDCBinNum,CDetTotMin,CDetTotMax,TDCBinNum,CDetTotMin,CDetTotMax);
  TH2D* h2CDetx2VsCDetx1 = new TH2D("h2CDetx2VsCDetx1", "CDet Layer 2 x vs CDet Layer 1 x;CDet Layer 1 x (m);CDet Layer 2 x (m)",600,-1.5,1.5,600,-1.5,1.5);
  TH1I* hNpairPerEvent = new TH1I("hNpairPerEvent", "CDet accepted pairs per event;N_{pairs};Events", 20, 0, 20);
  TH1D* hCDetBarLe1 = new TH1D("hCDetBarLe1", "CDet Bar 1L27 Good Time;Layer 1 LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDetBarLe2 = new TH1D("hCDetBarLe2", "CDet Bar 2L27 Good Time;Layer 2 LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet1BarLeVsTot = new TH2D("hCDet1BarLeVsTot", "CDet Bar 1L27 Le vs Tot;Tot (ns);LE (ns)", TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet2BarLeVsTot = new TH2D("hCDet2BarLeVsTot", "CDet Bar 2L27 Le vs Tot;Tot (ns);LE (ns)", TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet1LeVsTot = new TH2D("hCDet1LeVsTot", "CDet Layer 1 Le vs Tot;Tot (ns);LE (ns)", TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hCDet2LeVsTot = new TH2D("hCDet2LeVsTot", "CDet Layer 2 Le vs Tot;Tot (ns);LE (ns)", TDCBinNum,CDetTotMin,CDetTotMax, TDCBinNum, CDetMin, CDetMax);
  TH2D* hECalVsCDetDt = new TH2D("hECalVsCDetDt", "ECal Time vs CDet dt;ECal ADC Time (ns);CDet dt_12 (ns)", NADCBins, ECalMin, ECalMax,TDCBinNum,DiffMin,DiffMax);
  TH2D* hECalVsCDetDtSingle = new TH2D("hECalVsCDetDtSingle", "CDet Single dt vs ECal Time;ECal ADC Time (ns);CDet dt_12 (ns)", NADCBins, ECalMin, ECalMax, TDCBinNum,DiffMin,DiffMax);
  TH2D* hECalVsCDetT = new TH2D("hECalVsCDetT", "CDet t vs ECal Time;ECal ADC Time (ns);CDet t (ns)", NADCBins, ECalMin, ECalMax,TDCBinNum,CDetMin,CDetMax);
  TH2D* hECalVsCDetTSingle = new TH2D("hECalVsCDetTSingle", "CDet Single t vs ECal Time;ECal ADC Time (ns);CDet t (ns)", NADCBins, ECalMin, ECalMax, TDCBinNum,CDetMin,CDetMax);

  TH2D* hCDetTimeDiffvsx1 = new TH2D("hCDetTimeDiffvsx1", "CDet Time Diff vs x1 position; x1 pos (m);CDet dt (ns)", 600,-1.5,1.5,TDCBinNum,DiffMin,DiffMax);
  TH2D* hCDetTimeDiffvsx2 = new TH2D("hCDetTimeDiffvsx2", "CDet Time Diff vs x2 position; x2 pos (m);CDet dt (ns)", 600,-1.5,1.5,TDCBinNum,DiffMin,DiffMax);

  TH2D* hCDetTimeDiffvsy1 = new TH2D("hCDetTimeDiffvsy1", "CDet Time Diff vs y1 position; y1 pos (m);CDet dt (ns)", 600,-0.5,0.5,TDCBinNum,DiffMin,DiffMax);
  TH2D* hCDetTimeDiffvsy2 = new TH2D("hCDetTimeDiffvsy2", "CDet Time Diff vs y2 position; y2 pos (m);CDet dt (ns)", 600,-0.5,0.5,TDCBinNum,DiffMin,DiffMax);

  TH1D* hDtCDetECal = new TH1D("hDtCDetECal", "CDet t - ECal t;CDet t - ECal t (ns);Counts", TDCBinNum, tdiffECalCDetMin, tdiffECalCDetMax);
  TH2D* hDtvsDxCDetECal = new TH2D("hDtvsDxCDetECal", "CDet ECal dt vs dx;dx_ECalCDet (m);dt_ECalCDet (ns)", NXDiffBins, XDiffLow, XDiffHigh, TDCBinNum, tdiffECalCDetMin, tdiffECalCDetMax);

  TH1D* hCDet1800Le = new TH1D("hCDet1800Le", "CDet Pixel 1800 Good Time; LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);
  TH1D* hCDet1107Le = new TH1D("hCDet1107Le", "CDet Pixel 1107 Good Time; LE (ns);Counts", TDCBinNum, CDetMin, CDetMax);

  const size_t Nev = vGoodLe.size();
  pairs_CDet.clear();
  pairs_CDet.resize(Nev);

  // ------------------------------
  // Tuning parameters 
  // ------------------------------
  double dt0   = 0.0;   // ns; start 0, later set to peak of dt histogram
  double dtWin = diffMaxCut;  // ns; start wide (10–15)
  double dxWin = xdiffMaxCut;   // m; ~one bar

  double sigT  = 1.0;   // ns; timing scale for score
  double sigX  = 0.01;   // m; start ~1, later tighten toward 0.5

  // dx histogram for ALL candidate pairs (before unique-matching).
  // Use a symmetric window around 0 that matches your initial dxWin.
  const double dxMin = -dxWin;
  const double dxMax =  dxWin;
  // Choose a reasonable binning for x; 0.5 cm/bin is a good start.
  const double dxBinW = 0.5; // cm
  const int dxBins = std::max(1, (int)std::round((dxMax - dxMin)/dxBinW));
  
  int sumNhits = 0;
  int sumGoodHits1 = 0;
  int sumGoodHits2 = 0;
  //vectors for plotting good hits in layers 1 and 2
  // std::vector<int> vGoodHits1PerEvent(Nev, 0);
  // std::vector<int> vGoodHits2PerEvent(Nev, 0);
  // std::vector<int> vGoodHitsPerEvent(Nev, 0);
  
  for (size_t ev = 0; ev < Nev; ev++) { //iterate through events
    sumNhits += vnhits1[ev]+vnhits2[ev];
    int countGoodHits1 = 0;
    int countGoodHits2 = 0;
    std::vector<double> vCDet1Time;
    std::vector<double> vCDet2Time;
    std::vector<double> vCDet1Tot;
    std::vector<double> vCDet2Tot;
    std::vector<double> vCDet1x;
    std::vector<double> vCDet2x;
    std::vector<double> vCDet1y;
    std::vector<double> vCDet2y;
    std::vector<double> vCDet1z;
    std::vector<double> vCDet2z;
    std::vector<double> vCDet1ID;
    std::vector<double> vCDet2ID;
    double t_ECal = v_GoodECalAdcTime[ev];
    // double x_ECal_actal = v_GoodECalX[ev];

    // ------------------------------
    // First pass: split hits into layers
    // ------------------------------
    const size_t Nhits = std::min(vGoodLe[ev].size(), vGoodTot[ev].size());
    for (size_t ihit = 0; ihit < Nhits; ++ihit) {
      if (vGoodLe[ev][ihit] >= LeMin && vGoodLe[ev][ihit] <= LeMax && vGoodTot[ev][ihit] >= TotMinCut && vGoodTot[ev][ihit] <= TotMaxCut && t_ECal >= ECalMin && t_ECal <= ECalMax){
        if (vGoodID[ev][ihit] >= 0 && vGoodID[ev][ihit] <= 1343){ //layer 1 hits
          vCDet1Time.push_back(vGoodLe[ev][ihit]);
          vCDet1Tot.push_back(vGoodTot[ev][ihit]);
          vCDet1ID.push_back(vGoodID[ev][ihit]);
          vCDet1x.push_back(vCDetGoodX[ev][ihit]);
          vCDet1y.push_back(vCDetGoodY[ev][ihit]);
          vCDet1z.push_back(vCDetGoodZ[ev][ihit]);
          countGoodHits1++;
        }
        else if (vGoodID[ev][ihit] >= 1344 && vGoodID[ev][ihit] <= 2687){//layer 2 hits
          vCDet2Time.push_back(vGoodLe[ev][ihit]);
          vCDet2Tot.push_back(vGoodTot[ev][ihit]);
          vCDet2ID.push_back(vGoodID[ev][ihit]);
          vCDet2x.push_back(vCDetGoodX[ev][ihit]);
          vCDet2y.push_back(vCDetGoodY[ev][ihit]);
          vCDet2z.push_back(vCDetGoodZ[ev][ihit]);
          countGoodHits2++;
        }
      }
    }

    //require at least one hit in each layer to even try pairing:
    if (vCDet1Time.empty() || vCDet2Time.empty()) continue;

    // ------------------------------
    // Second pass: build candidate pairs (all pairs that pass dt/dx cuts)
    // ------------------------------
    std::vector<Cand> cands;
    cands.reserve(vCDet1Time.size() * vCDet2Time.size());

    for (int i1 = 0; i1 < (int)vCDet1Time.size(); i1++) {
      for (int j2 = 0; j2 < (int)vCDet2Time.size(); j2++) {

        const double dt = vCDet2Time[j2] - vCDet1Time[i1];
        const double dx = vCDet2x[j2] - vCDet1x[i1];

        if (fabs(dt - dt0) > dtWin) continue;
        if (fabs(dx) > dxWin) continue;

        const double score = (((dt - dt0)*(dt - dt0)) / (sigT*sigT)) + ((dx*dx) / (sigX*sigX));

        cands.push_back({i1, j2, dt, dx, score});
      }
    }

    if (cands.empty()) continue;

    // ------------------------------
    // Sort candidates by score (best first)
    // ------------------------------
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) {
                return a.score < b.score;
              });

    // ------------------------------
    // Select non-conflicting pairs (greedy one-to-one matching)
    //   - ensures you don’t reuse the same hit in L1 or L2 twice
    // ------------------------------
    std::vector<char> used1(vCDet1Time.size(), 0);
    std::vector<char> used2(vCDet2Time.size(), 0);

    pairs_CDet[ev].clear();
    pairs_CDet[ev].reserve(std::min(vCDet1Time.size(), vCDet2Time.size())); 

    int nPairs = 0;

    for (const auto &c : cands) {
      if (used1[c.i1] || used2[c.j2]) continue;

      used1[c.i1] = 1;
      used2[c.j2] = 1;
      PairHit p;
      p.t1 = vCDet1Time[c.i1];
      p.tot1 = vCDet1Tot[c.i1];
      p.x1 = vCDet1x[c.i1];
      p.y1 = vCDet1y[c.i1];
      p.z1 = vCDet1z[c.i1];
      p.id1 = vCDet1ID[c.i1];

      p.t2 = vCDet2Time[c.j2];
      p.tot2 = vCDet2Tot[c.j2];
      p.x2 = vCDet2x[c.j2];
      p.y2 = vCDet2y[c.j2];
      p.z2 = vCDet2z[c.j2];
      p.id2 = vCDet2ID[c.j2];

      p.dt = c.dt;
      p.dx = c.dx;
      p.score = c.score;

      pairs_CDet[ev].push_back(p);

      if (!allowMultiplePairs) break;
    }
    hNpairPerEvent->Fill((int)pairs_CDet[ev].size());

    // ------------------------------
    // Third pass: fill histograms using the ACCEPTED pairs
    // ------------------------------    
    //if (ev <= 20) {std::cout << "CDet1 Size= " << vCDet1Time.size() << " CDet2 Size= " << vCDet2Time.size() << " Npairs= " << pairs_CDet[ev].size() << "\n";
   // }

    for (size_t ip = 0; ip < pairs_CDet[ev].size(); ip++) {
      const auto &p = pairs_CDet[ev][ip];

      if (ev <= 20) {
        if (ip == 0) std::cout << "--------new event ----------\n";
        std::cout << "pair=" << ip
                  << " t1=" << p.t1 << " x1=" << p.x1 << " TOT1= " << p.tot1 << " ID1=" << p.id1 << "\n";
        std::cout << "pair=" << ip
                  << " t2=" << p.t2 << " x2=" << p.x2 << " TOT2= " << p.tot2 << " ID2=" << p.id2 << "\n";
        std::cout << "pair=" << ip
                  << " tdiff=" << p.dt << " xdiff=" << p.dx
                  << " score=" << p.score << "\n";
        std::cout << "===============\n";
      }

      // Apply final/tighter cuts (these can be narrower than dtWin/dxWin)
      if (p.dt >= diffMinCut && p.dt <= diffMaxCut && p.dx >= xdiffMinCut && p.dx <= xdiffMaxCut) {
        double t_pair = (p.t1 + p.t2) / 2;
        double dt_EC = t_pair - t_ECal;
        if (dt_EC >= tdiffECalCDetMin && dt_EC <= tdiffECalCDetMax){
          double x_pair = (p.x1 + p.x2) / 2;
          double y_pair = (p.y1 + p.y2) / 2;
          double z_pair = (p.z1 + p.z2) / 2;
          double x_ECal = v_GoodECalX[ev]*z_pair/ECal_dist;
          double dx_EC = x_pair - x_ECal;

          hCDetTimeDiff->Fill(p.dt);
          hCDetXDiff->Fill(p.dx);
          hCDetLe2vs1->Fill(p.t1, p.t2);
          hCDetTot2vs1->Fill(p.tot1, p.tot2);
          hCDetLe1->Fill(p.t1);
          hCDetLe2->Fill(p.t2);
          h2CDetx2VsCDetx1->Fill(p.x1,p.x2);
          hCDet1LeVsTot->Fill(p.tot1,p.t1);
          hCDet2LeVsTot->Fill(p.tot2,p.t2);
          hECalVsCDetDt->Fill(t_ECal,p.dt);
          hCDetTimeDiffvsx1->Fill(p.x1, p.dt);
          hCDetTimeDiffvsx2->Fill(p.x2, p.dt);
          hCDetTimeDiffvsy1->Fill(p.y1, p.dt);
          hCDetTimeDiffvsy2->Fill(p.y2, p.dt);

          hECalVsCDetT->Fill(t_ECal,t_pair);

          hDtCDetECal->Fill(dt_EC);
          hDtvsDxCDetECal->Fill(dx_EC, dt_EC);

          //if (p.id1 >= 432 && p.id1 <= 447){
          if (p.id1 == 429){//417 hot hot hot
            hCDetBarLe1->Fill(p.t1);
            hCDet1BarLeVsTot->Fill(p.tot1, p.t1);
            hCDetBarLe2->Fill(p.t2);
            hCDet2BarLeVsTot->Fill(p.tot2, p.t2);
            hECalVsCDetDtSingle->Fill(t_ECal,p.dt);
            hECalVsCDetTSingle->Fill(t_ECal,p.t1);
          }
          if (p.id1 == 1107){
            hCDet1107Le->Fill(p.t1);
          }
          if (p.id2 == 1800){
            hCDet1800Le->Fill(p.t1);
          }
          //if (p.id2 >= 1776 && p.id2 <= 1791){
            //hCDetBarLe2->Fill(p.t2);
            //hCDet2BarLeVsTot->Fill(p.tot2, p.t2);
          //}
        }//ecal cdet tdiff cut
	    }//fill histogram with cuts
    }// end pair hits loop

  }// end event loop 
  // ------------------------------
  // Make Plots
  // ------------------------------    
  TCanvas *cCDetLayerTimes = new TCanvas("cCDetLayerTimes", "CDet Layer 1 and 2 LE",900,700);
  cCDetLayerTimes->Divide(1,2);

  cCDetLayerTimes->cd(1);
  //gPad->SetLogz();
  hCDetLe1->Draw();
  
  cCDetLayerTimes->cd(2);
  //gPad->SetLogz();
  hCDetLe2->Draw();
  
  // TCanvas *cCDetLeVsTot = new TCanvas("cCDetLeVsTot", "CDet LE vs Tot",900,700);
  // cCDetLeVsTot->Divide(1,2);

  // cCDetLeVsTot->cd(1);
  // //gPad->SetLogz();
  // hCDet1LeVsTot->Draw();
  
  // cCDetLeVsTot->cd(2);
  // //gPad->SetLogz();
  // hCDet2LeVsTot->Draw();


  // ----- plots for 1 bar ----- 
  TCanvas *cCDetLayerTimes1Bar = new TCanvas("cCDetLayerTimes1Bar", "CDet Single Paddles",900,700);
  cCDetLayerTimes1Bar->Divide(1,3);

  cCDetLayerTimes1Bar->cd(1);
  //gPad->SetLogz();
  hCDetBarLe1->Draw();
  
  cCDetLayerTimes1Bar->cd(2);
  //gPad->SetLogz();
  hCDet1107Le->Draw();

  cCDetLayerTimes1Bar->cd(3);
  hCDet1800Le->Draw();

  // TCanvas *cCDetLeVsTotBar = new TCanvas("cCDetLeVsTotBar", "CDet LE vs Tot (1 Paddle)",900,700);
  // cCDetLeVsTotBar->Divide(1,2);

  // cCDetLeVsTotBar->cd(1);
  // //gPad->SetLogz();
  // hCDet1BarLeVsTot->Draw();
  
  // cCDetLeVsTotBar->cd(2);
  // //gPad->SetLogz();
  // hCDet2BarLeVsTot->Draw();
  // ----- End plots for 1 bar -----

  TCanvas *cCDetTDiff = new TCanvas("cCDetTDiff", "CDet Time Diff",900,700);
  // ---- Gaussian fit on the "NoCuts" histogram ----
  // TF1 *fGaus = new TF1("fGaus", "gaus", DiffMin, DiffMax);
  // int maxBin = hCDetTimeDiff->GetMaximumBin();
  // double peakX = hCDetTimeDiff->GetBinCenter(maxBin);
  // fGaus->SetParameters(hCDetTimeDiff->GetMaximum(), peakX, diffMaxCut/2); // amp, mean, sigma guess

  // fGaus->SetRange(diffMinCut, diffMaxCut);
  // hCDetTimeDiff->Fit(fGaus, "R");

  // draw the fit on top of the already-drawn histogram
  // fGaus->SetLineColor(kBlack);

  hCDetTimeDiff->Draw();
  // fGaus->Draw("SAME");

  // TCanvas *cCDetXDiff = new TCanvas("cCDetXDiff", "CDet X Diff",900,700);
  // hCDetXDiff->Draw();

  TCanvas *cCDetLayer2v1 = new TCanvas("cCDetLayer2v1", "CDet Layer 2 vs 1",900,700);
  hCDetLe2vs1->SetMinimum(20);
  hCDetLe2vs1->Draw("COLZ");

  TCanvas *cCDetTotLayer2v1 = new TCanvas("cCDetTotLayer2v1", "CDet Tot Layer 2 vs 1",900,700);
  hCDetTot2vs1->SetMinimum(40);
  hCDetTot2vs1->Draw("COLZ");

  TCanvas *cNpair = new TCanvas("cNpair", "CDet accepted pairs per event", 900, 700);
  hNpairPerEvent->Draw();

  // TCanvas *cXplot = new TCanvas("cXplot", "CDet layer 2 vs 1 xposition", 900, 700);
  // h2CDetx2VsCDetx1->Draw("COLZ");

  // TCanvas * cDTvsECalT = new TCanvas("cDTvsECalT", " CDet dt vs ECal t", 900,700);
  // hECalVsCDetDt->Draw("COLZ");

  // TCanvas * cDTvsECalTSingle = new TCanvas("cDTvsECalTSingle", " CDet Single dt vs ECal t", 900,700);
  // hECalVsCDetDtSingle->Draw("COLZ");

  TCanvas * cCDetTvsECalT = new TCanvas("cCDetTvsECalT", " CDet t vs ECal t", 900,700);
  hECalVsCDetT->SetMinimum(40);
  hECalVsCDetT->Draw("COLZ");

  // ---------------------------------------------
  // Linear parametrization: <CDet t> vs ECal ADC t
  // Use a profile along X (ECal time) and fit with a straight line.
  // This yields a robust mean-trend fit for the 2D distribution.
  // ---------------------------------------------
  TProfile* pCDetTvsECalT = hECalVsCDetT->ProfileX("pCDetTvsECalT");
  pCDetTvsECalT->SetMarkerStyle(20);
  pCDetTvsECalT->SetMarkerSize(0.6);

  TF1* fCDetTvsECalT_lin = new TF1("fCDetTvsECalT_lin", "pol1", ECalMin, ECalMax);
  // Quiet fit, respect range
  pCDetTvsECalT->Fit(fCDetTvsECalT_lin, "QR");

  // Overlay the profile points and fit on top of the COLZ plot
  pCDetTvsECalT->Draw("SAME");
  fCDetTvsECalT_lin->Draw("SAME");

  const double p0 = fCDetTvsECalT_lin->GetParameter(0);
  const double p1 = fCDetTvsECalT_lin->GetParameter(1);
  const double e0 = fCDetTvsECalT_lin->GetParError(0);
  const double e1 = fCDetTvsECalT_lin->GetParError(1);
  std::cout << "\n[plotCDetLayersTimeComp] Linear fit for <CDet t> vs ECal ADC t:\n"
            << "  <t_CDet> = p0 + p1 * t_ECal\n"
            << "  p0 = " << p0 << " +/- " << e0 << " ns\n"
            << "  p1 = " << p1 << " +/- " << e1 << " (ns/ns)\n\n";

  TPaveText* ptFit = new TPaveText(0.14, 0.80, 0.52, 0.92, "NDC");
  ptFit->SetFillColor(0);
  ptFit->SetTextAlign(12);
  ptFit->AddText("<t_{CDet}> = p0 + p1 t_{ECal}");
  ptFit->AddText(Form("p0 = %.3f #pm %.3f ns", p0, e0));
  ptFit->AddText(Form("p1 = %.5f #pm %.5f", p1, e1));
  ptFit->Draw("SAME");

  TCanvas * cCDetTvsECalTSingle = new TCanvas("cCDetTvsECalTSingle", " CDet Single t vs ECal t", 900,700);
  hECalVsCDetTSingle->Draw("COLZ");

  TCanvas *cDtCDetECal = new TCanvas("cDtCDetECal", "CDet ECal dt",900,700);
  hDtCDetECal->Draw();

  TCanvas *cDtvsDxCDetECal = new TCanvas("cDtvsDxCDetECal", "CDet ECal dt vs dx",900,700);
  hDtvsDxCDetECal->Draw("COLZ");

  TCanvas *cDtCDetLayersvsPos = new TCanvas("cDtCDetLayersvsPos", "CDet dt vs position", 900,700);
  cDtCDetLayersvsPos->Divide(2,2);
  cDtCDetLayersvsPos->cd(1);

  hCDetTimeDiffvsx1->Draw();

  cDtCDetLayersvsPos->cd(2);
  hCDetTimeDiffvsx2->Draw();

  cDtCDetLayersvsPos->cd(3);
  hCDetTimeDiffvsy1->Draw();

  cDtCDetLayersvsPos->cd(4);
  hCDetTimeDiffvsy2->Draw();
}

void plotECalCDetTimeComp(double Width = 1, double diffMinCut = 70, double diffMaxCut = 115, double LeMin = 0.02, double LeMax = 60, double TotMin = 0, double TotMax = 150, double DiffMin = 0, double DiffMax = 130, double CDetTotMin = 0, double CDetTotMax = 80, double CDetMin = 0, double CDetMax = 60,double ECalMin = 62, double ECalMax = 130){
  
  int NADCBins = (int)((ECalMax-ECalMin)/4); //4ns bins for ECal, since fADC 4ns resolution
  int TDCBinNum = (int)((DiffMax-DiffMin)/Width);

  TH1D* hECalMinusCDetTime = new TH1D("hECalMinusCDetTime", "ECal-CDet Time;Time Diff (ns);Counts", TDCBinNum, DiffMin, DiffMax);
  TH1D* hECalMinusCDetTimeNoCuts = new TH1D("hECalMinusCDetTimeNoCuts", "ECal-CDet Time;Time Diff (ns);Counts", TDCBinNum, DiffMin, DiffMax);
  TH2D* h2ECalMinusCDetTime = new TH2D("h2ECalMinusCDetTime", "ECal-CDet Time vs CDet Time;CDet 'Good' Time (ns);ECal-CDet Time (ns)", TDCBinNum, CDetMin, CDetMax, TDCBinNum, DiffMin, DiffMax);
  TH2D* h2ECalMinusCDetTot = new TH2D("h2ECalMinusCDetTot", "ECal-CDet Time vs CDet Tot;CDet 'Good' Tot (ns);ECal-CDet Time (ns)", TDCBinNum, CDetTotMin, CDetTotMax, TDCBinNum, DiffMin, DiffMax);
  TH2D* hECalVsCDet = new TH2D("hECalVsCDet", "ECal Time vs CDet Time;CDet LE Time (ns);ECal ADC Time (ns)",TDCBinNum,CDetMin,CDetMax, NADCBins, ECalMin, ECalMax);
  TH2D* h2ECalxVsCDetx = new TH2D("h2ECalxVsCDetx", "ECal Good x vs CDet Good x;CDet Good x (m);ECal Good x (m)",600,-1.5,1.5,200,-1.5,1.5);
  TH2D* h2ECalxVsCDetxNoProject = new TH2D("h2ECalxVsCDetxNoProject", "ECal Actual x vs CDet Good x;CDet Good x (m);ECal Actual x (m)",600,-1.5,1.5,200,-1.5,1.5);
  // TH1I* hGoodHitsPerEvent = new TH1I("hGoodHitsPerEvent", "Good hits per event;N_{good hits};Events", 100, 0, 100);

  const size_t Nev = std::min(vGoodLe.size(),v_GoodECalAdcTime.size());
  int sumNhits = 0;
  int sumGoodHits1 = 0;
  int sumGoodHits2 = 0;
  //vectors for plotting good hits in layers 1 and 2
  // std::vector<int> vGoodHits1PerEvent(Nev, 0);
  // std::vector<int> vGoodHits2PerEvent(Nev, 0);
  // std::vector<int> vGoodHitsPerEvent(Nev, 0);

  for (size_t ev = 0; ev < Nev; ev++) { //iterate through events
    sumNhits += vnhits1[ev]+vnhits2[ev];
    int countGoodHits1 = 0;
    int countGoodHits2 = 0;

    double t_ECal = v_GoodECalAdcTime[ev];
    double x_ECal_actal = v_GoodECalX[ev];
  
    const size_t Nhits = std::min(vGoodLe[ev].size(), vGoodTot[ev].size());
    for (size_t ihit = 0; ihit < Nhits; ++ihit) {
      double t_CDet = vGoodLe[ev][ihit];
      double tot = vGoodTot[ev][ihit];
      double t_diff = t_ECal-t_CDet;
      hECalMinusCDetTimeNoCuts->Fill(t_diff);
      
      if (t_CDet >= LeMin && t_CDet <= LeMax && tot >= TotMin && tot <= TotMax && t_diff >= diffMinCut && t_diff <= diffMaxCut){
        double x_CDet = vCDetGoodX[ev][ihit];
        double x_ECal = v_GoodECalX[ev]*vCDetGoodZ[ev][ihit]/ECal_dist;
        hECalMinusCDetTime->Fill(t_diff);
        h2ECalMinusCDetTime->Fill(t_CDet, t_diff);
        h2ECalMinusCDetTot->Fill(tot, t_diff);
        hECalVsCDet->Fill(t_CDet, t_ECal);
        h2ECalxVsCDetx->Fill(x_CDet,x_ECal);
        h2ECalxVsCDetxNoProject->Fill(x_CDet,x_ECal_actal);

        int sbselemid = (Int_t)vGoodID[ev][ihit];
        int sbsrown = sbselemid%672;
        int sbscoln = sbselemid/672;
        int mylayern = sbscoln/2;
        int mypaddlen = sbscoln*672 + sbsrown;

        if (mylayern == 0) countGoodHits1++;
        else countGoodHits2++;

	    }//fill histogram with cuts
    }//finished looking at hits
    // vGoodHitsPerEvent[ev] = countGoodHits1 + countGoodHits2;
    // vGoodHits1PerEvent[ev] = countGoodHits1;
    // vGoodHits2PerEvent[ev] = countGoodHits2;
    
    sumGoodHits1 += countGoodHits1; 
    sumGoodHits2 += countGoodHits2;

  } //finished looking at all events
  double aveNhits = (sumNhits / Nev);
  double aveGoodHits = ((sumGoodHits1+sumGoodHits2) / Nev);
  double aveGoodHits1 = (sumGoodHits1 / Nev);
  double aveGoodHits2 = (sumGoodHits2 / Nev);
  std::cout << "ave nhits/event w/o cuts= " << aveNhits << endl;
  std::cout << "ave nGoodHits/event= " << aveGoodHits << endl;
  std::cout << "ave nGoodHits Layer 1/event= " << aveGoodHits1 << endl;
  std::cout << "ave nGoodHits Layer 2/event= " << aveGoodHits2 << endl;

  //make profile to fit time diff vs CDet time with linear fit
  // TProfile *prof = h2ECalMinusCDetTime->ProfileX("prof");
  // prof->Fit("pol1");

  // TF1 *fit = prof->GetFunction("pol1");  // retrieve automatic fit

  //make canvas and draw hist
  TCanvas *cTimeDiff = new TCanvas("cTimeDiff", "ECal ADCtime Minus CDet Good LE",900,700);

  hECalMinusCDetTime->SetLineColor(kRed);
  hECalMinusCDetTimeNoCuts->SetLineColor(kBlue);
  hECalMinusCDetTimeNoCuts->Draw("HIST");
  hECalMinusCDetTime->Draw("SAME");

  // ---- Gaussian fit on the "NoCuts" histogram ----
  TF1 *fGausNoCuts = new TF1("fGausNoCuts", "gaus", DiffMin, DiffMax);
  int maxBin = hECalMinusCDetTimeNoCuts->GetMaximumBin();
  double peakX = hECalMinusCDetTimeNoCuts->GetBinCenter(maxBin);
  fGausNoCuts->SetParameters(hECalMinusCDetTimeNoCuts->GetMaximum(), peakX, 20.0); // amp, mean, sigma guess

  // (optional) restrict fit range around the peak so you fit the main bump, not the tails/background
  double fitLo = peakX - 40.0;
  double fitHi = peakX + 40.0;
  fGausNoCuts->SetRange(fitLo, fitHi);

  // do the fit ("R" uses the TF1 range, "0" suppresses ROOT fit printout if you want)
  hECalMinusCDetTimeNoCuts->Fit(fGausNoCuts, "R");

  // draw the fit on top of the already-drawn histogram
  fGausNoCuts->SetLineColor(kBlack);
  fGausNoCuts->Draw("SAME");

  std::cout << "NoCuts Gaussian: mean = " << fGausNoCuts->GetParameter(1) << "  sigma = " << fGausNoCuts->GetParameter(2) << std::endl;

  auto leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(hECalMinusCDetTime,"Cuts","l");
  leg->AddEntry(hECalMinusCDetTimeNoCuts,"NoCuts","l");
  leg->Draw();

  //gaussian subtracted hist
  TH1* hNoCutsSub = SubtractFitFromHist(hECalMinusCDetTimeNoCuts, fGausNoCuts);
  hNoCutsSub->SetLineColor(kMagenta);
  hNoCutsSub->Draw("HIST SAME");

  TCanvas *cXComp = new TCanvas("cXComp", "ECal x vs CDet x",900,700);
  cXComp->SetLogz();
  h2ECalxVsCDetx->Draw("COLZ");

  TCanvas *c2DtimeComps = new TCanvas("c2DtimeComps", "ECal-CDet Time Comparisons",900,700);
  c2DtimeComps->Divide(1,3);

  c2DtimeComps->cd(1);
  //gPad->SetLogz();
  hECalVsCDet->Draw("COLZ");

  c2DtimeComps->cd(2);
  //gPad->SetLogz();
  h2ECalMinusCDetTime->Draw("COLZ"); //heatmap
  
  c2DtimeComps->cd(3);
  //gPad->SetLogz();
  h2ECalMinusCDetTot->Draw("COLZ"); //heatmap

  TCanvas *cXCompActual = new TCanvas("cXCompActual", "ECal x vs CDet x w/o projection",900,700); //uses non projected ecal x
  cXCompActual->SetLogz();
  h2ECalxVsCDetx->Draw("COLZ");

} //end routine

void plotRawXCorrelation(double tDiffMin = 80, double tDiffMax = 100){

  TH2D* h2RawECalxVsCDetx = new TH2D("h2ECalxVsCDetx", "ECal Good x vs CDet Good x;CDet Good x (m);ECal Good x (m)",600,-1.5,1.5,200,-1.5,1.5);

  if (vCDetX.size() != v_ECalX.size()){
    std::cout << "vCDetX and vECalX not the same size" << std::endl;
  }
  const size_t Nev = vCDetX.size();
  for (size_t ev = 0; ev < Nev; ev++){
    double t_e = v_ECalAdcTime[ev];
    //std::cout << "event = " << ev << " " << "size of cdetX = " << vCDetX[ev].size() << std::endl;
    //std::cout << " " <<std::endl;
    //std::cout << "event = " << ev << " " << "size of cdetLE = " << vRawLe[ev].size() << std::endl;
    //std::cout << "event = " << ev << " " << "size of cdetLE = " << vRawTe[ev].size() << std::endl;
    //std::cout << "event = " << ev << " " << "size of cdetLE = " << vRawTot[ev].size() << std::endl;
    const size_t Nhits = vCDetX[ev].size();
    for (size_t ihit = 0; ihit < Nhits; ihit++){
      double t_c = vRawLe[ev][ihit];
      double t_diff = t_e - t_c;
      if (t_diff >= tDiffMin && t_diff <= tDiffMax){
        double x_c = vCDetX[ev][ihit];
        double z_c = vCDetZ[ev][ihit];
        double x_e = v_ECalX[ev]*(vCDetZ[ev][ihit]/ECal_dist);
        //std::cout << "CDetX = " << x_c << " ECalX = " << x_e << "CDetZ = " << z_c << std::endl;
        h2RawECalxVsCDetx->Fill(x_c,x_e);
      }//if statement for cuts
    }//end hit loop
  }//end event loop
  TCanvas *cRawXComp = new TCanvas("cRawXComp", "ECal x vs Raw CDet x",900,700);
  h2RawECalxVsCDetx->Draw("COLZ");
}

void plotTimeECalVsCDet(double Width = 0.0160167/2, double LeMin = 0.02, double LeMax = 60, double TotMin = 0, double TotMax = 150, double CDetMin = 0, double CDetMax = 60, double ECalMin = 0, double ECalMax = 250){
  int NADCBins = (int)((ECalMax-ECalMin)/4); //4ns bins for ECal, since fADC 4ns resolution
  int TDCBinNum = (int)((CDetMax-CDetMin)/Width);
  TH2D* hECalVsCDet = new TH2D("hECalVsCDet", "ECal Time vs CDet Time;CDet LE Time (ns);ECal ADC Time (ns)",TDCBinNum,CDetMin,CDetMax, NADCBins, ECalMin, ECalMax);
  const size_t Nev = std::min(vGoodLe.size(),v_GoodECalAdcTime.size());

  for (size_t ev = 0; ev < Nev; ev++) { //iterate through events
    double t_ECal = v_GoodECalAdcTime[ev];

    const size_t Nhits = std::min(vGoodLe[ev].size(), vGoodTot[ev].size());
    for (size_t ihit = 0; ihit < Nhits; ++ihit) {
      double t_CDet = vGoodLe[ev][ihit];
      double tot = vGoodTot[ev][ihit];
      if (t_CDet >= LeMin && t_CDet <= LeMax && tot >= TotMin && tot <= TotMax){
        hECalVsCDet->Fill(t_CDet, t_ECal);
      }
    }//finished looking at hits
  } //finished looking at all events

  //make canvas and draw hist
  TCanvas *cTimeComp = new TCanvas("cTimeComp", "ECal ADCtime vs CDet LE",900,700);
  hECalVsCDet->Draw("COLZ");
} //end routine

void plotXDiffSections(double le_min = 0, double le_max = 60, double tDiffMin = 80, double tDiffMax = 100){
  //define nonchanging histograms
  TH1D* h1 = new TH1D("h1", "xDiffLayer1;xdiff (m);Counts", NXDiffBins,XDiffLow,XDiffHigh);
  TH1D* h2 = new TH1D("h2", "xDiffLayer2;xdiff (m);Counts", NXDiffBins,XDiffLow,XDiffHigh);
  TH1D* h3 = new TH1D("h3", "good le;le (ns);Counts", NTDCBins, le_min, le_max);

  const size_t Nev = vCDetGoodX.size();
  for (size_t ev = 0; ev < Nev; ++ev) {
    double t_ECal = v_GoodECalAdcTime[ev];
    for (int n = 0; n < vCDetGoodX[ev].size(); n++){
      double t_CDet = vGoodLe[ev][n];
      double t_diff = t_ECal-t_CDet;
      if (vGoodLe[ev][n] <= le_max && vGoodLe[ev][n] >= le_min && t_diff >= tDiffMin && t_diff <= tDiffMax){
        if (vGoodLayer[ev][n]==0){ //layer 1
          h1->Fill(vCDetGoodX[ev][n]-(v_GoodECalX[ev]*vCDetGoodZ[ev][n]/ECal_dist));
        }
        else if (vGoodLayer[ev][n]==1){ //layer 2
          h2->Fill(vCDetGoodX[ev][n]-(v_GoodECalX[ev]*vCDetGoodZ[ev][n]/ECal_dist));
        }
        h3->Fill(vGoodLe[ev][n]);
      }
    }
  } // individual time step histograms filled

  TCanvas *cTimeComp = new TCanvas("cTimeComp", "ECal ADCtime vs CDet LE",900,700);
  TCanvas* c1 = new TCanvas("c1", "xDiff Layer 1",900,700);
  h1->Draw("HIST");

  TCanvas* c2 = new TCanvas("c2", "xDiff Layer 2",900,700);
  h2->Draw("HIST");

  TCanvas* c3 = new TCanvas("c3", "le",900,700);
  h3->Draw("HIST");

} //end plotXDiffSections

void plotPMTRates(Int_t mymodule=1, Int_t mylayer=1, Int_t choice = 1){

  Int_t offsetl = (mylayer-1)*1344 + (mymodule-1)*224;
  Int_t offsetr = (mylayer-1)*1344 + 672 + (mymodule-1)*224;
  Int_t xcpos = 50 + (choice-1)*400;

  TCanvas *caPMT = new TCanvas(TString::Format("PMT Rates Left %d",choice), TString::Format("PMT Rates Left %d",choice), xcpos,50,400,550);
  caPMT->Divide(2,7,0.01,0.01,0);
  TCanvas *caPMTT = new TCanvas(TString::Format("PMT Rates Right %d",choice), TString::Format("PMT Rates Right %d",choice), xcpos,650,400,550);
  caPMTT->Divide(2,7,0.01,0.01,0);

  Double_t histymax = 12500;
  Double_t rate_levelu = 1800; // Expect 600 for 50000 events @ 5uA
  Double_t rate_levell = 600; // Expect 600 for 50000 events @ 5uA
  TLine *lineu = new TLine(0, rate_levelu, 2688, rate_levelu);
  TLine *linel = new TLine(0, rate_levell, 2688, rate_levell);

  for (Int_t ii=0; ii<14; ii++) {
    caPMT->cd(ii+1);
    gPad->SetLogy();
    gPad->DrawFrame(offsetl+ii*16,1.,offsetl+ii*16+16,histymax);
    hAllRawPMT->Draw("sames");

    //draw unused pixels on each PMT plot
    double xmin = offsetl + ii * 16;
    double xmax = xmin + 16;
    for (double x : missingPixelBins) {
      if (x < xmin || x >= xmax) continue;

      int bin = hAllRawPMT->FindBin(x);
      double xlow = hAllRawPMT->GetBinLowEdge(bin);
      double xup = xlow + hAllRawPMT->GetBinWidth(bin);
      double yup = hAllRawPMT->GetBinContent(bin);
      
      TBox *box = new TBox(xlow, 0, xup, yup);
      box->SetFillColor(kBlack);
      box->SetFillStyle(1001);
      box->Draw("sames");
    }
	
    lineu->Draw("sames");
    linel->Draw("sames");
  }
  caPMT->Update();

  for (Int_t ii=0; ii<14; ii++) {
    caPMTT->cd(ii+1);
    gPad->SetLogy();
    gPad->DrawFrame(offsetr+ii*16,1.,offsetr+ii*16+16.,histymax);
    hAllRawPMT->Draw("sames");

    double xmin = offsetr + ii * 16;
    double xmax = xmin + 16;
    for (double x : missingPixelBins) {
      if (x < xmin || x >= xmax) continue;

      int bin = hAllRawPMT->FindBin(x);
      double xlow = hAllRawPMT->GetBinLowEdge(bin);
      double xup = xlow + hAllRawPMT->GetBinWidth(bin);
      double yup = hAllRawPMT->GetBinContent(bin);
      
      TBox *boxx = new TBox(xlow, 0, xup, yup);
      boxx->SetFillColor(kBlack);
      boxx->SetFillStyle(1001);
      boxx->Draw("sames");
    }
    
    lineu->Draw("sames");
    linel->Draw("sames");
  }
  caPMTT->Update();
  
 
  return;


}
TCanvas *plotBarRates(){

  TCanvas *cabar = new TCanvas("BarRates", "BarRates", 50,50,800,800);
  cabar->Divide(3,4,0.01,0.01,0);

  Double_t histymax = 200000.0;
  Double_t rate_level = 600; // Expect 600 for 50000 events @ 5uA
  TLine *line = new TLine(0, rate_level, 168, rate_level);

  cabar->cd(1);
  gPad->SetLogy();
  gPad->DrawFrame(0.,1.,14.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(2);
  gPad->SetLogy();
  gPad->DrawFrame(14.,1.,28.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(3);
  gPad->SetLogy();
  gPad->DrawFrame(28.,1.,42.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(4);
  gPad->SetLogy();
  gPad->DrawFrame(42.,1.,56.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(5);
  gPad->SetLogy();
  gPad->DrawFrame(56.,1.,70.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(6);
  gPad->SetLogy();
  gPad->DrawFrame(70.,1.,84.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(7);
  gPad->SetLogy();
  gPad->DrawFrame(84.,1.,98.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(8);
  gPad->SetLogy();
  gPad->DrawFrame(98.,1.,112.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(9);
  gPad->SetLogy();
  gPad->DrawFrame(112.,1.,126.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(10);
  gPad->SetLogy();
  gPad->DrawFrame(126.,1.,140.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(11);
  gPad->SetLogy();
  gPad->DrawFrame(140.,1.,154.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(12);
  gPad->SetLogy();
  gPad->DrawFrame(154.,1.,168.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  return cabar;


}

TCanvas *plotGoodTDC2D(){

  TCanvas *cac = new TCanvas("all2d", "all2d", 50,50,800,800);
  cac->Divide(2,2,0.01,0.01,0);
  
  cac->cd(1);
  gPad->SetLogz();
  h2AllGoodLe->Draw("colz");
  cac->cd(2);
  gPad->SetLogz();
  h2AllGoodTe->Draw("colz");
  cac->cd(3);
  gPad->SetLogz();
  h2AllGoodTot->Draw("colz");

  return cac;
}

TCanvas *plotRefTDC() {
  hRefRawLe = new TH1F(TString::Format("hRefRawLe"),
            TString::Format("hRefRawLe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefRawTe = new TH1F(TString::Format("hRefRawTe"),
            TString::Format("hRefRawTe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefRawTot = new TH1F(TString::Format("hRefRawTot"),
            TString::Format("hRefRawTot"),
            RefNTotBins, RefTotBinLow, RefTotBinHigh);
  hRefRawPMT = new TH1F(TString::Format("hRefRawPMT"),
            TString::Format("hRefRawPMT"),
            32, 2688, 2720);

  
  //fill histograms
  for (double val : vRefRawLe)   hRefRawLe->Fill(val);
  for (double val : vRefRawTe)   hRefRawTe->Fill(val);
  for (double val : vRefRawTot)  hRefRawTot->Fill(val);
  for (int    val : vRefRawPMT)  hRefRawPMT->Fill(val);

  //make canvas
  TCanvas *cbb = new TCanvas("ref", "ref", 850,50, 1200,800);
  cbb->Divide(2,2,0.01,0.01,0);

  cbb->cd(1);
  gPad->SetLogy();
  hRefRawLe->Draw();

  cbb->cd(2);
  gPad->SetLogy();
  hRefRawTe->Draw();

  cbb->cd(3);
  gPad->SetLogy();
  hRefRawTot->Draw();

  cbb->cd(4);
  gPad->SetLogy();
  hRefRawPMT->Draw();

  return cbb;
}


TCanvas *plotCDetTDC(){

  TCanvas *canvas[NumHalfModules];
  for (int cmodule=1;cmodule<=NumModules;cmodule++) {
   for (int layer=1;layer<=NumLayers;layer++) {
    for (int side=1;side<=NumSides;side++){
     int xposition1 = 100 + (layer-1)*1000 + (side-1)*500;
     int yposition1 = 100 + (NumModules-cmodule)*400;
     int xposition2 = 150 + (layer-1)*1000 + (side-1)*500;
     int yposition2 = 150 + (NumModules-cmodule)*400;

     int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
     int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;

     TString cname;
     cname.Form("c1_%d_%d_%d",cmodule,layer,side);
     TCanvas *c1 = new TCanvas(cname, cname, xposition1,yposition1,380,280);
     c1->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c1;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c1->cd(ii+1);
        if (ii == 1) {
		cout << "side_group_plot = " << side_group_plot << endl;
	}
	
	hRawTe[elemID_start + ii ]->Draw();
	
     }

     cname.Form("c2_%d_%d_%d",cmodule,layer,side);
     TCanvas *c2 = new TCanvas(cname, cname, xposition2,yposition2,380,280);
     c2->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c2;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c2->cd(ii+1);
        if (ii == 1) {
                cout << "side_group_plot = " << side_group_plot << endl;
        }

        hRawLe[elemID_start + ii + 1]->Draw();

     }


    }
   }
  }

  return canvas[0];

}

TCanvas *plotHalfModule(int cmodule = 1, int side = 1, int layer = 1){


  TCanvas *c4 = new TCanvas("c4", "c4", 50,50,1000,1200);
  c4->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
  TCanvas *c4b = new TCanvas("c4b", "c4b", 1050,50,1000,1200);
  c4b->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);

  int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;


  for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c4->cd(ii+1);
  	hRawTe[elemID_start + ii ]->Draw();
	c4b->cd(ii+1);
  	hRawLe[elemID_start + ii ]->Draw();
  }

  return c4;

}

TCanvas *plotBarTDC(int bar = 39, int side = 1, int layer = 1){

  int mymodule = (bar-1)/NumBars+1; //mymodule in this case represnts top, middle, or bottom as opposed todetector labels whe
  int paddle_start = (bar - (mymodule-1)*NumBars - 1)*NumPaddles;
  std::cout << "mymodule = " << mymodule << "  paddle_start = " << paddle_start << std::endl;

  int side_group = 6*(layer-1) + 2*(mymodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (mymodule-1)*NumBars*NumPaddles + paddle_start;


  std::cout << "side_group = " << side_group << std::endl;

  TCanvas *c3 = new TCanvas("c3", "c3", 150,150,600,450);
  c3->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333 = new TCanvas("c333", "c333", 750,150,600,450);
  c333->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c3a = new TCanvas("c3a", "c3a", 150,650,600,450);
  c3a->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333a = new TCanvas("c333a", "c333a", 750,650,600,450);
  c333a->Divide(4,4, 0.01, 0.01, 0);

  for (int ii = 0; ii < NumPaddles; ii++) {

        c3->cd(ii+1);
  	hRawLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c333->cd(ii+1);
  	hGoodLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c3a->cd(ii+1);
  	hRawTe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c333a->cd(ii+1);
  	hGoodTe[elemID_start + ii ]->Draw();
  }

  return c3;

}

TCanvas *plotTOTvsLE(){

  TCanvas *c123 = new TCanvas("c123", "c123", 50,50,1000,1000);

  c123->cd();
  h2TDCTOTvsLE->Draw("colz");

  return c123;

}

TCanvas *plotTOTvsXDiff(){

  TCanvas *c1234 = new TCanvas("c1234", "c1234", 50,50,1000,1000);
  c1234->Divide(2,2,0.01,0.01,0);

  c1234->cd(1);
  h2TOTvsXDiff1->Draw("colz");
  c1234->cd(2);
  h2TOTvsXDiff2->Draw("colz");
  c1234->cd(3);
  h2LEvsXDiff1->Draw("colz");
  c1234->cd(4);
  h2LEvsXDiff2->Draw("colz");

  return c1234;

}


TCanvas *plotRowColLayer(){


  TCanvas *c5 = new TCanvas("c5", "c5", 50,50,800,800);
  c5->Divide(4,2, 0.01, 0.01, 0);

  c5->cd(1);
  gPad->SetLogy();
  hRowLayer1Side1->Draw();
  c5->cd(2);
  gPad->SetLogy();
  hRowLayer1Side2->Draw();
  c5->cd(3);
  gPad->SetLogy();
  hRowLayer2Side1->Draw();
  c5->cd(4);
  gPad->SetLogy();
  hRowLayer2Side2->Draw();
  c5->cd(5);
  gPad->SetLogy();
  hRow->Draw();
  c5->cd(6);
  gPad->SetLogy();
  hCol->Draw();
  c5->cd(7);
  gPad->SetLogy();
  hLayer->Draw();
  c5->cd(8);
  gPad->SetLogy();
  hHitPMT->Draw();

  return c5;

}

TCanvas *plotNhits(){
  TCanvas *c55 = new TCanvas("c55", "c5", 50,50,800,800);
  c55->Divide(3,4, 0.01, 0.01, 0);

  c55->cd(1);
  hnhits1->Draw();
  c55->cd(2);
  hngoodhits1->Draw();
  c55->cd(3);
  hngoodTDChits1->Draw();
  c55->cd(4);
  hnhits2->Draw();
  c55->cd(5);
  hngoodhits2->Draw();
  c55->cd(6);
  hngoodTDChits2->Draw();
  
  c55->cd(7);
  hnpaddles->Draw();
  c55->cd(8);
  hngoodpaddles->Draw();
  c55->cd(9);
  hngoodTDCpaddles->Draw();
  
  c55->cd(10);
  hnhits_ev->Draw();
  c55->cd(11);
  hngoodhits_ev->Draw();
  c55->cd(12);
  hngoodTDChits_ev->Draw();

  return c55;

}

TCanvas *plotTDC2d(){

  h2d_RawLE  = new TH2F("h2d_RawLE","Raw LE vs PMT", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTE  = new TH2F("h2d_RawTE","Raw TE vs PMT", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTot = new TH2F("h2d_RawTot","Raw TOT vs PMT", NTotBins,TotBinLow,TotBinHigh,nTdc+1,0,nTdc+1);

  // TH2D *h2d_RawLE = new TH2D("h2d_RawLE", "Raw LE vs PMT", 400, 0, 200, 2700, 0, 2700);
  // TH2D *h2d_RawTE = new TH2D("h2d_RawTE", "Raw TE vs PMT", 400, 0, 200, 2700, 0, 2700);
  // TH2D *h2d_RawTot = new TH2D("h2d_RawTot", "Raw Tot vs PMT", 400, 0, 200, 2700, 0, 2700);

  for (size_t evt = 0; evt < vRawLe.size(); evt++) {
    for (size_t hit = 0; hit < vRawLe[evt].size(); hit++) {
      h2d_RawLE->Fill(vRawLe[evt][hit], vRawID[evt][hit]);
      h2d_RawTE->Fill(vRawTe[evt][hit], vRawID[evt][hit]);
      h2d_RawTot->Fill(vRawTot[evt][hit], vRawID[evt][hit]);
    }
  }



  TCanvas *c6 = new TCanvas("c6", "c6", 50,50,800,800);
  c6->Divide(2,2, 0.01, 0.01, 0);

  c6->cd(1);
  h2d_RawLE->Draw();
  c6->cd(2);
  h2d_RawTE->Draw();
  c6->cd(3);
  h2d_RawTot->Draw();
  c6->cd(4);
  h2d_Mult->Draw();

  return c6;

}

auto *plotEECalCDet() {

  TCanvas *c1717 = new TCanvas("c1717", "c1717", 50,50,800,800);
  c1717->Divide(1,2, 0.01, 0.01, 0);

  c1717->cd(1);
  hEECalCDet1->Draw();
  c1717->cd(2);
  hEECalCDet2->Draw();

  return c1717;
}

auto plotXYZ(){

  hHitX = new TH1F("HitXposition","HitXPosition",1000,-2.0,2.0);
  hHitY = new TH1F("HitYposition","HitYPosition",200,-0.5,0.5);
  hHitZ = new TH1F("HitZposition","HitZPosition",200,5.5,6.0);
   size_t N = vCDetGoodX.size();
    for (size_t ev = 0; ev < N; ev++){
      for (size_t hit = 0; hit < vCDetGoodX[ev].size(); hit++) {
      //if (vRefRawTot[i] >= cutTotMin && vRefRawTot[i] <= cutTotMax) {
        hHitX->Fill(vCDetGoodX[ev][hit]);
        hHitY->Fill(vCDetGoodY[ev][hit]);
        hHitZ->Fill(vCDetGoodZ[ev][hit]);
      //}
      }
    }

   TCanvas *c7 = new TCanvas("c7", "c7", 800,800);
   c7->Draw();
   TPad *p1 = new TPad("p1","p1",0.05,0.0,0.45,1.0);
   p1->Draw();
   p1->Divide(1,3);

   p1->cd(1);
   gPad->SetLogy();
   hHitX->Draw();

   p1->cd(2);
   gPad->SetLogy();
   hHitY->Draw();
   
   p1->cd(3);
   gPad->SetLogy();
   hHitZ->Draw();

   c7->cd(0);
   TPad *p2 = new TPad("p1","p1",0.55,0.0,0.95,1.0);
   p2->Draw();
   p2->Divide(2,1);

   p2->cd(1);
   gPad->SetLogz();
   hHitXY1->Draw("colz");
   p2->cd(2);
   gPad->SetLogz();
   hHitXY2->Draw("colz");


  return c7;

}

/*auto plotXYECalCDet(){

   TCanvas *c8 = new TCanvas("c8", "c7", 1200,1200);
   c8->Divide(3,3);

   c8->cd(1);
   gPad->SetLogz();
   hXECalCDet1->Draw("colz");

   c8->cd(2);
   gPad->SetLogz();
   hXECalCDet2->Draw("colz");
   
   c8->cd(3);
   hXECal->Draw();
   
   c8->cd(4);
   hYECal->Draw();
   
   c8->cd(5);
    // Define the Gaussian + constant background function
    TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

    // Set initial parameters:
    // [0] amplitude, [1] mean, [2] sigma, [3] constant background
    fitFunc->SetParameters(hXDiffECalCDet1->GetMaximum(), 0.02, 0.01, hXDiffECalCDet1->GetMinimum());
    fitFunc->SetParNames("Amplitude", "Mean", "Sigma", "Background");

    // Optional: set limits on the parameters if needed
    fitFunc->SetParLimits(1, -0.05, 0.05);  // constrain the mean near 0.02
    fitFunc->SetParLimits(2, 0.001, 0.05);  // positive sigma

    // Fit the histogram
    hXDiffECalCDet1->Fit(fitFunc, "R");  // "R" = use function range only

    // Draw the result
    hXDiffECalCDet1->Draw();
    fitFunc->Draw("same");
    // Extract Gaussian parameters
    double A = fitFunc->GetParameter(0); // Amplitude
    double mu = fitFunc->GetParameter(1); // Mean
    double sigma = fitFunc->GetParameter(2); // Sigma
    double bg = fitFunc->GetParameter(3); // Background level

    // Define integration limits
    double x_min = mu - 3*sigma;
    double x_max = mu + 3*sigma;

    // Signal: integral of the Gaussian part only over ±3σ
    TF1* gausOnly = new TF1("gausOnly", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
    gausOnly->SetParameters(A, mu, sigma);
    double signal = gausOnly->Integral(x_min, x_max);

    // Noise: integral of background over same range
    double noise = bg * (x_max - x_min);

    // Compute signal-to-noise ratio
    double snr = (noise > 0) ? signal / noise : 0;

    std::cout << "Signal (Gaussian, ±3σ): " << signal << std::endl;
    std::cout << "Noise (Background, ±3σ): " << noise << std::endl;
    std::cout << "Signal-to-Noise Ratio: " << snr << std::endl;
   

   c8->cd(6);
    // Define the Gaussian + constant background function
    TF1* fitFunc2 = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

    // Set initial parameters:
    // [0] amplitude, [1] mean, [2] sigma, [3] constant background
    fitFunc2->SetParameters(hXDiffECalCDet2->GetMaximum(), 0.02, 0.01, hXDiffECalCDet2->GetMinimum());
    fitFunc2->SetParNames("Amplitude", "Mean", "Sigma", "Background");

    // Optional: set limits on the parameters if needed
    fitFunc2->SetParLimits(1, -0.05, 0.05);  // constrain the mean near 0.02
    fitFunc2->SetParLimits(2, 0.001, 0.05);  // positive sigma

    // Fit the histogram
    hXDiffECalCDet2->Fit(fitFunc2, "R");  // "R" = use function range only

    // Draw the result
    hXDiffECalCDet2->Draw();
    fitFunc2->Draw("same");

    // Extract Gaussian parameters
    double A2 = fitFunc2->GetParameter(0); // Amplitude
    double mu2 = fitFunc2->GetParameter(1); // Mean
    double sigma2 = fitFunc2->GetParameter(2); // Sigma
    double bg2 = fitFunc2->GetParameter(3); // Background level

    // Define integration limits
    double x_min2 = mu2 - 3*sigma2;
    double x_max2 = mu2 + 3*sigma2;

    // Signal: integral of the Gaussian part only over ±3σ
    TF1* gausOnly2 = new TF1("gausOnly2", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
    gausOnly2->SetParameters(A2, mu2, sigma2);
    double signal2 = gausOnly2->Integral(x_min2, x_max2);

    // Noise: integral of background over same range
    double noise2 = bg2 * (x_max2 - x_min2);

    // Compute signal-to-noise ratio
    double snr2 = (noise2 > 0) ? signal2 / noise2 : 0;

    std::cout << "Signal (Gaussian, ±3σ): " << signal2 << std::endl;
    std::cout << "Noise (Background, ±3σ): " << noise2 << std::endl;
    std::cout << "Signal-to-Noise Ratio: " << snr2 << std::endl;

   c8->cd(7);
   hYECalCDet1->Draw();
   
   c8->cd(8);
   hYECalCDet2->Draw();

   
   c8->cd(9);
   hXYECal->Draw("colz");
   

  return c8;
}*/

auto plotXYECalCDet(){

   // 4x3 layout so we can add the two new "min-hit" histograms cleanly
   TCanvas *c8 = new TCanvas("c8", "plotXYECalCDet", 1600,1200);
   c8->Divide(4,3);

   // ---------------- Row 1: XECal vs XCDet (cut-based and min-hit-based) ----------------

   c8->cd(1);
   gPad->SetLogz();
   hXECalCDet1->SetMinimum(10);
   hXECalCDet1->Draw("colz");

   c8->cd(2);
   gPad->SetLogz();
   hXECalCDet2->SetMinimum(10);
   hXECalCDet2->Draw("colz");

   c8->cd(3);
   gPad->SetLogz();
   hXECalCDet1_min->SetMinimum(10);
   hXECalCDet1_min->Draw("colz");

   c8->cd(4);
   gPad->SetLogz();
   hXECalCDet2_min->SetMinimum(10);
   hXECalCDet2_min->Draw("colz");

   // ---------------- Row 2: ECal X/Y + Xdiff fits ----------------

   c8->cd(5);
   hXECal->Draw();

   c8->cd(6);
   hYECal->Draw();

   // ---- Fit XDiff layer 1 (same code as you had, moved to pad 7)
   c8->cd(7);
   {
     TF1* fitFunc = new TF1("fitFunc1", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

     fitFunc->SetParameters(hXDiffECalCDet1->GetMaximum(), 0.02, 0.01, hXDiffECalCDet1->GetMinimum());
     fitFunc->SetParNames("Amplitude", "Mean", "Sigma", "Background");

     fitFunc->SetParLimits(1, -0.05, 0.05);
     fitFunc->SetParLimits(2, 0.001, 0.05);

     hXDiffECalCDet1->Fit(fitFunc, "R");
     hXDiffECalCDet1->Draw();
     fitFunc->Draw("same");

     double A     = fitFunc->GetParameter(0);
     double mu    = fitFunc->GetParameter(1);
     double sigma = fitFunc->GetParameter(2);
     double bg    = fitFunc->GetParameter(3);

     double x_min = mu - 3*sigma;
     double x_max = mu + 3*sigma;

     TF1* gausOnly = new TF1("gausOnly1", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
     gausOnly->SetParameters(A, mu, sigma);

     double signal = gausOnly->Integral(x_min, x_max);
     double noise  = bg * (x_max - x_min);
     double snr    = (noise > 0) ? signal / noise : 0;

     std::cout << "Layer1: Signal (Gaussian, ±3σ): " << signal << std::endl;
     std::cout << "Layer1: Noise  (Background, ±3σ): " << noise  << std::endl;
     std::cout << "Layer1: Signal-to-Noise Ratio: " << snr << std::endl;
   }

   // ---- Fit XDiff layer 2 (same code as you had, moved to pad 8)
   c8->cd(8);
   {
     TF1* fitFunc2 = new TF1("fitFunc2", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

     fitFunc2->SetParameters(hXDiffECalCDet2->GetMaximum(), 0.02, 0.01, hXDiffECalCDet2->GetMinimum());
     fitFunc2->SetParNames("Amplitude", "Mean", "Sigma", "Background");

     fitFunc2->SetParLimits(1, -0.05, 0.05);
     fitFunc2->SetParLimits(2, 0.001, 0.05);

     hXDiffECalCDet2->Fit(fitFunc2, "R");
     hXDiffECalCDet2->Draw();
     fitFunc2->Draw("same");

     double A2     = fitFunc2->GetParameter(0);
     double mu2    = fitFunc2->GetParameter(1);
     double sigma2 = fitFunc2->GetParameter(2);
     double bg2    = fitFunc2->GetParameter(3);

     double x_min2 = mu2 - 3*sigma2;
     double x_max2 = mu2 + 3*sigma2;

     // IMPORTANT FIX: use (x_min2, x_max2) here (your current file uses x_min/x_max by accident)
     TF1* gausOnly2 = new TF1("gausOnly2", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min2, x_max2);
     gausOnly2->SetParameters(A2, mu2, sigma2);

     double signal2 = gausOnly2->Integral(x_min2, x_max2);
     double noise2  = bg2 * (x_max2 - x_min2);
     double snr2    = (noise2 > 0) ? signal2 / noise2 : 0;

     std::cout << "Layer2: Signal (Gaussian, ±3σ): " << signal2 << std::endl;
     std::cout << "Layer2: Noise  (Background, ±3σ): " << noise2  << std::endl;
     std::cout << "Layer2: Signal-to-Noise Ratio: " << snr2 << std::endl;
   }

   // ---------------- Row 3: YECal vs YCDet and XY ECal ----------------

   c8->cd(9);
   hYECalCDet1->Draw();

   c8->cd(10);
   hYECalCDet2->Draw();

   c8->cd(11);
   gPad->SetLogz();
   hXYECal->Draw("colz");

   // pad 12 left intentionally empty for now (room for future additions)

   return c8;
}


//------------Some routines when Ben getting familiar with branches
auto plotDpp(int nbins = 100, double xmin = -0.1, double xmax =  0.1)
{
    TH1D* h = new TH1D("hHeep_dpp", "heep_dpp;#delta p/p;Counts", nbins, xmin, xmax);
    for (double x : vheep_dpp) if (std::isfinite(x)) h->Fill(x);
    auto* c9 = new TCanvas(Form("c_%s","heep_dpp"), "heep_dpp", 900, 650); h->Draw();
    return c9;
}

//auto find hist range
static std::pair<double,double> MinMaxFlat(const std::vector<std::vector<double>>& vv){
  double mn =  std::numeric_limits<double>::infinity();
  double mx = -std::numeric_limits<double>::infinity();
  for (const auto& v : vv){
    for (double x : v){
      if (std::isfinite(x)) {
        if (x < mn) mn = x;
        if (x > mx) mx = x;
      }
    }
  }
  if (!std::isfinite(mn) || !std::isfinite(mx)) { // empty or all non-finite
    mn = 0.0; mx = 1.0;
  }
  if (mn == mx) { // collapse -> pad a bit
    mn -= 0.5; mx += 0.5;
  }
  return {mn, mx};
}

auto plotECalClusX(double xmin = -1.5, double xmax = 1.5, bool log=true)
{
  int nbin = std::ceil((xmax - xmin)/0.0425);
  TH1D* h_ECal_clus_x = new TH1D("h_ECal_clus_x", "hECalClusX;Clus x (m);Counts", nbin,-1.5,1.5);

  //fill hist from vectors
  for (const auto& vec : v_ECal_clus_x) {
    for (double val : vec){
      h_ECal_clus_x->Fill(val);
    }
  }
  TCanvas* c_ECal_clus_x = new TCanvas("c_ECal_clus_x", "ECal Cluster x",800,600);
  if (log==true){
    c_ECal_clus_x->SetLogy();
  }
  h_ECal_clus_x->Draw("HIST");
}
auto plotECalNclus()
{
  int nclus_max = *std::max_element(v_ECal_nclus.begin(), v_ECal_nclus.end());
  int nclus_min = *std::min_element(v_ECal_nclus.begin(),v_ECal_nclus.end());
  TH1D* h_ECal_nclus = new TH1D("h_ECal_nclus", "hECalNclus;Nclus;Counts", nclus_max+1, nclus_min-0.5, nclus_max + 0.5);
  for (auto val : v_ECal_nclus){
    h_ECal_nclus->Fill(val);
  }
  TCanvas* c_ECal_nclus = new TCanvas("c_ECal_nclus", "ECal Cluster Count",800,600);
  h_ECal_nclus->Draw("HIST");
}

auto plotECalClusE(int binlow = 0, int binhigh=12)
{
  TH1D* h_ECal_clus_e = new TH1D("h_ECal_clus_e", "hECalCluse;Clus E (GeV);Counts", (binhigh+binlow)*100,binlow-0.5,binhigh+0.5);

  //fill hist from vectors
  for (const auto& vec : v_ECal_clus_e) {
    for (double val : vec){
      h_ECal_clus_e->Fill(val);
    }
  }
  TCanvas* c_ECal_clus_e = new TCanvas("c_ECal_clus_e", "ECal Cluster e",800,600);
  h_ECal_clus_e->Draw("HIST");
}
auto plotECalClusAdcTime()
{
  double adctime_min = std::numeric_limits<double>::max();
  double adctime_max = std::numeric_limits<double>::lowest();

  for (const auto& subvec : v_ECal_clus_adctime) {
    if (subvec.empty()) continue; // skip events with no clusters

    // find min and max within this event
    double local_min = *std::min_element(subvec.begin(), subvec.end());
    double local_max = *std::max_element(subvec.begin(), subvec.end());

    // update global range
    if (local_min < adctime_min) adctime_min = local_min;
    if (local_max > adctime_max) adctime_max = local_max;
  }

  TH1D* h_ECal_clus_adctime = new TH1D("h_ECal_clus_adctime", "hECalClusAdctime;Clus adctime;Counts", 1000,adctime_min-10,adctime_max+10);

  //fill hist from vectors
  for (const auto& vec : v_ECal_clus_adctime) {
    for (double val : vec){
      h_ECal_clus_adctime->Fill(val);
    }
  }
  TCanvas* c_ECal_clus_adctime = new TCanvas("c_ECal_clus_adctime", "ECal Cluster adctime",800,600);
  h_ECal_clus_adctime->Draw("HIST");
}
