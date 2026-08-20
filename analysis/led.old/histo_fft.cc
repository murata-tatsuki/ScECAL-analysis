#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1S.h"
#include "TH2D.h"
#include "TH2S.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include "TCanvas.h"
// #include "langaus_newConvoluted.C"
#include "EBUdecode.h"
#include "EBUdecode.cxx"
#include "fft_fullrange_tempcorr_plusNphoto_v5.C"

using namespace std;

// ssa後にcutをかけて1p.e.のgainを求めるマクロ
// cut有り無しを両方表示するマクロ
// fftを用いて


const bool fitting = true;

double fixTemp = 25;



int main(int argc, char* argv[])
{ 
  const TString& runFile = argv[1];
  
  analyzeOneRun_full_(runFile);
}