#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include "math.h"
#include "sys/types.h"

#ifndef __CINT__ 
#include "TROOT.h"
#include "TChain.h"
#include "TClass.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TAttFill.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TSpectrum.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TFitResult.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TArrow.h"
#include "TArrayL.h"
#include "TArrayF.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TPDF.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"
using namespace std;
#endif

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Bool_t> BoolVec;
typedef vector<Int_t> IntVec;
typedef vector<Float_t> FloatVec;
typedef vector<Double_t> DoubleVec;
typedef vector<TLorentzVector> LorentzVec;
#else
typedef vector<Bool_t, allocator<Bool_t>> BoolVec;
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<Float_t, allocator<Float_t>> FloatVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
#endif


