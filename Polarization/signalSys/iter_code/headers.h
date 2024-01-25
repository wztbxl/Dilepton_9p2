#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "sys/types.h"
#include "dirent.h"
#include "math.h"
#include "string.h"

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "EmbedTree.h"
#include "THnSparse.h"
using namespace std;
#endif

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Float_t> FloatVec;
typedef vector<TLorentzVector> LorentzVec;
#else
typedef vector<Float_t, allocator<Float_t>> FloatVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
#endif

const Double_t PI = TMath::Pi();
const Double_t muMass = 0.105658367;

//-----cuts
const Double_t mSmearPtRes = 0.012;
const Double_t mMaxVz = 100.;
const Double_t mMaxVr = 2.;

const Double_t mMinTrkPt = 1.;
const Double_t mMaxTrkPt = 20.;
Int_t mMinNHitsDedx = 10;
Int_t mMinNHitsFit = 15;
const Double_t mMinFitHitsFaction = 0.52;
Double_t mMaxDca = 3.;
const Double_t mMinTrkEta = -0.8;
const Double_t mMaxTrkEta = 0.8;
const Double_t mMinTrkPhi = -2*PI;
const Double_t mMaxTrkPhi = 2*PI;

const Double_t mMuPtCut = 1.3;
Double_t mMinNSigmaPiCut = -2.;
const Double_t mMaxNSigmaPiCut = 3.;
const Double_t mMuDYSigCut = 3.;
const Double_t mMuDZSigCut = 3.;
const Double_t mMuMindTofCut = -1e+4;
const Double_t mMuMaxdTofCut = 1.;
const Double_t mPairYCut = 0.5;

int nIterate = 0;
int Frame = 0;
int nPt = 0;
Double_t lamThe = 0;
Double_t lamPhi = 0;
const char *frameName[2] = {"HX", "CS"};

TRandom3 *myRandom;

//---claim--function
void bookHistograms();
void writeHistograms(char* outFile);
Bool_t passEvent(EmbedTree *evt);
Bool_t passTrack(Double_t pt, Int_t nHitsDedx, Int_t nHitsFit, Double_t nHitsFrac, Double_t dca, Double_t eta, Double_t phi);
Double_t weightFunc(double pT);
Double_t weightFunc(Double_t pT, TH1D* hShape);
Double_t smearPt(double pT);
Bool_t checkMuCandidate(Double_t pt, Double_t nSigmaPion, Double_t dy, Double_t dz, Double_t dtof, Int_t matchFlag);
Bool_t checkMuDeltaY(Double_t pt, Double_t dy);
Bool_t checkMuDeltaZ(Double_t pt, Double_t dz);
Bool_t checkMuDeltaTof(Double_t dtof);
Bool_t calMtdTriggerEff(Double_t pt, TH1D* hEff);
Bool_t calMtdResponseEff(Double_t pt, TH1D* hCosmic, TH1D*hEmb );
Bool_t isActiveModule(Int_t mRunId, Int_t backleg, Int_t module);
void calPolarization(TLorentzVector iVec,TLorentzVector vec, THnSparse * hn, THnSparse *hnCS, double w, double pTheta, double pPhi );
Bool_t  calTpcRespEff(Double_t pt, Double_t eta, Double_t phi, TH2D* hEff);
Double_t weightPola( double ptheta, double pphi, double costheta, double phi );

// book histograms
TH1D* hEvent;
TH2D* hEtaVsPhi;
TH2D* hHitMap;
TH2D* hHitEtaVsPhi;

THnSparse        *mhMcMPtCosPhi;
THnSparse        *mhMcMPtCosPhiCS;
THnSparse        *mhRcMPtCosPhi;
THnSparse        *mhRcMPtCosPhiCS;

TH1D* hJpsiSpec[4][2];
TH1D* hNEmbedJpsi;

TLorentzVector mcMuon1, mcMuon2, mcJpsi;
TLorentzVector rcMuon1, rcMuon2, rcJpsi;

