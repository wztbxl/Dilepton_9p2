//-- define constant --//
const Double_t Mmuon = 0.105658367;//Gev
const double PI = TMath::Pi();

//-- cut --//
Double_t mVzCut  = 100.;
Double_t mVxyCut = 2.;
Double_t mTpcMinPtCut = 1.3;
Double_t mTpcMaxPtCut = 1e4;
Double_t mMuEtaCut = 0.8;
Double_t mNHitsFitMaxCut = 0.52;
Int_t mNHitsFitCut = 15;
Int_t mNHitsDedxCut = 10;
Double_t mDcaCut = 3.;

Double_t mMinNSigmaPiCut = -2.0;
Double_t mMaxNSigmaPiCut = 3.;

Double_t mMuDYSigCut = 3.;
Double_t mMuDZSigCut = 3.;
Double_t mMuMindTofCut = -1.0;
Double_t mMuMaxdTofCut = 1.0;

Double_t mPairYCut = 0.5;

//-- global variable --//
LorentzVec muPlus;
LorentzVec muMinus;
DoubleVec muPlusInfo;
DoubleVec muMinusInfo;

//-- define histograms --//
TH1D* hEvent;
TH2D* hEtaVsPhi;
TH2D* hHitMap;
TH2D* hHitEtaVsPhi;

THnSparse        *mhULMPtCosPhi;
THnSparse        *mhLSMPtCosPhi;
THnSparse        *mhULMPtCosPhiCS;
THnSparse        *mhLSMPtCosPhiCS;

TH3D* mhULPtCos;
TH3D* mhLSPtCos;
TH3D* mhULPtPhi;
TH3D* mhLSPtPhi;

TH3D* hMuonPtEtaPhi;

//-- LorentzVector --//
TLorentzVector muon1, muon2, Jpsi;
