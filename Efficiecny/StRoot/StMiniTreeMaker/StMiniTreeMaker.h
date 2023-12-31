#ifndef STMINITREEMAKER_HH
#define STMINITREEMAKER_HH

/***************************************************************************
 *
 * $Id: StMiniTreeMaker.h 2015/04/09  Exp $ 
 * StMiniTreeMaker - class to produce miniTree for mtd related analysis
 * Author: Shuai Yang
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/

#include "StMaker.h"
#include "StTreeStructure.h"

#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TRandom3.h"


#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#ifndef ST_NO_NAMESPACES
using std::vector;
using std::map;
#endif

class TH1D;
class TH2D;
class TH1I;
class TH3F;
class TH3D;
class TH2F;
class TString;
class TFile;

class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StRefMultCorr;

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Int_t> IntVec;
typedef vector<Double_t> DoubleVec;
#else
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
#endif

class StMiniTreeMaker : public StMaker
{
  public:
	StMiniTreeMaker(const Char_t *name = "StMiniTreeMaker");
	~StMiniTreeMaker();

	Int_t Init();
	Int_t InitRun(const Int_t runNumber);
	Int_t Make();
	Int_t Finish();

	void setTriggerIDs(const IntVec triggerids);
	void setUseDefaultVtx(const Bool_t flag);
	void setMaxVtxR(const Double_t max);
	void setMaxVtxZ(const Double_t max);
	void setMaxVzDiff(const Double_t max);
	void setMinTrackPt(const Double_t min);
	void setMaxTrackEta(const Double_t max);
	void setMinNHitsFit(const Int_t min);
	void setMinNHitsFitRatio(const Double_t min);
	void setMinNHitsDedx(const Int_t min);
	void setMaxDca(const Double_t max);
	void setMaxnSigmaE(const Double_t max);
	void setMaxBeta2TOF(const Double_t max);
	void setmPhotonicEMassCut(const double max);
	void setmElectronCutMode(const double mode);
	void setFillHisto(const Bool_t fill);
	void setFillTree(const Bool_t fill);
	void setOutFileName(const TString name);
	void setStreamName(const TString name);
	void setPrintMemory(const Bool_t pMem);
	void setPrintCpu(const Bool_t pCpu);
	void setPrintConfig(const Bool_t print);

  protected:
	void printConfig();
	void bookHistos();
	void writeHistos();
	Bool_t processPicoEvent();
	Bool_t isValidTrack(StPicoTrack *pTrack, TVector3 vtxPos) const;
	bool isElectron(StPicoTrack *pTrack);
	bool isElectronBetaCut(StPicoTrack *pTrack);
	bool isElectronSigmaECut(StPicoTrack *pTrack);
	double phiVangle(TLorentzVector e1, TLorentzVector e2, int q1, int q2);
	bool isPiKP_masscut(double msquare);

  private:

	ifstream indata_001;//for load bad run
	map<Int_t,Int_t> mBadRunId_001;

	StPicoDstMaker *mPicoDstMaker;
	StPicoDst *mPicoDst;
	StRefMultCorr *refMultCorr; //decide centrality

	const Double_t PI = TMath::Pi(); //const PI for draw histo

	Bool_t mPrintMemory;		// Flag to print out memory usage
	Bool_t mPrintCpu;			// Flag to print out CPU usage
	Bool_t mPrintConfig;		// Flag to print out task configuration
	TString mStreamName;		// Data stream name
	Bool_t mDefaultVtx;			// Use Default Vertex
	Bool_t mSelectVtxRank;		// Vertex ranking > 0
	Double_t mMaxVtxR;			// Maximum vertex r
	Double_t mMaxVtxZ;			// Maximum vertex z
	Double_t mMaxVzDiff;		// Maximum VpdVz-TpcVz
	Double_t mMinTrkPt;			// Minimum track pt
	Double_t mMaxTrkEta;		// Maximum track eta
	Int_t mMinNHitsFit;			// Minimum number of hits used for track fit
	Double_t mMinNHitsFitRatio; // Minimum ratio of hits used for track fit
	Int_t mMinNHitsDedx;		// Minimum number of hits used for de/dx
	Double_t mMaxDca;			// Maximum track dca
	Double_t mMaxnSigmaE;		// Maximum nSigmaE cut
	Double_t mMinnSigmaE;		// Maximum nSigmaE cut
	Double_t mMaxBeta2TOF;		// Maximum |1-1./beta| for TpcE
	Int_t mElectronCutMode;      //Switch the mode of Eff selection
	double mPhotonicEMassCut;   // Photonic elextron mass cut

	Bool_t mFillHisto;	// Flag of fill the histogram
	TFile *fOutFile;	  // Output file
	TString mOutFileName; // Name of the output file
	StEvtData mEvtData;
	TTree *mEvtTree; // Pointer to the event tree
	const double Melectron = 0.00051099907;
	const int MaxNElectron = 10000;
	TF1* PileupUplimit;
    TF1* PileupLowlimit;
    TF1* PileupLimit;

	IntVec mTriggerIDs;

	DoubleVec vEtaPlusCosPart;
	DoubleVec vEtaPlusSinPart;
	DoubleVec vEtaPlusPtWeight;
	DoubleVec vEtaMinusCosPart;
	DoubleVec vEtaMinusSinPart;
	DoubleVec vEtaMinusPtWeight;

	//define histograms ongoing...
	TH1D *hEvent;
	TH2D *hVtxYvsVtxX;
	TH2D *hVPDVzvsTPCVz;
	TH1D *hVzDiff;
	//some histo to check refMult
	TH1D* hRefMult;
	TH1D* hCentrality9;
	TH2D* hRefMultvsnTOFMatch;
	TH2D* hMsquraevsRefMult;
	TH3D* hRefMultvsnTOFMatchvsVz;
	TH2D* hRefMultvsnChargeParticle;
	TH2D* hnTOFMatchvsnChargePartile;

	TH2D *hdEdxvsP;
	TH2D *hdNdxvsP;
	TH2D *hnSigEvsP;
	TH2D *hBetavsP;
	TH2D* hBetavsPwonSigmaE;
	TH2D *hnSigEvsPWTOF;
	TH2D *hULMvsphiV;
	TH2D *hLPosMvsphiV;
	TH2D *hLNegMvsphiV;

	// define histogram for check untag and likesign method 
	TH1D *htrkMassdistrbution;
	TH2D* hULMvsPtLarge;
	TH2D *hLMPosvsPt;
	TH2D *hLMNegvsPt;
	TH2D *hULDcavsPt;
	TH2D *hLPosDcavsPt;
	TH2D *hLNegDcavsPt;

	

	//histo for cal eff
	TH3F *hDenPiPlusTofEff;
	TH3F *hNumPiPlusTofEff;
	TH3F *hDenPiMinusTofEff;
	TH3F *hNumPiMinusTofEff;
	TH3F *hDenPiPlusTofEffCen[9];
	TH3F *hNumPiPlusTofEffCen[9];
	TH3F *hDenPiMinusTofEffCen[9];
	TH3F *hNumPiMinusTofEffCen[9];
	TH2F *hULMvsPt;
	TH2F *hLPosMvsPt;
	TH2F *hLNegMvsPt;
	//beta Eff
	TH2F *hPEPlusBetavsP;
	TH2F *hPEMinusBetavsP;
	TH2F *hPEPlusBetavsPt;
	TH2F *hPEMinusBetavsPt;
	TH1D* hBeta;
	TH1D* hTofLocalY;
	TH1D* hTofMatchFlag;
	TH1D* hTofCellID;

	TH3D* hPureElectronNSigmaEvsPtvsCen;	
	TH3D* hPurePositronNSigmaEvsPtvsCen;
	TH3D* hPureElectronNSigmaEvsPvsCen;
	TH3D* hPurePositronNSigmaEvsPvsCen;
	TH3D* hPureElectronNSigmaEvsPhivsCen;
	TH3D* hPurePositronNSigmaEvsPhivsCen;
	TH3F *hPEPlusBetavsPCen;
	TH3F *hPEMinusBetavsPCen;
	TH3F *hPEPlusBetavsPtCen;
	TH3F *hPEMinusBetavsPtCen;
	TH3F *hDenPEPlusTofEffCen[9];
	TH3F *hNumPEPlusTofEffCen[9];
	TH3F *hDenPEMinusTofEffCen[9];
	TH3F *hNumPEMinusTofEffCen[9]; 


	// for check something that in the start of analysis
	TH1D *hMsquare;
	TH2D *hMsquarevsP;
	TH2D *hPurePionNSigmaEvsP;
	TH2D *hMergePionNSigmaEvsP;
	TH2D *hPureKaonNSigmaEvsP;
	TH2D *hPureProtonNSigmaEvsP;
	TH3D *hPurePionNSigmaEvsPCen;
	TH3D *hMergePionNSigmaEvsPCen;
	TH3D *hPureKaonNSigmaEvsPCen;
	TH3D *hPureProtonNSigmaEvsPCen;
	TH3D *hPureElectronNSigmaEvsPCen;
	
	TH2D *hPureProtonPhivsEta;
	TH2D *hPureElectronNSigmaEvsP;
	TH2D *hPEElectronnSigmaEvsP;
	TH2D *hPEPositronnSigmaEvsP;
	TH2D *hPurePionNSigmaEvsPt;
	TH2D *hMergePionNSigmaEvsPt;
	TH2D *hPureKaonNSigmaEvsPt;
	TH2D *hPureProtonNSigmaEvsPt;
	TH2D *hPureElectronNSigmaEvsPt;
	TH2D *hPurePositronNSigmaEvsPt;
	TH2D *hSignaEvsP;
	TH2D *hSigmaEvsPwithNSigmaE;
	TH2D *hSigmaEvsPwithNSigEandBeta;

	TH2D* hPurePionSigmaPionvsP;

	//check the vertical band 
	TH2D* hEtavsPhi_vband;
	TH2D* hTOFEtavsPhi; // electron TOF eta vs Phi
	TH2D* hTOFEtavsPhi_cellID;
	TH2D* hTOFEtavsPhi_vband;
	TH1D* hTOFCellID_vband;
	TH2D* nSigmaE_vband;
	TH2D* hBetavsP_Pion;
	TH2D* hBetavsP_Kaon;
	TH2D* hBetavsP_Proton;
	TH1D* hLocalY_vband;
	TH2D* hEtavsPhi_pT1;
	TH2D* hnHitsFitvsP_vBand;
	TH2D* hnHitsdEdxvsP_vBand;
	TH2D* hnHitsratiovsP_vBand;
	TH2D* hnHitsPossvsnHitsMax;


	ClassDef(StMiniTreeMaker, 1)
};

inline void StMiniTreeMaker::setTriggerIDs(const IntVec triggerids) { mTriggerIDs = triggerids; }
inline void StMiniTreeMaker::setUseDefaultVtx(const Bool_t flag) { mDefaultVtx = flag; }
inline void StMiniTreeMaker::setMaxVtxR(const Double_t max) { mMaxVtxR = max; }
inline void StMiniTreeMaker::setMaxVtxZ(const Double_t max) { mMaxVtxZ = max; }
inline void StMiniTreeMaker::setMaxVzDiff(const Double_t max) { mMaxVzDiff = max; }
inline void StMiniTreeMaker::setMinTrackPt(const Double_t min) { mMinTrkPt = min; }
inline void StMiniTreeMaker::setMaxTrackEta(const Double_t max) { mMaxTrkEta = max; }
inline void StMiniTreeMaker::setMinNHitsFit(const Int_t min) { mMinNHitsFit = min; }
inline void StMiniTreeMaker::setMinNHitsFitRatio(const Double_t min) { mMinNHitsFitRatio = min; }
inline void StMiniTreeMaker::setMinNHitsDedx(const Int_t min) { mMinNHitsDedx = min; }
inline void StMiniTreeMaker::setMaxDca(const Double_t max) { mMaxDca = max; }
inline void StMiniTreeMaker::setMaxnSigmaE(const Double_t max) { mMaxnSigmaE = max; }
inline void StMiniTreeMaker::setMaxBeta2TOF(const Double_t max) { mMaxBeta2TOF = max; }
inline void StMiniTreeMaker::setmPhotonicEMassCut(const Double_t max) { mPhotonicEMassCut = max; }
inline void StMiniTreeMaker::setmElectronCutMode(const Double_t mode) { mElectronCutMode = mode; }
inline void StMiniTreeMaker::setFillHisto(const Bool_t fill) { mFillHisto = fill; }
inline void StMiniTreeMaker::setOutFileName(const TString name) { mOutFileName = name; }
inline void StMiniTreeMaker::setStreamName(const TString name) { mStreamName = name; }
inline void StMiniTreeMaker::setPrintMemory(const Bool_t pMem) { mPrintMemory = pMem; }
inline void StMiniTreeMaker::setPrintCpu(const Bool_t pCpu) { mPrintCpu = pCpu; }
inline void StMiniTreeMaker::setPrintConfig(const Bool_t print) { mPrintConfig = print; }
#endif
