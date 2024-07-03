//////////////////////////////////////////////////////////
// This class is used for meson decay
// ver 1.0 2011/04/27 huangbc
//
//////////////////////////////////////////////////////////

#ifndef VECMESON_H
#define VECMESON_H

#include "TROOT.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "MesonConstant.h"

#define nSmearFac 3

class VecMeson
{
	public:
		VecMeson(ParticleTypes particle=omega,DecayMode dmode=twobody);
		~VecMeson();
		Int_t Init();

		void 		SetUseCocktailInput(Bool_t kFlag);
		void        SetUseScaleEff(Bool_t kFlag);
		void 		SetInputPtSpectra(Int_t kFlag);
		void 		SetNumberOfTracks(Int_t nTrks);
		void        SetCentralityIdx(Int_t idx);
		void        SetRapidityRange(Double_t min, Double_t max);
		void 		SetPtSmearPar(Double_t *par);
		Double_t 	GetSmear(Double_t pT);
		Double_t 	GetSmear2(Double_t pT);
		TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
		TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass);
		void Polarization(int icharge,int jcharge,TLorentzVector ivector,TLorentzVector jvector);
		void GetPtPhiCentBin(TLorentzVector pair,TLorentzVector Positron, int _mCentrality,float eventphi,int &ptindex,int &yindex,int &phiindex,int &CentIndex,double &costhe,Bool_t tangent, Int_t Flag );


		Double_t	SampleMesonMass();
		Double_t 	EvalEff3D(TLorentzVector electron,int ipttpc,int ietatpc,int iphitpc,int ipttof,int ietatof,int iphitof,int charge);
		Double_t 	EvalScaleEff3D(TLorentzVector electron,int ipttpc,int ietatpc,int iphitpc,int ipttof,int ietatof,int iphitof,int charge);

		void 		Rot(Double_t pin[3], Double_t pout[3], Double_t costheta, Double_t sintheta, Double_t cosphi, Double_t sinphi);
		void 		DalitzDecay(TLorentzVector parent);
		void 		GenerateDecay();

		int mDebug = 0;

		TF1         *funMeson;
		TH1D        *histMeson;

		TF1 		*fRapidity;
		TF1 		*massfunMeson;
		TF1 		*massfunrho;
		TF1 		*massfunee;
		TH1D 		*hmassfunee;


		TH2D        *hCocktail;

		TH1D *pTRes1D[19];
		TF1  *pTResFun[19];
		double pTResPar[19][7] =  {
									{ 1.6637, 1.9347, 1.53494, 6.03621, 0.000156225, 0.00720455, 20732.2},
									{ 1.74984, 2.03021, 1.33671, 4.28385, -0.000362419, 0.00728465, 20634.8},
									{ 1.81318, 2.04927, 1.24935, 4.09537, -0.000337762, 0.00747221, 20052.4},
									{ 1.82843, 2.14228, 1.22715, 3.43323, -0.000276705, 0.0077572, 18898.4},
									{ 1.84743, 2.09733, 1.21416, 3.74061, -0.000273966, 0.00802061, 17987.8},
									{ 1.87571, 2.03729, 1.18184, 4.29472, -0.000241213, 0.00828199, 17232.7},
									{ 1.86892, 2.09715, 1.20676, 3.65561, -0.000252983, 0.00859205, 16551.7},
									{ 1.89405, 2.0887, 1.17424, 3.73803, -0.000247621, 0.00891028, 15987.6},
									{ 1.89175, 2.11199, 1.19518, 3.76984, -0.000279682, 0.00932575, 15363.4},
									{ 1.89054, 2.05286, 1.19011, 3.96226, -0.000364476, 0.0096548, 14829.9},
									{ 1.88233, 2.03773, 1.20658, 4.23582, -0.000353858, 0.0100075, 14333.1},
									{ 1.87335, 1.97932, 1.22586, 4.68102, -0.000425858, 0.0103718, 13800.8},
									{ 1.87471, 1.95117, 1.23945, 5.14921, -0.000455547, 0.0107908, 13298.1},
									{ 1.89675, 1.89297, 1.20876, 5.56997, -0.000532793, 0.0111968, 12805.9},
									{ 1.90299, 1.94859, 1.20021, 5.1045, -0.000556637, 0.0117113, 12298.1},
									{ 1.86911, 1.8589, 1.25773, 6.118, -0.000591135, 0.0120378, 11911.8},
									{ 1.86789, 1.84631, 1.28047, 6.36404, -0.00070459, 0.0125204, 11485.3},
									{ 1.87603, 1.87682, 1.26887, 6.17248, -0.000704703, 0.012979, 11086.7},
									{ 1.86888, 1.83297, 1.28275, 6.73365, -0.000792004, 0.0134098, 10725.3},
									};
		TH2D *PtRes2D;
		TF1         *funSmearPt; // momentum smearing
		TF1         *funSmearPtEmb; //momemtum resolution from embedding
		TF1         *momShape;

		TRandom3	*myRandom;

		TH1D        *hSampledPt;
		TH2D        *hEPSingleTrkEffvsPt;
		TH2D        *hEMSingleTrkEffvsPt;
		TH2D        *hRCPairRapidityvsParentRapidity;
		TH2D        *hMCPairPtvsParentPt;
		TH2D        *hMCAcc0PairPtvsParentPt;
		TH2D        *hMCAcc1PairPtvsParentPt;
		TH2D		*hMCMvsPt;
		TH2D		*hMCAcc0MvsPt;
		TH2D		*hMCAcc1MvsPt;
		TH2D		*hMCAcc2MvsPt;
		TH2D		*hRCAcc1MvsPt3D;
		TH2D		*hRCAcc2MvsPt3D;
		TH1D		*hMCAcc1PairRapidity;
		TH2D        *hMCAcc1PairEMPtvsEPPt;
		TH2D*       hPlusSmearvsPt;
		TH2D*       hMinusSmearvsPt;
		TH2D*       hMCAcc1MvsPtwoSmear;

    	TH1D*       RapiditySTARAcc;
    	TH1D*       RapiditywoSTARAcc;
		TH1D*		RapidityOutSTARAcc;
    	TH2D*       MeePtRapidityOver1;
    	TH1D*       MeeFullRapidity;
		TH1D*		MeeWopTCut;
		TH1D*		MeeWopTEtaCut;
		TH1D*		MeeWoEtaCut;
		TH2D*		pTWithLargeMee;
		TH1D*		MeeWpTEtaWoRapidity;
		TH1D*		MeeWpTEtaWRapidity;
		TH2D*		PtMVsPtPLargeMass;
		TH2D*		EtaMVsEtaPLargeMass;
		TH1D*		CutRecorder;
		

		static const Int_t mPtBins = 10;
		static const Int_t mYBins = 20;
		static const Int_t mPhiBins= 7;
		static const Int_t mCenBins = 9; //16; //9;
		const Double_t mPairPtCut[11]= {0.3,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0};
		const Double_t mCentCut[10]= {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,8.5,9.5};
		TAxis *PtAxis;
		TAxis *YAxis;
		TAxis *PhiAxis;
		TAxis *CentAxis;

		// for the polarization, not only Jpsi, pair itself polarization
		// TTree *tree;
		Float_t positron_theta_hx=-99.,positron_theta_cs=-99.,positron_phi_hx=-99.,positron_phi_cs=-99.;
		Float_t electron_theta_hx=-99.,electron_theta_cs=-99.,electron_phi_hx=-99.,electron_phi_cs=-99.;
		Float_t pair_pt,pair_eta,pair_phi,pair_InvM;
		Float_t lepton1_pt,lepton1_eta,lepton1_phi,lepton1_InvM;
		Float_t lepton2_pt,lepton2_eta,lepton2_phi,lepton2_InvM;
		TLorentzVector lepton1,lepton2;

		//Efficiency histograms for the polarization calculations
		TH1D* hMCAcc0Mass[mCenBins][mPtBins][mPhiBins];
		TH1D* hMCAcc1Mass[mCenBins][mPtBins][mPhiBins];
		TH1D* hRCAcc1Mass[mCenBins][mPtBins][mPhiBins];
		TH2D* hMCAcc0_CosthetapT;
		TH2D* hMCAcc1_CosthetapT;
		TH2D* hRCAcc1_CosthetapT;
		// TH3D* hMCAcc0PairCosThetaInvMPt;
		// TH3D* hMCAcc1PairCosThetaInvMPt;
		// TH3D* hRCAcc0PairCosThetaInvMPt;
		// TH3D* hMCAcc0PairCosThetaInvMPt_CS;
		// TH3D* hMCAcc1PairCosThetaInvMPt_CS;
		// TH3D* hRCAcc0PairCosThetaInvMPt_CS;
		TH2D* hMCAcc0PairCosThetaPt_HX;
		TH2D* hMCAcc1PairCosThetaPt_HX;
		TH2D* hRCAcc1PairCosThetaPt_HX;
		TH2D* hMCAcc0PairCosThetaPt_CS;
		TH2D* hMCAcc1PairCosThetaPt_CS;
		TH2D* hRCAcc1PairCosThetaPt_CS;
		TH2D* hMCAcc0PairPhiPt_HX;
		TH2D* hMCAcc1PairPhiPt_HX;
		TH2D* hRCAcc1PairPhiPt_HX;
		TH2D* hMCAcc0PairPhiPt_CS;
		TH2D* hMCAcc1PairPhiPt_CS;
		TH2D* hRCAcc1PairPhiPt_CS;

		Char_t  	MesonType[256];


	private:
		TLorentzVector Parent;
		TLorentzVector fProducts[3];
		Bool_t         mUseCocktailInput;
		Bool_t         mUseScaleEff;
		Int_t 		   mUseTsaPtSpectra;
		Int_t          mParIndex;
		Int_t          mDmode;
		Double_t       mMass;
		Double_t       mWidth;
		Int_t          mNTrks;
		Int_t          mCenIdx;
		Double_t       mMinRap;
		Double_t       mMaxRap;
		Double_t       mPtSmearPar[nSmearFac];
};

inline void VecMeson::SetUseCocktailInput(Bool_t kFlag) {mUseCocktailInput=kFlag;}
inline void VecMeson::SetUseScaleEff(Bool_t kFlag) {mUseScaleEff=kFlag;}
inline void VecMeson::SetInputPtSpectra(Int_t kFlag) {mUseTsaPtSpectra=kFlag;}
inline void VecMeson::SetNumberOfTracks(Int_t nTrks) {mNTrks=nTrks;}
inline void VecMeson::SetCentralityIdx(Int_t idx) {mCenIdx=idx;}
inline void VecMeson::SetRapidityRange(Double_t min, Double_t max) { mMinRap = min; mMaxRap = max;}
inline void VecMeson::SetPtSmearPar(Double_t *par){for(int i=0;i<nSmearFac;i++) mPtSmearPar[i]=par[i];}
#endif
