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
#include "TH3.h"
#include "TH3D.h"
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
									{ 1.66268, 1.99131, 1.46197, 4.49114, 9.71842e-05, 0.00715602, 87798.4, }, 
									{ 1.76017, 2.04854, 1.27276, 3.72734, -0.000686213, 0.00720795, 88120.8, }, 
									{ 1.80722, 2.09845, 1.21865, 3.42635, -0.000751521, 0.00737833, 86201.4, }, 
									{ 1.83816, 2.09986, 1.1918, 3.44015, -0.000747686, 0.00755741, 83504.5, }, 
									{ 1.85328, 2.05284, 1.17883, 3.63617, -0.000727032, 0.00775252, 78975.7, }, 
									{ 1.86194, 2.02749, 1.17761, 3.70037, -0.000687092, 0.00794673, 75846.9, }, 
									{ 1.87446, 2.02995, 1.16724, 3.73425, -0.000659317, 0.00809439, 74103.8, }, 
									{ 1.89458, 2.00627, 1.14958, 4.01988, -0.000680513, 0.00824834, 72515.1, }, 
									{ 1.8972, 2.01839, 1.14874, 4.0307, -0.000664356, 0.00838704, 71066.7, }, 
									{ 1.91504, 2.03513, 1.13439, 3.90167, -0.000664157, 0.00856234, 69665.2, }, 
									{ 1.91352, 2.04179, 1.13677, 3.94066, -0.000668979, 0.00874685, 68018.1, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
									{ 1.91236, 2.00246, 1.14803, 4.23183, -0.000672078, 0.00886363, 33536, }, 
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
		//3D histograms
		TH3D* hMCAcc0PairCosThetaInvMPt_HX;
		TH3D* hMCAcc1PairCosThetaInvMPt_HX;
		TH3D* hRCAcc1PairCosThetaInvMPt_HX;
		TH3D* hMCAcc0PairCosThetaInvMPt_CS;
		TH3D* hMCAcc1PairCosThetaInvMPt_CS;
		TH3D* hRCAcc1PairCosThetaInvMPt_CS;
		TH3D* hMCAcc0PairPhiInvMPt_HX;
		TH3D* hMCAcc1PairPhiInvMPt_HX;
		TH3D* hRCAcc1PairPhiInvMPt_HX;
		TH3D* hMCAcc0PairPhiInvMPt_CS;
		TH3D* hMCAcc1PairPhiInvMPt_CS;
		TH3D* hRCAcc1PairPhiInvMPt_CS;
		//2D histograms
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
