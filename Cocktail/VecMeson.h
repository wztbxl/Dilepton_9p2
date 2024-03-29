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

		Double_t	SampleMesonMass();
		Double_t 	EvalEff3D(TLorentzVector electron,int ipttpc,int ietatpc,int iphitpc,int ipttof,int ietatof,int iphitof,int charge);
		Double_t 	EvalScaleEff3D(TLorentzVector electron,int ipttpc,int ietatpc,int iphitpc,int ipttof,int ietatof,int iphitof,int charge);

		void 		Rot(Double_t pin[3], Double_t pout[3], Double_t costheta, Double_t sintheta, Double_t cosphi, Double_t sinphi);
		void 		DalitzDecay(TLorentzVector parent);
		void 		GenerateDecay();

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
		// double pTResPar[19][7] =  {
		// 							{ 1.6637, 1.9347, 1.53494, 6.03621, 0.000156225, 0.00720455, 20732.2},
		// 							{ 1.74984, 2.03021, 1.33671, 4.28385, -0.000362419, 0.00728465, 20634.8},
		// 							{ 1.81318, 2.04927, 1.24935, 4.09537, -0.000337762, 0.00747221, 20052.4},
		// 							{ 1.82843, 2.14228, 1.22715, 3.43323, -0.000276705, 0.0077572, 18898.4},
		// 							{ 1.84743, 2.09733, 1.21416, 3.74061, -0.000273966, 0.00802061, 17987.8},
		// 							{ 1.87571, 2.03729, 1.18184, 4.29472, -0.000241213, 0.00828199, 17232.7},
		// 							{ 1.86892, 2.09715, 1.20676, 3.65561, -0.000252983, 0.00859205, 16551.7},
		// 							{ 1.89405, 2.0887, 1.17424, 3.73803, -0.000247621, 0.00891028, 15987.6},
		// 							{ 1.89175, 2.11199, 1.19518, 3.76984, -0.000279682, 0.00932575, 15363.4},
		// 							{ 1.89054, 2.05286, 1.19011, 3.96226, -0.000364476, 0.0096548, 14829.9},
		// 							{ 1.88233, 2.03773, 1.20658, 4.23582, -0.000353858, 0.0100075, 14333.1},
		// 							{ 1.87335, 1.97932, 1.22586, 4.68102, -0.000425858, 0.0103718, 13800.8},
		// 							{ 1.87471, 1.95117, 1.23945, 5.14921, -0.000455547, 0.0107908, 13298.1},
		// 							{ 1.89675, 1.89297, 1.20876, 5.56997, -0.000532793, 0.0111968, 12805.9},
		// 							{ 1.90299, 1.94859, 1.20021, 5.1045, -0.000556637, 0.0117113, 12298.1},
		// 							{ 1.86911, 1.8589, 1.25773, 6.118, -0.000591135, 0.0120378, 11911.8},
		// 							{ 1.86789, 1.84631, 1.28047, 6.36404, -0.00070459, 0.0125204, 11485.3},
		// 							{ 1.87603, 1.87682, 1.26887, 6.17248, -0.000704703, 0.012979, 11086.7},
		// 							{ 1.86888, 1.83297, 1.28275, 6.73365, -0.000792004, 0.0134098, 10725.3},
		// 							};
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
