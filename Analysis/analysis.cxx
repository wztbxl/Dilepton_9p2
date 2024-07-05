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
#include "TTimer.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "miniDst.h"
#include "cuts.h"
#include "RefMfun.h"
#include "StRefMultCorr.h"
#include "pileup.h"

using namespace std;
#endif

void bookHistograms();
void writeHistograms(char* outFile);
void makeTags();
void makeRealPairs();
void makeMixPairs();
void copyCurrentToBuffer();
Bool_t Init();
Bool_t passEvent(miniDst* event);
Bool_t passTrack(miniDst* event, Int_t i);
Int_t  getCentralityBin9(Int_t cenBin16);
Double_t calReWeight(Double_t refMultCorr);
Double_t calCosTheta(TLorentzVector eVec,TLorentzVector eeVec);
Double_t reCalEventPlane(miniDst* event, Bool_t rejElectron = kFALSE);
Double_t phiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2);
bool nPi_K_P_rejection(int refmult, int nPi_K_P );
void Polarization(int icharge,int jcharge,TLorentzVector ivector,TLorentzVector jvector);
void fillHistograms(std::string  unlikeOrlike,TLorentzVector JPSI);
void fill3DHistograms(std::string unlikeOrlike,TLorentzVector JPSI,int i,int j,int pairs);
void GetPtPhiCentBin(TLorentzVector pair,TLorentzVector Positron, int _mCentrality,float eventphi,int &ptindex,int &yindex,int &phiindex,int &CentIndex,double &costhe,Bool_t tangent, Int_t Flag );

int mDebug = 0;
// int mDebug = 1;

int nPi_K_P_tof = 0;//used for pile rejection
TF1* f_upper = new TF1("f_upper","pol5",0,350);
TF1* f_lower = new TF1("f_lower","pol5",0,350);

TTimer   *timer;
TRandom3 *myRandom;

//variables 
Int_t dayIndex;
Int_t runIndex;
Short_t mCentrality;
map<Int_t,Int_t> mTotalDayId;
map<Int_t,Int_t> mTotalRunId;
map<Int_t,Int_t> mBadRunId_001;
map<Int_t,Int_t> mBadRunId_021;

Float_t bField;
Float_t reWeight;
Int_t iran = 0;

//for the polarization, not only Jpsi
// TTree *tree;
Float_t positron_theta_hx=-99.,positron_theta_cs=-99.,positron_phi_hx=-99.,positron_phi_cs=-99.;
Float_t electron_theta_hx=-99.,electron_theta_cs=-99.,electron_phi_hx=-99.,electron_phi_cs=-99.;
Float_t pair_pt,pair_eta,pair_phi,pair_InvM;
Float_t lepton1_pt,lepton1_eta,lepton1_phi,lepton1_InvM;
Float_t lepton2_pt,lepton2_eta,lepton2_phi,lepton2_InvM;
TLorentzVector lepton1,lepton2;

const Float_t pairPtCut = 0.1;

//StRefMultCorr* refMultCorrUtil;
//default categories for mixevent
//mCenBins=16; mVzBins=10; mEveBins=24; mMaxEventsInBuffer=100;
const Int_t mCenBins = 9; //16; //9;
const Int_t mVzBins = 10; //10; //6;
const Int_t mEveBins = 12; //24; //12;
const Int_t mPtBins = 10;
const Int_t mYBins = 20;
const Int_t mPhiBins= 7;
const Int_t mMaxEventsInBuffer = 350; //100; //50;
const Int_t mMaxElectrons = 70;
const Float_t mPhiVCutMRange = 0.2;
Float_t current_EQx[mMaxElectrons],current_EQy[mMaxElectrons];
vector <Float_t> vCurrent_eQx; 
vector <Float_t> vCurrent_eQy; 
vector <Float_t> vCurrent_e;
int current_nE = 0;
int current_nEPlus = 0;
int current_nEMinus = 0;
TLorentzVector current_ePlus[mMaxElectrons];
TLorentzVector current_eMinus[mMaxElectrons];
TLorentzVector current_e[mMaxElectrons];
int current_ePlus_CellID[mMaxElectrons];
int current_eMinus_CellID[mMaxElectrons];
int current_ePlus_tags[mMaxElectrons];
int current_eMinus_tags[mMaxElectrons];
Int_t cenBufferPointer, vzBufferPointer, eveBufferPointer;
Int_t nEventsInBuffer[mCenBins][mVzBins][mEveBins];
Bool_t bufferFullFlag[mCenBins][mVzBins][mEveBins];
Int_t buffer_nEPlus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer];
Int_t buffer_nEMinus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer];
TLorentzVector buffer_ePlus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];
TLorentzVector buffer_eMinus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];
int buffer_ePlus_CellID[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];
int buffer_eMinus_CellID[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];

//***** constrain the bad dedx calibration geometry *****
TF1 *funPosHi;
TF1 *funPosLow;
TF1 *funNegHi;
TF1 *funNegLow;
Float_t par[4][4];
Float_t parErr[4][4];

//********* define function and histograms *********
TF1 *phiVcut;
TF1 *Pileuplimit;
TF1 *PileupUplimit;
TF1 *PileupLowlimit;
TF1 *Delta_Psi2;

//in Init function
TProfile2D *ShiftFactorcos[mArrayLength];
TProfile2D *ShiftFactorsin[mArrayLength];
TProfile2D *etapluszplusQx;
TProfile2D *etapluszminusQx;
TProfile2D *etaminuszplusQx;
TProfile2D *etaminuszminusQx;
TProfile2D *etapluszplusQy;
TProfile2D *etapluszminusQy;
TProfile2D *etaminuszplusQy;
TProfile2D *etaminuszminusQy;
//in passEvent function
TH1D* hnEvts;
TH1D* hRunID;
TH1D* hTriggerID;
TH1F *hCentrality9;
TH1F *hRefMult;
TH1F *hVertexZ;
TH1F *hVzDiff;
TH1D *hVr;
TH1F *hBField;
TH2D *hnTofHitsvsRefMult;
TH2D* hnTofHitsvsRefMult_noCut;
TH2D* hnTofHitsvsRefMult_Vz35;
TH2D* hVxvsVy;
//eventPlane
TH1F *hRawEventPlane;
TH1F *hNewEventPlane;
TH1D *hNewEventPlaneEast;//rejection electron
TH1D *hNewEventPlaneWest;
TH1F *hReCenterEventPlane;
TH1D *hReCenterEventPlaneWest;
TH1D *hReCenterEventPlaneEast;
TH1F *hFinalEventPlane;
TH1D *hFinalEventPlaneWest;
TH1D *hFinalEventPlaneEast;
TH1F *hFinalEventPlane_Fit;
TH2D *hEventPlaneWestvsEast;
TH1D *hCosthetastar;
TH1D *hMixCosthetastar;
TH2F *hDelta_Psi2;
TH2F *hDelta_Psi2_FitvsFactor;
TH1D* hLargeDiffEvt_Day;
TH1D* hLargeDiffEvt_vz;
TH1D* hLargeDiffEvt_vr;
TProfile *EventPlanRes;

TH2F *hInclusiveEPhivsPt;
TH2F *hExclusiveEPhivsPt;
TH2F *hCut3EPhivsPt;
TH2F *hnEMinusvsEPlus;
//angleV 
//without phiV cut
// TH2F *hULMvsPtwophiV;
// TH2F *hLPosMvsPtwophiV;
// TH2F *hLNegMvsPtwophiV;
// TH2F *hMixULMvsPtwophiV;
// TH2F *hMixLPosMvsPtwophiV;
// TH2F *hMixLNegMvsPtwophiV;
//with phiV cut
TH2F *hULMvsPt;
TH2F *hLPosMvsPt;
TH2F *hLNegMvsPt;
TH2F *hMixULMvsPt;
TH2F *hMixLPosMvsPt;
TH2F *hMixLNegMvsPt;
//**********************
//add centrality dimension
TH3F *hULMvsPtCen;
TH3F *hLPosMvsPtCen;
TH3F *hLNegMvsPtCen;
TH3F *hMixULMvsPtCen;
TH3F *hMixLPosMvsPtCen;
TH3F *hMixLNegMvsPtCen;

TH1D* hCellIDDiff;
TH2F* hnSigmaEvsP;
//QA plot to dig the extra E issue;
TH2D* hdEdxvsP;
TH1D* hPt_Electron;
TH1D* hPt_Positron;

//histograms for the polarization
TH2F* hPairPhiPt;
TH2F* hPairPhiPtBG;
TH2F* hPairCosThetaPt;
TH2F* hPairCosThetaPtBG;
TH2F* hPairPhiPtHX;
TH2F* hPairPhiPtHXBG;
TH2F* hPairCosThetaPtCS;
TH2F* hPairCosThetaPtCSBG;
TH2F* hPairPhiPtCS;
TH2F* hPairPhiPtCSBG;
TH3F* hPairCosThetaInvMPt;
TH3F* hPairCosThetaInvMPtBG;
TH3F* hPairCosThetaInvMPtCS;
TH3F* hPairCosThetaInvMPtCSBG;
TH3F* hPairPhiInvMPt;
TH3F* hPairPhiInvMPtBG;
TH3F* hPairPhiInvMPtCS;
TH3F* hPairPhiInvMPtCSBG;
TH3F* hPairCosThetaPhiPt;
TH3F* hPairCosThetaPhiPtBG;
TH3F* hPairCosThetaPhiPtCS;
TH3F* hPairCosThetaPhiPtCSBG;

//M vs Pt
TH1F *hRapdity;
TH2F *hULMvsPtT;
TH2F *hULMvsPtTW;
TH1F *hULM[mCenBins][mPtBins][mPhiBins];
TH1F *hULYM[mCenBins][mYBins][mPhiBins];

//LS event
TH2F *hLSMvsPtT;
TH2F *hLSMvsPtTW;
TH1F *hLSPlusM[mCenBins][mPtBins][mPhiBins];
TH1F *hLSPlusYM[mCenBins][mYBins][mPhiBins];
TH1F *hLSMinusM[mCenBins][mPtBins][mPhiBins];
TH1F *hLSMinusYM[mCenBins][mYBins][mPhiBins];

//mix event
TH2F *hMixULMvsPtT;
TH2F *hMixULMvsPtTW;
TH1F *hMixULM[mCenBins][mPtBins][mPhiBins];
TH1F *hMixLSPosM[mCenBins][mPtBins][mPhiBins];
TH1F *hMixLSNegM[mCenBins][mPtBins][mPhiBins];
TH1F *hMixULYM[mCenBins][mYBins][mPhiBins];
TH1F *hMixLSPosYM[mCenBins][mYBins][mPhiBins];
TH1F *hMixLSMinusYM[mCenBins][mYBins][mPhiBins];


TAxis *PtAxis;
TAxis *YAxis;
TAxis *PhiAxis;
TAxis *CentAxis;
const Double_t mPairPtCut[11]= {0.3,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0};
const Double_t mCentCut[10]= {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,8.5,9.5};
Int_t runId;
Int_t EvtID;
Double_t finalEventPlane;

int main(int argc, char** argv)
{
	if(argc!=1&&argc!=3) return -1;

	TString inFile="test.list";
	char outFile[1024];
	sprintf(outFile,"test/test");
	if(argc==3){
		inFile = argv[1];
		sprintf(outFile,"%s",argv[2]);
	}

	//+---------------------------------+
	//| open files and add to the chain |
	//+---------------------------------+
	TChain *chain = new TChain("miniDst");

	Int_t ifile=0;
	char filename[512];
	ifstream *inputStream = new ifstream;
	inputStream->open(inFile.Data());
	if (!(inputStream)) {
		printf("can not open list file\n");
		return 0;
	}
	for(;inputStream->good();){
		inputStream->getline(filename,512);
		if(inputStream->good()) {
			TFile *ftmp = new TFile(filename);
			if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
				cout<<"something wrong"<<endl;
			} else {
				cout << filename << endl;
				cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
				chain->Add(filename);
				ifile++;
				if(mDebug) ftmp->Print();
			}
			delete ftmp;
		}
	}
	delete inputStream;

	//intialization
	bookHistograms();

	if( !Init() ){
		cout<<"The initialization is failed !!!"<<endl;
		return 0;
	}else{
		timer = new TTimer();
		myRandom = new TRandom3();

		//+-------------------+
		//| initialize buffer |
		//+-------------------+
		if(mDebug) cout<<"before memset events"<<endl;
		memset(nEventsInBuffer,0,sizeof(nEventsInBuffer));
		memset(bufferFullFlag,0,sizeof(bufferFullFlag));
		memset(buffer_nEPlus,0,sizeof(buffer_nEPlus));
		memset(buffer_nEMinus,0,sizeof(buffer_nEMinus));
		if(mDebug) cout<<"after memset events"<<endl;
	}


	//+-------------+
	//| loop events |
	//+-------------+
	if(mDebug) chain->Print();
	miniDst *event = new miniDst(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;
	//refMultCorrUtil = new StRefMultCorr("refmult");
	for(int i=0;i<nEvts;i++){

		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;
		if(mDebug) cout << "begin " << i << "th entry...." << endl;
		// chain->GetEntry(i);
		event->GetEntry(i);
		 runId = event->mRunId;

		map<Int_t,Int_t>::iterator iter = mTotalDayId.find((runId/1000)%1000);
		if(iter != mTotalDayId.end())
			dayIndex = iter->second;
		else{
			dayIndex = -1;
			cout<<"Can not find the dayId in the mTotalDayId list"<<endl;
			continue;
		}
    // cout << "begin " << i << "th entry...." << endl;
		if(mDebug) cout << "Day index =  " << dayIndex << endl;

		if(dayIndex<0) continue;

		iter = mTotalRunId.find(runId);
		if(iter != mTotalRunId.end()){
			runIndex = iter->second;
		}
		else{
			cout<<"runNumber:"<<runId<<endl;
			cout<<"Can not find the runNumber in the mTotalRunId list"<<endl;
			continue;

		}

		if(i%1000==0){
			long long tmp = (long long)timer->GetAbsTime();
			UInt_t seed = tmp/myRandom->Rndm();
			myRandom->SetSeed(seed);
			//cout<<"random number:"<<myRandom->Uniform(-1,1)<<endl;
		}

		if(mDebug) cout << "Before pass Event = " << endl;
		if(!passEvent(event)) continue; 
		EvtID = event->mEventId;
		if(mDebug) cout << "EvtID = " << EvtID << endl;
		current_nE=0;
		current_nEPlus=0;
		current_nEMinus=0;
		Int_t npTrks = event->mNTrks;
		nPi_K_P_tof = event->mnChargeParticle;
		if(mDebug) cout << "npTrks = " << npTrks << endl;
		// nPi_K_P_tof = 0;
		for(int j=0;j<npTrks;j++) passTrack(event,j); //Trk loop
		if(mDebug) cout << "after passtrack" << endl;
		hnEMinusvsEPlus->Fill(current_nEPlus,current_nEMinus);
		// cout << "nPi_K_P_tof = " << nPi_K_P_tof << endl; 

		// if(!nPi_K_P_rejection(event->mRefMult,nPi_K_P_tof))
		// {
		// 	continue;
		// }
		

		finalEventPlane = reCalEventPlane(event, kTRUE);//reject the electron contribution
		if(mDebug) cout << "after recal Event Plane" << endl;
		if(finalEventPlane<0) continue;
		eveBufferPointer = (Int_t)(finalEventPlane/TMath::Pi()*mEveBins);
		// // cout<<"eveBufferPointer:"<<eveBufferPointer<<endl;
		// if(eveBufferPointer<0 || eveBufferPointer>=mEveBins) continue;
		// eveBufferPointer = 0;


		for (int i = 0; i < mMaxElectrons; i++)
		{
			current_ePlus_tags[i] = 1;
			current_eMinus_tags[i] = 1;
		}
		if(mDebug) cout << "beforetags " << endl;
		makeTags();
		if(mDebug) cout << "after tags " << endl;
		makeRealPairs();
		if(mDebug) cout << "after real pairs " << endl;
		makeMixPairs();
		if(mDebug) cout << "after mixed pairs " << endl;
		copyCurrentToBuffer();
		if(mDebug) cout << "after copy ro buffer " << endl;
	}

  cout << "start checking buffer full flag " << endl;
  for (int iCent = 0; iCent < mCenBins; iCent++)
  {
    for (int iVz = 0; iVz < mVzBins; iVz++)
    {
      for (int iEve = 0; iEve < mEveBins; iEve++)
      {
         cout << "bufferFullFlag = " << bufferFullFlag[iCent][iVz][iEve] << endl;
      }
    }

  } 

	writeHistograms(outFile);
	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
bool nPi_K_P_rejection(int refmult, int nPi_K_P )
{
	if ( nPi_K_P >= f_lower->Eval(refmult) && nPi_K_P < f_upper->Eval(refmult))
	{
		return kTRUE; // pass the cut
	} else return kFALSE; // did not pass the cut
}

//________________________________________________________________
Bool_t passEvent(miniDst* event)
{
	Int_t runId  = event->mRunId;

	Float_t vx = event->mVertexX;
	Float_t vy = event->mVertexY;
	Float_t vz = event->mVertexZ;
	Float_t vr = sqrt(vx*vx+vy*vy);
	Float_t vpdVz = event->mVpdVz;
	Float_t ZDCrate = event->mZDCRate;

	Float_t vzDiff = vz - vpdVz;
	Int_t mnTOFMatch = event->mnTOFMatch;
	Int_t refMult = event->mRefMult;
	if(mDebug) cout << "refMult = " << refMult <<endl;
	int  nTrigs = event->mNTrigs;
	bool fireTrigger = kFALSE;
	bool RefMVzCorFlag = kFALSE;
	Bool_t is001Trigger = kFALSE;
	Bool_t is021Trigger = kFALSE;
	for(int i=0; i< nTrigs; i++){
		int trigId = event->mTrigId[i];
		if(trigId == mTrigId[0] || trigId == mTrigId[2]){
			fireTrigger = kTRUE;
		}
		if(trigId == mTrigId[0])RefMVzCorFlag = kTRUE, is001Trigger = kTRUE;
		if(trigId == mTrigId[2])is021Trigger = kTRUE;
		hTriggerID->Fill(trigId);
	}
	if(!fireTrigger) return kFALSE;
	bField = event->mBField;
	// mCentrality = event->mCentrality;

	map<Int_t, Int_t>::iterator iter_001 = mBadRunId_001.find(runId);
	if(iter_001 != mBadRunId_001.end() && is001Trigger){
		if(mDebug) cout<<"bad run, continue"<<endl;
		return kFALSE;
	}

	map<Int_t, Int_t>::iterator iter_021 = mBadRunId_001.find(runId);
	if(iter_021 != mBadRunId_001.end() && is021Trigger){ // using same bad runlist for the test
		if(mDebug) cout<<"bad run, continue"<<endl;
		return kFALSE;
	}

	hnEvts->Fill(1);
	hRunID->Fill(runId);
	// reWeight = 1.;
	Double_t RefMultCorr = refMult;
	mCentrality = (int)event->mCentrality;
	if(mDebug) cout << "mCentrality = " << mCentrality <<endl;
  // mCentrality = mCentrality+1;
  cenBufferPointer = mCentrality;
  RefMultCorr = event->mGRefMultCorr;
  reWeight = event->mEvtWeight;

  // if(RefMVzCorFlag)RefMultCorr = GetRefMultCorr(refMult, vz);
	// reWeight = GetWeight(RefMultCorr);
	// mCentrality = GetCentrality(RefMultCorr);
	cenBufferPointer = mCentrality-1;
	if (cenBufferPointer <0 || cenBufferPointer >8) return kFALSE;// 0-8 for 70-80% - 0-5%
	
	
	// refMultCorrUtil->init(runId);
	// refMultCorrUtil->initEvent(refMult,vz,ZDCrate);
	// mCentrality = refMultCorrUtil->getCentralityBin9();
	// reWeight = refMultCorrUtil->getWeight();
	// double refMultCorr = refMultCorrUtil->getRefMultCorr();
	// cenBufferPointer = mCentrality;
	// if (cenBufferPointer <0 || cenBufferPointer >8) return kFALSE;

	//cout << cenBufferPointer<<endl;
	//if(refMult<300) cout<<"reWeight: "<<reWeight<<endl;

	hnTofHitsvsRefMult_noCut->Fill(refMult,mnTOFMatch);
	if(TMath::Abs(vx)<1.e-5 && TMath::Abs(vy)<1.e-5 && TMath::Abs(vz)<1.e-5) return kFALSE;
	if(!pileupRejection(vz, refMult, mnTOFMatch)) return kFALSE;
	hnEvts->Fill(5);
	if(TMath::Abs(vz)>=mVzCut) return kFALSE;//vz should also be in the range listed in the parameters file to do the refMult correction
	hnEvts->Fill(2);
	if(vr>=mVrCut) return kFALSE;
	hnEvts->Fill(3);
	// hnTofHitsvsRefMult_noCut->Fill(refMult,mnTOFMatch);
	// if(TMath::Abs(vzDiff)>=mVzDiffCut) return kFALSE;
	//pile up rejection
	hnTofHitsvsRefMult_Vz35->Fill(refMult,mnTOFMatch);
	hnEvts->Fill(4);
	// if (mnTOFMatch < Pileuplimit->Eval(refMult)) return kFALSE;
	if (mDebug) cout << "after pileup" << endl;

	hBField->Fill(bField);
	hVertexZ->Fill(vz);
	hVzDiff->Fill(vzDiff);
	hVr->Fill(vr);
	hVxvsVy->Fill(vx,vy);

	hRefMult->Fill(refMult,reWeight);
	hnTofHitsvsRefMult->Fill(refMult,mnTOFMatch);
	

	Int_t centrality9 = mCentrality;
	hCentrality9->Fill(centrality9,reWeight);

	vzBufferPointer = (Int_t)((vz+mVzCut)/(2*mVzCut)*mVzBins);
	if(vzBufferPointer<0 || vzBufferPointer>=mVzBins) return kFALSE;
	if (mDebug) cout << "after Vz" << endl;


	return kTRUE;
}
//______________________________________________________________
Bool_t passTrack(miniDst* event, Int_t i)
{
	Int_t charge = event->mCharge[i];
	Int_t nHitsFit = event->mNHitsFit[i];
	Int_t nHitsDedx = event->mNHitsDedx[i];
	Int_t nHitsPoss = event->mNHitsPoss[i];
	Float_t nSigmaE = event->mNSigmaE[i];
	Float_t dca = event->mDca[i];
	Float_t pt = event->mPt[i];
	Float_t eta = event->mEta[i];
	Float_t phi = event->mPhi[i];
	Float_t beta2TOF = event->mBeta2TOF[i];
	Float_t TOFLoaclY = event->mTOFLocalY[i];
	Float_t ratio = 1.0*nHitsFit/nHitsPoss;
	int CellID = event->mTOFCellID[i];
	TVector3 mom;
	mom.SetPtEtaPhi(pt,eta,phi);
	Float_t p = mom.Mag();
	double msquare =  -999;
	msquare = pow(p, 2) * (1 - pow(beta2TOF, 2)) / pow(beta2TOF, 2);

	// if(TMath::Abs(msquare-0.879)<0.020 || TMath::Abs(msquare-0.243)<0.005 || TMath::Abs(msquare-0.019)<0.003) nPi_K_P_tof = nPi_K_P_tof+1;

//   if(charge < 0) return kFALSE;
	// if(charge<0)cout<<"charge="<<charge<<endl;
	if(pt<mTpcePtCut[0] || pt>mTpcePtCut[1]) return kFALSE;
	// if(nHitsFit<15) return kFALSE;
	if(nHitsFit<mTpceNHitsFitCut) return kFALSE;
	if(ratio<mTpceNHitsFitRatioCut) return kFALSE;
	// if(nHitsDedx<20) return kFALSE;
	if(nHitsDedx<mTpceNHitsDedxCut) return kFALSE;
	// if(dca>0.8) return kFALSE;
	if(dca>mTpceDcaCut) return kFALSE;
	if(TMath::Abs(eta)>mTpceEtaCut) return kFALSE;
	hInclusiveEPhivsPt->Fill(charge*pt,phi);
	if(beta2TOF<=0. || TMath::Abs(1.-1./beta2TOF)>mTpceBeta2TOFCut) return kFALSE;
	if(abs(TOFLoaclY) > 1.8) return kFALSE;
	hnSigmaEvsP->Fill(p,nSigmaE);

	hExclusiveEPhivsPt->Fill(charge*pt,phi);
	Float_t mTpceNSigmaECutLow;
	if(p<.8){
		mTpceNSigmaECutLow = 3.0*p - 3.15; 
	}else{
		mTpceNSigmaECutLow = mTpceNSigmaECut[0];
	}
	if(nSigmaE<mTpceNSigmaECutLow+mNSigmaEShift || nSigmaE>mTpceNSigmaECut[1]+mNSigmaEShift) return kFALSE;
	hCut3EPhivsPt->Fill(charge*pt,phi);

	// if(current_nE != current_nEPlus+current_nEMinus) current_nE = current_nEPlus+current_nEMinus;
	// cout << "current_nE = " << current_nE << " current_nEPlus = " << current_nEPlus << " current_nEMinus = " << current_nEMinus << endl;
	vCurrent_eQx.push_back(pt*TMath::Cos(2*phi));
	vCurrent_eQy.push_back(pt*TMath::Sin(2*phi));
	// cout << "current_nE = " << current_nE << " current_nEPlus = " << current_nEPlus << " current_nEMinus = " << current_nEMinus << endl;
	current_EQx[current_nE] = pt*TMath::Cos(2*phi);
	current_EQy[current_nE] = pt*TMath::Sin(2*phi);
	current_e[current_nE].SetPtEtaPhiM(pt,eta,phi,Melectron);
	// current_nE = vCurrent_eQx.size();
	// cout << "after vector "<< endl;
	current_nE++;


	if(charge==1){
		current_ePlus[current_nEPlus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_ePlus_CellID[current_nEPlus] = CellID;
		hPt_Positron->Fill(pt);
		current_nEPlus++;
		// cout << "eP "<< endl;

	}
	else if(charge==-1){
		current_eMinus[current_nEMinus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_eMinus_CellID[current_nEMinus] = CellID;
		hPt_Electron->Fill(pt);
		current_nEMinus++;
		// cout << "eM "<< endl;
	}

	return kTRUE;
}

void makeTags()
{
	// cout << " nElectron = " << current_nEMinus << endl;
	// cout << " nPositron = " << current_nEPlus << endl;

	TLorentzVector pair(0,0,0,0);
	// e+e- and turn the electron under cuts tag to 0
	for (int i = 0; i < current_nEPlus; i++)
	{
		for ( int j = 0; j < current_nEMinus; j++)
		{
			pair = current_ePlus[i]+current_eMinus[j];
			// cout << "pair PseudoRapidity = " << pair.PseudoRapidity() << endl;
			if(TMath::Abs(pair.Rapidity())<=mPairYCut)
			{
				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_ePlus[i],current_eMinus[j],1,-1);
				// if( angleV < angleVcut && pair.M()<mPhiVCutMRange )
				// {
				// 	current_ePlus_tags[i] = 0;
				// 	current_eMinus_tags[j] = 0;
				// }
				// if( pair.M() < 0.055 )
				// {
				// 	current_ePlus_tags[i] = 0;
				// 	current_eMinus_tags[j] = 0;
				// }
			}
			
		}
	}
	// for (int i = 0; i < current_nEPlus; i++)
	// {
	// 	for (int j = i+1; j < current_nEPlus; j++)
	// 	{
	// 		pair = current_ePlus[i]+current_ePlus[j];
	// 		// cout << "pair PseudoRapidity = " << pair.PseudoRapidity() << endl;
	// 		if(TMath::Abs(pair.Rapidity())<=mPairYCut)
	// 		{

	// 			Double_t angleVcut = phiVcut->Eval(pair.M());
	// 			Double_t angleV = phiVAngle(current_ePlus[i],current_ePlus[j],1,1);
	// 			// if( angleV < angleVcut && pair.M()<mPhiVCutMRange )
	// 			// {
	// 			// 	current_ePlus_tags[i] = 1;
	// 			// 	current_ePlus_tags[j] = 1;
	// 			// }
	// 			if( pair.M() < mMassCut )
	// 			{
	// 				current_ePlus_tags[i] = 0;
	// 				current_ePlus_tags[j] = 0;
	// 			}
	// 		}
	// 	}
	// }
	// for (int i = 0; i < current_nEMinus; i++)
	// {	
	// 	for (int j = i+1; j < current_nEMinus; j++)
	// 	{
	// 		pair = current_eMinus[i]+current_eMinus[j];
	// 		// cout << "pair PseudoRapidity = " << pair.PseudoRapidity() << endl;
	// 		if(TMath::Abs(pair.Rapidity())<=mPairYCut)
	// 		{
	// 			Double_t angleVcut = phiVcut->Eval(pair.M());
	// 			Double_t angleV = phiVAngle(current_eMinus[i],current_eMinus[j],-1,-1);
	// 			// if( angleV < angleVcut && pair.M()<mPhiVCutMRange )
	// 			// {
	// 			// 	current_eMinus_tags[i] = 1;
	// 			// 	current_eMinus_tags[j] = 1;
	// 			// }
	// 			if( pair.M() < mMassCut )
	// 			{
	// 				current_eMinus_tags[i] = 0;
	// 				current_eMinus_tags[j] = 0;
	// 			}
	// 		}
			
	// 	}
	// }


}

void makeRealPairs()
{
	Int_t _PtIndex = -999;
	Int_t _YIndex = -999;
	Int_t _PhiIndex = -999;
	Int_t _CentIndex = -999;
	double costhetastar =-999.;
	//+--------------------------+
	//| current e+  + current e- |
	//+--------------------------+
	TLorentzVector pair(0,0,0,0);
	for(Int_t i=0;i<current_nEPlus;i++){
		if ( current_ePlus_tags[i] == 0) continue;
		for(Int_t j=0;j<current_nEMinus;j++){ 
			if ( current_eMinus_tags[j] == 0 ) continue;
			pair = current_ePlus[i]+current_eMinus[j];
			if(TMath::Abs(pair.Rapidity())<=mPairYCut){
				// hULMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
				//hULMvsPtwophiV->Fill(pair.Pt(),pair.M());

				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_ePlus[i],current_eMinus[j],1,-1);
				// if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hULMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				// hULMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				
				if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
					hULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
					//hULMvsPt->Fill(pair.Pt(),pair.M());
					hULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);

					if(mDebug) cout << "before polarization calcultion" << endl;
					//this for the pair polarization it self
					if (pair.M()>0.2 && pair.M() < 1.1 )
					{
						if(mDebug) cout << "before polarization calcultion" << endl;
						Polarization(1,-1,current_ePlus[i],current_eMinus[j]);
						if(mDebug) cout << "before Fill3D calcultion" << endl;
						fill3DHistograms("unlike",pair,i,j,0);
						if(mDebug) cout << "before Fillhis calcultion" << endl;
						fillHistograms("unlike",pair);
					}
					if(mDebug) cout << "after polarization calcultion" << endl;

					GetPtPhiCentBin(pair, current_ePlus[i], mCentrality, finalEventPlane, _PtIndex, _YIndex, _PhiIndex, _CentIndex,costhetastar, 0, 1);// 0 for how to calculate costheta* 1 for nonsense 
					hCosthetastar->Fill(costhetastar);
					if(_PtIndex > mPtBins-1 || _PtIndex<0 || _PhiIndex > mPhiBins-1 || _PhiIndex <0)continue; 

					hULM[_CentIndex][_PtIndex][_PhiIndex]->Fill(pair.M(),reWeight);
					hULYM[_CentIndex][_YIndex][_PhiIndex]->Fill(pair.M(),reWeight);

					if(pair.Pt()<pairPtCut){
						Double_t costheta = calCosTheta(current_ePlus[i], pair);
						// hULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						//hULMvsPt->Fill(pair.Pt(),pair.M());
						// hULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						//hULCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

						// hULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
						// hULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[j].Pt(), reWeight);
					}
				}
			}
		}//end of e- loop
	}//end of e+ loop

	//+--------------------------+
	//| current e+  + current e+ |
	//+--------------------------+
	for(Int_t i=0;i<current_nEPlus;i++){
		if ( current_ePlus_tags[i] == 0 ) continue;
		for(Int_t j=i+1;j<current_nEPlus;j++){
			if ( current_ePlus_tags[j] == 0 ) continue;
			pair = current_ePlus[i]+current_ePlus[j];
			if(TMath::Abs(pair.Rapidity())<=mPairYCut){
				// hLPosMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
				//hLPosMvsPtwophiV->Fill(pair.Pt(),pair.M());

				// double TOF1x = cos(current_ePlus[i].Phi())*220;
				// double TOF1y = sin(current_ePlus[i].Phi())*220;
				// double TOF2x = cos(current_ePlus[j].Phi())*220;
				// double TOF2y = sin(current_ePlus[j].Phi())*220;
				// if ( sqrt( (TOF1x-TOF2x)*(TOF1x-TOF2x) + (TOF1y-TOF2y)*(TOF1y-TOF2y) ) < 6 ) continue;

				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_ePlus[i],current_ePlus[j],1,1);
				// if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hLPosMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				// hLPosMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
				// if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
					hLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
					//hLPosMvsPt->Fill(pair.Pt(),pair.M());
					hLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					if (pair.M()>0.2 && pair.M() < 1.1 )
					{
						Polarization(1,1,current_ePlus[i],current_ePlus[j]);
						fill3DHistograms("like",pair,i,j,0);
						fillHistograms("like",pair);
					}
					GetPtPhiCentBin(pair,current_ePlus[i], mCentrality, finalEventPlane, _PtIndex, _YIndex, _PhiIndex, _CentIndex, costhetastar, 0, 1);
					hLSPlusM[_CentIndex][_PtIndex][_PhiIndex]->Fill(pair.M(), reWeight);
					hLSPlusYM[_CentIndex][_YIndex][_PhiIndex]->Fill(pair.M(), reWeight);
					// hLPosMvsPhiCen->Fill(pair.Phi(),cenBufferPointer,pair.M(),reWeight);

					if(pair.Pt()<pairPtCut){
						Double_t costheta = calCosTheta(current_ePlus[i], pair);
						//hLPosCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[j].Pt(), reWeight);
						// hLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						//hLPosMvsPt->Fill(pair.Pt(),pair.M());
						// hLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					}
				}
			}
		}//end of e+ loop
	}//end of e+ loop

	//+--------------------------+
	//| current e-  + current e- |
	//+--------------------------+
	for(Int_t i=0;i<current_nEMinus;i++){
		if ( current_eMinus_tags[i] == 0 ) continue;
		for(Int_t j=i+1;j<current_nEMinus;j++){
			if ( current_eMinus_tags[j] == 0 ) continue;
			pair = current_eMinus[i]+current_eMinus[j];
			if(TMath::Abs(pair.Rapidity())<=mPairYCut){
				// hLNegMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
				//hLNegMvsPtwophiV->Fill(pair.Pt(),pair.M());

				// double TOF1x = cos(current_eMinus[i].Phi())*220;
				// double TOF1y = sin(current_eMinus[i].Phi())*220;
				// double TOF2x = cos(current_eMinus[j].Phi())*220;
				// double TOF2y = sin(current_eMinus[j].Phi())*220;
				// if ( sqrt( (TOF1x-TOF2x)*(TOF1x-TOF2x) + (TOF1y-TOF2y)*(TOF1y-TOF2y) ) < 6 ) continue;

				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_eMinus[i],current_eMinus[j],-1,-1);
				// if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hLNegMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				// hLNegMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
					hLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
					//hLNegMvsPt->Fill(pair.Pt(),pair.M());
					hLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					// hLNegMvsPhiCen->Fill(pair.Phi(),cenBufferPointer,pair.M(),reWeight);

					if (pair.M()>0.2 && pair.M() < 1.1 )
					{
						Polarization(-1,-1,current_eMinus[i],current_eMinus[j]);
						fill3DHistograms("like",pair,i,j,0);
						fillHistograms("like",pair);
					}
					GetPtPhiCentBin(pair,current_eMinus[i], mCentrality, finalEventPlane, _PtIndex, _YIndex, _PhiIndex, _CentIndex, costhetastar, 0, 1);
					hLSMinusM[_CentIndex][_PtIndex][_PhiIndex]->Fill(pair.M(), reWeight);
					hLSMinusYM[_CentIndex][_YIndex][_PhiIndex]->Fill(pair.M(), reWeight);

					if(pair.Pt()<pairPtCut){
						Double_t costheta = calCosTheta(current_eMinus[i], pair);
						//hLNegCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[i].Pt(), reWeight);
						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[j].Pt(), reWeight);
						// hLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						//hLNegMvsPt->Fill(pair.Pt(),pair.M());
						// hLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					}
				}
			}
		}//end of e- loop
	}//end of e- loop
}
//_____________________________________________________________________________
void makeMixPairs()
{
	Int_t _PtIndex = -999;
	Int_t _YIndex = -999;
	Int_t _PhiIndex = -999;
	Int_t _CentIndex = -999;
	double costhetastar =-999.;
	TLorentzVector pair(0,0,0,0);
	for(Int_t iBufferEvent=0;iBufferEvent<nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer];iBufferEvent++){
		//+-------------------------+
		//| current e+  + buffer e- |
		//+-------------------------+
		for(Int_t i=0;i<current_nEPlus;i++){
			for(Int_t j=0;j<buffer_nEMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_ePlus_tags[i] == 0 ) continue;
				pair = current_ePlus[i] + buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
				if ( abs( current_ePlus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) == 0 ) continue;
				// if ( abs( current_ePlus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;
				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixULMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_ePlus[i],buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],1,-1);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);

						// if (pair.M()>0.2 && pair.M() < 1.1 )
						// {
							GetPtPhiCentBin(pair,current_ePlus[i], mCentrality, finalEventPlane, _PtIndex, _YIndex, _PhiIndex, _CentIndex, costhetastar, 0, 1);
							hMixULM[_CentIndex][_PtIndex][_PhiIndex]->Fill(pair.M(), reWeight);
							hMixULYM[_CentIndex][_YIndex][_PhiIndex]->Fill(pair.M(), reWeight);
						// }
						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(current_ePlus[i], pair);
							//hMixULCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}
					}
				}
			}//end of buffer e- loop
		}//end of current e+ loop

		//+-------------------------+
		//| current e-  + buffer e+ |
		//+-------------------------+
		for(Int_t i=0;i<current_nEMinus;i++){
			for(Int_t j=0;j<buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_eMinus_tags[i] == 0) continue;
				pair = current_eMinus[i] + buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
				if ( abs( current_eMinus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) == 0 ) continue;
				// if ( abs( current_eMinus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;

				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixULMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_eMinus[i],buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],-1,1);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						// if (pair.M()>0.2 && pair.M() < 1.1 )
						// {
							GetPtPhiCentBin(pair,buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j], mCentrality, finalEventPlane, _PtIndex, _YIndex, _PhiIndex, _CentIndex, costhetastar, 0, 1);
							hMixULM[_CentIndex][_PtIndex][_PhiIndex]->Fill(pair.M(), reWeight);
							hMixULYM[_CentIndex][_YIndex][_PhiIndex]->Fill(pair.M(), reWeight);
						// }

						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j], pair);
							//hMixULCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[i].Pt(), reWeight);
							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}
					}
				}
			}//end of buffer e+ loop
		}//end of current e- loop

		//+-------------------------+
		//| current e+  + buffer e+ |
		//+-------------------------+
		for(Int_t i=0;i<current_nEPlus;i++){
			for(Int_t j=0;j<buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_ePlus_tags[i] == 0 )  continue;
				pair = current_ePlus[i] + buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];

				if ( abs( current_ePlus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) == 0 ) continue;
				// if ( abs( current_ePlus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;
				hCellIDDiff->Fill(  current_ePlus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] );

				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixLPosMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_ePlus[i],buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],1,1);
					// if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hMixLPosMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);

					// if (pair.M()>0.2 && pair.M() < 1.1 )
					// {
						GetPtPhiCentBin(pair,current_ePlus[i], mCentrality, finalEventPlane, _PtIndex, _YIndex, _PhiIndex, _CentIndex, costhetastar, 0, 1);
						hMixLSPosM[_CentIndex][_PtIndex][_PhiIndex]->Fill(pair.M(), reWeight);
						hMixLSPosYM[_CentIndex][_YIndex][_PhiIndex]->Fill(pair.M(), reWeight);
					// }
						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(current_ePlus[i], pair);
							// hMixLPosCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}

					}
				}
			}//end of buffer e+ loop
		}//endl of current e+ loop

		//+-------------------------+
		//| current e-  + buffer e- |
		//+-------------------------+
		for(Int_t i=0;i<current_nEMinus;i++){
			for(Int_t j=0;j<buffer_nEMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_eMinus_tags[i] == 0 ) continue;
				pair = current_eMinus[i] + buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];

				if ( abs( current_eMinus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;
				hCellIDDiff->Fill(  current_eMinus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] );

				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixLNegMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_eMinus[i],buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],-1,-1);
					// if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hMixLNegMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					// if (pair.M()>0.2 && pair.M() < 1.1 )
					// {
						GetPtPhiCentBin(pair,current_eMinus[i], mCentrality, finalEventPlane, _PtIndex, _YIndex, _PhiIndex, _CentIndex, costhetastar, 0, 1);
						hMixLSNegM[_CentIndex][_PtIndex][_PhiIndex]->Fill(pair.M(), reWeight);
						hMixLSMinusYM[_CentIndex][_YIndex][_PhiIndex]->Fill(pair.M(), reWeight);
					// }

						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(current_eMinus[i], pair);
							//hMixLNegCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[i].Pt(), reWeight);
							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}
					}
				}
			}//end of buffer e- loop
		}//endl of current e- loop
	}
}
//_____________________________________________________________________________
void copyCurrentToBuffer()
{
	if(nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer]>=mMaxEventsInBuffer) bufferFullFlag[cenBufferPointer][vzBufferPointer][eveBufferPointer] = kTRUE;
	Int_t eventPointer = -1;
	if(bufferFullFlag[cenBufferPointer][vzBufferPointer][eveBufferPointer]){
		eventPointer = (Int_t)myRandom->Uniform(0,mMaxEventsInBuffer-1.e-6);
	}else{
		eventPointer = nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer];
	}

	buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nEPlus;
	buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nEPlus;
	int nTrks = 0;
	for(Int_t i=0;i<current_nEPlus; i++){
		// if (current_ePlus_tags[i] == 0 ) continue;
		buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_ePlus[i];
		buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_ePlus_CellID[i];
		nTrks++;
	}

	buffer_nEMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nEMinus;
	nTrks = 0;
	for(Int_t i=0;i<current_nEMinus;i++){
		// if (current_eMinus_tags[i] == 0 ) continue;
		buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_eMinus[i];
		buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_eMinus_CellID[i];
		nTrks++;
	}

	if(nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer]<mMaxEventsInBuffer){
		nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer]++;
	}
}
//_____________________________________________________________________________
Int_t getCentralityBin9(Int_t cenBin16)
{
	if(cenBin16<0 || cenBin16>15) return -1;
	else{
		if(cenBin16==15) return 8;
		else if(cenBin16==14) return 7;
		else return (Int_t)(0.5*cenBin16);
	}
}
//_____________________________________________________________________________
Double_t reCalEventPlane(miniDst* event, Bool_t rejElectron)
{
	Int_t runId  = event->mRunId;
	Float_t vz = event->mVertexZ;
	Float_t vy = event->mVertexY;
	Float_t vx = event->mVertexX;
	Float_t vr = sqrt(vx*vx+vy*vy);
	//west
	Float_t mPlusQx = event->mEtaPlusQx;
	Float_t mPlusQy = event->mEtaPlusQy;
	//east
	Float_t mMinusQx = event->mEtaMinusQx;
	Float_t mMinusQy = event->mEtaMinusQy;
	Int_t mEtaPlusNTrks = event->mEtaPlusNTrks;
	Int_t mEtaMinusNTrks = event->mEtaMinusNTrks;
	Float_t Qx = mPlusQx + mMinusQx; 
	Float_t Qy = mMinusQy + mPlusQy;
	int dayIndex = -99;
	map<Int_t,Int_t>::iterator iter = mTotalDayId.find((runId/1000)%1000);
		if(iter != mTotalDayId.end())
			dayIndex = iter->second;
		else{
			dayIndex = -1;
			cout<<"Can not find the dayId in the mTotalDayId list"<<endl;
		}

	if(mDebug) cout << "before get raw Q" << endl;
	TVector2 mRawQ(Qx,Qy);
	Double_t rawEP = 0.5*mRawQ.Phi();
	if(rawEP<0.) rawEP += TMath::Pi();
	hRawEventPlane->Fill(rawEP);

	//reject the electron for etaplus and etaminus
	if(rejElectron){ //reject the contribution of electron
		for(Int_t i=0;i<current_nE;i++){
			if (current_e[i].Eta() < 0)
			{
				mMinusQx -= current_EQx[i];
				mMinusQy -= current_EQy[i];
				mEtaMinusNTrks = mEtaMinusNTrks-1;

			} else if (current_e[i].Eta() > 0)
			{
				mPlusQx -= current_EQx[i];
				mPlusQy -= current_EQy[i];
				mEtaPlusNTrks = mEtaPlusNTrks-1;
			}
			Qx -= current_EQx[i];
			Qy -= current_EQy[i];
		}
	}

	Double_t eventPlane = -1;
	Double_t eventPlane_East = -1;
	Double_t eventPlane_West = -1;
	TVector2 Q(Qx,Qy);
	TVector2 QPlus(mPlusQx,mPlusQy);
	TVector2 QMinus(mMinusQx,mMinusQy);

	if((Q.Mod())!=0.){
		eventPlane = 0.5*Q.Phi();
		if(eventPlane<0.) eventPlane +=TMath::Pi();
	}
	if((QPlus.Mod())!=0.){
		eventPlane_West = 0.5*QPlus.Phi();
		if(eventPlane_West<0.) eventPlane_West +=TMath::Pi();
	}
	if((QMinus.Mod())!=0.){
		eventPlane_East = 0.5*QMinus.Phi();
		if(eventPlane_East<0.) eventPlane_East +=TMath::Pi();
	}
	hNewEventPlane->Fill(eventPlane);// with electron rejection
	hNewEventPlaneEast->Fill(eventPlane_East);
	hNewEventPlaneWest->Fill(eventPlane_West);

	if(eventPlane<0.) return eventPlane;

	if(mDebug) cout << "before recenter" << endl;
	//********* get recenter number and recenter *********
	Double_t mReCenterQx, mReCenterQy;
	Double_t mReCenterQxEast, mReCenterQyEast;
	Double_t mReCenterQxWest, mReCenterQyWest;
	if(vz>0){
		mReCenterQx = Qx - mEtaPlusNTrks*etapluszplusQx->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszplusQx->GetBinContent(runIndex+1, mCentrality);
		mReCenterQy = Qy - mEtaPlusNTrks*etapluszplusQy->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszplusQy->GetBinContent(runIndex+1, mCentrality);

		mReCenterQxEast = mMinusQx - mEtaMinusNTrks*etaminuszplusQx->GetBinContent(runIndex+1, mCentrality);
		mReCenterQyEast = mMinusQy - mEtaMinusNTrks*etaminuszplusQy->GetBinContent(runIndex+1, mCentrality);

		mReCenterQxWest = mPlusQx - mEtaPlusNTrks*etapluszplusQx->GetBinContent(runIndex+1, mCentrality);
		mReCenterQyWest = mPlusQy - mEtaPlusNTrks*etapluszplusQy->GetBinContent(runIndex+1, mCentrality);
	}
	else{
		mReCenterQx = Qx - mEtaPlusNTrks*etapluszminusQx->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszminusQx->GetBinContent(runIndex+1, mCentrality);
		mReCenterQy = Qy - mEtaPlusNTrks*etapluszminusQy->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszminusQy->GetBinContent(runIndex+1, mCentrality);

		mReCenterQxEast = mMinusQx - mEtaMinusNTrks*etaminuszminusQx->GetBinContent(runIndex+1, mCentrality);
		mReCenterQyEast = mMinusQy - mEtaMinusNTrks*etaminuszminusQy->GetBinContent(runIndex+1, mCentrality);

		mReCenterQxWest = mPlusQx - mEtaPlusNTrks*etapluszminusQx->GetBinContent(runIndex+1, mCentrality);
		mReCenterQyWest = mPlusQy - mEtaPlusNTrks*etapluszminusQy->GetBinContent(runIndex+1, mCentrality);
	}
	if(mDebug) cout << "before cal recentred EP" << endl;

	Double_t recenterEP_noFlat;
	Double_t recenterEP;
	Double_t recenterEPEast;
	Double_t recenterEPWest;
	Double_t recenterEP_2;
	TVector2 *mReCenterQ = new TVector2(mReCenterQx, mReCenterQy);
	TVector2 *mReCenterQWest = new TVector2(mReCenterQxWest, mReCenterQyWest);
    TVector2 *mReCenterQEast = new TVector2(mReCenterQxEast, mReCenterQyEast);

	if(mReCenterQ->Mod() > 0){
		recenterEP = 0.5*mReCenterQ->Phi();
		if(recenterEP<0.) recenterEP += TMath::Pi();
		hReCenterEventPlane->Fill(recenterEP);
	}
	if(mReCenterQWest->Mod() > 0){
		recenterEPWest = 0.5*mReCenterQWest->Phi();
		if(recenterEPWest<0.) recenterEPWest += TMath::Pi();
		hReCenterEventPlaneWest->Fill(recenterEPWest);
	}
	if(mReCenterQEast->Mod() > 0){
		recenterEPEast = 0.5*mReCenterQEast->Phi();
		if(recenterEPEast<0.) recenterEPEast += TMath::Pi();
		hReCenterEventPlaneEast->Fill(recenterEPEast);
	}
	recenterEP_noFlat = recenterEP;

	//for now just using the flat event plane to do the calculation
	hEventPlaneWestvsEast->Fill(recenterEPEast,recenterEPWest);
	EventPlanRes->Fill(mCentrality, cos(2*(recenterEPEast-recenterEPWest)));
	if(mDebug) cout << "before shift" << endl;
	
	//*********  get shift factor and add shift deltaPhi *********
	Float_t shiftCorrcos[mArrayLength];
	Float_t shiftCorrsin[mArrayLength];
	for(Int_t i=0;i<mArrayLength;i++){
		shiftCorrcos[i] = ShiftFactorcos[i]->GetBinContent(dayIndex+1,mCentrality);
		shiftCorrsin[i] = ShiftFactorsin[i]->GetBinContent(dayIndex+1,mCentrality);
	}
	recenterEP_2 = recenterEP;
	Double_t deltaPhi=0;
	Double_t deltaPhi_2=0;
	for(Int_t i=0;i<mArrayLength;i++){
		deltaPhi += 1./(i+1)*(-1.*shiftCorrsin[i]*cos(2.*(i+1)*recenterEP) + shiftCorrcos[i]*sin(2.*(i+1)*recenterEP));
	}
	deltaPhi = deltaPhi/2.;
	if(deltaPhi<0.) deltaPhi += TMath::Pi();
	if(deltaPhi>=TMath::Pi()) deltaPhi -= TMath::Pi();
	recenterEP += deltaPhi;
	if(recenterEP<0.) recenterEP += TMath::Pi();
	if(recenterEP>=TMath::Pi()) recenterEP -= TMath::Pi();
	// hDelta_Psi2->Fill(recenterEP,deltaPhi);
	hFinalEventPlane->Fill(recenterEP);

	deltaPhi_2 = Delta_Psi2->Eval(recenterEP_2);
	if(deltaPhi_2<0.) deltaPhi_2 += TMath::Pi();
	if(deltaPhi_2>=TMath::Pi()) deltaPhi_2 -= TMath::Pi();
	recenterEP_2 += deltaPhi_2;
	if(recenterEP_2<0.) recenterEP_2 += TMath::Pi();
	if(recenterEP_2>=TMath::Pi()) recenterEP_2 -= TMath::Pi();
	hDelta_Psi2->Fill(recenterEP_2,deltaPhi_2);


	if(mDebug) cout << "before final" << endl;
	hDelta_Psi2_FitvsFactor->Fill(recenterEP_2,recenterEP);
	hFinalEventPlane_Fit->Fill(recenterEP_2);
	if (abs(recenterEP_2-recenterEP) > 0.02)
	{
		hLargeDiffEvt_Day->Fill(dayIndex);
		hLargeDiffEvt_vz->Fill(vz);
		hLargeDiffEvt_vr->Fill(vr);
	}
	
	return recenterEP_noFlat;
	// return recenterEP_2;
	// return recenterEP;
}
//____________________________________________________________
//get the costheta* of electron
//Phi index is the costheta* bin
void GetPtPhiCentBin(TLorentzVector pair,TLorentzVector Positron, int _mCentrality,float eventphi,int &ptindex,int &yindex,int &phiindex,int &CentIndex,double &costhe,Bool_t tangent, Int_t Flag ){
	if(mDebug){
		cout<<"in get pT phi bin"<<endl;
	}
	Int_t _PtIndex = -999;
	Int_t _YIndex = -999;
	Int_t _PhiIndex = -999;

	_PtIndex = PtAxis->FindBin(pair.Pt()) - 1;

	//get y
	
	_YIndex = YAxis->FindBin(pair.Rapidity()) -1;

	//get phi
	TVector3 vBetaPhi = -1.0*pair.BoostVector(); 
	Positron.Boost(vBetaPhi);
	TVector3 vPositronRest = Positron.Vect().Unit(); // positron momentum direction in J/psi rest frame	
	if(mDebug){
		cout<<"positron phi:"<<vPositronRest.Phi()<<endl;
	}

	Double_t CosThetaStar = -999.;
	if(!tangent){
		TVector3 Q2(TMath::Sin(eventphi),-1.0*TMath::Cos(eventphi),0.0);
		TVector3 Q2_Unit = Q2.Unit();
		CosThetaStar =  vPositronRest.Dot(Q2_Unit);
	}else{
		TVector3 Q2(TMath::Sin(eventphi),1.0*TMath::Cos(eventphi),0.0);
		TVector3 Q2_Unit = Q2.Unit();
		CosThetaStar =  vPositronRest.Dot(Q2_Unit);
	}

	_PhiIndex = PhiAxis->FindBin(TMath::Abs(CosThetaStar))-1;

	if(mDebug){
		cout<<"Centality:"<<_mCentrality<<endl;
	}
	float _mCent = _mCentrality;
	CentIndex = CentAxis->FindBin(_mCent) -1;
	if(mDebug){
		cout<<"CentIndex:"<<CentIndex<<endl;
	}

	costhe = CosThetaStar;
	ptindex = _PtIndex;
	yindex = _YIndex;
	phiindex = _PhiIndex;
}
//____________________________________________________________
Double_t phiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2)
{
	Double_t pt1 = e1.Pt();
	Double_t eta1 = e1.Eta();
	Double_t phi1 = e1.Phi();

	Double_t pt2 = e2.Pt();
	Double_t eta2 = e2.Eta();
	Double_t phi2 = e2.Phi();

	TVector3 e1Mom,e2Mom;
	if(q1>0&&q2<0){
		e2Mom.SetPtEtaPhi(pt1,eta1,phi1);//e+
		e1Mom.SetPtEtaPhi(pt2,eta2,phi2);//e-
	}else if(q1<0&&q2>0){
		e2Mom.SetPtEtaPhi(pt2,eta2,phi2);//e+
		e1Mom.SetPtEtaPhi(pt1,eta1,phi1);//e-
	}else if(q1==q2&&TMath::Abs(q1)==1){
		Double_t ran = myRandom->Uniform(-1,1);
		if(ran>0){
			e2Mom.SetPtEtaPhi(pt1,eta1,phi1);
			e1Mom.SetPtEtaPhi(pt2,eta2,phi2);
		}
		else{
			e2Mom.SetPtEtaPhi(pt2,eta2,phi2);
			e1Mom.SetPtEtaPhi(pt1,eta1,phi1);
		}
	}else return -1;
	Double_t mN = 0.;
	if(bField<0.) mN = -1.;
	if(bField>0.) mN = 1.;

	TVector3 pu=e1Mom+e2Mom;
	TVector3 pv=e1Mom.Cross(e2Mom);
	TVector3 pw=pu.Cross(pv);
	TVector3 pnz(0.,0.,mN);
	TVector3 pwc=pu.Cross(pnz);

	Double_t angleV = pw.Angle(pwc);

	return angleV;
}
//____________________________________________________________
Double_t calCosTheta(TLorentzVector eVec,TLorentzVector eeVec)
{
	//eVec: positron TLorentzVector  eeVec: ee pair LorentzVector
	TLorentzVector positron(eVec); //positron
	TLorentzVector beam(0., 0., sqrt(pow(96.5,2)-pow(Mproton,2)), 96.5); // UU@193 GeV

	TVector3 dir = eeVec.BoostVector();

	positron.Boost(-1*dir);
	beam.Boost(-1*dir);

	Float_t theta = positron.Angle(beam.Vect());
	return TMath::Cos(theta);
}
//____________________________________________________________
void Polarization(int icharge,int jcharge,TLorentzVector ivector,TLorentzVector jvector){
	TLorentzVector mpositron,melectron,JPSI;
	TLorentzVector Proton1(0.,0.,100.,100.),Proton2(0.,0.,-100.,100.);
	TVector3 XX,YY,ZZ;
	JPSI = ivector+jvector;
	mpositron = icharge>=jcharge? ivector:jvector;
	melectron = jcharge<=icharge? jvector:ivector;
	Int_t NFRAME =4;
	Double_t theta[NFRAME];
	Double_t phi[NFRAME];
	TVector3 XXHX,YYHX,ZZHX;
	ZZHX = JPSI.Vect();
	YYHX = JPSI.Vect().Cross(Proton1.Vect());
	XXHX = YYHX.Cross(ZZHX);

	mpositron.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());
	melectron.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());
	Proton1.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());
	Proton2.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());

	theta[0]= mpositron.Angle(JPSI.Vect());
	phi[0] = TMath::ATan2(mpositron.Vect().Dot(YYHX.Unit()),mpositron.Vect().Dot(XXHX.Unit()));
	electron_theta_hx = melectron.Angle(JPSI.Vect());
	electron_phi_hx = TMath::ATan2(melectron.Vect().Dot(YYHX.Unit()),(melectron.Vect().Dot(XXHX.Unit())));

	ZZ = Proton1.Vect()*(1/(Proton1.Vect()).Mag())-Proton2.Vect()*(1/(Proton2.Vect()).Mag());
	YY = Proton1.Vect().Cross(Proton2.Vect());
	XX = Proton1.Vect()*(1/(Proton1.Vect()).Mag())+Proton2.Vect()*(1/(Proton2.Vect()).Mag());

	theta[1] = mpositron.Angle(ZZ);
	phi[1] = TMath::ATan2(mpositron.Vect().Dot(YY.Unit()),mpositron.Vect().Dot(XX.Unit()));

	positron_theta_hx = theta[0];
	positron_theta_cs = theta[1];
	positron_phi_hx = phi[0];
	positron_phi_cs = phi[1];

	electron_theta_cs = melectron.Angle(ZZ);
	electron_phi_cs = TMath::ATan2(melectron.Vect().Dot(YY.Unit()),melectron.Vect().Dot(XX.Unit()));

	pair_pt = JPSI.Pt();
	pair_eta = JPSI.Eta();
	pair_phi = JPSI.Phi();
	pair_InvM = JPSI.M();

	lepton1_pt = mpositron.Pt();
	lepton1_eta = mpositron.Eta();
	lepton1_phi = mpositron.Phi();
	lepton1_InvM = mpositron.M();

	lepton2_pt = melectron.Pt();
	lepton2_eta = melectron.Eta();
	lepton2_phi = melectron.Phi();
	lepton2_InvM = melectron.M();
}
//____________________________________________________________
void fillHistograms(std::string unlikeOrlike, TLorentzVector JPSI)
{
	if(unlikeOrlike.compare("unlike")==0){
		hPairCosThetaPt->Fill(TMath::Cos(positron_theta_hx),JPSI.Pt());
		hPairPhiPtHX->Fill(positron_phi_hx,JPSI.Pt());
		hPairCosThetaPtCS->Fill(TMath::Cos(positron_theta_cs),JPSI.Pt());
		hPairPhiPtCS->Fill(positron_phi_cs,JPSI.Pt());
		hPairCosThetaPhiPt->Fill(TMath::Cos(positron_theta_hx),positron_phi_hx,JPSI.Pt());	
		hPairCosThetaPhiPtCS->Fill(TMath::Cos(positron_theta_cs),positron_phi_cs,JPSI.Pt());	
	}
	else{
		hPairCosThetaPtBG->Fill(TMath::Cos(positron_theta_hx),JPSI.Pt(),0.5);
		hPairPhiPtHXBG->Fill(positron_phi_hx,JPSI.Pt(),0.5);
		hPairCosThetaPtCSBG->Fill(TMath::Cos(positron_theta_cs),JPSI.Pt(),0.5);
		hPairPhiPtCSBG->Fill(positron_phi_cs,JPSI.Pt(),0.5);
		hPairCosThetaPtBG->Fill(TMath::Cos(electron_theta_hx),JPSI.Pt(),0.5);
		hPairPhiPtHXBG->Fill(electron_phi_hx,JPSI.Pt(),0.5);
		hPairCosThetaPtCSBG->Fill(TMath::Cos(electron_theta_cs),JPSI.Pt(),0.5);
		hPairPhiPtCSBG->Fill(electron_phi_cs,JPSI.Pt(),0.5);
		hPairCosThetaPhiPtBG->Fill(TMath::Cos(positron_theta_hx),positron_phi_hx,JPSI.Pt(),0.5);
		hPairCosThetaPhiPtBG->Fill(TMath::Cos(electron_theta_hx),electron_phi_hx,JPSI.Pt(),0.5);
		hPairCosThetaPhiPtCSBG->Fill(TMath::Cos(positron_theta_cs),positron_phi_cs,JPSI.Pt(),0.5);
		hPairCosThetaPhiPtCSBG->Fill(TMath::Cos(electron_theta_cs),electron_phi_cs,JPSI.Pt(),0.5);
	}
}
//____________________________________________________________
void fill3DHistograms(std::string unlikeOrlike, TLorentzVector JPSI,int i,int j,int pairs){
	if(unlikeOrlike.compare("unlike")==0){
		hPairCosThetaInvMPt->Fill(TMath::Cos(positron_theta_hx),JPSI.M(),JPSI.Pt());
		hPairCosThetaInvMPtCS->Fill(TMath::Cos(positron_theta_cs),JPSI.M(),JPSI.Pt());
		hPairPhiInvMPt->Fill(positron_phi_hx,JPSI.M(),JPSI.Pt());
		hPairPhiInvMPtCS->Fill(positron_phi_cs,JPSI.M(),JPSI.Pt());
	}
	else{	
		hPairCosThetaInvMPtBG->Fill(TMath::Cos(positron_theta_hx),JPSI.M(),JPSI.Pt(),0.5);
		hPairCosThetaInvMPtBG->Fill(TMath::Cos(electron_theta_hx),JPSI.M(),JPSI.Pt(),0.5);
		hPairCosThetaInvMPtCSBG->Fill(TMath::Cos(positron_theta_cs),JPSI.M(),JPSI.Pt(),0.5);
		hPairCosThetaInvMPtCSBG->Fill(TMath::Cos(electron_theta_cs),JPSI.M(),JPSI.Pt(),0.5);
		hPairPhiInvMPtBG->Fill(positron_phi_hx,JPSI.M(),JPSI.Pt(),0.5);
		hPairPhiInvMPtBG->Fill(electron_phi_hx,JPSI.M(),JPSI.Pt(),0.5);
		hPairPhiInvMPtCSBG->Fill(positron_phi_cs,JPSI.M(),JPSI.Pt(),0.5);
		hPairPhiInvMPtCSBG->Fill(electron_phi_cs,JPSI.M(),JPSI.Pt(),0.5);
	}
}
//____________________________________________________________
void bookHistograms()
{
	char buf[500];
	for(int i=0;i<mArrayLength;i++){
		sprintf(buf,"ShiftFactorcos_%d",i);
		ShiftFactorcos[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"ShiftFactorsin_%d",i);
		ShiftFactorsin[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
	}

	hnEvts = new TH1D("hnEvts","hnEvts",6,0.5,6.5);
	hnEvts->GetXaxis()->SetBinLabel(1,"nPicoEvents");
	hnEvts->GetXaxis()->SetBinLabel(2,"|Vz| < 35cm");
	hnEvts->GetXaxis()->SetBinLabel(3,"Vr < 2cm");
	hnEvts->GetXaxis()->SetBinLabel(4,"|TPC_{vz} - VPD_{Vz}| < 10 cm");
	hnEvts->GetXaxis()->SetBinLabel(5,"after nTOFHits rejection");
	hCentrality9 = new TH1F("hCentrality9","hCentrality9;Centrality;Counts",16,0,16);
	hRefMult = new TH1F("hRefMult","hRefMult;dN_{ch}/d#eta;Counts",1000,0,1000);
	hVertexZ = new TH1F("hVertexZ","hVertexZ;TPC VertexZ (cm);Counts",2000,-100,100);
	hVzDiff = new TH1F("hVzDiff","hVzDiff;Vz_{TPC} - Vz_{VPD} (cm);Counts",200,-10,10);
	hVr = new TH1D("hVr","hVr;V_{r} (cm);Counts",500,0,5);
	hBField = new TH1F("hBField","hBField;Magnetic Filed (KiloGauss);Counts",400,-10,10);
	hnTofHitsvsRefMult = new TH2D("hnTofHitsvsRefMult",";RefMult;nTofHits",500,0,500,500,0,500);
	hnTofHitsvsRefMult_noCut = new TH2D("hnTofHitsvsRefMult_noCut",";RefMult;nTofHits",500,0,500,500,0,500); 
  	hnTofHitsvsRefMult_Vz35 = new TH2D("hnTofHitsvsRefMult_Vz35",";RefMult;nTofHits",500,0,500,500,0,500);
  	hVxvsVy = new TH2D("hVxvsVy",";Vx;Vy",100,0,10,100,0,10);
	hRunID = new TH1D("hRunID",";RunID;nCounts",214990,21030025,21245015);
	hTriggerID = new TH1D("hTriggerID",";Trigger ID;nCounts",4,780000-1,780040-1);

	const Int_t    nPtBins   = 100;
	const Double_t ptLow     = 0;
	const Double_t ptHi      = 5;
	const Int_t    nMassBins = 400;
	const Double_t massLow   = 0;
	const Double_t massHi    = 4;

	//eventPlane
	hRawEventPlane = new TH1F("hRawEventPlane","hRawEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hNewEventPlane = new TH1F("hNewEventPlane","hNewEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hNewEventPlaneEast = new TH1D("hNewEventPlaneEast","hNewEventPlaneEast;Reaction Plane East (rad); Counts",300,0,TMath::Pi());
	hNewEventPlaneWest = new TH1D("hNewEventPlaneWest","hNewEventPlaneWest;Reaction Plane West (rad); Counts",300,0,TMath::Pi());
	hReCenterEventPlane = new TH1F("hReCenterEventPlane","hReCenterEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hReCenterEventPlaneEast = new TH1D("hReCenterEventPlaneEast","hReCenterEventPlaneEast;Reaction Plane East (rad); Counts",300,0,TMath::Pi());
	hReCenterEventPlaneWest = new TH1D("hReCenterEventPlaneWest","hReCenterEventPlaneWest;Reaction Plane West (rad); Counts",300,0,TMath::Pi());
	hEventPlaneWestvsEast = new TH2D("hEventPlaneWestvsEast","hEventPlaneWestvsEast; EP east; EP west",300,0,TMath::Pi(),300,0,TMath::Pi());
	hFinalEventPlane = new TH1F("hFinalEventPlane","hFinalEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hFinalEventPlaneEast = new TH1D("hFinalEventPlaneEast","hFinalEventPlaneEast;Reaction Plane East (rad); Counts",300,0,TMath::Pi());
	hFinalEventPlaneWest = new TH1D("hFinalEventPlaneWest","hFinalEventPlaneWest;Reaction Plane West (rad); Counts",300,0,TMath::Pi());
	hFinalEventPlane_Fit = new TH1F("hFinalEventPlane_Fit","hFinalEventPlane_Fit;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hDelta_Psi2 = new TH2F("hDelta_Psi2","hDelta_Psi2;recenter #Psi_{2};#Delta#Psi_{2}",300,0,TMath::Pi(),600,-TMath::Pi()-0.1,TMath::Pi()+0.1);
	hDelta_Psi2_FitvsFactor = new TH2F("hDelta_Psi2_FitvsFactor","hDelta_Psi2_FitvsFactor;Fit #Delta#Psi_{2};Factor #Delta#Psi_{2}",300,0-0.2,TMath::Pi()+0.2,300,0-0.2,TMath::Pi()+0.2);
	hLargeDiffEvt_Day = new TH1D("hLargeDiffEvt_Day","hLargeDiffEvt_Day;Day Index; nCounts",mTotalDay+3,0,mTotalDay+3);
	hLargeDiffEvt_vz = new TH1D("hLargeDiffEvt_vz","hLargeDiffEvt_vz;",1200,-60,60);
	hLargeDiffEvt_vr = new TH1D("hLargeDiffEvt_vr","hLargeDiffEvt_vr;",500,0,5);
	EventPlanRes = new TProfile("EventPlanRes","EventPlanRes",10,-0.5,9.5);
	


	hInclusiveEPhivsPt = new TH2F("hInclusiveEPhivsPt","hInclusiveEPhivsPt;q*p_{T} (GeV/c); #phi",200,-10,10,600,-TMath::Pi(),TMath::Pi());
	hExclusiveEPhivsPt = new TH2F("hExclusiveEPhivsPt","hExclusiveEPhivsPt;q*p_{T} (GeV/c); #phi",200,-10,10,600,-TMath::Pi(),TMath::Pi());
	hCut3EPhivsPt = new TH2F("hCut3EPhivsPt","hCut3EPhivsPt;q*p_{T} (GeV/c); #phi",2000,-10,10,1800,-TMath::Pi(),TMath::Pi());

	hnEMinusvsEPlus = new TH2F("hnEMinusvsEPlus","hnEMinusvsEPlus;# e^{+};# e^{-}",30,0,30,30,0,30);
	hnSigmaEvsP = new TH2F("hnSigmaEvsP","; P (GeV/c); n#sigma_{e}",600,0,6,1000,-5,5);

	//with phiV cut
	hULMvsPt = new TH2F("hULMvsPt","hULMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hLPosMvsPt = new TH2F("hLPosMvsPt","hLPosMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hLNegMvsPt = new TH2F("hLNegMvsPt","hLNegMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hMixULMvsPt = new TH2F("hMixULMvsPt","hMixULMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hMixLPosMvsPt = new TH2F("hMixLPosMvsPt","hMixLPosMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hMixLNegMvsPt = new TH2F("hMixLNegMvsPt","hMixLNegMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);

	//add centrality dimension
	hULMvsPtCen = new TH3F("hULMvsPtCen","hULMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hLPosMvsPtCen = new TH3F("hLPosMvsPtCen","hLPosMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hLNegMvsPtCen = new TH3F("hLNegMvsPtCen","hLNegMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixULMvsPtCen = new TH3F("hMixULMvsPtCen","hMixULMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixLPosMvsPtCen = new TH3F("hMixLPosMvsPtCen","hMixLPosMvsPtCen;p_{T} (GeV/c);Centrality;Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixLNegMvsPtCen = new TH3F("hMixLNegMvsPtCen","hMixLNegMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	// hULMvsPhiCen = new TH3F("hULMvsPhiCen","hULMvsPhiCen;Phi;Centrality;M_{ee} (GeV/c^{2})",60,-3.14-0.1,3.14+0.1,16,0,16,500,massLow,massHi);
	// hLPosMvsPhiCen = new TH3F("hLPosMvsPhiCen","hLPosMvsPhiCen;Phi;Centrality;M_{ee} (GeV/c^{2})",60,-3.14-0.1,3.14+0.1,16,0,16,500,massLow,massHi);
	// hLNegMvsPhiCen = new TH3F("hLNegMvsPhiCen","hLNegMvsPhiCen;Phi;Centrality;M_{ee} (GeV/c^{2})",60,-3.14-0.1,3.14+0.1,16,0,16,500,massLow,massHi);
	hCellIDDiff = new TH1D("hCellIDDiff","hCellIDDiff;ID Diff ;counts;",32,-16,16);
	hdEdxvsP = new TH2D();
	hPt_Electron = new TH1D("hPt_Electron",";p_{T};Counts",1000,0,10);
	hPt_Positron = new TH1D("hPt_Positron",";p_{T};Counts",1000,0,10);

	hULMvsPtCen->Sumw2();
	hLPosMvsPtCen->Sumw2();
	hLNegMvsPtCen->Sumw2();
	hMixULMvsPtCen->Sumw2();
	hMixLPosMvsPtCen->Sumw2();
	hMixLNegMvsPtCen->Sumw2();

	//add the histograms for the dilepton polarization, in the mass region 0.2-1.1
	hPairPhiPt = new TH2F("hPairPhiPt","Pair Phi vs #Phi;#Phi;P_{T} GeV/c",360,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairPhiPt->Sumw2();
	hPairPhiPtBG = new TH2F("hPairPhiPtBG","Pair Phi vs #Phi;#Phi;P_{T} GeV/c",360,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairPhiPtBG->Sumw2();

	hPairCosThetaPt = new TH2F("hPairCosThetaPt","Pair Pt vs Cos(#theta); Cos(#theta); Pair Pt;",10,-1,1,120,0,30);
	hPairCosThetaPt->Sumw2();
	hPairCosThetaPtBG = new TH2F("hPairCosThetaPtBG","Pair Pt vs Cos(#theta); Cos(#theta); Pair Pt;",10,-1,1,120,0,30);
	hPairCosThetaPtBG->Sumw2();

	hPairPhiPtHX = new TH2F("hPairPhiPtHX","Pair Pt vs #phi;#phi;Pair Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairPhiPtHX->Sumw2();
	hPairPhiPtHXBG = new TH2F("hPairPhiPtHXBG","Pair Pt vs #phi;#phi;Pair Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairPhiPtHXBG->Sumw2();

	hPairCosThetaPtCS = new TH2F("hPairCosThetaPtCS","Pair Pt vs Cos(#theta);Cos(#theta);Pair Pt",10,-1,1,120,0,30);
	hPairCosThetaPtCS->Sumw2();
	hPairCosThetaPtCSBG = new TH2F("hPairCosThetaPtCSBG","Pair Pt vs Cos(#theta);Cos(#theta);Pair Pt",10,-1,1,120,0,30);
	hPairCosThetaPtCSBG->Sumw2();

	hPairPhiPtCS = new TH2F("hPairPhiPtCS","Pair Pt vs #phi;#phi;Pair Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairPhiPtCS->Sumw2();
	hPairPhiPtCSBG = new TH2F("hPairPhiPtCSBG","Pair Pt vs #phi;#phi;Pair Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairPhiPtCSBG->Sumw2();

	hPairCosThetaInvMPt = new TH3F("hPairCosThetaInvMPt","hPairCosThetaInvMPt",20,-1,1,250,0,2.5,50,0,5);
	hPairCosThetaInvMPtBG= new TH3F("hPairCosThetaInvMPtBG","hPairCosThetaInvMPtBG",20,-1,1,250,0,2.5,50,0,5);
	hPairCosThetaInvMPtCS= new TH3F("hPairCosThetaInvMPtCS","hPairCosThetaInvMPtCS",20,-1,1,250,0,2.5,50,0,5);
	hPairCosThetaInvMPtCSBG= new TH3F("hPairCosThetaInvMPtCSBG","hPairCosThetaInvMPtCSBG",20,-1,1,250,0,2.5,50,0,5);

	hPairPhiInvMPt = new TH3F("hPairPhiInvMPt","hPairPhiInvMPt",40,-TMath::Pi(),TMath::Pi(),250,0,2.5,50,0,5);
	hPairPhiInvMPtBG = new TH3F("hPairPhiInvMPtBG","hPairPhiInvMPtBG",40,-TMath::Pi(),TMath::Pi(),250,0,2.5,50,0,5);
	hPairPhiInvMPtCS = new TH3F("hPairPhiInvMPtCS","hPairPhiInvMPtCS",40,-TMath::Pi(),TMath::Pi(),250,0,2.5,50,0,5);
	hPairPhiInvMPtCSBG = new TH3F("hPairPhiInvMPtCSBG","hPairPhiInvMPtCSBG",40,-TMath::Pi(),TMath::Pi(),250,0,2.5,50,0,5);

	hPairCosThetaInvMPt->Sumw2();
	hPairCosThetaInvMPtBG->Sumw2();
	hPairCosThetaInvMPtCS->Sumw2();
	hPairCosThetaInvMPtCSBG->Sumw2();

	hPairPhiInvMPt->Sumw2();
	hPairPhiInvMPtBG->Sumw2();
	hPairPhiInvMPtCS->Sumw2();
	hPairPhiInvMPtCSBG->Sumw2();

	hPairCosThetaPhiPt = new TH3F("hPairCosThetaPhiPt","hPairCosThetaPhiPt;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairCosThetaPhiPtBG = new TH3F("hPairCosThetaPhiPtBG","hPairCosThetaPhiPtBG;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairCosThetaPhiPt->Sumw2();
	hPairCosThetaPhiPtBG->Sumw2();

	hPairCosThetaPhiPtCS = new TH3F("hPairCosThetaPhiPtCS","hPairCosThetaPhiPtCS;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairCosThetaPhiPtCSBG = new TH3F("hPairCosThetaPhiPtCSBG","hPairCosThetaPhiPtCSBG;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hPairCosThetaPhiPtCS->Sumw2();
	hPairCosThetaPhiPtCSBG->Sumw2();

	hCosthetastar = new TH1D("hCosthetastar","costhetastar; cos(#theta^{*});coutns",200,-1,1);
	hMixCosthetastar = new TH1D("hMixCosthetastar","costhetastar; cos(#theta^{*});coutns",200,-1,1);

	//Mass vs Pt
	hULMvsPtT = new TH2F("hULMvsPtT","hULMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hULMvsPtTW = new TH2F("hULMvsPtTW","hULMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hRapdity = new TH1F("hRapdity","Rapdity",100,-2,2);
	for(int i=0; i<mCenBins;i++){
		for(int j=0; j<mPtBins;j++){
			for(int k=0; k<mPhiBins; k++){
				hULM[i][j][k] = new TH1F(Form("hULM_%d_%d_%d",i,j,k),"hULM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);// cent, pT, costheta
				hLSPlusM[i][j][k] = new TH1F(Form("hLSPlusM_%d_%d_%d",i, j, k), "hLSPlus;M_{ee} (GeV/c^{2})", nMassBins, massLow, massHi);
				hLSMinusM[i][j][k] = new TH1F(Form("hLSMinusM_%d_%d_%d",i,j,k),"hLSMinusM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);
				hMixULM[i][j][k] = new TH1F(Form("hMixULM_%d_%d_%d",i,j,k),"hMixULM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);
				hMixLSPosM[i][j][k] = new TH1F(Form("hMixLSPosM_%d_%d_%d",i,j,k),"hMixLSPosM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);
				hMixLSNegM[i][j][k] = new TH1F(Form("hMixLSNegM_%d_%d_%d",i,j,k),"hMixLSNegM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);
			}
		}
		for(int j=0; j<mYBins;j++){
			for(int k=0; k<mPhiBins; k++){
				hULYM[i][j][k] = new TH1F(Form("hULYM_%d_%d_%d",i,j,k),"hULM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);// cent, rapidity, costheta
				hLSPlusYM[i][j][k] = new TH1F(Form("hLSPlusYM_%d_%d_%d",i,j,k),"hLSPlusYM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);// cent, rapidity, costheta
				hLSMinusYM[i][j][k] = new TH1F(Form("hLSMinusYM_%d_%d_%d",i,j,k),"hLSMinusYM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);// cent, rapidity, costheta
				hMixULYM[i][j][k] = new TH1F(Form("hMixULYM_%d_%d_%d",i,j,k),"hMixULYM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);// cent, rapidity, costheta
				hMixLSPosYM[i][j][k] = new TH1F(Form("hMixLSPosYM_%d_%d_%d",i,j,k),"hMixLSPosYM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);// cent, rapidity, costheta
				hMixLSMinusYM[i][j][k] = new TH1F(Form("hMixLSMinusYM_%d_%d_%d",i,j,k),"hMixLSMinusYM;M_{ee} (GeV/c^{2})",nMassBins,massLow,massHi);// cent, rapidity, costheta


			}
		}
	}



}
//=======================================================================================
void writeHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.histo.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	TFile *mFile = new TFile(buf,"recreate");
	mFile->cd();

	//in passEvent function
	hnEvts->Write();
	hCentrality9->Write();
	hRefMult->Write();
	hVertexZ->Write();
	hVzDiff->Write();
	hVr->Write();
	hBField->Write();
	hnTofHitsvsRefMult_noCut->Write();
    hnTofHitsvsRefMult_Vz35->Write();
	hnTofHitsvsRefMult->Write();
	hRunID->Write();
	hTriggerID->Write();
	cout << "writing passevent done" << endl;

	//eventPlane
	hRawEventPlane->Write();
	hNewEventPlane->Write();
	hFinalEventPlane->Write();
	hFinalEventPlane_Fit->Write();
	hReCenterEventPlane->Write();
	hDelta_Psi2->Write();
	hDelta_Psi2_FitvsFactor->Write();
	hLargeDiffEvt_Day->Write();
	hLargeDiffEvt_vz->Write();
	hLargeDiffEvt_vr->Write();
	cout << "writing event plane done" << endl;
	
	hInclusiveEPhivsPt->Write();
	hExclusiveEPhivsPt->Write();
	hCut3EPhivsPt->Write();

	hnEMinusvsEPlus->Write();
	cout << "writing electron done" << endl;

	//angleV 

	//without phiV cut
	// hULMvsPtwophiV->Write();
	// hLPosMvsPtwophiV->Write();
	// hLNegMvsPtwophiV->Write();
	// hMixULMvsPtwophiV->Write();
	// hMixLPosMvsPtwophiV->Write();
	// hMixLNegMvsPtwophiV->Write();
	cout << "writing angle V done" << endl;

	//with phiV cut
	hULMvsPt->Write();
	hLPosMvsPt->Write();
	hLNegMvsPt->Write();
	hMixULMvsPt->Write();
	hMixLPosMvsPt->Write();
	hMixLNegMvsPt->Write();

	cout << "writing phiV done" << endl;
	
	//add centrality dimension
	hULMvsPtCen->Write();
	hLPosMvsPtCen->Write();
	hLNegMvsPtCen->Write();
	hMixULMvsPtCen->Write();
	hMixLPosMvsPtCen->Write();
	hMixLNegMvsPtCen->Write();
	// hULMvsPhiCen->Write();
  	// hLPosMvsPhiCen->Write();
  	// hLNegMvsPhiCen->Write();
	hCellIDDiff->Write();
  	hnSigmaEvsP->Write();
	hVxvsVy->Write();
	hdEdxvsP->Write();
	hPt_Electron->Write();
	hPt_Positron->Write();
	cout << "writing 3D done" << endl;

	hPairPhiPt->Write();
	hPairPhiPtBG->Write();
	hPairCosThetaPt->Write();
	hPairCosThetaPtBG->Write();
	hPairPhiPtHX->Write();
	hPairPhiPtHXBG->Write();
	hPairCosThetaPtCS->Write();
	hPairCosThetaPtCSBG->Write();
	hPairPhiPtCS->Write();
	hPairPhiPtCSBG->Write();
	hPairCosThetaInvMPt->Write();
	hPairCosThetaInvMPtBG->Write();
	hPairCosThetaInvMPtCS->Write();
	hPairCosThetaInvMPtCSBG->Write();
	hPairPhiInvMPt->Write();
	hPairPhiInvMPtBG->Write();
	hPairPhiInvMPtCS->Write();
	hPairPhiInvMPtCSBG->Write();
	hPairCosThetaPhiPt->Write();
	hPairCosThetaPhiPtBG->Write();
	hPairCosThetaPhiPtCS->Write();
	hPairCosThetaPhiPtCSBG->Write();
	cout << "writing loacl Polarization done" << endl;

	hNewEventPlaneEast->Write();
	hNewEventPlaneWest->Write();
	hReCenterEventPlaneWest->Write();
	hReCenterEventPlaneEast->Write();
	hFinalEventPlaneWest->Write();
	hFinalEventPlaneEast->Write();
	hEventPlaneWestvsEast->Write();
	hCosthetastar->Write();
	hMixCosthetastar->Write();
	EventPlanRes->Write();

	cout << "writing EP done" << endl;
	for(int i=0; i<mCenBins;i++){
		for(int j=0; j<mPtBins;j++){
			for(int k=0; k<mPhiBins; k++){
				hULM[i][j][k]->Write();
				hLSPlusM[i][j][k]->Write();
				hLSMinusM[i][j][k]->Write();
				hMixULM[i][j][k]->Write();
				hMixLSPosM[i][j][k]->Write();
				hMixLSNegM[i][j][k]->Write();
			}
		}
		for(int j=0; j<mYBins;j++){
			for(int k=0; k<mPhiBins; k++){
				hULYM[i][j][k]->Write();
				hLSPlusYM[i][j][k]->Write();
				hLSMinusYM[i][j][k]->Write();
				hMixULYM[i][j][k]->Write();
				hMixLSPosYM[i][j][k]->Write();
				hMixLSMinusYM[i][j][k]->Write();
			}
		}
	}

	// hULCosThetavsMvsCen->Write();
	// hLPosCosThetavsMvsCen->Write();
	// hLNegCosThetavsMvsCen->Write();
	// hMixULCosThetavsMvsCen->Write();
	// hMixLPosCosThetavsMvsCen->Write();
	// hMixLNegCosThetavsMvsCen->Write();

	// hULePtvsMeevsCen->Write();
	// hLSePtvsMeevsCen->Write();
	// hMixULePtvsMeevsCen->Write();
	// hMixLSePtvsMeevsCen->Write();
}
//==============================================================================================
Bool_t Init()
{
	cout<<endl;

	ifstream indata;

	indata.open("/star/u/wangzhen/run20/Dielectron/DataQA/mTotalDayList.dat");
	mTotalDayId.clear();
	if(indata.is_open()){
		cout<<"read in day number list and recode day number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata>>oldId){
			mTotalDayId[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the day number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	indata.open("/star/u/wangzhen/run20/Dielectron/DataQA/mTotalRunList.dat");
	mTotalRunId.clear();
	if(indata.is_open()){
		cout<<"read in total run number list and recode run number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata>>oldId){
			mTotalRunId[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	//read in bad run for 580001 and 580021
	ifstream indata_001;
	indata_001.open("/star/u/wangzhen/run20/Dielectron/BadRunList/BadRunList.dat");
	mBadRunId_001.clear();
	if(indata_001.is_open()){
		cout<<"read in bad run list for 9.2 GeV Au+Au GeV ";
		Int_t oldId;
		Int_t newId=0;
		while(indata_001>>oldId){
			mBadRunId_001[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total bad run  list !!!"<<endl;
		return kFALSE;
	}
	indata_001.close();

	PtAxis = new TAxis(mPtBins, mPairPtCut);
	YAxis = new TAxis(mYBins, -1,1);
	PhiAxis = new TAxis(mPhiBins, 0, 1);
	CentAxis = new TAxis(mCenBins,mCentCut);


	cout<<"bad run for trigger 580001"<<endl;
	for(map<Int_t,Int_t>::iterator iter=mBadRunId_001.begin();iter!=mBadRunId_001.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;


	TFile *fReCenter = TFile::Open("/star/u/wangzhen/run20/Dielectron/FlatEvtPlane/reCenter/output_all/reCenter.root");
	if(fReCenter->IsOpen()){
		cout<<"read in re-center root file ...";
		etapluszplusQx   = (TProfile2D *)fReCenter->Get("etapluszplusQx");
		etapluszminusQx  = (TProfile2D *)fReCenter->Get("etapluszminusQx");
		etaminuszplusQx  = (TProfile2D *)fReCenter->Get("etaminuszplusQx");
		etaminuszminusQx = (TProfile2D *)fReCenter->Get("etaminuszminusQx");
		etapluszplusQy   = (TProfile2D *)fReCenter->Get("etapluszplusQy");
		etapluszminusQy  = (TProfile2D *)fReCenter->Get("etapluszminusQy");
		etaminuszplusQy  = (TProfile2D *)fReCenter->Get("etaminuszplusQy");
		etaminuszminusQy = (TProfile2D *)fReCenter->Get("etaminuszminusQy");
	}

	TFile *fShift = TFile::Open("/star/u/wangzhen/run20/Dielectron/FlatEvtPlane/shift/output_all/shift.histo.root");
	if(fShift->IsOpen()){
		cout<<"read in shiftfactor root file ...";
		for(int i=0;i<mArrayLength;i++){
			ShiftFactorcos[i] = (TProfile2D*)fShift->Get(Form("shiftfactorcos_%d",i));
			ShiftFactorsin[i] = (TProfile2D*)fShift->Get(Form("shiftfactorsin_%d",i));
		}
		cout<<" [OK]"<<endl;
	}

	Pileuplimit = new TF1("Pileuplimit","0.7*x-10",0,1000);
	phiVcut=new TF1("phiVcut","0.84326*exp((-49.4819)*x)+(-0.996609)*x+(0.19801)",0.,1.0); //jie's cut
	PileupUplimit = new TF1("PileupUplimit","pol6",0,340);
	PileupUplimit->SetParameters(7.14109e+00,3.24086e+00,-5.75451e-02,8.41265e-04,-5.77820e-06,1.82626e-08,-2.17213e-11);
	PileupLowlimit = new TF1("PileupLowlimit","pol6",0,340);
	PileupLowlimit->SetParameters(-6.33246e+00,7.90568e-02,3.03279e-02,-5.03738e-04,3.82206e-06,-1.30813e-08,1.64832e-11);
	Delta_Psi2 = new TF1("Delta_Psi2","0.5*( 2*[0]*sin(2*x)-2*[1]*cos(2*x)+[3]*sin(4*x)-[2]*cos(4*x) )",-TMath::Pi(),TMath::Pi());
	Delta_Psi2->SetParNames("<cos2#Psi_{2}>","<sin2#Psi_{2}>","<cos4#Psi_{2}>","<sin4#Psi_{2}>");
	Delta_Psi2->SetParameters(0.001461,0.000840,0.002069,0.002289);
	f_upper->SetParameters(6.32816,0.689232,-0.00185181,6.31563e-06,-8.29481e-09);
	f_lower->SetParameters(-5.20165,0.144438,0.00186397,-1.28471e-05,4.28608e-08);


	cout<<"Initialization DONE !!!"<<endl;
	cout<<endl;

	return kTRUE;
}
