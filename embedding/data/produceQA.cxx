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

#include "EVENT.h"
#include "./StRefMultCorr/StRefMultCorr.h"
#include "./StRefMultCorr/CentralityMaker.h"
#include "cuts.h"

using namespace std;
#endif

StRefMultCorr *refMultCorrUtil;

map<Int_t,Int_t> mBadRunId;
Float_t par[4][4];
Float_t parErr[4][4];

const Int_t mMaxElectrons = 100;
Int_t current_nEPlus;
Int_t current_nEMinus;
TLorentzVector current_ePlus[mMaxElectrons];
Int_t current_ePlusTrkId[mMaxElectrons];
Bool_t current_phePlusTag[mMaxElectrons];
TLorentzVector current_eMinus[mMaxElectrons];
Int_t current_eMinusTrkId[mMaxElectrons];
Bool_t current_pheMinusTag[mMaxElectrons];

bool Init();
void bookHistograms();
bool passEvent(EVENT* event);
bool passTrack(EVENT* event, Int_t i);
void photonE(EVENT* event);
void writeHistograms(char* outFile);

//define functions and histograms
TF1 *funPosHi;
TF1 *funPosLow;
TF1 *funNegHi;
TF1 *funNegLow;

TH2F *hVyvsVx;
TH2F *hVyvsVz;
TH2F *hRefMultvsRefMultCorr;

TH2F *hULMvsPt;
TH2F *hLPosMvsPt;
TH2F *hLNegMvsPt;

TH3F *hEtavsPtQ;
TH3F *hPhivsPtQ;
TH3F *hYvsPtQ;
TH3F *hNHitsFitvsPtQ;
TH3F *hNHitsFitvsEtaQ;
TH3F *hNHitsFitvsPhiQ;
TH3F *hNHitsPossvsPtQ;
TH3F *hNHitsPossvsEtaQ;
TH3F *hNHitsPossvsPhiQ;
TH3F *hNHitsDedxvsPtQ;
TH3F *hDcavsPtQ;
TH3F *hDedxvsPQ;
TH3F *hDedxvsEtaQ;
TH3F *hDedxvsPhiQ;
TH3F *hNSigmaEvsPQ;
TH3F *hNSigmaEvsEtaQ;
TH3F *hNSigmaEvsPhiQ;
TH3F *hNSigmaPivsPQ;
TH3F *hNSigmaKvsPQ;
TH3F *hNSigmaPvsPQ;

TH3F *hPiEtavsPtQ;
TH3F *hPiPhivsPtQ;
TH3F *hPiYvsPtQ;
TH3F *hPiNHitsFitvsPtQ;
TH3F *hPiNHitsFitvsEtaQ;
TH3F *hPiNHitsFitvsPhiQ;
TH3F *hPiNHitsPossvsPtQ;
TH3F *hPiNHitsPossvsEtaQ;
TH3F *hPiNHitsPossvsPhiQ;
TH3F *hPiNHitsDedxvsPtQ;
TH3F *hPiDcavsPtQ;
TH3F *hPiDedxvsPQ;
TH3F *hPiDedxvsEtaQ;
TH3F *hPiDedxvsPhiQ;
TH3F *hPiNSigmaEvsPQ;
TH3F *hPiNSigmaPivsPQ;
TH3F *hPiNSigmaPivsEtaQ;
TH3F *hPiNSigmaPivsPhiQ;
TH3F *hPiNSigmaKvsPQ;
TH3F *hPiNSigmaPvsPQ;

// wangzhen add it for check in 3D
TH3F *hPENHitsFitsvsPtvsEtaMinus;
TH3F *hPENHitsFitsvsPtvsPhiMinus;
TH3F *hPENHitsFitsvsPtvsEtaPlus;
TH3F *hPENHitsFitsvsPtvsPhiPlus;
TH3F *hPENHitsDedxvsPtvsEtaMinus;
TH3F *hPENHitsDedxvsPtvsPhiMinus;
TH3F *hPENHitsDedxvsPtvsEtaPlus;
TH3F *hPENHitsDedxvsPtvsPhiPlus;
TH3F *hPEDcavsPtvsEtaMinus;
TH3F *hPEDcavsPtvsPhiMinus;
TH3F *hPEDcavsPtvsEtaPlus;
TH3F *hPEDcavsPtvsPhiPlus;


//to check track origin Z and phi distribution
TH1D *hMinusOriginZ;
TH1D *hMinusOriginPhi;
TH2D *hMinusOriginZvsPhi;
TH1D *hPlusOriginZ;
TH1D *hPlusOriginPhi;
TH2D *hPlusOriginZvsPhi;
TH2D* hPhotonEtavsPhi;


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
				cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
				chain->Add(filename);
				ifile++;
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
	}
	//refMultCorrUtil = new StRefMultCorr();

	//+-------------+
	//| loop events |
	//+-------------+
	EVENT *event = new EVENT(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;
	for(int i=0;i<nEvts;i++){

		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;
		if (event==NULL) continue;
		event->GetEntry(i);

		Int_t runId = event->mRunId;
		map<Int_t,Int_t>::iterator iter = mBadRunId.find(runId);
		if(iter != mBadRunId.end()) continue;

		if(!passEvent(event)) continue;

		current_nEPlus = 0;
		current_nEMinus = 0;
		Int_t nTrks = event->mNTrks;
		for(int j=0;j<nTrks;j++) passTrack(event,j);

		photonE(event);

	}

	writeHistograms(outFile);
	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
bool passEvent(EVENT* event)
{
	bool eventflag=kFALSE;
	for(Int_t i=0;i<event->mNTrigs;i++){
		if(event->mTrigId[i] == 780010) eventflag=kTRUE; //vpd-zdce-tac-protected
		if(event->mTrigId[i] == 780020) eventflag=kTRUE; //vpd-zdce-tac-protected
		//if(event->mTrigId[i] == 810030) eventflag=kTRUE; //vpd-zdce-tac-protected
		//if(event->mTrigId[i] == 810040) eventflag=kTRUE; //vpd-zdce-tac-protected
	}
	if(!eventflag) return kFALSE;

	Float_t runId = event->mRunId;
	Float_t zdcRate = event->mZDCRate;
	Int_t   refMult = event->mRefMult;
	Int_t mnTOFMatch = event->mnTOFMatch;
  Float_t vx = event->mVertexX;
	Float_t vy = event->mVertexY;
	Float_t vz = event->mVertexZ;
	Float_t vr = sqrt(vx*vx+vy*vy);
	Float_t vpdVz = event->mVpdVz;
	Float_t vzDiff = vz - vpdVz;
  int mCentrality = event->mCentrality;

	if(TMath::Abs(vx)<1.e-5 && TMath::Abs(vy)<1.e-5 && TMath::Abs(vz)<1.e-5) return kFALSE;
	if(vr>=mVrCut) return kFALSE;
	if(TMath::Abs(vz)>=35) return kFALSE;//vz should also be in the range listed in the parameters file to do the refMult correction
	// if(TMath::Abs(vzDiff)>=mVzDiffCut) return kFALSE;
	

  refMultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
	refMultCorrUtil->init(runId);
  Bool_t isPileUpEvt_Cen = !refMultCorrUtil->passnTofMatchRefmultCut(refMult*1.0,mnTOFMatch*1.0);
  if(isPileUpEvt_Cen) return kFALSE;
	refMultCorrUtil->initEvent(refMult,vz,zdcRate);
	Double_t refMultCor = refMultCorrUtil->getRefMultCorr();//if you want to call the getRefMultCorr() with no argument, it must be called after initEvent()

	hVyvsVx->Fill(vx,vy);
	hVyvsVz->Fill(vz,vy);
	//hRefMultvsRefMultCorr->Fill(refMultCor,refMult);

	return kTRUE;
}
//______________________________________________________________
bool passTrack(EVENT* event, Int_t i)
{
	Int_t charge = event->mCharge[i];
	Int_t nHitsFit = event->mNHitsFit[i];
	Int_t nHitsPoss = event->mNHitsPoss[i];
	Int_t nHitsDedx = event->mNHitsDedx[i];
	Float_t ratio = 1.0*nHitsFit/nHitsPoss;
	Float_t dedx = event->mDedx[i];
	Float_t nSigmaE = event->mNSigmaE[i];
	/*	Float_t nSigmaPi = event->mNSigmaPi[i];
		Float_t nSigmaK = event->mNSigmaK[i];
		Float_t nSigmaP = event->mNSigmaP[i];*/
	Float_t dca = event->mDca[i];
	Float_t pt = event->mPt[i];
	Float_t eta = event->mEta[i];
	Float_t phi = event->mPhi[i];
	//Float_t gPt = event->mgPt[i];
	//Float_t gEta = event->mgEta[i];
	//Float_t gPhi = event->mgPhi[i];
	//Float_t gOriginX = event->mgOriginX[i]/100.;
//Float_t gOriginY = event->mgOriginY[i]/100.;
	Float_t gOriginZ = event->mgOriginZ[i];
	Float_t beta2TOF = event->mBeta2TOF[i];

	TVector3 mom;
	mom.SetPtEtaPhi(pt,eta,phi);
	Float_t p = mom.Mag();

	if(pt<mTpcePtCut[0] || pt>mTpcePtCut[1]) return kFALSE;
	if(pt<mTpcePtCut[0]) return kFALSE;
	// if(nHitsFit<15) return kFALSE;
	if(nHitsFit<15) return kFALSE;
	if(ratio<mTpceNHitsFitRatioCut) return kFALSE;
	if(nHitsDedx<10) return kFALSE;
	// if(nHitsDedx<mTpceNHitsDedxCut) return kFALSE;
	if(dca>3) return kFALSE;
	//if(dca>mTpceDcaCut) return kFALSE;
	if(TMath::Abs(eta)>mTpceEtaCut) return kFALSE;
	// if(TMath::Abs(eta)>mTpceEtaCut) return kFALSE;

	/*	if(charge*pt>0. && eta>0.
		&& phi<=funPosHi->Eval(charge*pt)
		&& phi>=funPosLow->Eval(charge*pt)
		) return kFALSE;

		if(charge*pt<0. && eta>0.
		&& phi<=funNegHi->Eval(charge*pt)
		&& phi>=funNegLow->Eval(charge*pt)
		) return kFALSE;
		*/
	//Float_t mTpceNSigmaECutLow;
	//if(p<1.){
	//	mTpceNSigmaECutLow = (mTpceNSigmaECut[0]+2)/(1.-mTpcePtCut[0])*(p-mTpcePtCut[0]) - 2;
	//}else{
	//	mTpceNSigmaECutLow = mTpceNSigmaECut[0];
	//}
	//if(nSigmaE<mTpceNSigmaECutLow || nSigmaE>mTpceNSigmaECut[1]) return kFALSE;

	//just using pico track quality cuts and nSigmaPi cut to select Pion

	Bool_t piFlag = kFALSE;
	Float_t mSquare;
	if(beta2TOF>0.) mSquare = p*p*(1-beta2TOF*beta2TOF)/beta2TOF/beta2TOF;
	else mSquare = -999;

	if( TMath::Abs(mSquare-0.019)<0.003) piFlag = kTRUE;

	/*if(TMath::Abs(nSigmaPi)<=2.0 && piFlag){
	  TLorentzVector pion(0,0,0,0);
	  pion.SetPtEtaPhiM(pt,eta,phi,Mpion);
	  if(charge>0){
	  hPiEtavsPtQ->Fill(0.5,pt,eta);
	  hPiPhivsPtQ->Fill(0.5,pt,phi);
	  hPiYvsPtQ->Fill(0.5,pt,pion.Rapidity());;
	  hPiNHitsFitvsPtQ->Fill(0.5,pt,nHitsFit);
	  hPiNHitsFitvsEtaQ->Fill(0.5,eta,nHitsFit);
	  hPiNHitsFitvsPhiQ->Fill(0.5,phi,nHitsFit);
	  hPiNHitsPossvsPtQ->Fill(0.5,pt,nHitsPoss);
	  hPiNHitsPossvsEtaQ->Fill(0.5,eta,nHitsPoss);
	  hPiNHitsPossvsPhiQ->Fill(0.5,phi,nHitsPoss);
	  hPiNHitsDedxvsPtQ->Fill(0.5,pt,nHitsDedx);
	  hPiDcavsPtQ->Fill(0.5,pt,dca);
	  hPiDedxvsPQ->Fill(0.5,p,dedx);
	  hPiDedxvsEtaQ->Fill(0.5,eta,dedx);
	  hPiDedxvsPhiQ->Fill(0.5,phi,dedx);
	  hPiNSigmaEvsPQ->Fill(0.5,p,nSigmaE);
	  hPiNSigmaPivsPQ->Fill(0.5,p,nSigmaPi);
	  hPiNSigmaPivsEtaQ->Fill(0.5,eta,nSigmaPi);
	  hPiNSigmaPivsPhiQ->Fill(0.5,phi,nSigmaPi);
	  hPiNSigmaKvsPQ->Fill(0.5,p,nSigmaK);
	  hPiNSigmaPvsPQ->Fill(0.5,p,nSigmaP);
	  }else if(charge<0){
	  hPiEtavsPtQ->Fill(-0.5,pt,eta);
	  hPiPhivsPtQ->Fill(-0.5,pt,phi);
	  hPiYvsPtQ->Fill(-0.5,pt,pion.Rapidity());;
	  hPiNHitsFitvsPtQ->Fill(-0.5,pt,nHitsFit);
	  hPiNHitsFitvsEtaQ->Fill(-0.5,eta,nHitsFit);
	  hPiNHitsFitvsPhiQ->Fill(-0.5,phi,nHitsFit);
	  hPiNHitsPossvsPtQ->Fill(-0.5,pt,nHitsPoss);
	  hPiNHitsPossvsEtaQ->Fill(-0.5,eta,nHitsPoss);
	  hPiNHitsPossvsPhiQ->Fill(-0.5,phi,nHitsPoss);
	  hPiNHitsDedxvsPtQ->Fill(-0.5,pt,nHitsDedx);
	  hPiDcavsPtQ->Fill(-0.5,pt,dca);
	  hPiDedxvsPQ->Fill(-0.5,p,dedx);
	  hPiDedxvsEtaQ->Fill(-0.5,eta,dedx);
	  hPiDedxvsPhiQ->Fill(-0.5,phi,dedx);
	  hPiNSigmaEvsPQ->Fill(-0.5,p,nSigmaE);
	  hPiNSigmaPivsPQ->Fill(-0.5,p,nSigmaPi);
	  hPiNSigmaPivsEtaQ->Fill(-0.5,eta,nSigmaPi);
	  hPiNSigmaPivsPhiQ->Fill(-0.5,phi,nSigmaPi);
	  hPiNSigmaKvsPQ->Fill(-0.5,p,nSigmaK);
	  hPiNSigmaPvsPQ->Fill(-0.5,p,nSigmaP);
	  }
	  }*/	

	//use pico track quality cuts,tof,nSigmaE to select electon&positron
	//if(dca>1.0) return kFALSE;
	if(beta2TOF<=0. || TMath::Abs(1.-1./beta2TOF)>0.025) return kFALSE;
	//if(TMath::Abs(nSigmaE)>2.0) return kFALSE;
	if (p > 0.2 && p < 1.0)
	{
		if (nSigmaE < 1.5625 * (p - 0.2) - 2 - 0.34 || nSigmaE > 2 - 0.34)
			return kFALSE;
	}
	if (p >= 1.0)
	{
		if (nSigmaE < (-1.09) || nSigmaE > 1.66)
			return kFALSE;
	}//analysis cut condition

	/*if(nHitsFit<mTpceNHitsFitCut) cout<< "nHitsFit = " << nHitsFit<<endl;
	if(nHitsDedx<mTpceNHitsDedxCut) cout<< "nHitsDedx = " << nHitsDedx <<endl;
	if(dca>2.5) cout<< "dca = " << dca <<endl;*/


	if(charge>0){
		current_ePlus[current_nEPlus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_ePlusTrkId[current_nEPlus] = i;
		current_nEPlus++;
	}else if(charge<0){
		current_eMinus[current_nEMinus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_eMinusTrkId[current_nEMinus] = i;
		current_nEMinus++;
	}

	return kTRUE;
}
//____________________________________________________________
void photonE(EVENT *event)
{
	for(Int_t i=0;i<current_nEPlus;i++) current_phePlusTag[i] = kFALSE;
	for(Int_t i=0;i<current_nEMinus;i++) current_pheMinusTag[i] = kFALSE;

	TLorentzVector pair(0,0,0,0);

	for(Int_t i=0;i<current_nEPlus;i++){
		for(Int_t j=0;j<current_nEMinus;j++){
			pair = current_ePlus[i] + current_eMinus[j];
			hULMvsPt->Fill(pair.Pt(),pair.M());
			if(pair.M()>0. && pair.M()<mPheMassCutWTof){
				if (TMath::Abs(pair.Rapidity())<=mPairYCut)
				{
					current_phePlusTag[i] = kTRUE;
					current_pheMinusTag[j] = kTRUE;
					hPhotonEtavsPhi->Fill(pair.Eta(),pair.Phi());
				}
			}
		}
	}

	for(Int_t i=0;i<current_nEPlus;i++){
		for(Int_t j=i+1;j<current_nEPlus;j++){
			pair = current_ePlus[i] + current_ePlus[j];
			hLPosMvsPt->Fill(pair.Pt(),pair.M());
			if(pair.M()>0. && pair.M()<mPheMassCutWTof){
				if (TMath::Abs(pair.Rapidity())<=mPairYCut)
				{
					current_phePlusTag[i] = kFALSE;
					current_phePlusTag[j] = kFALSE;
				}
			}
		}
	}

	for(Int_t i=0;i<current_nEMinus;i++){
		for(Int_t j=i+1;j<current_nEMinus;j++){
			pair = current_eMinus[i] + current_eMinus[j];
			hLNegMvsPt->Fill(pair.Pt(),pair.M());
			if(pair.M()>0. && pair.M()<mPheMassCutWTof){
				if (TMath::Abs(pair.Rapidity())<=mPairYCut)
				{
					current_pheMinusTag[i] = kFALSE;
					current_pheMinusTag[j] = kFALSE;
				}
			}
		}
	}

	for(Int_t i=0;i<current_nEPlus;i++){
		if(!current_phePlusTag[i]) continue;

		Float_t pt = current_ePlus[i].Pt();
		Float_t eta = current_ePlus[i].Eta();
		Float_t phi = current_ePlus[i].Phi();
		Float_t p = current_ePlus[i].P();
		Float_t Y = current_ePlus[i].Rapidity();
		//cout<<"Eta:"<<eta<<"  Y:"<<Y<<endl;
		Int_t trkId = current_ePlusTrkId[i];
		Int_t nHitsFit = event->mNHitsFit[trkId];
		Int_t nHitsPoss = event->mNHitsPoss[trkId];
		Int_t nHitsDedx = event->mNHitsDedx[trkId];
		Float_t dca = event->mDca[trkId];
		Float_t dEdx = event->mDedx[trkId];
		Float_t nSigmaE = event->mNSigmaE[trkId];
		Float_t originZ = event->mgOriginZ[trkId];
		Float_t originPhi = event->mgPhi[trkId];
		/*		Float_t nSigmaPi = event->mNSigmaPi[trkId]/1000.;
				Float_t nSigmaK = event->mNSigmaK[trkId]/1000.;
				Float_t nSigmaP = event->mNSigmaP[trkId]/1000.;
				*/
		
		if(nHitsFit<mTpceNHitsFitCut) cout<< "nHitsFit = " << nHitsFit<<endl;
		if(nHitsDedx<mTpceNHitsDedxCut) cout<< "nHitsDedx = " << nHitsDedx <<endl;
		if(dca>2.5) cout<< "dca = " << dca <<endl;

		hEtavsPtQ->Fill(0.5,pt,eta);
		hPhivsPtQ->Fill(0.5,pt,phi);
		hYvsPtQ->Fill(0.5,pt,Y);
		hNHitsFitvsPtQ->Fill(0.5,pt,nHitsFit);
		hNHitsFitvsEtaQ->Fill(0.5,eta,nHitsFit);
		hNHitsFitvsPhiQ->Fill(0.5,phi,nHitsFit);
		hNHitsPossvsPtQ->Fill(0.5,pt,nHitsPoss);
		hNHitsPossvsEtaQ->Fill(0.5,eta,nHitsPoss);
		hNHitsPossvsPhiQ->Fill(0.5,phi,nHitsPoss);
		hNHitsDedxvsPtQ->Fill(0.5,pt,nHitsDedx);
		hDcavsPtQ->Fill(0.5,pt,dca);
		hDedxvsPQ->Fill(0.5,p,dEdx);
		hDedxvsEtaQ->Fill(0.5,eta,dEdx);
		hDedxvsPhiQ->Fill(0.5,phi,dEdx);
		hNSigmaEvsPQ->Fill(0.5,p,nSigmaE);
		hNSigmaEvsEtaQ->Fill(0.5,eta,nSigmaE);
		hNSigmaEvsPhiQ->Fill(0.5,phi,nSigmaE);
		//		hNSigmaPivsPQ->Fill(0.5,p,nSigmaPi);
		//		hNSigmaKvsPQ->Fill(0.5,p,nSigmaK);
		//		hNSigmaPvsPQ->Fill(0.5,p,nSigmaP);
		hPENHitsFitsvsPtvsEtaPlus->Fill(pt,eta,nHitsFit);
		hPENHitsFitsvsPtvsPhiPlus->Fill(pt,phi,nHitsFit);
		hPENHitsDedxvsPtvsEtaPlus->Fill(pt,eta,nHitsDedx);
		hPENHitsDedxvsPtvsPhiPlus->Fill(pt,phi,nHitsDedx);
		hPEDcavsPtvsEtaPlus->Fill(pt,eta,dca);
		hPEDcavsPtvsPhiPlus->Fill(pt,phi,dca);
		//Zhen add it
		hPlusOriginZ->Fill(originZ);
		hPlusOriginPhi->Fill(originPhi);
		hPlusOriginZvsPhi->Fill(originZ,originPhi);

	}

	for(Int_t i=0;i<current_nEMinus;i++){
		if(!current_pheMinusTag[i]) continue;

		Float_t pt = current_eMinus[i].Pt();
		Float_t eta = current_eMinus[i].Eta();
		Float_t phi = current_eMinus[i].Phi();
		Float_t p = current_eMinus[i].P();
		Float_t Y = current_eMinus[i].Rapidity();
		Int_t trkId = current_eMinusTrkId[i];
		Int_t nHitsFit = event->mNHitsFit[trkId];
		Int_t nHitsPoss = event->mNHitsPoss[trkId];
		Int_t nHitsDedx = event->mNHitsDedx[trkId];
		Float_t dca = event->mDca[trkId];
		Float_t dEdx = event->mDedx[trkId];
		Float_t nSigmaE = event->mNSigmaE[trkId];
		Float_t originZ = event->mgOriginZ[trkId];
		Float_t originPhi = event->mgPhi[trkId];
		//		Float_t nSigmaPi = event->mNSigmaPi[trkId]/1000.;
		//		Float_t nSigmaK = event->mNSigmaK[trkId]/1000.;
		//		Float_t nSigmaP = event->mNSigmaP[trkId]/1000.;
		if(nHitsFit<mTpceNHitsFitCut) cout<< "nHitsFit = " << nHitsFit<<endl;
		if(nHitsDedx<mTpceNHitsDedxCut) cout<< "nHitsDedx = " << nHitsDedx <<endl;
		if(dca>2.5) cout<< "dca = " << dca <<endl;

		hEtavsPtQ->Fill(-0.5,pt,eta);
		hPhivsPtQ->Fill(-0.5,pt,phi);
		hYvsPtQ->Fill(-0.5,pt,Y);
		hNHitsFitvsPtQ->Fill(-0.5,pt,nHitsFit);
		hNHitsFitvsEtaQ->Fill(-0.5,eta,nHitsFit);
		hNHitsFitvsPhiQ->Fill(-0.5,phi,nHitsFit);
		hNHitsPossvsPtQ->Fill(-0.5,pt,nHitsPoss);
		hNHitsPossvsEtaQ->Fill(-0.5,eta,nHitsPoss);
		hNHitsPossvsPhiQ->Fill(-0.5,phi,nHitsPoss);
		hNHitsDedxvsPtQ->Fill(-0.5,pt,nHitsDedx);
		hDcavsPtQ->Fill(-0.5,pt,dca);
		hDedxvsPQ->Fill(-0.5,p,dEdx);
		hDedxvsEtaQ->Fill(-0.5,eta,dEdx);
		hDedxvsPhiQ->Fill(-0.5,phi,dEdx);
		hNSigmaEvsPQ->Fill(-0.5,p,nSigmaE);
		hNSigmaEvsEtaQ->Fill(-0.5,eta,nSigmaE);
		hNSigmaEvsPhiQ->Fill(-0.5,phi,nSigmaE);
		hPENHitsFitsvsPtvsEtaMinus->Fill(pt,eta,nHitsFit);
		hPENHitsFitsvsPtvsPhiMinus->Fill(pt,phi,nHitsFit);
		hPENHitsDedxvsPtvsEtaMinus->Fill(pt,eta,nHitsDedx);
		hPENHitsDedxvsPtvsPhiMinus->Fill(pt,phi,nHitsDedx);
		hPEDcavsPtvsEtaMinus->Fill(pt,eta,dca);
		hPEDcavsPtvsPhiMinus->Fill(pt,phi,dca);
		//zhen add it
		hMinusOriginZ->Fill(originZ);
		hMinusOriginPhi->Fill(originPhi);
		hMinusOriginZvsPhi->Fill(originZ,originPhi);
		//		hNSigmaPivsPQ->Fill(-0.5,p,nSigmaPi);
		//		hNSigmaKvsPQ->Fill(-0.5,p,nSigmaK);
		//		hNSigmaPvsPQ->Fill(-0.5,p,nSigmaP);
	}

}
//____________________________________________________________
void bookHistograms()
{
	hVyvsVx = new TH2F("hVyvsVx","hVyvsVx;Vx (cm); Vy (cm)",150,-1.5,1.5,150,-1.5,1.5);
	hVyvsVz = new TH2F("hVyvsVz","hVyvsVz;Vz (cm); Vy (cm)",150,-1.5,1.5,150,-1.5,1.5);
	hRefMultvsRefMultCorr = new TH2F("hRefMultvsRefMultCorr","hRefMultvsRefMultCorr; refMultCorr; refMult",1000,0,1000,1000,0,1000);

	hULMvsPt = new TH2F("hULMvsPt","hULMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})",100,0,10,3000,0,0.3);
	hLPosMvsPt = new TH2F("hLPosMvsPt","hLPosMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})",100,0,10,3000,0,0.3);
	hLNegMvsPt = new TH2F("hLNegMvsPt","hLNegMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})",100,0,10,3000,0,0.3);

	hEtavsPtQ = new TH3F("hEtavsPtQ","hEtavsPtQ;q;p_{T} (GeV/c);#eta",2,-1,1,500,0,10,200,-2,2);
	hPhivsPtQ = new TH3F("hPhivsPtQ","hPhivsPtQ;q;p_{T} (GeV/c);#phi",2,-1,1,500,0,10,360,-TMath::Pi(),TMath::Pi());
	hYvsPtQ = new TH3F("hYvsPtQ","hYvsPtQ;q;p_{T} (GeV/c);Rapidity",2,-1,1,500,0,10,200,-2,2);
	hNHitsFitvsPtQ = new TH3F("hNHitsFitvsPtQ","hNHitsFitvsPtQ;q;p_{T} (GeV/c);nHitsFit",2,-1,1,500,0,10,50,0,50);
	hNHitsFitvsEtaQ = new TH3F("hNHitsFitvsEtaQ","hNHitsFitvsEtaQ;q;#eta;nHitsFit",2,-1,1,200,-2,2,50,0,50);
	hNHitsFitvsPhiQ = new TH3F("hNHitsFitvsPhiQ","hNHitsFitvsPhiQ;q;#phi;nHitsFit",2,-1,1,360,-TMath::Pi(),TMath::Pi(),50,0,50);
	hNHitsPossvsPtQ = new TH3F("hNHitsPossvsPtQ","hNHitsPossvsPtQ;q;p_{T} (GeV/c);nHitsPoss",2,-1,1,500,0,10,50,0,50);
	hNHitsPossvsEtaQ = new TH3F("hNHitsPossvsEtaQ","hNHitsPossvsEtaQ;q;#eta;nHitsPoss",2,-1,1,200,-2,2,50,0,50);
	hNHitsPossvsPhiQ = new TH3F("hNHitsPossvsPhiQ","hNHitsPossvsPhiQ;q;#phi;nHitsPoss",2,-1,1,360,-TMath::Pi(),TMath::Pi(),50,0,50);
	hNHitsDedxvsPtQ = new TH3F("hNHitsDedxvsPtQ","hNHitsDedxvsPtQ;q;p_{T} (GeV/c);nHitsDedx",2,-1,1,500,0,10,50,0,50);
	hDcavsPtQ = new TH3F("hDcavsPtQ","hDcavsPtQ;q;p_{T} (GeV/c);dca (cm)",2,-1,1,500,0,10,2000,-1.e-6,4-1.e-6);
	hDedxvsPQ = new TH3F("hDedxvsPQ","hDedxvsPQ;q;p (GeV/c);dE/dx (KeV/cm)",2,-1,1,500,0,10,300,-1.e-6,15-1.e-6);
	hDedxvsEtaQ = new TH3F("hDedxvsEtaQ","hDedxvsEtaQ;q;#eta;dE/dx (KeV/cm)",2,-1,1,200,-2,2,300,-1.e-6,15-1.e-6);
	hDedxvsPhiQ = new TH3F("hDedxvsPhiQ","hDedxvsPhiQ;q;#phi;dE/dx (KeV/cm)",2,-1,1,360,-TMath::Pi(),TMath::Pi(),300,-1.e-6,15-1.e-6);
	hNSigmaEvsPQ = new TH3F("hNSigmaEvsPQ","hNSigmaEvsPQ;q;p (GeV/c);n#sigma_{e}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hNSigmaEvsEtaQ = new TH3F("hNSigmaEvsEtaQ","hNSigmaEvsEtaQ;q;#eta;n#sigma_{e}",2,-1,1,200,-2,2,600,-15-1.e-6,15-1.e-6);
	hNSigmaEvsPhiQ = new TH3F("hNSigmaEvsPhiQ","hNSigmaEvsPhiQ;q;#phi;n#sigma_{e}",2,-1,1,360,-TMath::Pi(),TMath::Pi(),600,-15-1.e-6,15-1.e-6);
	hNSigmaPivsPQ = new TH3F("hNSigmaPivsPQ","hNSigmaPivsPQ;q;p (GeV/c);n#sigma_{#pi}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hNSigmaKvsPQ = new TH3F("hNSigmaKvsPQ","hNSigmaKvsPQ;q;p (GeV/c);n#sigma_{k}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hNSigmaPvsPQ = new TH3F("hNSigmaPvsPQ","hNSigmaPvsPQ;q;p (GeV/c);n#sigma_{p}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);

	hPiEtavsPtQ = new TH3F("hPiEtavsPtQ","hPiEtavsPtQ;q;p_{T} (GeV/c);#eta",2,-1,1,500,0,10,200,-2,2);
	hPiPhivsPtQ = new TH3F("hPiPhivsPtQ","hPiPhivsPtQ;q;p_{T} (GeV/c);#phi",2,-1,1,500,0,10,360,-TMath::Pi(),TMath::Pi());
	hPiYvsPtQ = new TH3F("hPiYvsPtQ","hPiYvsPtQ;q;p_{T} (GeV/c);Rapidity",2,-1,1,500,0,10,200,-2,2);
	hPiNHitsFitvsPtQ = new TH3F("hPiNHitsFitvsPtQ","hPiNHitsFitvsPtQ;q;p_{T} (GeV/c);nHitsFit",2,-1,1,500,0,10,50,0,50);
	hPiNHitsFitvsEtaQ = new TH3F("hPiNHitsFitvsEtaQ","hPiNHitsFitvsEtaQ;q;#eta;nHitsFit",2,-1,1,200,-2,2,50,0,50);
	hPiNHitsFitvsPhiQ = new TH3F("hPiNHitsFitvsPhiQ","hPiNHitsFitvsPhiQ;q;#phi;nHitsFit",2,-1,1,360,-TMath::Pi(),TMath::Pi(),50,0,50);
	hPiNHitsPossvsPtQ = new TH3F("hPiNHitsPossvsPtQ","hPiNHitsPossvsPtQ;q;p_{T} (GeV/c);nHitsPoss",2,-1,1,500,0,10,50,0,50);
	hPiNHitsPossvsEtaQ = new TH3F("hPiNHitsPossvsEtaQ","hPiNHitsPossvsEtaQ;q;#eta;nHitsPoss",2,-1,1,200,-2,2,50,0,50);
	hPiNHitsPossvsPhiQ = new TH3F("hPiNHitsPossvsPhiQ","hPiNHitsPossvsPhiQ;q;#phi;nHitsPoss",2,-1,1,360,-TMath::Pi(),TMath::Pi(),50,0,50);
	hPiNHitsDedxvsPtQ = new TH3F("hPiNHitsDedxvsPtQ","hPiNHitsDedxvsPtQ;q;p_{T} (GeV/c);nHitsDedx",2,-1,1,500,0,10,50,0,50);
	hPiDcavsPtQ = new TH3F("hPiDcavsPtQ","hPiDcavsPtQ;q;p_{T} (GeV/c);dca (cm)",2,-1,1,500,0,10,2000,-1.e-6,4-1.e-6);
	hPiDedxvsPQ = new TH3F("hPiDedxvsPQ","hPiDedxvsPQ;q;p (GeV/c);dE/dx (KeV/cm)",2,-1,1,500,0,10,300,-1.e-6,15-1.e-6);
	hPiDedxvsEtaQ = new TH3F("hPiDedxvsEtaQ","hPiDedxvsEtaQ;q;#eta;dE/dx (KeV/cm)",2,-1,1,200,-2,2,300,-1.e-6,15-1.e-6);
	hPiDedxvsPhiQ = new TH3F("hPiDedxvsPhiQ","hPiDedxvsPhiQ;q;#phi;dE/dx (KeV/cm)",2,-1,1,360,-TMath::Pi(),TMath::Pi(),300,-1.e-6,15-1.e-6);
	hPiNSigmaEvsPQ = new TH3F("hPiNSigmaEvsPQ","hPiNSigmaEvsPQ;q;p (GeV/c);n#sigma_{e}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hPiNSigmaPivsPQ = new TH3F("hPiNSigmaPivsPQ","hPiNSigmaPivsPQ;q;p (GeV/c);n#sigma_{#pi}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hPiNSigmaPivsEtaQ = new TH3F("hPiNSigmaPivsEtaQ","hPiNSigmaPivsEtaQ;q;#eta;n#sigma_{e}",2,-1,1,200,-2,2,600,-15-1.e-6,15-1.e-6);
	hPiNSigmaPivsPhiQ = new TH3F("hPiNSigmaPivsPhiQ","hPiNSigmaPivsPhiQ;q;#phi;n#sigma_{e}",2,-1,1,360,-TMath::Pi(),TMath::Pi(),600,-15-1.e-6,15-1.e-6);
	hPiNSigmaKvsPQ = new TH3F("hPiNSigmaKvsPQ","hPiNSigmaKvsPQ;q;p (GeV/c);n#sigma_{k}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hPiNSigmaPvsPQ = new TH3F("hPiNSigmaPvsPQ","hPiNSigmaPvsPQ;q;p (GeV/c);n#sigma_{p}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);

	// wangzhen add it for 3D check
	hPENHitsFitsvsPtvsEtaMinus = new TH3F("hPENitsFitsvsPtvsEtaMinus","hPENitsFitsvsPtvsEtaMinus;p_{T};#eta;nHitsFit",500,0,10,200,-2,2,50,0,50);
	hPENHitsFitsvsPtvsEtaPlus = new TH3F("hPENitsFitsvsPtvsEtaPlus","hPENitsFitsvsPtvsEtaPlus;p_{T};#eta;nHitsFit",500,0,10,200,-2,2,50,0,50);
	hPENHitsFitsvsPtvsPhiMinus = new TH3F("hPENitsFitsvsPtvsPhiMinus","hPENitsFitsvsPtvsPhiMinus;p_{T};#phi;nHitsFit",500,0,10,360,-TMath::Pi(),TMath::Pi(),50,0,50);
	hPENHitsFitsvsPtvsPhiPlus = new TH3F("hPENitsFitsvsPtvsPhiPlus","hPENitsFitsvsPtvsPhiPlus;p_{T};#phi;nHitsFit",500,0,10,360,-TMath::Pi(),TMath::Pi(),50,0,50);
	hPENHitsDedxvsPtvsEtaMinus = new TH3F("hPENHitsDedxvsPtvsEtaMinus","hPENHitsDedxvsPtvsEtaMinus;p_{T};#eta;nHitsDedx",500,0,10,200,-2,2,50,0,50);
    hPENHitsDedxvsPtvsEtaPlus = new TH3F("hPENHitsDedxvsPtvsEtaPlus","hPENHitsDedxvsPtvsEtaPlus;p_{T};#eta;nHitsDedx",500,0,10,200,-2,2,50,0,50);
    hPENHitsDedxvsPtvsPhiMinus = new TH3F("hPENHitsDedxvsPtvsPhiMinus","hPENHitsDedxvsPtvsPhiMinus;p_{T};#phi;nHitsDedx",500,0,10,360,-TMath::Pi(),TMath::Pi(),50,0,50);
    hPENHitsDedxvsPtvsPhiPlus = new TH3F("hPENHitsDedxvsPtvsPhiPlus","hPENHitsDedxvsPtvsPhiPlus;p_{T};#phi;nHitsDedx",500,0,10,360,-TMath::Pi(),TMath::Pi(),50,0,50);
	hPEDcavsPtvsEtaMinus = new TH3F("hPEDcavsPtvsEtaMinus","hPEDcavsPtvsEtaMinus;p_{T};#eta;Dca(cm)",500,0,10,360,-2,2,400,0,10);
    hPEDcavsPtvsEtaPlus = new TH3F("hPEDcavsPtvsEtaPlus","hPEDcavsPtvsEtaPlus;p_{T};#eta;Dca(cm)",500,0,10,360,-2,2,400,0,10);
    hPEDcavsPtvsPhiMinus = new TH3F("hPEDcavsPtvsPhiMinus","hPEDcavsPtvsPhiMinus;p_{T};#phi;Dca(cm)",500,0,10,360,-TMath::Pi(),TMath::Pi(),400,0,10);
    hPEDcavsPtvsPhiPlus = new TH3F("hPEDcavsPtvsPhiPlus","hPEDcavsPtvsPhiPlus;p_{T};#phi;Dca(cm)",500,0,10,360,-TMath::Pi(),TMath::Pi(),400,0,10);

	// to chenk origin Z and Phi distribution
	hMinusOriginZ = new TH1D("hMinusOriginZ","hMinusOriginZ",600,-30,30);
	hMinusOriginPhi = new TH1D("hMinusOriginPhi","hMinusOriginPhi",360,-TMath::Pi(),TMath::Pi());
	hMinusOriginZvsPhi = new TH2D("hMinusOriginZvsPhi","hMinusOriginZvsPhi",600,-30,30,360,-TMath::Pi(),TMath::Pi());
	hPlusOriginZ = new TH1D("hPlusOriginZ","hPlusOriginZ",600,-30,30);
	hPlusOriginPhi = new TH1D("hPlusOriginPhi","PlushOriginPhi",360,-TMath::Pi(),TMath::Pi());
	hPlusOriginZvsPhi = new TH2D("hPlusOriginZvsPhi","hPlusOriginZvsPhi",600,-30,30,360,-TMath::Pi(),TMath::Pi());
	hPhotonEtavsPhi = new TH2D("hPhotonEtavsPhi","hPhotonEtavsPhi",480,-1.2,1.2,360,-TMath::Pi(),TMath::Pi());

}
//____________________________________________________________
void writeHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.QAhisto.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	TFile *mFile = new TFile(buf,"recreate");
	mFile->cd();

	hVyvsVx->Write();
	hVyvsVz->Write();
	hRefMultvsRefMultCorr->Write();

	hULMvsPt->Write();
	hLPosMvsPt->Write();
	hLNegMvsPt->Write();

	hEtavsPtQ->Write();
	hPhivsPtQ->Write();
	hYvsPtQ->Write();
	hNHitsFitvsPtQ->Write();
	hNHitsFitvsEtaQ->Write();
	hNHitsFitvsPhiQ->Write();
	hNHitsPossvsPtQ->Write();
	hNHitsPossvsEtaQ->Write();
	hNHitsPossvsPhiQ->Write();
	hNHitsDedxvsPtQ->Write();
	hDcavsPtQ->Write();
	hDedxvsPQ->Write();
	hDedxvsEtaQ->Write();
	hDedxvsPhiQ->Write();
	hNSigmaEvsPQ->Write();
	hNSigmaEvsEtaQ->Write();
	hNSigmaEvsPhiQ->Write();
	hNSigmaPivsPQ->Write();
	hNSigmaKvsPQ->Write();
	hNSigmaPvsPQ->Write();

	hPiEtavsPtQ->Write();
	hPiPhivsPtQ->Write();
	hPiYvsPtQ->Write();
	hPiNHitsFitvsPtQ->Write();
	hPiNHitsFitvsEtaQ->Write();
	hPiNHitsFitvsPhiQ->Write();
	hPiNHitsPossvsPtQ->Write();
	hPiNHitsPossvsEtaQ->Write();
	hPiNHitsPossvsPhiQ->Write();
	hPiNHitsDedxvsPtQ->Write();
	hPiDcavsPtQ->Write();
	hPiDedxvsPQ->Write();
	hPiDedxvsEtaQ->Write();
	hPiDedxvsPhiQ->Write();
	hPiNSigmaEvsPQ->Write();
	hPiNSigmaPivsPQ->Write();
	hPiNSigmaPivsEtaQ->Write();
	hPiNSigmaPivsPhiQ->Write();
	hPiNSigmaKvsPQ->Write();
	hPiNSigmaPvsPQ->Write();

	hPENHitsFitsvsPtvsEtaMinus->Write();
	hPENHitsFitsvsPtvsEtaPlus ->Write();
	hPENHitsFitsvsPtvsPhiMinus->Write();
	hPENHitsFitsvsPtvsPhiPlus ->Write();
	hPENHitsDedxvsPtvsEtaMinus->Write();
	hPENHitsDedxvsPtvsPhiMinus->Write();
	hPENHitsDedxvsPtvsEtaPlus->Write();
	hPENHitsDedxvsPtvsPhiPlus->Write();
	hPEDcavsPtvsEtaMinus->Write();
	hPEDcavsPtvsPhiMinus->Write();
	hPEDcavsPtvsEtaPlus->Write();
	hPEDcavsPtvsPhiPlus->Write();

	hMinusOriginZ->Write();
	hMinusOriginPhi->Write();
	hMinusOriginZvsPhi->Write();
	hPlusOriginZ->Write();
	hPlusOriginPhi->Write();
	hPlusOriginZvsPhi->Write();
	hPhotonEtavsPhi->Write();

}
//____________________________________________________________
bool Init()
{
	cout<<endl;

	ifstream indata;

	indata.open("mBadRunList.dat");
	mBadRunId.clear();
	if(indata.is_open()){
		cout<<"read in the bad run number list ...";
		Int_t id;
		while(indata>>id){
			mBadRunId[id] = id;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the bad run number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();
	for(map<Int_t,Int_t>::iterator iter=mBadRunId.begin();iter!=mBadRunId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;
	cout<<endl;

	indata.open("/star/u/syang/run12/uu/minibias/badTPCSectorGeo/fitPars.dat");
	if(!indata.is_open()){                                                      
		cout<<"Failed to load the bad TPC sector geometry parameters !!!"<<endl;
		return kFALSE;
	}else{
		cout<<"Load bad TPC sector geometry parameters ...";
		Int_t idx = 0;
		Double_t tmp0,tmp1;
		while(indata>>tmp0>>tmp1){
			par[idx/4][idx%4] = tmp0;
			parErr[idx/4][idx%4] = tmp1;
			idx++;
		}
		cout<<" [OK]"<<endl;
	}
	indata.close();
	for(Int_t i=0;i<4;i++){
		for(Int_t j=0;j<4;j++){
			cout<<par[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	funPosHi = new TF1("funPosHi","[0]*exp(-([1]/x)**[2])+[3]",0.1,60);
	funPosHi->SetParameters(par[0][0],par[0][1],par[0][2],par[0][3]);
	funPosLow = new TF1("funPosLow","[0]*exp(-([1]/x)**[2])+[3]",0.1,60);
	funPosLow->SetParameters(par[1][0],par[1][1],par[1][2],par[1][3]);
	funNegHi = new TF1("funNegHi","[0]*exp(-([1]/x)**[2])+[3]",-60,-0.1);
	funNegHi->SetParameters(par[2][0],par[2][1],par[2][2],par[2][3]);
	funNegLow = new TF1("funNegLow","[0]*exp(-([1]/x)**[2])+[3]",-60,-0.1);
	funNegLow->SetParameters(par[3][0],par[3][1],par[3][2],par[3][3]);

	cout<<"Initialization DONE !!!"<<endl;

	return kTRUE;
}
