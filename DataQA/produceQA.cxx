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

#include "miniDst.h"
#include "StRefMultCorr.h"
#include "cuts.h"

using namespace std;
#endif


StRefMultCorr *refMultCorrUtil;

Int_t runIndex;
Int_t randomId;
map<Int_t,Int_t> mTotalRunId;

bool Init();
void bookHistograms();
bool passEvent(miniDst* event);
bool passTrack(miniDst* event, Int_t i);
void writeHistograms(char* outFile);

//define histograms
//inclusive QA
TH1D *hnEvts;
TH2F *hVyvsVx;
TH2F *hVyvsVz;
TH1F *hVr;
TH2F *hVPDVzvsTPCVz;
TH2F *hDeltaZvsTPCVz;
TH2F *hDeltaZvsVPDVz;
TH2F *hTPCVzvsRefMult;
TH2F *hVPDVzvsRefMult;
TH2F *hDeltaZvsRefMult;
TH2F *hZDCXvsRefMult;
TH2F *hBBCXvsRefMult;
TH1F *hNHitsFit;
TH1F *hNHitsPoss;
TH1F *hNHitsDedx;
TH1F *hNHitsFitoverNHitsPoss;//new
TH2F *hDcavsPt;
TH2F *hBetavsP;
TH2F *hBetavsEta;
TH2F *hBetavsPhi;
TH2F *hDedxvsP;
TH2F *hDedxvsEta;
TH2F *hDedxvsPhi;
TH2F *hNSigmaEvsP_woTOF;
TH2F *hNSigmaEvsP_wTOF;
TH2F *hNSigmaEvsP;
TH2F *hNSigmaEvsEta;
TH2F *hNSigmaEvsPhi;
TH2F *hMSquarevsP;
TH2F *hEtavsPhi;
TH2F *hEEtavsPhi;
TH2F *hEtavsPt;
TH2F *hEEtavsPt;
TH2F *hPhivsPt;
TH2F *hEPhivsPt;
TH2F *hEVxvsVy;
TH2F *hEVyvsVz;
TH2F *hEVxvsVz;

TH2D* hnTofMatchvsRefMult;
TH2D* hnTofMatchvsRefMult_Rejection;

//run by run QA
TProfile *hNHitsFitoverNHitsPossvsRunIndex;//new
TProfile *hVrvsRunIndex;//new
TProfile *hNHitsFitvsRunIndex;
TProfile *hZDCXvsRunIndex;
TProfile *hBBCXvsRunIndex;
TProfile *hTPCVzvsRunIndex;
TProfile *hVPDVzvsRunIndex;
TProfile *hDeltaZvsRunIndex;
TProfile *hRefMultvsRunIndex;
TProfile *hRefMultCorrvsRunIndex;
TProfile *hPtvsRunIndex;
TProfile *hgPtvsRunIndex;//new
TProfile *hEtavsRunIndex;
TProfile *hPhivsRunIndex;
TProfile *hDcavsRunIndex;
TProfile *hBetavsRunIndex;
TProfile *hEBetaDiffvsRunIndex;
TProfile *hKBetaDiffvsRunIndex;
TProfile *hDedxvsRunIndex;
TProfile *hNSigmaEvsRunIndex;
TProfile* hNHitsdEdxvsRunIndex;
TProfile* hnTrkvsRunIndex;

TF1* Pileuplimit;

int mTotalRun;
int mTotalDay;

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
        //ofstream out;
        //out.open("./outlog.txt");
	//intialization
	if( !Init() ){
		cout<<"The initialization is failed !!!"<<endl;
		return 0;
	}
	bookHistograms();
	//refMultCorrUtil = new StRefMultCorr("refmult");

	//+-------------+
	//+-------------+
	//| loop events |
	//+-------------+
	miniDst *event = new miniDst(chain);
	// Int_t nEvts  = 20;
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;
	int i_Trk = 0;
	for(int i=0;i<nEvts;i++){

		//initialize the struct

		// cout << "begin " << i << "th entry...." << endl;
		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;
		event->GetEntry(i);
		// event->Show(i);

		Int_t runId = event->mRunId;
		// cout << "runId = " << runId << endl;
		map<Int_t,Int_t>::iterator iter = mTotalRunId.find(runId);
		if(iter != mTotalRunId.end())
			runIndex = iter->second;
		else{
			runIndex = -1;
			cout << "Run numner = " <<iter->second<< endl;
			cout<<"Can not find the runNumber in the runNumber list"<<endl;
		}

		if(runIndex<0) continue;

		// cout << "passEvt = " << runId << endl;
		if(!passEvent(event)) continue; 
		// hnEvts->Fill(1);

		Int_t nTrks = event->mNTrks;
		if(nTrks>0 && i_Trk == 0) 
		{
			i_Trk++;
			event->Show(i);
		}
		// cout << "end passEvt = " << i << "nTrks = " << nTrks << endl;
		for(int j=0;j<nTrks;j++) 
		{
			passTrack(event,j);
		}
	}

	writeHistograms(outFile);
	delete chain;
        //out.close();
	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
bool passEvent(miniDst* event)
{
	bool eventflag=kFALSE;
	for(Int_t i=0;i<event->mNTrigs;i++){
		// Trigger ID for AuAu @ 9.2 GeV
		if(event->mTrigId[i] == 780010) eventflag = kTRUE; 
		if(event->mTrigId[i] == 780020) eventflag = kTRUE; 
	}
	if(!eventflag) return kFALSE;
	hnEvts->Fill(1);// minbias trigger events

	Int_t runId = event->mRunId;
	Int_t zdcRate =event->mZDCRate;
	Int_t bbcRate =event->mBBCRate;
	Int_t refMult = event->mRefMult;
	Int_t mnTOFMatch = event->mnTOFMatch;
	Float_t vx = event->mVertexX;
	Float_t vy = event->mVertexY;
	Float_t vz = event->mVertexZ;
	Float_t vpdVz = event->mVpdVz;
	Float_t vzDiff = vz - vpdVz;
	Int_t nTrks = event->mNTrks;
        Float_t vr = sqrt(pow(vx,2)+pow(vy,2));

	//already have event cuts in pico production

	hVyvsVx->Fill(vx,vy);
	hVyvsVz->Fill(vz,vy);
        hVr->Fill(vr);
	hVPDVzvsTPCVz->Fill(vz,vpdVz);
	hDeltaZvsTPCVz->Fill(vz,vzDiff);
	hDeltaZvsVPDVz->Fill(vpdVz,vzDiff);
	hTPCVzvsRefMult->Fill(refMult,vz);
	hVPDVzvsRefMult->Fill(refMult,vpdVz);
	hDeltaZvsRefMult->Fill(refMult,vzDiff);
	hZDCXvsRefMult->Fill(refMult,zdcRate/1000.);
	hBBCXvsRefMult->Fill(refMult,bbcRate/1000.);

	hZDCXvsRunIndex->Fill(runIndex,zdcRate/1000.);
	hBBCXvsRunIndex->Fill(runIndex,bbcRate/1000.);
	hTPCVzvsRunIndex->Fill(runIndex,vz);
        //vr hist
        hVrvsRunIndex->Fill(runIndex,vr);//new
	hVPDVzvsRunIndex->Fill(runIndex,vpdVz);
	hDeltaZvsRunIndex->Fill(runIndex,vzDiff);
	hRefMultvsRunIndex->Fill(runIndex,refMult);
	hnTofMatchvsRefMult->Fill(refMult,mnTOFMatch);
	hnTrkvsRunIndex->Fill(runIndex,nTrks);


	if (mnTOFMatch > Pileuplimit->Eval(refMult))
	{
		hnEvts->Fill(2); // with pile up rejection
		hnTofMatchvsRefMult_Rejection->Fill(refMult,mnTOFMatch);
		if (abs(vz) < 35)
		{
			hnEvts->Fill(3); // |Vz| < 35cm
			if(vr < 2)
			{
				hnEvts->Fill(4); // Vr < 2cm
			}
		}

	}
//	if(TMath::Abs(vz)<mVzCut){//make sure vz is in the range listed in the parameters file
//		refMultCorrUtil->init(runId);
//		refMultCorrUtil->initEvent(refMult,vz,zdcRate);
//		Double_t refMultCor = refMultCorrUtil->getRefMultCorr();//if you want to call the getRefMultCorr() with no argument, it must be called after initEvent()
//		hRefMultCorrvsRunIndex->Fill(runIndex,refMultCor);
//	}
	return kTRUE;
}
//______________________________________________________________
bool passTrack(miniDst* event, Int_t i)
{
	//debug 
	// cout << "passing "<< i << "th track" << endl;
	//debug
	Int_t charge = event->mCharge[i];
	Int_t nHitsFit = event->mNHitsFit[i];
	Int_t nHitsDedx = event->mNHitsDedx[i];
	Int_t nHitsPoss = event->mNHitsPoss[i];
        Float_t FitoverPoss = float (nHitsFit)/float (nHitsPoss);//new 
	Float_t dedx = event->mDedx[i];
	Float_t nSigmaE = event->mNSigmaE[i];
	Float_t dca = event->mDca[i];
	Float_t pt = event->mPt[i];
	Float_t eta = event->mEta[i];
	Float_t phi = event->mPhi[i];
	Float_t beta2TOF = event->mBeta2TOF[i];
	Float_t mgOriginX = event->mgOriginX[i]; 
	Float_t mgOriginY = event->mgOriginY[i];
	Float_t mgOriginZ = event->mgOriginZ[i];
	TVector3 mom;
	mom.SetPtEtaPhi(pt,eta,phi);
	Float_t p = mom.Mag();
	//for reading data test 
	// cout << "in the tree pT = " << event->mPt[i] << " eta = " << event->mEta[i] << " phi = " << event->mPhi[i] << endl;  
	// cout << "pT = " << pt << " eta = " << eta << " phi = " << phi << endl;
	//already applied some track quality cuts in pico production


	hNHitsFit->Fill(charge*nHitsFit);
        //ofstream out;
        //out.open("./hitsfitlog.txt");
        //cout<<nHitsFit<<"   "<<i<<"   "<<endl;
        //out.close();
        hNHitsFitvsRunIndex->Fill(runIndex,nHitsFit);
		hNHitsdEdxvsRunIndex->Fill(runIndex,nHitsDedx);
	hNHitsPoss->Fill(charge*nHitsPoss);
        hNHitsFitoverNHitsPoss->Fill(charge*FitoverPoss);//new
        hNHitsFitoverNHitsPossvsRunIndex->Fill(runIndex,FitoverPoss);//new
        //fit/poss int to float
	hNHitsDedx->Fill(charge*nHitsDedx);
	hEtavsPhi->Fill(phi,eta);
	hEtavsPt->Fill(charge*pt,eta);
	hPhivsPt->Fill(charge*pt,phi);
	hDcavsPt->Fill(charge*pt, dca);
	if(beta2TOF>0.){
		hBetavsP->Fill(charge*p,1./beta2TOF);
		hBetavsEta->Fill(eta,1./beta2TOF);
		hBetavsPhi->Fill(phi,1./beta2TOF);
	}
	hDedxvsP->Fill(charge*p,dedx);
	hDedxvsEta->Fill(eta,dedx);
	hDedxvsPhi->Fill(phi,dedx);
	hNSigmaEvsP->Fill(charge*p,nSigmaE);
	hNSigmaEvsEta->Fill(eta,nSigmaE);
	hNSigmaEvsPhi->Fill(phi,nSigmaE);

	hPtvsRunIndex->Fill(runIndex,pt);
	hEtavsRunIndex->Fill(runIndex,eta);
	hPhivsRunIndex->Fill(runIndex,phi);
	hDcavsRunIndex->Fill(runIndex,dca);
	if(beta2TOF>0.) hBetavsRunIndex->Fill(runIndex,1./beta2TOF);
	hDedxvsRunIndex->Fill(runIndex,dedx);
	hNSigmaEvsRunIndex->Fill(runIndex,nSigmaE);
	// cout<< "reading track " << i << endl;

	Float_t msquare = -999.;
	if(beta2TOF>0.) msquare = pow(p,2)*(1-pow(beta2TOF,2))/pow(beta2TOF,2);
	if(msquare>-990.) hMSquarevsP->Fill(p,msquare);
	// cout<< "Filling track " << i << endl;


	Float_t expBeta2TOF = -999;

	Float_t mTpceNSigmaECutLow;
	if(p<.8){                                                                                          
		mTpceNSigmaECutLow = 3.0*p - 3.15; 
	}else{
		mTpceNSigmaECutLow = mTpceNSigmaECut[0];
	}

	hNSigmaEvsP_woTOF->Fill(p,nSigmaE);

	if(beta2TOF>0. && TMath::Abs(1.-1./beta2TOF)<=mTpceBeta2TOFCut)
	{
		hNSigmaEvsP_wTOF->Fill(p,nSigmaE);
	}

	if(beta2TOF>0.
			&& TMath::Abs(1.-1./beta2TOF)<=mTpceBeta2TOFCut 
			&& nSigmaE>=mTpceNSigmaECutLow && nSigmaE<=mTpceNSigmaECut[1]){
		hEEtavsPhi->Fill(phi,eta);
		hEEtavsPt->Fill(charge*pt,eta);
		hEPhivsPt->Fill(charge*pt,phi);
		expBeta2TOF = p/sqrt(pow(Melectron,2)+pow(p,2));
		hEBetaDiffvsRunIndex->Fill(runIndex,1./beta2TOF-1./expBeta2TOF);
		hEVxvsVy->Fill(mgOriginX,mgOriginY);
		hEVyvsVz->Fill(mgOriginY,mgOriginZ);
		hEVxvsVz->Fill(mgOriginX,mgOriginZ);
		
	}
	// cout<< "end fill track " << i << endl;

	return kTRUE;

}
//____________________________________________________________
void bookHistograms()
{
	//inclusive QA
	hnEvts = new TH1D("hnEvts","hnEvts",5,0.5,5.5);
	hnEvts->GetXaxis()->SetBinLabel(1,"nPicoEvents");
	hnEvts->GetXaxis()->SetBinLabel(2,"after Pileup rejection");
	hnEvts->GetXaxis()->SetBinLabel(3,"|Vz| < 35cm");
	hnEvts->GetXaxis()->SetBinLabel(4,"Vr < 2cm");
	hVyvsVx = new TH2F("hVyvsVx","hVyvsVx;Vx (cm); Vy (cm)",1000,-5,5,1000,-5,5);
	hVyvsVz = new TH2F("hVyvsVz","hVyvsVz;Vz (cm); Vy (cm)",1000,-5,5,1000,-5,5);
	hVr = new TH1F("hVr","hVr;Vr (cm);Counts",500,0,10);
        hVPDVzvsTPCVz = new TH2F("hVPDVzvsTPCVz","hVPDVzvsTPCVz; TPC Vz (cm); VPD Vz (cm)",400,-200,200,400,-200,200);
	hDeltaZvsTPCVz = new TH2F("hDeltaZvsTPCVz","hDeltaZvsTPCVz; TPC Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",400,-200,200,200,-10,10);
	hDeltaZvsVPDVz = new TH2F("hDeltaZvsVPDVz","hDeltaZvsVPDVz; VPD Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",400,-200,200,200,-10,10);
	hTPCVzvsRefMult = new TH2F("hTPCVzvsRefMult","hTPCVzvsRefMult; refMult; TPC Vz (cm)",200,0,1000,400,-200,200);
	hVPDVzvsRefMult = new TH2F("hVPDVzvsRefMult","hVPDVzvsRefMult; refMult; VPD Vz (cm)",200,0,1000,400,-200,200);
	hDeltaZvsRefMult = new TH2F("hDeltaZvsRefMult","hDeltaZvsRefMult; refMult; Vz_{TPC} - Vz_{VPD} (cm)",200,0,1000,200,-10,10);
	hZDCXvsRefMult = new TH2F("hZDCXvsRefMult","hZDCXvsRefMult; refMult; zdcRate (kHz)",1000,0,1000,1000,0,100);
	hBBCXvsRefMult = new TH2F("hBBCXvsRefMult","hBBCXvsRefMult; refMult; bbcRate (kHz)",1000,0,1000,1000,0,100);
	hnTofMatchvsRefMult_Rejection = new TH2D("hnTofMatchvsRefMult_Rejection","hnTofMatchvsRefMult_Rejection; refmult; nTOFmatch",600,0,600,600,0,600);
	hnTofMatchvsRefMult = new TH2D("hnTofMatchvsRefMult","hnTofMatchvsRefMult; refmult; nTOFmatch",600,0,600,600,0,600);
	hNHitsFit = new TH1F("hNHitsFit","hNHitsFit;nHitsFit;Counts",200,-100,100);
	hNHitsPoss = new TH1F("hNHitsPoss","hNHitsPoss;nHitsPoss;Counts",200,-100,100);
	hNHitsDedx = new TH1F("hNHitsDedx","hNHitsDedx;nHitsDedx;Counts",200,-100,100);
	hDcavsPt = new TH2F("hDcavsPt","hDcavsPt;q*p_{T} (GeV/c);dca (cm)",200,-10,10,300,-1.e-6,10.-1.e-6);
	hBetavsP = new TH2F("hBetavsP","hBetavsP;q*p (GeV/c);1/#beta",200,-10,10,500,0,5);
	hBetavsEta = new TH2F("hBetavsEta","hBetavsEta;#eta;1/#beta",400,-2,2,500,0,5);
	hBetavsPhi = new TH2F("hBetavsPhi","hBetavsPhi;#phi;1/#beta",360,-TMath::Pi(),TMath::Pi(),500,0,5);
	hDedxvsP = new TH2F("hDedxvsP","hDedxvsP;q*p (GeV/c);dEdx (KeV/cm)",200,-10,10,300,-1.e-6,30-1.e-6);
	hDedxvsEta = new TH2F("hDedxvsEta","hDedxvsEta;#eta;dEdx (KeV/cm)",400,-2,2,300,-1.e-6,30-1.e-6);
	hDedxvsPhi = new TH2F("hDedxvsPhi","hDedxvsPhi;#phi;dEdx (KeV/cm)",360,-TMath::Pi(),TMath::Pi(),300,-1.e-6,30-1.e-6);
	hNSigmaEvsP_woTOF = new TH2F("hNSigmaEvsP_woTOF",";q (GeV/c); n#sigma_{e}",1000,0,5,2000,-10+1.e-6,10+1.e-6);
	hNSigmaEvsP_wTOF = new TH2F("hNSigmaEvsP_wTOF",";q (GeV/c); n#sigma_{e}",1000,0,5,2000,-10+1.e-6,10+1.e-6);
	hNSigmaEvsP = new TH2F("hNSigmaEvsP","hNSigmaEvsP;q*p (GeV/c);n#sigma_{e}",2000,-10,10,4000,-20-1.e-6,20-1.e-6);
	hNSigmaEvsEta = new TH2F("hNSigmaEvsEta","hNSigmaEvsEta;#eta;n#sigma_{e}",400,-2,2,4000,-20-1.e-6,20-1.e-6);
	hNSigmaEvsPhi = new TH2F("hNSigmaEvsPhi","hNSigmaEvsPhi;#phi;n#sigma_{e}",360,-TMath::Pi(),TMath::Pi(),4000,-20-1.e-6,20-1.e-6);
	hMSquarevsP = new TH2F("hMSquarevsP","hMSquarevsP; p (GeV/c);m^{2} ( (GeV/c^{2})^{2} )",100,0,10,1200,-0.2,1);
	hEtavsPhi = new TH2F("hEtavsPhi","hEtavsPhi;#phi;#eta",360,-TMath::Pi(),TMath::Pi(),400,-2,2);
	hEEtavsPhi = new TH2F("hEEtavsPhi","hEEtavsPhi;#phi;#eta",360,-TMath::Pi(),TMath::Pi(),400,-2,2);
	hEtavsPt = new TH2F("hEtavsPt","hEtavsPt; q*p_{T} (GeV/c); #eta",200,-10,10,400,-2,2);
	hEEtavsPt = new TH2F("hEEtavsPt","hEEtavsPt; q*p_{T} (GeV/c); #eta",200,-10,10,400,-2,2);
	hPhivsPt = new TH2F("hPhivsPt","hPhivsPt; q*p_{T} (GeV/c); #phi",2000,-10,10,1080,-TMath::Pi(),TMath::Pi());
	hEPhivsPt = new TH2F("hEPhivsPt","hEPhivsPt; q*p_{T} (GeV/c); #phi",4000,-5,5,1080,-TMath::Pi(),TMath::Pi());
    hNHitsFitoverNHitsPoss = new TH1F("hNHitsFitoverNHitsPoss","hNHitsFitoverNHitsPoss;#frac{NHitsFit}{NHitsPoss};Counts",200,-1,1);//new
	hEVxvsVy = new TH2F("hEVxvsVy","hEVxvsVy; Vx (cm); Vy (cm)", 200,-10,10, 200, -10, 10);
	hEVyvsVz = new TH2F("hEVyvsVz","hEVyvsVz; Vz (cm); Vx (cm)", 6000, -300,300 , 200, -10,10);
	hEVxvsVz = new TH2F("hEVxvsVz","hEVxvsVz; Vz (cm); Vy (cm)", 6000, -300,300 , 200, -10,10);
	
	//run by run QA
	hZDCXvsRunIndex = new TProfile("hZDCXvsRunIndex","hZDCXvsRunIndex;runIndex;zdcRate (KHz)",mTotalRun,0,mTotalRun,0,100);
	hBBCXvsRunIndex = new TProfile("hBBCXvsRunIndex","hBBCXvsRunIndex;runIndex;bbcRate (KHz)",mTotalRun,0,mTotalRun,0,250);
        //bin set to 0.1cm 5.23
	hTPCVzvsRunIndex = new TProfile("hTPCVzvsRunIndex","hTPCVzvsRunIndex;runIndex;TPC Vz (cm)",mTotalRun,0,mTotalRun,-250,250);
	hVPDVzvsRunIndex = new TProfile("hVPDVzvsRunIndex","hVPDVzvsRunIndex;runIndex;VPD Vz (cm)",mTotalRun,0,mTotalRun,-250,250);
        hVrvsRunIndex = new TProfile("hVrvsRunIndex","hVrvsRunIndex;runIndex;Vr (cm))",mTotalRun,0,mTotalRun,0,5);//new
        hNHitsFitoverNHitsPossvsRunIndex = new TProfile("hNHitsFitoverNHitsPossvsRunIndex","hNHitsFitoverNHitsPossvsRunIndex;runIndex;#frac{NHitsFit}{NHitsPoss}",mTotalRun,0,mTotalRun,0,1.1);//new
	hDeltaZvsRunIndex = new TProfile("hDeltaZvsRunIndex","hDeltaZvsRunIndex;runIndex; Vz_{TPC} - Vz_{VPD} (cm)",mTotalRun,0,mTotalRun,-50,50);
	hRefMultvsRunIndex = new TProfile("hRefMultvsRunIndex","hRefMultvsRunIndex;runIndex; refMult",mTotalRun,0,mTotalRun,0,1000);
	hRefMultCorrvsRunIndex = new TProfile("hRefMultCorrvsRunIndex","hRefMultCorrvsRunIndex;runIndex; refMultCorr",mTotalRun,0,mTotalRun,0,1000);
	hPtvsRunIndex = new TProfile("hPtvsRunIndex","hPtvsRunIndex;runIndex;p_{T} (GeV/c)",mTotalRun,0,mTotalRun,0,10);//when I profile it ,it mean I add a cut form 0 to 10
	hEtavsRunIndex = new TProfile("hEtavsRunIndex","hEtavsRunIndex;runIndex;#eta",mTotalRun,0,mTotalRun,-2,2);
	hPhivsRunIndex = new TProfile("hPhivsRunIndex","hPhivsRunIndex;runIndex;#phi",mTotalRun,0,mTotalRun,-TMath::Pi(),TMath::Pi());
	hDcavsRunIndex = new TProfile("hDcavsRunIndex","hDcavsRunIndex;runIndex;dca (cm)",mTotalRun,0,mTotalRun,-1.e-6,5-1.e-6);//when I profile it ,it mean I add a cut form 0 to 3
	hBetavsRunIndex = new TProfile("hBetavsRunIndex","hBetavsRunIndex;runIndex;1/#beta",mTotalRun,0,mTotalRun,0,5);
	hDedxvsRunIndex = new TProfile("hDedxvsRunIndex","hDedxvsRunIndex;runIndex;dE/dX [KeV/cm]",mTotalRun,0,mTotalRun,0,15);
	hEBetaDiffvsRunIndex = new TProfile("hEBetaDiffvsRunIndex","hEBetaDiffvsRunIndex;runIndex;1/#beta - 1/#beta_{exp}",mTotalRun,0,mTotalRun,-0.05,0.05);
	hKBetaDiffvsRunIndex = new TProfile("hKBetaDiffvsRunIndex","hKBetaDiffvsRunIndex;runIndex;1/#beta - 1/#beta_{exp}",mTotalRun,0,mTotalRun,-0.05,0.05);
	hNSigmaEvsRunIndex = new TProfile("hNSigmaEvsRunIndex","hNSigmaEvsRunIndex;runIndex;n#sigma_{e}",mTotalRun,0,mTotalRun,-30-1.e-6,30-1.e-6);
    hNHitsFitvsRunIndex = new TProfile("hNHitsFitvsRunIndex","hNHitsFitvsRunIndex;runIndex;NHitsFit",mTotalRun,0,mTotalRun,0,100);
    hNHitsdEdxvsRunIndex = new TProfile("hNHitsDedxvsRunIndex","hNHitsDedxvsRunIndex;runIndex;NHitsDedx",mTotalRun,0,mTotalRun,0,100);
	hnTrkvsRunIndex = new TProfile("hnTrkvsRunIndex","hnTrkvsRunIndex;runIndex;nTrk",mTotalRun,0,mTotalRun,0,100);
}
//=======================================================================================
void writeHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.QAhisto.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	TFile *mFile = new TFile(buf,"recreate");
	mFile->cd();

	//inclusive QA
	hnEvts->Write();
	hVyvsVx->Write();
	hVyvsVz->Write();
        hVr->Write();
	hVPDVzvsTPCVz->Write();
	hDeltaZvsTPCVz->Write();
	hDeltaZvsVPDVz->Write();
	hTPCVzvsRefMult->Write();
	hVPDVzvsRefMult->Write();
	hDeltaZvsRefMult->Write();
	hZDCXvsRefMult->Write();
	hBBCXvsRefMult->Write();
	hnTofMatchvsRefMult->Write();
	hnTofMatchvsRefMult_Rejection->Write();

	hNHitsFit->Write();
	hNHitsPoss->Write();
	hNHitsDedx->Write();
	hDcavsPt->Write();
	hBetavsP->Write();
	hBetavsEta->Write();
	hBetavsPhi->Write();
	hDedxvsP->Write();
	hDedxvsEta->Write();
	hDedxvsPhi->Write();
	hNSigmaEvsP_woTOF->Write();
	hNSigmaEvsP_wTOF->Write();
	hNSigmaEvsP->Write();
	hNSigmaEvsEta->Write();
	hNSigmaEvsPhi->Write();
	hMSquarevsP->Write();
	hEtavsPhi->Write();
	hEEtavsPhi->Write();
	hEtavsPt->Write();
	hEEtavsPt->Write();
	hPhivsPt->Write();
	hEPhivsPt->Write();
        hNHitsFitoverNHitsPoss->Write();
	hEVxvsVy->Write();
	hEVyvsVz->Write();
	hEVxvsVz->Write();
	//run by run QA
	hZDCXvsRunIndex->Write();
	hBBCXvsRunIndex->Write();
        hVrvsRunIndex->Write();
        hNHitsFitoverNHitsPossvsRunIndex->Write();
	hTPCVzvsRunIndex->Write();
	hVPDVzvsRunIndex->Write();
        hNHitsFitvsRunIndex->Write();
	hDeltaZvsRunIndex->Write();
	hRefMultvsRunIndex->Write();
	hRefMultCorrvsRunIndex->Write();
	hPtvsRunIndex->Write();
	hEtavsRunIndex->Write();
	hPhivsRunIndex->Write();
	hDcavsRunIndex->Write();
	hBetavsRunIndex->Write();
	hEBetaDiffvsRunIndex->Write();
	hKBetaDiffvsRunIndex->Write();
	hDedxvsRunIndex->Write();
	hNSigmaEvsRunIndex->Write();
	hNHitsdEdxvsRunIndex->Write();
	hnTrkvsRunIndex->Write();
}
//==============================================================================================
bool Init()
{
	ifstream indata;
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
	mTotalRun = mTotalRunId.size();
	cout << "Totally we have" << mTotalRun << " Runs!" << endl;
	for(map<Int_t,Int_t>::iterator iter=mTotalRunId.begin();iter!=mTotalRunId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;
	cout<<endl;

	Pileuplimit = new TF1("Pileuplimit","0.7*x-10",0,1000);

	return kTRUE;
}
