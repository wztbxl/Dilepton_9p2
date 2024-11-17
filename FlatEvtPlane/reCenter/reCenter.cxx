#include "/star/u/wangzhen/run20/Dielectron/func_macro/headers.h"
#include "/star/u/wangzhen/run20/Dielectron/func_macro/function.C"
#include "/star/u/wangzhen/run20/Dielectron/Analysis/cuts.h"
#include "/star/u/wangzhen/run20/Dielectron/Analysis/RefMfun.h"

#include "CentralityMaker.h"
#include "StRefMultCorr.h"
#include "miniDst.h"

TFile *fOutFile;

map<Int_t,Int_t> mTotalRunId;
map<Int_t,Int_t> mBadRunId_001;
map<Int_t,Int_t> mBadRunId_021;

bool Init();
// Double_t GetRefMultCorr(const Int_t RefMult, const Double_t z);
// Double_t GetWeight(Double_t refMultCorr);
// Int_t GetCentrality(Double_t refMultCorr);
void bookHistograms(char* outFile);
bool passEvent(miniDst const* const event);

//define histograms
TProfile2D *etapluszplusQx;
TProfile2D *etapluszminusQx;
TProfile2D *etaminuszplusQx;
TProfile2D *etaminuszminusQx;
TProfile2D *etapluszplusQy;
TProfile2D *etapluszminusQy;
TProfile2D *etaminuszplusQy;
TProfile2D *etaminuszminusQy;
TProfile2D *etaplusQx;
TProfile2D *etaminusQx;
TProfile2D *etaplusQy;
TProfile2D *etaminusQy;
TProfile *etaplusQx_cent;
TProfile *etaminusQx_cent;
TProfile *etaplusQy_cent;
TProfile *etaminusQy_cent;
TProfile *etaplusQx_cent_RejectE;
TProfile *etaminusQx_cent_RejectE;
TProfile *etaplusQy_cent_RejectE;
TProfile *etaminusQy_cent_RejectE;
TH1F *hRefMult;
TH1F *hRefMultCorZ;
TH1F *hRefMultCor;
TH1F *hCentrality;
TH1F *hCentralityCor;
TH3F *hQXvsQYvsRunIndex;
// TH3F *hQXvsQYvsRunIndex_runindex;
TH3F* hQXvsQYvsCentrality_east;
TH3F* hQXvsQYvsCentrality_west;
// TH3F *hQXvsQYvsRunIndex_east;
// TH3F *hQXvsQYvsRunIndex_west;
TH3F *hQXvsQYvsCent_east;
TH3F *hQXvsQYvsCent_west;

TH2D *hnTofHitsvsRefMult;
TH2D* hnTofHitsvsRefMult_noCut;
TH2D* hnTofHitsvsRefMult_Vz35;
TH3D* hnPlusTrkvsnMinusTrkvsCent;
TF1* Pileuplimit;
TF1* PileupUplimit;
TF1* PileupLowlimit;


int main(int argc, char** argv)
{
	if(argc!=1&&argc!=3) return -1;

	TString inFile="test.list";
	char outFile[1024];
	sprintf(outFile,"test");
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
	if( !Init() ){
		cout<<"The initialization is failed !!!"<<endl;
		return 0;
	}
	bookHistograms(outFile);

	//+-------------+
	//| loop events |
	//+-------------+
	miniDst *event = new miniDst(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;

	for(int i=0;i<nEvts;i++){

		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;

		event->GetEntry(i);

		if(!passEvent(event)) continue;
	}

	if(fOutFile){
		fOutFile->cd();
		fOutFile->Write();
		fOutFile->Close();
	}

	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
bool passEvent(miniDst const* const event)
{
	Bool_t validTrig = kFALSE;
	Bool_t RefMVzCorFlag = kFALSE;
	Bool_t is001Trigger = kFALSE;
	Bool_t is021Trigger = kFALSE;
	for(Int_t i=0; i<event->mNTrigs; i++){
		if(event->mTrigId[i] == 780010) validTrig = kTRUE, is001Trigger = kTRUE, RefMVzCorFlag= kTRUE;
		//if(event->mTrigId[i] == 580011) validTrig = kTRUE;
		if(event->mTrigId[i] == 780020) validTrig = kTRUE, is021Trigger = kTRUE;
	}
	if(!validTrig){
		return kFALSE;
	}

	Int_t runId  = event->mRunId;

	Int_t runIndex;
	map<Int_t, Int_t>::iterator iter = mTotalRunId.find(runId);
	if(iter != mTotalRunId.end()){
		runIndex = iter->second;
	}
	else{
		cout<<"Can not find the runNumber in the mTotalRunId list"<<endl;
		return kFALSE;
	}

	map<Int_t, Int_t>::iterator iter_001 = mBadRunId_001.find(runId);
	if(iter_001 != mBadRunId_001.end() && validTrig){
		//cout<<"bad run, continue"<<endl;
		return kFALSE;
	}
	Int_t   mCentrality  = event->mCentrality;
	Int_t	refMult 	 = event->mRefMult;
  Int_t mnTOFMatch = event->mnTOFMatch;
  Float_t zdcRate = event->mZDCRate;
	Float_t vx           = event->mVertexX;
	Float_t vy           = event->mVertexY;
	Float_t vz           = event->mVertexZ;
	Float_t vpdVz        = event->mVpdVz;
	Float_t vr           = sqrt(vx*vx + vy*vy);    
	Float_t vzDiff       = vz - vpdVz;

	//for  the official centrality defination
	StRefMultCorr* mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
	// cout << "after refMultCorr defination" << endl;
	//using offical badrun list
	mRefMultCorr->init((Int_t)runId);
	// cout << "after refMultCorr init" << endl;
	mRefMultCorr->initEvent(refMult,vz,zdcRate);
	// cout << "after refMultCorr initEvent" << endl;
	if (mRefMultCorr->isBadRun(runId))
	{
		return kFALSE;
	}
	// cout << "after refMultCorr isBadRun" << endl;
	Double_t RefMultCorr  = mRefMultCorr->getRefMultCorr();
	// cout << "after refMultCorr getRefMultCorr " << endl;
	Double_t reweight  = mRefMultCorr->getWeight();
	// cout << "after refMultCorr getWeight" << endl;
	mCentrality = mRefMultCorr->getCentralityBin9();//9 Centrality bin
	// cout << "after refMultCorr getCentralityBin9" << endl;
	//offical pile up pileupRejection
	if  ( mRefMultCorr->isPileUpEvent(refMult,mnTOFMatch,vz ) ) return kFALSE;
    // Double_t RefMultCorr = refMult;
	// if(RefMVzCorFlag)RefMultCorr = GetRefMultCorr(refMult, vz);	
    // Double_t reweight = 1.;// no RefMultCorr!!!!	
    // Double_t reweight = GetWeight(RefMultCorr);	
	// mCentrality = GetCentrality(RefMultCorr);


	Float_t  mEtaPlusQx        = event->mEtaPlusQx;
	Float_t  mEtaPlusQy        = event->mEtaPlusQy;
	Float_t  mEtaPlusPtWeight  = event->mEtaPlusPtWeight;
	Short_t  mEtaPlusNTrks     = event->mEtaPlusNTrks;
	Float_t  mEtaMinusQx       = event->mEtaMinusQx;
	Float_t  mEtaMinusQy       = event->mEtaMinusQy;
	Float_t  mEtaMinusPtWeight = event->mEtaMinusPtWeight;
	Short_t  mEtaMinusNTrks    = event->mEtaMinusNTrks;

	//for the electron rejection
	Float_t  mEtaPlusQx_rejectE        = event->mEtaPlusQx;
	Float_t  mEtaPlusQy_rejectE        = event->mEtaPlusQy;
	Float_t  mEtaPlusPtWeight_rejectE  = event->mEtaPlusPtWeight;
	Float_t  mEtaMinusQx_rejectE       = event->mEtaMinusQx;
	Float_t  mEtaMinusQy_rejectE       = event->mEtaMinusQy;
	Float_t  mEtaMinusPtWeight_rejectE = event->mEtaMinusPtWeight;


  hnTofHitsvsRefMult_noCut->Fill(refMult,mnTOFMatch);
	if(vr>mVrCut)                     return kFALSE;
	if(TMath::Abs(vz)>mVzCut)         return kFALSE;
	// if(TMath::Abs(vzDiff)>mVzDiffCut) return kFALSE;
	// if (mnTOFMatch < Pileuplimit->Eval(refMult)) return kFALSE;
 	 hnTofHitsvsRefMult->Fill(refMult,mnTOFMatch);
	  if(mCentrality<0)                 return kFALSE;
	hRefMult->Fill(refMult);
	hRefMultCorZ->Fill(RefMultCorr);
	hRefMultCor->Fill(RefMultCorr, reweight);
	hCentrality->Fill(mCentrality);
	hCentralityCor->Fill(mCentrality,reweight);
	// if(mCentrality < 1 || mCentrality > 9 ) return kFALSE;

	// if(vz>0){
	// 	if(mEtaPlusNTrks>0){
	// 		etapluszplusQx->Fill(runIndex, mCentrality, mEtaPlusQx/mEtaPlusNTrks, mEtaPlusNTrks);
	// 		etapluszplusQy->Fill(runIndex, mCentrality, mEtaPlusQy/mEtaPlusNTrks, mEtaPlusNTrks);
	// 	}

	// 	if(mEtaMinusNTrks>0){
	// 		etaminuszplusQx->Fill(runIndex, mCentrality, mEtaMinusQx/mEtaMinusNTrks, mEtaMinusNTrks);
	// 		etaminuszplusQy->Fill(runIndex, mCentrality, mEtaMinusQy/mEtaMinusNTrks, mEtaMinusNTrks);
	// 	}
	// }
	// else{
	// 	if(mEtaPlusNTrks>0){
	// 		etapluszminusQx->Fill(runIndex, mCentrality, mEtaPlusQx/mEtaPlusNTrks, mEtaPlusNTrks);
	// 		etapluszminusQy->Fill(runIndex, mCentrality, mEtaPlusQy/mEtaPlusNTrks, mEtaPlusNTrks);
	// 	}

	// 	if(mEtaMinusNTrks>0){
	// 		etaminuszminusQx->Fill(runIndex, mCentrality, mEtaMinusQx/mEtaMinusNTrks, mEtaMinusNTrks);
	// 		etaminuszminusQy->Fill(runIndex, mCentrality, mEtaMinusQy/mEtaMinusNTrks, mEtaMinusNTrks);
	// 	}
	// }

	//add the electron selection cut to reject the electrons
	Int_t npTrks = event->mNTrks;
	for(int j=0;j<npTrks;j++)
	{
		Int_t charge = event->mCharge[j];
		Int_t nHitsFit = event->mNHitsFit[j];
		Int_t nHitsDedx = event->mNHitsDedx[j];
		Int_t nHitsPoss = event->mNHitsPoss[j];
		Float_t nSigmaE = event->mNSigmaE[j];
		Float_t dca = event->mDca[j];
		Float_t pt = event->mPt[j];
		Float_t eta = event->mEta[j];
		Float_t phi = event->mPhi[j];
		Float_t beta2TOF = event->mBeta2TOF[j];
		Float_t TOFLoaclY = event->mTOFLocalY[j];
		Float_t ratio = 1.0*nHitsFit/nHitsPoss;
		int CellID = event->mTOFCellID[j];
		TVector3 mom;
		mom.SetPtEtaPhi(pt,eta,phi);
		Float_t p = mom.Mag();
		double msquare =  -999;
		msquare = pow(p, 2) * (1 - pow(beta2TOF, 2)) / pow(beta2TOF, 2);

		if(pt<0.2 || pt>30.) continue;;
		// if(nHitsFit<15) continue;;
		if(nHitsFit<20) continue;;
		if(ratio<0.52) continue;;
		// if(nHitsDedx<20) continue;;
		if(nHitsDedx<15) continue;;
		// if(dca>0.8) continue;;
		if(dca>1.) continue;;
		if(TMath::Abs(eta)>1.) continue;;
		if(beta2TOF<=0. || TMath::Abs(1.-1./beta2TOF)>0.025) continue;;
		if(abs(TOFLoaclY) > 1.8) continue;;
		Float_t mTpceNSigmaECutLow;
		if(p<.8){
			mTpceNSigmaECutLow = 3.0*p - 3.15; 
		}else{
			mTpceNSigmaECutLow = -0.75;
		}
		if(nSigmaE<mTpceNSigmaECutLow+mNSigmaEShift || nSigmaE>2.0+mNSigmaEShift) continue;;
		
		if (eta < 0)
			{
				mEtaMinusQx_rejectE -= pt*TMath::Cos(2*phi);
				mEtaMinusQy_rejectE -= pt*TMath::Sin(2*phi);

			} else if (eta > 0)
			{
				mEtaPlusQx_rejectE -= pt*TMath::Cos(2*phi);
				mEtaPlusQy_rejectE -= pt*TMath::Sin(2*phi);
			}

	}

	if(mEtaPlusNTrks > 0)
	{
		etaplusQx_cent_RejectE->Fill(mCentrality,mEtaPlusQx_rejectE);
		etaplusQy_cent_RejectE->Fill(mCentrality,mEtaPlusQy_rejectE);
	} 
	if(mEtaMinusNTrks > 0)
	{
		etaminusQx_cent_RejectE->Fill(mCentrality,mEtaMinusQx_rejectE);
		etaminusQy_cent_RejectE->Fill(mCentrality,mEtaMinusQy_rejectE);
	}
	hnPlusTrkvsnMinusTrkvsCent->Fill(mEtaPlusNTrks,mEtaMinusNTrks,mCentrality);

	Float_t Qx = mEtaPlusQx_rejectE + mEtaMinusQx_rejectE; 
	Float_t Qy = mEtaPlusQy_rejectE + mEtaMinusQy_rejectE;
	TVector2 Q(Qx,Qy);
	TVector2 mRawQWest(mEtaPlusQx_rejectE,mEtaPlusQy_rejectE);
	TVector2 mRawQEast(mEtaMinusQx_rejectE,mEtaMinusQy_rejectE);
	if (mRawQEast.Mod() <= 0 || mRawQWest.Mod() <= 0 ) return kFALSE;

	//pT weight
	if(vz>0){
		if(mEtaPlusNTrks>0){
			etapluszplusQx->Fill(runIndex, mCentrality, mEtaPlusQx);
			etapluszplusQy->Fill(runIndex, mCentrality, mEtaPlusQy);
			// hQXvsQYvsRunIndex_west->Fill(mEtaPlusQx,mEtaPlusQy,runIndex);
			etaplusQx->Fill(runIndex,mCentrality,mEtaPlusQx);
			etaplusQy->Fill(runIndex,mCentrality,mEtaPlusQy);
			etaplusQx_cent->Fill(mCentrality,mEtaPlusQx);
			etaplusQy_cent->Fill(mCentrality,mEtaPlusQy);
			hQXvsQYvsCent_west->Fill(mEtaPlusQx,mEtaPlusQy,mCentrality);
			
		}

		if(mEtaMinusNTrks>0){
			etaminuszplusQx->Fill(runIndex, mCentrality, mEtaMinusQx);
			etaminuszplusQy->Fill(runIndex, mCentrality, mEtaMinusQy);
			// hQXvsQYvsRunIndex_east->Fill(mEtaMinusQx,mEtaMinusQy,runIndex);
			etaminusQx->Fill(runIndex,mCentrality,mEtaMinusQx);
			etaminusQy->Fill(runIndex,mCentrality,mEtaMinusQy);
			etaminusQx_cent->Fill(mCentrality,mEtaMinusQx);
			etaminusQy_cent->Fill(mCentrality,mEtaMinusQy);
			hQXvsQYvsCent_east->Fill(mEtaMinusQx,mEtaMinusQy,mCentrality);
			
		}
	}
	else{
		if(mEtaPlusNTrks>0){
			etapluszminusQx->Fill(runIndex, mCentrality, mEtaPlusQx);
			etapluszminusQy->Fill(runIndex, mCentrality, mEtaPlusQy);
			// hQXvsQYvsRunIndex_west->Fill(mEtaPlusQx,mEtaPlusQy,runIndex);
			etaplusQx->Fill(runIndex,mCentrality,mEtaPlusQx);
			etaplusQy->Fill(runIndex,mCentrality,mEtaPlusQy);
			etaplusQx_cent->Fill(mCentrality,mEtaPlusQx);
			etaplusQy_cent->Fill(mCentrality,mEtaPlusQy);
			hQXvsQYvsCent_west->Fill(mEtaPlusQx,mEtaPlusQy,mCentrality);
		}

		if(mEtaMinusNTrks>0){
			etaminuszminusQx->Fill(runIndex, mCentrality, mEtaMinusQx);
			etaminuszminusQy->Fill(runIndex, mCentrality, mEtaMinusQy);
			// hQXvsQYvsRunIndex_east->Fill(mEtaMinusQx,mEtaMinusQy,runIndex);
			etaminusQx->Fill(runIndex,mCentrality,mEtaMinusQx);
			etaminusQy->Fill(runIndex,mCentrality,mEtaMinusQy);
			etaminusQx_cent->Fill(mCentrality,mEtaMinusQx);
			etaminusQy_cent->Fill(mCentrality,mEtaMinusQy);
			hQXvsQYvsCent_east->Fill(mEtaMinusQx,mEtaMinusQy,mCentrality);
		}
	}
	
	Qx = mEtaPlusQx+mEtaMinusQx;
	Qy = mEtaPlusQy+mEtaMinusQy;

	hQXvsQYvsRunIndex->Fill(Qx,Qy,mCentrality);
	hQXvsQYvsCentrality_east->Fill(mEtaMinusQx,mEtaMinusQy,mCentrality);
	hQXvsQYvsCentrality_west->Fill(mEtaPlusQx,mEtaPlusQy,mCentrality);
	// hQXvsQYvsRunIndex_runindex->Fill(Qx,Qy,runIndex);

	return kTRUE;
}
//____________________________________________________________
void bookHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.histo.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	fOutFile = new TFile(buf,"recreate");

	cout << "nRuns = " << mTotalRunId.size() << endl;
	const Int_t mTotalRun = mTotalRunId.size();
	const Int_t mTotalCentrality = 12;

	etapluszplusQx = new TProfile2D("etapluszplusQx","etapluszplusQx;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etapluszminusQx = new TProfile2D("etapluszminusQx","etapluszminusQx;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etaminuszplusQx = new TProfile2D("etaminuszplusQx","etaminuszplusQx;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etaminuszminusQx = new TProfile2D("etaminuszminusQx","etaminuszminusQx;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality); 

	etapluszplusQy = new TProfile2D("etapluszplusQy","etapluszplusQy;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etapluszminusQy = new TProfile2D("etapluszminusQy","etapluszminusQy;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etaminuszplusQy = new TProfile2D("etaminuszplusQy","etaminuszplusQy;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etaminuszminusQy = new TProfile2D("etaminuszminusQy","etaminuszminusQy;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);

	etaplusQx = new TProfile2D("etaplusQx","etaplusQx;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etaminusQx = new TProfile2D("etaminusQx","etaminusQx;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etaplusQy = new TProfile2D("etaplusQy","etaplusQy;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);
	etaminusQy = new TProfile2D("etaminusQy","etaminusQy;runIndex;Centrality",mTotalRun,0,mTotalRun,mTotalCentrality,0,mTotalCentrality);

	etaplusQx_cent = new TProfile("etaplusQx_cent","etaplusQx_cent;Centrality;Q_{x}^{#eta>0} ",mTotalCentrality,0,mTotalCentrality);
	etaminusQx_cent = new TProfile("etaminusQx_cent","etaminusQx_cent;Centrality;Q_{x}^{#eta<0} ",mTotalCentrality,0,mTotalCentrality);
	etaplusQy_cent = new TProfile("etaplusQy_cent","etaplusQy_cent;Centrality;Q_{y}^{#eta>0} ",mTotalCentrality,0,mTotalCentrality);
	etaminusQy_cent = new TProfile("etaminusQy_cent","etaminusQy_cent;Centrality;Q_{y}^{#eta<0} ",mTotalCentrality,0,mTotalCentrality);
	etaplusQx_cent_RejectE = new TProfile("etaplusQx_cent_RejectE","etaplusQx_cent_rejectE;Centrality;Q_{x}^{#eta>0} ",mTotalCentrality,0,mTotalCentrality);
	etaminusQx_cent_RejectE = new TProfile("etaminusQx_cent_RejectE","etaminusQx_cent_rejectE;Centrality;Q_{x}^{#eta<0} ",mTotalCentrality,0,mTotalCentrality);
	etaplusQy_cent_RejectE = new TProfile("etaplusQy_cent_RejectE","etaplusQy_cent_rejectE;Centrality;Q_{y}^{#eta>0} ",mTotalCentrality,0,mTotalCentrality);
	etaminusQy_cent_RejectE = new TProfile("etaminusQy_cent_RejectE","etaminusQy_cent_rejectE;Centrality;Q_{y}^{#eta<0} ",mTotalCentrality,0,mTotalCentrality);

	hQXvsQYvsRunIndex = new TH3F("hQXvsQYvsRunIndex","; Qx; Qy; Centrality",200,-20,20,200,-20,20,10,0,10);
	hQXvsQYvsCentrality_east = new TH3F("hQXvsQYvsCentrality_east","; Qx; Qy; Centrality",200,-20,20,200,-20,20,10,0,10);
	hQXvsQYvsCentrality_west = new TH3F("hQXvsQYvsCentrality_west","; Qx; Qy; Centrality",200,-20,20,200,-20,20,10,0,10);
	// hQXvsQYvsRunIndex_runindex = new TH3F("hQXvsQYvsRunIndex_runindex","; Qx; Qy; runindex",200,-20,20,200,-20,20,mTotalRun,0,mTotalRun);
	// hQXvsQYvsRunIndex_east = new TH3F("hQXvsQYvsRunIndex_east","; Qx; Qy; runindex",200,-20,20,200,-20,20,mTotalRun,0,mTotalRun);
	// hQXvsQYvsRunIndex_west = new TH3F("hQXvsQYvsRunIndex_west","; Qx; Qy; runindex",200,-20,20,200,-20,20,mTotalRun,0,mTotalRun);
	hQXvsQYvsCent_east = new TH3F("hQXvsQYvsCent_east","; Qx; Qy; Centrality",200,-20,20,200,-20,20,mTotalCentrality,0,mTotalCentrality);
	hQXvsQYvsCent_west = new TH3F("hQXvsQYvsCent_west","; Qx; Qy; Centrality",200,-20,20,200,-20,20,mTotalCentrality,0,mTotalCentrality);

	hRefMult  = new TH1F("RefMult","RefMult;RefMult",600,0,600);
	hRefMultCorZ  = new TH1F("RefMultCorZ","RefMult corrected for z dependence ;RefMult",600,0,600);
	hRefMultCor  = new TH1F("RefMultCor","RefMult corrected for z and weight ;Mult",600,0,600);
	hCentrality = new TH1F("Centrality","Centrality; centrality",16,0,16);
	hCentralityCor = new TH1F("CentralityCor","Centrality with weight; centrality",16,0,16);
  hnTofHitsvsRefMult = new TH2D("hnTofHitsvsRefMult",";RefMult;nTofHits",500,0,500,500,0,500);
  hnTofHitsvsRefMult_noCut = new TH2D("hnTofHitsvsRefMult_noCut",";RefMult;nTofHits",500,0,500,500,0,500); 
  hnTofHitsvsRefMult_Vz35 = new TH2D("hnTofHitsvsRefMult_Vz35",";RefMult;nTofHits",500,0,500,500,0,500);
  hnPlusTrkvsnMinusTrkvsCent = new TH3D("hnPlusTrkvsnMinusTrkvsCent","nPlusTrk;nMinusTrk;Cent",500,0,500,500,0,500,10,0,10);
  Pileuplimit = new TF1("Pileuplimit","0.7*x-10",0,1000);

}
//____________________________________________________________
bool Init()
{
	cout<<endl;

	ifstream indata;
	indata.open("/star/u/wangzhen/run20/Dielectron/DataQA/mTotalRunList.dat");
	mTotalRunId.clear();
	if(indata.is_open()){
		cout<<"read in total run number list and recode run number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata>>oldId){
			mTotalRunId[oldId] = newId;
			cout << "Run " << oldId << " index is " << newId << endl;
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
		cout<<"read in total run number list and recode run number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata_001>>oldId){
			mBadRunId_001[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total bad run number list !!!"<<endl;
		return kFALSE;
	}
	indata_001.close();
	
	cout<< "run list: " <<endl;
	for(map<Int_t,Int_t>::iterator iter=mTotalRunId.begin();iter!=mTotalRunId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;

	cout<<"bad run list:"<<endl;
	for(map<Int_t,Int_t>::iterator iter=mBadRunId_001.begin();iter!=mBadRunId_001.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;

  return kTRUE;
}
//---------------------------------------------------
