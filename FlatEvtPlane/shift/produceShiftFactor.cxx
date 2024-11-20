#include "/star/u/wangzhen/run20/Dielectron/func_macro/headers.h"                                                
#include "/star/u/wangzhen/run20/Dielectron/func_macro/function.C"
#include "/star/u/wangzhen/run20/Dielectron/Analysis/cuts.h"
#include "/star/u/wangzhen/run20/Dielectron/Analysis/RefMfun.h"

#include "CentralityMaker.h"
#include "StRefMultCorr.h"
#include "miniDst.h"

TFile *fOutFile;

map<Int_t,Int_t> mTotalRunId;
map<Int_t,Int_t> mTotalDayId;
map<Int_t,Int_t> mBadRunId_001;
map<Int_t,Int_t> mBadRunId_021;
TF1* Pileuplimit;

bool Init();
// Double_t GetRefMultCorr(const Int_t RefMult, const Double_t z);
// Double_t GetWeight(Double_t refMultCorr);
// Int_t GetCentrality(Double_t refMultCorr);
void bookHistograms(char* outFile);
bool passEvent(miniDst const* const event);

TFile *fReCenter;
TProfile2D *etapluszplusQx;
TProfile2D *etapluszminusQx;
TProfile2D *etaminuszplusQx;
TProfile2D *etaminuszminusQx;
TProfile2D *etapluszplusQy;
TProfile2D *etapluszminusQy;
TProfile2D *etaminuszplusQy;
TProfile2D *etaminuszminusQy;
TProfile *etaplusQx_cent;
TProfile *etaminusQx_cent;
TProfile *etaplusQy_cent;
TProfile *etaminusQy_cent;
TProfile *etaplusQx_cent_rejectE;
TProfile *etaminusQx_cent_rejectE;
TProfile *etaplusQy_cent_rejectE;
TProfile *etaminusQy_cent_rejectE;

const int nCent = 9;
//define histograms
TH1D *hRawEventPlane;
TH1D *hReCenterEventPlane;
TH1D* hReCenterEventPlane_rejectE;
TH2D* hRawQxQy[nCent];
TH2D* hRecenterQxQy[nCent];
TH1D* hReCenterEventPlaneEast;
TH1D* hReCenterEventPlaneWest;
TProfile2D *shiftfactorcos[mArrayLength];
TProfile2D *shiftfactorsin[mArrayLength];
TProfile *shiftfactorcos_cent[mArrayLength];
TProfile *shiftfactorsin_cent[mArrayLength];
TProfile *shiftfactorcos_cent_rejectE[mArrayLength];
TProfile *shiftfactorsin_cent_rejectE[mArrayLength];
TProfile *shiftfactorcos_cent_east_rejectE[mArrayLength];
TProfile *shiftfactorcos_cent_west_rejectE[mArrayLength];
TProfile *shiftfactorsin_cent_east_rejectE[mArrayLength];
TProfile *shiftfactorsin_cent_west_rejectE[mArrayLength];


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
	if(runId<mMinRunId) return kFALSE;

    Int_t runIndex;
	map<Int_t, Int_t>::iterator iter = mTotalRunId.find(runId);
	if(iter != mTotalRunId.end()){
		runIndex = iter->second;
	}
	else{
		cout<<"Can not find the runNumber in the mTotalRunId list"<<endl;
		return kFALSE;
	}

	Int_t dayIndex;
	iter = mTotalDayId.find((runId/1000)%1000);
	if(iter != mTotalDayId.end()){
		dayIndex = iter->second;
	}
	else{
		cout<<"Can not find the dayNumber in the mTotalDayId list"<<endl;
		return kFALSE;
	}

	map<Int_t, Int_t>::iterator iter_001 = mBadRunId_001.find(runId);
	if(iter_001 != mBadRunId_001.end() && is001Trigger){
		//cout<<"bad run, continue"<<endl;
		return kFALSE;
	}


	Int_t   mCentrality  = event->mCentrality;
	Double_t vx           = event->mVertexX;
	Double_t vy           = event->mVertexY;
	Double_t vz           = event->mVertexZ;
	Int_t  refMult       = event->mRefMult;
	Int_t mnTOFMatch = event->mnTOFMatch;
	Float_t zdcRate = event->mZDCRate;
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
	// cout << "centrality = " << mCentrality << endl;
	// cout << "after refMultCorr getCentralityBin9" << endl;
	//offical pile up pileupRejection
	if  ( mRefMultCorr->isPileUpEvent(refMult,mnTOFMatch,vz ) ) return kFALSE;

	if(vr>mVrCut)                     return kFALSE;
	if(TMath::Abs(vz)>mVzCut)         return kFALSE;
	// if(TMath::Abs(vzDiff)>mVzDiffCut) return kFALSE;
	if(mCentrality<0)                 return kFALSE;

    // Double_t RefMultCorr = refMult;
	// if(RefMVzCorFlag)RefMultCorr = GetRefMultCorr(refMult, vz);	
    // Double_t reweight = GetWeight(RefMultCorr);	
	// hnTofHitsvsRefMult_noCut->Fill(refMult,mnTOFMatch);
	if(vr>mVrCut)                     return kFALSE;
	if(TMath::Abs(vz)>mVzCut)         return kFALSE;
	// if(TMath::Abs(vzDiff)>mVzDiffCut) return kFALSE;
	// if (mnTOFMatch < Pileuplimit->Eval(refMult)) return kFALSE;
 	//  hnTofHitsvsRefMult->Fill(refMult,mnTOFMatch);
	  if(mCentrality<0)                 return kFALSE;
	// mCentrality = GetCentrality(RefMultCorr);
	// if(mCentrality < 1 || mCentrality > 9 ) return kFALSE;

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
	double mPlusQx = mEtaPlusQx;
	double mPlusQy = mEtaPlusQy;
	double mMinusQx = mEtaMinusQx;
	double mMinusQy = mEtaMinusQy;

	Double_t mRawQx = mPlusQx + mMinusQx; 
	Double_t mRawQy = mPlusQy + mMinusQy;
	hRawQxQy[mCentrality]->Fill(mRawQx,mRawQy);

	TVector2 *mRawQ = new TVector2(mRawQx, mRawQy);
	if(mRawQ->Mod() > 0){
		Double_t mRawEventPlane = 0.5*TMath::ATan2(mRawQy,mRawQx);
		// Double_t mRawEventPlane = 0.5*mRawQ->Phi();
		if(mRawEventPlane<0.) mRawEventPlane += TMath::Pi();
		hRawEventPlane->Fill(mRawEventPlane);
	}

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
	double mPlusQx_rejectE 	= mEtaPlusQx_rejectE;
	double mPlusQy_rejectE 	= mEtaPlusQy_rejectE;
	double mMinusQx_rejectE = mEtaMinusQx_rejectE;
	double mMinusQy_rejectE = mEtaMinusQy_rejectE;
	Float_t Qx = mEtaPlusQx_rejectE + mEtaMinusQx_rejectE; 
	Float_t Qy = mEtaPlusQy_rejectE + mEtaMinusQy_rejectE;
	TVector2 Q(Qx,Qy);
	TVector2 mRawQWest(mEtaPlusQx_rejectE,mEtaPlusQy_rejectE);
	TVector2 mRawQEast(mEtaMinusQx_rejectE,mEtaMinusQy_rejectE);
	if (mRawQEast.Mod() <= 0 || mRawQWest.Mod() <= 0 ) return kFALSE;
	
	// reCenter process
	Double_t mReCenterQx, mReCenterQy;
	Double_t mReCenterQx_rejectE, mReCenterQy_rejectE;
	Double_t mReCenterQx_east, mReCenterQy_east; 
	Double_t mReCenterQx_west, mReCenterQy_west; 
	// if(vz>0){
	// 	mReCenterQx = mRawQx - mEtaPlusNTrks*etapluszplusQx->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszplusQx->GetBinContent(runIndex+1, mCentrality);
	// 	mReCenterQy = mRawQy - mEtaPlusNTrks*etapluszplusQy->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszplusQy->GetBinContent(runIndex+1, mCentrality);
	// }
	// else{
	// 	mReCenterQx = mRawQx - mEtaPlusNTrks*etapluszminusQx->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszminusQx->GetBinContent(runIndex+1, mCentrality);
	// 	mReCenterQy = mRawQy - mEtaPlusNTrks*etapluszminusQy->GetBinContent(runIndex+1, mCentrality) - mEtaMinusNTrks*etaminuszminusQy->GetBinContent(runIndex+1, mCentrality);
	// }
	mPlusQx     = mPlusQx-etaplusQx_cent->GetBinContent(mCentrality+1);
	mPlusQy     = mPlusQy-etaplusQy_cent->GetBinContent(mCentrality+1);
	mMinusQx    = mMinusQx-etaminusQx_cent->GetBinContent(mCentrality+1);
	mMinusQy    = mMinusQy-etaminusQy_cent->GetBinContent(mCentrality+1);
	mPlusQx_rejectE     = mPlusQx_rejectE-etaplusQx_cent_rejectE->GetBinContent(mCentrality+1);
	mPlusQy_rejectE     = mPlusQy_rejectE-etaplusQy_cent_rejectE->GetBinContent(mCentrality+1);
	mMinusQx_rejectE    = mMinusQx_rejectE-etaminusQx_cent_rejectE->GetBinContent(mCentrality+1);
	mMinusQy_rejectE    = mMinusQy_rejectE-etaminusQy_cent_rejectE->GetBinContent(mCentrality+1);
	mReCenterQx = mPlusQx + mMinusQx;
	mReCenterQy = mPlusQy + mMinusQy;
	mReCenterQx_rejectE = mPlusQx_rejectE+mPlusQy_rejectE;
	mReCenterQy_rejectE = mPlusQy_rejectE+mMinusQy_rejectE;
	mReCenterQx_east = mMinusQx_rejectE;
	mReCenterQy_east = mMinusQy_rejectE;
	mReCenterQx_west = mPlusQx_rejectE;
	mReCenterQy_west = mPlusQy_rejectE;
	hRecenterQxQy[mCentrality]->Fill(mReCenterQx,mReCenterQy);

	// reCenter process
	Double_t mReCenterQx_run, mReCenterQy_run;
	if(vz>0){
		mReCenterQx_run = mRawQx - etapluszplusQx->GetBinContent(runIndex+1, mCentrality+1) - etaminuszplusQx->GetBinContent(runIndex+1, mCentrality+1);
		mReCenterQy_run = mRawQy - etapluszplusQy->GetBinContent(runIndex+1, mCentrality+1) - etaminuszplusQy->GetBinContent(runIndex+1, mCentrality+1);
	}
	else{
		mReCenterQx_run = mRawQx - etapluszminusQx->GetBinContent(runIndex+1, mCentrality+1) - etaminuszminusQx->GetBinContent(runIndex+1, mCentrality+1);
		mReCenterQy_run = mRawQy - etapluszminusQy->GetBinContent(runIndex+1, mCentrality+1) - etaminuszminusQy->GetBinContent(runIndex+1, mCentrality+1);
	}

    TVector2 *mReCenterQ = new TVector2(mReCenterQx, mReCenterQy);
	TVector2 *mReCenterQEast = new TVector2(mReCenterQx_east, mReCenterQy_east);
    TVector2 *mReCenterQWest = new TVector2(mReCenterQx_west, mReCenterQy_west);
    TVector2 *mReCenterQ_rejectE = new TVector2(mReCenterQx_rejectE, mReCenterQy_rejectE);
	Double_t mReCenterEventPlane = -1;
	Double_t mReCenterEventPlane_rejectE = -1;
	Double_t recenterEPEast;
	Double_t recenterEPWest;
	if(mReCenterQ->Mod() > 0){
		mReCenterEventPlane = 0.5*TMath::ATan2(mReCenterQy,mReCenterQx);
		// mReCenterEventPlane = 0.5*mReCenterQ->Phi();
		if(mReCenterEventPlane<0.) mReCenterEventPlane += TMath::Pi();
		hReCenterEventPlane->Fill(mReCenterEventPlane);
	}
	if(mReCenterQ->Mod() > 0){
		mReCenterEventPlane_rejectE = 0.5*TMath::ATan2(mReCenterQy,mReCenterQx);
		// mReCenterEventPlane = 0.5*mReCenterQ->Phi();
		if(mReCenterEventPlane_rejectE<0.) mReCenterEventPlane_rejectE += TMath::Pi();
		hReCenterEventPlane_rejectE->Fill(mReCenterEventPlane_rejectE);
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

	for(Int_t j=0; j<mArrayLength; j++){

			if (mEtaPlusNTrks>0)
			{
				shiftfactorcos_cent_west_rejectE[j]->Fill(mCentrality,cos(2*(j+1)*recenterEPWest));
				shiftfactorsin_cent_west_rejectE[j]->Fill(mCentrality,sin(2*(j+1)*recenterEPWest));
			}
			if(mEtaMinusNTrks>0)
			{
				shiftfactorcos_cent_east_rejectE[j]->Fill(mCentrality,cos(2*(j+1)*recenterEPEast));
				shiftfactorsin_cent_east_rejectE[j]->Fill(mCentrality,sin(2*(j+1)*recenterEPEast));
			}
			// if(mEtaPlusNTrks>0 &&  mEtaPlusNTrks > 0)
			// {
			// 	shiftfactorcos[j]->Fill(dayIndex,mCentrality,cos(2*(j+1)*mReCenterEventPlane));
			// 	shiftfactorsin[j]->Fill(dayIndex,mCentrality,sin(2*(j+1)*mReCenterEventPlane));
			// 	shiftfactorcos_cent[j]->Fill(mCentrality,cos(2*(j+1)*mReCenterEventPlane));
			// 	shiftfactorsin_cent[j]->Fill(mCentrality,sin(2*(j+1)*mReCenterEventPlane));
			// 	shiftfactorcos_cent_rejectE[j]->Fill(mCentrality,cos(2*(j+1)*mReCenterEventPlane_rejectE));
			// 	shiftfactorsin_cent_rejectE[j]->Fill(mCentrality,sin(2*(j+1)*mReCenterEventPlane_rejectE));
			// }
			shiftfactorcos[j]->Fill(dayIndex,mCentrality,cos(2*(j+1)*mReCenterEventPlane));
			shiftfactorsin[j]->Fill(dayIndex,mCentrality,sin(2*(j+1)*mReCenterEventPlane));
			shiftfactorcos_cent[j]->Fill(mCentrality,cos(2*(j+1)*mReCenterEventPlane));
			shiftfactorsin_cent[j]->Fill(mCentrality,sin(2*(j+1)*mReCenterEventPlane));
			shiftfactorcos_cent_rejectE[j]->Fill(mCentrality,cos(2*(j+1)*mReCenterEventPlane_rejectE));
			shiftfactorsin_cent_rejectE[j]->Fill(mCentrality,sin(2*(j+1)*mReCenterEventPlane_rejectE));
			
	}

	return kTRUE;
}
//____________________________________________________________
void bookHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.histo.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	fOutFile = new TFile(buf,"recreate");

	for(Int_t i=0;i<mArrayLength;i++){
		sprintf(buf,"shiftfactorcos_%d",i);
		shiftfactorcos[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorsin_%d",i);
		shiftfactorsin[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorcos_cent_%d",i);
		shiftfactorcos_cent[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorsin_cent_%d",i);
		shiftfactorsin_cent[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorcos_cent_rejectE_%d",i);
		shiftfactorcos_cent_rejectE[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorsin_cent_rejectE_%d",i);
		shiftfactorsin_cent_rejectE[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorcos_cent_east_rejectE_%d",i);
		shiftfactorcos_cent_east_rejectE[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorcos_cent_west_rejectE_%d",i);
		shiftfactorcos_cent_west_rejectE[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorsin_cent_east_rejectE_%d",i);
		shiftfactorsin_cent_east_rejectE[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"shiftfactorsin_cent_west_rejectE_%d",i);
		shiftfactorsin_cent_west_rejectE[i] = new TProfile(buf,buf,mTotalCentrality,0,mTotalCentrality);
	}
	hReCenterEventPlaneEast = new TH1D("hReCenterEventPlaneEast","hReCenterEventPlaneEast;Reaction Plane East (rad); Counts",300,0,TMath::Pi());
	hReCenterEventPlaneWest = new TH1D("hReCenterEventPlaneWest","hReCenterEventPlaneWest;Reaction Plane West (rad); Counts",300,0,TMath::Pi());

	for(int i = 0; i < nCent; i++)
	{
		hRawQxQy[i] = new TH2D(Form("hRawQxQy_icent%d",i),Form("hRawQxQy_icent%d; Qx, Qy",i),400,-10,10,400,-10,10);
		hRecenterQxQy[i] = new TH2D(Form("hRecenterQxQy_icent%d",i),Form("hRecenterQxQy_icent%d; Qx, Qy",i),400,-10,10,400,-10,10);
	}

	hRawEventPlane = new TH1D("hRawEventPlane","hRawEventPlane,Event Plane",360,0,TMath::Pi());
	hReCenterEventPlane = new TH1D("hReCenterEventPlane","hReCenterEventPlane;Event Plane",360,0,TMath::Pi());
	hReCenterEventPlane_rejectE = new TH1D("hReCenterEventPlane_rejectE","hReCenterEventPlane_rejectE;Event Plane",360,0,TMath::Pi());
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
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();
	for(map<Int_t,Int_t>::iterator iter=mTotalRunId.begin();iter!=mTotalRunId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;
	cout<<endl;

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
	for(map<Int_t,Int_t>::iterator iter=mTotalDayId.begin();iter!=mTotalDayId.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;
	cout<<endl;


	//read in bad run for 580001 and 580021
	ifstream indata_001;
	indata_001.open("/star/u/wangzhen/run20/Dielectron/BadRunList/BadRunList.dat");
	mBadRunId_001.clear();
	if(indata_001.is_open()){
		cout<<"read in total bad run number list and recode bad run number ...";
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
	
  cout<<"bad run for trigger Au+Au 9.2 GeV"<<endl;
	for(map<Int_t,Int_t>::iterator iter=mBadRunId_001.begin();iter!=mBadRunId_001.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;



    fReCenter = TFile::Open("/star/u/wangzhen/run20/Dielectron/FlatEvtPlane/reCenter/output_all/reCenter.root");
	etapluszplusQx   = (TProfile2D *)fReCenter->Get("etapluszplusQx");
	etapluszminusQx  = (TProfile2D *)fReCenter->Get("etapluszminusQx");
	etaminuszplusQx  = (TProfile2D *)fReCenter->Get("etaminuszplusQx");
	etaminuszminusQx = (TProfile2D *)fReCenter->Get("etaminuszminusQx");
	etapluszplusQy   = (TProfile2D *)fReCenter->Get("etapluszplusQy");
	etapluszminusQy  = (TProfile2D *)fReCenter->Get("etapluszminusQy");
	etaminuszplusQy  = (TProfile2D *)fReCenter->Get("etaminuszplusQy");
	etaminuszminusQy = (TProfile2D *)fReCenter->Get("etaminuszminusQy");
	etaplusQx_cent   = (TProfile*)fReCenter->Get("etaplusQx_cent");
	etaminusQx_cent  = (TProfile*)fReCenter->Get("etaminusQx_cent");
	etaplusQy_cent   = (TProfile*)fReCenter->Get("etaplusQy_cent");
	etaminusQy_cent  = (TProfile*)fReCenter->Get("etaminusQy_cent");
	etaplusQx_cent_rejectE   = (TProfile*)fReCenter->Get("etaplusQx_cent_RejectE");
	etaminusQx_cent_rejectE  = (TProfile*)fReCenter->Get("etaminusQx_cent_RejectE");
	etaplusQy_cent_rejectE   = (TProfile*)fReCenter->Get("etaplusQy_cent_RejectE");
	etaminusQy_cent_rejectE  = (TProfile*)fReCenter->Get("etaminusQy_cent_RejectE");

	cout<<"Initialization DONE !!!"<<endl;
	cout<<endl;

	return kTRUE;
}

