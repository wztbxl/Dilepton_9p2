#include "/star/u/wangzhen/run20/Dielectron/func_macro/headers.h"
#include "/star/u/wangzhen/run20/Dielectron/func_macro/function.C"
#include "/star/u/wangzhen/run20/Dielectron/Analysis/cuts.h"
#include "/star/u/wangzhen/run20/Dielectron/Analysis/RefMfun.h"

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
TH1F *hRefMult;
TH1F *hRefMultCorZ;
TH1F *hRefMultCor;
TH1F *hCentrality;
TH1F *hCentralityCor;
TH2D *hnTofHitsvsRefMult;
TH2D* hnTofHitsvsRefMult_noCut;
TH2D* hnTofHitsvsRefMult_Vz35;
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
	Float_t vx           = event->mVertexX;
	Float_t vy           = event->mVertexY;
	Float_t vz           = event->mVertexZ;
	Float_t vpdVz        = event->mVpdVz;
	Float_t vr           = sqrt(vx*vx + vy*vy);    
	Float_t vzDiff       = vz - vpdVz;

    Double_t RefMultCorr = refMult;
	// if(RefMVzCorFlag)RefMultCorr = GetRefMultCorr(refMult, vz);	
    Double_t reweight = 1.;// no RefMultCorr!!!!	
    // Double_t reweight = GetWeight(RefMultCorr);	
	mCentrality = GetCentrality(RefMultCorr);


	Float_t  mEtaPlusQx        = event->mEtaPlusQx;
	Float_t  mEtaPlusQy        = event->mEtaPlusQy;
	//Float_t  mEtaPlusPtWeight  = event->mEtaPlusPtWeight;
	Short_t  mEtaPlusNTrks     = event->mEtaPlusNTrks;
	Float_t  mEtaMinusQx       = event->mEtaMinusQx;
	Float_t  mEtaMinusQy       = event->mEtaMinusQy;
	//Float_t  mEtaMinusPtWeight = event->mEtaMinusPtWeight;
	Short_t  mEtaMinusNTrks    = event->mEtaMinusNTrks;

  hnTofHitsvsRefMult_noCut->Fill(refMult,mnTOFMatch);
	if(vr>mVrCut)                     return kFALSE;
	if(TMath::Abs(vz)>mVzCut)         return kFALSE;
	// if(TMath::Abs(vzDiff)>mVzDiffCut) return kFALSE;
	if (mnTOFMatch < Pileuplimit->Eval(refMult)) return kFALSE;
  hnTofHitsvsRefMult->Fill(refMult,mnTOFMatch);
  if(mCentrality<0)                 return kFALSE;
	hRefMult->Fill(refMult);
	hRefMultCorZ->Fill(RefMultCorr);
	hRefMultCor->Fill(RefMultCorr, reweight);
	hCentrality->Fill(mCentrality);
	hCentralityCor->Fill(mCentrality,reweight);

	if(vz>0){
		if(mEtaPlusNTrks>0){
			etapluszplusQx->Fill(runIndex, mCentrality, mEtaPlusQx/mEtaPlusNTrks, mEtaPlusNTrks);
			etapluszplusQy->Fill(runIndex, mCentrality, mEtaPlusQy/mEtaPlusNTrks, mEtaPlusNTrks);
		}

		if(mEtaMinusNTrks>0){
			etaminuszplusQx->Fill(runIndex, mCentrality, mEtaMinusQx/mEtaMinusNTrks, mEtaMinusNTrks);
			etaminuszplusQy->Fill(runIndex, mCentrality, mEtaMinusQy/mEtaMinusNTrks, mEtaMinusNTrks);
		}
	}
	else{
		if(mEtaPlusNTrks>0){
			etapluszminusQx->Fill(runIndex, mCentrality, mEtaPlusQx/mEtaPlusNTrks, mEtaPlusNTrks);
			etapluszminusQy->Fill(runIndex, mCentrality, mEtaPlusQy/mEtaPlusNTrks, mEtaPlusNTrks);
		}

		if(mEtaMinusNTrks>0){
			etaminuszminusQx->Fill(runIndex, mCentrality, mEtaMinusQx/mEtaMinusNTrks, mEtaMinusNTrks);
			etaminuszminusQy->Fill(runIndex, mCentrality, mEtaMinusQy/mEtaMinusNTrks, mEtaMinusNTrks);
		}
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
	hRefMult  = new TH1F("RefMult","RefMult;RefMult",600,0,600);
	hRefMultCorZ  = new TH1F("RefMultCorZ","RefMult corrected for z dependence ;RefMult",600,0,600);
	hRefMultCor  = new TH1F("RefMultCor","RefMult corrected for z and weight ;Mult",600,0,600);
	hCentrality = new TH1F("Centrality","Centrality; centrality",16,0,16);
	hCentralityCor = new TH1F("CentralityCor","Centrality with weight; centrality",16,0,16);
  hnTofHitsvsRefMult = new TH2D("hnTofHitsvsRefMult",";RefMult;nTofHits",500,0,500,500,0,500);
  hnTofHitsvsRefMult_noCut = new TH2D("hnTofHitsvsRefMult_noCut",";RefMult;nTofHits",500,0,500,500,0,500); 
  hnTofHitsvsRefMult_Vz35 = new TH2D("hnTofHitsvsRefMult_Vz35",";RefMult;nTofHits",500,0,500,500,0,500);

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
