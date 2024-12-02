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
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TH3D.h"
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
#include "TAxis.h"

#include "miniDst.h"
#include "StRefMultCorr.h"
#include "cuts.h"

using namespace std;
#endif


StRefMultCorr *refMultCorrUtil;

Int_t runIndex;
Int_t randomId;
map<Int_t,Int_t> mTotalRunId;
map<TString, TH2*> mHisto_2D;
map<TString, TH1*> mHisto_1D;

bool Init();
void bookHistograms();
bool passEvent(miniDst* event);
bool passTrack(miniDst* event, Int_t i);
void writeHistograms(char* outFile);

//define histograms


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

	// hEvent->GetXaxis()->SetBinLabel(1, "All events");
	// hEvent->GetXaxis()->SetBinLabel(3, "minbias");
	// hEvent->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
	// hEvent->GetXaxis()->SetBinLabel(9,  Form("|V_{r}|<%1.2f cm",10));
	// hEvent->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm",100));
	// hEvent->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm",13));
	//+---------------------------------+
	//| open files and add to the chain |
	//+---------------------------------+
	TFile* inputFile = new TFile("/star/u/wangzhen/run20/Dielectron_Common/minitree/output/9p2_Phiweight/DD0E67C29473FCB5A275A2F2CDC9A2A0_602.root");
	//TH1
	mHisto_1D["hEvent"]                  = (TH1D*)inputFile->Get("hEvent");
	mHisto_1D["hVzDiff"]                 = (TH1D*)inputFile->Get("hVzDiff");
	mHisto_1D["hCentrality"]             = (TH1D*)inputFile->Get("hCentrality");
	mHisto_1D["hRawQx"]                  = (TH1D*)inputFile->Get("hRawQx");
	mHisto_1D["hRawQy"]                  = (TH1D*)inputFile->Get("hRawQy");
	mHisto_1D["hRawEventPlane"]          = (TH1D*)inputFile->Get("hRawEventPlane");
	mHisto_1D["hVz"]                     = (TH1D*)inputFile->Get("hVz");
	mHisto_1D["hNHitsFit"]               = (TH1D*)inputFile->Get("hNHitsFit");
	mHisto_1D["hNHitsPoss"]              = (TH1D*)inputFile->Get("hNHitsPoss");
	mHisto_1D["hNHitsdEdx"]              = (TH1D*)inputFile->Get("hNHitsdEdx");

	//TH2 
	mHisto_2D["hVtxYvsVtxX"]             = (TH2D*)inputFile->Get("hVtxYvsVtxX");
	mHisto_2D["hVPDVzvsTPCVz"]           = (TH2D*)inputFile->Get("hVPDVzvsTPCVz");
	mHisto_2D["hGRefMultvsGRefMultCorr"] = (TH2D*)inputFile->Get("hGRefMultvsGRefMultCorr");
	mHisto_2D["hdEdxvsP"]                = (TH2D*)inputFile->Get("hdEdxvsP");
	mHisto_2D["hdNdxvsP"]                = (TH2D*)inputFile->Get("hdNdxvsP");
	mHisto_2D["hnSigEvsP"]               = (TH2D*)inputFile->Get("hnSigEvsP");
	mHisto_2D["hBetavsP"]                = (TH2D*)inputFile->Get("hBetavsP");
	mHisto_2D["hDcavsPt"]                = (TH2D*)inputFile->Get("hDcavsPt");
	mHisto_2D["hBetavsEta"]              = (TH2D*)inputFile->Get("hBetavsEta");
	mHisto_2D["hBetavsPhi"]              = (TH2D*)inputFile->Get("hBetavsPhi");
	mHisto_2D["dEdxvsEta"]               = (TH2D*)inputFile->Get("dEdxvsEta");
	mHisto_2D["dEdxvsPhi"]               = (TH2D*)inputFile->Get("dEdxvsPhi");
	mHisto_2D["hNSigmaEvsEta"]           = (TH2D*)inputFile->Get("hNSigmaEvsEta");
	mHisto_2D["hNSigmaEvsPhi"]           = (TH2D*)inputFile->Get("hNSigmaEvsPhi");
	mHisto_2D["hMSquarevsP"]             = (TH2D*)inputFile->Get("hMSquarevsP");
	mHisto_2D["hEtavsPhi"]               = (TH2D*)inputFile->Get("hEtavsPhi");
	mHisto_2D["hEtavsPt"]                = (TH2D*)inputFile->Get("hEtavsPt");
	mHisto_2D["hPhivsPt"]                = (TH2D*)inputFile->Get("hPhivsPt");
	mHisto_2D["hEEtavsPhi"]              = (TH2D*)inputFile->Get("hEEtavsPhi");
	mHisto_2D["hEEtavsPt"]               = (TH2D*)inputFile->Get("hEEtavsPt");
	mHisto_2D["hEPhivsPt"]               = (TH2D*)inputFile->Get("hEPhivsPt");
	mHisto_2D["hnTOFMatchvsRefmult"]     = (TH2D*)inputFile->Get("hnTOFMatchvsRefmult");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent0"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent0");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent1"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent1");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent2"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent2");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent3"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent3");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent4"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent4");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent5"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent5");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent6"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent6");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent7"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent7");
	mHisto_2D["hPrimaryTrackPhiVsEta_Cent8"]     = (TH2D*)inputFile->Get("hPrimaryTrackPhiVsEta_Cent8");

	//try to load the histograms without any hard code
	//you can reference that https://root-forum.cern.ch/t/accessing-all-histograms-inside-a-rootfile/36070/24
	//https://root-forum.cern.ch/t/reading-all-histograms-from-a-root-file/11607/2

	for( auto his : mHisto_1D )
	{
		his.second->Reset("ICES");
	}
	for( auto his : mHisto_2D)
	{
		his.second->Reset("ICES");
	}

	vector <TString>histo_name1D = { "hEvent", "hVzDiff", "hCentrality", "hRawQx", "hRawQy", "hRawEventPlane", "hVz", "hNHitsFit", "hNHitsPoss", "hNHitsdEdx"};
	vector <TString>histo_name2D = {"hVtxYvsVtxX", "hVPDVzvsTPCVz", "hGRefMultvsGRefMultCorr", "hdEdxvsP", "hdNdxvsP", "hnSigEvsP", "hBetavsP", "hDcavsPt", "hBetavsEta", "hBetavsPhi", "dEdxvsEta", "dEdxvsPhi", "hNSigmaEvsEta", "hNSigmaEvsPhi", "hMSquarevsP", "hEtavsPhi", "hEtavsPt", "hPhivsPt", "hEEtavsPhi", "hEEtavsPt", "hEPhivsPt", "hnTOFMatchvsRefmult", "hPrimaryTrackPhiVsEta_Cent0", "hPrimaryTrackPhiVsEta_Cent1", "hPrimaryTrackPhiVsEta_Cent2", "hPrimaryTrackPhiVsEta_Cent3", "hPrimaryTrackPhiVsEta_Cent4", "hPrimaryTrackPhiVsEta_Cent5", "hPrimaryTrackPhiVsEta_Cent6", "hPrimaryTrackPhiVsEta_Cent7", "hPrimaryTrackPhiVsEta_Cent8"};
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
				for(auto his : histo_name1D)
				{
					TH1D* histem = (TH1D*)ftmp->Get(his);
					mHisto_1D[his]->Add(histem);
				}
				for(auto his : histo_name2D)
				{
					TH2D* histem = (TH2D*)ftmp->Get(his);
					mHisto_2D[his]->Add(histem);
				}

				ifile++;
			}
			delete ftmp;
		}
	}
	delete inputStream;
        //ofstream out;
        //out.open("./outlog.txt");
	//intialization



	char buf[1024];
	sprintf(buf,"%s.QAhisto.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	TFile *mFile = new TFile(buf,"recreate");
	mFile->cd();
	for(auto his : mHisto_2D)
	{
		his.second->Write();
	}
	for(auto his : mHisto_1D)
	{
		his.second->Write();
	}
    //     //out.close();
	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________

