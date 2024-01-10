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

#include "MCEVENT.h"
// #include "StRefMultCorr/StRefMultCorr.h"
// #include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr.h"
#include "CentralityMaker.h"
#include "cuts.h"
//#include "RefMfun.h"

using namespace std;
#endif

StRefMultCorr *refMultCorrUtil;
Int_t mCentrality;

const Double_t PI = TMath::Pi();

void bookHistograms();
void writeHistograms(char* outFile);
bool Init();

//***** constrain the bad dedx calibration TPC sector geometry *****
TF1 *funPosHi;
TF1 *funPosLow;
TF1 *funNegHi;
TF1 *funNegLow;
Float_t par[4][4];
Float_t parErr[4][4];

//define histograms
TH2F *hVyvsVx;
TH2F *hVyvsVz;
TH1F *hVzDiff;
TH2F *hRcMcVxDiffvsMcVx;
TH2F *hRcMcVyDiffvsMcVy;
TH2F *hRcMcVzDiffvsMcVz;
TH2F *hRefMultvsRefMultCorr;
TH2F *hNMatchTrksvsInputTrks;

TH3F *hPtResvsPtCen;
TH3F *hMcEtavsPtQ;
TH3F *hMcPhivsPtQ;
TH3F *hMcYvsPtQ;
TH3F *hRcEtavsPtQ;
TH3F *hRcPhivsPtQ;
TH3F *hRcYvsPtQ;
TH3F *hRcPtvsMcPtQ;
TH3F *hRcEtavsMcEtaQ;
TH3F *hRcPhivsMcPhiQ;
TH3F *hRcYvsMcYQ;
TH3F *hPtResvsPtQ;
TH3F *hPResvsPQ;
TH3F *hRcDcavsPtQ;
TH3F *hRcNHitsFitvsPtQ;
TH3F *hRcNHitsFitvsEtaQ;
TH3F *hRcNHitsFitvsPhiQ;
TH3F *hRcNHitsPossvsPtQ;
TH3F *hRcNHitsPossvsEtaQ;
TH3F *hRcNHitsPossvsPhiQ;
TH3F *hRcNHitsDedxvsPtQ;
TH3F *hRcNHitsCommonvsPtQ;
TH3F *hRcNHitsCommonvsEtaQ;
TH3F *hRcNHitsCommonvsPhiQ;
TH3F *hRcDedxvsPQ;
TH3F *hRcDedxvsEtaQ;
TH3F *hRcDedxvsPhiQ;
TH3F *hRcNSigmaEvsPQ;
TH3F *hRcNSigmaEvsEtaQ;
TH3F *hRcNSigmaEvsPhiQ;
TH3F *hRcNSigmaPivsPQ;
TH3F *hRcNSigmaKvsPQ;
TH3F *hRcNSigmaPvsPQ;

//Zhen add it for 3D check
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

//TPC Tracking Efficiency
TH3F *hDenEPlusTpcEff;
TH3F *hNumEPlusTpcEff;
TH3F *hDenEMinusTpcEff;
TH3F *hNumEMinusTpcEff;

TH3F *hDenEPlusTpcEffCen[9];
TH3F *hNumEPlusTpcEffCen[9];
TH3F *hDenEMinusTpcEffCen[9];
TH3F *hNumEMinusTpcEffCen[9];
//check origin phi and Z
TH1F *hMinusOriginZ;
TH1F *hMinusOriginPhi;
TH1F *hPlusOriginZ;
TH1F *hPlusOriginPhi;
//Check Centrality 
TH1F *hCentrality9;
TH1F *hRefMult;
TH1F* hRefMultCorr;
//for pTEtaWeight
TH2D* pTEtaWeight[2];


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
	TChain *chain = new TChain("mcT");

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

	refMultCorrUtil = new StRefMultCorr();

	//+-------------+
	//| loop events |
	//+-------------+
	MCEVENT *event = new MCEVENT(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;
	for(int i=0;i<nEvts;i++){

		//initialize the struct
		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;
		event->GetEntry(i);

		Float_t mcVertexX = event->muPriVertexX;
		Float_t mcVertexY = event->muPriVertexY;
		Float_t mcVertexZ = event->muPriVertexZ;
		Float_t rcVertexX = event->muPriVertexX;
		Float_t rcVertexY = event->muPriVertexY;
		Float_t rcVertexZ = event->muPriVertexZ;
		Float_t rcVertexR = sqrt(rcVertexX*rcVertexX+rcVertexY*rcVertexY);
		Float_t rcVpdVz = event->muVpdVz; 
		Float_t vzDiff = rcVertexZ - rcVpdVz;
		Float_t rcRefMult = event->muRefMult;
		Float_t rcRefMultCorr = event->refmult_corr;
		Int_t triggerID ;
		if(event->triggerId[0] > 0) triggerID = event->triggerId[0];
		if(event->triggerId[1] > 0) triggerID = event->triggerId[1];
		if(event->triggerId[2] > 0) triggerID = event->triggerId[2];
		if(event->triggerId[3] > 0) triggerID = event->triggerId[3];
		if (triggerID != 810010 && triggerID != 810020 && triggerID != 810030 && triggerID != 810040) continue;

		hVyvsVx->Fill(rcVertexX,rcVertexY);
		hVyvsVz->Fill(rcVertexZ,rcVertexY);
		hVzDiff->Fill(vzDiff);
		hRcMcVxDiffvsMcVx->Fill(mcVertexX, rcVertexX-mcVertexX);
		hRcMcVyDiffvsMcVy->Fill(mcVertexY, rcVertexY-mcVertexY);
		hRcMcVzDiffvsMcVz->Fill(mcVertexZ, rcVertexZ-mcVertexZ);
		hRefMultvsRefMultCorr->Fill(rcRefMultCorr,rcRefMult);

		
		if( TMath::Abs(rcVertexX)<1.e-5 
				&& TMath::Abs(rcVertexY)<1.e-5
				&& TMath::Abs(rcVertexZ)<1.e-5 ) continue;
		if(rcVertexR>=mVrCut) continue;
		if(TMath::Abs(rcVertexZ)>=mVzCut) continue;
		if(TMath::Abs(vzDiff)>=mVzDiffCut) continue;

		Int_t runId = event->runId;
		Int_t zdcRate = event->zdcX;
		Int_t refMult = event->muRefMult;
		refMultCorrUtil->init(runId);
		refMultCorrUtil->initEvent(refMult,rcVertexZ,zdcRate);
		//mCentrality = refMultCorrUtil->getCentralityBin16();
		hRefMult->Fill(refMult);
		Double_t RefMultCorr = refMultCorrUtil->getRefMultCorr();
    	mCentrality = refMultCorrUtil->getCentralityBin9();//Centrality defined by offical, 0 is 70-80%,8 is 0-5%,0 is refMult<7 
    	Double_t weight = refMultCorrUtil->getWeight();
    // mCentrality = mCentrality + 1;
		hRefMultCorr->Fill(RefMultCorr);
		hCentrality9->Fill(mCentrality,weight);
		// mCentrality = mCentrality - 1;
		printf("Centrality is %d \n",mCentrality);
		//cout<<refMult<<"   "<<rcRefMultCorr<<"   read:"<<refMultCorrUtil->getRefMultCorr()<<endl;

		if(mCentrality<0 || mCentrality>8) continue;

		Int_t nMatchedE = 0;
		for(Int_t i=0;i<event->nMcE;i++){
			Int_t geantId = event->geantId[i];
			Float_t q;
			if(geantId==2) q=0.5;
			else if(geantId==3) q=-0.5;
			else cout<<"The geantId is wrong !"<<endl;
			Float_t mcPt = event->mcPt[i];
			Float_t mcEta = event->mcEta[i];
			Float_t mcPhi = event->mcPhi[i];
			TLorentzVector mcFourMom(0.,0.,0.,0.);
			mcFourMom.SetPtEtaPhiM(mcPt,mcEta,mcPhi,Melectron);
			Float_t mcY = mcFourMom.Rapidity();
			Float_t mcP = mcFourMom.P();
			Float_t rcPt = event->rcPt[i];
			Float_t rcEta = event->rcEta[i];
			Float_t rcPhi = event->rcPhi[i];
			TLorentzVector rcFourMom(0.,0.,0.,0.);
			rcFourMom.SetPtEtaPhiM(rcPt,rcEta,rcPhi,Melectron);
			Float_t rcY = rcFourMom.Rapidity();
			Float_t rcP = rcFourMom.P();
			Int_t rcNHitsFit = event->rcNHitsFit[i];
			Int_t rcNHitsPoss = event->rcNHitsPoss[i];
			Int_t rcNHitsDedx = event->rcNHitsDedx[i];
			// Int_t rcNHitsCommon = event->rcNHitsCommon[i];
			Float_t rcDedx = event->rcDedx[i];
			Float_t rcNSigmaE = event->rcNSigmaE[i];
			Float_t rcNSigmaPi = event->rcNSigmaPi[i];
			Float_t rcNSigmaK = event->rcNSigmaK[i];
			Float_t rcNSigmaP = event->rcNSigmaP[i];
			Float_t rcDca = event->rcDca[i];
			Float_t rcPhiFirst = event->rcPhi[i];

			hMcEtavsPtQ->Fill(q,mcPt,mcEta);
			hMcPhivsPtQ->Fill(q,mcPt,mcPhi);
			hMcYvsPtQ->Fill(q,mcPt,mcY);

			/*if(q>0. && mcEta>0. 
					&& mcPhi<=funPosHi->Eval(mcPt)
					&& mcPhi>=funPosLow->Eval(mcPt)
			  ){
				continue;
			}else if(q<0. && mcEta>0.
					&& mcPhi<=funNegHi->Eval(-1*mcPt)
					&& mcPhi>=funNegLow->Eval(-1*mcPt)
					){
				continue;
			}**/

			if(mcPhi<-PI) mcPhi += 2*PI;

			Float_t ratio = 0.;
			if(rcDedx>0) ratio = rcNHitsFit*1./rcNHitsPoss;

			// if(TMath::Abs(mcEta)<=1.0){
			// 	if(q>0.){
			// 	   	hDenEPlusTpcEff->Fill(mcPt,mcEta,mcPhi);
			// 		hDenEPlusTpcEffCen[mCentrality]->Fill(mcPt,mcEta,mcPhi);
			// 	}else{
			// 	   	hDenEMinusTpcEff->Fill(mcPt,mcEta,mcPhi); 
			// 		hDenEMinusTpcEffCen[mCentrality]->Fill(mcPt,mcEta,mcPhi);
			// 	}

			// 	if(rcNHitsFit>=mTpceNHitsFitCut//for systemiac uncentraity change this +/-5
			// 	//if(rcNHitsFit>=20//for systemiac uncentraity change this +/-5
			// 			&& ratio>=mTpceNHitsFitRatioCut
			// 			&& rcDca<=mTpceDcaCut
			// 			// && rcDca<=0.8
			// 			&& rcNHitsDedx>=15){
			// 			// && rcNHitsDedx>=20){
			// 		if(q>0){
			// 		   	hNumEPlusTpcEff->Fill(mcPt,mcEta,mcPhi);
			// 			hNumEPlusTpcEffCen[mCentrality]->Fill(mcPt,mcEta,mcPhi);
			// 		}else{
			// 		   	hNumEMinusTpcEff->Fill(mcPt,mcEta,mcPhi);
			// 			hNumEMinusTpcEffCen[mCentrality]->Fill(mcPt,mcEta,mcPhi);
			// 		}
			// 	}
			// }
			if(TMath::Abs(mcEta)<=1.0){
				if(q>0.){
				   	hDenEPlusTpcEff->Fill(mcPt,mcEta,mcPhi);
					hDenEPlusTpcEffCen[mCentrality]->Fill(mcPt,mcEta,mcPhi,weight);
				}else{
				   	hDenEMinusTpcEff->Fill(mcPt,mcEta,mcPhi); 
					hDenEMinusTpcEffCen[mCentrality]->Fill(mcPt,mcEta,mcPhi,weight);
				}

				if(rcNHitsFit>=mTpceNHitsFitCut//for systemiac uncentraity change this +/-5, origin is 20
				 // if(rcNHitsFit>=1//for systemiac uncentraity change this +/-5
						&& ratio>=mTpceNHitsFitRatioCut
						&& rcDca<=mTpceDcaCut // for sys. uncent. is 1.2, origin is 1
						// && rcDca<=3
						// ){
						&& rcNHitsDedx>mTpceNHitsDedxCut){ //for sys. uncent. is 20, origin is 15
						//&& rcNHitsDedx>=1){
					if(q>0){
					   	hNumEPlusTpcEff->Fill(rcPt,rcEta,rcPhi);
						hNumEPlusTpcEffCen[mCentrality]->Fill(rcPt,rcEta,rcPhi,weight);
					}else{
					   	hNumEMinusTpcEff->Fill(rcPt,rcEta,rcPhi);
						hNumEMinusTpcEffCen[mCentrality]->Fill(rcPt,rcEta,rcPhi,weight);
					}
				}
			}

			if(mcPhi>PI) mcPhi -= 2*PI;

			//The cuts below are applied in my pico production macro
			if(rcPt<0.2) continue;
			if(rcNHitsFit<10) continue;
			if(ratio<0.52) continue;
			if(rcNHitsDedx<10) continue;
			if(rcDca>3) continue;
			if(TMath::Abs(rcEta)> 1.3) continue;

			nMatchedE++;

			int pTEtaBinMinus = pTEtaWeight[0]->FindBin(rcPt,rcEta);
			int pTEtaBinPlus = pTEtaWeight[1]->FindBin(rcPt,rcEta);
			// double weightMinus = pTEtaWeight[0]->GetBinContent(pTEtaBinMinus);
			// double weightPlus = pTEtaWeight[1]->GetBinContent(pTEtaBinPlus);
			double weightMinus = 1; // this version do not add the weight, will check the difference with and without weight
			double weightPlus = 1;

    	hPtResvsPtCen->Fill(mCentrality,mcPt,(rcPt-mcPt)/mcPt,weight);

			if (q > 0)
			{
				hRcEtavsPtQ->Fill(q,rcPt,rcEta,weightPlus);
				hRcPhivsPtQ->Fill(q,rcPt,rcPhi,weightPlus);
				hRcYvsPtQ->Fill(q,rcPt,rcY,weightPlus);
				hRcPtvsMcPtQ->Fill(q,mcPt,rcPt,weightPlus);
				hRcEtavsMcEtaQ->Fill(q,mcEta,rcEta,weightPlus);
				hRcPhivsMcPhiQ->Fill(q,mcPhi,rcPhi,weightPlus);
				hRcYvsMcYQ->Fill(q,mcY,rcY,weightPlus);
				hPtResvsPtQ->Fill(q,mcPt,(rcPt-mcPt)/mcPt,weightPlus);
				hPResvsPQ->Fill(q,mcP,(rcP-mcP)/mcP,weightPlus);

				hRcDcavsPtQ->Fill(q,rcPt,rcDca,weightPlus);
				hRcNHitsFitvsPtQ->Fill(q,rcPt,rcNHitsFit,weightPlus);
				hRcNHitsFitvsEtaQ->Fill(q,rcEta,rcNHitsFit,weightPlus);
				hRcNHitsFitvsPhiQ->Fill(q,rcPhi,rcNHitsFit,weightPlus);
				hRcNHitsPossvsPtQ->Fill(q,rcPt,rcNHitsPoss,weightPlus);
				hRcNHitsPossvsEtaQ->Fill(q,rcEta,rcNHitsPoss,weightPlus);
				hRcNHitsPossvsPhiQ->Fill(q,rcPhi,rcNHitsPoss,weightPlus);
				hRcNHitsDedxvsPtQ->Fill(q,rcPt,rcNHitsDedx,weightPlus);
				// hRcNHitsCommonvsPtQ->Fill(q,rcPt,rcNHitsCommon,weightPlus);
				// hRcNHitsCommonvsEtaQ->Fill(q,rcEta,rcNHitsCommon,weightPlus);
				// hRcNHitsCommonvsPhiQ->Fill(q,rcPhi,rcNHitsCommon,weightPlus);
			}
			if (q < 0)
			{
				hRcEtavsPtQ->Fill(q,rcPt,rcEta,weightMinus);
				hRcPhivsPtQ->Fill(q,rcPt,rcPhi,weightMinus);
				hRcYvsPtQ->Fill(q,rcPt,rcY,weightMinus);
				hRcPtvsMcPtQ->Fill(q,mcPt,rcPt,weightMinus);
				hRcEtavsMcEtaQ->Fill(q,mcEta,rcEta,weightMinus);
				hRcPhivsMcPhiQ->Fill(q,mcPhi,rcPhi,weightMinus);
				hRcYvsMcYQ->Fill(q,mcY,rcY,weightMinus);
				hPtResvsPtQ->Fill(q,mcPt,(rcPt-mcPt)/mcPt,weightMinus);
				hPResvsPQ->Fill(q,mcP,(rcP-mcP)/mcP,weightMinus);

				hRcDcavsPtQ->Fill(q,rcPt,rcDca,weightMinus);
				hRcNHitsFitvsPtQ->Fill(q,rcPt,rcNHitsFit,weightMinus);
				hRcNHitsFitvsEtaQ->Fill(q,rcEta,rcNHitsFit,weightMinus);
				hRcNHitsFitvsPhiQ->Fill(q,rcPhi,rcNHitsFit,weightMinus);
				hRcNHitsPossvsPtQ->Fill(q,rcPt,rcNHitsPoss,weightMinus);
				hRcNHitsPossvsEtaQ->Fill(q,rcEta,rcNHitsPoss,weightMinus);
				hRcNHitsPossvsPhiQ->Fill(q,rcPhi,rcNHitsPoss,weightMinus);
				hRcNHitsDedxvsPtQ->Fill(q,rcPt,rcNHitsDedx,weightMinus);
				// hRcNHitsCommonvsPtQ->Fill(q,rcPt,rcNHitsCommon,weightMinus);
				// hRcNHitsCommonvsEtaQ->Fill(q,rcEta,rcNHitsCommon,weightMinus);
				// hRcNHitsCommonvsPhiQ->Fill(q,rcPhi,rcNHitsCommon,weightMinus);
			}
			
			if (q<0)
			{
				hPENHitsFitsvsPtvsEtaMinus->Fill(rcPt,rcEta,rcNHitsFit,weightMinus);
				hPENHitsFitsvsPtvsPhiMinus->Fill(rcPt,rcPhi,rcNHitsFit,weightMinus);
				hPENHitsDedxvsPtvsEtaMinus->Fill(rcPt,rcEta,rcNHitsDedx,weightMinus);
				hPENHitsDedxvsPtvsPhiMinus->Fill(rcPt,rcPhi,rcNHitsDedx,weightMinus);
				hPEDcavsPtvsEtaMinus->Fill(rcPt,rcEta,rcDca,weightMinus);
				hPEDcavsPtvsPhiMinus->Fill(rcPt,rcPhi,rcDca,weightMinus);
				hMinusOriginPhi->Fill(rcPhiFirst,weightMinus);
			}
			if (q>0)
			{
				hPENHitsFitsvsPtvsEtaPlus->Fill(rcPt,rcEta,rcNHitsFit,weightPlus);
            	hPENHitsFitsvsPtvsPhiPlus->Fill(rcPt,rcPhi,rcNHitsFit,weightPlus);
				hPENHitsDedxvsPtvsEtaPlus->Fill(rcPt,rcEta,rcNHitsDedx,weightPlus);
				hPENHitsDedxvsPtvsPhiPlus->Fill(rcPt,rcPhi,rcNHitsDedx,weightPlus);
				hPEDcavsPtvsEtaPlus->Fill(rcPt,rcEta,rcDca,weightPlus);
				hPEDcavsPtvsPhiPlus->Fill(rcPt,rcPhi,rcDca,weightPlus);
				hPlusOriginPhi->Fill(rcPhiFirst,weightPlus);
			}

			if ( q > 0)
			{
				hRcDedxvsPQ->Fill(q,rcP,rcDedx,weightPlus);
				hRcDedxvsEtaQ->Fill(q,rcEta,rcDedx,weightPlus);
				hRcDedxvsPhiQ->Fill(q,rcPhi,rcDedx,weightPlus);
				hRcNSigmaEvsPQ->Fill(q,rcP,rcNSigmaE,weightPlus);
				hRcNSigmaEvsEtaQ->Fill(q,rcEta,rcNSigmaE,weightPlus);
				hRcNSigmaEvsPhiQ->Fill(q,rcPhi,rcNSigmaE,weightPlus);
				hRcNSigmaPivsPQ->Fill(q,rcP,rcNSigmaPi,weightPlus);
				hRcNSigmaKvsPQ->Fill(q,rcP,rcNSigmaK,weightPlus);
				hRcNSigmaPvsPQ->Fill(q,rcP,rcNSigmaP,weightPlus);
			}
			if (q < 0)
			{
				hRcDedxvsPQ->Fill(q,rcP,rcDedx,weightMinus);
				hRcDedxvsEtaQ->Fill(q,rcEta,rcDedx,weightMinus);
				hRcDedxvsPhiQ->Fill(q,rcPhi,rcDedx,weightMinus);
				hRcNSigmaEvsPQ->Fill(q,rcP,rcNSigmaE,weightMinus);
				hRcNSigmaEvsEtaQ->Fill(q,rcEta,rcNSigmaE,weightMinus);
				hRcNSigmaEvsPhiQ->Fill(q,rcPhi,rcNSigmaE,weightMinus);
				hRcNSigmaPivsPQ->Fill(q,rcP,rcNSigmaPi,weightMinus);
				hRcNSigmaKvsPQ->Fill(q,rcP,rcNSigmaK,weightMinus);
				hRcNSigmaPvsPQ->Fill(q,rcP,rcNSigmaP,weightMinus);
			}
		}

		//Int_t nMcTrks = event->nMcTrks;
		Int_t nMcTrks = event->nMcE;
		hNMatchTrksvsInputTrks->Fill(nMcTrks,nMatchedE);
	}

	writeHistograms(outFile);
	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//____________________________________________________________
void bookHistograms()
{
	hVyvsVx = new TH2F("hVyvsVx","hVyvsVx;Vx (cm); Vy (cm)",150,-1.5,1.5,150,-1.5,1.5);
	hVyvsVz = new TH2F("hVyvsVz","hVyvsVz;Vz (cm); Vy (cm)",240,-60,60,150,-1.5,1.5);
	hVzDiff = new TH1F("hVzDiff","hVzDiff;Vz_{TPC} - Vz_{VPD} (cm);Counts",200,-10,10);
	hRcMcVxDiffvsMcVx = new TH2F("hRcMcVxDiffvsMcVx","hRcMcVxDiffvsMcVx; McV_{x}; RcV_{x} - McV_{x}",150,-1.5,1.5,200,-1.,1.);
	hRcMcVyDiffvsMcVy = new TH2F("hRcMcVyDiffvsMcVy","hRcMcVyDiffvsMcVy; McV_{y}; RcV_{y} - McV_{y}",150,-1.5,1.5,200,-1.,1.);
	hRcMcVzDiffvsMcVz = new TH2F("hRcMcVzDiffvsMcVz","hRcMcVzDiffvsMcVz; McV_{z}; RcV_{z} - McV_{z}",240,-60,60,500,-5.,5.);
	hRefMultvsRefMultCorr = new TH2F("hRefMultvsRefMultCorr","hRefMultvsRefMultCorr; refMultCorr; refMult",1000,0,1000,1000,0,1000);
	hNMatchTrksvsInputTrks = new TH2F("hNMatchTrksvsInputTrks","hNMatchTrksvsInputTrks; # of InputTrks(pGeantID==0); # of MatchedTrks",50,0,50,50,0,50);

	hMcEtavsPtQ = new TH3F("hMcEtavsPtQ","hMcEtavsPtQ;q;p_{T} (GeV/c);#eta",2,-1,1,500,0,10,200,-2,2);
	hMcPhivsPtQ = new TH3F("hMcPhivsPtQ","hMcPhivsPtQ;q;p_{T} (GeV/c);#phi",2,-1,1,500,0,10,360,-PI,PI);
	hMcYvsPtQ = new TH3F("hMcYvsPtQ","hMcYvsPtQ;q;p_{T} (GeV/c);Rapidity",2,-1,1,500,0,10,200,-2,2);
	hRcEtavsPtQ = new TH3F("hRcEtavsPtQ","hRcEtavsPtQ;q;p_{T} (GeV/c);#eta",2,-1,1,500,0,10,200,-2,2);
	hRcPhivsPtQ = new TH3F("hRcPhivsPtQ","hRcPhivsPtQ;q;p_{T} (GeV/c);#phi",2,-1,1,500,0,10,360,-PI,PI);
	hRcYvsPtQ = new TH3F("hRcYvsPtQ","hRcYvsPtQ;q;p_{T} (GeV/c);Rapidity",2,-1,1,500,0,10,200,-2,2);
	hRcPtvsMcPtQ = new TH3F("hRcPtvsMcPtQ","hRcPtvsMcPtQ;q;Mc p_{T} (GeV/c);Rc p_{T} (GeV/c)",2,-1,1,500,0,10,500,0,10);
	hRcEtavsMcEtaQ = new TH3F("hRcEtavsMcEtaQ","hRcEtavsMcEtaQ;q;Mc #eta;Rc #eta",2,-1,1,200,-2,2,200,-2,2);
	hRcPhivsMcPhiQ = new TH3F("hRcPhivsMcPhiQ","hRcPhivsMcPhiQ;q;Mc #phi;Rc #phi",2,-1,1,360,-PI,PI,360,-PI,PI);
	hRcYvsMcYQ = new TH3F("hRcYvsMcYQ","hRcYvsMcYQ;q;Mc Rapidity;Rc Rapidity",2,-1,1,200,-2,2,200,-2,2);
	hPtResvsPtQ = new TH3F("hPtResvsPtQ","hPtResvsPtQ;q;p_{T}^{MC} (GeV/c);(p_{T}^{RC}-p_{T}^{MC})/p_{T}^{MC}",2,-1,1,500,0,10,4000,-2,2);
	hPResvsPQ = new TH3F("hPResvsPQ","hPResvsPQ;q;p_{MC} (GeV/c);(p_{RC}-p_{MC})/p_{MC}",2,-1,1,500,0,10,4000,-2,2);
  hPtResvsPtCen = new TH3F("hPtResvsPtCen","hPtResvsPtCen;Centrality;p_{T}^{MC} (GeV/c);(p_{T}^{RC}-p_{T}^{MC})/p_{T}^{MC}",13,-1,12,500,0,10,4000,-2,2);

	hRcNHitsFitvsPtQ = new TH3F("hRcNHitsFitvsPtQ","hRcNHitsFitvsPtQ;q;p_{T} (GeV/c);nHitsFit",2,-1,1,500,0,10,50,0,50);
	hRcNHitsFitvsEtaQ = new TH3F("hRcNHitsFitvsEtaQ","hRcNHitsFitvsEtaQ;q;#eta;nHitsFit",2,-1,1,200,-2,2,50,0,50);
	hRcNHitsFitvsPhiQ = new TH3F("hRcNHitsFitvsPhiQ","hRcNHitsFitvsPhiQ;q;#phi;nHitsFit",2,-1,1,360,-PI,PI,50,0,50);
	hRcNHitsPossvsPtQ = new TH3F("hRcNHitsPossvsPtQ","hRcNHitsPossvsPtQ;q;p_{T} (GeV/c);nHitsPoss",2,-1,1,500,0,10,50,0,50);
	hRcNHitsPossvsEtaQ = new TH3F("hRcNHitsPossvsEtaQ","hRcNHitsPossvsEtaQ;q;#eta;nHitsPoss",2,-1,1,200,-2,2,50,0,50);
	hRcNHitsPossvsPhiQ = new TH3F("hRcNHitsPossvsPhiQ","hRcNHitsPossvsPhiQ;q;#phi;nHitsPoss",2,-1,1,360,-PI,PI,50,0,50);
	hRcNHitsDedxvsPtQ = new TH3F("hRcNHitsDedxvsPtQ","hRcNHitsDedxvsPtQ;q;p_{T} (GeV/c);nHitsDedx",2,-1,1,500,0,10,50,0,50);
	hRcNHitsCommonvsPtQ = new TH3F("hRcNHitsCommonvsPtQ","hRcNHitsCommonvsPtQ;q;p_{T} (GeV/c);nHitsCommon",2,-1,1,500,0,10,50,0,50);
	hRcNHitsCommonvsEtaQ = new TH3F("hRcNHitsCommonvsEtaQ","hRcNHitsCommonvsEtaQ;q;#eta;nHitsCommon",2,-1,1,200,-2,2,50,0,50);
	hRcNHitsCommonvsPhiQ = new TH3F("hRcNHitsCommonvsPhiQ","hRcNHitsCommonvsPhiQ;q;#phi;nHitsCommon",2,-1,1,360,-PI,PI,50,0,50);
	hRcDcavsPtQ = new TH3F("hRcDcavsPtQ","hRcDcavsPtQ;q;p_{T} (GeV/c);dca (cm)",2,-1,1,500,0,10,2000,-1.e-6,4-1.e-6);
	hRcDedxvsPQ = new TH3F("hRcDedxvsPQ","hRcDedxvsPQ;q;p (GeV/c);dE/dx (KeV/cm)",2,-1,1,500,0,10,300,-1.e-6,15-1.e-6);
	hRcDedxvsEtaQ = new TH3F("hRcDedxvsEtaQ","hRcDedxvsEtaQ;q;#eta;dE/dx (KeV/cm)",2,-1,1,200,-2,2,300,-1.e-6,15-1.e-6);
	hRcDedxvsPhiQ = new TH3F("hRcDedxvsPhiQ","hRcDedxvsPhiQ;q;#phi;dE/dx (KeV/cm)",2,-1,1,360,-PI,PI,300,-1.e-6,15-1.e-6);
	hRcNSigmaEvsPQ = new TH3F("hRcNSigmaEvsPQ","hRcNSigmaEvsPQ;q;p (GeV/c);n#sigma_{e}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hRcNSigmaEvsEtaQ = new TH3F("hRcNSigmaEvsEtaQ","hRcNSigmaEvsEtaQ;q;#eta;n#sigma_{e}",2,-1,1,200,-2,2,600,-15-1.e-6,15-1.e-6);
	hRcNSigmaEvsPhiQ = new TH3F("hRcNSigmaEvsPhiQ","hRcNSigmaEvsPhiQ;q;#phi;n#sigma_{e}",2,-1,1,360,-PI,PI,600,-15-1.e-6,15-1.e-6);
	hRcNSigmaPivsPQ = new TH3F("hRcNSigmaPivsPQ","hRcNSigmaPivsPQ;q;p (GeV/c);n#sigma_{#pi}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hRcNSigmaKvsPQ = new TH3F("hRcNSigmaKvsPQ","hRcNSigmaKvsPQ;q;p (GeV/c);n#sigma_{k}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);
	hRcNSigmaPvsPQ = new TH3F("hRcNSigmaPvsPQ","hRcNSigmaPvsPQ;q;p (GeV/c);n#sigma_{p}",2,-1,1,500,0,10,600,-15-1.e-6,15-1.e-6);

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

	hDenEPlusTpcEff = new TH3F("hDenEPlusTpcEff","hDenEPlusTpcEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);
	hNumEPlusTpcEff = new TH3F("hNumEPlusTpcEff","hNumEPlusTpcEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);
	hDenEMinusTpcEff = new TH3F("hDenEMinusTpcEff","hDenEMinusTpcEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);
	hNumEMinusTpcEff = new TH3F("hNumEMinusTpcEff","hNumEMinusTpcEff;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);

	for(Int_t i=0;i<9;i++){
		hDenEPlusTpcEffCen[i] = new TH3F(Form("hDenEPlusTpcEffCenBin%d",i),"hDenEPlusTpcEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);
		hNumEPlusTpcEffCen[i] = new TH3F(Form("hNumEPlusTpcEffCenBin%d",i),"hNumEPlusTpcEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);
		hDenEMinusTpcEffCen[i] = new TH3F(Form("hDenEMinusTpcEffCenBin%d",i),"hDenEMinusTpcEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);
		hNumEMinusTpcEffCen[i] = new TH3F(Form("hNumEMinusTpcEffCenBin%d",i),"hNumEMinusTpcEffCen;p_{T} (GeV/c); #eta; #phi",200,0,10,100,-1,1,120,-PI-PI/12,PI+PI/12);
	}

	hMinusOriginZ = new TH1F("hMinusOriginZ","hMinusOriginZ",600,-30,30);
	hMinusOriginPhi = new TH1F("hMinusOriginPhi","hMinusOriginPhi",360,-PI,PI);
	hPlusOriginZ = new TH1F("hPlusOriginZ","hPlusOriginZ",600,-30,30);
	hPlusOriginPhi = new TH1F("hPlusOriginPhi","hPlusOriginPhi",360,-PI,PI);

	hRefMult = new TH1F("hRefMult","hRefMult;refmult;counts",600,0,600);
    hRefMultCorr = new TH1F("hRefMultCorr","hRefMultCorr;refmultcorr;counts",600,0,600);
    hCentrality9 = new TH1F("hCentrality9","hCentrality9",10,0,10);
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
	hVzDiff->Write();
	hRcMcVxDiffvsMcVx->Write();
	hRcMcVyDiffvsMcVy->Write();
	hRcMcVzDiffvsMcVz->Write();
	hRefMultvsRefMultCorr->Write();
	hNMatchTrksvsInputTrks->Write();

	hMcEtavsPtQ->Write();
	hMcPhivsPtQ->Write();
	hMcYvsPtQ->Write();
	hRcEtavsPtQ->Write();
	hRcPhivsPtQ->Write();
	hRcYvsPtQ->Write();
	hRcPtvsMcPtQ->Write();
	hRcEtavsMcEtaQ->Write();
	hRcPhivsMcPhiQ->Write();
	hRcYvsMcYQ->Write();
	hPtResvsPtQ->Write();
	hPResvsPQ->Write();
	hRcDcavsPtQ->Write();
	hRcNHitsFitvsPtQ->Write();
	hRcNHitsFitvsEtaQ->Write();
	hRcNHitsFitvsPhiQ->Write();
	hRcNHitsPossvsPtQ->Write();
	hRcNHitsPossvsEtaQ->Write();
	hRcNHitsPossvsPhiQ->Write();
	hRcNHitsDedxvsPtQ->Write();
	hRcNHitsCommonvsPtQ->Write();
	hRcNHitsCommonvsEtaQ->Write();
	hRcNHitsCommonvsPhiQ->Write();
	hRcDedxvsPQ->Write();
	hRcDedxvsEtaQ->Write();
	hRcDedxvsPhiQ->Write();
	hRcNSigmaEvsPQ->Write();
	hRcNSigmaEvsEtaQ->Write();
	hRcNSigmaEvsPhiQ->Write();
	hRcNSigmaPivsPQ->Write();
	hRcNSigmaKvsPQ->Write();
	hRcNSigmaPvsPQ->Write();

	hPENHitsFitsvsPtvsEtaMinus->Write();
	hPENHitsFitsvsPtvsPhiMinus->Write();
	hPENHitsFitsvsPtvsEtaPlus->Write();
	hPENHitsFitsvsPtvsPhiPlus->Write();
	hPENHitsDedxvsPtvsEtaMinus->Write();
	hPENHitsDedxvsPtvsPhiMinus->Write();
	hPENHitsDedxvsPtvsEtaPlus->Write();
	hPENHitsDedxvsPtvsPhiPlus->Write();
	hPEDcavsPtvsEtaMinus->Write();
	hPEDcavsPtvsPhiMinus->Write();
	hPEDcavsPtvsEtaPlus->Write();
	hPEDcavsPtvsPhiPlus->Write();


	for(Int_t i=0;i<9;i++){
		hDenEPlusTpcEffCen[i]->Write();
		hNumEPlusTpcEffCen[i]->Write();
		hDenEMinusTpcEffCen[i]->Write();
		hNumEMinusTpcEffCen[i]->Write();
	}

	hDenEPlusTpcEff->Write();
	hNumEPlusTpcEff->Write();
	hDenEMinusTpcEff->Write();
	hNumEMinusTpcEff->Write();

	hMinusOriginZ ->Write();
	hMinusOriginPhi ->Write();
	hPlusOriginZ ->Write();
	hPlusOriginPhi ->Write();

	hCentrality9->Write();
    hRefMult->Write();
    hRefMultCorr->Write();
    hPtResvsPtCen->Write();
}
//____________________________________________________________
bool Init()
{
	ifstream indata;

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

	TFile* infile = new TFile("/star/u/wangzhen/QA/wangzhen/embedding/myEmbedding/basicQA/mc/weight.root");
	pTEtaWeight[0] = (TH2D*)infile->Get("pTEtaWeightMinus");
	pTEtaWeight[1] = (TH2D*)infile->Get("pTEtaWeightPlus");

	cout<<"Initialization DONE !!!"<<endl;

	return kTRUE;
}

