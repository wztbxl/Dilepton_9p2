#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "sys/types.h"
#include "dirent.h"

#include "math.h"
#include "string.h"
//Add the data structure
#include <iomanip>
//#include <stdio.h> 

using std::cout;
using std::endl;

#ifndef __CINT__  
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h" 
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "VecMeson.h"
#endif

//using std::cout;
//using std::endl;

int main(int argc, char** argv){

	system("mkdir -p output");

	VecMeson *vm;
	char outputName[256];
	double ptSmearParameters[2] = {0.0};

	if(argc==1){
		vm = new VecMeson(virtualphoton,twobody);
		//vm = new VecMeson(pi0,dalitz);
		vm->SetNumberOfTracks(1.e4);
		vm->SetCentralityIdx(0);
		//vm->SetUseCocktailInput(0);
		cout<<"Input meson:"<<vm->MesonType<<endl;
		system("mkdir -p output/test");
		sprintf(outputName,"output/test/test.root");
	}
	else if(argc==6){
		ParticleTypes ptype = ParticleTypes(atoi(argv[2]));
		cout << ptype <<endl;
		DecayMode dmode = DecayMode(atoi(argv[3]));
		//from 19.6 GeV data
		ptSmearParameters[0] = 5.713e-3;
		// ptSmearParameters[1] = 0.007473;
		ptSmearParameters[1] = 7.92e-3;
		int step = atoi(argv[5]);
		//ptSmearParameters[0] = 0.008000+step*0.00004;//used to produce the pt smear scan
		//ptSmearParameters[1] = 0.007473;

		vm = new VecMeson(ptype,dmode);
		//vm->SetNumberOfTracks(2.e7);
		//vm->SetRapidityRange(-1., 1.);
		//vm->SetNumberOfTracks(2.4e8);
     	vm->SetNumberOfTracks(5.e7);
		vm->SetRapidityRange(-1.2, 1.2);
		//vm->SetUseCocktailInput(0); //default is 1 -- use 2-D cocktail to sample virtual photon pt and mass
		//vm->SetUseScaleEff(1);      //default is 1; 0-use specialized 3-D TPC&TOF efficiency for each centrality bin, 1-scale 3-D MB TPC&TOF efficiencies for each centrality bin
		//vm->SetInputPtSpectra(1);   //default is 1; 0-PHENIX(mT scaling), 1-STAR(Tsallis fit)
		vm->SetPtSmearPar(ptSmearParameters);//set the pTSmearParameters to get the find parameters

		if (strcmp(argv[1],"080")==0) {
			vm->SetCentralityIdx(0);
		}
		else if (strcmp(argv[1],"010")==0) {
			vm->SetCentralityIdx(1);
		}
		else if (strcmp(argv[1],"1040")==0) {
			vm->SetCentralityIdx(2);
		}
		else if (strcmp(argv[1],"4080")==0) {
			vm->SetCentralityIdx(3);
		}
		else if (strcmp(argv[1],"4060")==0) {
			vm->SetCentralityIdx(4);
		}
		else if (strcmp(argv[1],"6080")==0) {
			vm->SetCentralityIdx(5);
		}
		else if (strcmp(argv[1],"6070")==0) {
			vm->SetCentralityIdx(6);
		}
		else if (strcmp(argv[1],"7080")==0) {
			vm->SetCentralityIdx(7);
		}
		else{
			cout<<"Centrality string is wrong !"<<endl;
			return -1;
		}

		cout<<"Meson is "<<vm->MesonType<<endl;
		if(strcmp(vm->MesonType, argv[4])!=0){
			cout<<"Meson's name is wrong !"<<endl;
			return -1;
		}

		system(Form("mkdir -p output/Cen%s",argv[1]));
		system(Form("mkdir -p output/Cen%s/%s",argv[1],argv[4]));
		if(dmode==2) sprintf(outputName,"output/Cen%s/%s/%s2ee_%s.root",argv[1],argv[4],argv[4],argv[5]);
		if(dmode==3) sprintf(outputName,"output/Cen%s/%s/%sdalitz_%s.root",argv[1],argv[4],argv[4],argv[5]);
	}
	else{
		cout<<"The # of arguments is wrong !"<<endl;
		return -1;
	}

	vm->Init();

	cout<<"------------ processing decay ----------"<<endl;
	vm->GenerateDecay();

	cout<<"------------ Writing file ----------"<<endl;
	TFile *fout = new TFile(outputName,"recreate");;
	fout->cd();
	vm->hSampledPt->Write();
	vm->hEPSingleTrkEffvsPt->Write();
	vm->hEMSingleTrkEffvsPt->Write();
	vm->hRCPairRapidityvsParentRapidity->Write();
	vm->hMCPairPtvsParentPt->Write();
	vm->hMCAcc0PairPtvsParentPt->Write();
	vm->hMCAcc1PairPtvsParentPt->Write();
	vm->hMCMvsPt->Write();
	vm->hMCAcc0MvsPt->Write();
	vm->hMCAcc1MvsPt->Write();
	vm->hMCAcc2MvsPt->Write();
	vm->hRCAcc1MvsPt3D->Write();
	vm->hRCAcc2MvsPt3D->Write();
	vm->hMCAcc1PairRapidity->Write();
	vm->hMCAcc1PairEMPtvsEPPt->Write();
	vm->hPlusSmearvsPt->Write();
	vm->hMinusSmearvsPt->Write();
  	vm->RapiditySTARAcc->Write();
  	vm->RapiditywoSTARAcc->Write();
  	vm->MeePtRapidityOver1->Write();
  	vm->MeeFullRapidity->Write();
	vm->MeeWopTCut->Write();
	vm->MeeWopTEtaCut->Write();
	vm->MeeWoEtaCut->Write();
	vm->pTWithLargeMee->Write();
	vm->RapidityOutSTARAcc->Write();
	vm->MeeWpTEtaWoRapidity->Write();
	vm->MeeWpTEtaWRapidity->Write();
	vm->PtMVsPtPLargeMass->Write();
	vm->EtaMVsEtaPLargeMass->Write();
	vm->CutRecorder->Write();
	vm->hMCAcc1MvsPtwoSmear->Write();
	static const Int_t mPtBins = 10;
	static const Int_t mYBins = 20;
	static const Int_t mPhiBins= 7;
	static const Int_t mCenBins = 9; //16; //9;
	for (int i = 0; i < mCenBins; i++)
	{
		for( int j = 0; j < mPtBins; j++)
		{
			for (int k = 0; k < mPhiBins; k++)
			{
				vm->hMCAcc0Mass[i][j][k]->Write();
				vm->hMCAcc1Mass[i][j][k]->Write();
				vm->hRCAcc1Mass[i][j][k]->Write();
			}
		}
	}
	vm->hMCAcc0_CosthetapT->Write();
	vm->hMCAcc1_CosthetapT->Write();
	vm->hRCAcc1_CosthetapT->Write();
	vm->hMCAcc0PairCosThetaPt_HX->Write();
	vm->hMCAcc1PairCosThetaPt_HX->Write();
	vm->hRCAcc1PairCosThetaPt_HX->Write();
	vm->hMCAcc0PairCosThetaPt_CS->Write();
	vm->hMCAcc1PairCosThetaPt_CS->Write();
	vm->hRCAcc1PairCosThetaPt_CS->Write();
	vm->hMCAcc0PairPhiPt_HX->Write();
	vm->hMCAcc1PairPhiPt_HX->Write();
	vm->hRCAcc1PairPhiPt_HX->Write();
	vm->hMCAcc0PairPhiPt_CS->Write();
	vm->hMCAcc1PairPhiPt_CS->Write();
	vm->hRCAcc1PairPhiPt_CS->Write();
	vm->hMCAcc0PairCosThetaInvMPt_HX->Write();
	vm->hMCAcc1PairCosThetaInvMPt_HX->Write();
	vm->hRCAcc1PairCosThetaInvMPt_HX->Write();
	vm->hMCAcc0PairCosThetaInvMPt_CS->Write();
	vm->hMCAcc1PairCosThetaInvMPt_CS->Write();
	vm->hRCAcc1PairCosThetaInvMPt_CS->Write();
	vm->hMCAcc0PairPhiInvMPt_HX->Write();
	vm->hMCAcc1PairPhiInvMPt_HX->Write();
	vm->hRCAcc1PairPhiInvMPt_HX->Write();
	vm->hMCAcc0PairPhiInvMPt_CS->Write();
	vm->hMCAcc1PairPhiInvMPt_CS->Write();
	vm->hRCAcc1PairPhiInvMPt_CS->Write();
	

	//fout->Write();
	fout->Close();

	cout<<"output file: "<<outputName<<" has been written!"<<endl;
	cout << "------------end of program-----------------" << endl;
	return 0;
}
