#ifndef EFFUTIL_H
#define EFFUTIL_H

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

/*
 * Centrality Sequence:
 * 0 - 0-80%
 * 1 - 0-10%
 * 2 - 10-40%
 * 3 - 40-80%
 * 4 - 40-60%
 * 5 - 60-80%
 * 6 - 60-70%
 * 7 - 70-80%
 */
const Int_t nCenBins = 4; 
const Int_t mCenBinLow[nCenBins]    = {0,  14, 8,  0}; // 16 Centrality , need to change
const Int_t mCenBinHi[nCenBins]     = {15, 15, 13, 7};
const Int_t CentralityLow[nCenBins] = {0,  0,  10, 40};  // 
const Int_t CentralityHi[nCenBins]  = {80, 10, 40, 80};
const Int_t NPt_TPC[nCenBins]  = {40, 20, 20, 20}; //0~4 GeV/c
const Int_t NEta_TPC[nCenBins] = {10, 10, 10, 10};
const Int_t NPhi_TPC[nCenBins] = {36, 36, 36, 36};
const Int_t NPt_TOF[nCenBins]  = {80, 40, 40, 40}; //0~4 GeV/c
const Int_t NEta_TOF[nCenBins] = {10, 10, 10, 10};
const Int_t NPhi_TOF[nCenBins] = {36, 36, 36, 36};
int CenIdx = -999;

TF1    *funndEdxEffPt;
TF1    *funnSigEffP;
TF1    *funBetaEffP;
TF1    *fEPiEff_Tof;        
TF1    *fMBEPiEff_Tof;        
TF1    *fTpcEffRatioToMB;
TF1    *fTofEffRatioToMB;

const double PI = 3.1415926;
const Double_t EpsilonVal = 1.e-9;
const Int_t    maxEtaBins = 40;
const Int_t    maxPhiBins = 120;
const Int_t nEta_TPC = 10;//3
const Int_t nPhi_TPC = 36;//4
const Int_t nEta_TOF = 10;//3
const Int_t nPhi_TOF = 36;//4
const Int_t mPHENIX = 0;
const double ptl_Tpc = 0.2;
const double pth_Tpc = 2.5;
const double ptl_Tof = 0.2;
const double pth_Tof = 3;//pt limit in tof eff
TH1D   *hEff_Tpc_Pos[maxEtaBins][maxPhiBins];
TH1D   *hEff_Tpc_Neg[maxEtaBins][maxPhiBins];
TH1D   *hEff_Tof_Pos[maxEtaBins][maxPhiBins];
TH1D   *hEff_Tof_Neg[maxEtaBins][maxPhiBins];
TH1D   *hMBEff_Tpc_Pos[maxEtaBins][maxPhiBins];
TH1D   *hMBEff_Tpc_Neg[maxEtaBins][maxPhiBins];
TH1D   *hMBEff_Tof_Pos[maxEtaBins][maxPhiBins];
TH1D   *hMBEff_Tof_Neg[maxEtaBins][maxPhiBins];

TH1D* hTOFMatchEff[2];
TH1D* hElecPionRatio[2];
TF1*  f_ElecPoionRatio;
TF1*  f_betaCutEff;
TF1*  f_nSigmaEEff_lowpT;
TF1*  f_nSigmaEEff_HigpT;
TH1D* hTPCTrackingEffPlus[nEta_TPC];
TH1D* hTPCTrackingEffMinus[nEta_TPC];
TH1D* hBetaCutEff[2];
TH1D* hNSigmaECutEff[2];

//-----------------------Zhen add this for BES-II Dielectron analysis--------------------------//
//overall : add centrality selection part

void setCentrality(int icent)
{
	CenIdx = icent;
}

//Fit function for the e/pi ratio
double BES_II_EffRatio(double *x, double *par)
{
    double temp;
    double Num = 1+par[3]*x[0]+par[4]*x[0]*x[0];
    double Den = par[0]+TMath::Exp((x[0]-par[1])/par[2]);
    temp = Num/Den+par[5];
    return temp;
}

//Zhen add it to read the efficiency histograms
void InitialzeEffHist(int CenIdx = -999)
{
	//TOF efficiency : load histograms for pion efficiency of each eta and phi bin, and set the function of e/pi ratio
	TString Flag[2] = {"Plus","Minus"};
	//read TOF match efficiency
	// TFile *f1 = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/Effinput/TOFMatchEffhis.root");
	TFile *f1 = new TFile(Form("./Effinput/TOFMatchEffhis%d%d.root",CentralityLow[CenIdx],CentralityHi[CenIdx]));
	if(f1->IsOpen()) cout<<"TOF match efficiency file is open "<<endl;
	else cout<<"Fail to read TOF match efficiency file"<<endl;
	TString name;
	for (int i = 0; i < nEta_TOF; i++)
	{
		for (int j = 0; j < nPhi_TOF; j++)
		{
			name = Form("PionEff%sEta%dPhi%d",Flag[0].Data(),i,j);
			hMBEff_Tof_Pos[i][j] = (TH1D*)f1->Get(name);
			name = Form("PionEff%sEta%dPhi%d",Flag[1].Data(),i,j);
			hMBEff_Tof_Neg[i][j] = (TH1D*)f1->Get(name);
		}
	}
	cout<<"read TOF match efficiency ok"<<endl;
	f_ElecPoionRatio = new TF1("f_ElecPoionRatio",BES_II_EffRatio,0.2,5,6);
	switch (CenIdx)
	{
	case 0:
		f_ElecPoionRatio->SetParameters( 3.42168, 0.830217, 0.153163,-3.09743, 2.68287, 1.03397);//0-80%
		break;
	case 1:
		f_ElecPoionRatio->SetParameters( 3.42168, 0.830217, 0.153163,-3.09743, 2.68287, 1.03397);//0-10%
		break;
	case 2:
		f_ElecPoionRatio->SetParameters( 3.42168, 0.830217, 0.153163,-3.09743, 2.68287, 1.03397);//10-40%
		break;
	case 3:
		f_ElecPoionRatio->SetParameters( 3.42168, 0.830217, 0.153163,-3.09743, 2.68287, 1.03397);//40-80%
		break;
	
	default: cout << "you need to select a centraility index!!!!!!!!" << endl;
		break;
	}	
	//next step: reading a text file to get the parameters

	//get TPC tracking 
	TFile* f2 = new TFile(Form("./Effinput/TPCEffHisto%d%d.root",CentralityLow[CenIdx],CentralityHi[CenIdx]));
	if(f1->IsOpen()) cout<<"TPC match efficiency file is open "<<endl;
	else cout<<"Fail to read TPC match efficiency file"<<endl;
	for (int i = 0;i<nEta_TPC;i++)
	{
		for (int j = 0; j < nPhi_TPC; j++)
		{
			hMBEff_Tpc_Pos[i][j] = (TH1D*)f2->Get(Form("TPCEFFPlusEta%dPhi%d",i,j));
			hMBEff_Tpc_Neg[i][j] = (TH1D*)f2->Get(Form("TPCEFFMinusEta%dPhi%d",i,j));
		}
	}
	cout<<"read TPC match efficiency ok"<<endl;

	//get beta cut Eff
	//beta cut eff using a const 
	f_betaCutEff = new TF1("f_betaCutEff","pol0",0.2,6);
	switch (CenIdx)
	{
	case 0:
		f_betaCutEff->SetParameter(0,0.9689);//0-80%
		break;
	case 1:
		f_betaCutEff->SetParameter(0,0.9688);//0-10%
		break;
	case 2:
		f_betaCutEff->SetParameter(0,0.9690);//10-40%
		break;
	case 3:
		f_betaCutEff->SetParameter(0,0.9681);//40-80%
		break;
	
	default: cout << "you need to select a centraility index!!!!!!!!"<< endl;
		break;
	}	
	cout<<"read beta cut efficiency ok"<<endl;

	//get nSigmaE cut Eff
	f_nSigmaEEff_lowpT = new TF1("f_nSigmaEEff_lowpT","pol3",0.2,0.8);
	f_nSigmaEEff_HigpT = new TF1("f_nSigmaEEff_HigpT","pol1",0.8,2);
	f_betaCutEff = new TF1("f_betaCutEff","pol0",0.2,6);
	//for different centrality, now do not have result
	switch (CenIdx)
	{
	case 0:
		f_nSigmaEEff_lowpT->SetParameters(0.979,0.042,0.148,-0.617);//0-80%
		f_nSigmaEEff_HigpT->SetParameters(0.700,0.068);//0-80%
		break;
	case 1:
		f_nSigmaEEff_lowpT->SetParameters(0,0,0,0);//0-10%
		f_nSigmaEEff_HigpT->SetParameters(0,0);//0-10%
		break;
	case 2:
		f_nSigmaEEff_lowpT->SetParameters(0,0,0,0);//10-40%
		f_nSigmaEEff_HigpT->SetParameters(0,0);//10-40%
		break;
	case 3:
		f_nSigmaEEff_lowpT->SetParameters(0,0,0,0);//40-80%
		f_nSigmaEEff_HigpT->SetParameters(0,0);//40-80%
		break;
	
	default: cout << "you need to select a centraility index!!!!!!!!"<< endl;
		break;
	}	
	cout<<"read nSigmaE cut efficiency ok"<<endl;
}

double getEff(double pt, TH1D* histo , double ptLLimit, double ptHLimit)
{
	if (pt<ptLLimit) return 0;
	if (pt>ptHLimit) pt = ptHLimit;
	int Bin = histo->GetXaxis()->FindBin(pt+1.e-8);
	if (Bin<0) return 0;
	if (Bin>histo->GetNbinsX()) return histo->GetBinContent(histo->GetNbinsX()) ;
	return histo->GetBinContent(Bin);
}

void getEtaPhiBin_TPC(double Eta, double Phi, int *EtaBin, int *PhiBin)//phi range is 0~2pi
{
    *EtaBin = 0; *PhiBin = 0;
    *EtaBin = (int)((Eta+1)/(2.0/nEta_TPC));
    *PhiBin = (int)((Phi+PI)/(2*PI/nPhi_TPC));
    if (*EtaBin<0.) *EtaBin = 0;
    if (*EtaBin>nEta_TPC) *EtaBin = nEta_TPC-1;
    if (*PhiBin<0.) *PhiBin = 0;
    if (*PhiBin>nEta_TPC) *PhiBin = nPhi_TPC-1;
    return;
}

void getEtaPhiBin_TOF(double Eta, double Phi, int *EtaBin, int *PhiBin)//phi range is 0~2pi
{
    *EtaBin = 0; *PhiBin = 0;
    *EtaBin = (int)((Eta+1)/(2.0/nEta_TOF));
    *PhiBin = (int)((Phi+PI)/(2*PI/nPhi_TOF));
    if (*EtaBin<0.) *EtaBin = 0;
    if (*EtaBin>nEta_TOF) *EtaBin = nEta_TOF-1;
    if (*PhiBin<0.) *PhiBin = 0;
    if (*PhiBin>nEta_TOF) *PhiBin = nPhi_TOF-1;
    return;
}
//-----------------------Zhen add this for BES-II Dielectron analysis--------------------------//

#endif
