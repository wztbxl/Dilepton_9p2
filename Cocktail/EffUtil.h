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
const Int_t nCenBins = 8; 
const Int_t mCenBinLow[nCenBins]    = {0,  14, 8,  0, 4, 0, 2, 0}; // 16 Centrality Bins
const Int_t mCenBinHi[nCenBins]     = {15, 15, 13, 7, 7, 3, 3, 1};
const Int_t CentralityLow[nCenBins] = {0,  0,  10, 40, 40, 60, 60, 70};  // %
const Int_t CentralityHi[nCenBins]  = {80, 10, 40, 80, 60, 80, 70, 80};
const Int_t NPt_TPC[nCenBins]  = {40, 20, 20, 20, 20, 20, 20, 20}; //0~4 GeV/c
const Int_t NEta_TPC[nCenBins] = {20, 10, 10, 10, 2,  2, 2, 2};
const Int_t NPhi_TPC[nCenBins] = {60, 60, 60, 60, 60, 60, 24, 24};
const Int_t NPt_TOF[nCenBins]  = {80, 40, 40, 40, 40, 40, 40, 40}; //0~4 GeV/c
const Int_t NEta_TOF[nCenBins] = {10, 10, 10, 10, 5,  5, 4, 4};
const Int_t NPhi_TOF[nCenBins] = {60, 60, 60, 60, 60, 60, 24, 24};

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
const Int_t nEta_TPC = 4;//3
const Int_t nPhi_TPC = 3;//4
const Int_t mPHENIX = 0;
const double ptl_Tpc = 0.2;
const double pth_Tpc = 1.96;
const double ptl_Tof = 0.2;
const double pth_Tof = 4;//pt limit in tof eff
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
TH1D* hTPCTrackingEffPlus[nEta_TPC];
TH1D* hTPCTrackingEffMinus[nEta_TPC];
TH1D* hBetaCutEff[2];
TH1D* hNSigmaECutEff[2];

Double_t betaEff(Double_t *x, Double_t *par)
{
	if (TMath::Abs(par[0]-0)<EpsilonVal) {
		return 0.9776; //run12 UU 0-80% 
	}
	else if (TMath::Abs(par[0]-1)<EpsilonVal) {
		return 0.9616; //run12 UU 0-10%
	}
	else if (TMath::Abs(par[0]-2)<EpsilonVal) {
		return 0.9766; //run12 UU 10-40%
	}
	else if (TMath::Abs(par[0]-3)<EpsilonVal) {
		return 0.9852; //run12 UU 40-80%
	}
	else if (TMath::Abs(par[0]-4)<EpsilonVal) {
		return 0.9856; //run12 UU 40-60%
	}
	else if (TMath::Abs(par[0]-5)<EpsilonVal) {
		return 0.985; //run12 UU 60-80%
	}
	else if (TMath::Abs(par[0]-6)<EpsilonVal) {
		return 0.9856; //run12 UU 60-70%
	}
	else if (TMath::Abs(par[0]-7)<EpsilonVal) {
		return 0.9831; //run12 UU 70-80%
	}
	else {
		return 0;
	}
}

Double_t ndEdxEff(Double_t *x, Double_t *par)
{
	if (TMath::Abs(par[0]-0)<EpsilonVal) {
		if(x[0]<1.) return 0.5631+3.569*pow(x[0],1)-14.29*pow(x[0],2)+28.92*pow(x[0],3)-31.69*pow(x[0],4)+17.99*pow(x[0],5)-4.144*pow(x[0],6);
		else return 0.9348*exp(-pow(0.4422/x[0],4.328)); //run12 UU 0-80% 
	}
	else if (TMath::Abs(par[0]-1)<EpsilonVal) {
		if(x[0]<1.) return 0.3003+4.587*pow(x[0],1)-16.19*pow(x[0],2)+29.62*pow(x[0],3)-30.02*pow(x[0],4)+16.23*pow(x[0],5)-3.68*pow(x[0],6);
		else return 0.9994*exp(-pow(0.00003019/x[0],0.1751)); //run12 UU 0-10%
	}
	else if (TMath::Abs(par[0]-2)<EpsilonVal) {
		if(x[0]<1.) return 0.5712+3.511*pow(x[0],1)-14.17*pow(x[0],2)+28.87*pow(x[0],3)-31.75*pow(x[0],4)+17.96*pow(x[0],5)-4.071*pow(x[0],6);
		else return 0.9322*exp(-pow(0.674/x[0],9.207)); //run12 UU 10%-40%
	}
	else if (TMath::Abs(par[0]-3)<EpsilonVal) {
		if(x[0]<1.) return 0.7238+2.744*pow(x[0],1)-12.14*pow(x[0],2)+27.29*pow(x[0],3)-33.67*pow(x[0],4)+21.75*pow(x[0],5)-5.763*pow(x[0],6);
		else return 0.9718*exp(-pow(0.3868/x[0],3.831)); //run12 UU 40%-80%
	}
	else if (TMath::Abs(par[0]-4)<EpsilonVal) {
		if(x[0]<1.) return 0.7339+2.627*pow(x[0],1)-11.81*pow(x[0],2)+27.06*pow(x[0],3)-33.89*pow(x[0],4)+22.07*pow(x[0],5)-5.839*pow(x[0],6);
		else return 1.052*exp(-pow(0.006286/x[0],0.4487)); //run12 UU 40%-60%
	}
	else if (TMath::Abs(par[0]-5)<EpsilonVal) {
		if(x[0]<1.) return 0.7161+2.829*pow(x[0],1)-12.38*pow(x[0],2)+27.49*pow(x[0],3)-33.48*pow(x[0],4)+21.42*pow(x[0],5)-5.643*pow(x[0],6);
		else return 0.9622*exp(-pow(0.955/x[0],60.55)); //run12 UU 60%-80%
	}
	else if (TMath::Abs(par[0]-6)<EpsilonVal) {
		if(x[0]<1.) return 0.8237+1.375*pow(x[0],1)-4.678*pow(x[0],2)+6.849*pow(x[0],3)-3.973*pow(x[0],4)-0.00626*pow(x[0],5)+0.5556*pow(x[0],6);
		else return 1.028*exp(-pow(0.0008233/x[0],0.3632)); //run12 UU 60%-70%
	}
	else if (TMath::Abs(par[0]-7)<EpsilonVal) {
		if(x[0]<1.) return 0.5831+4.728*pow(x[0],1)-22.92*pow(x[0],2)+56.88*pow(x[0],3)-77.01*pow(x[0],4)+53.98*pow(x[0],5)-15.3*pow(x[0],6);
		else return 0.9728*exp(-pow(0.4783/x[0],4.913)); //run12 UU 70%-80%
	}
	else {
		return 0;
	}
}

Double_t EPiTofEffRatio(Double_t *x, Double_t *par)
{
	Double_t ratio;

	if (TMath::Abs(par[0]-0)<EpsilonVal) { //run12 UU 0-80% 
		ratio = (1.+3.231*x[0]*x[0]-3.422*x[0])/(5.858+exp((x[0]-0.4633)/0.1283))+1.089;
	}
	else if (TMath::Abs(par[0]-1)<EpsilonVal) { //run12 UU 0-10%
		ratio = (1.+2.89*x[0]*x[0]-3.345*x[0])/(5.909+exp((x[0]-0.7508)/0.04018))+1.069;
	}
	else if (TMath::Abs(par[0]-2)<EpsilonVal) { //run12 UU 10-40%
		ratio = (1.+3.307*x[0]*x[0]-3.162*x[0])/(3.465+exp((x[0]-0.1997)/0.2408))+1.059;
	}
	else if (TMath::Abs(par[0]-3)<EpsilonVal) { //run12 UU 40-80%
		ratio = (1.+3.336*x[0]*x[0]-3.136*x[0])/(3.95+exp((x[0]-0.08137)/0.2626))+1.067;
	}
	else if (TMath::Abs(par[0]-4)<EpsilonVal) { //run12 UU 40-60%
		ratio = (1.+3.235*x[0]*x[0]-3.136*x[0])/(-1.275+exp((x[0]+0.4121)/0.3385))+1.083;
	}
	else if (TMath::Abs(par[0]-5)<EpsilonVal) { //run12 UU 60-80%
		ratio = (1.+3.012*x[0]*x[0]-2.909*x[0])/(4.745+exp((x[0]-0.2186)/0.2546))+1.051;
	}
	else if (TMath::Abs(par[0]-6)<EpsilonVal) { //run12 UU 60-70%
		ratio = (1.+3.066*x[0]*x[0]-2.942*x[0])/(4.696+exp((x[0]-0.2152)/0.2503))+1.051;
	}
	else if (TMath::Abs(par[0]-7)<EpsilonVal) { //run12 UU 70-80%
		ratio = (1.+2.918*x[0]*x[0]-2.753*x[0])/(4.467+exp((x[0]-0.1639)/0.2938))+1.044;
	}
	else {
		ratio = 0;
	}

	return ratio;
}

Double_t TpcEffRatioToMB(Double_t *x, Double_t *par)
{
	Double_t ratio;

	//[0]/(1+[1]*exp(-x/[2]))
	if (TMath::Abs(par[0]-0)<EpsilonVal) { //run12 UU 0-80% 
		ratio = 1.;
	}
	else if (TMath::Abs(par[0]-1)<EpsilonVal) { //run12 UU 0-10%
		ratio = 0.869/(1+0.1183*exp(-x[0]/0.4292));
	}
	else if (TMath::Abs(par[0]-2)<EpsilonVal) { //run12 UU 10-40%
		ratio = 1.041/(1-0.0236*exp(-x[0]/0.364));
	}
	else if (TMath::Abs(par[0]-3)<EpsilonVal) { //run12 UU 40-80%
		ratio = 1.153/(1-0.1123*exp(-x[0]/0.5039));
	}
	else if (TMath::Abs(par[0]-4)<EpsilonVal) { //run12 UU 40-60%
		ratio = 1.149/(1-0.1107*exp(-x[0]/0.48));
	}
	else if (TMath::Abs(par[0]-5)<EpsilonVal) { //run12 UU 60-80%
		ratio = 1.171/(1-0.1229*exp(-x[0]/0.5754));
	}
	else if (TMath::Abs(par[0]-6)<EpsilonVal) { //run12 UU 60-70%
		ratio = 1.168/(1-0.1213*exp(-x[0]/0.5851));
	}
	else if (TMath::Abs(par[0]-7)<EpsilonVal) { //run12 UU 70-80%
		ratio = 1.177/(1-0.1258*exp(-x[0]/0.5559));
	}
	else {
		ratio = 0;
	}

	return ratio;
}

Double_t TofEffRatioToMB(Double_t *x, Double_t *par)
{
	Double_t ratio;

	//[0]/(1+[1]*exp(-x/[2]))
	if (TMath::Abs(par[0]-0)<EpsilonVal) { //run12 UU 0-80% 
		ratio = 1.;
	}
	else if (TMath::Abs(par[0]-1)<EpsilonVal) { //run12 UU 0-10%
		ratio = 0.955/(1-0.3216*exp(-x[0]/0.08479));
	}
	else if (TMath::Abs(par[0]-2)<EpsilonVal) { //run12 UU 10-40%
		ratio = 0.9941/(1-0.1276*exp(-x[0]/0.08344));
	}
	else if (TMath::Abs(par[0]-3)<EpsilonVal) { //run12 UU 40-80%
		ratio = 1.03/(1+0.2079*exp(-x[0]/0.1053));
	}
	else if (TMath::Abs(par[0]-4)<EpsilonVal) { //run12 UU 40-60%
		ratio = 1.03/(1+0.02641*exp(-x[0]/0.1824));
	}
	else if (TMath::Abs(par[0]-5)<EpsilonVal) { //run12 UU 60-80%
		ratio = 1.03/(1+0.3557*exp(-x[0]/0.09554));
	}
	else if (TMath::Abs(par[0]-6)<EpsilonVal) { //run12 UU 60-70%
		ratio = 1.029/(1+0.3589*exp(-x[0]/0.09342));
	}
	else if (TMath::Abs(par[0]-7)<EpsilonVal) { //run12 UU 70-80%
		ratio = 1.034/(1+0.3442*exp(-x[0]/0.1024));
	}
	else {
		ratio = 0;
	}

	return ratio;
}

Double_t nSigEEff(Double_t *x, Double_t *par)
{
	if (TMath::Abs(par[0]-0)<EpsilonVal) {
		return 1./(exp((x[0]-0.5237)/0.143)+5.372)+0.7845; //run12 UU 0-80% 
	}
	else if (TMath::Abs(par[0]-1)<EpsilonVal) {
		return 1./(exp((x[0]-0.5301)/0.1374)+5.063)+0.7648; //run12 UU 0-10%
	}
	else if (TMath::Abs(par[0]-2)<EpsilonVal) {
		return 1./(exp((x[0]-0.5284)/0.1431)+5.313)+0.7822; //run12 UU 10-40%
	}
	else if (TMath::Abs(par[0]-3)<EpsilonVal) {
		return 1./(exp((x[0]-0.5073)/0.1465)+5.96)+0.808; //run12 UU 40-80%
	}
	else if (TMath::Abs(par[0]-4)<EpsilonVal) {
		return 1./(exp((x[0]-0.5093)/0.1455)+5.898)+0.8059; //run12 UU 40-60%
	}
	else if (TMath::Abs(par[0]-5)<EpsilonVal) {
		return 1./(exp((x[0]-0.5192)/0.1338)+6.401)+0.819; //run12 UU 60-80%
	}
	else if (TMath::Abs(par[0]-6)<EpsilonVal) {
		return 1./(exp((x[0]-0.5214)/0.1317)+6.416)+0.8191; //run12 UU 60-70%
	}
	else if (TMath::Abs(par[0]-7)<EpsilonVal) {
		return 1./(exp((x[0]-0.5029)/0.1368)+6.621)+0.8252; //run12 UU 70-80%
	}
	else {
		return 0;
	}
}

//Zhen add it to use bin by bin efficiency
void InitialzeEffHist()
{
	TString Flag[2] = {"Plus","Minus"};
	//read TOF match efficiency
	TFile *f1 = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/Effinput/TOFMatchEffhis.root");
	if(f1->IsOpen()) cout<<"TOF match efficiency file is open "<<endl;
	else cout<<"Fail to read TOF match efficiency file"<<endl;
	TString name;
	for (int i = 0;i<2;i++)
	{
		name = "PionTofMatch1D"+Flag[i]+"_reb";
		hTOFMatchEff[i] = (TH1D*)f1->Get(name);
		name = "EoverPi"+Flag[i];
		hElecPionRatio[i] = (TH1D*)f1->Get(name);
	}
	cout<<"read TOF match efficiency ok"<<endl;
	//get TPC tracking 
	TFile* f2 = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/Effinput/TPCEffHisto.root");
	if(f1->IsOpen()) cout<<"TPC match efficiency file is open "<<endl;
	else cout<<"Fail to read TPC match efficiency file"<<endl;
	for (int j = 0;j<nEta_TPC;j++)
	{
		hTPCTrackingEffPlus[j] = (TH1D*)f2->Get(Form("TPCEffEta%d",j));
		hTPCTrackingEffMinus[j] = (TH1D*)f2->Get(Form("TPCEffEta%d",j));
	}
	cout<<"read TPC match efficiency ok"<<endl;
	//get beta cut Eff
	TFile* f3 = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/Effinput/BetaCutEff.root");
	if(f1->IsOpen()) cout<<"beta cut efficiency file is open "<<endl;
	else cout<<"Fail to read beta cut efficiency file"<<endl;
	for (int i = 0;i<2;i++)
	{
		name = "betacutEffpT"+Flag[i];
		hBetaCutEff[i] = (TH1D*)f3->Get(name);
	}
	cout<<"read beta cut efficiency ok"<<endl;
	//get nSigmaE cut Eff
	TFile* f4 = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/Effinput/nSigmaECutEffhisto.root");
	if(f1->IsOpen()) 
	{
		cout<<"nSigmaE cut efficiency file is open "<<endl;
		hNSigmaECutEff[0] = (TH1D*)f4->Get("Cen0-80nSigmaEcutEffpTPlus");
		hNSigmaECutEff[1] = (TH1D*)f4->Get("Cen0-80nSigmaEcutEffpTMinus");
	}
	else cout<<"Fail to read nSgimaE cut efficiency file"<<endl;
	/*for (int i = 0;i<2;i++)
	{
		name = "nSigmaEcutEffpT"+Flag[i];
		hNSigmaECutEff[i] = (TH1D*)f4->Get(name);
	}*/
	//hNSigmaECutEff[0]->Draw();
	cout<<"read nSigmaE cut efficiency ok"<<endl;
}

double getEff(double pt, TH1D* histo , double ptLLimit, double ptHLimit)
{
	if (pt<ptLLimit) return 0;
	if (pt>ptHLimit) pt = ptHLimit;
	int Bin = histo->GetXaxis()->FindBin(pt+1.e-8);
	if (Bin<0) return 0;
	return histo->GetBinContent(Bin);
}
void getEtaPhiBin(double Eta, double Phi, int *EtaBin, int *PhiBin)//phi range is 0~2pi
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
////////

void InitializeEffFun(Int_t cenIdx)
{
	//****** Efficiency for single track ******
	funndEdxEffPt = new TF1("funndEdxEffPt",ndEdxEff,0.,10,1);
	funndEdxEffPt->SetParameter(0, cenIdx);
	funndEdxEffPt->SetNpx(1000);

	funnSigEffP = new TF1("funnSigEffP",nSigEEff,0.,10,1);
	funnSigEffP->SetParameter(0, cenIdx);
	funnSigEffP->SetNpx(1000);

	funBetaEffP = new TF1("funBetaEffP",betaEff,0.,10,1);
	funBetaEffP->SetParameter(0, cenIdx);
	funBetaEffP->SetNpx(1000);

	fEPiEff_Tof = new TF1("fEPiEff_Tof",EPiTofEffRatio,0.,10,1);
	fEPiEff_Tof->SetParameter(0, cenIdx);
	fEPiEff_Tof->SetNpx(1000);

    fMBEPiEff_Tof = new TF1("fMBEPiEff_Tof",EPiTofEffRatio,0.,10,1);
	fMBEPiEff_Tof->SetParameter(0, 0);
	fMBEPiEff_Tof->SetNpx(1000);

    fTpcEffRatioToMB = new TF1("fTpcEffRatioToMB",TpcEffRatioToMB,0.,10,1);
	fTpcEffRatioToMB->SetParameter(0, cenIdx);
	fTpcEffRatioToMB->SetNpx(1000);

    fTofEffRatioToMB = new TF1("fTofEffRatioToMB",TofEffRatioToMB,0.,10,1);
	fTofEffRatioToMB->SetParameter(0, cenIdx);
	fTofEffRatioToMB->SetNpx(1000);

	char plusname[100];
	char minusname[100];

	TFile *fTpcTrkEff = TFile::Open("/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/inputRootFiles/ElectronTpcTrackEff4AllCens.root");	
	for(int i=0;i<NEta_TPC[cenIdx];i++) {
		for(int j=0;j<NPhi_TPC[cenIdx];j++) {
			sprintf(plusname, "hEPlusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[cenIdx], CentralityHi[cenIdx], i, j);
			sprintf(minusname, "hEMinusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[cenIdx], CentralityHi[cenIdx], i, j);

			hEff_Tpc_Pos[i][j] = (TH1D *)fTpcTrkEff->Get(plusname);
			hEff_Tpc_Neg[i][j] = (TH1D *)fTpcTrkEff->Get(minusname);
		}
	}
	//MB TPC efficiency
	for(int i=0;i<NEta_TPC[0];i++) {
		for(int j=0;j<NPhi_TPC[0];j++) {
			sprintf(plusname, "hEPlusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[0], CentralityHi[0], i, j);
			sprintf(minusname, "hEMinusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[0], CentralityHi[0], i, j);

			hMBEff_Tpc_Pos[i][j] = (TH1D *)fTpcTrkEff->Get(plusname);
			hMBEff_Tpc_Neg[i][j] = (TH1D *)fTpcTrkEff->Get(minusname);
		}
	}

	TFile *fTofMatchEff = TFile::Open("/star/u/syang/run12/uu/minibias/eff/pairEff/cocktails_new/inputRootFiles/PionTofMatchEff4AllCens.root"); 
	for(int i=0;i<NEta_TOF[cenIdx];i++) {
		for(int j=0;j<NPhi_TOF[cenIdx];j++) {
			sprintf(plusname, "hPiPlusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[cenIdx], CentralityHi[cenIdx], i, j);
			sprintf(minusname, "hPiMinusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[cenIdx], CentralityHi[cenIdx], i, j);

			hEff_Tof_Pos[i][j] = (TH1D *)fTofMatchEff->Get(plusname);
			hEff_Tof_Neg[i][j] = (TH1D *)fTofMatchEff->Get(minusname);
		}
	}
	//MB TOF efficiency
	for(int i=0;i<NEta_TOF[0];i++) {
		for(int j=0;j<NPhi_TOF[0];j++) {
			sprintf(plusname, "hPiPlusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[0], CentralityHi[0], i, j);
			sprintf(minusname, "hPiMinusEff_Cen%d_%d_Eta%dPhi%d", CentralityLow[0], CentralityHi[0], i, j);

			hMBEff_Tof_Pos[i][j] = (TH1D *)fTofMatchEff->Get(plusname);
			hMBEff_Tof_Neg[i][j] = (TH1D *)fTofMatchEff->Get(minusname);
		}
	}
}


void tpcPtEtaPhi2Bin(int icen, double pt, double eta, double phi, int *ipt, int *ieta, int *iphi)
{   
	*ipt = 0; *ieta = 0; *iphi = 0;

	*ipt = (int)(pt/(4./NPt_TPC[icen]));
	if(*ipt<0) *ipt = 0;
	if(*ipt>=NPt_TPC[icen]) *ipt = NPt_TPC[icen]-1;

	*ieta = (int)((eta+1.)/(2./NEta_TPC[icen]));
	if(*ieta<0) *ieta = 0;
	if(*ieta>=NEta_TPC[icen]) *ieta = NEta_TPC[icen]-1;

	if(phi<-TMath::Pi()+TMath::Pi()/12) phi += 2*TMath::Pi();
	*iphi = (int)((phi+TMath::Pi()-TMath::Pi()/12)/(2*TMath::Pi()/NPhi_TPC[icen]));
	if(*iphi<0) *iphi = 0;
	if(*iphi>=NPhi_TPC[icen]) *iphi = NPhi_TPC[icen]-1;
	return;
}

void tofPtEtaPhi2Bin(int icen, double pt, double eta, double phi, int *ipt, int *ieta, int *iphi)
{   
	*ipt = 0; *ieta = 0; *iphi = 0;

	*ipt = (int)(pt/(4./NPt_TOF[icen]));
	if(*ipt<0) *ipt = 0;
	if(*ipt>=NPt_TOF[icen]) *ipt = NPt_TOF[icen]-1;

	*ieta = (int)((eta+1.)/(2./NEta_TOF[icen]));
	if(*ieta<0) *ieta = 0;
	if(*ieta>=NEta_TOF[icen]) *ieta = NEta_TOF[icen]-1;

	if(phi<-TMath::Pi()+TMath::Pi()/12) phi += 2*TMath::Pi();
	*iphi = (int)((phi+TMath::Pi()-TMath::Pi()/12)/(2*TMath::Pi()/NPhi_TOF[icen]));
	if(*iphi<0) *iphi = 0;
	if(*iphi>=NPhi_TOF[icen]) *iphi = NPhi_TOF[icen]-1;
	return;
}

//-------------------------------------------------------
// return values
// Within Acceptance
//   0 - east arm
//   1 - west arm
// Outside Acceptance
//  -1 - invalid charge
//  -2 - fails minimum pT cut (200 MeV)
//  -3 - fails polar angle cut (this is an event vertex dependent rapdity cut)
//  -4 - fails azimuthal angle cut
int PHENIXFilter(int charge, double px, double py, double pz, double evt_zvtx)
{
	// Phi acceptance west arm
	double phi_bot_w = -0.589;  // ~ -3*pi/16
	double phi_top_w =  0.982;  // ~ 5*pi/16
	// Phi acceptance east arm
	double phi_bot_e =  3.731;  // ~ 11*pi/16
	double phi_top_e =  2.160;  // ~ 19*pi/16

	// Theta acceptance at evt_zvtx=0
	const double theta_min = 1.23;
	const double theta_max = 1.92;

	// slopes for each line tuned ++ field (Run-4/5/6)
	double sl1 =  0.309;    // CRK cut
	double sl2 =  0.206;    // DCH cut
	double slz =  0.004224; // Z_theta slope

	if (charge==0) return -1;

	const double pt = hypot(px, py);
	if (pt < 0.2) return -2;  //Lower pt cut


	// --------------- ACCEPTANCE CHECK -------------------------

	// --------------------- Theta acceptance ----------------------
	const double theta  = atan2(pt,pz);
	const double th1 = theta-evt_zvtx*slz;
	if (th1>theta_max || th1<theta_min) return -3;

	// --------------------- Phi acceptance ----------------------
	const double q_pt = charge/pt;
	const double phi  = atan2(py,px);
	const double phi1 = phi-q_pt*sl1;
	const double phi2 = phi-q_pt*sl2;

	if (phi2>phi_top_e && phi2<phi_bot_e &&
			phi1>phi_top_e && phi1<phi_bot_e) // Fall into East phi acceptance
	{
		return 0;
	}
	else if (phi2>phi_bot_w && phi2<phi_top_w &&
			phi1>phi_bot_w && phi1<phi_top_w)  // Fall into West phi acceptance
	{
		return 1;
	}

	// bad track
	return -4;
}

#endif
