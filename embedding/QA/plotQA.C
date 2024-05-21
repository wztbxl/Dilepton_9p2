#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TAttFill.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TLegend.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPDF.h"

#include "cuts.h"

using namespace::std;

//***** constrain TPC sector7 geometry *****
TF1 *funPosHi;
TF1 *funPosLow;
TF1 *funNegHi;
TF1 *funNegLow;
TF1 *funNegHiMirror;
TF1 *funNegLowMirror;
Double_t par[4][4];
Double_t parErr[4][4];

bool Init();
TH1F *calRatio(TH1F *hNum, TH1F *hDen, TString name);
TH2D* histo(TString name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, TString xTitle, TString yTitle);
TLatex* drawLatex(Double_t x, Double_t y, TString text, Int_t textFont, Double_t textSize, Int_t colorIndex);
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth,Int_t lineStyle , Int_t lineColor);
void drawLines(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth,Int_t lineStyle, Int_t lineColor);
void setPad(float left, float right, float top, float bottom);
void clearPad(TCanvas *c, Int_t nPads);
void setHisto(TH1F *h,Int_t MarkerStyle, Float_t MarkerSize, Int_t MarkerColor,Int_t LineColor);
void setHisto(TH1D *h,Int_t MarkerStyle, Float_t MarkerSize, Int_t MarkerColor,Int_t LineColor);
void pdfAction(TCanvas *c, TPDF *ps);

void plotQA(TString inFile="minibias"){

	cout<<"Start to draw the embedding QA plots ..."<<endl;
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(0);

	Int_t binLow,binHi;
	Int_t binLow1,binHi1;
	Float_t textX = 0.55;
	Float_t textY = 0.75;
	Int_t textFont = 42;
	Float_t textSize = 0.08;
	Int_t textColor = 1;

	//const Int_t nMomBins = 18;
	const int nMomBins = 13;
	Float_t mom[nMomBins+1] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.5,2.5,4.0};
	//const Int_t nPiMomBins = 27;
	const int nPiMomBins = 27;
	Float_t pimom[nPiMomBins+1] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,4.0};
	//const Int_t nMomBins = 8;
	//Float_t mom[nMomBins+1] = {0.2,0.5,0.8,1.1,1.4,1.7,2.0,2.5,4.0};
	//const Int_t nPiMomBins = 9;
	//Float_t pimom[nPiMomBins+1] = {0.2,0.5,0.8,1.1,1.4,1.7,2.0,2.5,3.0,4.0};
    const Int_t NumEtaBins = 20;
    Float_t eta[NumEtaBins + 1] = {-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	const Int_t NumPhiBins = 20;
    double PhiRange = TMath::Pi()/10;

	Int_t momBinLow,momBinHi;

	//TFile *fMc = new TFile( "./mc/test/test.QAhisto.root" );
	TFile *fMc = new TFile( Form("./%sMC.QAhistoFullsample.root",inFile.Data()) );
	if(!fMc->IsOpen()){
		cout<<"Fail to open the MC QA root file!";
		return;
	}
	TFile *fData = new TFile( Form("./%sData.QAhistoFullsample.root",inFile.Data()) );
	if(!fData->IsOpen()){
		cout<<"Fail to open the DATA QA root file!";
		return;
	}

	if( !Init() ){
		cout<<"Fail to initialize !"<<endl;
		return;
	}

	//embedding 
	TH2F *hVyvsVx = (TH2F *)fMc->Get("hVyvsVx");
	TH2F *hVyvsVz = (TH2F *)fMc->Get("hVyvsVz");
	TH1F *hVzDiff = (TH1F *)fMc->Get("hVzDiff");
	TH2F *hRcMcVxDiffvsMcVx = (TH2F *)fMc->Get("hRcMcVxDiffvsMcVx");
	TH2F *hRcMcVyDiffvsMcVy = (TH2F *)fMc->Get("hRcMcVyDiffvsMcVy");
	TH2F *hRcMcVzDiffvsMcVz = (TH2F *)fMc->Get("hRcMcVzDiffvsMcVz");
	TH2F *hRefMultvsRefMultCorr = (TH2F *)fMc->Get("hRefMultvsRefMultCorr");
	TH2F *hNMatchTrksvsInputTrk = (TH2F *)fMc->Get("hNMatchTrksvsInputTrks");
	TH3F *hRcPtvsMcPtQ = (TH3F *)fMc->Get("hRcPtvsMcPtQ");
	TH3F *hRcEtavsMcEtaQ = (TH3F *)fMc->Get("hRcEtavsMcEtaQ");
	TH3F *hRcPhivsMcPhiQ = (TH3F *)fMc->Get("hRcPhivsMcPhiQ");
	TH3F *hMcEtavsPtQ = (TH3F *)fMc->Get("hMcEtavsPtQ");
	hMcEtavsPtQ->Sumw2();
	TH3F *hMcPhivsPtQ = (TH3F *)fMc->Get("hMcPhivsPtQ");
	hMcPhivsPtQ->Sumw2();
	TH3F *hRcEtavsPtQ = (TH3F *)fMc->Get("hRcEtavsPtQ");
	hRcEtavsPtQ->Sumw2();
	TH3F *hRcPhivsPtQ = (TH3F *)fMc->Get("hRcPhivsPtQ");
	hRcPhivsPtQ->Sumw2();
	TH3F *hPtResvsPtQ = (TH3F *)fMc->Get("hPtResvsPtQ");
	hPtResvsPtQ->Sumw2();
	TH3F *hPResvsPQ = (TH3F *)fMc->Get("hPResvsPQ");
	hPResvsPQ->Sumw2();
	TH3F *hRcNHitsFitvsPtQ = (TH3F *)fMc->Get("hRcNHitsFitvsPtQ");
	hRcNHitsFitvsPtQ->Sumw2();
	TH3F *hRcNHitsFitvsEtaQ = (TH3F *)fMc->Get("hRcNHitsFitvsEtaQ");
	hRcNHitsFitvsEtaQ->Sumw2();
	TH3F *hRcNHitsFitvsPhiQ = (TH3F *)fMc->Get("hRcNHitsFitvsPhiQ");
	hRcNHitsFitvsPhiQ->Sumw2();
	TH3F *hRcNHitsPossvsPtQ = (TH3F *)fMc->Get("hRcNHitsPossvsPtQ");
	hRcNHitsPossvsPtQ->Sumw2();
	TH3F *hRcNHitsPossvsEtaQ = (TH3F *)fMc->Get("hRcNHitsPossvsEtaQ");
	hRcNHitsPossvsEtaQ->Sumw2();
	TH3F *hRcNHitsPossvsPhiQ = (TH3F *)fMc->Get("hRcNHitsPossvsPhiQ");
	hRcNHitsPossvsPhiQ->Sumw2();
	TH3F *hRcNHitsDedxvsPtQ = (TH3F *)fMc->Get("hRcNHitsDedxvsPtQ");
	hRcNHitsDedxvsPtQ->Sumw2();
	TH3F *hRcNHitsCommonvsPtQ = (TH3F *)fMc->Get("hRcNHitsCommonvsPtQ");
	hRcNHitsCommonvsPtQ->Sumw2();
	TH3F *hRcNHitsCommonvsEtaQ = (TH3F *)fMc->Get("hRcNHitsCommonvsEtaQ");
	hRcNHitsCommonvsEtaQ->Sumw2();
	TH3F *hRcNHitsCommonvsPhiQ = (TH3F *)fMc->Get("hRcNHitsCommonvsPhiQ");
	hRcNHitsCommonvsPhiQ->Sumw2();
	TH3F *hRcDcavsPtQ = (TH3F *)fMc->Get("hRcDcavsPtQ");
	hRcDcavsPtQ->Sumw2();
	TH3F *hRcDedxvsPQ = (TH3F *)fMc->Get("hRcDedxvsPQ");
	hRcDedxvsPQ->Sumw2();
	TH3F *hRcNSigmaEvsPQ = (TH3F *)fMc->Get("hRcNSigmaEvsPQ");
	hRcNSigmaEvsPQ->Sumw2();
	
	cout<<1<<endl;
	TH3F *hRcNHitsFitPtEtaMinus = (TH3F* )fMc->Get("hPENitsFitsvsPtvsEtaMinus");
	hRcNHitsFitPtEtaMinus->Sumw2();
	TH3F *hRcNHitsFitPtEtaPlus = (TH3F* )fMc->Get("hPENitsFitsvsPtvsEtaPlus");
	hRcNHitsFitPtEtaPlus->Sumw2();
	TH3F *hRcNHitsFitPtPhiMinus = (TH3F* )fMc->Get("hPENitsFitsvsPtvsPhiMinus");
	hRcNHitsFitPtPhiMinus->Sumw2();
	TH3F *hRcNHitsFitPtPhiPlus = (TH3F* )fMc->Get("hPENitsFitsvsPtvsPhiPlus");
	hRcNHitsFitPtPhiPlus->Sumw2();
	TH3F *hRcNHitsDedxPtEtaMinus = (TH3F* )fMc->Get("hPENHitsDedxvsPtvsEtaMinus");
	hRcNHitsDedxPtEtaMinus->Sumw2();
	TH3F *hRcNHitsDedxPtEtaPlus = (TH3F* )fMc->Get("hPENHitsDedxvsPtvsEtaPlus");
	hRcNHitsDedxPtEtaPlus->Sumw2();
	TH3F *hRcDcaPtEtaMinus = (TH3F* )fMc->Get("hPEDcavsPtvsEtaMinus");
	hRcDcaPtEtaMinus->Sumw2();
	TH3F *hRcDcaPtEtaPlus = (TH3F* )fMc->Get("hPEDcavsPtvsEtaPlus");
	hRcDcaPtEtaPlus->Sumw2();
	
	//real data
	TH2F *hULMvsPt = (TH2F *)fData->Get("hULMvsPt");
	hULMvsPt->Sumw2();
	TH2F *hLPosMvsPt = (TH2F *)fData->Get("hLPosMvsPt");
	hLPosMvsPt->Sumw2();
	TH2F *hLNegMvsPt = (TH2F *)fData->Get("hLNegMvsPt");
	hLNegMvsPt->Sumw2();
	TH3F *hDataEtavsPtQ = (TH3F *)fData->Get("hEtavsPtQ");
	hDataEtavsPtQ->Sumw2();
	TH3F *hDataPhivsPtQ = (TH3F *)fData->Get("hPhivsPtQ");
	hDataPhivsPtQ->Sumw2();
	TH3F *hDataNHitsFitvsPtQ = (TH3F *)fData->Get("hNHitsFitvsPtQ");
	hDataNHitsFitvsPtQ->Sumw2();
	TH3F *hDataNHitsFitvsEtaQ = (TH3F *)fData->Get("hNHitsFitvsEtaQ");
	hDataNHitsFitvsEtaQ->Sumw2();
	TH3F *hDataNHitsFitvsPhiQ = (TH3F *)fData->Get("hNHitsFitvsPhiQ");
	hDataNHitsFitvsPhiQ->Sumw2();
	TH3F *hDataNHitsPossvsPtQ = (TH3F *)fData->Get("hNHitsPossvsPtQ");
	hDataNHitsPossvsPtQ->Sumw2();
	TH3F *hDataNHitsPossvsEtaQ = (TH3F *)fData->Get("hNHitsPossvsEtaQ");
	hDataNHitsPossvsEtaQ->Sumw2();
	TH3F *hDataNHitsPossvsPhiQ = (TH3F *)fData->Get("hNHitsPossvsPhiQ");
	hDataNHitsPossvsPhiQ->Sumw2();
	TH3F *hDataNHitsDedxvsPtQ = (TH3F *)fData->Get("hNHitsDedxvsPtQ");
	hDataNHitsDedxvsPtQ->Sumw2();
	TH3F *hDataDcavsPtQ = (TH3F *)fData->Get("hDcavsPtQ");
	hDataDcavsPtQ->Sumw2();
	TH3F *hDataDedxvsPQ = (TH3F *)fData->Get("hDedxvsPQ");
	hDataDedxvsPQ->Sumw2();
	TH3F *hDataNSigmaEvsPQ = (TH3F *)fData->Get("hNSigmaEvsPQ");
	hDataNSigmaEvsPQ->Sumw2();

	TH3F *hDataNHitsFitPtEtaMinus = (TH3F* )fData->Get("hPENitsFitsvsPtvsEtaMinus");
	hDataNHitsFitPtEtaMinus->Sumw2();
	TH3F *hDataNHitsFitPtEtaPlus = (TH3F* )fData->Get("hPENitsFitsvsPtvsEtaPlus");
	hDataNHitsFitPtEtaPlus->Sumw2();
	TH3F *hDataNHitsFitPtPhiMinus = (TH3F* )fData->Get("hPENitsFitsvsPtvsPhiMinus");
	hDataNHitsFitPtPhiMinus->Sumw2();
	TH3F *hDataNHitsFitPtPhiPlus = (TH3F* )fData->Get("hPENitsFitsvsPtvsPhiPlus");
	hDataNHitsFitPtPhiPlus->Sumw2();

	TH3F *hDataPiEtavsPtQ = (TH3F *)fData->Get("hPiEtavsPtQ");
	hDataPiEtavsPtQ->Sumw2();
	TH3F *hDataPiPhivsPtQ = (TH3F *)fData->Get("hPiPhivsPtQ");
	hDataPiPhivsPtQ->Sumw2();
	TH3F *hDataPiNHitsFitvsPtQ = (TH3F *)fData->Get("hPiNHitsFitvsPtQ");
	hDataPiNHitsFitvsPtQ->Sumw2();
	TH3F *hDataPiNHitsFitvsEtaQ = (TH3F *)fData->Get("hPiNHitsFitvsEtaQ");
	hDataPiNHitsFitvsEtaQ->Sumw2();
	TH3F *hDataPiNHitsFitvsPhiQ = (TH3F *)fData->Get("hPiNHitsFitvsPhiQ");
	hDataPiNHitsFitvsPhiQ->Sumw2();
	TH3F *hDataPiNHitsPossvsPtQ = (TH3F *)fData->Get("hPiNHitsPossvsPtQ");
	hDataPiNHitsPossvsPtQ->Sumw2();
	TH3F *hDataPiNHitsPossvsEtaQ = (TH3F *)fData->Get("hPiNHitsPossvsEtaQ");
	hDataPiNHitsPossvsEtaQ->Sumw2();
	TH3F *hDataPiNHitsPossvsPhiQ = (TH3F *)fData->Get("hPiNHitsPossvsPhiQ");
	hDataPiNHitsPossvsPhiQ->Sumw2();
	TH3F *hDataPiNHitsDedxvsPtQ = (TH3F *)fData->Get("hPiNHitsDedxvsPtQ");
	hDataPiNHitsDedxvsPtQ->Sumw2();
	TH3F *hDataPiDcavsPtQ = (TH3F *)fData->Get("hPiDcavsPtQ");
	hDataPiDcavsPtQ->Sumw2();
	TH3F *hDataPiDedxvsPQ = (TH3F *)fData->Get("hPiDedxvsPQ");
	hDataPiDedxvsPQ->Sumw2();
	TH3F *hDataPiNSigmaEvsPQ = (TH3F *)fData->Get("hPiNSigmaEvsPQ");
	hDataPiNSigmaEvsPQ->Sumw2();

	TCanvas *c = new TCanvas("c", "c",0,0,800,600);

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,1000,750);
	TPDF *ps = new TPDF(Form("plots/%s.embeddingQA.pdf",inFile.Data()),111);
	ps->Off();

	Int_t nColumns = 2;
	Int_t nRaws = 2;
	Int_t nPads = nColumns*nRaws;

	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hVyvsVx->Draw("colz");
	c1->cd(2);
	gPad->SetLogz(1);
	hVyvsVz->Draw("colz");
	c1->cd(3);
	gPad->SetLogy(0);
	TH1F *hVz = (TH1F *)hVyvsVz->ProjectionX("hVz");
	hVz->Draw();
	c1->cd(4);
	gPad->SetLogy(0);
	hVzDiff->GetXaxis()->SetRangeUser(-5,5);
	hVzDiff->Draw();
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	TH1F *hRcMcVxDiff = (TH1F *)hRcMcVxDiffvsMcVx->ProjectionY("hRcMcVxDiff");
	hRcMcVxDiff->GetXaxis()->SetRangeUser(-1,1);
	hRcMcVxDiff->Draw();
	c1->cd(2);
	TH1F *hRcMcVyDiff = (TH1F *)hRcMcVyDiffvsMcVy->ProjectionY("hRcMcVyDiff");
	hRcMcVyDiff->GetXaxis()->SetRangeUser(-1,1);
	hRcMcVyDiff->Draw();
	c1->cd(3);
	TH1F *hRcMcVzDiff = (TH1F *)hRcMcVzDiffvsMcVz->ProjectionY("hRcMcVzDiff");
	hRcMcVzDiff->GetXaxis()->SetRangeUser(-1,1);
	hRcMcVzDiff->Draw();
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRefMultvsRefMultCorr->Draw("colz");
	c1->cd(2);
	gPad->SetLogz(1);
	hNMatchTrksvsInputTrk->GetXaxis()->SetTitle("#McTracks (pGeantID==0)");
	hNMatchTrksvsInputTrk->GetXaxis()->SetRangeUser(0,40);
	hNMatchTrksvsInputTrk->GetYaxis()->SetTitle("#RcTracks (primary && pGeantID==0)");
	hNMatchTrksvsInputTrk->GetYaxis()->SetRangeUser(0,40);
	hNMatchTrksvsInputTrk->Draw("colz");
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	hMcEtavsPtQ->GetYaxis()->SetRangeUser(0,5.0);
	hMcEtavsPtQ->GetZaxis()->SetRangeUser(-1.5,1.5);
	hMcPhivsPtQ->GetYaxis()->SetRangeUser(0,5.0);
	c1->cd(1);
	gPad->SetLogz(1);
	hMcEtavsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hMcEtavsPtE = (TH2F *)hMcEtavsPtQ->Project3D("ZY"); 
	hMcEtavsPtE->SetNameTitle("hMcEtavsPtE","hMcEtavsPtE");
	hMcEtavsPtE->Draw("colz");
	drawLatex(textX,textY,"Emb. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hMcEtavsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hMcEtavsPtP = (TH2F *)hMcEtavsPtQ->Project3D("ZY");
	hMcEtavsPtP->SetNameTitle("hMcEtavsPtP","hMcEtavsPtP");
	hMcEtavsPtP->Draw("colz");
	drawLatex(textX,textY,"Emb. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hMcPhivsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hMcPhivsPtE = (TH2F *)hMcPhivsPtQ->Project3D("ZY"); 
	hMcPhivsPtE->SetNameTitle("hMcPhivsPtE","hMcPhivsPtE");
	hMcPhivsPtE->Draw("colz");
	drawLatex(textX,textY,"Emb. MC e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hMcPhivsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hMcPhivsPtP = (TH2F *)hMcPhivsPtQ->Project3D("ZY");
	hMcPhivsPtP->SetNameTitle("hMcPhivsPtP","hMcPhivsPtP");
	hMcPhivsPtP->Draw("colz");
	drawLatex(textX,textY,"Emb. MC e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
	leg->SetFillColor(10);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);

	clearPad(c1,nPads);
	c1->cd(1);
	TH1F *hMcPtE = (TH1F *)hMcEtavsPtE->ProjectionX("hMcPtE");
	TH1F *hMcPtP = (TH1F *)hMcEtavsPtP->ProjectionX("hMcPtP");
	hMcPtE->SetLineColor(1);
	hMcPtP->SetLineColor(2);
	hMcPtP->SetMinimum(1.);
	hMcPtE->Draw("histe");
	hMcPtP->Draw("histesame");
	leg->AddEntry(hMcPtE,"Emb. MC e^{-}","pl");
	leg->AddEntry(hMcPtP,"Emb. MC e^{+}","pl");
	leg->Draw("same");
	c1->cd(2);
	TH1F *hMcEtaE = (TH1F *)hMcEtavsPtE->ProjectionY("hMcEtaE");
	TH1F *hMcEtaP = (TH1F *)hMcEtavsPtP->ProjectionY("hMcEtaP");
	hMcEtaE->SetLineColor(1);
	hMcEtaP->SetLineColor(2);
	hMcEtaP->SetMinimum(1.);
	hMcEtaE->Draw("histe");
	hMcEtaP->Draw("histesame");
	leg->Draw("same");
	c1->cd(3);
	TH1F *hMcPhiE = (TH1F *)hMcPhivsPtE->ProjectionY("hMcPhiE");
	TH1F *hMcPhiP = (TH1F *)hMcPhivsPtP->ProjectionY("hMcPhiP");
	hMcPhiE->SetLineColor(1);
	hMcPhiP->SetLineColor(2);
	hMcPhiP->SetMinimum(1.);
	hMcPhiE->Draw("histe");
	hMcPhiP->Draw("histesame");
	leg->Draw("same");
	pdfAction(c1,ps);

	c1->Clear();
	c1->cd();
	drawLatex(0.08,0.9,"Compare Embedded Electron and Positron with Data",62,0.05,1);
	drawLatex(0.08,0.8,"Event Level Cuts:",62,0.05,1);
	drawLatex(0.15,0.72,Form("|V_{x}| > 0  ||  |V_{y}| > 0  ||  |V_{z}|> 0"),42,0.05,1);
	drawLatex(0.15,0.62,Form("#sqrt{V_{x}^{2}+V_{y}^{2}} < %4.1f cm",mVrCut),42,0.05,1);
	drawLatex(0.15,0.55,Form("|V_{z}| < %4.1f cm",mVzCut),42,0.05,1);
	drawLatex(0.15,0.48,Form("|VpdV_{z}-TpcV_{z}| < %4.1f cm",mVzDiffCut),42,0.05,1);
	drawLatex(0.08,0.40,"Track Level Cuts:",62,0.05,1);
	drawLatex(0.15,0.32,Form("p_{T} >= %4.1f GeV/c",mTpcePtCut[0]),42,0.05,1);
	drawLatex(0.15,0.27,Form("|#eta| <= %4.1f",mTpceEtaCut),42,0.05,1);
	drawLatex(0.42,0.32,Form("nHitsFit >= %d",15),42,0.05,1);
	drawLatex(0.42,0.27,Form("ratio >= %4.2f",mTpceNHitsFitRatioCut),42,0.05,1);
	drawLatex(0.69,0.32,Form("nHitsDedx >= %d",mTpceNHitsDedxCut),42,0.05,1);
	drawLatex(0.69,0.27,Form("dca <= %4.1f cm",2.5),42,0.05,1);
	drawLatex(0.10,0.20,"Data:",62,0.05,1);
	drawLatex(0.15,0.14,"For e:",62,0.05,1);
	drawLatex(0.28,0.14,"1.5625*p-2.5625<n#sigma_{e}<1.66(p<1 GeV/c) -1.09<n#sigma_{e}<1.66 (p>=1 GeV/c)",42,0.02,1);
	drawLatex(0.75,0.14,Form("|1-1/#beta| <= %4.3f",mTpceBeta2TOFCut),42,0.05,1);
	drawLatex(0.15,0.09,"For #pi:",62,0.05,1);
	drawLatex(0.35,0.09,"|n#sigma_{#pi}| <= 2.0",42,0.05,1);
	drawLatex(0.60,0.09,"|m^{2}_{#pi} - 0.019| <= 0.003",42,0.05,1);
	drawLatex(0.10,0.02,"Embedding: ",62,0.05,1);
	drawLatex(0.35,0.02,"parent GeantId == 0",42,0.05,1);
	drawLatex(0.75,0.02,Form("|#eta| <= %4.1f",mTpceEtaCut),42,0.05,1);
	pdfAction(c1,ps);

	hRcPtvsMcPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcPtvsMcPtE = (TH2F *)hRcPtvsMcPtQ->Project3D("ZY");
	hRcPtvsMcPtE->SetNameTitle("hRcPtvsMcPtE","hRcPtvsMcPtE");
	hRcPtvsMcPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcPtvsMcPtP = (TH2F *)hRcPtvsMcPtQ->Project3D("ZY");
	hRcPtvsMcPtP->SetNameTitle("hRcPtvsMcPtP","hRcPtvsMcPtP");
	hRcEtavsMcEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcEtavsMcEtaE = (TH2F *)hRcEtavsMcEtaQ->Project3D("ZY");
	hRcEtavsMcEtaE->SetNameTitle("hRcEtavsMcEtaE","hRcEtavsMcEtaE");
	hRcEtavsMcEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcEtavsMcEtaP = (TH2F *)hRcEtavsMcEtaQ->Project3D("ZY");
	hRcEtavsMcEtaP->SetNameTitle("hRcEtavsMcEtaP","hRcEtavsMcEtaP");
	hRcPhivsMcPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcPhivsMcPhiE = (TH2F *)hRcPhivsMcPhiQ->Project3D("ZY");
	hRcPhivsMcPhiE->SetNameTitle("hRcPhivsMcPhiE","hRcPhivsMcPhiE");
	hRcPhivsMcPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcPhivsMcPhiP = (TH2F *)hRcPhivsMcPhiQ->Project3D("ZY");
	hRcPhivsMcPhiP->SetNameTitle("hRcPhivsMcPhiP","hRcPhivsMcPhiP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcPtvsMcPtE->GetXaxis()->SetRangeUser(0,5);
	hRcPtvsMcPtE->GetYaxis()->SetRangeUser(0,5);
	hRcPtvsMcPtE->Draw("colz");
	drawLatex(0.20,0.75,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcPtvsMcPtP->GetXaxis()->SetRangeUser(0,5);
	hRcPtvsMcPtP->GetYaxis()->SetRangeUser(0,5);
	hRcPtvsMcPtP->Draw("colz");
	drawLatex(0.20,0.75,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hRcEtavsMcEtaE->GetXaxis()->SetRangeUser(-1.5,1.5);
	hRcEtavsMcEtaE->GetYaxis()->SetRangeUser(-1.5,1.5);
	hRcEtavsMcEtaE->Draw("colz");
	drawLatex(0.20,0.75,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hRcEtavsMcEtaP->GetXaxis()->SetRangeUser(-1.5,1.5);
	hRcEtavsMcEtaP->GetYaxis()->SetRangeUser(-1.5,1.5);
	hRcEtavsMcEtaP->Draw("colz");
	drawLatex(0.20,0.75,"Rec. MC e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcPhivsMcPhiE->Draw("colz");
	drawLatex(0.20,0.75,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcPhivsMcPhiP->Draw("colz");
	drawLatex(0.20,0.75,"Rec. MC e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	//phe quality histograms
	c1->cd();
	TH1F *hULM = (TH1F *)hULMvsPt->ProjectionY("hULM");
	TH1F *hLPosM = (TH1F *)hLPosMvsPt->ProjectionY("hLPosM");
	TH1F *hLNegM = (TH1F *)hLNegMvsPt->ProjectionY("hLNegM");
	TH1F *hBg = (TH1F *)hLPosM->Clone("hBg");
	hBg->Add(hLNegM);
	TH1F *hSig = (TH1F *)hULM->Clone("hSig");
	hSig->Add(hBg,-1);
	setHisto(hULM,20,0.5,1,1);
	setHisto(hBg,24,0.5,2,2);
	setHisto(hSig,20,0.5,4,4);
	//Int_t xBinLow = 1;
	Int_t xBinLow = hULM->GetXaxis()->FindBin(0.0+1.e-4);
	Int_t xBinHi = hULM->GetXaxis()->FindBin(0.015-1.e-4);
	Float_t Num = hSig->Integral(xBinLow,xBinHi);
	Float_t Den = hBg->Integral(xBinLow,xBinHi);
	hULM->GetXaxis()->SetRangeUser(0,0.07);
	hULM->Draw("pe");
	hBg->Draw("pesame");
	hSig->Draw("pesame");
	//drawLine(0.006,hULM->GetMinimum(),0.006,hSig->GetMaximum(),2,2,kViolet);
	drawLine(0.015,hULM->GetMinimum(),0.015,hSig->GetMaximum(),2,2,kViolet);
	drawLatex(0.3,0.75,Form("S/B = %4.2f:1",Num/Den),42,0.05,1);
	pdfAction(c1,ps);

	hPtResvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hPtResvsPtE = (TH2F *)hPtResvsPtQ->Project3D("ZY");
	hPtResvsPtE->SetNameTitle("hPtResvsPtE","hPtResvsPtE");
	hPtResvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hPtResvsPtP = (TH2F *)hPtResvsPtQ->Project3D("ZY");
	hPtResvsPtP->SetNameTitle("hPtResvsPtP","hPtResvsPtP");

	TF1 *funPtRes = new TF1("funPtRes","sqrt([0]*[0]*x*x+[1]*[1])",0.,10.);
	funPtRes->SetParameters(0.0036,0.008);

	c1->Clear();
	c1->Divide(nColumns,nRaws);
	nPads = nColumns*nRaws;
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hPtResvsPtE->SetNameTitle("hPtResvsPtE","hPtResvsPtE");
	hPtResvsPtE->RebinX(10);
	hPtResvsPtE->GetXaxis()->SetRangeUser(0,4.5);
	hPtResvsPtE->GetYaxis()->SetRangeUser(-0.2,0.2);
	hPtResvsPtE->DrawClone("colz");
	hPtResvsPtE->FitSlicesY();
	TH1F *hMeanPtResvsPtE = (TH1F *)gDirectory->Get("hPtResvsPtE_1");
	setHisto(hMeanPtResvsPtE,20,0.8,1,1);
	hMeanPtResvsPtE->Draw("same");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	TH1F *hSigmaPtResvsPtE = (TH1F *)gDirectory->Get("hPtResvsPtE_2");
	hSigmaPtResvsPtE->GetXaxis()->SetRangeUser(0,4.5);
	hSigmaPtResvsPtE->GetYaxis()->SetRangeUser(0,0.025);
	hSigmaPtResvsPtE->GetYaxis()->SetTitle("#sigma_{p_{T}}/p_{T}^{MC}");
	setHisto(hSigmaPtResvsPtE,20,0.8,1,1);
	hSigmaPtResvsPtE->Fit("funPtRes","R","",0.2,4.0);
	hSigmaPtResvsPtE->Draw("pesame");
	drawLatex(textX-0.2,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hPtResvsPtP->SetNameTitle("hPtResvsPtP","hPtResvsPtP");
	hPtResvsPtP->RebinX(10);
	hPtResvsPtP->GetXaxis()->SetRangeUser(0,4.5);
	hPtResvsPtP->GetYaxis()->SetRangeUser(-0.2,0.2);
	hPtResvsPtP->DrawClone("colz");
	hPtResvsPtP->FitSlicesY();
	TH1F *hMeanPtResvsPtP = (TH1F *)gDirectory->Get("hPtResvsPtP_1");
	setHisto(hMeanPtResvsPtP,20,0.8,1,1);
	hMeanPtResvsPtP->Draw("same");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(4);
	TH1F *hSigmaPtResvsPtP = (TH1F *)gDirectory->Get("hPtResvsPtP_2");
	hSigmaPtResvsPtP->GetXaxis()->SetRangeUser(0,4.5);
	hSigmaPtResvsPtP->GetYaxis()->SetRangeUser(0,0.025);
	hSigmaPtResvsPtP->GetYaxis()->SetTitle("#sigma_{p_{T}}/p_{T}^{MC}");
	setHisto(hSigmaPtResvsPtP,20,0.8,1,1);
	hSigmaPtResvsPtP->Fit("funPtRes","R","",0.2,4.0);
	hSigmaPtResvsPtP->Draw("pesame");
	drawLatex(textX-0.2,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	c->cd();
	setPad(0.13, 0.07, 0.08, 0.12);
	hSigmaPtResvsPtP->GetYaxis()->SetTitleOffset(1.2);
	hSigmaPtResvsPtP->Draw("p");
	drawLatex(0.2,0.75,"Rec. MC e^{+}",22,0.06,1);
	c->SaveAs("positron_momRes.png");
	c->SaveAs("positron_momRes.pdf");
	c->SaveAs("positron_momRes.eps");

	hPResvsPQ->GetXaxis()->SetRange(1,1);
	TH2F *hPResvsPE = (TH2F *)hPResvsPQ->Project3D("ZY");
	hPResvsPE->SetNameTitle("hPResvsPE","hPResvsPE");
	hPResvsPQ->GetXaxis()->SetRange(2,2);
	TH2F *hPResvsPP = (TH2F *)hPResvsPQ->Project3D("ZY");
	hPResvsPP->SetNameTitle("hPResvsPP","hPResvsPP");

	c1->Clear();
	c1->Divide(nColumns,nRaws);
	nPads = nColumns*nRaws;
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hPResvsPE->SetNameTitle("hPResvsPE","hPResvsPE");
	hPResvsPE->RebinX(10);
	hPResvsPE->GetXaxis()->SetRangeUser(0,7.);
	hPResvsPE->GetYaxis()->SetRangeUser(-0.2,0.2);
	hPResvsPE->DrawClone("colz");
	hPResvsPE->FitSlicesY();
	TH1F *hMeanPResvsPE = (TH1F *)gDirectory->Get("hPResvsPE_1");
	setHisto(hMeanPResvsPE,20,0.8,1,1);
	hMeanPResvsPE->Draw("same");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	TH1F *hSigmaPResvsPE = (TH1F *)gDirectory->Get("hPResvsPE_2");
	hSigmaPResvsPE->GetXaxis()->SetRangeUser(0,7.);
	hSigmaPResvsPE->GetYaxis()->SetRangeUser(0,0.025);
	hSigmaPResvsPE->GetYaxis()->SetTitle("Momentum Resolution");
	setHisto(hSigmaPResvsPE,20,0.8,1,1);
	hSigmaPResvsPE->Draw("pe");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hPResvsPP->SetNameTitle("hPResvsPP","hPResvsPP");
	hPResvsPP->RebinX(10);
	hPResvsPP->GetXaxis()->SetRangeUser(0,7.);
	hPResvsPP->GetYaxis()->SetRangeUser(-0.2,0.2);
	hPResvsPP->DrawClone("colz");
	hPResvsPP->FitSlicesY();
	TH1F *hMeanPResvsPP = (TH1F *)gDirectory->Get("hPResvsPP_1");
	setHisto(hMeanPResvsPP,20,0.8,1,1);
	hMeanPResvsPP->Draw("same");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(4);
	TH1F *hSigmaPResvsPP = (TH1F *)gDirectory->Get("hPResvsPP_2");
	hSigmaPResvsPP->GetXaxis()->SetRangeUser(0,7.);
	hSigmaPResvsPP->GetYaxis()->SetRangeUser(0,0.025);
	hSigmaPResvsPP->GetYaxis()->SetTitle("Momentum Resolution");
	setHisto(hSigmaPResvsPP,20,0.8,1,1);
	hSigmaPResvsPP->Draw("pe");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	hRcEtavsPtQ->GetYaxis()->SetRangeUser(0,5.0);
	hRcEtavsPtQ->GetZaxis()->SetRangeUser(-1.5,1.5);
	hRcEtavsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcEtavsPtE = (TH2F *)hRcEtavsPtQ->Project3D("ZY"); 
	hRcEtavsPtE->SetNameTitle("hRcEtavsPtE","hRcEtavsPtE");
	hRcEtavsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcEtavsPtP = (TH2F *)hRcEtavsPtQ->Project3D("ZY");
	hRcEtavsPtP->SetNameTitle("hRcEtavsPtP","hRcEtavsPtP");
	hDataEtavsPtQ->GetYaxis()->SetRangeUser(0.,5.);
	hDataEtavsPtQ->GetZaxis()->SetRangeUser(-1.5,1.5);
	hDataEtavsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataEtavsPtE = (TH2F *)hDataEtavsPtQ->Project3D("ZY"); 
	hDataEtavsPtE->SetNameTitle("hDataEtavsPtE","hDataEtavsPtE");
	hDataEtavsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataEtavsPtP = (TH2F *)hDataEtavsPtQ->Project3D("ZY");
	hDataEtavsPtP->SetNameTitle("hDataEtavsPtP","hDataEtavsPtP");
	hDataPiEtavsPtQ->GetYaxis()->SetRangeUser(0.,5.);
	hDataPiEtavsPtQ->GetZaxis()->SetRangeUser(-1.5,1.5);
	hDataPiEtavsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiEtavsPtM = (TH2F *)hDataPiEtavsPtQ->Project3D("ZY"); 
	hDataPiEtavsPtM->SetNameTitle("hDataPiEtavsPtM","hDataPiEtavsPtM");
	hDataPiEtavsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiEtavsPtP = (TH2F *)hDataPiEtavsPtQ->Project3D("ZY");
	hDataPiEtavsPtP->SetNameTitle("hDataPiEtavsPtP","hDataPiEtavsPtP");
	cout<<"# of embedded eMinus:"<<hRcEtavsPtE->GetEntries()<<endl;;
	cout<<"# of embedded ePlus:"<<hRcEtavsPtP->GetEntries()<<endl;;
	cout<<"# of data eMinus: "<<hDataEtavsPtE->GetEntries()<<endl;
	cout<<"# of data ePlus: "<<hDataEtavsPtP->GetEntries()<<endl;
	//cout<<"# of data piMinus: "<<hDataPiEtavsPtM->GetEntries()<<endl;
	//cout<<"# of data piPlus: "<<hDataPiEtavsPtP->GetEntries()<<endl;

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcEtavsPtE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcEtavsPtP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataEtavsPtE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataEtavsPtP->Draw("colz");
	drawLatex(textX,textY,"Data. e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiEtavsPtM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiEtavsPtP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns=3;
	nRaws=3;
	nPads=nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	TH1F *hRcEtaE[nMomBins];
	TH1F *hRcEtaP[nMomBins];
	TH1F *hDataEtaE[nMomBins];
	TH1F *hDataEtaP[nMomBins];
	Int_t rebY = 5;
	hRcEtavsPtE->RebinY(rebY);
	hRcEtavsPtP->RebinY(rebY);
	hDataEtavsPtE->RebinY(rebY);
	hDataEtavsPtP->RebinY(rebY);
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcEtavsPtE->GetXaxis()->FindBin(mom[i]);
		momBinHi = hRcEtavsPtE->GetXaxis()->FindBin(mom[i+1]);
		hRcEtaE[i] = (TH1F *)hRcEtavsPtE->ProjectionY(Form("hRcEtaE_PtBin%d",i+1),momBinLow,momBinHi);
		hRcEtaP[i] = (TH1F *)hRcEtavsPtP->ProjectionY(Form("hRcEtaP_PtBin%d",i+1),momBinLow,momBinHi);
		hDataEtaE[i] = (TH1F *)hDataEtavsPtE->ProjectionY(Form("hDataEtaE_PtBin%d",i+1),momBinLow,momBinHi);
		hDataEtaP[i] = (TH1F *)hDataEtavsPtP->ProjectionY(Form("hDataEtaP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcEtaE[i]->SetLineColor(1);
		hRcEtaP[i]->SetLineColor(2);
		setHisto(hDataEtaE[i],20,0.8,4,4);
		setHisto(hDataEtaP[i],20,0.8,kViolet,kViolet);
		binLow = hDataEtaE[i]->GetXaxis()->FindBin(-0.7);
		binHi = hDataEtaE[i]->GetXaxis()->FindBin(0.7);
		hDataEtaE[i]->Scale(hRcEtaE[i]->Integral(binLow,binHi)*1./hDataEtaE[i]->Integral(binLow,binHi));
		hDataEtaP[i]->Scale(hRcEtaP[i]->Integral(binLow,binHi)*1./hDataEtaP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		hRcEtaP[i]->SetMinimum(1.);
		hRcEtaE[i]->Draw("histe");
		hRcEtaP[i]->Draw("histesame");
		hDataEtaE[i]->Draw("pesame");
		hDataEtaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.35);
		leg->SetX2NDC(0.65);
		leg->SetY1NDC(0.2);
		leg->SetY2NDC(0.4);
		leg->AddEntry(hRcEtaE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcEtaP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataEtaE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataEtaP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	TH1F *hRcEtaE1[nPiMomBins];
	TH1F *hRcEtaP1[nPiMomBins];
	TH1F *hDataPiEtaM[nPiMomBins];
	TH1F *hDataPiEtaP[nPiMomBins];
	hDataPiEtavsPtM->RebinY(rebY);
	hDataPiEtavsPtP->RebinY(rebY);
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcEtavsPtE->GetXaxis()->FindBin(pimom[i]);
		momBinHi = hRcEtavsPtE->GetXaxis()->FindBin(pimom[i+1]);
		hRcEtaE1[i] = (TH1F *)hRcEtavsPtE->ProjectionY(Form("hRcEtaE1_PtBin%d",i+1),momBinLow,momBinHi);
		hRcEtaP1[i] = (TH1F *)hRcEtavsPtP->ProjectionY(Form("hRcEtaP1_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiEtaM[i] = (TH1F *)hDataPiEtavsPtM->ProjectionY(Form("hDataPiEtaM_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiEtaP[i] = (TH1F *)hDataPiEtavsPtP->ProjectionY(Form("hDataPiEtaP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcEtaE1[i]->SetLineColor(1);
		hRcEtaP1[i]->SetLineColor(2);
		setHisto(hDataPiEtaM[i],22,0.8,4,4);
		setHisto(hDataPiEtaP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiEtaM[i]->GetXaxis()->FindBin(-0.7);
		binHi = hDataPiEtaM[i]->GetXaxis()->FindBin(0.7);
		hDataPiEtaM[i]->Scale(hRcEtaE1[i]->Integral(binLow,binHi)*1./hDataPiEtaM[i]->Integral(binLow,binHi));
		hDataPiEtaP[i]->Scale(hRcEtaP1[i]->Integral(binLow,binHi)*1./hDataPiEtaP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		hRcEtaP1[i]->SetMinimum(1.);
		hRcEtaE1[i]->Draw("histe");
		hRcEtaP1[i]->Draw("histesame");
		hDataPiEtaM[i]->Draw("pesame");
		hDataPiEtaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.35);
		leg->SetX2NDC(0.65);
		leg->SetY1NDC(0.2);
		leg->SetY2NDC(0.4);
		leg->AddEntry(hRcEtaE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcEtaP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiEtaM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiEtaP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		/*if(i%nPads==nPads-1)
			pdfAction(c1,ps);*/
	}
	/*if(nPiMomBins%nPads!=0)
		pdfAction(c1,ps);*/

	hRcPhivsPtQ->GetYaxis()->SetRangeUser(0,5.0);
	hRcPhivsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcPhivsPtE = (TH2F *)hRcPhivsPtQ->Project3D("ZY"); 
	hRcPhivsPtE->SetNameTitle("hRcPhivsPtE","hRcPhivsPtE");
	hRcPhivsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcPhivsPtP = (TH2F *)hRcPhivsPtQ->Project3D("ZY");
	hRcPhivsPtP->SetNameTitle("hRcPhivsPtP","hRcPhivsPtP");
	hDataPhivsPtQ->GetYaxis()->SetRangeUser(0.,5.);
	hDataPhivsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPhivsPtE = (TH2F *)hDataPhivsPtQ->Project3D("ZY"); 
	hDataPhivsPtE->SetNameTitle("hDataPhivsPtE","hDataPhivsPtE");
	hDataPhivsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPhivsPtP = (TH2F *)hDataPhivsPtQ->Project3D("ZY");
	hDataPhivsPtP->SetNameTitle("hDataPhivsPtP","hDataPhivsPtP");
	hDataPiPhivsPtQ->GetYaxis()->SetRangeUser(0.,5.);
	hDataPiPhivsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiPhivsPtM = (TH2F *)hDataPiPhivsPtQ->Project3D("ZY"); 
	hDataPiPhivsPtM->SetNameTitle("hDataPiPhivsPtM","hDataPiPhivsPtM");
	hDataPiPhivsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiPhivsPtP = (TH2F *)hDataPiPhivsPtQ->Project3D("ZY");
	hDataPiPhivsPtP->SetNameTitle("hDataPiPhivsPtP","hDataPiPhivsPtP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcPhivsPtE->Draw("colz");
	funNegHiMirror->SetRange(0.18,5.);
	//funNegHiMirror->Draw("same");
	funNegLowMirror->SetRange(0.18,5.);
	//funNegLowMirror->Draw("same");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcPhivsPtP->Draw("colz");
	funPosHi->SetRange(0.18,5.);
	//funPosHi->Draw("same");
	funPosLow->SetRange(0.18,5.);
	//funPosLow->Draw("same");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataPhivsPtE->Draw("colz");
	//funNegHiMirror->Draw("same");
	//funNegLowMirror->Draw("same");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataPhivsPtP->Draw("colz");
	//funPosHi->Draw("same");
	//funPosLow->Draw("same");
	drawLatex(textX,textY,"Data. e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiPhivsPtM->Draw("colz");
	//funNegHiMirror->Draw("same");
	//funNegLowMirror->Draw("same");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiPhivsPtP->Draw("colz");
	//funPosHi->Draw("same");
	//funPosLow->Draw("same");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	TH1F *hRcPhiE[nMomBins];
	TH1F *hRcPhiP[nMomBins];
	TH1F *hDataPhiE[nMomBins];
	TH1F *hDataPhiP[nMomBins];
	rebY = 15;
	hRcPhivsPtE->RebinY(rebY);
	hRcPhivsPtP->RebinY(rebY);
	hDataPhivsPtE->RebinY(3);
	hDataPhivsPtP->RebinY(3);
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcPhivsPtE->GetXaxis()->FindBin(mom[i]+1.e-3);
		momBinHi = hRcPhivsPtE->GetXaxis()->FindBin(mom[i+1]-1.e-3);
		hRcPhiE[i] = (TH1F *)hRcPhivsPtE->ProjectionY(Form("hRcPhiE_PtBin%d",i+1),momBinLow,momBinHi);
		hRcPhiP[i] = (TH1F *)hRcPhivsPtP->ProjectionY(Form("hRcPhiP_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPhiE[i] = (TH1F *)hDataPhivsPtE->ProjectionY(Form("hDataPhiE_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPhiP[i] = (TH1F *)hDataPhivsPtP->ProjectionY(Form("hDataPhiP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcPhiE[i]->SetLineColor(1);
		hRcPhiP[i]->SetLineColor(2);
		//hDataPhiE[i]->SetLineColor(4);
		//hDataPhiE[i]->SetLineColor(kViolet);
		setHisto(hDataPhiE[i],20,0.8,4,4);
		setHisto(hDataPhiP[i],20,0.8,kViolet,kViolet);
		binLow = hRcPhiE[i]->GetXaxis()->FindBin(1.);
		binHi = hRcPhiE[i]->GetXaxis()->FindBin(3.);
		binLow1 = hDataPhiE[i]->GetXaxis()->FindBin(1.);
		binHi1 = hDataPhiE[i]->GetXaxis()->FindBin(3.);
		//cout<<binLow<<"  "<<binHi<<endl;
		//cout<<binLow1<<"  "<<binHi1<<endl;
		hDataPhiE[i]->Scale(rebY/3.*hRcPhiE[i]->Integral(binLow,binHi)*1./hDataPhiE[i]->Integral(binLow1,binHi1));
		hDataPhiP[i]->Scale(rebY/3.*hRcPhiP[i]->Integral(binLow,binHi)*1./hDataPhiP[i]->Integral(binLow1,binHi1));
		c1->cd(i%nPads+1);
		hDataPhiP[i]->SetMinimum(1.);
		hDataPhiE[i]->Draw("pe");
		hDataPhiE[i]->Draw("histsame");
		hDataPhiP[i]->Draw("pesame");
		hDataPhiP[i]->Draw("histsame");
		hRcPhiE[i]->Draw("histesame");
		hRcPhiP[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.55);
		leg->SetX2NDC(0.75);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcPhiE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcPhiP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPhiE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataPhiP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	TH1F *hRcPhiE1[nPiMomBins];
	TH1F *hRcPhiP1[nPiMomBins];
	TH1F *hDataPiPhiM[nPiMomBins];
	TH1F *hDataPiPhiP[nPiMomBins];
	hDataPiPhivsPtM->RebinY(3);
	hDataPiPhivsPtP->RebinY(3);
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcPhivsPtE->GetXaxis()->FindBin(pimom[i]+1.e-3);
		momBinHi = hRcPhivsPtE->GetXaxis()->FindBin(pimom[i+1]-1.e-3);
		hRcPhiE1[i] = (TH1F *)hRcPhivsPtE->ProjectionY(Form("hRcPhiE1_PtBin%d",i+1),momBinLow,momBinHi);
		hRcPhiP1[i] = (TH1F *)hRcPhivsPtP->ProjectionY(Form("hRcPhiP1_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiPhiM[i] = (TH1F *)hDataPiPhivsPtM->ProjectionY(Form("hDataPiPhiM_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiPhiP[i] = (TH1F *)hDataPiPhivsPtP->ProjectionY(Form("hDataPiPhiP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcPhiE1[i]->SetLineColor(1);
		hRcPhiP1[i]->SetLineColor(2);
		setHisto(hDataPiPhiM[i],22,0.8,4,4);
		setHisto(hDataPiPhiP[i],22,0.8,kViolet,kViolet);
		binLow = hRcPhiE1[i]->GetXaxis()->FindBin(0.);
		binHi = hRcPhiP1[i]->GetXaxis()->FindBin(3.);
		binLow1 = hDataPiPhiM[i]->GetXaxis()->FindBin(0.);
		binHi1 = hDataPiPhiM[i]->GetXaxis()->FindBin(3.);
		hDataPiPhiM[i]->Scale(rebY/3.*hRcPhiE1[i]->Integral(binLow,binHi)*1./hDataPiPhiM[i]->Integral(binLow1,binHi1));
		hDataPiPhiP[i]->Scale(rebY/3.*hRcPhiP1[i]->Integral(binLow,binHi)*1./hDataPiPhiP[i]->Integral(binLow1,binHi1));
		c1->cd(i%nPads+1);
		hDataPiPhiP[i]->SetMinimum(1.);
		Float_t max = hDataPiPhiP[i]->GetMaximum();
		hDataPiPhiP[i]->SetMaximum(max*1.2);
		hDataPiPhiM[i]->Draw("pe");
		hDataPiPhiM[i]->Draw("histsame");
		hDataPiPhiP[i]->Draw("pesame");
		hDataPiPhiP[i]->Draw("histsame");
		hRcPhiE1[i]->Draw("histesame");
		hRcPhiP1[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.55);
		leg->SetX2NDC(0.75);
		leg->SetY1NDC(0.20);
		leg->SetY2NDC(0.40);
		leg->AddEntry(hRcPhiE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcPhiP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiPhiM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiPhiP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		//if(i%nPads==nPads-1)
			//pdfAction(c1,ps);
	}
	//if(nPiMomBins%nPads!=0)
		//pdfAction(c1,ps);

	hRcNHitsCommonvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsCommonvsPtE = (TH2F *)hRcNHitsCommonvsPtQ->Project3D("ZY");
	hRcNHitsCommonvsPtE->SetNameTitle("hRcNHitsCommonvsPtE","hRcNHitsCommonvsPtE");
	hRcNHitsCommonvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsCommonvsPtP = (TH2F *)hRcNHitsCommonvsPtQ->Project3D("ZY");
	hRcNHitsCommonvsPtP->SetNameTitle("hRcNHitsCommonvsPtP","hRcNHitsCommonvsPtP");
	hRcNHitsCommonvsEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsCommonvsEtaE = (TH2F *)hRcNHitsCommonvsEtaQ->Project3D("ZY");
	hRcNHitsCommonvsEtaE->SetNameTitle("hRcNHitsCommonvsEtaE","hRcNHitsCommonvsEtaE");
	hRcNHitsCommonvsEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsCommonvsEtaP = (TH2F *)hRcNHitsCommonvsEtaQ->Project3D("ZY");
	hRcNHitsCommonvsEtaP->SetNameTitle("hRcNHitsCommonvsEtaP","hRcNHitsCommonvsEtaP");
	hRcNHitsCommonvsPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsCommonvsPhiE = (TH2F *)hRcNHitsCommonvsPhiQ->Project3D("ZY");
	hRcNHitsCommonvsPhiE->SetNameTitle("hRcNHitsCommonvsPhiE","hRcNHitsCommonvsPhiE");
	hRcNHitsCommonvsPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsCommonvsPhiP = (TH2F *)hRcNHitsCommonvsPhiQ->Project3D("ZY");
	hRcNHitsCommonvsPhiP->SetNameTitle("hRcNHitsCommonvsPhiP","hRcNHitsCommonvsPhiP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	c1->Divide(nColumns,nRaws);
	nPads = nColumns*nRaws;
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsCommonvsPtE->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsCommonvsPtE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsCommonvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsCommonvsPtP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hRcNHitsCommonvsEtaE->GetXaxis()->SetRangeUser(-1.25,1.25);
	hRcNHitsCommonvsEtaE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hRcNHitsCommonvsEtaP->GetXaxis()->SetRangeUser(-1.25,1.25);
	hRcNHitsCommonvsEtaP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsCommonvsPhiE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsCommonvsPhiP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	TH1F *hRcNHitsCommonPtE[nMomBins];
	TH1F *hRcNHitsCommonPtP[nMomBins];
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcNHitsCommonvsPtE->GetXaxis()->FindBin(mom[i]+1.e-3);
		momBinHi = hRcNHitsCommonvsPtE->GetXaxis()->FindBin(mom[i+1]-1.e-3);
		hRcNHitsCommonPtE[i] = (TH1F *)hRcNHitsCommonvsPtE->ProjectionY(Form("hRcNHitsCommonPtE_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsCommonPtP[i] = (TH1F *)hRcNHitsCommonvsPtP->ProjectionY(Form("hRcNHitsCommonPtP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsCommonPtE[i]->SetLineColor(1);
		hRcNHitsCommonPtP[i]->SetLineColor(2);
		c1->cd(i%nPads+1);
		hRcNHitsCommonPtE[i]->Draw("histe");
		hRcNHitsCommonPtP[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsCommonPtE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsCommonPtP[i],"Rec. MC e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.45,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.04,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	const Int_t rebEta = 10;
	clearPad(c1,nPads);
	hRcNHitsCommonvsEtaE->RebinX(rebEta);
	hRcNHitsCommonvsEtaP->RebinX(rebEta);
	const Int_t nEtaBins = hRcNHitsCommonvsEtaE->GetNbinsX();
	TH1F *hRcNHitsCommonEtaE[nEtaBins];
	TH1F *hRcNHitsCommonEtaP[nEtaBins];
	Int_t nMinEtaBin = hRcNHitsCommonvsEtaE->GetXaxis()->FindBin(-1.+1.e-3);
	Int_t nMaxEtaBin = hRcNHitsCommonvsEtaE->GetXaxis()->FindBin(1-1.e-3);
	for(Int_t i=nMinEtaBin;i<=nMaxEtaBin;i++){
		hRcNHitsCommonEtaE[i] = (TH1F *)hRcNHitsCommonvsEtaE->ProjectionY(Form("hRcNHitsCommonEtaE_EtaBin%d",i+1-nMinEtaBin),i,i);
		hRcNHitsCommonEtaP[i] = (TH1F *)hRcNHitsCommonvsEtaP->ProjectionY(Form("hRcNHitsCommonEtaP_EtaBin%d",i+1-nMinEtaBin),i,i);
		hRcNHitsCommonEtaE[i]->SetLineColor(1);
		hRcNHitsCommonEtaP[i]->SetLineColor(2);
		c1->cd((i-nMinEtaBin)%nPads+1);
		hRcNHitsCommonEtaP[i]->SetMinimum(1.);;
		hRcNHitsCommonEtaE[i]->Draw("histe");
		hRcNHitsCommonEtaP[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsCommonEtaE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsCommonEtaP[i],"Rec. MC e^{+}","pl");
		leg->DrawClone("same");
		Float_t EtaLow = hRcNHitsCommonvsEtaE->GetXaxis()->GetBinLowEdge(i);
		Float_t EtaHi = hRcNHitsCommonvsEtaE->GetXaxis()->GetBinLowEdge(i+1);
		drawLatex(0.40,0.95,Form("%3.1f<#eta<%3.1f",EtaLow,EtaHi),42,0.06,1);
		if((i-nMinEtaBin)%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if((nMaxEtaBin-nMinEtaBin+1)%nPads!=0)
		pdfAction(c1,ps);

	clearPad(c1,nPads);
	const Int_t rebPhi = 30;
	hRcNHitsCommonvsPhiE->RebinX(rebPhi);
	hRcNHitsCommonvsPhiP->RebinX(rebPhi);
	const Int_t nPhiBins = hRcNHitsCommonvsPhiE->GetNbinsX();
	TH1F *hRcNHitsCommonPhiE[nPhiBins];
	TH1F *hRcNHitsCommonPhiP[nPhiBins];
	for(Int_t i=0;i<nPhiBins;i++){
		hRcNHitsCommonPhiE[i] = (TH1F *)hRcNHitsCommonvsPhiE->ProjectionY(Form("hRcNHitsCommonPhiE_PhiBin%d",i+1),i+1,i+1);
		hRcNHitsCommonPhiP[i] = (TH1F *)hRcNHitsCommonvsPhiP->ProjectionY(Form("hRcNHitsCommonPhiP_PhiBin%d",i+1),i+1,i+1);
		hRcNHitsCommonPhiE[i]->SetLineColor(1);
		hRcNHitsCommonPhiP[i]->SetLineColor(2);
		c1->cd(i%nPads+1);
		hRcNHitsCommonPhiP[i]->SetMinimum(1.);;
		hRcNHitsCommonPhiE[i]->Draw("histe");
		hRcNHitsCommonPhiP[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsCommonPhiE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsCommonPhiP[i],"Rec. MC e^{+}","pl");
		leg->DrawClone("same");
		Float_t PhiLow = hRcNHitsCommonvsPhiE->GetXaxis()->GetBinLowEdge(i+1);
		Float_t PhiHi = hRcNHitsCommonvsPhiE->GetXaxis()->GetBinLowEdge(i+2);
		drawLatex(0.40,0.95,Form("%3.2f<#phi<%3.2f",PhiLow,PhiHi),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nPhiBins%nPads!=0)
		pdfAction(c1,ps);

	hRcNHitsFitvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsFitvsPtE = (TH2F *)hRcNHitsFitvsPtQ->Project3D("ZY");
	hRcNHitsFitvsPtE->SetNameTitle("hRcNHitsFitvsPtE","hRcNHitsFitvsPtE");
	hRcNHitsFitvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsFitvsPtP = (TH2F *)hRcNHitsFitvsPtQ->Project3D("ZY");
	hRcNHitsFitvsPtP->SetNameTitle("hRcNHitsFitvsPtP","hRcNHitsFitvsPtP");
	hDataNHitsFitvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNHitsFitvsPtE = (TH2F *)hDataNHitsFitvsPtQ->Project3D("ZY");
	hDataNHitsFitvsPtE->SetNameTitle("hDataNHitsFitvsPtE","hDataNHitsFitvsPtE");
	hDataNHitsFitvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNHitsFitvsPtP = (TH2F *)hDataNHitsFitvsPtQ->Project3D("ZY");
	hDataNHitsFitvsPtP->SetNameTitle("hDataNHitsFitvsPtP","hDataNHitsFitvsPtP");
	hDataPiNHitsFitvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNHitsFitvsPtM = (TH2F *)hDataPiNHitsFitvsPtQ->Project3D("ZY");
	hDataPiNHitsFitvsPtM->SetNameTitle("hDataPiNHitsFitvsPtM","hDataPiNHitsFitvsPtM");
	hDataPiNHitsFitvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNHitsFitvsPtP = (TH2F *)hDataPiNHitsFitvsPtQ->Project3D("ZY");
	hDataPiNHitsFitvsPtP->SetNameTitle("hDataPiNHitsFitvsPtP","hDataPiNHitsFitvsPtP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsFitvsPtE->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsFitvsPtE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsFitvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsFitvsPtP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNHitsFitvsPtE->GetXaxis()->SetRangeUser(0,5.);
	hDataNHitsFitvsPtE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNHitsFitvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataNHitsFitvsPtP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNHitsFitvsPtM->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNHitsFitvsPtM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNHitsFitvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNHitsFitvsPtP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	TH1F *hRcNHitsFitE[nMomBins];
	TH1F *hRcNHitsFitP[nMomBins];
	TH1F *hDataNHitsFitE[nMomBins];
	TH1F *hDataNHitsFitP[nMomBins];
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcNHitsFitvsPtE->GetXaxis()->FindBin(mom[i]);
		momBinHi = hRcNHitsFitvsPtE->GetXaxis()->FindBin(mom[i+1]);
		hRcNHitsFitE[i] = (TH1F *)hRcNHitsFitvsPtE->ProjectionY(Form("hRcNHitsFitE_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsFitP[i] = (TH1F *)hRcNHitsFitvsPtP->ProjectionY(Form("hRcNHitsFitP_PtBin%d",i+1),momBinLow,momBinHi);
		hDataNHitsFitE[i] = (TH1F *)hDataNHitsFitvsPtE->ProjectionY(Form("hDataNHitsFitE_PtBin%d",i+1),momBinLow,momBinHi);
		hDataNHitsFitP[i] = (TH1F *)hDataNHitsFitvsPtP->ProjectionY(Form("hDataNHitsFitP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsFitE[i]->SetLineColor(1);
		hRcNHitsFitP[i]->SetLineColor(2);
		setHisto(hDataNHitsFitE[i],20,0.8,4,4);
		setHisto(hDataNHitsFitP[i],20,0.8,kViolet,kViolet);
		binLow = hDataNHitsFitE[i]->GetXaxis()->FindBin(30);
		binHi = hDataNHitsFitE[i]->GetXaxis()->FindBin(40);
		// hDataNHitsFitE[i]->Scale(hRcNHitsFitE[i]->Integral(binLow,binHi)*1./hDataNHitsFitE[i]->Integral(binLow,binHi));
		// hDataNHitsFitP[i]->Scale(hRcNHitsFitP[i]->Integral(binLow,binHi)*1./hDataNHitsFitP[i]->Integral(binLow,binHi));
		hDataNHitsFitE[i]->Scale(1./hDataNHitsFitE[i]->Integral());
		hDataNHitsFitP[i]->Scale(1./hDataNHitsFitP[i]->Integral());
		hRcNHitsFitE[i]->Scale(1./hRcNHitsFitE[i]->Integral());
		hRcNHitsFitP[i]->Scale(1./hRcNHitsFitP[i]->Integral());
		c1->cd(i%nPads+1);
		hRcNHitsFitP[i]->SetMinimum(1.);;
		hRcNHitsFitE[i]->Draw("histe");
		hRcNHitsFitP[i]->Draw("histesame");
		hDataNHitsFitE[i]->Draw("pesame");
		hDataNHitsFitP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsFitE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsFitP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNHitsFitE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNHitsFitP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	clearPad(c1,nPads);
	TH1F *hRcNHitsFitE1[nPiMomBins];
	TH1F *hRcNHitsFitP1[nPiMomBins];
	TH1F *hDataPiNHitsFitM[nPiMomBins];
	TH1F *hDataPiNHitsFitP[nPiMomBins];
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcNHitsFitvsPtE->GetXaxis()->FindBin(pimom[i]);
		momBinHi = hRcNHitsFitvsPtE->GetXaxis()->FindBin(pimom[i+1]);
		hRcNHitsFitE1[i] = (TH1F *)hRcNHitsFitvsPtE->ProjectionY(Form("hRcNHitsFitE1_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsFitP1[i] = (TH1F *)hRcNHitsFitvsPtP->ProjectionY(Form("hRcNHitsFitP1_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiNHitsFitM[i] = (TH1F *)hDataPiNHitsFitvsPtM->ProjectionY(Form("hDataPiNHitsFitM_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiNHitsFitP[i] = (TH1F *)hDataPiNHitsFitvsPtP->ProjectionY(Form("hDataPiNHitsFitP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsFitE1[i]->SetLineColor(1);
		hRcNHitsFitP1[i]->SetLineColor(2);
		setHisto(hDataPiNHitsFitM[i],22,0.8,4,4);
		setHisto(hDataPiNHitsFitP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiNHitsFitM[i]->GetXaxis()->FindBin(30);
		binHi = hDataPiNHitsFitM[i]->GetXaxis()->FindBin(40);
		hDataPiNHitsFitM[i]->Scale(hRcNHitsFitE1[i]->Integral(binLow,binHi)*1./hDataPiNHitsFitM[i]->Integral(binLow,binHi));
		hDataPiNHitsFitP[i]->Scale(hRcNHitsFitP1[i]->Integral(binLow,binHi)*1./hDataPiNHitsFitP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		hRcNHitsFitP1[i]->SetMinimum(1.);
		hRcNHitsFitE1[i]->Draw("histe");
		hRcNHitsFitP1[i]->Draw("histesame");
		hDataPiNHitsFitM[i]->Draw("pesame");
		hDataPiNHitsFitP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsFitE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsFitP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNHitsFitM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNHitsFitP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		/*if(i%nPads==nPads-1)
			pdfAction(c1,ps);*/
	}
	/*if(nPiMomBins%nPads!=0)
		pdfAction(c1,ps);*/
	
	cout<<"before fit pt eta"<<endl;
	clearPad(c1,nPads);
	TH1D *hDataNHitsFitPtEta[nMomBins * NumEtaBins];
    TH1D *hMCNHitsFitPtEta[nMomBins * NumEtaBins];
	for (int i = 0; i < nMomBins; i++)
    {
        for (int j = 0; j < NumEtaBins; j++)
        {

            int PtBinLow = hRcNHitsFitPtEtaMinus->GetXaxis()->FindBin(mom[i] + 1.e-3);
            int PtBinHigh = hRcNHitsFitPtEtaMinus->GetXaxis()->FindBin(mom[i + 1] - 1.e-3);
            int EtaBinLow = hRcNHitsFitPtEtaMinus->GetYaxis()->FindBin(eta[j] + 1.e-3);
            int EtaBinHigh = hRcNHitsFitPtEtaMinus->GetYaxis()->FindBin(eta[j + 1] - 1.e-3);
            hDataNHitsFitPtEta[i * NumEtaBins + j] = hDataNHitsFitPtEtaMinus->ProjectionZ(Form("p_{T}: %1.2f~%1.2f #eta: %1.2f~%1.2f data",mom[i],mom[i+1],eta[j],eta[j+1]),PtBinLow, PtBinHigh, EtaBinLow, EtaBinHigh);
            hMCNHitsFitPtEta[i * NumEtaBins + j] = hRcNHitsFitPtEtaMinus->ProjectionZ(Form("p_{T}: %1.2f~%1.2f #eta: %1.2f~%1.2f MC",mom[i],mom[i+1],eta[j],eta[j+1]),PtBinLow, PtBinHigh, EtaBinLow, EtaBinHigh);
            setHisto(hDataNHitsFitPtEta[i * NumEtaBins + j],20,0.8,1,1);
            setHisto(hMCNHitsFitPtEta[i * NumEtaBins + j],20,0.8,4,4);
            hMCNHitsFitPtEta[i * NumEtaBins + j]->SetTitle(Form("p_{T}: %1.2f~%1.2f #eta: %1.2f~%1.2f",mom[i],mom[i+1],eta[j],eta[j+1]));
            c1->cd((i * NumEtaBins + j)%nPads+1);
            int binLow = hDataNHitsFitPtEta[i * NumEtaBins + j]->GetXaxis()->FindBin(30);
		    int binHi = hDataNHitsFitPtEta[i * NumEtaBins + j]->GetXaxis()->FindBin(40);
            hDataNHitsFitPtEta[i * NumEtaBins + j]->Scale(hMCNHitsFitPtEta[i * NumEtaBins + j]->Integral(binLow,binHi)*1.0/hDataNHitsFitPtEta[i * NumEtaBins + j]->Integral(binLow,binHi));
            hMCNHitsFitPtEta[i * NumEtaBins + j]->Draw("pe");
            hDataNHitsFitPtEta[i * NumEtaBins + j]->Draw("pesame");

            leg->Clear();
            leg->SetX1NDC(0.2);
            leg->SetX2NDC(0.45);
            leg->SetY1NDC(0.65);
            leg->SetY2NDC(0.85);
            leg->AddEntry(hDataNHitsFitPtEta[i * NumEtaBins + j], "Data e^{-}", "pl");
            leg->AddEntry(hMCNHitsFitPtEta[i * NumEtaBins + j], "Rec. MC e^{-}", "pl");
            leg->Draw("same");
            /*if((i * NumEtaBins + j)%nPads==nPads-1) 
            {
                //TString name = 
                pdfAction(c1,ps);
                clearPad(c1,nPads);
            }*/
        }
    }
	if(nMomBins * NumEtaBins %nPads!=0)
		pdfAction(c1,ps);

	cout<<"before fit pT phi"<<endl;
	clearPad(c1,nPads);
	TH1D *hDataNHitsFitPtPhi[nMomBins * NumPhiBins];
    TH1D *hMCNHitsFitPtPhi[nMomBins * NumPhiBins];
	for (int i = 0; i < nMomBins; i++)
    {
        for (int j = 0; j < NumPhiBins; j++)
        {
			int PtBinLow = hRcNHitsFitPtPhiMinus->GetXaxis()->FindBin(mom[i] + 1.e-3);
            int PtBinHigh = hRcNHitsFitPtPhiMinus->GetXaxis()->FindBin(mom[i + 1] - 1.e-3);
            int PhiBinLow = hRcNHitsFitPtPhiMinus->GetYaxis()->FindBin(-TMath::Pi()+j*PhiRange + 1.e-3);
            int PhiBinHigh = hRcNHitsFitPtPhiMinus->GetYaxis()->FindBin(-TMath::Pi()+(j+1)*PhiRange - 1.e-3);
            hDataNHitsFitPtPhi[i * NumPhiBins + j] = hDataNHitsFitPtPhiMinus->ProjectionZ(Form("p_{T}: %1.2f~%1.2f #Phi: %1.2f~%1.2f data",mom[i],mom[i+1],-TMath::Pi()+j*PhiRange + 1.e-3,-TMath::Pi()+(j+1)*PhiRange - 1.e-3),PtBinLow, PtBinHigh, PhiBinLow, PhiBinHigh);
            hMCNHitsFitPtPhi[i * NumPhiBins + j] = hRcNHitsFitPtPhiMinus->ProjectionZ(Form("p_{T}: %1.2f~%1.2f #Phi: %1.2f~%1.2f MC",mom[i],mom[i+1],-TMath::Pi()+j*PhiRange + 1.e-3,-TMath::Pi()+(j+1)*PhiRange - 1.e-3),PtBinLow, PtBinHigh, PhiBinLow, PhiBinHigh);
            setHisto(hDataNHitsFitPtPhi[i * NumPhiBins + j],20,0.8,1,1);
            setHisto(hMCNHitsFitPtPhi[i * NumPhiBins + j],20,0.8,4,4);
			hMCNHitsFitPtPhi[i * NumPhiBins + j]->SetTitle(Form("p_{T}: %1.2f~%1.2f #Phi: %1.2f~%1.2f",mom[i],mom[i+1],-TMath::Pi()+j*PhiRange,-TMath::Pi()+(j+1)*PhiRange));
            hMCNHitsFitPtPhi[i * NumPhiBins + j]->SetTitleSize(0.015);
            c1->cd((i * NumPhiBins + j)%nPads+1);
            int binLow = hDataNHitsFitPtPhi[i * NumPhiBins + j]->GetXaxis()->FindBin(30);
		    int binHi = hDataNHitsFitPtPhi[i * NumPhiBins + j]->GetXaxis()->FindBin(40);
            hDataNHitsFitPtPhi[i * NumPhiBins + j]->Scale(hMCNHitsFitPtPhi[i * NumPhiBins + j]->Integral(binLow,binHi)*1.0/hDataNHitsFitPtPhi[i * NumPhiBins + j]->Integral(binLow,binHi));
			hMCNHitsFitPtPhi[i * NumPhiBins + j]->Draw("pe");
            hDataNHitsFitPtPhi[i * NumPhiBins + j]->Draw("pesame");

			leg->Clear();
            leg->SetX1NDC(0.2);
            leg->SetX2NDC(0.45);
            leg->SetY1NDC(0.65);
            leg->SetY2NDC(0.85);
            leg->AddEntry(hDataNHitsFitPtPhi[i * NumPhiBins + j], "Data e^{-}", "pl");
            leg->AddEntry(hMCNHitsFitPtPhi[i * NumPhiBins + j], "Rec. MC e^{-}", "pl");
            leg->Draw("same");
            /*if((i * NumPhiBins + j)%nPads==nPads-1) 
            {
                //TString name = 
                pdfAction(c1,ps);
                //c1->SaveAs(Form("pic_phi%d.png",k));
                //k++;
                clearPad(c1,nPads);
            }*/

		}
	}
	if(nMomBins * NumPhiBins %nPads!=0)
		pdfAction(c1,ps);

	hRcNHitsFitvsEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsFitvsEtaE = (TH2F *)hRcNHitsFitvsEtaQ->Project3D("ZY");
	hRcNHitsFitvsEtaE->SetNameTitle("hRcNHitsFitvsEtaE","hRcNHitsFitvsEtaE");
	hRcNHitsFitvsEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsFitvsEtaP = (TH2F *)hRcNHitsFitvsEtaQ->Project3D("ZY");
	hRcNHitsFitvsEtaP->SetNameTitle("hRcNHitsFitvsEtaP","hRcNHitsFitvsEtaP");
	hDataNHitsFitvsEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNHitsFitvsEtaE = (TH2F *)hDataNHitsFitvsEtaQ->Project3D("ZY");
	hDataNHitsFitvsEtaE->SetNameTitle("hDataNHitsFitvsEtaE","hDataNHitsFitvsEtaE");
	hDataNHitsFitvsEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNHitsFitvsEtaP = (TH2F *)hDataNHitsFitvsEtaQ->Project3D("ZY");
	hDataNHitsFitvsEtaP->SetNameTitle("hDataNHitsFitvsEtaP","hDataNHitsFitvsEtaP");
	hDataPiNHitsFitvsEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNHitsFitvsEtaM = (TH2F *)hDataPiNHitsFitvsEtaQ->Project3D("ZY");
	hDataPiNHitsFitvsEtaM->SetNameTitle("hDataPiNHitsFitvsEtaM","hDataPiNHitsFitvsEtaM");
	hDataPiNHitsFitvsEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNHitsFitvsEtaP = (TH2F *)hDataPiNHitsFitvsEtaQ->Project3D("ZY");
	hDataPiNHitsFitvsEtaP->SetNameTitle("hDataPiNHitsFitvsEtaP","hDataPiNHitsFitvsEtaP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsFitvsEtaE->GetXaxis()->SetRangeUser(-1.25,1.25);
	hRcNHitsFitvsEtaE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsFitvsEtaP->GetXaxis()->SetRangeUser(-1.25,1.25);
	hRcNHitsFitvsEtaP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNHitsFitvsEtaE->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataNHitsFitvsEtaE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNHitsFitvsEtaP->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataNHitsFitvsEtaP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNHitsFitvsEtaM->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataPiNHitsFitvsEtaM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNHitsFitvsEtaP->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataPiNHitsFitvsEtaP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);

	clearPad(c1,nPads);
	hRcNHitsFitvsEtaE->RebinX(rebEta);
	hRcNHitsFitvsEtaP->RebinX(rebEta);
	hDataNHitsFitvsEtaE->RebinX(rebEta);
	hDataNHitsFitvsEtaP->RebinX(rebEta);
	TH1F *hRcNHitsFitEtaE[nEtaBins];
	TH1F *hRcNHitsFitEtaP[nEtaBins];
	TH1F *hDataNHitsFitEtaE[nEtaBins];
	TH1F *hDataNHitsFitEtaP[nEtaBins];
	for(Int_t i=nMinEtaBin;i<=nMaxEtaBin;i++){
		hRcNHitsFitEtaE[i] = (TH1F *)hRcNHitsFitvsEtaE->ProjectionY(Form("hRcNHitsFitE_EtaBin%d",i-nMinEtaBin+1),i,i);
		hRcNHitsFitEtaP[i] = (TH1F *)hRcNHitsFitvsEtaP->ProjectionY(Form("hRcNHitsFitP_EtaBin%d",i-nMinEtaBin+1),i,i);
		hDataNHitsFitEtaE[i] = (TH1F *)hDataNHitsFitvsEtaE->ProjectionY(Form("hDataNHitsFitEtaE_EtaBin%d",i-nMinEtaBin+1),i,i);
		hDataNHitsFitEtaP[i] = (TH1F *)hDataNHitsFitvsEtaP->ProjectionY(Form("hDataNHitsFitEtaP_EtaBin%d",i-nMinEtaBin+1),i,i);
		hRcNHitsFitEtaE[i]->SetLineColor(1);
		hRcNHitsFitEtaP[i]->SetLineColor(2);
		setHisto(hDataNHitsFitEtaE[i],20,0.8,4,4);
		setHisto(hDataNHitsFitEtaP[i],20,0.8,kViolet,kViolet);
		binLow = hDataNHitsFitEtaE[i]->GetXaxis()->FindBin(30);
		binHi = hDataNHitsFitEtaE[i]->GetXaxis()->FindBin(40);
		hDataNHitsFitEtaE[i]->Scale(hRcNHitsFitEtaE[i]->Integral(binLow,binHi)*1./hDataNHitsFitEtaE[i]->Integral(binLow,binHi));
		hDataNHitsFitEtaP[i]->Scale(hRcNHitsFitEtaP[i]->Integral(binLow,binHi)*1./hDataNHitsFitEtaP[i]->Integral(binLow,binHi));
		c1->cd((i-nMinEtaBin)%nPads+1);
		hRcNHitsFitEtaP[i]->SetMinimum(1.);;
		hRcNHitsFitEtaE[i]->Draw("histe");
		hRcNHitsFitEtaP[i]->Draw("histesame");
		hDataNHitsFitEtaE[i]->Draw("pesame");
		hDataNHitsFitEtaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsFitEtaE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsFitEtaP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNHitsFitEtaE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNHitsFitEtaP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		Float_t EtaLow = hRcNHitsFitvsEtaE->GetXaxis()->GetBinLowEdge(i);
		Float_t EtaHi = hRcNHitsFitvsEtaE->GetXaxis()->GetBinLowEdge(i+1);
		drawLatex(0.40,0.95,Form("%3.1f<#eta<%3.1f",EtaLow,EtaHi),42,0.06,1);
		//if((i-nMinEtaBin)%nPads==nPads-1)
		//	pdfAction(c1,ps);
	}
	//if((nMaxEtaBin-nMinEtaBin+1)%nPads!=0)
		//pdfAction(c1,ps);

	clearPad(c1,nPads);
	hDataPiNHitsFitvsEtaM->RebinX(rebEta);
	hDataPiNHitsFitvsEtaP->RebinX(rebEta);
	TH1F *hDataPiNHitsFitEtaM[nEtaBins];
	TH1F *hDataPiNHitsFitEtaP[nEtaBins];
	for(Int_t i=nMinEtaBin;i<=nMaxEtaBin;i++){
		hDataPiNHitsFitEtaM[i] = (TH1F *)hDataPiNHitsFitvsEtaM->ProjectionY(Form("hDataPiNHitsFitEtaM_EtaBin%d",i-nMinEtaBin+1),i,i);
		hDataPiNHitsFitEtaP[i] = (TH1F *)hDataPiNHitsFitvsEtaP->ProjectionY(Form("hDataPiNHitsFitEtaP_EtaBin%d",i-nMinEtaBin+1),i,i);
		setHisto(hDataPiNHitsFitEtaM[i],22,0.8,4,4);
		setHisto(hDataPiNHitsFitEtaP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiNHitsFitEtaM[i]->GetXaxis()->FindBin(30);
		binHi = hDataPiNHitsFitEtaM[i]->GetXaxis()->FindBin(40);
		hDataPiNHitsFitEtaM[i]->Scale(hRcNHitsFitEtaE[i]->Integral(binLow,binHi)*1./hDataPiNHitsFitEtaM[i]->Integral(binLow,binHi));
		hDataPiNHitsFitEtaP[i]->Scale(hRcNHitsFitEtaP[i]->Integral(binLow,binHi)*1./hDataPiNHitsFitEtaP[i]->Integral(binLow,binHi));
		c1->cd((i-nMinEtaBin)%nPads+1);
		hRcNHitsFitEtaE[i]->Draw("histe");
		hRcNHitsFitEtaP[i]->Draw("histesame");
		hDataPiNHitsFitEtaM[i]->Draw("pesame");
		hDataPiNHitsFitEtaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsFitEtaE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsFitEtaP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNHitsFitEtaM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNHitsFitEtaP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		Float_t EtaLow = hRcNHitsFitvsEtaE->GetXaxis()->GetBinLowEdge(i);
		Float_t EtaHi = hRcNHitsFitvsEtaE->GetXaxis()->GetBinLowEdge(i+1);
		drawLatex(0.40,0.95,Form("%3.1f<#eta<%3.1f",EtaLow,EtaHi),42,0.06,1);
		//if((i-nMinEtaBin)%nPads==nPads-1)
			//pdfAction(c1,ps);
	}
	//if((nMaxEtaBin-nMinEtaBin+1)%nPads!=0)
	//	pdfAction(c1,ps);

	hRcNHitsFitvsPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsFitvsPhiE = (TH2F *)hRcNHitsFitvsPhiQ->Project3D("ZY");
	hRcNHitsFitvsPhiE->SetNameTitle("hRcNHitsFitvsPhiE","hRcNHitsFitvsPhiE");
	hRcNHitsFitvsPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsFitvsPhiP = (TH2F *)hRcNHitsFitvsPhiQ->Project3D("ZY");
	hRcNHitsFitvsPhiP->SetNameTitle("hRcNHitsFitvsPhiP","hRcNHitsFitvsPhiP");
	hDataNHitsFitvsPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNHitsFitvsPhiE = (TH2F *)hDataNHitsFitvsPhiQ->Project3D("ZY");
	hDataNHitsFitvsPhiE->SetNameTitle("hDataNHitsFitvsPhiE","hDataNHitsFitvsPhiE");
	hDataNHitsFitvsPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNHitsFitvsPhiP = (TH2F *)hDataNHitsFitvsPhiQ->Project3D("ZY");
	hDataNHitsFitvsPhiP->SetNameTitle("hDataNHitsFitvsPhiP","hDataNHitsFitvsPhiP");
	hDataPiNHitsFitvsPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNHitsFitvsPhiM = (TH2F *)hDataPiNHitsFitvsPhiQ->Project3D("ZY");
	hDataPiNHitsFitvsPhiM->SetNameTitle("hDataPiNHitsFitvsPhiM","hDataPiNHitsFitvsPhiM");
	hDataPiNHitsFitvsPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNHitsFitvsPhiP = (TH2F *)hDataPiNHitsFitvsPhiQ->Project3D("ZY");
	hDataPiNHitsFitvsPhiP->SetNameTitle("hDataPiNHitsFitvsPhiP","hDataPiNHitsFitvsPhiP");


	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsFitvsPhiE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsFitvsPhiP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNHitsFitvsPhiE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNHitsFitvsPhiP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNHitsFitvsPhiM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNHitsFitvsPhiP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);

	clearPad(c1,nPads);
	hRcNHitsFitvsPhiE->RebinX(rebPhi);
	hRcNHitsFitvsPhiP->RebinX(rebPhi);
	hDataNHitsFitvsPhiE->RebinX(rebPhi);
	hDataNHitsFitvsPhiP->RebinX(rebPhi);
	TH1F *hRcNHitsFitPhiE[nPhiBins];
	TH1F *hRcNHitsFitPhiP[nPhiBins];
	TH1F *hDataNHitsFitPhiE[nPhiBins];
	TH1F *hDataNHitsFitPhiP[nPhiBins];
	for(Int_t i=0;i<nPhiBins;i++){
		hRcNHitsFitPhiE[i] = (TH1F *)hRcNHitsFitvsPhiE->ProjectionY(Form("hRcNHitsFitE_PhiBin%d",i+1),i+1,i+1);
		hRcNHitsFitPhiP[i] = (TH1F *)hRcNHitsFitvsPhiP->ProjectionY(Form("hRcNHitsFitP_PhiBin%d",i+1),i+1,i+1);
		hDataNHitsFitPhiE[i] = (TH1F *)hDataNHitsFitvsPhiE->ProjectionY(Form("hDataNHitsFitPhiE_PhiBin%d",i+1),i+1,i+1);
		hDataNHitsFitPhiP[i] = (TH1F *)hDataNHitsFitvsPhiP->ProjectionY(Form("hDataNHitsFitPhiP_PhiBin%d",i+1),i+1,i+1);
		hRcNHitsFitPhiE[i]->SetLineColor(1);
		hRcNHitsFitPhiP[i]->SetLineColor(2);
		setHisto(hDataNHitsFitPhiE[i],20,0.8,4,4);
		setHisto(hDataNHitsFitPhiP[i],20,0.8,kViolet,kViolet);
		binLow = hDataNHitsFitPhiE[i]->GetXaxis()->FindBin(30);
		binHi = hDataNHitsFitPhiE[i]->GetXaxis()->FindBin(40);
		// hDataNHitsFitPhiE[i]->Scale(hRcNHitsFitPhiE[i]->Integral(binLow,binHi)*1./hDataNHitsFitPhiE[i]->Integral(binLow,binHi));
		// hDataNHitsFitPhiP[i]->Scale(hRcNHitsFitPhiP[i]->Integral(binLow,binHi)*1./hDataNHitsFitPhiP[i]->Integral(binLow,binHi));
		hDataNHitsFitPhiP[i]->Scale(1./hDataNHitsFitPhiP[i]->Integral());
		hDataNHitsFitPhiE[i]->Scale(1./hDataNHitsFitPhiE[i]->Integral());
		hRcNHitsFitPhiE[i]->Scale(1./hRcNHitsFitPhiE[i]->Integral());
		hRcNHitsFitPhiP[i]->Scale(1./hRcNHitsFitPhiP[i]->Integral());
		c1->cd(i%nPads+1);
		hRcNHitsFitPhiP[i]->SetMinimum(1.);;
		hRcNHitsFitPhiE[i]->Draw("histe");
		hRcNHitsFitPhiP[i]->Draw("histesame");
		hDataNHitsFitPhiE[i]->Draw("pesame");
		hDataNHitsFitPhiP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsFitPhiE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsFitPhiP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNHitsFitPhiE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNHitsFitPhiP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		Float_t PhiLow = hRcNHitsFitvsPhiE->GetXaxis()->GetBinLowEdge(i+1);
		Float_t PhiHi = hRcNHitsFitvsPhiE->GetXaxis()->GetBinLowEdge(i+2);
		drawLatex(0.40,0.95,Form("%3.2f<#phi<%3.2f",PhiLow,PhiHi),42,0.06,1);
		//if(i%nPads==nPads-1)
		//	pdfAction(c1,ps);
	}
	//if(nPhiBins%nPads!=0)
	//	pdfAction(c1,ps);

	clearPad(c1,nPads);
	hDataPiNHitsFitvsPhiM->RebinX(rebPhi);
	hDataPiNHitsFitvsPhiP->RebinX(rebPhi);
	TH1F *hDataPiNHitsFitPhiM[nPhiBins];
	TH1F *hDataPiNHitsFitPhiP[nPhiBins];
	for(Int_t i=0;i<nPhiBins;i++){
		hDataPiNHitsFitPhiM[i] = (TH1F *)hDataPiNHitsFitvsPhiM->ProjectionY(Form("hDataPiNHitsFitPhiM_PhiBin%d",i+1),i+1,i+1);
		hDataPiNHitsFitPhiP[i] = (TH1F *)hDataPiNHitsFitvsPhiP->ProjectionY(Form("hDataPiNHitsFitPhiP_PhiBin%d",i+1),i+1,i+1);
		setHisto(hDataPiNHitsFitPhiM[i],22,0.8,4,4);
		setHisto(hDataPiNHitsFitPhiP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiNHitsFitPhiM[i]->GetXaxis()->FindBin(30);
		binHi = hDataPiNHitsFitPhiM[i]->GetXaxis()->FindBin(40);
		hDataPiNHitsFitPhiM[i]->Scale(hRcNHitsFitPhiE[i]->Integral(binLow,binHi)*1./hDataPiNHitsFitPhiM[i]->Integral(binLow,binHi));
		hDataPiNHitsFitPhiP[i]->Scale(hRcNHitsFitPhiP[i]->Integral(binLow,binHi)*1./hDataPiNHitsFitPhiP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		hRcNHitsFitPhiE[i]->Draw("histe");
		hRcNHitsFitPhiP[i]->Draw("histesame");
		hDataPiNHitsFitPhiM[i]->Draw("pesame");
		hDataPiNHitsFitPhiP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsFitPhiE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsFitPhiP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNHitsFitPhiM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNHitsFitPhiP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		Float_t PhiLow = hRcNHitsFitvsPhiE->GetXaxis()->GetBinLowEdge(i+1);
		Float_t PhiHi = hRcNHitsFitvsPhiE->GetXaxis()->GetBinLowEdge(i+2);
		drawLatex(0.40,0.95,Form("%3.2f<#phi<%3.2f",PhiLow,PhiHi),42,0.06,1);
		//if(i%nPads==nPads-1)
		//	pdfAction(c1,ps);
	}
	//if(nPhiBins%nPads!=0)
	//	pdfAction(c1,ps);


	hRcNHitsPossvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsPossvsPtE = (TH2F *)hRcNHitsPossvsPtQ->Project3D("ZY");
	hRcNHitsPossvsPtE->SetNameTitle("hRcNHitsPossvsPtE","hRcNHitsPossvsPtE");
	hRcNHitsPossvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsPossvsPtP = (TH2F *)hRcNHitsPossvsPtQ->Project3D("ZY");
	hRcNHitsPossvsPtP->SetNameTitle("hRcNHitsPossvsPtP","hRcNHitsPossvsPtP");
	hDataNHitsPossvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNHitsPossvsPtE = (TH2F *)hDataNHitsPossvsPtQ->Project3D("ZY");
	hDataNHitsPossvsPtE->SetNameTitle("hDataNHitsPossvsPtE","hDataNHitsPossvsPtE");
	hDataNHitsPossvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNHitsPossvsPtP = (TH2F *)hDataNHitsPossvsPtQ->Project3D("ZY");
	hDataNHitsPossvsPtP->SetNameTitle("hDataNHitsPossvsPtP","hDataNHitsPossvsPtP");
	hDataPiNHitsPossvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNHitsPossvsPtM = (TH2F *)hDataPiNHitsPossvsPtQ->Project3D("ZY");
	hDataPiNHitsPossvsPtM->SetNameTitle("hDataPiNHitsPossvsPtM","hDataPiNHitsPossvsPtM");
	hDataPiNHitsPossvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNHitsPossvsPtP = (TH2F *)hDataPiNHitsPossvsPtQ->Project3D("ZY");
	hDataPiNHitsPossvsPtP->SetNameTitle("hDataPiNHitsPossvsPtP","hDataPiNHitsPossvsPtP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsPossvsPtE->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsPossvsPtE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsPossvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsPossvsPtP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNHitsPossvsPtE->GetXaxis()->SetRangeUser(0,5.);
	hDataNHitsPossvsPtE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNHitsPossvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataNHitsPossvsPtP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNHitsPossvsPtM->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNHitsPossvsPtM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNHitsPossvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNHitsPossvsPtP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	TH1F *hRcNHitsPossE[nMomBins];
	TH1F *hRcNHitsPossP[nMomBins];
	TH1F *hDataNHitsPossE[nMomBins];
	TH1F *hDataNHitsPossP[nMomBins];
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcNHitsPossvsPtE->GetXaxis()->FindBin(mom[i]);
		momBinHi = hRcNHitsPossvsPtE->GetXaxis()->FindBin(mom[i+1]);
		hRcNHitsPossE[i] = (TH1F *)hRcNHitsPossvsPtE->ProjectionY(Form("hRcNHitsPossE_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsPossP[i] = (TH1F *)hRcNHitsPossvsPtP->ProjectionY(Form("hRcNHitsPossP_PtBin%d",i+1),momBinLow,momBinHi);
		hDataNHitsPossE[i] = (TH1F *)hDataNHitsPossvsPtE->ProjectionY(Form("hDataNHitsPossE_PtBin%d",i+1),momBinLow,momBinHi);
		hDataNHitsPossP[i] = (TH1F *)hDataNHitsPossvsPtP->ProjectionY(Form("hDataNHitsPossP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsPossE[i]->SetLineColor(1);
		hRcNHitsPossP[i]->SetLineColor(2);
		setHisto(hDataNHitsPossE[i],20,0.8,4,4);
		setHisto(hDataNHitsPossP[i],20,0.8,kViolet,kViolet);
		binLow = hDataNHitsPossE[i]->GetXaxis()->FindBin(40);
		binHi = hDataNHitsPossE[i]->GetXaxis()->FindBin(45);
		// hDataNHitsPossE[i]->Scale(hRcNHitsPossE[i]->Integral(binLow,binHi)*1./hDataNHitsPossE[i]->Integral(binLow,binHi));
		// hDataNHitsPossP[i]->Scale(hRcNHitsPossP[i]->Integral(binLow,binHi)*1./hDataNHitsPossP[i]->Integral(binLow,binHi));
		hDataNHitsPossE[i]->Scale(1./hDataNHitsPossE[i]->Integral());
		hDataNHitsPossP[i]->Scale(1./hDataNHitsPossP[i]->Integral());
		hRcNHitsPossE[i]->Scale(1./hRcNHitsPossE[i]->Integral());
		hRcNHitsPossP[i]->Scale(1./hRcNHitsPossP[i]->Integral());
		c1->cd(i%nPads+1);
		gPad->SetLogy(1);
		hDataNHitsPossP[i]->SetMinimum(1.);
		hDataNHitsPossE[i]->Draw("pe");
		hDataNHitsPossP[i]->Draw("pesame");
		hRcNHitsPossE[i]->Draw("histesame");
		hRcNHitsPossP[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsPossE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsPossP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNHitsPossE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNHitsPossP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	clearPad(c1,nPads);
	TH1F *hRcNHitsPossE1[nPiMomBins];
	TH1F *hRcNHitsPossP1[nPiMomBins];
	TH1F *hDataPiNHitsPossM[nPiMomBins];
	TH1F *hDataPiNHitsPossP[nPiMomBins];
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcNHitsPossvsPtE->GetXaxis()->FindBin(pimom[i]);
		momBinHi = hRcNHitsPossvsPtE->GetXaxis()->FindBin(pimom[i+1]);
		hRcNHitsPossE1[i] = (TH1F *)hRcNHitsPossvsPtE->ProjectionY(Form("hRcNHitsPossE1_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsPossP1[i] = (TH1F *)hRcNHitsPossvsPtP->ProjectionY(Form("hRcNHitsPossP1_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiNHitsPossM[i] = (TH1F *)hDataPiNHitsPossvsPtM->ProjectionY(Form("hDataPiNHitsPossM_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiNHitsPossP[i] = (TH1F *)hDataPiNHitsPossvsPtP->ProjectionY(Form("hDataPiNHitsPossP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsPossE1[i]->SetLineColor(1);
		hRcNHitsPossP1[i]->SetLineColor(2);
		setHisto(hDataPiNHitsPossM[i],22,0.8,4,4);
		setHisto(hDataPiNHitsPossP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiNHitsPossM[i]->GetXaxis()->FindBin(40);
		binHi = hDataPiNHitsPossM[i]->GetXaxis()->FindBin(45);
		// hDataPiNHitsPossM[i]->Scale(hRcNHitsPossE1[i]->Integral(binLow,binHi)*1./hDataPiNHitsPossM[i]->Integral(binLow,binHi));
		// hDataPiNHitsPossP[i]->Scale(hRcNHitsPossP1[i]->Integral(binLow,binHi)*1./hDataPiNHitsPossP[i]->Integral(binLow,binHi));
		hDataPiNHitsPossM[i]->Scale(1./hDataPiNHitsPossM[i]->Integral());
		hDataPiNHitsPossP[i]->Scale(1./hDataPiNHitsPossP[i]->Integral());
		hRcNHitsPossE1[i]->Scale(1./hRcNHitsPossE1[i]->Integral());
		hRcNHitsPossP1[i]->Scale(1./hRcNHitsPossP1[i]->Integral());
		c1->cd(i%nPads+1);
		gPad->SetLogy(1);
		hRcNHitsPossP1[i]->SetMinimum(1.);
		hRcNHitsPossE1[i]->Draw("histe");
		hRcNHitsPossP1[i]->Draw("histesame");
		hDataPiNHitsPossM[i]->Draw("pesame");
		hDataPiNHitsPossP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsPossE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsPossP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNHitsPossM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNHitsPossP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nPiMomBins%nPads!=0)
		pdfAction(c1,ps);

	hRcNHitsPossvsEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsPossvsEtaE = (TH2F *)hRcNHitsPossvsEtaQ->Project3D("ZY");
	hRcNHitsPossvsEtaE->SetNameTitle("hRcNHitsPossvsEtaE","hRcNHitsPossvsEtaE");
	hRcNHitsPossvsEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsPossvsEtaP = (TH2F *)hRcNHitsPossvsEtaQ->Project3D("ZY");
	hRcNHitsPossvsEtaP->SetNameTitle("hRcNHitsPossvsEtaP","hRcNHitsPossvsEtaP");
	hDataNHitsPossvsEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNHitsPossvsEtaE = (TH2F *)hDataNHitsPossvsEtaQ->Project3D("ZY");
	hDataNHitsPossvsEtaE->SetNameTitle("hDataNHitsPossvsEtaE","hDataNHitsPossvsEtaE");
	hDataNHitsPossvsEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNHitsPossvsEtaP = (TH2F *)hDataNHitsPossvsEtaQ->Project3D("ZY");
	hDataNHitsPossvsEtaP->SetNameTitle("hDataNHitsPossvsEtaP","hDataNHitsPossvsEtaP");
	hDataPiNHitsPossvsEtaQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNHitsPossvsEtaM = (TH2F *)hDataPiNHitsPossvsEtaQ->Project3D("ZY");
	hDataPiNHitsPossvsEtaM->SetNameTitle("hDataPiNHitsPossvsEtaM","hDataPiNHitsPossvsEtaM");
	hDataPiNHitsPossvsEtaQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNHitsPossvsEtaP = (TH2F *)hDataPiNHitsPossvsEtaQ->Project3D("ZY");
	hDataPiNHitsPossvsEtaP->SetNameTitle("hDataPiNHitsPossvsEtaP","hDataPiNHitsPossvsEtaP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsPossvsEtaE->GetXaxis()->SetRangeUser(-1.25,1.25);
	hRcNHitsPossvsEtaE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsPossvsEtaP->GetXaxis()->SetRangeUser(-1.25,1.25);
	hRcNHitsPossvsEtaP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNHitsPossvsEtaE->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataNHitsPossvsEtaE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNHitsPossvsEtaP->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataNHitsPossvsEtaP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNHitsPossvsEtaM->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataPiNHitsPossvsEtaM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNHitsPossvsEtaP->GetXaxis()->SetRangeUser(-1.25,1.25);
	hDataPiNHitsPossvsEtaP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//dfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);

	clearPad(c1,nPads);
	hRcNHitsPossvsEtaE->RebinX(rebEta);
	hRcNHitsPossvsEtaP->RebinX(rebEta);
	hDataNHitsPossvsEtaE->RebinX(rebEta);
	hDataNHitsPossvsEtaP->RebinX(rebEta);
	TH1F *hRcNHitsPossEtaE[nEtaBins];
	TH1F *hRcNHitsPossEtaP[nEtaBins];
	TH1F *hDataNHitsPossEtaE[nEtaBins];
	TH1F *hDataNHitsPossEtaP[nEtaBins];
	for(Int_t i=nMinEtaBin;i<=nMaxEtaBin;i++){
		hRcNHitsPossEtaE[i] = (TH1F *)hRcNHitsPossvsEtaE->ProjectionY(Form("hRcNHitsPossE_EtaBin%d",i-nMinEtaBin+1),i,i);
		hRcNHitsPossEtaP[i] = (TH1F *)hRcNHitsPossvsEtaP->ProjectionY(Form("hRcNHitsPossP_EtaBin%d",i-nMinEtaBin+1),i,i);
		hDataNHitsPossEtaE[i] = (TH1F *)hDataNHitsPossvsEtaE->ProjectionY(Form("hDataNHitsPossEtaE_EtaBin%d",i-nMinEtaBin+1),i,i);
		hDataNHitsPossEtaP[i] = (TH1F *)hDataNHitsPossvsEtaP->ProjectionY(Form("hDataNHitsPossEtaP_EtaBin%d",i-nMinEtaBin+1),i,i);
		hRcNHitsPossEtaE[i]->SetLineColor(1);
		hRcNHitsPossEtaP[i]->SetLineColor(2);
		setHisto(hDataNHitsPossEtaE[i],20,0.8,4,4);
		setHisto(hDataNHitsPossEtaP[i],20,0.8,kViolet,kViolet);
		binLow = hDataNHitsPossEtaE[i]->GetXaxis()->FindBin(40);
		binHi = hDataNHitsPossEtaE[i]->GetXaxis()->FindBin(45);
		// hDataNHitsPossEtaE[i]->Scale(hRcNHitsPossEtaE[i]->Integral(binLow,binHi)*1./hDataNHitsPossEtaE[i]->Integral(binLow,binHi));
		// hDataNHitsPossEtaP[i]->Scale(hRcNHitsPossEtaP[i]->Integral(binLow,binHi)*1./hDataNHitsPossEtaP[i]->Integral(binLow,binHi));
		hDataNHitsPossEtaP[i]->Scale(1./hDataNHitsPossEtaP[i]->Integral());
		hDataNHitsPossEtaE[i]->Scale(1./hDataNHitsPossEtaE[i]->Integral());
		hRcNHitsPossEtaE[i]->Scale(1./hRcNHitsPossEtaE[i]->Integral());
		hRcNHitsPossEtaP[i]->Scale(1./hRcNHitsPossEtaP[i]->Integral());
		c1->cd((i-nMinEtaBin)%nPads+1);
		gPad->SetLogy(1);
		hRcNHitsPossEtaP[i]->SetMinimum(1.);;
		hRcNHitsPossEtaE[i]->Draw("histe");
		hRcNHitsPossEtaP[i]->Draw("histesame");
		hDataNHitsPossEtaE[i]->Draw("pesame");
		hDataNHitsPossEtaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsPossEtaE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsPossEtaP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNHitsPossEtaE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNHitsPossEtaP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		Float_t EtaLow = hRcNHitsPossvsEtaE->GetXaxis()->GetBinLowEdge(i);
		Float_t EtaHi = hRcNHitsPossvsEtaE->GetXaxis()->GetBinLowEdge(i+1);
		drawLatex(0.40,0.95,Form("%3.1f<#eta<%3.1f",EtaLow,EtaHi),42,0.06,1);
		if((i-nMinEtaBin)%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if((nMaxEtaBin-nMinEtaBin+1)%nPads!=0)
		pdfAction(c1,ps);

	clearPad(c1,nPads);
	hDataPiNHitsPossvsEtaM->RebinX(rebEta);
	hDataPiNHitsPossvsEtaP->RebinX(rebEta);
	TH1F *hDataPiNHitsPossEtaM[nEtaBins];
	TH1F *hDataPiNHitsPossEtaP[nEtaBins];
	for(Int_t i=nMinEtaBin;i<=nMaxEtaBin;i++){
		hDataPiNHitsPossEtaM[i] = (TH1F *)hDataPiNHitsPossvsEtaM->ProjectionY(Form("hDataPiNHitsPossEtaM_EtaBin%d",i-nMinEtaBin+1),i,i);
		hDataPiNHitsPossEtaP[i] = (TH1F *)hDataPiNHitsPossvsEtaP->ProjectionY(Form("hDataPiNHitsPossEtaP_EtaBin%d",i-nMinEtaBin+1),i,i);
		setHisto(hDataPiNHitsPossEtaM[i],22,0.8,4,4);
		setHisto(hDataPiNHitsPossEtaP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiNHitsPossEtaM[i]->GetXaxis()->FindBin(40);
		binHi = hDataPiNHitsPossEtaM[i]->GetXaxis()->FindBin(45);
		hDataPiNHitsPossEtaM[i]->Scale(hRcNHitsPossEtaE[i]->Integral(binLow,binHi)*1./hDataPiNHitsPossEtaM[i]->Integral(binLow,binHi));
		hDataPiNHitsPossEtaP[i]->Scale(hRcNHitsPossEtaP[i]->Integral(binLow,binHi)*1./hDataPiNHitsPossEtaP[i]->Integral(binLow,binHi));
		c1->cd((i-nMinEtaBin)%nPads+1);
		gPad->SetLogy(1);
		hRcNHitsPossEtaE[i]->Draw("histe");
		hRcNHitsPossEtaP[i]->Draw("histesame");
		hDataPiNHitsPossEtaM[i]->Draw("pesame");
		hDataPiNHitsPossEtaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsPossEtaE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsPossEtaP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNHitsPossEtaM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNHitsPossEtaP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		Float_t EtaLow = hRcNHitsPossvsEtaE->GetXaxis()->GetBinLowEdge(i);
		Float_t EtaHi = hRcNHitsPossvsEtaE->GetXaxis()->GetBinLowEdge(i+1);
		drawLatex(0.40,0.95,Form("%3.1f<#eta<%3.1f",EtaLow,EtaHi),42,0.06,1);
		//if((i-nMinEtaBin)%nPads==nPads-1)
		//	pdfAction(c1,ps);
	}
	//if((nMaxEtaBin-nMinEtaBin+1)%nPads!=0)
	//	pdfAction(c1,ps);

	hRcNHitsPossvsPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsPossvsPhiE = (TH2F *)hRcNHitsPossvsPhiQ->Project3D("ZY");
	hRcNHitsPossvsPhiE->SetNameTitle("hRcNHitsPossvsPhiE","hRcNHitsPossvsPhiE");
	hRcNHitsPossvsPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsPossvsPhiP = (TH2F *)hRcNHitsPossvsPhiQ->Project3D("ZY");
	hRcNHitsPossvsPhiP->SetNameTitle("hRcNHitsPossvsPhiP","hRcNHitsPossvsPhiP");
	hDataNHitsPossvsPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNHitsPossvsPhiE = (TH2F *)hDataNHitsPossvsPhiQ->Project3D("ZY");
	hDataNHitsPossvsPhiE->SetNameTitle("hDataNHitsPossvsPhiE","hDataNHitsPossvsPhiE");
	hDataNHitsPossvsPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNHitsPossvsPhiP = (TH2F *)hDataNHitsPossvsPhiQ->Project3D("ZY");
	hDataNHitsPossvsPhiP->SetNameTitle("hDataNHitsPossvsPhiP","hDataNHitsPossvsPhiP");
	hDataPiNHitsPossvsPhiQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNHitsPossvsPhiM = (TH2F *)hDataPiNHitsPossvsPhiQ->Project3D("ZY");
	hDataPiNHitsPossvsPhiM->SetNameTitle("hDataPiNHitsPossvsPhiM","hDataPiNHitsPossvsPhiM");
	hDataPiNHitsPossvsPhiQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNHitsPossvsPhiP = (TH2F *)hDataPiNHitsPossvsPhiQ->Project3D("ZY");
	hDataPiNHitsPossvsPhiP->SetNameTitle("hDataPiNHitsPossvsPhiP","hDataPiNHitsPossvsPhiP");


	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsPossvsPhiE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsPossvsPhiP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNHitsPossvsPhiE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNHitsPossvsPhiP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNHitsPossvsPhiM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNHitsPossvsPhiP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);

	clearPad(c1,nPads);
	hRcNHitsPossvsPhiE->RebinX(rebPhi);
	hRcNHitsPossvsPhiP->RebinX(rebPhi);
	hDataNHitsPossvsPhiE->RebinX(rebPhi);
	hDataNHitsPossvsPhiP->RebinX(rebPhi);
	TH1F *hRcNHitsPossPhiE[nPhiBins];
	TH1F *hRcNHitsPossPhiP[nPhiBins];
	TH1F *hDataNHitsPossPhiE[nPhiBins];
	TH1F *hDataNHitsPossPhiP[nPhiBins];
	for(Int_t i=0;i<nPhiBins;i++){
		hRcNHitsPossPhiE[i] = (TH1F *)hRcNHitsPossvsPhiE->ProjectionY(Form("hRcNHitsPossE_PhiBin%d",i+1),i+1,i+1);
		hRcNHitsPossPhiP[i] = (TH1F *)hRcNHitsPossvsPhiP->ProjectionY(Form("hRcNHitsPossP_PhiBin%d",i+1),i+1,i+1);
		hDataNHitsPossPhiE[i] = (TH1F *)hDataNHitsPossvsPhiE->ProjectionY(Form("hDataNHitsPossPhiE_PhiBin%d",i+1),i+1,i+1);
		hDataNHitsPossPhiP[i] = (TH1F *)hDataNHitsPossvsPhiP->ProjectionY(Form("hDataNHitsPossPhiP_PhiBin%d",i+1),i+1,i+1);
		hRcNHitsPossPhiE[i]->SetLineColor(1);
		hRcNHitsPossPhiP[i]->SetLineColor(2);
		setHisto(hDataNHitsPossPhiE[i],20,0.8,4,4);
		setHisto(hDataNHitsPossPhiP[i],20,0.8,kViolet,kViolet);
		binLow = hDataNHitsPossPhiE[i]->GetXaxis()->FindBin(40);
		binHi = hDataNHitsPossPhiE[i]->GetXaxis()->FindBin(45);
		hDataNHitsPossPhiE[i]->Scale(hRcNHitsPossPhiE[i]->Integral(binLow,binHi)*1./hDataNHitsPossPhiE[i]->Integral(binLow,binHi));
		hDataNHitsPossPhiP[i]->Scale(hRcNHitsPossPhiP[i]->Integral(binLow,binHi)*1./hDataNHitsPossPhiP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		gPad->SetLogy(1);
		hRcNHitsPossPhiP[i]->SetMinimum(1.);;
		hRcNHitsPossPhiE[i]->Draw("histe");
		hRcNHitsPossPhiP[i]->Draw("histesame");
		hDataNHitsPossPhiE[i]->Draw("pesame");
		hDataNHitsPossPhiP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsPossPhiE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsPossPhiP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNHitsPossPhiE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNHitsPossPhiP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		Float_t PhiLow = hRcNHitsPossvsPhiE->GetXaxis()->GetBinLowEdge(i+1);
		Float_t PhiHi = hRcNHitsPossvsPhiE->GetXaxis()->GetBinLowEdge(i+2);
		drawLatex(0.40,0.95,Form("%3.2f<#phi<%3.2f",PhiLow,PhiHi),42,0.06,1);
		//if(i%nPads==nPads-1)
		//	pdfAction(c1,ps);
	}
	//if(nPhiBins%nPads!=0)
	//	pdfAction(c1,ps);

	clearPad(c1,nPads);
	hDataPiNHitsPossvsPhiM->RebinX(rebPhi);
	hDataPiNHitsPossvsPhiP->RebinX(rebPhi);
	TH1F *hDataPiNHitsPossPhiM[nPhiBins];
	TH1F *hDataPiNHitsPossPhiP[nPhiBins];
	for(Int_t i=0;i<nPhiBins;i++){
		hDataPiNHitsPossPhiM[i] = (TH1F *)hDataPiNHitsPossvsPhiM->ProjectionY(Form("hDataPiNHitsPossPhiM_PhiBin%d",i+1),i+1,i+1);
		hDataPiNHitsPossPhiP[i] = (TH1F *)hDataPiNHitsPossvsPhiP->ProjectionY(Form("hDataPiNHitsPossPhiP_PhiBin%d",i+1),i+1,i+1);
		setHisto(hDataPiNHitsPossPhiM[i],22,0.8,4,4);
		setHisto(hDataPiNHitsPossPhiP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiNHitsPossPhiM[i]->GetXaxis()->FindBin(40);
		binHi = hDataPiNHitsPossPhiM[i]->GetXaxis()->FindBin(45);
		hDataPiNHitsPossPhiM[i]->Scale(hRcNHitsPossPhiE[i]->Integral(binLow,binHi)*1./hDataPiNHitsPossPhiM[i]->Integral(binLow,binHi));
		hDataPiNHitsPossPhiP[i]->Scale(hRcNHitsPossPhiP[i]->Integral(binLow,binHi)*1./hDataPiNHitsPossPhiP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		gPad->SetLogy(1);
		hRcNHitsPossPhiE[i]->Draw("histe");
		hRcNHitsPossPhiP[i]->Draw("histesame");
		hDataPiNHitsPossPhiM[i]->Draw("pesame");
		hDataPiNHitsPossPhiP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsPossPhiE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsPossPhiP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNHitsPossPhiM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNHitsPossPhiP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		Float_t PhiLow = hRcNHitsPossvsPhiE->GetXaxis()->GetBinLowEdge(i+1);
		Float_t PhiHi = hRcNHitsPossvsPhiE->GetXaxis()->GetBinLowEdge(i+2);
		drawLatex(0.40,0.95,Form("%3.2f<#phi<%3.2f",PhiLow,PhiHi),42,0.06,1);
		//if(i%nPads==nPads-1)
		//	pdfAction(c1,ps);
	}
	//if(nPhiBins%nPads!=0)
	//	pdfAction(c1,ps);


	hRcNHitsDedxvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNHitsDedxvsPtE = (TH2F *)hRcNHitsDedxvsPtQ->Project3D("ZY");
	hRcNHitsDedxvsPtE->SetNameTitle("hRcNHitsDedxvsPtE","hRcNHitsDedxvsPtE");
	hRcNHitsDedxvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNHitsDedxvsPtP = (TH2F *)hRcNHitsDedxvsPtQ->Project3D("ZY");
	hRcNHitsDedxvsPtP->SetNameTitle("hRcNHitsDedxvsPtP","hRcNHitsDedxvsPtP");
	hDataNHitsDedxvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNHitsDedxvsPtE = (TH2F *)hDataNHitsDedxvsPtQ->Project3D("ZY");
	hDataNHitsDedxvsPtE->SetNameTitle("hDataNHitsDedxvsPtE","hDataNHitsDedxvsPtE");
	hDataNHitsDedxvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNHitsDedxvsPtP = (TH2F *)hDataNHitsDedxvsPtQ->Project3D("ZY");
	hDataNHitsDedxvsPtP->SetNameTitle("hDataNHitsDedxvsPtP","hDataNHitsDedxvsPtP");
	hDataPiNHitsDedxvsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNHitsDedxvsPtM = (TH2F *)hDataPiNHitsDedxvsPtQ->Project3D("ZY");
	hDataPiNHitsDedxvsPtM->SetNameTitle("hDataPiNHitsDedxvsPtM","hDataPiNHitsDedxvsPtM");
	hDataPiNHitsDedxvsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNHitsDedxvsPtP = (TH2F *)hDataPiNHitsDedxvsPtQ->Project3D("ZY");
	hDataPiNHitsDedxvsPtP->SetNameTitle("hDataPiNHitsDedxvsPtP","hDataPiNHitsDedxvsPtP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNHitsDedxvsPtE->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsDedxvsPtE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNHitsDedxvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hRcNHitsDedxvsPtP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNHitsDedxvsPtE->GetXaxis()->SetRangeUser(0,5.);
	hDataNHitsDedxvsPtE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNHitsDedxvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataNHitsDedxvsPtP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNHitsDedxvsPtM->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNHitsDedxvsPtM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNHitsDedxvsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNHitsDedxvsPtP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	TH1F *hRcNHitsDedxE[nMomBins];
	TH1F *hRcNHitsDedxP[nMomBins];
	TH1F *hDataNHitsDedxE[nMomBins];
	TH1F *hDataNHitsDedxP[nMomBins];
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcNHitsDedxvsPtE->GetXaxis()->FindBin(mom[i]);
		momBinHi = hRcNHitsDedxvsPtE->GetXaxis()->FindBin(mom[i+1]);
		hRcNHitsDedxE[i] = (TH1F *)hRcNHitsDedxvsPtE->ProjectionY(Form("hRcNHitsDedxE_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsDedxP[i] = (TH1F *)hRcNHitsDedxvsPtP->ProjectionY(Form("hRcNHitsDedxP_PtBin%d",i+1),momBinLow,momBinHi);
		hDataNHitsDedxE[i] = (TH1F *)hDataNHitsDedxvsPtE->ProjectionY(Form("hDataNHitsDedxE_PtBin%d",i+1),momBinLow,momBinHi);
		hDataNHitsDedxP[i] = (TH1F *)hDataNHitsDedxvsPtP->ProjectionY(Form("hDataNHitsDedxP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsDedxE[i]->SetLineColor(1);
		hRcNHitsDedxP[i]->SetLineColor(2);
		setHisto(hDataNHitsDedxE[i],20,0.8,4,4);
		setHisto(hDataNHitsDedxP[i],20,0.8,kViolet,kViolet);
		binLow = hDataNHitsDedxE[i]->GetXaxis()->FindBin(20);
		binHi = hDataNHitsDedxE[i]->GetXaxis()->FindBin(30);
		// hDataNHitsDedxE[i]->Scale(hRcNHitsDedxE[i]->Integral(binLow,binHi)*1./hDataNHitsDedxE[i]->Integral(binLow,binHi));
		// hDataNHitsDedxP[i]->Scale(hRcNHitsDedxP[i]->Integral(binLow,binHi)*1./hDataNHitsDedxP[i]->Integral(binLow,binHi));
		hDataNHitsDedxE[i]->Scale(1./hDataNHitsDedxE[i]->Integral());
		hDataNHitsDedxP[i]->Scale(1./hDataNHitsDedxP[i]->Integral());
		hRcNHitsDedxE[i]->Scale(1./hRcNHitsDedxE[i]->Integral());
		hRcNHitsDedxP[i]->Scale(1./hRcNHitsDedxP[i]->Integral());
		c1->cd(i%nPads+1);
		gPad->SetLogy(0);
		hDataNHitsDedxP[i]->SetMinimum(1.);
		hDataNHitsDedxE[i]->Draw("pe");
		hDataNHitsDedxP[i]->Draw("pesame");
		hRcNHitsDedxE[i]->Draw("histesame");
		hRcNHitsDedxP[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.2);
		leg->SetX2NDC(0.45);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsDedxE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsDedxP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNHitsDedxE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNHitsDedxP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.15,0.85,Form("%3.1f<p_{T}<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	TH1F *hRcNHitsDedxE1[nPiMomBins];
	TH1F *hRcNHitsDedxP1[nPiMomBins];
	TH1F *hDataPiNHitsDedxM[nPiMomBins];
	TH1F *hDataPiNHitsDedxP[nPiMomBins];
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcNHitsDedxvsPtE->GetXaxis()->FindBin(pimom[i]);
		momBinHi = hRcNHitsDedxvsPtE->GetXaxis()->FindBin(pimom[i+1]);
		hRcNHitsDedxE1[i] = (TH1F *)hRcNHitsDedxvsPtE->ProjectionY(Form("hRcNHitsDedxE1_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsDedxP1[i] = (TH1F *)hRcNHitsDedxvsPtP->ProjectionY(Form("hRcNHitsDedxP1_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiNHitsDedxM[i] = (TH1F *)hDataPiNHitsDedxvsPtM->ProjectionY(Form("hDataPiNHitsDedxM_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiNHitsDedxP[i] = (TH1F *)hDataPiNHitsDedxvsPtP->ProjectionY(Form("hDataPiNHitsDedxP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcNHitsDedxE1[i]->SetLineColor(1);
		hRcNHitsDedxP1[i]->SetLineColor(2);
		setHisto(hDataPiNHitsDedxM[i],22,0.8,4,4);
		setHisto(hDataPiNHitsDedxP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiNHitsDedxM[i]->GetXaxis()->FindBin(20);
		binHi = hDataPiNHitsDedxM[i]->GetXaxis()->FindBin(30);
		hDataPiNHitsDedxM[i]->Scale(hRcNHitsDedxE1[i]->Integral(binLow,binHi)*1./hDataPiNHitsDedxM[i]->Integral(binLow,binHi));
		hDataPiNHitsDedxP[i]->Scale(hRcNHitsDedxP1[i]->Integral(binLow,binHi)*1./hDataPiNHitsDedxP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		hDataPiNHitsDedxP[i]->SetMinimum(1.);
		hDataPiNHitsDedxM[i]->Draw("pe");
		hDataPiNHitsDedxP[i]->Draw("pesame");
		hRcNHitsDedxE1[i]->Draw("histesame");
		hRcNHitsDedxP1[i]->Draw("histesame");
		leg->Clear();
		leg->SetX1NDC(0.6);
		leg->SetX2NDC(0.85);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNHitsDedxE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNHitsDedxP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNHitsDedxM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNHitsDedxP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		//if(i%nPads==nPads-1)
			//pdfAction(c1,ps);
	}
	//if(nPiMomBins%nPads!=0)
		//pdfAction(c1,ps);

	hRcDcavsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcDcavsPtE = (TH2F *)hRcDcavsPtQ->Project3D("ZY");
	hRcDcavsPtE->SetNameTitle("hRcDcavsPtE","hRcDcavsPtE");
	hRcDcavsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcDcavsPtP = (TH2F *)hRcDcavsPtQ->Project3D("ZY");
	hRcDcavsPtP->SetNameTitle("hRcDcavsPtP","hRcDcavsPtP");
	hDataDcavsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataDcavsPtE = (TH2F *)hDataDcavsPtQ->Project3D("ZY");
	hDataDcavsPtE->SetNameTitle("hDataDcavsPtE","hDataDcavsPtE");
	hDataDcavsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataDcavsPtP = (TH2F *)hDataDcavsPtQ->Project3D("ZY");
	hDataDcavsPtP->SetNameTitle("hDataDcavsPtP","hDataDcavsPtP");
	hDataPiDcavsPtQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiDcavsPtM = (TH2F *)hDataPiDcavsPtQ->Project3D("ZY");
	hDataPiDcavsPtM->SetNameTitle("hDataPiDcavsPtM","hDataPiDcavsPtM");
	hDataPiDcavsPtQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiDcavsPtP = (TH2F *)hDataPiDcavsPtQ->Project3D("ZY");
	hDataPiDcavsPtP->SetNameTitle("hDataPiDcavsPtP","hDataPiDcavsPtP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcDcavsPtE->GetXaxis()->SetRangeUser(0,5.);
	hRcDcavsPtE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcDcavsPtP->GetXaxis()->SetRangeUser(0,5.);
	hRcDcavsPtP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataDcavsPtE->GetXaxis()->SetRangeUser(0,5.);
	hDataDcavsPtE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataDcavsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataDcavsPtP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiDcavsPtM->GetXaxis()->SetRangeUser(0,5.);
	hDataPiDcavsPtM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiDcavsPtP->GetXaxis()->SetRangeUser(0,5.);
	hDataPiDcavsPtP->Draw("colz");
	drawLatex(textX,textY,"Data. #pi^{+}",textFont,textSize,textColor);
	//pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	TH1F *hRcDcaE[nMomBins];
	TH1F *hRcDcaP[nMomBins];
	TH1F *hDataDcaE[nMomBins];
	TH1F *hDataDcaP[nMomBins];
	rebY = 10;
	hRcDcavsPtE->RebinY(rebY);
	hRcDcavsPtP->RebinY(rebY);
	hDataDcavsPtE->RebinY(rebY);
	hDataDcavsPtP->RebinY(rebY);
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcDcavsPtE->GetXaxis()->FindBin(mom[i]);
		momBinHi = hRcDcavsPtE->GetXaxis()->FindBin(mom[i+1]);
		hRcDcaE[i] = (TH1F *)hRcDcavsPtE->ProjectionY(Form("hRcDcaE_PtBin%d",i+1),momBinLow,momBinHi);
		hRcDcaP[i] = (TH1F *)hRcDcavsPtP->ProjectionY(Form("hRcDcaP_PtBin%d",i+1),momBinLow,momBinHi);
		hDataDcaE[i] = (TH1F *)hDataDcavsPtE->ProjectionY(Form("hDataDcaE_PtBin%d",i+1),momBinLow,momBinHi);
		hDataDcaP[i] = (TH1F *)hDataDcavsPtP->ProjectionY(Form("hDataDcaP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcDcaE[i]->SetLineColor(1);
		hRcDcaP[i]->SetLineColor(2);
		setHisto(hDataDcaE[i],20,0.8,4,4);
		setHisto(hDataDcaP[i],20,0.8,kViolet,kViolet);
		binLow = hDataDcaE[i]->GetXaxis()->FindBin(0.1);
		binHi = hDataDcaE[i]->GetXaxis()->FindBin(1.0);
		// hDataDcaE[i]->Scale(hRcDcaE[i]->Integral(binLow,binHi)*1./hDataDcaE[i]->Integral(binLow,binHi));
		// hDataDcaP[i]->Scale(hRcDcaP[i]->Integral(binLow,binHi)*1./hDataDcaP[i]->Integral(binLow,binHi));
		hDataDcaE[i]->Scale(1./hDataDcaE[i]->Integral());
		hDataDcaP[i]->Scale(1./hDataDcaP[i]->Integral());
		hRcDcaE[i]->Scale(1./hRcDcaE[i]->Integral());
		hRcDcaP[i]->Scale(1./hRcDcaP[i]->Integral());
		c1->cd(i%nPads+1);
		gPad->SetLogy(0);
		hRcDcaP[i]->SetMinimum(1.);
		hRcDcaP[i]->GetXaxis()->SetRangeUser(0,1.);
		hRcDcaE[i]->Draw("histe");
		hRcDcaP[i]->Draw("histesame");
		hDataDcaE[i]->Draw("pesame");
		hDataDcaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.6);
		leg->SetX2NDC(0.85);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcDcaE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcDcaP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataDcaE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataDcaP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	TH1F *hRcDcaE1[nPiMomBins];
	TH1F *hRcDcaP1[nPiMomBins];
	TH1F *hDataPiDcaM[nPiMomBins];
	TH1F *hDataPiDcaP[nPiMomBins];
	hDataPiDcavsPtM->RebinY(rebY);
	hDataPiDcavsPtP->RebinY(rebY);
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcDcavsPtE->GetXaxis()->FindBin(pimom[i]);
		momBinHi = hRcDcavsPtE->GetXaxis()->FindBin(pimom[i+1]);
		hRcDcaE1[i] = (TH1F *)hRcDcavsPtE->ProjectionY(Form("hRcDcaE1_PtBin%d",i+1),momBinLow,momBinHi);
		hRcDcaP1[i] = (TH1F *)hRcDcavsPtP->ProjectionY(Form("hRcDcaP1_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiDcaM[i] = (TH1F *)hDataPiDcavsPtM->ProjectionY(Form("hDataPiDcaM_PtBin%d",i+1),momBinLow,momBinHi);
		hDataPiDcaP[i] = (TH1F *)hDataPiDcavsPtP->ProjectionY(Form("hDataPiDcaP_PtBin%d",i+1),momBinLow,momBinHi);
		hRcDcaE1[i]->SetLineColor(1);
		hRcDcaP1[i]->SetLineColor(2);
		setHisto(hDataPiDcaM[i],22,0.8,4,4);
		setHisto(hDataPiDcaP[i],22,0.8,kViolet,kViolet);
		binLow = hDataPiDcaM[i]->GetXaxis()->FindBin(0.1);
		binHi = hDataPiDcaM[i]->GetXaxis()->FindBin(1.0);
		hDataPiDcaM[i]->Scale(hRcDcaE1[i]->Integral(binLow,binHi)*1./hDataPiDcaM[i]->Integral(binLow,binHi));
		hDataPiDcaP[i]->Scale(hRcDcaP1[i]->Integral(binLow,binHi)*1./hDataPiDcaP[i]->Integral(binLow,binHi));
		c1->cd(i%nPads+1);
		hRcDcaP1[i]->SetMinimum(1.);
		hRcDcaP1[i]->GetXaxis()->SetRangeUser(0,1.);
		hRcDcaE1[i]->Draw("histe");
		hRcDcaP1[i]->Draw("histesame");
		hDataPiDcaM[i]->Draw("pesame");
		hDataPiDcaP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.6);
		leg->SetX2NDC(0.85);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcDcaE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcDcaP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiDcaM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiDcaP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p_{T}<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		/*if(i%nPads==nPads-1)
			pdfAction(c1,ps);*/
	}
	/*if(nPiMomBins%nPads!=0)
		pdfAction(c1,ps);*/

	hRcDedxvsPQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcDedxvsPE = (TH2F *)hRcDedxvsPQ->Project3D("ZY");
	hRcDedxvsPE->SetNameTitle("hRcDedxvsPE","hRcDedxvsPE");
	hRcDedxvsPQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcDedxvsPP = (TH2F *)hRcDedxvsPQ->Project3D("ZY");
	hRcDedxvsPP->SetNameTitle("hRcDedxvsPP","hRcDedxvsPP");
	hDataDedxvsPQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataDedxvsPE = (TH2F *)hDataDedxvsPQ->Project3D("ZY");
	hDataDedxvsPE->SetNameTitle("hDataDedxvsPE","hDataDedxvsPE");
	hDataDedxvsPQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataDedxvsPP = (TH2F *)hDataDedxvsPQ->Project3D("ZY");
	hDataDedxvsPP->SetNameTitle("hDataDedxvsPP","hDataDedxvsPP");
	hDataPiDedxvsPQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiDedxvsPM = (TH2F *)hDataPiDedxvsPQ->Project3D("ZY");
	hDataPiDedxvsPM->SetNameTitle("hDataPiDedxvsPM","hDataPiDedxvsPM");
	hDataPiDedxvsPQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiDedxvsPP = (TH2F *)hDataPiDedxvsPQ->Project3D("ZY");
	hDataPiDedxvsPP->SetNameTitle("hDataPiDedxvsPP","hDataPiDedxvsPP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcDedxvsPE->GetXaxis()->SetRangeUser(0,5.);
	hRcDedxvsPE->GetYaxis()->SetRangeUser(2.5,5.5);
	hRcDedxvsPE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcDedxvsPP->GetXaxis()->SetRangeUser(0,5.);
	hRcDedxvsPP->GetYaxis()->SetRangeUser(2.5,5.5);
	hRcDedxvsPP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataDedxvsPE->GetXaxis()->SetRangeUser(0,5.);
	hDataDedxvsPE->GetYaxis()->SetRangeUser(2.5,5.5);
	hDataDedxvsPE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataDedxvsPP->GetXaxis()->SetRangeUser(0,5.);
	hDataDedxvsPP->GetYaxis()->SetRangeUser(2.5,5.5);
	hDataDedxvsPP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiDedxvsPM->GetXaxis()->SetRangeUser(0,5.);
	hDataPiDedxvsPM->GetYaxis()->SetRangeUser(1.,5.5);
	hDataPiDedxvsPM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiDedxvsPP->GetXaxis()->SetRangeUser(0,5.);
	hDataPiDedxvsPP->GetYaxis()->SetRangeUser(1.,5.5);
	hDataPiDedxvsPP->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	TH1F *hRcDedxE[nMomBins];
	TH1F *hRcDedxP[nMomBins];
	TH1F *hDataDedxE[nMomBins];
	TH1F *hDataDedxP[nMomBins];
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcDedxvsPE->GetXaxis()->FindBin(mom[i]);
		momBinHi = hRcDedxvsPE->GetXaxis()->FindBin(mom[i+1]);
		hRcDedxE[i] = (TH1F *)hRcDedxvsPE->ProjectionY(Form("hRcDedxE_PBin%d",i+1),momBinLow,momBinHi);
		hRcDedxP[i] = (TH1F *)hRcDedxvsPP->ProjectionY(Form("hRcDedxP_PBin%d",i+1),momBinLow,momBinHi);
		hDataDedxE[i] = (TH1F *)hDataDedxvsPE->ProjectionY(Form("hDataDedxE_PBin%d",i+1),momBinLow,momBinHi);
		hDataDedxP[i] = (TH1F *)hDataDedxvsPP->ProjectionY(Form("hDataDedxP_PBin%d",i+1),momBinLow,momBinHi);
		hRcDedxE[i]->SetLineColor(1);
		hRcDedxP[i]->SetLineColor(2);
		setHisto(hDataDedxE[i],20,0.8,4,4);
		setHisto(hDataDedxP[i],20,0.8,kViolet,kViolet);
		hDataDedxE[i]->Scale(hRcDedxE[i]->GetMaximum()*1./hDataDedxE[i]->GetMaximum());
		hDataDedxP[i]->Scale(hRcDedxP[i]->GetMaximum()*1./hDataDedxP[i]->GetMaximum());
		c1->cd(i%nPads+1);
		gPad->SetLogy(0);
		hRcDedxP[i]->SetMinimum(1.);
		hRcDedxE[i]->Draw("histe");
		hRcDedxP[i]->Draw("histesame");
		hDataDedxE[i]->Draw("pesame");
		hDataDedxP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.6);
		leg->SetX2NDC(0.85);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcDedxE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcDedxP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataDedxE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataDedxP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	TH1F *hRcDedxE1[nPiMomBins];
	TH1F *hRcDedxP1[nPiMomBins];
	TH1F *hDataPiDedxM[nPiMomBins];
	TH1F *hDataPiDedxP[nPiMomBins];
	hRcDedxvsPE->GetYaxis()->SetRangeUser(1.,5.5);
	hRcDedxvsPP->GetYaxis()->SetRangeUser(1.,5.5);
	hDataPiDedxvsPM->GetYaxis()->SetRangeUser(1.,5.5);
	hDataPiDedxvsPP->GetYaxis()->SetRangeUser(1.,5.5);
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcDedxvsPE->GetXaxis()->FindBin(pimom[i]);
		momBinHi = hRcDedxvsPE->GetXaxis()->FindBin(pimom[i+1]);
		hRcDedxE1[i] = (TH1F *)hRcDedxvsPE->ProjectionY(Form("hRcDedxE1_PBin%d",i+1),momBinLow,momBinHi);
		hRcDedxP1[i] = (TH1F *)hRcDedxvsPP->ProjectionY(Form("hRcDedxP1_PBin%d",i+1),momBinLow,momBinHi);
		hDataPiDedxM[i] = (TH1F *)hDataPiDedxvsPM->ProjectionY(Form("hDataPiDedxM_PBin%d",i+1),momBinLow,momBinHi);
		hDataPiDedxP[i] = (TH1F *)hDataPiDedxvsPP->ProjectionY(Form("hDataPiDedxP_PBin%d",i+1),momBinLow,momBinHi);
		hRcDedxE1[i]->SetLineColor(1);
		hRcDedxP1[i]->SetLineColor(2);
		setHisto(hDataPiDedxM[i],22,0.8,4,4);
		setHisto(hDataPiDedxP[i],22,0.8,kViolet,kViolet);
		hDataPiDedxM[i]->Scale(hRcDedxE1[i]->GetMaximum()*1./hDataPiDedxM[i]->GetMaximum());
		hDataPiDedxP[i]->Scale(hRcDedxP1[i]->GetMaximum()*1./hDataPiDedxP[i]->GetMaximum());
		c1->cd(i%nPads+1);
		gPad->SetLogy(0);
		hRcDedxP1[i]->SetMinimum(1.);
		hRcDedxE1[i]->Draw("histe");
		hRcDedxP1[i]->Draw("histesame");
		hDataPiDedxM[i]->Draw("pesame");
		hDataPiDedxP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.6);
		leg->SetX2NDC(0.85);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcDedxE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcDedxP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiDedxM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiDedxP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		/*if(i%nPads==nPads-1)
			pdfAction(c1,ps);*/
	}
	/*if(nPiMomBins%nPads!=0)
		pdfAction(c1,ps);*/

	hRcNSigmaEvsPQ->GetXaxis()->SetRange(1,1);
	TH2F *hRcNSigmaEvsPE = (TH2F *)hRcNSigmaEvsPQ->Project3D("ZY");
	hRcNSigmaEvsPE->SetNameTitle("hRcNSigmaEvsPE","hRcNSigmaEvsPE");
	hRcNSigmaEvsPQ->GetXaxis()->SetRange(2,2);
	TH2F *hRcNSigmaEvsPP = (TH2F *)hRcNSigmaEvsPQ->Project3D("ZY");
	hRcNSigmaEvsPP->SetNameTitle("hRcNSigmaEvsPP","hRcNSigmaEvsPP");
	hDataNSigmaEvsPQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataNSigmaEvsPE = (TH2F *)hDataNSigmaEvsPQ->Project3D("ZY");
	hDataNSigmaEvsPE->SetNameTitle("hDataNSigmaEvsPE","hDataNSigmaEvsPE");
	hDataNSigmaEvsPQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataNSigmaEvsPP = (TH2F *)hDataNSigmaEvsPQ->Project3D("ZY");
	hDataNSigmaEvsPP->SetNameTitle("hDataNSigmaEvsPP","hDataNSigmaEvsPP");
	hDataPiNSigmaEvsPQ->GetXaxis()->SetRange(1,1);
	TH2F *hDataPiNSigmaEvsPM = (TH2F *)hDataPiNSigmaEvsPQ->Project3D("ZY");
	hDataPiNSigmaEvsPM->SetNameTitle("hDataPiNSigmaEvsPM","hDataPiNSigmaEvsPM");
	hDataPiNSigmaEvsPQ->GetXaxis()->SetRange(2,2);
	TH2F *hDataPiNSigmaEvsPP = (TH2F *)hDataPiNSigmaEvsPQ->Project3D("ZY");
	hDataPiNSigmaEvsPP->SetNameTitle("hDataPiNSigmaEvsPP","hDataPiNSigmaEvsPP");

	c1->Clear();
	nColumns = 2;
	nRaws = 2;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hRcNSigmaEvsPE->GetXaxis()->SetRangeUser(0,5.);
	hRcNSigmaEvsPE->GetYaxis()->SetRangeUser(-4.,4.);
	hRcNSigmaEvsPE->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hRcNSigmaEvsPP->GetXaxis()->SetRangeUser(0,5.);
	hRcNSigmaEvsPP->GetYaxis()->SetRangeUser(-4.,4.);
	hRcNSigmaEvsPP->Draw("colz");
	drawLatex(textX,textY,"Rec. MC e^{+}",textFont,textSize,textColor);
	c1->cd(3);
	gPad->SetLogz(1);
	hDataNSigmaEvsPE->GetXaxis()->SetRangeUser(0,5.);
	hDataNSigmaEvsPE->GetYaxis()->SetRangeUser(-4.,4.);
	hDataNSigmaEvsPE->Draw("colz");
	drawLatex(textX,textY,"Data e^{-}",textFont,textSize,textColor);
	c1->cd(4);
	gPad->SetLogz(1);
	hDataNSigmaEvsPP->GetXaxis()->SetRangeUser(0,5.);
	hDataNSigmaEvsPP->GetYaxis()->SetRangeUser(-4.,4.);
	hDataNSigmaEvsPP->Draw("colz");
	drawLatex(textX,textY,"Data e^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	clearPad(c1,nPads);
	c1->cd(1);
	gPad->SetLogz(1);
	hDataPiNSigmaEvsPM->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNSigmaEvsPM->GetYaxis()->SetRangeUser(-10,2.);
	hDataPiNSigmaEvsPM->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{-}",textFont,textSize,textColor);
	c1->cd(2);
	gPad->SetLogz(1);
	hDataPiNSigmaEvsPP->GetXaxis()->SetRangeUser(0,5.);
	hDataPiNSigmaEvsPP->GetYaxis()->SetRangeUser(-10,2.);
	hDataPiNSigmaEvsPP->Draw("colz");
	drawLatex(textX,textY,"Data #pi^{+}",textFont,textSize,textColor);
	pdfAction(c1,ps);

	c1->Clear();
	nColumns = 3;
	nRaws = 3;
	nPads = nColumns*nRaws;
	c1->Divide(nColumns,nRaws);
	clearPad(c1,nPads);
	TH1F *hRcNSigmaEE[nMomBins];
	TH1F *hRcNSigmaEP[nMomBins];
	TH1F *hDataNSigmaEE[nMomBins];
	TH1F *hDataNSigmaEP[nMomBins];
	rebY = 2;
	hRcNSigmaEvsPE->RebinY(rebY);
	hRcNSigmaEvsPP->RebinY(rebY);
	hDataNSigmaEvsPE->RebinY(rebY);
	hDataNSigmaEvsPP->RebinY(rebY);
	hRcNSigmaEvsPE->GetYaxis()->SetRangeUser(-4.,4.);
	hRcNSigmaEvsPP->GetYaxis()->SetRangeUser(-4.,4.);
	hDataNSigmaEvsPE->GetYaxis()->SetRangeUser(-4.,4.);
	hDataNSigmaEvsPP->GetYaxis()->SetRangeUser(-4.,4.);
	for(Int_t i=0;i<nMomBins;i++){
		momBinLow = hRcNSigmaEvsPE->GetXaxis()->FindBin(mom[i]);
		momBinHi = hRcNSigmaEvsPE->GetXaxis()->FindBin(mom[i+1]);
		hRcNSigmaEE[i] = (TH1F *)hRcNSigmaEvsPE->ProjectionY(Form("hRcNSigmaEE_PBin%d",i+1),momBinLow,momBinHi);
		hRcNSigmaEP[i] = (TH1F *)hRcNSigmaEvsPP->ProjectionY(Form("hRcNSigmaEP_PBin%d",i+1),momBinLow,momBinHi);
		hDataNSigmaEE[i] = (TH1F *)hDataNSigmaEvsPE->ProjectionY(Form("hDataNSigmaEE_PBin%d",i+1),momBinLow,momBinHi);
		hDataNSigmaEP[i] = (TH1F *)hDataNSigmaEvsPP->ProjectionY(Form("hDataNSigmaEP_PBin%d",i+1),momBinLow,momBinHi);
		hRcNSigmaEE[i]->SetLineColor(1);
		hRcNSigmaEP[i]->SetLineColor(2);
		setHisto(hDataNSigmaEE[i],20,0.8,4,4);
		setHisto(hDataNSigmaEP[i],20,0.8,kViolet,kViolet);
		hDataNSigmaEE[i]->Scale(hRcNSigmaEE[i]->GetMaximum()*1./hDataNSigmaEE[i]->GetMaximum());
		hDataNSigmaEP[i]->Scale(hRcNSigmaEP[i]->GetMaximum()*1./hDataNSigmaEP[i]->GetMaximum());
		c1->cd(i%nPads+1);
		gPad->SetLogy(0);
		hRcNSigmaEP[i]->SetMinimum(1.);
		hRcNSigmaEE[i]->Draw("histe");
		hRcNSigmaEP[i]->Draw("histesame");
		hDataNSigmaEE[i]->Draw("pesame");
		hDataNSigmaEP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.6);
		leg->SetX2NDC(0.85);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNSigmaEE[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNSigmaEP[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataNSigmaEE[i],"Data e^{-}","pl");
		leg->AddEntry(hDataNSigmaEP[i],"Data e^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p<%3.1f (GeV/c)",mom[i],mom[i+1]),42,0.06,1);
		if(i%nPads==nPads-1)
			pdfAction(c1,ps);
	}
	if(nMomBins%nPads!=0)
		pdfAction(c1,ps);

	TH1F *hRcNSigmaEE1[nPiMomBins];
	TH1F *hRcNSigmaEP1[nPiMomBins];
	TH1F *hDataPiNSigmaEM[nPiMomBins];
	TH1F *hDataPiNSigmaEP[nPiMomBins];
	hDataPiNSigmaEvsPM->RebinY(rebY);
	hDataPiNSigmaEvsPP->RebinY(rebY);
	hRcNSigmaEvsPE->GetYaxis()->SetRangeUser(-10.,2.);
	hRcNSigmaEvsPP->GetYaxis()->SetRangeUser(-10.,2.);
	hDataPiNSigmaEvsPM->GetYaxis()->SetRangeUser(-10.,2.);
	hDataPiNSigmaEvsPP->GetYaxis()->SetRangeUser(-10.,2.);
	for(Int_t i=0;i<nPiMomBins;i++){
		momBinLow = hRcNSigmaEvsPE->GetXaxis()->FindBin(pimom[i]);
		momBinHi = hRcNSigmaEvsPE->GetXaxis()->FindBin(pimom[i+1]);
		hRcNSigmaEE1[i] = (TH1F *)hRcNSigmaEvsPE->ProjectionY(Form("hRcNSigmaEE1_PBin%d",i+1),momBinLow,momBinHi);
		hRcNSigmaEP1[i] = (TH1F *)hRcNSigmaEvsPP->ProjectionY(Form("hRcNSigmaEP1_PBin%d",i+1),momBinLow,momBinHi);
		hDataPiNSigmaEM[i] = (TH1F *)hDataPiNSigmaEvsPM->ProjectionY(Form("hDataPiNSigmaEM_PBin%d",i+1),momBinLow,momBinHi);
		hDataPiNSigmaEP[i] = (TH1F *)hDataPiNSigmaEvsPP->ProjectionY(Form("hDataPiNSigmaEP_PBin%d",i+1),momBinLow,momBinHi);
		hRcNSigmaEE1[i]->SetLineColor(1);
		hRcNSigmaEP1[i]->SetLineColor(2);
		setHisto(hDataPiNSigmaEM[i],22,0.8,4,4);
		setHisto(hDataPiNSigmaEP[i],22,0.8,kViolet,kViolet);
		hDataPiNSigmaEM[i]->Scale(hRcNSigmaEE1[i]->GetMaximum()*1./hDataPiNSigmaEM[i]->GetMaximum());
		hDataPiNSigmaEP[i]->Scale(hRcNSigmaEP1[i]->GetMaximum()*1./hDataPiNSigmaEP[i]->GetMaximum());
		c1->cd(i%nPads+1);
		gPad->SetLogy(0);
		hRcNSigmaEP1[i]->SetMinimum(1.);
		hRcNSigmaEE1[i]->Draw("histe");
		hRcNSigmaEP1[i]->Draw("histesame");
		hDataPiNSigmaEM[i]->Draw("pesame");
		hDataPiNSigmaEP[i]->Draw("pesame");
		leg->Clear();
		leg->SetX1NDC(0.18);
		leg->SetX2NDC(0.38);
		leg->SetY1NDC(0.65);
		leg->SetY2NDC(0.85);
		leg->AddEntry(hRcNSigmaEE1[i],"Rec. MC e^{-}","pl");
		leg->AddEntry(hRcNSigmaEP1[i],"Rec. MC e^{+}","pl");
		leg->AddEntry(hDataPiNSigmaEM[i],"Data #pi^{-}","pl");
		leg->AddEntry(hDataPiNSigmaEP[i],"Data #pi^{+}","pl");
		leg->DrawClone("same");
		drawLatex(0.35,0.95,Form("%3.1f<p<%3.1f (GeV/c)",pimom[i],pimom[i+1]),42,0.06,1);
		/*if(i%nPads==nPads-1)
			pdfAction(c1,ps);*/
	}
	/*if(nPiMomBins%nPads!=0)
		pdfAction(c1,ps);*/

	ps->On();
	ps->Close();

	cout<<"End of program !"<<endl;
	return;
}
//___________________________________________________________________
bool Init()
{
	ifstream indata;

	indata.open("./fitPars.dat");
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
	funPosHi = new TF1("funPosHi","[0]*exp(-([1]/x)**[2])+[3]",0.2,60);
	funPosHi->SetParameters(par[0][0],par[0][1],par[0][2],par[0][3]);
	funPosHi->SetLineWidth(1);
	funPosLow = new TF1("funPosLow","[0]*exp(-([1]/x)**[2])+[3]",0.2,60);
	funPosLow->SetParameters(par[1][0],par[1][1],par[1][2],par[1][3]);
	funPosLow->SetLineWidth(1);
	funNegHi = new TF1("funNegHi","[0]*exp(-([1]/x)**[2])+[3]",-60,-0.2);
	funNegHi->SetParameters(par[2][0],par[2][1],par[2][2],par[2][3]);
	funNegHi->SetLineWidth(1);
	funNegLow = new TF1("funNegLow","[0]*exp(-([1]/x)**[2])+[3]",-60,-0.2);
	funNegLow->SetParameters(par[3][0],par[3][1],par[3][2],par[3][3]);
	funNegLow->SetLineWidth(1);

	funNegHiMirror = new TF1("funNegHiMirror","[0]*exp(-([1]/(-1*x))**[2])+[3]",0.2,60);
	funNegHiMirror->SetParameters(par[2][0],par[2][1],par[2][2],par[2][3]);
	funNegHiMirror->SetLineWidth(1);
	funNegLowMirror = new TF1("funNegLowMirror","[0]*exp(-([1]/(-1*x))**[2])+[3]",0.2,60);
	funNegLowMirror->SetParameters(par[3][0],par[3][1],par[3][2],par[3][3]);
	funNegLowMirror->SetLineWidth(1);

	cout<<"Initialization DONE !!!"<<endl;

	return kTRUE;
}
//___________________________________________________________________
TH1F *calRatio(TH1F *hNum, TH1F *hDen, TString name)
{
	if(hNum->GetNbinsX() != hDen->GetNbinsX()){
		cout<<"The bin numbers of \""<<hNum->GetName()<<"\" and \""<<hDen->GetName()<<"\""<<" are not the same !"<<endl;
		return NULL;
	}

	TH1F *hRatio = (TH1F *)hNum->Clone(name.Data());
	hRatio->SetTitle(name.Data());
	hRatio->GetXaxis()->SetTitle(hNum->GetXaxis()->GetTitle());
	hRatio->Reset();
	Int_t nBinsX = hNum->GetNbinsX();
	for(Int_t i=0;i<nBinsX;i++){
		Float_t nNum = hNum->GetBinContent(i+1);
		Float_t nNumErr = hNum->GetBinError(i+1);
		Float_t nDen = hDen->GetBinContent(i+1);
		Float_t nDenErr = hDen->GetBinError(i+1);
		Float_t ratio = 0.;
		Float_t ratioErr = 0.;
		if(nNum>0 && nDen>0){
			ratio = nNum/nDen;
			ratioErr = ratio*sqrt(pow(nNumErr,2)/pow(nNum,2)+pow(nDenErr,2)/pow(nDen,2));
		}
		hRatio->SetBinContent(i+1,ratio);
		hRatio->SetBinError(i+1,ratioErr);
	}
	return hRatio;
};
//___________________________________________________________________
TH2D* histo(TString name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, TString xTitle, TString yTitle)
{
	TH2D *dd = new TH2D(name.Data(),"",100,xlow,xup,100,ylow,yup);
	dd->GetXaxis()->SetTitle(xTitle.Data());
	dd->GetYaxis()->SetTitle(yTitle.Data());

	dd->GetXaxis()->SetTitleSize(0.055);
	dd->GetXaxis()->SetTitleOffset(0.9);
	dd->GetXaxis()->SetLabelSize(0.045);
	dd->GetYaxis()->SetTitleSize(0.055);
	dd->GetYaxis()->SetTitleOffset(1);
	dd->GetYaxis()->SetLabelSize(0.045);
	dd->GetXaxis()->CenterTitle(kTRUE);
	dd->GetYaxis()->CenterTitle(kTRUE);
	dd->GetXaxis()->SetNdivisions(512);
	return dd;
}
//___________________________________________________________________
TLatex* drawLatex(Double_t x, Double_t y, TString text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
	TLatex *latex = new TLatex(x,y,text.Data());
	latex->SetNDC();
	latex->SetTextFont(textFont);
	latex->SetTextSize(textSize);
	latex->SetTextColor(colorIndex);
	latex->Draw("same");
	return latex;
}
//___________________________________________________________________
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
	TLine *l1 = new TLine(xlow,ylow,xup,yup);
	l1->SetLineWidth(lineWidth);
	l1->SetLineColor(lineColor);
	l1->SetLineStyle(lineStyle);
	l1->Draw("same");
	return l1;
}
//___________________________________________________________________
void drawLines(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
	drawLine(xlow,ylow,xup,ylow,lineWidth,lineStyle,lineColor);
	drawLine(xlow,yup,xup,yup,lineWidth,lineStyle,lineColor);
	drawLine(xlow,ylow,xlow,yup,lineWidth,lineStyle,lineColor);
	drawLine(xup,ylow,xup,yup,lineWidth,lineStyle,lineColor);
}
//___________________________________________________________________
void setPad(float left, float right, float top, float bottom)
{
	gPad->SetFillColor(10);
	gPad->SetBorderMode(0);
	gPad->SetBorderSize(0);
	gPad->SetFrameFillColor(10);
	gPad->SetFrameBorderMode(0);
	gPad->SetFrameBorderSize(0);
	gPad->SetLeftMargin(left);
	gPad->SetRightMargin(right);
	gPad->SetTopMargin(top);
	gPad->SetBottomMargin(bottom);
}
//___________________________________________________________________
void clearPad(TCanvas *c, Int_t nPads)
{
	for(Int_t i=0;i<nPads;i++){
		c->cd(i+1);
		gPad->Clear();
	}
}
//___________________________________________________________________
void setHisto(TH1F *h,Int_t MarkerStyle, Float_t MarkerSize, Int_t MarkerColor,Int_t LineColor)
{
	h->SetMarkerStyle(MarkerStyle);
	h->SetMarkerSize(MarkerSize);
	h->SetMarkerColor(MarkerColor);
	h->SetLineColor(LineColor);
};
//___________________________________________________________________
void setHisto(TH1D *h,Int_t MarkerStyle, Float_t MarkerSize, Int_t MarkerColor,Int_t LineColor)
{
	h->SetMarkerStyle(MarkerStyle);
	h->SetMarkerSize(MarkerSize);
	h->SetMarkerColor(MarkerColor);
	h->SetLineColor(LineColor);
};
//___________________________________________________________________
void pdfAction(TCanvas *c, TPDF *ps)
{
	ps->On();
	c->Update();
	c->cd();
	ps->NewPage();
	ps->Off();
};

