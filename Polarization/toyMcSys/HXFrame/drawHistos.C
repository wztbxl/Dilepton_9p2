#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TVirtualPad.h"

//----------------------------------------------------------
TCanvas * drawHistos(TList *list, TString title, TString hTitle, Bool_t setXRange, Double_t rangeMin, Double_t rangeMax, Bool_t setYRange, Double_t min, Double_t max, Bool_t setLogy, Bool_t useLeg, TString *legName, Bool_t drawLeg, TString legTitle, Float_t llx, Float_t lhx, Float_t lly, Float_t lhy,Bool_t drawP, const Float_t size, const Float_t legSize, Bool_t shift, Double_t xshift, Bool_t setMarker, Bool_t setColor, const TString drawOption, const Int_t titleFont)
{
  TCanvas *c = new TCanvas(title.Data(),title.Data(),800,600);
  //SetPadMargin(gPad,0.11,0.11);
  
  TLegend *leg = new TLegend(llx,lly,lhx,lhy);
  if(legTitle.Length()>0)
    leg->SetHeader(legTitle.Data());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(legSize);
  
  Int_t nHistos=list->GetEntries();
  printf("Input %d histograms\n",nHistos);
  title="";
  for(Int_t i=0; i<nHistos; i++)
    {
      TH1 *h = (TH1*)list->At(i);
      if(drawP && setMarker)
	{
	  h->SetMarkerStyle(20);
	  //h->SetMarkerSize(0.9);
	}
      if(setColor) h->SetMarkerColor(gColor[i]);
      if(setColor) h->SetLineColor(gColor[i]);
      
      c->cd();
      if(setLogy) gPad->SetLogy();
      if(!shift)
	{
	  if(i==0)
	    {
	      if(setYRange)  h->GetYaxis()->SetRangeUser(min, max);
	      if(setXRange)  h->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
	      if(hTitle.Length()>0) h->SetTitle(hTitle.Data());
	      title = h->GetTitle();
	      h->SetTitle("");
	      if(drawOption.Length()>0)
		h->DrawCopy(drawOption.Data());
	      else
		{
		  if(drawP)	h->DrawCopy("PE");
		  else	h->DrawCopy("HIST");
		}
	    }
	  else
	    {
	      if(drawOption.Length()>0)
		h->DrawCopy(Form("sames %s",drawOption.Data()));
	      else
		{
		  if(drawP)	h->DrawCopy("sames PE");
		  else	h->DrawCopy("sames HIST");
		}
	    }
	}
      else
	{
	  TH1F *htmp = (TH1F*)h->Clone(Form("%s_tmp",h->GetName()));
	  TGraphErrors *hGrTmp = new TGraphErrors(htmp);
	  hGrTmp->GetXaxis()->SetBinLabel(1,"test");
	  Double_t sign = 0;
	  if(i%2==1) sign = (i+1)/2;
	  else       sign = -1*i/2;
	  offset_x(hGrTmp,xshift*sign);
	  if(i==0)
	    {
	      if(setYRange)  hGrTmp->GetYaxis()->SetRangeUser(min, max);
	      if(setXRange)  hGrTmp->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
	      if(hTitle.Length()>0) hGrTmp->SetTitle(hTitle.Data());
	      title = hGrTmp->GetHistogram()->GetTitle();
	      hGrTmp->SetTitle("");
	      hGrTmp->GetYaxis()->SetTitleOffset(1.1);
	      hGrTmp->Draw("APZ");
	    }
	  else
	    hGrTmp->Draw("SAMES PZ");
	}
      char *legOption = Form("L");
      if(drawP) legOption= Form("PL");
      const char *legLabel = h->GetName();
      if(useLeg)
	legLabel = legName[i].Data();
      leg->AddEntry(h,legLabel,legOption);
    }
  if(drawLeg)
    leg->Draw();

  if(title.Length()>0)
    {
      TPaveText *t1 = GetTitleText(title,size,titleFont);
      t1->Draw();
    }
  return c;
}

//----------------------------------------------------------
TCanvas * drawGraphs(TList *list, TString title, TString hTitle, Bool_t setXRange=kTRUE, Double_t rangeMin=0, Double_t rangeMax=200, Bool_t setYRange=kTRUE, Double_t min=10, Double_t max=1e-10, Bool_t setLogy=kTRUE, Bool_t useLeg=kFALSE, TString *legName=0x0, Bool_t drawLeg=kTRUE, TString legTitle = "", Float_t llx=0.3, Float_t lhx=0.5, Float_t lly=0.6, Float_t lhy=0.8,Bool_t drawP=kTRUE, const Float_t size = 0.04, const Float_t legSize = 0.04, Bool_t shift=kFALSE, Double_t xshift=1, Bool_t setMarker = kTRUE, Bool_t setColor = kTRUE, char *drawOption = "P", Int_t font = 62)
{
  TCanvas *c = new TCanvas(title.Data(),title.Data(),800,600);
  //SetPadMargin(gPad,0.11,0.11);
  
  TLegend *leg = new TLegend(llx,lly,lhx,lhy);
  if(legTitle.Length()>0)
    leg->SetHeader(legTitle.Data());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(legSize);
  
  Int_t nHistos=list->GetEntries();
  printf("Input %d graphs\n",nHistos);
  title="";
  for(Int_t i=0; i<nHistos; i++)
    {
      TGraph *h = (TGraph*)list->At(i);
      if(setMarker) h->SetMarkerStyle(21);
      if(setColor) 
	{
	  h->SetMarkerColor(gColor[i]);
	  h->SetLineColor(gColor[i]);
	}
      
      c->cd();
      if(setLogy) gPad->SetLogy();
      if(!shift)
	{
	  if(i==0)
	    {
	      if(setYRange)  h->GetYaxis()->SetRangeUser(min, max);
	      if(setXRange)  h->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
	      if(hTitle.Length()>0) h->SetTitle(hTitle.Data());
	      title = h->GetHistogram()->GetTitle();
	      h->SetTitle("");
	      h->Draw(Form("A%s",drawOption));
	    }
	  else
	    {
	      h->Draw(Form("sames %s",drawOption));
	    }
	}
      else
	{
	  TGraphErrors *hGrTmp = new TGraphErrors(*(TGraphErrors*)h);
	  Double_t sign = 0;
	  if(i%2==1) sign = (i+1)/2;
	  else       sign = -1*i/2;
	  offset_x(hGrTmp,xshift*sign);
	  if(i==0)
	    {
	      if(setYRange)  hGrTmp->GetYaxis()->SetRangeUser(min, max);
	      if(setXRange)  hGrTmp->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
	      if(hTitle.Length()>0) hGrTmp->SetTitle(hTitle.Data());
	      title = hGrTmp->GetHistogram()->GetTitle();
	      hGrTmp->SetTitle("");
	      hGrTmp->Draw(Form("A%s",drawOption));
	    }
	  else
	    hGrTmp->Draw(Form("sames %s",drawOption));
	}
      const char *legLabel = h->GetName();
      if(useLeg)
	legLabel = legName[i].Data();
      leg->AddEntry(h,legLabel,"PLE");
    }
  if(drawLeg)
    leg->Draw();

  if(title.Length()>0)
    {
      TPaveText *t1 = GetTitleText(title,size,font);
      t1->Draw();
    }
  return c;
}

//----------------------------------------------------------
TCanvas * drawRatios(TList *list, TString title, TString hTitle, Bool_t setXRange=kTRUE, Double_t rangeMin=0, Double_t rangeMax=200, Bool_t setYRange=kTRUE, Double_t min=10, Double_t max=1e-10, Bool_t setLogy=kTRUE, Bool_t useLeg=kFALSE, TString *legName=0x0, Bool_t drawLeg=kTRUE, TString legTitle = "", Float_t llx=0.3, Float_t lhx=0.5, Float_t lly=0.6, Float_t lhy=0.8,Bool_t drawP=kTRUE, const Float_t size = 0.04, const Float_t legSize = 0.04, Bool_t setMarker = kTRUE, Bool_t setColor = kTRUE)
{
  TCanvas *c = new TCanvas(title.Data(),title.Data(),800,600);
  SetPadMargin(gPad,0.11,0.11);
  
  TLegend *leg = new TLegend(llx,lly,lhx,lhy);
  if(legTitle.Length()>0)
    leg->SetHeader(legTitle.Data());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(legSize);
  
  Int_t nHistos=list->GetEntries();
  printf("Input %d histograms\n",nHistos);
  title="";
  for(Int_t i=0; i<nHistos; i++)
    {
      TGraph *h = (TGraph*)list->At(i);
      if(drawP && setMarker)
	{
	  h->SetMarkerStyle(20);
	  //h->SetMarkerSize(0.9);
	  if(setColor) h->SetMarkerColor(gColor[i]);
	}
      if(setMarker && setColor) h->SetLineColor(gColor[i]);
      
      c->cd();
      if(setLogy) gPad->SetLogy();
      if(setYRange)  h->GetYaxis()->SetRangeUser(min, max);
      if(setXRange)  h->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
      if(hTitle.Length()>0) h->SetTitle(hTitle.Data());
      title = h->GetTitle();
      h->SetTitle("");
      //h->SetTitleOffset(1.2,"Y");
      if(i==0)
	{
	  h->Draw("aP");
	}
      else
	{
	  h->Draw("sames PE");
	}

      char *legOption = "L";
      if(drawP) legOption="PL";
      const char *legLabel = h->GetName();
      if(useLeg)
	legLabel = legName[i].Data();
      leg->AddEntry(h,legLabel,legOption);
    }
  if(drawLeg)
    leg->Draw();

  if(title.Length()>0)
    {
      TPaveText *t1 = GetTitleText(title,size);
      t1->Draw();
    }
  return c;
}

//------------------------------------------------
TCanvas *drawSystematics(TList *list, TString title, TString hTitle1,  TString hTitle2, Bool_t setXRange=kTRUE, Double_t rangeMin=0, Double_t rangeMax=200, Bool_t setYRange=kTRUE, Double_t min=10, Double_t max=1e-10, Bool_t setYRange1=kTRUE, Double_t min1=10, Double_t max1=1e-10, Bool_t setLogy=kTRUE, Bool_t useLeg=kFALSE, TString *legName=0x0, Bool_t drawLeg=kTRUE, TString legTitle = "", Float_t llx=0.3, Float_t lhx=0.5, Float_t lly=0.6, Float_t lhy=0.8, Double_t legSize=0.04, Double_t size=0.04)
{

  TCanvas *c = new TCanvas(title.Data(),title.Data(),1100,580);
  c->Divide(2,1);
  c->cd(1);
  if(setLogy) gPad->SetLogy();
  TLegend *leg = new TLegend(llx,lly,lhx,lhy);
  if(legTitle.Length()>0)
    leg->SetHeader(legTitle.Data());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(legSize);

  TH1 *h = (TH1F*)list->At(0);
  h->SetMarkerStyle(21);
  h->SetMarkerColor(1);
  h->SetLineColor(1);
  leg->AddEntry(h,legName[0],"PLE");

  Int_t nHistos=list->GetEntries();
  printf("Input %d histograms\n",nHistos);
  TString title1="";
  for(Int_t i=1; i<nHistos; i++)
    {
      TH1 *htmp = (TH1*)list->At(i);
      if(hTitle1.Length()>0) htmp->SetTitle(hTitle1.Data());
      title1 = htmp->GetTitle();
      htmp->SetTitle("");
      TGraphAsymmErrors *gr = new TGraphAsymmErrors(htmp);
      gr->GetXaxis()->SetTitle(htmp->GetXaxis()->GetTitle());
      //gr->GetXaxis()->SetTitleOffset(1.2);
      //gr->GetYaxis()->SetTitleOffset(1.5);
      gr->GetYaxis()->SetTitle(htmp->GetYaxis()->GetTitle());
      if(setYRange)  gr->GetYaxis()->SetRangeUser(min, max);
      if(setXRange)  gr->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
      gr->SetMarkerSize(0);
      gr->SetFillStyle(0);
      gr->SetLineColor(gColor[i]);
      if(i==1)
	gr->Draw("ape5");
      else
	gr->Draw("samepe5");
      leg->AddEntry(gr,legName[i],"F");
    }
  h->DrawCopy("sames");
  
  if(drawLeg)
    leg->Draw();
  if(title1.Length()>0)
    {
      TPaveText *t1 = GetTitleText(title1,size);
      t1->Draw();
    }
  
  c->cd(2);
  TString title2="";
  TH1 *h1 = (TH1*)h->Clone(Form("%s_clone",h->GetName()));
  h1->SetMarkerStyle(21);
  h1->SetMarkerColor(1);
  h1->SetLineColor(1);
  for(Int_t ibin=1; ibin<=h1->GetNbinsX(); ibin++)
    {
      h1->SetBinError(ibin,h1->GetBinError(ibin)/h1->GetBinContent(ibin));
      h1->SetBinContent(ibin,1);
    }
  for(Int_t i=1; i<nHistos; i++)
    {
      TH1* htmp = (TH1*)list->At(i);
      if(hTitle2.Length()>0) htmp->SetTitle(hTitle2.Data());
      title2 = htmp->GetTitle();
      htmp->SetTitle("");
      TGraphAsymmErrors *gr = new TGraphAsymmErrors(htmp);
      gr->GetXaxis()->SetTitle(htmp->GetXaxis()->GetTitle());
      gr->GetYaxis()->SetTitle(htmp->GetYaxis()->GetTitle());
      Int_t npoint = gr->GetN();
      Double_t x,y;
      for(Int_t ipoint=0; ipoint<npoint; ipoint++)
	{
	  gr->GetPoint(ipoint,x,y);
	  gr->SetPoint(ipoint,x,1);
	  gr->SetPointEYhigh(ipoint,gr->GetErrorYhigh(ipoint)/y);
	  gr->SetPointEYlow(ipoint,gr->GetErrorYlow(ipoint)/y);
	}
      gr->SetTitle(Form(";p_{T,jet} (GeV/c); Relative uncertainty"));
      if(setYRange1)  gr->GetYaxis()->SetRangeUser(min1, max1);
      if(setXRange)  gr->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
      gr->SetMarkerSize(0);
      gr->SetFillStyle(0);
      gr->SetLineColor(gColor[i]);
      //gr->GetXaxis()->SetTitleOffset(1.2);
      //gr->GetYaxis()->SetTitleOffset(1.5);
      if(i==1)
	gr->Draw("ape5");
      else
	gr->Draw("samepe5");
    }
  h1->DrawCopy("sames");
  if(title2.Length()>0)
    {
      TPaveText *t1 = GetTitleText(title2,size);
      t1->Draw();
    }

  return c;
}


//-----------------------------------------
TCanvas *drawLists(TList **list, Int_t nList, TString cTitle, TString hTitle="", Bool_t setLogy=kTRUE, Bool_t drawLeg=kFALSE, TString *legName=0x0, TString legTitle = "", Float_t llx=0.3, Float_t lhx=0.5, Float_t lly=0.6, Float_t lhy=0.8,Bool_t drawP=kTRUE, const Float_t size = 0.04, const Float_t legSize = 0.04, const Double_t bottomMargin=-1, const Double_t leftMargin = -1, const Double_t rightMargin = -1)
{
  TCanvas *c = 0x0;
  if(nList==3)
    c = new TCanvas(cTitle.Data(),cTitle.Data(),1250,450);
  else if(nList==1)
    c = new TCanvas(cTitle.Data(),cTitle.Data(),800,600);
  else
    c = new TCanvas(cTitle.Data(),cTitle.Data(),1150,700);

  if(nList>1)
    {  
      if(nList==2)
	c->Divide(2,1);
      else if(nList==3)
	c->Divide(3,1);
      else if(nList==9)
	c->Divide(3,3);
      else if(nList%2==0)
	c->Divide(nList/2,2);
      else
	c->Divide((nList+1)/2,2);  
    }

  TLegend *leg = new TLegend(llx,lly,lhx,lhy);
  if(legTitle.Length()>0)
    leg->SetHeader(legTitle.Data());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(legSize);
  
  for(Int_t il = 0; il<nList; il++)
    {
      TList *iList = list[il];
      Int_t nHistos=iList->GetEntries();
      c->cd(il+1);
      if(setLogy) gPad->SetLogy();
      if(bottomMargin>0) SetPadMargin(gPad,bottomMargin,leftMargin, rightMargin);
      else
	{
	  if(nList==4) SetPadMargin(gPad,0.15,0.1,0.05);
	}
      for(Int_t i=0; i<nHistos; i++)
	{
	  TH1 *h = (TH1*)iList->At(i);
	  if(hTitle.Length()>0)
	    h->SetTitle(hTitle.Data());
	  if(drawP)
	    {
	      h->SetMarkerStyle(20);
	      h->SetMarkerSize(0.9);
	    }

	  h->SetMarkerColor(color[i]);
	  h->SetLineColor(color[i]);
	  
	  if(i==0)
	    {
	      TPaveText *t1 = GetTitleText(h->GetTitle(),size);
	      h->SetTitle("");
	      if(drawP)
		h->DrawCopy("P");
	      else
		h->DrawCopy("HIST");
	      t1->Draw();
	    }
	  else
	    {
	      if(drawP)
		h->DrawCopy("samesP");
	      else
		h->DrawCopy("HIST sames");
	    }

	  if(il==0)
	    {
	      if(drawP)
		{
		  if(legName)
		    leg->AddEntry(h,legName[i].Data(),"P");
		  else
		    leg->AddEntry(h,h->GetName(),"P");
		}
	      else
		{
		  if(legName)
		    leg->AddEntry(h,legName[i].Data(),"l");
		  else
		    leg->AddEntry(h,h->GetName(),"l");
		}
	    }
	}
    }
  if(drawLeg)
    {
      c->cd(1);
      leg->Draw();
    }

  return c;
}

//-----------------------------------------
TCanvas *drawGraphLists(TList **list, Int_t nList, TString cTitle, TString hTitle="", Bool_t setLogy=kTRUE, Bool_t drawLeg=kFALSE, TString *legName=0x0, TString legTitle = "", Float_t llx=0.3, Float_t lhx=0.5, Float_t lly=0.6, Float_t lhy=0.8,Bool_t drawP=kTRUE, const Float_t size = 0.04, const Float_t legSize = 0.04, const Double_t bottomMargin=0.1, const Double_t leftMargin = 0.1, const Double_t rightMargin = 0.1)
{
  TCanvas *c = 0x0;
  if(nList==3)
    c = new TCanvas(cTitle.Data(),cTitle.Data(),1250,450);
  else if(nList==1)
    c = new TCanvas(cTitle.Data(),cTitle.Data(),800,600);
  else
    c = new TCanvas(cTitle.Data(),cTitle.Data(),1150,700);

  if(nList>1)
    {  
      if(nList==2)
	c->Divide(2,1);
      else if(nList==3)
	c->Divide(3,1);
      else if(nList==9)
	c->Divide(3,3);
      else if(nList%2==0)
	c->Divide(nList/2,2);
      else
	c->Divide((nList+1)/2,2);  
    }

  TLegend *leg = new TLegend(llx,lly,lhx,lhy);
  if(legTitle.Length()>0)
    leg->SetHeader(legTitle.Data());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(legSize);
  
  for(Int_t il = 0; il<nList; il++)
    {
      TList *iList = list[il];
      if(!iList) continue;
      Int_t nHistos=iList->GetEntries();
      c->cd(il+1);
      gPad->SetLogy((Int_t)setLogy);
      SetPadMargin(gPad,bottomMargin,leftMargin, rightMargin);

      for(Int_t i=0; i<nHistos; i++)
	{
	  TGraph *h = (TGraph*)iList->At(i);
	  const char *histo_title = h->GetHistogram()->GetTitle();
	  if(hTitle.Length()>0)
	    h->SetTitle(hTitle.Data());
	  h->SetMarkerColor(color[i]);
	  h->SetLineColor(color[i]);
	  
	  if(i==0)
	    {
	      TPaveText *t1 = GetTitleText(histo_title,size);
	      h->SetTitle("");
	      h->Draw("AP");
	      t1->Draw();
	    }
	  else
	    {
	      h->Draw("samesP");
	    }

	  if(il==0)
	    {
	      if(legName)
		leg->AddEntry(h,legName[i].Data(),"PLE");
	      else
		leg->AddEntry(h,h->GetName(),"PLE");
	    }
	}
    }

  if(drawLeg)
    {
      c->cd(1);
      leg->Draw();
    }

  return c;
}

//-----------------------------------------
TCanvas *draw1D(TH1 *h, TString hTitle, Bool_t setLog, Bool_t drawP, const Float_t size, const TString drawOpt, const Int_t titleFont,
		const Int_t wh, const Int_t ww)
{
  TCanvas *c = new TCanvas(h->GetName(),h->GetName(),wh,ww);
  if(setLog) gPad->SetLogy();
  if(hTitle.Length()>0) h->SetTitle(hTitle);
  TPaveText *t1 = GetTitleText(h->GetTitle(),size,titleFont);
  h->SetTitle("");
  if(drawOpt.Length()>0)
    h->DrawCopy(drawOpt.Data());
  else
    {
      if(drawP)
	h->DrawCopy("PE");
      else
	h->DrawCopy("HIST");
    }
  t1->Draw();
  return c;
}

//-----------------------------------------
TCanvas *draw2D(TH2 *h, const TString hTitle, const Float_t size, const Bool_t logz, const char *drawOption)
{
  TCanvas *c = new TCanvas(h->GetName(),h->GetName(),800,600);
  if(hTitle.Length()>0)
    h->SetTitle(hTitle);
  if(logz) gPad->SetLogz();
  TPaveText *t1 = GetTitleText(h->GetTitle(),size);
  h->SetTitle("");
  h->DrawCopy(drawOption);
  t1->Draw();
  return c;
}

//-----------------------------------------
TCanvas *draw2DList(TList *list, const Float_t size = 0.04, const Bool_t logy = kFALSE, const TString hTitle = "", const Double_t bottomMargin=-1, const Double_t leftMargin = -1, const Double_t rightMargin = -1, const Bool_t logz = kTRUE, const Bool_t logx = kFALSE, const char *drawOption = "colz")
{
  TH2 *h = (TH2*)list->At(0);
  if(!h) return 0x0;
  Int_t nHisto = list->GetEntries();

  TCanvas *c = 0x0;
  if(nHisto==3)
    c = new TCanvas(h->GetName(),h->GetName(),1150,450);
  else if(nHisto==1)
    c = new TCanvas(h->GetName(),h->GetName(),800,600);
  else if(nHisto==2)
    c = new TCanvas(h->GetName(),h->GetName(),1100,500);
  else
    c = new TCanvas(h->GetName(),h->GetName(),1100,600);

  if(nHisto>1)
    {  
      if(nHisto==2)
	c->Divide(2,1);
      else if(nHisto==3)
	c->Divide(3,1);
      else if(nHisto%2==0)
	c->Divide(nHisto/2,2);
      else
	c->Divide((nHisto+1)/2,2);  
    }

  for(Int_t i=0; i<nHisto; i++)
    {
      TH2 *h = (TH2*)list->At(i);
      // if(hTitle.Length()>0)
      // 	h->SetTitle(hTitle);
      c->cd(i+1);
      if(bottomMargin>0) SetPadMargin(gPad,bottomMargin,leftMargin, rightMargin);
      else
	{
	  if(nHisto==4) SetPadMargin(gPad,0.15,0.1,0.05);
	}
      if(logz) gPad->SetLogz();
      if(logy) gPad->SetLogy();
      if(logx) gPad->SetLogx();
      TPaveText *t1 = GetTitleText(h->GetTitle(),size);
      h->SetTitle("");
      h->DrawCopy(drawOption);
      t1->Draw();
    }
  return c;
}

//-----------------------------------------
TCanvas *drawGraph(TGraph *h, const TString hTitle,Bool_t setLog, const Float_t size, const char *drawOption, const Int_t font)
{
  TCanvas *c = new TCanvas(h->GetName(),h->GetName(),800,600);
  if(setLog)
    gPad->SetLogy();
  if(hTitle.Length()>0)
    h->SetTitle(hTitle);
  TPaveText *t1 = GetTitleText(h->GetHistogram()->GetTitle(),size,font);
  h->SetTitle("");
  h->Draw(drawOption);
  t1->Draw("sames");
  return c;
}

//-----------------------------------------
void SetPadMargin(TVirtualPad *pad, const Double_t bottomMargin, const Double_t leftMargin, const Double_t rightMargin, const Double_t topMargin)
{
  if(!pad) return;
  pad->SetLeftMargin(leftMargin);
  pad->SetBottomMargin(bottomMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetTopMargin(topMargin);
}

//-----------------------------------------
void ScaleHistoTitle(TH1 *h, 
		     const Double_t xTitleSize, const Double_t xTitleOffset, const Double_t xLabelSize,
		     const Double_t yTitleSize, const Double_t yTitleOffset, const Double_t yLabelSize,
		     const Int_t font)
{
  if(!h) return;
  h->GetXaxis()->SetTitleFont(font);
  h->GetXaxis()->SetLabelFont(font);
  h->GetYaxis()->SetTitleFont(font);
  h->GetYaxis()->SetLabelFont(font);

  h->GetXaxis()->SetTitleSize(xTitleSize);
  h->GetXaxis()->SetTitleOffset(xTitleOffset);
  h->GetXaxis()->SetLabelSize(xLabelSize);
  h->GetYaxis()->SetTitleSize(yTitleSize);
  h->GetYaxis()->SetTitleOffset(yTitleOffset);
  h->GetYaxis()->SetLabelSize(yLabelSize);
}

//__________________________________________________________________________
void scaleHisto(TH1 *h, Double_t nEvents=1, Double_t acceptance=1, Bool_t binWidth=kTRUE, Bool_t rebin=kFALSE, Bool_t Sumw2=kTRUE, const Int_t type=0)
{
  Double_t entry = h->GetEntries();
  if(Sumw2)
    h->Sumw2();
  if(rebin)
    h->Rebin(2);
  h->Scale(1./nEvents);
  h->Scale(1./acceptance);
  if(binWidth)
    {
      for(Int_t ibin=1; ibin<h->GetNbinsX()+1; ibin++)
	{
	  h->SetBinContent(ibin, (h->GetBinContent(ibin))/(h->GetBinWidth(ibin)));
	  if(type==0)
	    h->SetBinError(ibin, (h->GetBinError(ibin))/(h->GetBinWidth(ibin)));
	}
    }
  h->SetEntries(entry);
}

//-----------------------------------------
TPaveText *GetTitleText(TString title, const Float_t size, const Int_t font)
{
  TPaveText* t1=new TPaveText(0.3530151,0.8968531,0.6532663,0.9965035,"brNDC");
  //TPaveText* t1=new TPaveText(0.3530151,0.7968531,0.6532663,0.8965035,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,title.Data());
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
}

//-----------------------------------------
TPaveText *GetJpsiPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, const TString runtype = "Run13_pp500", Double_t size = 0.04)
{
  TPaveText* t1 = GetPaveText(xl,xh,yl,yh,size);
  t1->SetFillStyle(1);
  t1->SetFillColor(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(size);
  t1->SetTextFont(42);
  t1->SetTextAlign(11);
  if(runtype=="Run13_pp500")
    {
      t1->AddText("#it{Run13 pp #sqrt{s} = 500 GeV}");
      t1->AddText("#it{di-muon trigger}");
    }
  else if(runtype=="Run14_AuAu200")
    {
      t1->AddText("#it{Run14 AuAu #sqrt{s_{NN}} = 200 GeV}");
      t1->AddText("#it{di-muon trigger}");
    }
  return t1;
}

//-----------------------------------------
TPaveText *GetHJetANPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, const TString system = "PbPb", Double_t size = 0.04, Int_t font = 42)
{
  TPaveText *t1 = GetPaveText(xl,xh,yl,yh,size,font);
  if(system=="PbPb") t1->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, 0-10%");
  else if (system=="Pythia_PL") t1->AddText("PYTHIA particle-level #sqrt{#it{s}} = 2.76 TeV");
  else if (system=="Pythia_DL") t1->AddText("PYTHIA detector-level #sqrt{#it{s}} = 2.76 TeV");
  else t1->AddText("#sqrt{#it{s}_{NN}} = 2.76 TeV");
  t1->AddText(Form("anti-k_{T} charged jets, #it{R} = 0.4"));
  if(system=="PbPb") t1->AddText(Form("40 < #it{p}_{T,jet}^{reco,ch} < 60 GeV/#it{c}"));
  else if (system=="Pythia_PL") t1->AddText(Form("40 < #it{p}_{T,jet}^{part,ch} < 60 GeV/#it{c}"));
  else if (system=="Pythia_DL") t1->AddText(Form("40 < #it{p}_{T,jet}^{det,ch} < 60 GeV/#it{c}"));
  t1->AddText("TT{20,50}-{8,9}");
  t1->SetTextAlign(11);
  return t1;
}

//-----------------------------------------
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size, const Int_t font)
{
  TPaveText* t1=new TPaveText(xl,yl,xh,yh,"brNDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->SetTextSize(size);
  t1->SetTextFont(font);
  return t1;
}



//--------------------------------------------
TLine *GetLine(Double_t xl, Double_t yl, Double_t xh, Double_t yh, Color_t color, Width_t width, Style_t style)
{
  TLine *line = new TLine(xl,yl,xh,yh);
  line->SetLineColor(color);
  line->SetLineWidth(width);
  line->SetLineStyle(style);
  return line;
}

//--------------------------------------------
void offset_x(TGraph* g, Double_t xoff)
{
  Int_t npoints = g->GetN();
  Double_t* x = g->GetX();
  Double_t* y = g->GetY();
  for (Int_t j=0; j<npoints; j++)
    {
      g->SetPoint(j, x[j] + xoff, y[j]);
    }
}


//--------------------------------------------
void offset_x_with_error(TGraphErrors* g, Double_t xoff)
{
  Int_t npoints = g->GetN();
  Double_t* x = g->GetX();
  Double_t* y = g->GetY();
  Double_t* ex = g->GetEX();
  Double_t* ey = g->GetEY();

  for (Int_t j=0; j<npoints; j++)
    {
      g->SetPoint(j, x[j] + xoff, y[j]);
      g->SetPointError(j, ex[j], ey[j]);
    }
} 

//--------------------------------------------
void offset_x_with_asym_error(TGraphAsymmErrors* g, Double_t xoff)
{
  Int_t npoints = g->GetN();
  Double_t* x = g->GetX();
  Double_t* y = g->GetY();
  Double_t* exl = g->GetEXlow();
  Double_t* exh = g->GetEXhigh();
  Double_t* eyl = g->GetEYlow();
  Double_t* eyh = g->GetEYhigh();

  for (Int_t j=0; j<npoints; j++)
    {
      g->SetPoint(j, x[j] + xoff, y[j]);
      g->SetPointError(j, exl[j], exh[j], eyl[j], eyh[j]);
    }
}   

//--------------------------------------------
void printStat(TFile *f)
{
  const char *triggerName[2] = {"MB","EMC"};

  printf("\n===== Event statistics =====\n");
  TH1F *h1 = (TH1F*)f->Get("fhEventStat");
  if(h1)
    {
      printf("EMCal QA: \n");
      printf("MB:  %5.3e -> %5.3e -> %5.3e\n",h1->GetBinContent(1), h1->GetBinContent(2), h1->GetBinContent(3));
      printf("EMC: %5.3e -> %5.3e -> %5.3e\n",h1->GetBinContent(1), h1->GetBinContent(4), h1->GetBinContent(5));
      printf("\n");
    }
  TH1F *h2 = (TH1F*)f->Get("fhJetQAEventStat");
  if(h2)
    {
      printf("Jet QA: \n");
      printf("MB:  %5.3e -> %5.6e -> %5.6e\n",h2->GetBinContent(1), h2->GetBinContent(2), h2->GetBinContent(3));
      printf("EMC: %5.3e -> %5.6e -> %5.6e\n",h2->GetBinContent(1), h2->GetBinContent(4), h2->GetBinContent(5));
      printf("\n");
    }

  h2 = (TH1F*)f->Get("fhJetEventStat_Baseline");
  if(h2)
    {
      printf("Jet finding: \n");
      printf("MB:  %5.3e -> %5.6e -> %5.6e \n",h2->GetBinContent(1), h2->GetBinContent(2), h2->GetBinContent(3));
      printf("EMC: %5.3e -> %5.6e -> %5.6e \n",h2->GetBinContent(1), h2->GetBinContent(4), h2->GetBinContent(5));
      printf("\n");
    }

  printf("=============================\n\n");
}

//================================================
TH1F *DivideTH1ForEff(TH1F *hMatch, TH1F *hBase, const char *effName)
{
  if(!hMatch || !hBase) return 0;
  TGraphAsymmErrors *graph = new TGraphAsymmErrors(hMatch, hBase,"cl=0.683 b(1,1) mode");
  graph->SetName(Form("graph_%s",hMatch->GetName()));
  TH1F *hEff = (TH1F*)hMatch->Clone(effName);
  hEff->Divide(hBase);
  double x,y;
  for(int ipoint=0; ipoint<graph->GetN(); ipoint++)
    {
      graph->GetPoint(ipoint,x,y);
      int bin = hEff->FindFixBin(x);
      double err1 = graph->GetErrorYhigh(ipoint);
      double err2 = graph->GetErrorYlow(ipoint);
      double err = (err1>err2) ? err1 : err2;
      hEff->SetBinError(bin,err);
    }
  delete graph;
  return hEff;
}

//================================================
double CrystalBall(double *x, double *par)
{
  double t = (x[0]-par[1])/par[2];
  if(par[0]<0) t = -t;

  double absAlpha = TMath::Abs(par[0]);
  double cb;
  if(t>-absAlpha)
    {
      cb = par[4]*exp(-0.5*t*t);
    }
  else
    {
      double aa = TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
      double bb = par[3]/absAlpha-absAlpha;
      cb = par[4]*aa/TMath::Power(bb-t,par[3]);
    }

  double f = cb + par[5]*exp(-0.5*TMath::Power((x[0]-par[6])/par[7],2));
  return cb;
}

//================================================
Double_t DoubleCrystalBall(Double_t *x, Double_t *par)
{
  Double_t N = par[0];
  Double_t mu = par[1];
  Double_t s = par[2];
  Double_t n = par[3];
  Double_t alpha2 = par[4];
  Double_t m = par[5];
  Double_t beta = par[6];

  Double_t A = TMath::Power(n/fabs(alpha2), n) * TMath::Exp(-alpha2*alpha2/2.);
  Double_t B = n/fabs(alpha2) - fabs(alpha2);

  Double_t C = TMath::Power(m/fabs(beta), m) * TMath::Exp(-beta*beta/2.);
  Double_t D = m/fabs(beta) - fabs(beta);

  Double_t norm = (x[0]-mu)/s;

  if(norm < -alpha2) {
    return N * A * TMath::Power(B-norm, -n);
  } else if(norm < beta) {
    return N * TMath::Exp(-0.5*norm*norm);
  } else {
    return N * C * TMath::Power(D+norm, -m);
  }
}


//================================================
double CrystalBallPlusPol1(double *x, double *par)
{
  double t = (x[0]-par[1])/par[2];
  if(par[0]<0) t = -t;

  double absAlpha = TMath::Abs(par[0]);
  double cb;
  if(t>-absAlpha)
    {
      cb = par[4]*exp(-0.5*t*t);
    }
  else
    {
      double aa = TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
      double bb = par[3]/absAlpha-absAlpha;
      cb = par[4]*aa/TMath::Power(bb-t,par[3]);
    }

  double f = cb + par[5] + par[6]*x[0]; 
  return f; 
}

//================================================
double CrystalBallPlusPol3(double *x, double *par)
{
  double t = (x[0]-par[1])/par[2];
  if(par[0]<0) t = -t;

  double absAlpha = TMath::Abs(par[0]);
  double cb;
  if(t>-absAlpha)
    {
      cb = par[4]*exp(-0.5*t*t);
    }
  else
    {
      double aa = TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
      double bb = par[3]/absAlpha-absAlpha;
      cb = par[4]*aa/TMath::Power(bb-t,par[3]);
    }

  double f = cb + par[5] + par[6]*x[0] + par[7]*x[0]*x[0] + par[8]*x[0]*x[0]*x[0]; 
  return f; 
}
