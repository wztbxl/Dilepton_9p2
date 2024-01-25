#include "TStyle.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TVirtualFitter.h"

#include <cmath>

#include "defs.h"
#include "drawHistos.C"

using namespace std;
TRandom3 *myRandom;
const double pi = TMath::Pi();
TH1F *hgCosTheta;
TH1F *hgPhi;
Int_t npfits;
const TString parName[4] = {"norm1","lambdaTheta","lambdaPhi","norm2"};
double iniPar[4] = { 1, 0.15, -0.05, 1 };
double iniVlow[4] = {0, 0, 0, 0};
double iniVhigh[4] = {0, 0, 0, 0};
double stp[4] = {0.01, 0.001, 0.001, 0.01}; 


//================================================
void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{

    TAxis *xaxis1  = hgCosTheta->GetXaxis();
    TAxis *xaxis2  = hgPhi->GetXaxis();

    int nbinX1 = hgCosTheta->GetNbinsX();
    int nbinX2 = hgPhi->GetNbinsX();

    double chi2 = 0; 
    double x;
    double tmp;
    double nstep = 20;
    npfits = 0;
    for (int ix = 1; ix <= nbinX1; ++ix) {
        x = xaxis1->GetBinCenter(ix);
        double xmin = xaxis1->GetBinLowEdge(ix);
        double xmax = xaxis1->GetBinUpEdge(ix);
        double size = (xmax-xmin)/nstep;
        if ( hgCosTheta->GetBinError(ix) > 0 ) {
            double data = hgCosTheta->GetBinContent(ix);
            double fit  = 0;
            for(int istep=0; istep<nstep; istep++)
            {
                double xtmp = xmin + size * (istep+0.5);
                fit += (p[0]*(3*(1+p[1]*xtmp*xtmp))/(2*(p[1]+3))) * size;
            }
            fit /= (xmax-xmin);
            tmp = (data-fit)/hgCosTheta->GetBinError(ix);
            //tmp = ( hgCosTheta->GetBinContent(ix) - p[0]*(3*(1+p[1]*x*x))/(2*(p[1]+3)))/hgCosTheta->GetBinError(ix);
            chi2 += tmp*tmp; 
            npfits++;
        }
    }
    for (int ix = 1; ix <= nbinX2; ++ix) {
        x = xaxis2->GetBinCenter(ix);
        double xmin = xaxis2->GetBinLowEdge(ix);
        double xmax = xaxis2->GetBinUpEdge(ix);
        double size = (xmax-xmin)/nstep;
        if ( hgPhi->GetBinError(ix) > 0 ) {
            double data = hgPhi->GetBinContent(ix);
            double fit  = 0;
            for(int istep=0; istep<nstep; istep++)
            {
                double xtmp = xmin + size * (istep+0.5);
                fit += (p[3]*((1+2*p[2]*cos(2*xtmp)/(3+p[1]))/(pi))) * size;
            }
            fit /= (xmax-xmin);
            tmp = (data-fit)/hgPhi->GetBinError(ix);
            //tmp = (hgPhi->GetBinContent(ix) - p[3]*((1+2*p[2]*cos(2*x)/(3+p[1]))/(pi)))/hgPhi->GetBinError(ix);
            chi2 += tmp*tmp; 
            npfits++;
        }
    }
    fval = chi2; 
}

//================================================
void ana_JpsiPol()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatY(0.9);                
    gStyle->SetStatX(0.9);  
    gStyle->SetStatW(0.2);                
    gStyle->SetStatH(0.2);

    myRandom = new TRandom3();
    myRandom->SetSeed(0);

    time_t t_start, t_end;
    t_start = time(NULL) ;

    fit_test_EffCorr();

    t_end = time(NULL) ;
    printf("time: %.0f s\n", difftime(t_end,t_start)) ;
}


//================================================
void fit_test_EffCorr(const int savePlot = 1, const int saveHisto = 1, const int savePDF = 1)
{
    TFile *feff[nthetaphi];
    TH2F *hJpsiPolVsPt[nthetaphi][2][nExpr+1][2];
    TH1F *hPolDis[nthetaphi][nPtBins][2][nExpr+1][2];
    TH1F *hEffPol[nthetaphi][nPtBins][2][nExpr+1];

    TCanvas *c[2];
    for(int it=0; it<nthetaphi; it++)
    {
        feff[it] = TFile::Open(Form("output/Iter%d.Loop1.Eff_%3.3f_%3.3f_%1.1f.root",gIter,ptheta[it],pphi[it],pthephi),"read");
        for(int i=0; i<2; i++)
        {
            for(int j=0; j<nExpr; j++) //previous nExpr+1
            {
                if((gIter==10 || gIter==1) && j>0) continue;
                for(int t=0; t<2; t++)
                {
                    if(gIter==10 || gIter==1) hJpsiPolVsPt[it][t][j][i] = (TH2F*)feff[it]->Get(Form("mh%s%sVsPt",hname[i],typeName[t]));
                    else  hJpsiPolVsPt[it][t][j][i] = (TH2F*)feff[it]->Get(Form("mh%s%sVsPt_Expr%d",hname[i],typeName[t],j));

                    for(int p=0; p<nPtBins; p++)
                    {
                        hPolDis[it][p][t][j][i] = (TH1F*)hJpsiPolVsPt[it][t][j][i]->ProjectionY(Form("hPolDis_Pt%d_%s_%s_%3.3f_%3.3f_Expr%d",p,typeName[t],hname[i],ptheta[it],pphi[it],j), p+1, p+1);
                    }
                }
            }
        }

        for(int j=0; j<nExpr; j++) 
        {
            if((gIter==10 || gIter==1) && j>0) continue;
            for(int t=0; t<2; t++)
            {
                for(int p=0; p<nPtBins; p++)
                {
                    hEffPol[it][p][t][j] = (TH1F*)hPolDis[it][p][t][j][1]->Clone(Form("Iter%d_hEffPol_Pt%d_%s_%3.3f_%3.3f_Expr%d",gIter,p,typeName[t],ptheta[it],pphi[it],j));
                    hEffPol[it][p][t][j]->Divide(hPolDis[it][p][t][j][0]);
                }
            }
        }
    }

    // Get Pseudo data
    TFile *fdata[nthetaphi];
    TH2F *hRcPolDisVsPt[nthetaphi][2][nExpr];
    TH1F *hRcPolDis[nthetaphi][nPtBins][2][nExpr];
    TH2F *hRcPolDisRelErr[nPtBins][2];
    for(int it=0; it<nthetaphi; it++)
    {
        fdata[it] = TFile::Open(Form("Output/ToyMC_Data_Iter0_%3.3f_%3.3f_%1.1f.root",ptheta[it],pphi[it],pthephi),"read");// for 500
        for(int i=0; i<nExpr; i++)
        {
            for(int t=0; t<2; t++)
            {
                hRcPolDisVsPt[it][t][i] = (TH2F*)fdata[it]->Get(Form("mhRc%sVsPt_Loop%d",typeName[t],i+1));
                hRcPolDisVsPt[it][t][i]->SetName(Form("%s_%3.3f_%3.3f",hRcPolDisVsPt[it][t][i]->GetName(),ptheta[it],pphi[it]));
                for(int p=0; p<nPtBins; p++)
                {
                    hRcPolDis[it][p][t][i] = (TH1F*)hRcPolDisVsPt[it][t][i]->ProjectionY(Form("Iter%d_hRcPolDis_Pt%d_%s_%3.3f_%3.3f_Expr%d",gIter,p,typeName[t],ptheta[it],pphi[it],i+1),p+1,p+1);
                    if(ptheta[it]==0.0 && pphi[it]==0.0)
                    {
                        if(i==0)
                        {
                            TAxis *axis = hRcPolDis[it][p][t][i]->GetXaxis();
                            int nbinsx = axis->GetNbins();
                            double xmin = axis->GetXmin();
                            double xmax = axis->GetXmax();
                            hRcPolDisRelErr[p][t] = new TH2F(Form("Iter%d_hRcPolDisRelErr_Pt%d_%s_%3.3f_%3.3f",gIter,p,typeName[t],ptheta[it],pphi[it]),"",nbinsx,xmin,xmax,nExpr,0,nExpr);
                        }
                        int nbins = hRcPolDis[it][p][t][i]->GetNbinsX();
                        for(int ibin=1; ibin<=nbins; ibin++)
                        {
                            if(hRcPolDis[it][p][t][i]->GetBinContent(ibin)>0)
                            {
                                double rel_err = hRcPolDis[it][p][t][i]->GetBinError(ibin)/hRcPolDis[it][p][t][i]->GetBinContent(ibin);
                                hRcPolDisRelErr[p][t]->SetBinContent(ibin, i+1, rel_err);
                                hRcPolDisRelErr[p][t]->SetBinError(ibin, i+1, 0);
                            }
                        }
                    }
                    if(gIter==10 || gIter==1) hRcPolDis[it][p][t][i]->Divide(hEffPol[it][p][t][0]);
                    else hRcPolDis[it][p][t][i]->Divide(hEffPol[it][p][t][i]);
                    hRcPolDis[it][p][t][i]->Scale(1./hRcPolDis[it][p][t][i]->Integral());
                }
            }
        }
    }

    // fit cos(theta) && phi simultaneously
    TF1 *fitJpsiPol[nthetaphi][nPtBins][2][nExpr];
    TH1F *hFitJpsiPol[nthetaphi][nPtBins][3];
    double myPar[nthetaphi][nPtBins][3][nExpr];
    double myParErr[nthetaphi][nPtBins][3][nExpr];
    double myMatrix[nthetaphi][nPtBins][nExpr];
    for(int it=0; it<nthetaphi; it++)
    {
        for(int p=0; p<nPtBins; p++)
        {
            for(int t=0; t<3; t++)
            {
                hFitJpsiPol[it][p][t] = new TH1F(Form("Iter%d_hFitJpsiPol_%s_Pt%d_%3.3f_%3.3f",gIter,typeName[t],p,ptheta[it],pphi[it]),"",nExpr,0,nExpr);
            }
        }
    }

    if(nExpr<50)
    {
        const int nExprSave = nExpr;
    }
    else
    {
        const int nExprSave = 5;
    }

    TCanvas *cFit[nthetaphi][nExprSave];
    for(int it=0; it<nthetaphi; it++)
    {
        //if(!(ptheta[it]==-1.0 && pphi[it]==0.0)) continue;
        for(int i=0; i<nExpr; i++)
        {
            if(savePDF && i<nExprSave)
            {
                cFit[it][i] = new TCanvas(Form("FitCosTheta_Expr%d_%3.3f_%3.3f",i,ptheta[it],pphi[it]), Form("FitCosTheta_Expr%d_%3.3f_%3.3f",i,ptheta[it],pphi[it]), 1200, 700);
                cFit[it][i]->Divide(4,2);
            }
            for(int p=0; p<nPtBins; p++)
            {
                int nbins_theta = hRcPolDis[it][p][0][i]->GetNbinsX();
                for(int ibin=1; ibin<=nbins_theta; ibin++)
                {
                    if( (p==0) && (ibin<=3 || ibin>=8))
                    {
                        hRcPolDis[it][p][0][i]->SetBinContent(ibin, 0);
                        hRcPolDis[it][p][0][i]->SetBinError(ibin, 0);
                    }
                    else if( (p==1 || p==2) && (ibin<=2 || ibin>=9))
                    {
                        hRcPolDis[it][p][0][i]->SetBinContent(ibin, 0);
                        hRcPolDis[it][p][0][i]->SetBinError(ibin, 0);
                    }
                    else if( p==3  && (ibin<2 || ibin>9))
                    {
                        hRcPolDis[it][p][0][i]->SetBinContent(ibin, 0);
                        hRcPolDis[it][p][0][i]->SetBinError(ibin, 0);
                    }
                }

                int nbins_phi = hRcPolDis[it][p][1][i]->GetNbinsX();
                for(int ibin=1; ibin<=nbins_phi; ibin++)
                {
                    if( (p==1 || p==2 || p==3) && (ibin<=4 || ibin>=12))
                    {
                        hRcPolDis[it][p][1][i]->SetBinContent(ibin, 0);
                        hRcPolDis[it][p][1][i]->SetBinError(ibin, 0);
                    }
                }

                hgCosTheta = (TH1F*)hRcPolDis[it][p][0][i]->Clone(Form("%s_fit",hRcPolDis[it][p][0][i]->GetName()));
                hgPhi = (TH1F*)hRcPolDis[it][p][1][i]->Clone(Form("%s_fit",hRcPolDis[it][p][1][i]->GetName()));

                //The default minimizer is Minuit
                TVirtualFitter::SetDefaultFitter("Minuit");
                TVirtualFitter * minuit = TVirtualFitter::Fitter(0,4);
                if(p==0 && ptheta[it]==-1.0 && pphi[it]==0.0)
                {
                    iniPar[1] = ptheta[it];
                    iniPar[2] = pphi[it];
                }
                for (int ipar=0; ipar<4; ++ipar)
                {
                    minuit->SetParameter(ipar, parName[ipar].Data(), iniPar[ipar], stp[ipar], iniVlow[ipar], iniVhigh[ipar]);
                }
                minuit->SetFCN(myFcn);
                double arglist[10];
                arglist[0] = 0;
                // set print level
                minuit->ExecuteCommand("SET PRINT",arglist,3);
                // minimize
                arglist[0] = 5000; // number of function calls
                arglist[1] = 0.01; // tolerance
                minuit->ExecuteCommand("MIGRAD",arglist,2);
                double chi2, edm, errdef; 
                int nvpar, nparx;
                minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
                int ndf = npfits-nvpar;
                //--- get error matrix ---//
                Double_t matrix[4][4];
                gMinuit->mnemat(&matrix[0][0],4);
                myMatrix[it][p][i] = matrix[1][2];
                //cout<<matrix[1][2]<<endl;
                //cout<<matrix[2][1]<<endl;

                // get result
                fitJpsiPol[it][p][0][i] = new TF1(Form("FitTheta_Pt%d_%3.3f_%3.3f_Expr%d",p,ptheta[it],pphi[it],i),"[0]*(3*(1+[1]*x*x))/(2*([1]+3))",-1,1);
                for(int ipar=0; ipar<2; ipar++)
                {
                    fitJpsiPol[it][p][0][i]->SetParameter(ipar, minuit->GetParameter(ipar));
                    fitJpsiPol[it][p][0][i]->SetParError(ipar, minuit->GetParError(ipar));
                }	
                fitJpsiPol[it][p][1][i] = new TF1(Form("FitPhi_Pt%d_%1.1f_Expr%d",p,ptheta[it],pphi[it],i),"[0]*((1+2*[1]*cos(2*x)/(3+[2]))/(pi))",0, pi);
                for(int ipar=0; ipar<3; ipar++)
                {
                    int index = ipar;
                    if(ipar==0) index = 3;
                    if(ipar==1) index = 2;
                    if(ipar==2) index = 1;
                    fitJpsiPol[it][p][1][i]->SetParameter(ipar, minuit->GetParameter(index));
                    fitJpsiPol[it][p][1][i]->SetParError(ipar, minuit->GetParError(index));
                }

                for(int t=0; t<2; t++)
                {
                    hFitJpsiPol[it][p][t]->SetBinContent(i+1, fitJpsiPol[it][p][t][i]->GetParameter(1));
                    hFitJpsiPol[it][p][t]->SetBinError(i+1, fitJpsiPol[it][p][t][i]->GetParError(1));
                    myPar[it][p][t][i] = fitJpsiPol[it][p][t][i]->GetParameter(1);
                    myParErr[it][p][t][i] = fitJpsiPol[it][p][t][i]->GetParError(1);
                    if(savePDF && i < nExprSave)
                    {
                        cFit[it][i]->cd(t*nPtBins + p + 1);
                        ScaleHistoTitle(hRcPolDis[it][p][t][i],0.05,0.9,0.04,0.05,1,0.04,62);
                        hRcPolDis[it][p][t][i]->SetTitle(Form(";%s;",typeTitle[t]));
                        hRcPolDis[it][p][t][i]->SetMarkerStyle(21);
                        hRcPolDis[it][p][t][i]->SetMaximum(1.5*hRcPolDis[it][p][t][i]->GetMaximum());
                        hRcPolDis[it][p][t][i]->SetMinimum(0.8*hRcPolDis[it][p][t][i]->GetMinimum());
                        hRcPolDis[it][p][t][i]->Draw();
                        fitJpsiPol[it][p][t][i]->SetLineStyle(2);
                        fitJpsiPol[it][p][t][i]->SetLineColor(4);
                        fitJpsiPol[it][p][t][i]->Draw("sames");
                        TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T,J/#psi} < %1.0f GeV/c",lowPtBins[p], highPtBins[p]),0.055);
                        t1->Draw();
                        t1 = GetPaveText(0.15,0.3,0.7,0.85,0.05);
                        t1->SetTextAlign(11);
                        t1->SetTextColor(kBlue);
                        if(t==0)
                        {
                            t1->AddText(Form("#lambda_{#theta}^{input } = %3.3f",ptheta[it]));
                            t1->AddText(Form("#lambda_{#theta}^{output} = %4.2f #pm %4.2f",fitJpsiPol[it][p][t][i]->GetParameter(1),fitJpsiPol[it][p][t][i]->GetParError(1)));
                        }
                        else
                        {
                            t1->AddText(Form("#lambda_{#varphi}^{input } = %3.3f",pphi[it]));
                            t1->AddText(Form("#lambda_{#varphi}^{output} = %4.2f #pm %4.2f",fitJpsiPol[it][p][t][i]->GetParameter(1),fitJpsiPol[it][p][t][i]->GetParError(1)));
                        }
                        t1->Draw();
                    }
                }
                calInvariant( myPar[it][p][0][i], myParErr[it][p][0][i], myPar[it][p][1][i], myParErr[it][p][1][i], myMatrix[it][p][i], myPar[it][p][2][i], myParErr[it][p][2][i]);
                hFitJpsiPol[it][p][2]->SetBinContent(i+1, myPar[it][p][2][i]);
                hFitJpsiPol[it][p][2]->SetBinError(i+1, myParErr[it][p][2][i]);
            }
        }
    }

    if(saveHisto)
    {
        TFile *fout1 = TFile::Open(Form("Resultfiles/%s.FitJpsiPol.root",run_type),"recreate");//update
        for(int it=0; it<nthetaphi; it++)
        {
            //if(!(ptheta[it]==-1.0 && pphi[it]==0.0)) continue;
            for(int t=0; t<2; t++)
            {
                for(int p=0; p<nPtBins; p++)
                {
                    hEffPol[it][p][t][0]->Write("", TObject::kOverwrite);
                }
            }
        }

        for(int it=0; it<nthetaphi; it++)
        {
            for(int t=0; t<2; t++)
            {
                for(int p=0; p<nPtBins; p++)
                {
                    if(ptheta[it]==0.0 && pphi[it]==0.0)
                    {
                        hRcPolDisRelErr[p][t]->Write("", TObject::kOverwrite);
                    }
                }
            }
        }
        fout1->Close();


        TFile *fout = TFile::Open(Form("Resultfiles/%s.JpsiPolPar.root",run_type),"update");//update
        for(int it=0; it<nthetaphi; it++)
        {
            //if(!(ptheta[it]==-1.0 && pphi[it]==0.0)) continue;
            for(int p=0; p<nPtBins; p++)
            {
                for(int t=0; t<3; t++)
                {
                    hFitJpsiPol[it][p][t]->Write("", TObject::kOverwrite);
                }
            }
        }
        fout->Close();
    }

    if(savePDF)
    {
        for(int it=0; it<nthetaphi; it++)
        {
            //if(!(ptheta[it]==-1.0 && pphi[it]==0.0)) continue;
            TString outPdfName = Form("PDF/FitJpsiPolDis_Iter%d_%3.3f_%3.3f.pdf",gIter,ptheta[it],pphi[it],pphi);
            for(int i=0; i<nExprSave; i++)
            {
                if(nExprSave!=1){
                    if(i==0) cFit[it][i]->Print(Form("%s(",outPdfName.Data()),Form("Title: Expr%d",i+1));
                    else if(i<nExprSave-1) cFit[it][i]->Print(Form("%s",outPdfName.Data()),Form("Title: Expr%d",i+1));
                    else cFit[it][i]->Print(Form("%s)",outPdfName.Data()),Form("Title: Expr%d",i+1));
                }
                else if(nExprSave==1)  cFit[it][i]->SaveAs(Form("%s",outPdfName.Data()));
                //else {
                //    cFit[it][i]->SaveAs(Form("PDF/%s.pdf", outPdfName.Data()));
                //    cout<<outPdfName.Data()<<endl;
                //}
            }
        }
    }

    for(int it=0; it<nthetaphi; it++)
    {
        //if(!(ptheta[it]==-1.0 && pphi[it]==0.0)) continue;
        feff[it]->Close();
        fdata[it]->Close();
    }
}

