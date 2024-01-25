const char* typeName[3] = {"Theta", "Phi", "inv"};
const char* typeTitle[3] = {"cos#theta", "#varphi", "inv"};
const char* polTitle[3] = {"#theta", "#varphi", "inv"};

const int nPtBins = 4;
const double lowPtBins[nPtBins] = {0,1,2,4};
const double highPtBins[nPtBins] = {1,2,4,10};

// Collins-Soper frame //
const int nthetaphi = 4;
double pPar[3][nthetaphi] = {{-0.26, 0.480, 0.630, -0.047}, {-0.034, 0.008, -0.124, -0.3}, {-0.351, 0.463, 0.237, -0.715}};//theta; phi; inv vs. pT bins
double pParErr[3][nthetaphi] = {{0.708, 0.581, 0.362, 0.281}, {0.058, 0.244, 0.253, 0.323}, {0.679, 1.03, 0.784, 0.593}}; // stat. err of theta; phi; inv vs. pT bins

//const double gxbins[nthetaphi] = {0.65, 1.47, 2.68, 4.92};//preliminary
const double gxbins[nPtBins] = {0.62, 1.46, 2.69, 4.92}; //published
const double gxbinsErr[nthetaphi] = {0};

int drawIter(int nIteration = 8, int savePDF=1, int saveHisto=1)
{
    gStyle->SetOptStat("erm");
    gStyle->SetOptFit(111);
    if(nIteration>=8){
        gStyle->SetStatY(0.9);
        gStyle->SetStatX(0.9);
    }
    
    TFile *f = TFile::Open(Form("./Resultfiles/Run15_pp200.JpsiPolPar.root"));
    
    TCanvas *c[3];
    for(int i=0; i<3; i++){
        c[i] = new TCanvas(Form("c_%d",i), Form("c_%s",typeName[i]), 1300, 700);
        c[i]->Divide(4,2);
    }
    
    TH1D *hJpsiPolParAll[nPtBins][3];
    TH1D *hJpsiPolPar[nPtBins][3];
    TH1D *hJpsiPolParErr[nPtBins][3];
    
    for(int t=0; t<3; t++){
        for(int p=0; p<nPtBins; p++){
            hJpsiPolParAll[p][t] = (TH1D*)f->Get(Form("Iter%d_hFitJpsiPol_%s_Pt%d_%3.3f_%3.3f",nIteration, typeName[t], p, pPar[0][p], pPar[1][p]));
            int total = hJpsiPolParAll[p][t]->GetNbinsX();
            if(t<2){
                hJpsiPolPar[p][t] = new TH1D(Form("Iter%d_pt%d_%s_mean", nIteration, p, typeName[t]),Form("; #lambda_{%s}; ", polTitle[t]), 200, -2, 2);
                hJpsiPolParErr[p][t] = new TH1D(Form("Iter%d_pt%d_%s_Err", nIteration, p, typeName[t]),Form("; Err(#lambda_{%s}); ", polTitle[t]), 100, 0, 1);
            }
            else if(t==2) {
                if(p==0) hJpsiPolPar[p][t] = new TH1D(Form("Iter%d_pt%d_%s_mean", nIteration, p, typeName[t]),Form("; #lambda_{%s}; ", polTitle[t]), 400, -4, 6);
                else hJpsiPolPar[p][t] = new TH1D(Form("Iter%d_pt%d_%s_mean", nIteration, p, typeName[t]),Form("; #lambda_{%s}; ", polTitle[t]), 200, -4, 6);
                
                hJpsiPolParErr[p][t] = new TH1D(Form("Iter%d_pt%d_%s_Err", nIteration, p, typeName[t]),Form("; Err(#lambda_{%s}); ", polTitle[t]), 60, 0, 3);
            }
            for(int n=1; n<=total; n++){
                double binContent = hJpsiPolParAll[p][t]->GetBinContent(n);
                double binErr = hJpsiPolParAll[p][t]->GetBinError(n);
                hJpsiPolPar[p][t]->Fill(binContent);
                hJpsiPolParErr[p][t]->Fill(binErr);
            }
            hJpsiPolPar[p][t]->Rebin(4);
            if(p!=0&&t==1)hJpsiPolParErr[p][t]->Rebin(2);
            if(p==0&&t==1){
                hJpsiPolPar[p][t]->GetXaxis()->SetRangeUser(-0.5, 0.5);
                hJpsiPolParErr[p][t]->GetXaxis()->SetRangeUser(0,0.2);
            }
        }
    }
    
    gStyle->SetOptStat("erm");
    // fit combined statistics //
    TF1 *fgPar[nPtBins][3];
    TF1 *fgParErr[nPtBins][3];
    for(int t=0; t<3; t++){
        for(int p=0; p<nPtBins; p++){
            c[t]->cd(1+p);
            ScaleHistoTitle(hJpsiPolPar[p][t],0.05,0.9,0.04,0.05,1,0.04,22);
            hJpsiPolPar[p][t]->SetTitle(Form("#lambda_{%s} distribution (%1.0f < p_{T} < %1.0f GeV/c)", polTitle[t], lowPtBins[p], highPtBins[p]));
            hJpsiPolPar[p][t]->SetMaximum(1.5*hJpsiPolPar[p][t]->GetMaximum());
            hJpsiPolPar[p][t]->SetMinimum(0.8*hJpsiPolPar[p][t]->GetMinimum());
            hJpsiPolPar[p][t]->Draw();
            fgPar[p][t] = new TF1(Form("fgParPt%d%s", p, typeName[t]), "gaus", -1, 1);
            fgPar[p][t]->SetLineStyle(2);
            fgPar[p][t]->SetLineColor(4);
            double mean = hJpsiPolPar[p][t]->GetMean();
            double rms = hJpsiPolPar[p][t]->GetRMS();
            if(t==1 && p==3)fgPar[p][t]->SetRange(mean-1.*rms, mean+2.5*rms);
            else fgPar[p][t]->SetRange(mean-2.5*rms, mean+2.5*rms);
            if(t!=2)hJpsiPolPar[p][t]->Fit(fgPar[p][t],"R");
            TPaveText *t1 = GetTitleText(Form("%1.0f < p_{T,J/#psi} < %1.0f GeV/c",lowPtBins[p], highPtBins[p]),0.055);
            t1 = GetPaveText(0.15,0.3,0.7,0.85,0.05);
            t1->SetTextAlign(11);
            t1->SetTextColor(kBlue);
            t1->AddText(Form("#lambda_{%s}^{input } = %3.2f", polTitle[t], pPar[t][p]));
            if(t!=2) t1->AddText(Form("#lambda_{%s}^{output} = %4.2f #pm %4.2f", polTitle[t],fgPar[p][t]->GetParameter(1),fgPar[p][t]->GetParError(1)));
            else t1->AddText(Form("#lambda_{%s}^{output} = %4.2f #pm %4.2f", polTitle[t], hJpsiPolPar[p][t]->GetMean(), hJpsiPolPar[p][t]->GetMeanError()));
            t1->Draw();
            
            c[t]->cd(p+5);
            ScaleHistoTitle(hJpsiPolParErr[p][t],0.05,0.9,0.04,0.05,1,0.04,22);
            hJpsiPolParErr[p][t]->SetTitle(Form("#lambda_{%s} error (%1.0f < p_{T} < %1.0f GeV/c)", polTitle[t], lowPtBins[p], highPtBins[p]));
            hJpsiPolParErr[p][t]->SetMaximum(1.5*hJpsiPolParErr[p][t]->GetMaximum());
            hJpsiPolParErr[p][t]->SetMinimum(0.8*hJpsiPolParErr[p][t]->GetMinimum());
            hJpsiPolParErr[p][t]->Draw();
            if(nIteration>=8) drawLine(pParErr[t][p],0, pParErr[t][p], 0.5*hJpsiPolParErr[p][t]->GetMaximum(), 2, 1, 2);
            fgParErr[p][t] = new TF1(Form("fgParErrPt%d%s", p, typeName[t]), "gaus", -1, 1);
            fgParErr[p][t]->SetLineStyle(2);
            fgParErr[p][t]->SetLineColor(4);
            mean = hJpsiPolParErr[p][t]->GetMean();
            rms = hJpsiPolParErr[p][t]->GetRMS();
            if((t==1 || t==2) && p>0) fgParErr[p][t]->SetRange(mean-2.5*rms, mean+1.*rms);
            else fgParErr[p][t]->SetRange(mean-2.5*rms, mean+2.5*rms);
            if(nIteration<8){
                if(t!=2)hJpsiPolParErr[p][t]->Fit(fgParErr[p][t],"R");
            TPaveText *t2 = GetTitleText(Form("%1.0f < p_{T,J/#psi} < %1.0f GeV/c",lowPtBins[p], highPtBins[p]),0.055);
            //t2->Draw();
            t2 = GetPaveText(0.15,0.3,0.7,0.85,0.05);
            t2->SetTextAlign(11);
            t2->SetTextColor(kBlue);
            if(t!=2){
                t2->AddText(Form("#sigma(#lambda_{%s}) = %3.2f", polTitle[t], fgPar[p][t]->GetParameter(2)));
                t2->AddText(Form("<Err(#lambda_{%s})> = %4.2f #pm %4.3f", polTitle[t],fgParErr[p][t]->GetParameter(1),fgParErr[p][t]->GetParError(1)));
            }
            else {
                t2->AddText(Form("RMS(#lambda_{%s}) = %3.2f", polTitle[t], hJpsiPolPar[p][t]->GetRMS()));
                t2->AddText(Form("<Err(#lambda_{%s})> = %4.2f #pm %4.3f", polTitle[t],hJpsiPolParErr[p][t]->GetMean(), hJpsiPolParErr[p][t]->GetMeanError()));
            }
                t2->Draw();}
        }
        
        if(!savePDF) continue;
        //c[t]->SaveAs(Form("Toy_getMean_CS_%s.pdf",typeName[t]));
        if(t==0)c[t]->Print(Form("getMean_CS_iter%d.pdf(", nIteration),Form("Title: %s",typeName[t]));
        else if(t==2) c[t]->Print(Form("getMean_CS_iter%d.pdf)", nIteration),Form("Title: %s",typeName[t]));
        else c[t]->Print(Form("getMean_CS_iter%d.pdf", nIteration),Form("Title: %s",typeName[t]));
    }
    
    double mPar[3][nPtBins];
    double mParErr[3][nPtBins];
    
    for(int p=0; p<nPtBins; p++){
        for(int t=0; t<2; t++){
            mPar[t][p] = fgPar[p][t]->GetParameter(1);
            mParErr[t][p] = fgParErr[p][t]->GetParameter(1);
        }
        mPar[2][p] = hJpsiPolPar[p][2]->GetMean();
        mParErr[2][p] = hJpsiPolParErr[p][2]->GetMean();
    }

    //----
    TGraph *gPolParCom[3][4];
       TGraph *gPolParComErr[3][4];
    double x[1], y[1];
       for(int p=0; p<nPtBins; p++){
           for(int t=0; t<3; t++){
               x[0] = pPar[t][p];
               y[0] = mPar[t][p];
               gPolParCom[t][p] = new TGraph(1, x, y);
               gPolParCom[t][p]->SetName(Form("gPolParCom%sPt%d", typeName[t], p));
               SetGraphMarker(gPolParCom[t][p], 29,1);
               if(t!=2) {
                   x[0] = fgPar[p][t]->GetParameter(2);
                   y[0] = mParErr[t][p];
                   gPolParComErr[t][p] = new TGraph(1, x, y);
                   gPolParComErr[t][p]->SetName(Form("gPolParErr%sPt%d", typeName[t], p));
               }
               else {
                   x[0] = hJpsiPolPar[p][t]->GetRMS();
                   y[0] = mParErr[t][p];
                   gPolParComErr[t][p] = new TGraph(1, x, y);
                   gPolParComErr[t][p]->SetName(Form("gPolParErr%sPt%d", typeName[t], p));
               }
               SetGraphMarker(gPolParComErr[t][p], 30,1);
           }}
    
       gStyle->SetOptStat(0);
       TCanvas *c2[2];
       for(int i=0; i<2; i++){
           c2[i] = new TCanvas(Form("c2_%d", i), Form("c2_%d", i), 800, 600);
           c2[i]->Divide(4,3);
       }
    
       for(int p=0; p<nPtBins; p++){
       for(int t=0; t<3; t++){
           c2[0]->cd(t*4+p+1);
           gPolParCom[t][p]->Draw();
           c2[1]->cd(t*4+p+1);
           gPolParComErr[t][p]->Draw();
       }}
    
       if(nIteration==8&&saveHisto){
           TFile *fout2 = new TFile(Form("CS_compare_%d.root", nIteration),"recreate");
           fout2->cd();
           for(int p=0; p<nPtBins; p++){
               for(int t=0; t<3; t++){
                   gPolParCom[t][p]->Write();
                   gPolParComErr[t][p]->Write();
               }}
       }

    //------- diff -------
    double diff[3][nPtBins];
    double diff_stat[3][nPtBins];
    double diff_err[3][nPtBins];
    double diff_err_err[3][nPtBins];
    
    for(int p=0; p<nPtBins; p++){
        for(int t=0; t<2; t++){
            diff[t][p] = mPar[t][p] - pPar[t][p];
            diff_stat[t][p] = fgPar[p][t]->GetParError(1);
            diff_err[t][p] = mParErr[t][p] - fgPar[p][t]->GetParameter(2);
            double a = fgParErr[p][t]->GetParError(1);
            double b = fgPar[p][t]->GetParError(2);
            diff_err_err[t][p] = sqrt(a*a + b*b);
        }
        for(int t=2; t<3; t++){
            diff[t][p] = mPar[t][p] - pPar[t][p];
            diff_stat[t][p] = hJpsiPolPar[p][t]->GetMeanError();
            diff_err[t][p] = mParErr[t][p] - hJpsiPolPar[p][t]->GetRMS();
            double a = hJpsiPolParErr[p][t]->GetMeanError();
            double b = hJpsiPolPar[p][t]->GetRMSError();
            diff_err_err[t][p] = sqrt(a*a + b*b);
        }
    }
    
    TCanvas *c1 = new TCanvas("c1", "c1", 1100, 300);
    c1->Divide(3,1);
    
    gStyle->SetOptStat(0);
    TGraphErrors *gPolParDiff[3];
    TGraphErrors *gPolParErrDiff[3];
        for(int t=0; t<3; t++){
            gPolParDiff[t] = new TGraphErrors(4, gxbins, diff[t], gxbinsErr, diff_stat[t]);
            gPolParDiff[t]->SetName(Form("gPolParDiff%s",typeName[t]));
            gPolParDiff[t]->SetMarkerStyle(21);
            gPolParDiff[t]->SetMarkerColor(1);
            gPolParDiff[t]->SetLineColor(1);
            gPolParErrDiff[t] = new TGraphErrors(4, gxbins, diff_err[t], gxbinsErr, diff_err_err[t]);
            gPolParErrDiff[t]->SetName(Form("gPolParErrDiff%s",typeName[t]));
            gPolParErrDiff[t]->SetMarkerStyle(25);
            gPolParErrDiff[t]->SetMarkerColor(2);
            gPolParErrDiff[t]->SetLineColor(2);
        }
    
    // Drawing results //
    TH2D *hLambda[3];
    TLegend *leg[3];
    for(int t=0; t<3; t++){
        hLambda[t] = new TH2D(Form("hLambda%s",typeName[t]), Form(";p_{T} (GeV/c);"), 6, 0, 6, 6, -0.4, 0.4);
        ScaleHistoTitle(hLambda[t],0.05,0.9,0.04,0.05,1,0.04,22);
        leg[t] = new TLegend(0.55, 0.69, 0.8, 0.89);
    }
    
    for(int t=0; t<3; t++){
        c1->cd(t+1);
        gPad->SetGrid(1,1);
        hLambda[t]->Draw();
        gPolParDiff[t]->Draw("psame");
        gPolParErrDiff[t]->Draw("psame");
        drawLine(0, 0, 6, 0, 1, 2, 1);
        leg[t]->SetTextFont(22);
        leg[t]->SetTextSize(0.05);
        leg[t]->SetBorderSize(0);
        leg[t]->SetFillColor(0);
        leg[t]->AddEntry(gPolParDiff[t], Form("#lambda_{%s}^{output} - #lambda_{%s}^{input}",polTitle[t], polTitle[t]), "p");
        if(t!=2)leg[t]->AddEntry(gPolParErrDiff[t], Form("<Err(#lambda_{%s})> - #sigma(#lambda_{%s})", polTitle[t], polTitle[t]), "p");
        else leg[t]->AddEntry(gPolParErrDiff[t], Form("<Err(#lambda_{%s})> - RMS(#lambda_{%s})", polTitle[t], polTitle[t]), "p");
        leg[t]->Draw("same");
        if(t==0)drawLatex(0.15, 0.8, "Collins-Soper frame", 22, 0.05, 2);
    }
    
    if(savePDF)c1->SaveAs(Form("CSFinal_iter%d.pdf", nIteration));
    
    if(saveHisto){
        TFile *fout = new TFile("CS_real.root","update");
        fout->cd();
        int nIter;
        if(nIteration==10) nIter = 0;
        else if(nIteration>0 && nIteration<10) nIter = nIteration;
        for(int t=0; t<3; t++){
            gPolParDiff[t]->SetName(Form("gPolParCSIter%d%s", nIter, typeName[t]));
            gPolParErrDiff[t]->SetName(Form("gPolParErrCSIter%d%s", nIter, typeName[t]));
            gPolParDiff[t]->Write();
            if(nIter==0) gPolParErrDiff[t]->Write();
        }
        fout->Close();
    }
    
    for(int t=0; t<3; t++){
        cout<<"********t:  "<<t<<endl;
        for(int p=0; p<nPtBins; p++){
            cout<<diff[t][p]<<endl;
        }
    }
    
    return 0;
}

//--------------------------------------------------------
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
    TLine *l1 = new TLine(xlow,ylow,xup,yup);
    l1->SetLineWidth(lineWidth);
    l1->SetLineColor(lineColor);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
    return l1;
}

//-----------------------------------------
void ScaleHistoTitle(TH1 *h,
                     const Double_t xTitleSize, const Double_t xTitleOffset, const Double_t xLabelSize,
                     const Double_t yTitleSize, const Double_t yTitleOffset, const Double_t yLabelSize,
                     const Int_t font)
{
    if(!h) return;
    h->SetLineColor(1);
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

//-----------------------------------------
void SetHistoMarker(TH2 *h, int marker, int color)
{
    if(!h) return;
    h->SetMarkerStyle(marker);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
}

//-----------------------------------------
void SetGraphMarker(TGraph *g, int marker, int color)
{
    if(!g) return;
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
}

//-----------------------------------------
TPaveText *GetTitleText(TString title, const Float_t size, const Int_t font=22)
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
TPaveText *GetPaveText(Double_t xl, Double_t xh, Double_t yl, Double_t yh, Double_t size, const Int_t font=22)
{
    TPaveText* t1=new TPaveText(xl,yl,xh,yh,"brNDC");
    t1->SetFillStyle(0);
    t1->SetBorderSize(0);
    t1->SetTextSize(size);
    t1->SetTextFont(font);
    return t1;
}

//__________________________________________________
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
//-------------------------------------------------
void calInvariant(double theta, double thetaErr, double phi, double phiErr, double covErr, double &par, double &parErr)
{
    par = (theta + 3.* phi)/(1 - phi);
    
    double err1 = thetaErr * thetaErr;
    double err2 = phiErr * phiErr;
    double err3 = 2 * covErr;
    double dErr1 = (1./(1-phi))*(1./(1-phi));
    double dErr2 = ((3+theta)/(1-phi)/(1-phi))*((3+theta)/(1-phi)/(1-phi));
    double dErr3 = (1./(1-phi))*((3+theta)/(1-phi)/(1-phi));
    parErr = TMath::Sqrt(err1*dErr1+err2*dErr2+err3*dErr3);
}

