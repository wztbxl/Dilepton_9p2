#include "constant.h"

const int nFiles = 12;
int color[nFiles] = { 2,   9,  8,  6, 38,  4,  2,  2,  9,  8,  6,  8};
int marker[nFiles] = {21, 22, 22, 20, 28, 28, 25, 22, 20, 28, 28, 25};

TH1D *hResult[2][nPtBins][nFiles][2];//frame & nPt & angle & files

int inFileName[2][nPtBins][nFiles-1][4] = {
    {
        {{4,4,0,1}, {4,3,0,1}, {4,3,1,1}, {4,3,2,1}, {2,3,0,1}, {4,5,0,1}, {2,4,0,1}, {2,4,1,1}, {2,4,2,1}, {4,4,2,1}, {4,4,1,1}},
        {{4,4,0,1}, {4,3,0,1}, {4,3,1,1}, {4,3,2,1}, {2,3,0,1}, {4,5,0,1}, {2,4,0,1}, {2,4,1,1}, {2,4,2,1}, {4,4,2,1}, {4,4,1,1}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,1,0}, {4,3,2,0}, {2,3,0,0}, {4,3,0,0}, {2,2,0,0}, {2,2,1,0}, {2,2,2,0}, {4,2,2,0}, {4,2,1,0}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,1,0}, {4,3,2,0}, {2,3,0,0}, {4,3,0,0}, {2,2,0,0}, {2,2,1,0}, {2,2,2,0}, {4,2,2,0}, {4,2,1,0}}
    },//HX
    {
        {{4,4,0,1}, {4,3,0,1}, {4,3,3,1}, {4,3,4,1}, {2,3,0,1}, {4,5,0,1}, {2,4,0,1}, {2,4,3,1}, {2,4,4,1}, {4,4,4,1}, {4,4,3,1}},
        {{4,4,0,1}, {4,3,0,1}, {4,3,3,1}, {4,3,4,1}, {2,3,0,1}, {4,5,0,1}, {2,4,0,1}, {2,4,3,1}, {2,4,4,1}, {4,4,4,1}, {4,4,3,1}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,3,0}, {4,3,4,0}, {2,3,0,0}, {4,3,0,0}, {2,2,0,0}, {2,2,3,0}, {2,2,4,0}, {4,2,4,0}, {4,2,3,0}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,3,0}, {4,3,4,0}, {2,3,0,0}, {4,3,0,0}, {2,2,0,0}, {2,2,3,0}, {2,2,4,0}, {4,2,4,0}, {4,2,3,0}}
    }//CS
};


const char *legName[nFiles] = {
    "Gaus+Pol4/Pol2; bin:40 MeV/c^{2}; fixed/free",
    "Gaus+Pol3/Pol2; bin:40 MeV/c^{2}; fixed/fixed",
    "Gaus+Pol3; bin:40 MeV/c^{2}; fixed/free; [x-0.1, 3.9]",
    "Gaus+Pol3; bin:40 MeV/c^{2}; fixed/free; [x+0.1, 4.1]",
    "Gaus+Pol3; bin:20 MeV/c^{2}; fixed/free",
    "Gaus+Pol5/Pol3; bin:40 MeV/c^{2}; fixed/free",
    "Gaus+Pol4; bin:20 MeV/c^{2}; fixed/free",
    "Gaus+Pol4; bin:20 MeV/c^{2}; fixed/free; [x-0.1, 3.9]",
    "Gaus+Pol4; bin:20 MeV/c^{2}; fixed/free; [x+0.1, 4.1]",
    "Gaus+Pol4; bin:40 MeV/c^{2}; fixed; [x-0.1, 3.9]",
    "Gaus+Pol4; bin:40 MeV/c^{2}; fixed; [x+0.1, 4.1]",
    "Gaus+Pol4; bin:40 MeV/c^{2}; fixed/free; counting "
};

int debug = 0;

int getSigFile(double dca=varDca[0], double nHitsFit=varNhitsFit[0], double nHitsDedx=varNhitsDedx[0], double nSigmaPi=varNsigmaPi[1])
{
    char buf[1024] = Form("_%1.0f_%2.0f_%2.0f_%2.1f",dca, nHitsFit, nHitsDedx, nSigmaPi);
    TFile *fin = TFile::Open(Form("./Data%s/All%s.root", buf, buf));

    for(int nF=0; nF<nFiles; nF++){
        TFile *fout = new TFile(Form("File%d/systematic_Signal.root", nF), "recreate");
        fout->cd();
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int an=0; an<2; an++){
                    if(nF<=nFiles-2){
                        if(debug) cout<<f<<p<<nF<<an<<" "<<Form("h%sPt%d%sFitRebin%dPol%dShift%dfix%d", frameName[f], p, angleName[an], inFileName[f][p][nF][0], inFileName[f][p][nF][1], inFileName[f][p][nF][2], inFileName[f][p][nF][3])<<" * "<<endl;
                        hResult[f][p][nF][an] = (TH1D*)fin->Get(Form("h%sPt%d%sFitRebin%dPol%dShift%dfix%d", frameName[f], p, angleName[an], inFileName[f][p][nF][0], inFileName[f][p][nF][1], inFileName[f][p][nF][2], inFileName[f][p][nF][3]));
                        setHistPoint(hResult[f][p][nF][an], marker[nF], color[nF]);
                        if(debug) cout<< hResult[f][p][nF][an]<<endl;
                    }
                    else{
                        if(debug) cout<<f<<p<<nF<<an<<" "<<Form("h%sPt%d%sCountRebin%dPol%dShift%dfix%d", frameName[f], p, angleName[an], inFileName[f][p][nF][0], inFileName[f][p][0][1], inFileName[f][p][0][2], inFileName[f][p][0][3])<<" * "<<endl;
                        hResult[f][p][nF][an] = (TH1D*)fin->Get(Form("h%sPt%d%sCountRebin%dPol%dShift%dfix%d", frameName[f], p, angleName[an], inFileName[f][p][0][0], inFileName[f][p][0][1], inFileName[f][p][0][2], inFileName[f][p][0][3]));
                        setHistPoint(hResult[f][p][nF][an], marker[nF], color[nF]);
                        if(debug) cout<<hResult[f][p][nF][an]<<endl;
                    }
                hResult[f][p][nF][an]->SetName(Form("h%s%sPt%d", frameName[f], angleName[an], p));
                hResult[f][p][nF][an]->Write();
                }
            }
        }
        fout->Close();
    }
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(4, 4);
    int nPad = 1;
    for(int f=0; f<2; f++){
        for(int an=0; an<2; an++){
        for(int p=0; p<nPtBins; p++){
                c1->cd(nPad);
            gPad->SetLeftMargin(0.2);
            if(an==0) {
                hResult[f][p][0][an]->SetTitle(Form("%s: %2.0f < p_{T} < %2.0f GeV/c", frameName[f], xPtBins[p], xPtBins[p+1]));
            }
            hResult[f][p][0][an]->SetTitleOffset(1.5,"y");
            hResult[f][p][0][an]->GetYaxis()->SetTitleSize(0.05);
            hResult[f][p][0][an]->GetXaxis()->SetTitleSize(0.05);
            if(f==1 && p==0 && an==0){
                hResult[f][p][0][an]->SetBinContent(3, 0);
                hResult[f][p][0][an]->SetBinError(3, 0);
                hResult[f][p][0][an]->SetBinContent(8, 0);
                hResult[f][p][0][an]->SetBinError(8, 0);
            }
                hResult[f][p][0][an]->Draw();
                for(int nF=1; nF<nFiles; nF++){
                    if(!hResult[f][p][nF][an]){
                        cout<<f<<", "<<p<<", "<<an<<", "<<nF<<", "<<endl;
                        continue;
                    }
                    if(f==1 && p==0 && an==0){
                        hResult[f][p][nF][an]->SetBinContent(3, 0);
                        hResult[f][p][nF][an]->SetBinError(3, 0);
                        hResult[f][p][nF][an]->SetBinContent(8, 0);
                        hResult[f][p][nF][an]->SetBinError(8, 0);
                    }
                    hResult[f][p][nF][an]->Draw("same");
                }
                nPad++;
            }
        }
    }
    c1->SaveAs(Form("sys_signal_counts.pdf"));
    
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    TLegend *leg = new TLegend(0.1,0.1,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    for(int nF=0; nF<nFiles; nF++){
        leg->AddEntry(hResult[0][0][nF][0], Form("%s",legName[nF]),"lp");
    }
    c2->cd();
    leg->Draw();
    c2->SaveAs(Form("sys_signal_legend.pdf"));
    return 0;
}

//------------------------------
void setHistPoint(TH1D *h, int marker, int color)
{
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
}

//------------------------------
double getMean(double *p, const int n, int debug=0)
{
    double mean = 0;
    if(debug) cout<<"*******"<<endl;
    int bn = 0;
    for(int i=0; i<n; i++){
        if(debug)cout<<p[i]<<endl;
        if(p[i]<0) {debug=1;continue;}
        mean += p[i];
        bn++;
    }
    if(bn!=0)mean = mean/bn;
    if(debug){
        cout<<"mean:"<<mean<<" bn:"<<bn<<endl;
        cout<<"*******"<<endl;}
    return mean;
}







