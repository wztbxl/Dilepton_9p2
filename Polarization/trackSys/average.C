#include "constant.h"

const int nFiles = 7;
int color[nFiles] = {2, 9, 8, 6, 3, 4, 2};
int marker[nFiles] = {21, 22, 22, 20, 28, 28, 25};

TH1D *hResult[2][nPtBins][nFiles][2];//frame & nPt & angle & files

int inFileName[2][nPtBins][6][4] = {
    {
        {{4,4,0,1}, {4,3,0,1}, {4,5,0,1}, {2,4,0,1}, {4,4,2,1}, {4,4,1,1}},
        {{4,4,0,1}, {4,3,0,1}, {4,5,0,1}, {2,4,0,1}, {4,4,2,1}, {4,4,1,1}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,0,0}, {2,2,0,0}, {4,2,2,0}, {4,2,1,0}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,0,0}, {2,2,0,0}, {4,2,2,0}, {4,2,1,0}}
    },//HX
    {
        {{4,4,0,1}, {4,3,0,1}, {4,5,0,1}, {2,4,0,1}, {4,4,4,1}, {4,4,3,1}},
        {{4,4,0,1}, {4,3,0,1}, {4,5,0,1}, {2,4,0,1}, {4,4,4,1}, {4,4,3,1}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,0,0}, {2,2,0,0}, {4,2,4,0}, {4,2,3,0}},
        {{4,2,0,0}, {4,2,0,1}, {4,3,0,0}, {2,2,0,0}, {4,2,4,0}, {4,2,3,0}}
    }//CS
};

int average(double dca=varDca[1], double nHitsFit=varNhitsFit[0], double nHitsDedx=varNhitsDedx[1], double nSigmaPi=varNsigmaPi[0])
{
    char buf[1024] = Form("_%1.0f_%2.0f_%2.0f_%2.1f",dca, nHitsFit, nHitsDedx, nSigmaPi);
    TFile *fin = new TFile(Form("./Data%s/All%s.root", buf, buf),"read");
    
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int an=0; an<2; an++){
            for(int nF=0; nF<nFiles; nF++){
                    if(nF<=nFiles-2){
                        hResult[f][p][nF][an] = (TH1D*)fin->Get(Form("h%sPt%d%sFitRebin%dPol%dShift%dfix%d", frameName[f], p, angleName[an], inFileName[f][p][nF][0], inFileName[f][p][nF][1], inFileName[f][p][nF][2], inFileName[f][p][nF][3]));
                        setHistPoint(hResult[f][p][nF][an], marker[nF], color[nF]);
                    }
                    else{
                        hResult[f][p][nF][an] = (TH1D*)fin->Get(Form("h%sPt%d%sCountRebin%dPol%dShift%dfix%d", frameName[f], p, angleName[an], inFileName[f][p][0][0], inFileName[f][p][0][1], inFileName[f][p][0][2], inFileName[f][p][0][3]));
                        setHistPoint(hResult[f][p][nF][an], marker[6], color[6]);
                    }
                }
            }
        }
    }
    
    TH1D *hNew[2][nPtBins][2];// frame & nPt & angue

    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int an=0; an<2; an++){
                hNew[f][p][an]=(TH1D*)hResult[f][p][0][an]->Clone(Form("%s_clone", hResult[f][p][0][an]->GetName()));
                hNew[f][p][an]->SetName(Form("h%s%sPt%d", frameName[f], angleName[an], p));
                setHistPoint(hNew[f][p][an], 29, 1);
                int nBins = hResult[f][p][0][an]->GetNbinsX();
                double **counts = new double*[nBins];
                double **countsErr = new double*[nBins];
                double *mean = new double[nBins];
                double *meanErr = new double[nBins];
                for(int i=0; i<nBins; i++){
                    counts[i] = new double[nFiles];
                    countsErr[i] = new double[nFiles];
                    for(int nF=0; nF<nFiles; nF++){
                        counts[i][nF] = hResult[f][p][nF][an]->GetBinContent(i+1);
                        countsErr[i][nF] = hResult[f][p][nF][an]->GetBinError(i+1);
                       
                        if( (i<4 && an==1) && (counts[i][nF]>300 || counts[i][nF]<-0.1) ){
                            counts[i][nF] = -1;
                            countsErr[i][nF] = -1;
                            hResult[f][p][nF][an]->SetBinContent(i+1, 0);
                            hResult[f][p][nF][an]->SetBinError(i+1, 0);
                        }
                        if(countsErr[i][nF]/counts[i][nF]>2.5){
                            counts[i][nF] = -1;
                            countsErr[i][nF] = -1;
                            hResult[f][p][nF][an]->SetBinContent(i+1, 0);
                            hResult[f][p][nF][an]->SetBinError(i+1, 0);
                        }
                        if(TMath::IsNaN(countsErr[i][nF])){
                            counts[i][nF] = -1;
                            countsErr[i][nF] = -1;
                            hResult[f][p][nF][an]->SetBinContent(i+1, 0);
                            hResult[f][p][nF][an]->SetBinError(i+1, 0);
                        }
                    }
                    mean[i] = (double)getMean(counts[i], nFiles);
                    meanErr[i] = (double)getMean(countsErr[i], nFiles);
                }
                
                for(int i=0; i<nBins; i++){
                    if(mean[i]-1<0) continue;
                    hNew[f][p][an]->SetBinContent(i+1, mean[i]);
                    hNew[f][p][an]->SetBinError(i+1, meanErr[i]);
                }
            }
        }
    }
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(4, 4);
    int nPad = 1;
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int an=0; an<2; an++){
                //cout<<"c1->cd("<<nPad<<"): "<<endl;
                c1->cd(nPad);
                hResult[f][p][0][an]->Draw();
                for(int nF=1; nF<nFiles; nF++){
                    if(!hResult[f][p][nF][an]){
                        cout<<f<<", "<<p<<", "<<an<<", "<<nF<<", "<<endl;
                        continue;
                    }
                    hResult[f][p][nF][an]->Draw("same");
                }
                hNew[f][p][an]->Draw("same");
                nPad++;
            }
        }
    }
    c1->SaveAs(Form("Data%s/Result%s.pdf", buf, buf));
    
    TFile *fout = new TFile(Form("Data%s/Result%s.root", buf, buf), "recreate");
    fout->cd();
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int an=0; an<2; an++){
                hNew[f][p][an]->Write();
            }
        }
    }
    
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
double getMean(double *p, const int n)
{
    double mean = 0;
    int bn = 0;
    for(int i=0; i<n; i++){
        if(p[i]<=0) {continue;}
        mean += p[i];
        bn++;
    }
    if(bn!=0)mean = mean/bn;
    return mean;
}







