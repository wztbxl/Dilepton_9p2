const char *frameName[2] = {"HX", "CS"};
const char *typeName[3] = {"Cos", "Phi", "inv"};
const char *labelName[3] = {"#lambda_{#theta}", "#lambda_{#varphi}", "#lambda_{inv}"};

int color[12] = { 2,   9,  8,  6, 38,  4,  2,  2,  9,  8,  6,  8};
int marker[12] = {21, 22, 22, 20, 28, 28, 25, 22, 20, 28, 28, 25};

TH1D *hRes[2][4][4]; //frame & pt & type
double pRes[2][4][4][12];// frame & type & pt
double pResErr[2][4][4][12];// frame & type & pt
TH1D *hFinalRes[2][4][12];// frame & type & number

const int nPtBins = 4;
double xPtBins[nPtBins+1] = {0, 1, 2, 4, 10};

int draw_Final_signal(int drawChi=1)
{
    int nFile = 0;
    for(int i=0; i<12; i++){
        TFile *fin = TFile::Open(Form("File%d/Iter_f%d/final_%d.root", i, i, i));
        for(int p=0; p<4; p++){
            for(int fr=0; fr<2; fr++){
                for(int t=0; t<3; t++){
                    hRes[fr][p][t] = (TH1D*)fin->Get(Form("hIter%s%sPt%d", frameName[fr], typeName[t], p));
                }
                hRes[fr][p][3] = (TH1D*)fin->Get(Form("hIterChi%sPt%d", frameName[fr], p));
                int binTheta = getConverBin(hRes[fr][p][0]);
                int binPhi = getConverBin(hRes[fr][p][1]);
                int binChose = max(binTheta, binPhi);
                for(int t=0; t<3; t++){
                    pRes[fr][t][p][nFile] = hRes[fr][p][t]->GetBinContent(binChose);
                    pResErr[fr][t][p][nFile] = hRes[fr][p][t]->GetBinError(binChose);
                }
                pRes[fr][3][p][nFile] = hRes[fr][p][3]->GetBinContent(binChose);
            }
        }

        for(int fr=0; fr<2; fr++){
            for(int t=0; t<3; t++){
                hFinalRes[fr][t][nFile] = new TH1D(Form("hFinal%s%s_%d", frameName[fr], typeName[t], nFile), Form("; p_{T} GeV/c; "), nPtBins, xPtBins);
                for(int p=0; p<4; p++){
                    hFinalRes[fr][t][nFile]->SetBinContent(p+1, pRes[fr][t][p][nFile]);
                    hFinalRes[fr][t][nFile]->SetBinError(p+1, pResErr[fr][t][p][nFile]);
                }
            }
            hFinalRes[fr][3][nFile] = new TH1D(Form("hFinal%sChi_%d", frameName[fr], nFile), Form("%s; p_{T} GeV/c; #Chi^{2}/ndf", frameName[fr]), nPtBins, xPtBins);
            for(int p=0; p<4; p++){
                hFinalRes[fr][3][nFile]->SetBinContent(p+1, pRes[fr][3][p][nFile]);
            }
        }
        nFile++;
    }

    //--- for draw together ---//
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "c1", 1300, 700);
    c1->Divide(3,2);
    for(int fr=0; fr<2; fr++){
        for(int t=0; t<3; t++){
            c1->cd(t+1+fr*3);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            for(int i=0; i<nFile; i++){
                if(t<2){
                    hFinalRes[fr][t][i]->SetMinimum(-1.5);
                    hFinalRes[fr][t][i]->SetMaximum(1.5);
                }
                else if(t==2){
                    hFinalRes[fr][t][i]->SetMinimum(-2);
                    hFinalRes[fr][t][i]->SetMaximum(3);
                }
                hFinalRes[fr][t][i]->GetYaxis()->SetTitle(Form("%s", labelName[t]));
                setHistoMarker(hFinalRes[fr][t][i], color[i], marker[i], 1);
                hFinalRes[fr][t][i]->Draw("same");
            }
        }
    }
    c1->cd(2); drawLatex(0.2, 0.8, "Helicity frame", 22, 0.055, 2);
    c1->cd(5); drawLatex(0.2, 0.8, "Collins-Soper frame", 22, 0.055, 2);
    c1->SaveAs("sys_signal_12.pdf");
    //----------------------//

    TH1D *hMean[2][3];
    double x[nPtBins]={0};
    double xShift[nPtBins]={0};
    double xErr[nPtBins]={0};
    double xErrShift[nPtBins]={0};
    for(int i=0; i<4; i++){
        double ptBinWidth = (xPtBins[i+1]-xPtBins[i])/2.;
        x[i]=xPtBins[i]+ptBinWidth;
        xErr[i] = 0.35;
        xShift[i] = x[i] + 0.2;
        //cout<<x[i]<<endl;
    }
    TGraphErrors *grRMS[2][3];
    TGraphErrors *grMeanRMS[2][3];
    TGraphErrors *grMean[2][3];

    for(int f=0; f<2; f++){
        for(int i=0; i<3; i++){
            hMean[f][i] = new TH1D(Form("hMean%s%s", frameName[f], typeName[i]), Form(""), nPtBins, xPtBins);
            if(i==2){
                hMean[f][i]->SetMinimum(-2);
                hMean[f][i]->SetMaximum(3);
            }
            else{
                hMean[f][i]->SetMinimum(-1);
                hMean[f][i]->SetMaximum(1);
            }
            for(int p=0; p<nPtBins; p++){
                double parMean = getMean(pRes[f][i][p], 12);
                double parMeanErr = getMean(pResErr[f][i][p], 12);
                hMean[f][i]->SetBinContent(p+1, parMean);
                hMean[f][i]->SetBinError(p+1, parMeanErr);
                setHistoMarker(hMean[f][i], 2, 21);
            }
        }
    }

    double y[4] = {0};
    double yErr[4] = {0};
    for(int f=0; f<2; f++){
        for(int i=0; i<3; i++){
            for(int p=0; p<4; p++){
                double parMean = getMean(pRes[f][i][p], 12);
                yErr[p] = getRMS(pRes[f][i][p], parMean, 12);
                y[p] = parMean;
            }
            grRMS[f][i] = new TGraphErrors(nPtBins, x, y, xErr, yErr);
            grRMS[f][i]->SetName(Form("grRms%s%s",frameName[f], typeName[i]));
            grRMS[f][i]->SetFillStyle(1001);
            grRMS[f][i]->SetFillColorAlpha(kRed, 0.35);
        }
    }

    gStyle->SetOptStat(0);
    //----- for chi^2/ndf -----//
    if(drawChi){
        TCanvas *c3 = new TCanvas("c3","c3", 800, 300);
        c3->Divide(2,1);
        for(int fr=0; fr<2; fr++){
            c3->cd(fr+1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            for(int i=0; i<nFile; i++){
                hFinalRes[fr][3][i]->SetTitleSize(0.08);
                hFinalRes[fr][3][i]->GetXaxis()->SetLabelSize(0.06);
                hFinalRes[fr][3][i]->GetYaxis()->SetLabelSize(0.06);
                hFinalRes[fr][3][i]->GetXaxis()->SetTitleSize(0.06);
                hFinalRes[fr][3][i]->GetYaxis()->SetTitleSize(0.06);
                hFinalRes[fr][3][i]->GetXaxis()->CenterTitle();
                hFinalRes[fr][3][i]->GetYaxis()->CenterTitle();
                hFinalRes[fr][3][i]->SetMinimum(0);
                hFinalRes[fr][3][i]->SetMaximum(3);
                hFinalRes[fr][3][i]->Draw("same");
            }
        }
        c3->SaveAs("sys_cuts_chi.pdf");
    }
    //------------------------//


    ofstream outfile("result_signal_fitting.txt");
    for(int fr=0; fr<2; fr++){
        outfile<<Form("%s:",frameName[fr])<<endl;
        for(int t=0; t<3; t++){
            outfile<<Form("%s:",typeName[t])<<endl;
            for(int p=0; p<4; p++){
                double parMean = getMean(pRes[fr][t][p], 12);
                double parErr = getMean(pResErr[fr][t][p], 12);
                double parMeanRMS = getRMS(pRes[fr][t][p], parMean, 12);
                outfile<< setprecision(3) << parMeanRMS<<";    ";
            }
            outfile<<endl;
        }
    }
    outfile.close();

    return 0;
}

//-------------------------
int getConverBin(TH1D *h1, int debug=0)
{
    int bin = 0;
    int nBins = h1->GetNbinsX();
    for(int i=1; i<=nBins-1; i++){
        double b1 = h1->GetBinContent(i);
        double b2 = h1->GetBinContent(i+1);
        if((b2-b1)<0.01 && (b2-b1)>-0.01) {
            if(debug)cout<<b2-b1<<endl;
            bin = i+1;
            break;
        }
    }
    return bin;
}

//----------------------------
int max(int a, int b)
{
    if(a>b || a==b) return a;
    else return b;
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

//__________________________________________________
void setPad(Double_t left, Double_t right, Double_t top, Double_t bottom, int color=10)
{
    gPad->SetFillColor(color);
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

//-----------------------------------
double getMean(double *p, const int n, int debug=0)
{
    double mean = 0;
    if(debug) cout<<"*******"<<endl;
    int bn = 0;
    for(int i=0; i<n; i++){
        if(debug)cout<<p[i]<<endl;
        mean += p[i];
        bn++;
    }
    if(bn!=0)mean = mean/bn;
    if(debug){
        cout<<"mean:"<<mean<<" bn:"<<bn<<endl;
        cout<<"*******"<<endl;}
    return mean;
}

//------------------------------------
double getRMS(double *a, double mean, const int loop, bool debug=false)
{
    double max=0;
    double b = 0;
    for(int i=0; i<loop; i++){
        b = a[i] - mean;
        max = max + b * b;
        if(debug) cout<<a[i]<<"; "<<mean<<"; "<<b<<"; "<<b*b<<max<<endl;
    }
    max = max/loop;
    if(debug) cout<<"RMS:"<<sqrt(max)<<endl;
    return sqrt(max);
}

//------------------------------------
void setHistoMarker(TH1D *dd, int color, int marker, int center=0)
{
    dd->SetLineColor(color);
    dd->SetMarkerColor(color);
    dd->SetMarkerStyle(marker);

    dd->GetXaxis()->SetTitleSize(0.055);
    dd->GetXaxis()->SetTitleFont(22);
    dd->GetXaxis()->SetTitleOffset(0.9);
    dd->GetXaxis()->SetLabelSize(0.045);
    dd->GetXaxis()->SetLabelFont(22);

    dd->GetYaxis()->SetTitleSize(0.06);
    dd->GetYaxis()->SetTitleFont(22);
    dd->GetYaxis()->SetTitleOffset(1.0);
    dd->GetYaxis()->SetLabelSize(0.045);
    dd->GetYaxis()->SetLabelFont(22);
    
    if(center){
        dd->GetXaxis()->CenterTitle();
        dd->GetYaxis()->CenterTitle();
    }
}

//--------------------------------------------------------
TLegend* drawLegend(TLegend* leg, Int_t textFont=22, Double_t textSize=0.04)
{
    leg->SetTextFont(textFont);
    //leg->SetTextAlign(11);
    //align = 10*HorizontalAlign + VerticalAlign (1=left adjusted, 2=centered, 3=right adjusted)
    leg->SetTextSize(textSize);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    return leg;
}
