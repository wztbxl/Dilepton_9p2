const int pDca[2] = {3, 1};
const int pNfit[2] = {15, 25};
const int pNdedx[2] = {10, 15};
const double pNsigmaPi[3] = {-1.5, -2.0, -2.5};

const char *frameName[2] = {"HX", "CS"};
const char *typeName[3] = {"Cos", "Phi", "inv"};
const char *labelName[3] = {"#lambda_{#theta}", "#lambda_{#varphi}", "#lambda_{inv}"};

TH1D *hRes[2][4][4]; //frame & pt & type
double pRes[2][4][4][24];// frame & type & pt
double pResErr[2][4][4][24];// frame & type & pt
TH1D *hFinalRes[2][4][24];// frame & type & number

const int nPtBins = 4;
double xPtBins[nPtBins+1] = {0, 1, 2, 4, 10};

int draw_Final_Track(int drawChi=1, int drawCheck=0)
{
    int nFile = 0;
    for(int pd=0; pd<2; pd++){
        for(int pNf=0; pNf<2; pNf++){
            for(int pNd=0; pNd<2; pNd++){
                for(int pNs=0; pNs<3; pNs++){
                    TFile *fin = TFile::Open(Form("./Data_%d_%d_%d_%2.1f/Iter_%d_%d_%d_%2.1f/final_%d_%d_%d_%2.1f.root", pDca[pd], pNfit[pNf], pNdedx[pNd], pNsigmaPi[pNs], pDca[pd], pNfit[pNf], pNdedx[pNd], pNsigmaPi[pNs], pDca[pd], pNfit[pNf], pNdedx[pNd], pNsigmaPi[pNs]));
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
                            if(pRes[fr][3][p][nFile]>3) continue;
                            hFinalRes[fr][3][nFile]->SetBinContent(p+1, pRes[fr][3][p][nFile]);
                        }
                    }
                    nFile++;
                }}}}

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
                hMean[f][i]->SetMinimum(-1.5);
                hMean[f][i]->SetMaximum(1.5);
            }
            for(int p=0; p<nPtBins; p++){
                double parMean = getMean(pRes[f][i][p], 24);
                double parMeanErr = getMean(pResErr[f][i][p], 24);
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
                double parMean = getMean(pRes[f][i][p], 24);
                yErr[p] = getRMS(pRes[f][i][p], parMean, 24);
                y[p] = parMean;
            }
            grRMS[f][i] = new TGraphErrors(nPtBins, x, y, xErr, yErr);
            grRMS[f][i]->SetName(Form("grRms%s%s",frameName[f], typeName[i]));
            grRMS[f][i]->SetFillStyle(1001);
            grRMS[f][i]->SetFillColorAlpha(kRed, 0.35);
        }
    }

    //--- for quick test ---//
    if(drawCheck){
        TCanvas *c1 = new TCanvas("c1", "c1", 1600, 600);
        c1->Divide(4,2);
        for(int fr=0; fr<2; fr++){
            for(int t=0; t<4; t++){
                c1->cd(t+1+fr*4);
                for(int i=0; i<nFile; i++){
                    hFinalRes[fr][t][i]->Draw("same");
                }
                if(t!=3) {
                    hMean[fr][t]->Draw("same");
                    grRMS[fr][t]->Draw("e2same");
                }
            }
        }
    }

    //return 0;

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

    //----- for Mean & RMS results -----//
    gStyle->SetOptStat(0);
    TCanvas *c[2];
    int nPad = 1;
    for(int f=0; f<2; f++){
        c[f] = new TCanvas(Form("c_%d",f),Form("%s",frameName[f]),1100,700);
        c[f]->Divide(3,3);
        double padL[3] = {0.01, 0.34, 0.66};
        double padR[3] = {0.31, 0.64, 0.96};
        for(int i=0; i<3; i++){
            double titlePadW = 0.08;
            double padW = (0.98-titlePadW)/2.;
            c[f]->cd(i+1);
            gPad->SetPad(padL[i],titlePadW+padW,padR[i],0.98);
            setPad(0.15,0.05,0.1,0);
            if(nPad==2) drawLatex(0.2, 0.8, "Helicity frame", 22, 0.055, 2);
            if(nPad==5) drawLatex(0.2, 0.8, "Collins-Soper frame", 22, 0.055, 2);
            for(int nF=0; nF<24; nF++){
                //hFinalRes[f][i][nF]->SetTitle(Form("%s;;", labelName[i]));
                setHistoMarker(hFinalRes[f][i][nF], 4, 1);
                if(i!=2) {
                    hFinalRes[f][i][nF]->SetMinimum(-1.5);
                    hFinalRes[f][i][nF]->SetMaximum(1.5);
                }
                else if(i==2) {
                    hFinalRes[f][i][nF]->SetMinimum(-2);
                    hFinalRes[f][i][nF]->SetMaximum(3);
                }
                hFinalRes[f][i][nF]->Draw("same");
            }
            drawLatex(0.55, 0.93, Form("%s", labelName[i]), 22, 0.08, 1);
            TLegend *leg;
            if(i==2){leg = new TLegend(0.6,0.7,0.85,0.8);
                drawLegend(leg, 22, 0.05);
                leg->AddEntry(hFinalRes[f][i][0],"varied cuts","lp");
                leg->Draw("same");}
            gStyle->SetOptStat(0);
            c[f]->cd(i+4);
            gPad->SetPad(padL[i],titlePadW,padR[i],titlePadW+padW);
            setPad(0.15,0.05,0,0.1);
            hMean[f][i]->Draw();
            grRMS[f][i]->Draw("e2same");
            if((i+4)==5 && f==0) drawLatex(0.3, 0.2, "Helicity frame", 22, 0.055, 2);
            if((i+4)==5 && f==1) drawLatex(0.3, 0.2, "Collins-Soper frame", 22, 0.055, 2);
            TLegend *leg2;
            if(i==2){leg2 = new TLegend(0.3,0.8,0.85,0.9);
                drawLegend(leg2, 22, 0.05);
                leg2->AddEntry(hMean[f][i],"mean of all the variations","lp");
                leg2->Draw("same");
                //drawLatex(0.55, 0.7, "(default)", 22, 0.04, 1);
            }
            c[f]->cd(i+7);
            gPad->SetPad(padL[i],0,padR[i],titlePadW);
            setPad(0.15,0.05,0,0.1);
            drawLatex(0.4,0.79,Form("p_{T} (GeV/c)"),22,0.3,1);
        }
        c[f]->SaveAs(Form("sys_cuts_%s.pdf",frameName[f]));
    }
    //---------------------------//

    ofstream outfile("result_track_cuts.txt");
    for(int fr=0; fr<2; fr++){
        outfile<<Form("%s:",frameName[fr])<<endl;
        for(int t=0; t<3; t++){
            outfile<<Form("%s:",typeName[t])<<endl;
            for(int p=0; p<4; p++){
                double parMean = getMean(pRes[fr][t][p], 24);
                double parErr = getMean(pResErr[fr][t][p], 24);
                double parMeanRMS = getRMS(pRes[fr][t][p], parMean, 24);
                outfile<< setprecision(3) <<parMean <<" "<< parErr <<" " << parMeanRMS<<";    ";
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
void setHistoMarker(TH1D *dd, int color, int marker)
{
    dd->SetLineColor(color);
    dd->SetMarkerColor(color);
    dd->SetMarkerStyle(marker);

    dd->GetXaxis()->SetTitleSize(0.055);
    dd->GetXaxis()->SetTitleFont(22);
    dd->GetXaxis()->SetTitleOffset(0.9);
    dd->GetXaxis()->SetLabelSize(0.045);
    dd->GetXaxis()->SetLabelFont(22);

    dd->GetYaxis()->SetTitleSize(0.055);
    dd->GetYaxis()->SetTitleFont(22);
    dd->GetYaxis()->SetTitleOffset(0.8);
    dd->GetYaxis()->SetLabelSize(0.045);
    dd->GetYaxis()->SetLabelFont(22);
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
