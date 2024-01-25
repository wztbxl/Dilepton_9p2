//--- for analysis note version v9

const int debug=0;
const char *frameName[2] ={"Hx", "Cs"};
const char *typeName[3] = {"Th", "Ph", "Inv1"};

int saveDimuonGraph()
{
    double result[2][3][4];// frame type pT
    double statErr[2][3][4];
    double sysTrac[2][3][4];
    double sysFit[2][3][4];
    double sysTotal[2][3][4];

    //---- read data ----//
    string tmpS;
    for(int fr=0; fr<2; fr++){
        ifstream fRes(Form("correct_factor_%s.txt", frameName[fr]));
        for(int i=0; i<6; i++){
            getline(fRes,tmpS);
        }
        for(int t=0; t<3; t++){
            getline(fRes,tmpS); fRes>>result[fr][t][0]>>statErr[fr][t][0]>>tmpS>>result[fr][t][1]>>statErr[fr][t][1]>>tmpS>>result[fr][t][2]>>statErr[fr][t][2]>>tmpS>>result[fr][t][3]>>statErr[fr][t][3]>>tmpS>>tmpS;
        }
    }
    if(debug){
        for(int fr=0; fr<2; fr++){
            for(int t=0; t<3; t++){
                for(int p=0; p<4; p++){
                    cout<<result[fr][t][p]<<"  "<<statErr[fr][t][p]<<"  "<<endl;
                }
            }
        }}
    
    ifstream fTrk("result_track_cuts.txt");
    for(int fr=0; fr<2; fr++){
        getline(fTrk,tmpS);
        for(int t=0; t<3; t++){
            getline(fTrk,tmpS); fTrk>>tmpS>>tmpS>>sysTrac[fr][t][0]>>tmpS>>tmpS>>tmpS>>sysTrac[fr][t][1]>>tmpS>>tmpS>>tmpS>>sysTrac[fr][t][2]>>tmpS>>tmpS>>tmpS>>sysTrac[fr][t][3]>>tmpS>>tmpS;
        }
    }
    if(debug){
        for(int fr=0; fr<2; fr++){
            for(int t=0; t<3; t++){
                for(int p=0; p<4; p++){
                    cout<<sysTrac[fr][t][p]<<endl;
                }
            }
        }}

    ifstream fFit("result_signal_fitting.txt");
    for(int fr=0; fr<2; fr++){
        getline(fFit,tmpS);
        for(int t=0; t<3; t++){
            getline(fFit,tmpS);
            for(int p=0; p<4; p++){
                fFit>>sysFit[fr][t][p]>>tmpS;
            }
            fFit>>tmpS;
        }
    }
    if(debug){
        for(int fr=0; fr<2; fr++){
            for(int t=0; t<3; t++){
                for(int p=0; p<4; p++){
                    cout<<sysFit[fr][t][p]<<"  ";
                }
                cout<<endl;
            }
        }}

    //--- save results: iteraction ----//
    ofstream outfile("result_sys_total.txt");
    for(int fr=0; fr<2; fr++){
        outfile<<Form("%s:",frameName[fr])<<endl;
        for(int t=0; t<3; t++){
            outfile<<Form("%s:",typeName[t])<<endl;
            for(int p=0; p<4; p++){
                sysTotal[fr][t][p] = sqrt(pow(sysTrac[fr][t][p],2) + pow(sysFit[fr][t][p],2));
                outfile<< setprecision(3) << sysTotal[fr][t][p] <<";    ";
            }
            outfile<<endl;
        }
    }
    outfile.close();
    
    
    //-- save data graph using TGraphAsymmErrors--//
    const int nPtBins = 4;
    double xPtBins[nPtBins+1] = {0, 1, 2, 4, 10};
    //double xMean[nPtBins] = {0.65, 1.47, 2.68, 4.92}; // preliminary
    double xMean[nPtBins] = {0.62, 1.46, 2.69, 4.92}; //published
    double exl[nPtBins], exh[nPtBins];
    for(int i=0; i<nPtBins; i++){
        exl[i] = xMean[i] - xPtBins[i];
        exh[i] = xPtBins[i+1] - xMean[i];
    }

    TGraphAsymmErrors *gErrAsy_dimuon_stat[2][3], *gErrAsy_dimuon_sys[2][3], *gErrAsy_dimuon_tot[2][3];
    for(int fr=0; fr<2; fr++){
        for(int t=0; t<3; t++){
            gErrAsy_dimuon_stat[fr][t] = new TGraphAsymmErrors(nPtBins, xMean, result[fr][t], exl, exh, statErr[fr][t], statErr[fr][t]);
            gErrAsy_dimuon_stat[fr][t]->SetName(Form("gErrAsy_dimuon_stat_%s_%s", frameName[fr], typeName[t]));
            gErrAsy_dimuon_sys[fr][t] = new TGraphAsymmErrors(nPtBins, xMean, result[fr][t], exl, exh, sysTotal[fr][t], sysTotal[fr][t]);
            gErrAsy_dimuon_sys[fr][t]->SetName(Form("gErrAsy_dimuon_sys_%s_%s", frameName[fr], typeName[t]));
            double total[4];
            for(int p=0; p<4; p++){
                total[p] = sqrt(pow(statErr[fr][t][p],2) + pow(sysTotal[fr][t][p],2));
            }
            gErrAsy_dimuon_tot[fr][t] = new TGraphAsymmErrors(nPtBins, xMean, result[fr][t], exl, exh, total, total);
            gErrAsy_dimuon_tot[fr][t]->SetName(Form("gErrAsy_dimuon_tot_%s_%s", frameName[fr], typeName[t]));
        }
    }
    
    //----- save final result -----//
    TFile *fout = new TFile("dimuonResult.root", "recreate");
    fout->cd();
    for(int fr=0; fr<2; fr++){
        for(int t=0; t<3; t++){
            gErrAsy_dimuon_stat[fr][t]->Write();
            gErrAsy_dimuon_sys[fr][t]->Write();
            gErrAsy_dimuon_tot[fr][t]->Write();
        }
    }


    return 0;
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

void setGraphStyle(TGraphErrors* gr, int color, int marker, int mSize, int lSize)
{
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(mSize);
    gr->SetLineWidth(lSize);
}

//__________________________________________________
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
    TLine *l1 = new TLine(xlow,ylow,xup,yup);
    l1->SetLineWidth(lineWidth);
    l1->SetLineColor(lineColor);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
    return l1;
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
