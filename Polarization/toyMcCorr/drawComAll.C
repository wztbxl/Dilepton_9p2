//-- for draw the results of the corrected result
#include "drawComFunc.h"

const int nPtBins = 4;
const double lowPtBins[nPtBins] = {0,1,2,4};
const double highPtBins[nPtBins] = {1,2,4,10};

const int debug = 0;

int drawComAll()
{
    double result[2][3][4];// frame type pT
    double statErr[2][3][4];

    //---- read data ----//
    string tmpS;
    ifstream fRes("result_track_cuts.txt");
    for(int fr=0; fr<2; fr++){
        getline(fRes,tmpS);
        for(int t=0; t<3; t++){
            getline(fRes,tmpS); fRes>>result[fr][t][0]>>statErr[fr][t][0]>>tmpS>>result[fr][t][1]>>statErr[fr][t][1]>>tmpS>>result[fr][t][2]>>statErr[fr][t][2]>>tmpS>>result[fr][t][3]>>statErr[fr][t][3]>>tmpS>>tmpS;
        }
    }
    if(debug){
        for(int fr=0; fr<2; fr++){
            for(int t=0; t<3; t++){
                for(int p=0; p<4; p++){
                    cout<<result[fr][t][p]<<"  "<<statErr[fr][t][p]<<endl;
                }
            }
        }}

    for(int fr=0; fr<2; fr++){ //frames

        gStyle->SetOptStat(0);
        TCanvas *c[2];
        for(int i=0; i<2; i++){
            c[i] = new TCanvas(Form("c%d", i), Form("c%d", i), 1100, 500);
            c[i]->Divide(4, 2);
        }
        
        TH2D *hPol[2];
        for(int t=0; t<2; t++){
            if(t==0)hPol[t] = new TH2D(Form("hPol%s", typeName[t]), Form("; input #lambda_{%s}; measured #lambda_{%s}", polTitle[t], polTitle[t]), 1, -1., 1., 1, -1., 1.);
            else hPol[t] = new TH2D(Form("hPol%s", typeName[t]), Form("; input #lambda_{%s}; measured #lambda_{%s}", polTitle[t], polTitle[t]), 1, -0.6, 0.6, 1, -0.6, 0.6);
            SetHistoMarker(hPol[t], 20, 1);
        }
        

        TH2D *hErr[2];
        for(int t=0; t<2; t++){
            hErr[t] = new TH2D(Form("hErr%s", typeName[t]), Form("; input #sigma_{#lambda_{%s}}; measured #sigma_{#lambda_{%s}}", polTitle[t], polTitle[t]), 1, 0, 1.5, 1, 0, 1.5);
            SetHistoMarker(hErr[t], 20, 1);
        }

        TGraph *gPolParCom[2][4];
        TGraph *gPolParComErr[2][4];

        TFile *fin[2];
        for(int t=0; t<2; t++){
            fin[t] = TFile::Open(Form("%s_%s_iter.root", frameName[fr], typeName[t]));
        }

        TF1 *flinear[2][4];
        TF1 *flinearErr[2][4];
        TF1 *f1 = new TF1("f1", "[0]*x", -1, 2);
        f1->SetParameter(0., 1.);
        f1->SetLineStyle(2);
        f1->SetLineColor(1);
        f1->SetLineWidth(1);

        gStyle->SetOptFit(111);
        gStyle->SetStatY(0.9);
        gStyle->SetStatX(0.53);
        //gStyle->SetStatW(0.4);
        //gStyle->SetStatH(0.2);
        for(int p=0; p<nPtBins; p++){
            for(int t=0; t<2; t++){
                flinear[t][p] = new TF1(Form("flinear%sPt%d",typeName[t], p), "[0]*x+[1]", -1, 1);
                flinear[t][p]->SetLineColor(4);
                flinear[t][p]->SetLineWidth(1);
                flinearErr[t][p] = new TF1(Form("flinearErr%sPt%d",typeName[t], p), "[0]*x", -2, 2);
                flinearErr[t][p]->SetLineColor(4);
                flinearErr[t][p]->SetLineWidth(1);
                //-- iter
                gPolParCom[t][p] = (TGraph*)fin[t]->Get(Form("gPolParCom%sPt%d", typeName[t], p));
                gPolParCom[t][p]->SetMarkerSize(2);
                gPolParComErr[t][p] = (TGraph*)fin[t]->Get(Form("gPolParErr%sPt%d", typeName[t], p));
                gPolParComErr[t][p]->SetMarkerSize(2);

                c[0]->cd(t*4+p+1);
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.15);
                hPol[t]->Draw();
                f1->Draw("same");

                if(t==0) drawLatex(0.3, 0.94, Form(" %2.0f < p_{T} < %2.0f GeV/c", lowPtBins[p], highPtBins[p]), 22, 0.06, 1);
                gPolParCom[t][p]->Fit(flinear[t][p]);
                flinear[t][p]->Draw("same");

                gPolParCom[t][p]->Draw("psame");

                //drawLatex(0.35, 0.2, Form(" From fitting: p0 = %3.2f #pm %3.2f", flinear[t][p]->GetParameter(0), flinear[t][p]->GetParError(0)), 132, 0.05, 4);
                //drawLatex(0.37, 0.15, Form(" p1 = %3.2f #pm %3.2f", flinear[t][p]->GetParameter(1), flinear[t][p]->GetParError(1)), 132, 0.05, 4);

                c[1]->cd(t*4+p+1);
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.15);
                hErr[t]->Draw();
                f1->Draw("same");
                gPolParComErr[t][p]->Fit(flinearErr[t][p]);
                flinearErr[t][p]->Draw("same");
                gPolParComErr[t][p]->Draw("psame");
                if(t==0)drawLatex(0.3, 0.94, Form(" %2.0f < p_{T} < %2.0f GeV/c", lowPtBins[p], highPtBins[p]), 22, 0.06, 1);
            }}

        double sysPar[2][4];
        double sysParErr[2][4];
        for(int p=0; p<4; p++){
            getCorrFactor(flinear[0][p], flinear[1][p], result[fr][0][p], result[fr][1][p], sysPar[0][p], sysPar[1][p]);
            getCorrFactor(flinearErr[0][p], flinearErr[1][p], statErr[fr][0][p], statErr[fr][1][p], sysParErr[0][p], sysParErr[1][p]);
            //cout<<sysPar[0][p]<<", "<<sysPar[1][p]<<endl;
            //cout<<sysParErr[0][p]<<", "<<sysParErr[1][p]<<endl;
        }

        if(fr==1)ofstream outfile(Form("correct_factor_Cs.txt"));
        else if(fr==0)ofstream outfile(Form("correct_factor_Hx.txt"));
        for(int t=0; t<2; t++){
            outfile<<Form("%s:",typeName[t])<<endl;
            for(int p=0; p<4; p++){
                outfile<< setprecision(3)<<sysPar[t][p] <<" "<< sysParErr[t][p] <<";    ";
            }
            outfile<<endl;
        }
        outfile<<endl;
        outfile<<"Final results:"<<endl;
        for(int t=0; t<3; t++){
            outfile<<Form("%s:",typeName[t])<<endl;
            for(int p=0; p<4; p++){
                if(t==0)outfile<< setprecision(3)<<result[fr][t][p]+sysPar[t][p] <<" "<< statErr[fr][t][p]+sysParErr[t][p] <<";    ";
                else if(t==1)outfile<< setprecision(3)<<result[fr][t][p]+sysPar[t][p] <<" "<< statErr[fr][t][p]+sysParErr[t][p] <<";    ";
                else if(t==2){
                    double inv, invErr;
                    calInvariant(result[fr][0][p]+sysPar[0][p], statErr[fr][0][p]+sysParErr[0][p], result[fr][1][p]+sysPar[1][p], statErr[fr][1][p]+sysParErr[1][p], 0., inv, invErr);
                    outfile<< setprecision(3)<<inv<<" "<< invErr <<";    ";
                }
            }
            outfile<<endl;
        }

        TLegend *leg2 = new TLegend(0.4, 0.2, 0.85, 0.35);
        leg2->SetHeader("Use y = p0 #bullet x + p1");
        leg2->SetTextFont(132);
        leg2->SetTextSize(0.05);
        leg2->SetBorderSize(0);
        leg2->SetFillColor(0);
        leg2->AddEntry(f1, "p0 = 1, p1 = 0", "l");
        leg2->AddEntry(flinear[0][0], "fit Last iteration", "l");
        c[0]->cd(4);
        leg2->Draw("same");

        TLegend *leg = new TLegend(0.4, 0.2, 0.85, 0.32);
        leg->SetTextFont(132);
        leg->SetTextSize(0.05);
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->AddEntry(gPolParCom[0][0], "Last Iteraction", "p");
        c[0]->cd(8);
        leg->Draw("same");

        TLegend *leg3 = new TLegend(0.4, 0.2, 0.85, 0.32);
        leg3->SetTextFont(132);
        leg3->SetTextSize(0.05);
        leg3->SetBorderSize(0);
        leg3->SetFillColor(0);
        leg3->AddEntry(gPolParComErr[0][0], "Last Iteraction", "p");
        c[1]->cd(8);
        leg3->Draw("same");

        for(int i=0; i<2; i++){
            for(t=0; t<2; t++){
                c[i]->cd(1+4*t);
                if(i==0) drawLatex(0.2, 0.6, Form("#lambda_{%s}", polTitle[t]), 22, 0.1, 1);
                else drawLatex(0.2, 0.6, Form("#sigma_{#lambda_{%s}}", polTitle[t]), 22, 0.1, 1);
            }
        }

        c[0]->SaveAs(Form("Toy_CorrFactor_%s.pdf", frameName[fr]));
        c[1]->SaveAs(Form("Toy_CorrFactor_%s_err.pdf", frameName[fr]));
    }
    return 0;
}
