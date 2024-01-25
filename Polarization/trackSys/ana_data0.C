#include "./constant.h"

int ana_data0(double dca=varDca[0], double nHitsFit=varNhitsFit[0], double nHitsDedx=varNhitsDedx[0], double nSigmaPi=varNsigmaPi[1])
{
    char buf[1024] = Form("_%1.0f_%2.0f_%2.0f_%2.1f",dca, nHitsFit, nHitsDedx, nSigmaPi);
    // pol4 //
    get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    drawHistrom(buf,4,4,0,1,1); //rebin & nPol & sLR & fix & savePDF
    closeFile();
    get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    drawHistrom(buf,4,4,1,1,1); //rebin & nPol & sLR & fix & savePDF
    closeFile();
    get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    drawHistrom(buf,4,4,2,1,1); //rebin & nPol & sLR & fix & savePDF
    closeFile();
    get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    drawHistrom(buf,4,4,3,1,1); //rebin & nPol & sLR & fix & savePDF
    closeFile();
    get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    drawHistrom(buf,4,4,4,1,1); //rebin & nPol & sLR & fix & savePDF
    closeFile();

    get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    drawHistrom(buf,2,4,0,1,1); //rebin & nPol & sLR & fix & savePDF
    closeFile();//rebin2&pol4

    // pol2 //
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,2,0,0,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,2,1,0,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,2,2,0,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,2,3,0,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,2,4,0,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();

    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,2,0,1,1);
    //closeFile();//fix&pol2
    
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,2,2,0,0,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();//rebin8&pol4

    //// pol3 //
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,3,0,1,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();
    //
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,3,0,0,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();
    //
    //// pol5 //
    //get_Angular_Histo(dca, nHitsFit, nHitsDedx, nSigmaPi);
    //drawHistrom(buf,4,5,0,1,1); //rebin & nPol & sLR & fix & savePDF
    //closeFile();

    return 0;
}

//------------- get invariant mass distribution -------------------------------------------//
void get_Angular_Histo(double dca, double nHitsFit, double nHitsDedx, double nSigmaPi)
{
    fin = TFile::Open(Form("./Data_%1.0f_%2.0f_%2.0f_%2.1f/Data_%1.0f_%2.0f_%2.0f_%2.1f.root", dca, nHitsFit, nHitsDedx, nSigmaPi, dca, nHitsFit, nHitsDedx, nSigmaPi),"read");

    //-- Get histogram --//
    for(int s=0; s<2; s++){
        hnMPtCosPhi[0][s] = (THnSparseF*)fin->Get(Form("mh%sMPtCosPhi",signName[s]));
        hnMPtCosPhi[1][s] = (THnSparseF*)fin->Get(Form("mh%sMPtCosPhiCS",signName[s]));
    }

    // for cos //
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int s=0; s<2; s++){
                hInvMass[f][s][p][0] = new TH1D[nCosBinsPt];
                for(int an=0; an<nCosBinsPt; an++){
                    if(an+1<=cutBin[f][p][0] || an+1 >= nCosBinsPt+1-cutBin[f][p][0]) continue;
                    hInvMass[f][s][p][0][an]=(TH1D*)getInvMass(hnMPtCosPhi[f][s], 2, p, an, Form("h%s%sPt%dCos%d",frameName[f],signName[s],p,an));
                }
            }
        }
    }
    // for phi //
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int s=0; s<2; s++){
                hInvMass[f][s][p][1] = new TH1D[nPhiBinsPt];
                for(int an=0; an<nPhiBinsPt; an++){
                    if(an+1<=cutBin[f][p][1] || an+1 >= nPhiBinsPt+1-cutBin[f][p][1]) continue;
                    hInvMass[f][s][p][1][an]=(TH1D*)getInvMass(hnMPtCosPhi[f][s], 3, p, an, Form("h%s%sPt%dPhi%d",frameName[f],signName[s],p,an));
                }
            }
        }
    }
    return ;
}
//--------------------------------------------------------
TH1D* getInvMass(THnSparseF* hMPtCosPhi, int axis, int nPt, int nAn, const char *name)
{
    THnSparseF *hClone = (THnSparseF*) hMPtCosPhi->Clone("hClone");
    hClone->GetAxis(0)->SetRangeUser(2.2, 4.2);
    hClone->GetAxis(1)->SetRangeUser(xPtBins[nPt], xPtBins[nPt+1]);
    if(axis==2)hClone->GetAxis(axis)->SetRangeUser(xCosBinsPt[nAn], xCosBinsPt[nAn+1]);
    else if(axis==3)hClone->GetAxis(axis)->SetRange(nAn+1, nAn+1);
    TH1D* h1 = (TH1D*) hClone->Projection(0);
    h1->SetName(name);
    return h1;
}

//------------- fit histograme -------------------------------------------//
TF1 *fitHisto(
        TH1D *h1,
        TH1D *h2,
        double *matrixPram,
        int rebin=4,
        int nPol=4,
        int nPt=0,
        int angle=0,
        int nSub=0,
        double rangeLow=fitCosRangeLow[0][0][0],
        double rangeUp=fitCosRangeUp[0][0][0],
        bool fix=true,
        const int nLoop=4)
{
    TH1D *hUL = (TH1D*)h1->Clone(Form("%s_clone",h1->GetName()));
    TH1D *hLS = (TH1D*)h2->Clone(Form("%s_clone",h1->GetName()));
    TH1D *hErr = new TH1D("hErr", "hErr", 40, 0, 40);

    hUL->Rebin(rebin);
    hLS->Rebin(rebin);

    TF1* fLS;
    if(nPol==2)
        fLS = new TF1("fLS", "pol2", rangeLow, rangeUp);
    else if(nPol==3)
        fLS = new TF1("fLS", "pol3", rangeLow, rangeUp);
    else if(nPol==4)
        fLS = new TF1("fLS", "pol4", rangeLow, rangeUp);
    else if(nPol==5)
        fLS = new TF1("fLS", "pol5", rangeLow, rangeUp);
    else
        cout<<"**** fitUL: fLS ERR only Pol2-Pol5"<<endl;

    //-- fit LS first --//
    fLS->SetNpx(10000);
    hLS->Fit(fLS,"R");
    double *pLS = fLS->GetParameters();

    //-- Loop fit UL --//
    const char *angleName[2] = {"Cos","Phi"};
    TF1 *fUL[nLoop];
    for(int i=0; i<nLoop; i++){
        if(nPol==2)
            fUL[i] = new TF1(Form("fUL_Pt%d_%s%d_Pol%d_Loop%d",nPt,angleName[angle], nSub, nPol,i), flinePol2, rangeLow, rangeUp, 3+nPol+1);
        else if(nPol==3)
            fUL[i] = new TF1(Form("fUL_Pt%d_%s%d_Pol%d_Loop%d",nPt,angleName[angle], nSub, nPol,i), flinePol3, rangeLow, rangeUp, 3+nPol+1);
        else if(nPol==4)
            fUL[i] = new TF1(Form("fUL_Pt%d_%s%d_Pol%d_Loop%d",nPt,angleName[angle], nSub, nPol,i), flinePol4, rangeLow, rangeUp, 3+nPol+1);
        else if(nPol==5)
            fUL[i] = new TF1(Form("fUL_Pt%d_%s%d_Pol%d_Loop%d",nPt,angleName[angle], nSub, nPol,i), flinePol5, rangeLow, rangeUp, 3+nPol+1);
        else
            cout<<"**** fitUL: fUL ERR only Pol2-Pol5"<<endl;

        double mean=3.1;
        double sigma=sigmaPt[nPt];
        if(fix){
            fUL[i] -> FixParameter(1,mean);
            fUL[i] -> FixParameter(2,sigma);
        }
        fUL[i]->SetLineColor(2);
    }//book fUL
    if(nPol==2)
        fUL[0] -> SetParameters(100, mean, sigma, pLS[0], pLS[1], pLS[2]);
    else if(nPol==3)
        fUL[0] -> SetParameters(100, mean, sigma, pLS[0], pLS[1], pLS[2], pLS[3]);
    else if(nPol==4)
        fUL[0] -> SetParameters(100, mean, sigma, pLS[0], pLS[1], pLS[2], pLS[3], pLS[4]);
    else if(nPol==5)
        fUL[0] -> SetParameters(100, mean, sigma, pLS[0], pLS[1], pLS[2], pLS[3], pLS[4], pLS[5]);
    else
        cout<<"**** fitUL: setParErr ERR only Pol2-Pol5"<<endl;
    fUL[0] -> FixParameter(1,mean);
    fUL[0] -> SetParLimits(0,0,10000);
    fUL[0] -> SetParLimits(2,0,1);
    fUL[0] -> FixParameter(2,sigma);

    TFitResultPtr r[nLoop];
    for(int i=0; i<nLoop; i++){
        fUL[i]->SetNpx(10000);
        r[i] = hUL->Fit(fUL[i],"RLVSN+","same");
        double *p = fUL[i] -> GetParameters();
        if(i < nLoop-1){
            fUL[i+1]->SetParameters(p);
        }
    }
    cout<<"TFitResultPtr:Status: "<< r[nLoop-1] -> Status() <<endl;
    cout<<"TFitResultPtr:CovMatrixStatus: "<< r[nLoop-1] -> CovMatrixStatus() <<endl;

    Double_t *ppMat =r[nLoop-1]->GetCovarianceMatrix().GetMatrixArray();

    if(nPol==2){
        for(int i=0; i<3; i++){
            matrixPram[0+i*3] = ppMat[21+i*6];
            matrixPram[1+i*3] = ppMat[22+i*6];
            matrixPram[2+i*3] = ppMat[23+i*6];
        }
    }
    else if(nPol==3){
        for(int i=0; i<4; i++){
            matrixPram[0+i*4] = ppMat[24+i*7];
            matrixPram[1+i*4] = ppMat[25+i*7];
            matrixPram[2+i*4] = ppMat[26+i*7];
            matrixPram[3+i*4] = ppMat[27+i*7];
        }
    }
    else if(nPol==4){
        for(int i=0; i<5; i++){
            matrixPram[0+i*5] = ppMat[27+i*8];
            matrixPram[1+i*5] = ppMat[28+i*8];
            matrixPram[2+i*5] = ppMat[29+i*8];
            matrixPram[3+i*5] = ppMat[30+i*8];
            matrixPram[4+i*5] = ppMat[31+i*8];
        }
    }
    else if(nPol==5){
        for(int i=0; i<6; i++){
            matrixPram[0+i*6] = ppMat[30+i*9];
            matrixPram[1+i*6] = ppMat[31+i*9];
            matrixPram[2+i*6] = ppMat[32+i*9];
            matrixPram[3+i*6] = ppMat[33+i*9];
            matrixPram[4+i*6] = ppMat[34+i*9];
            matrixPram[5+i*6] = ppMat[35+i*9];
        }
    }
    else
        cout<<"**** fitUL: ERR only Pol2-Pol5"<<endl;

    return fUL[nLoop-1];
}

//--------------------------------------
double flinePol2( double *x, double *par)
{
    if(reject && x[0] > 3.6 && x[0] < 3.8 ){
        TF1::RejectPoint();
        return 0;
    }
    return (par[0]/TMath::Sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-TMath::Power((x[0]-par[1])/TMath::Sqrt(2)/par[2],2))
            + par[3] + par[4] * x[0] + par[5] * x[0] * x[0]);
}

//--------------------------------------
double flinePol3( double *x, double *par)
{
    if(reject && x[0] > 3.6 && x[0] < 3.8 ){
        TF1::RejectPoint();
        return 0;
    }
    return (par[0]/TMath::Sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-TMath::Power((x[0]-par[1])/TMath::Sqrt(2)/par[2],2))
            + par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0]);
}

//--------------------------------------
double flinePol4( double *x, double *par)
{
    if(reject && x[0] > 3.6 && x[0] < 3.8 ){
        TF1::RejectPoint();
        return 0;
    }
    return (par[0]/TMath::Sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-TMath::Power((x[0]-par[1])/TMath::Sqrt(2)/par[2],2))
            + par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0]* x[0]);
}

//--------------------------------------
double flinePol5( double *x, double *par)
{
    if(reject && x[0] > 3.6 && x[0] < 3.8 ){
        TF1::RejectPoint();
        return 0;
    }
    return (par[0]/TMath::Sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-TMath::Power((x[0]-par[1])/TMath::Sqrt(2)/par[2],2))
            + par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0]* x[0] + par[8] * x[0] * x[0] * x[0]* x[0]* x[0]);
}


//------------- draw on canvas -------------------------------------------//
void drawHistrom(
        char *outFile,
        int rebin=4,
        int nPol=4,
        int sLR=0,//shift[sLR] {0, 0.1, -0.1, 0.05, -0.05}
        int fix=1,
        const int savePDF = 1)
{
    gStyle->SetOptStat(0);
    gStyle->SetGridWidth(0);
    gStyle->SetCanvasColor(10);

    // book canvas //
    TCanvas *cInv[2][nPtBins][2]; //frame & pt & angle
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            for(int an=0; an<2; an++){
                cInv[f][p][an] = new TCanvas(Form("cInv%sPt%d%s", frameName[f], p, angleName[an]), Form("cInv%sPt%d%s", frameName[f], p, angleName[an]), 800, 600);
                if(an==0) cInv[f][p][an]->Divide(4,3);
                else if(an==1) cInv[f][p][an]->Divide(4,4);
            }
        }
    }

    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            // for cos_theta distribution //
            hCosFit[f][p] = new TH1D(Form("h%sPt%dCosFitRebin%dPol%dShift%dfix%d", frameName[f], p, rebin, nPol, sLR, fix),Form(";cos#theta;# of J/#psi"), nCosBinsPt, xCosBinsPt);
            hCosCount[f][p] = new TH1D(Form("h%sPt%dCosCountRebin%dPol%dShift%dfix%d", frameName[f], p, rebin, nPol, sLR, fix),Form(";cos#theta;# of J/#psi"), nCosBinsPt, xCosBinsPt);

            for(int an=0; an<nCosBinsPt; an++){
                if(an+1<=cutBin[f][p][0] || an+1 >= nCosBinsPt+1-cutBin[f][p][0]) continue;
                if(hInvMass[f][0][p][0][an]){
                    cInv[f][p][0]->cd(an+1);

                    matrixCosPram[f][p][an] = new double[(nPol+1)*(nPol+1)];
                    TF1 *fULMCos = (TF1*)fitHisto(hInvMass[f][0][p][0][an], hInvMass[f][1][p][0][an], matrixCosPram[f][p][an], rebin, nPol, p, 0, an, fitCosRangeLow[f][p][an]+shift[sLR], fitCosRangeUp[f][p][an]+shift[sLR], fix);

                    double backBoundary = hInvMass[f][0][p][0][an]->GetBinCenter( hInvMass[f][0][p][0][an]->FindFirstBinAbove(0) );
                    double mLowLimit = fULMCos->GetParameter(1) - 3*fULMCos->GetParameter(2);
                    double mHighLimit = fULMCos->GetParameter(1) + 3*fULMCos->GetParameter(2);
                    if(backBoundary>massLowPt[p]) {
                        if(backBoundary<2.9 && massLowPt[p]<2.9) mLowLimit = 2.9;
                        else if(backBoundary<3.0 && massLowPt[p]<3.0) mLowLimit = 3.0;
                        else if(backBoundary>3.0 && backBoundary<3.1 && massLowPt[p]<3.0) mLowLimit = 3.0;
                    }
                    double tTotal_count = hInvMass[f][0][p][0][an]->Integral(hInvMass[f][0][p][0][an]->GetXaxis()->FindBin(mLowLimit+1e-6), hInvMass[f][0][p][0][an]->GetXaxis()->FindBin(mHighLimit-1e-6));

                    hInvMass[f][0][p][0][an]->Rebin(rebin);
                    hInvMass[f][0][p][0][an]->Draw();
                    hInvMass[f][1][p][0][an]->Rebin(rebin);
                    hInvMass[f][1][p][0][an]->Draw("samehist");
                    fULMCos->Draw("same");
                    TF1 *fSigCos = (TF1*)getFsig(fULMCos, fitCosRangeLow[f][p][an]+shift[sLR], fitCosRangeUp[f][p][an]+shift[sLR]);
                    fSigCos->Draw("same");
                    TF1 *fBkgCos = (TF1*)getFbkg(fULMCos, nPol, fitCosRangeLow[f][p][an]+shift[sLR], fitCosRangeUp[f][p][an]+shift[sLR]);
                    fBkgCos->Draw("same");
                    TLegend *legFit = new TLegend(0.65,0.7,0.89,0.89);
                    legFit->SetBorderSize(0);
                    legFit->SetTextSize(0.04);
                    legFit->AddEntry(hInvMass[f][0][p][0][an],"unlike-sign pairs","lp");
                    legFit->AddEntry(hInvMass[f][1][p][0][an],"like-sign pairs","l");
                    legFit->AddEntry(fSigCos,"signal","f");
                    legFit->AddEntry(fBkgCos,"background","l");
                    legFit->Draw("same");

                    double *par = fULMCos->GetParameters();
                    double *parErr = fULMCos->GetParErrors();
                    Double_t tChi = fULMCos->GetChisquare();
                    Double_t tNdf = fULMCos->GetNDF();
                    double binWidth = hInvMass[f][0][p][0][an]->GetBinWidth(1);
                    Double_t tTotal = (fULMCos->Integral(par[1]-5*par[2], par[1]+5*par[2]))/binWidth;
                    Double_t tSig = (fSigCos->Integral(2,4))/binWidth;
                    Double_t tBkg = tTotal - tSig;
                    Double_t tSignificance = tSig/sqrt(tSig+tBkg);

                    Double_t tBkg_count = fBkgCos -> Integral(mLowLimit, mHighLimit)/binWidth;
                    Double_t *tBkg_count_par = fBkgCos->GetParameters();
                    Double_t tBkg_count_err = (fBkgCos->IntegralError(mLowLimit, mHighLimit, tBkg_count_par, matrixCosPram[f][p][an]))/binWidth;
                    Double_t tSig_count = tTotal_count - tBkg_count;
                    Double_t tSig_count_err = sqrt(tTotal_count + pow(tBkg_count_err,2));

                    drawLatex(0.15,0.85,Form("#chi^{2}/ndf= %5.2f/%2.0f",tChi,tNdf),22,0.045,1);
                    drawLatex(0.15,0.8,Form("N_{J/#psi}= %5.0f#pm%3.0f",par[0]/binWidth,parErr[0]/binWidth),22,0.045,1);
                    drawLatex(0.15,0.75,Form("Mean= %3.2f#pm%1.4f",par[1],parErr[1]),22,0.045,1);
                    drawLatex(0.15,0.7,Form("#sigma= %1.2f#pm%1.2f",par[2],parErr[2]),22,0.045,1);
                    drawLatex(0.15,0.65,Form("S/B= %2.2f",tSig/tBkg),22,0.045,1);
                    drawLatex(0.15,0.6,Form("S/#sqrt{S+B}= %2.2f",tSignificance),22,0.045,1);
                    drawLatex(0.15,0.55,Form("N_{J/#psi counting}= %5.0f#pm%3.0f",tSig_count,tSig_count_err),22,0.045,1);
                    drawLatex(0.15,0.5,Form("bin width: %1.2f",binWidth),22,0.045,1);
                    
                    Double_t tLineHigh = hInvMass[f][0][p][0][an]->GetBinContent(hInvMass[f][0][p][0][an]->GetMaximumBin());
                    drawLine(mLowLimit,0,mLowLimit,tLineHigh,2,2,8);
                    drawLine(mHighLimit,0,mHighLimit,tLineHigh,2,2,8);
                    
                    //-- reject data points with significance < 3 and bad fitting (nan)
                    if(tSignificance!=tSignificance) continue;
                    if(tSignificance<3) continue;

                    hCosFit[f][p]->SetBinContent(an+1, par[0]/binWidth);
                    hCosFit[f][p]->SetBinError(an+1, parErr[0]/binWidth);
                    hCosCount[f][p]->SetBinContent(an+1, tSig_count);
                    hCosCount[f][p]->SetBinError(an+1, tSig_count_err);
                }
            }

            // for phi distribution //
            hPhiFit[f][p] = new TH1D(Form("h%sPt%dPhiFitRebin%dPol%dShift%dfix%d", frameName[f], p, rebin, nPol, sLR, fix),Form(";#varphi;# of J/#psi"), nPhiBinsPt, xPhiBinsPt);
            hPhiCount[f][p] = new TH1D(Form("h%sPt%dPhiCountRebin%dPol%dShift%dfix%d", frameName[f], p, rebin, nPol, sLR, fix),Form(";#varphi;# of J/#psi"), nPhiBinsPt, xPhiBinsPt);

            for(int an=0; an<nPhiBinsPt; an++){
                if(an+1<=cutBin[f][p][1] || an+1 >= nPhiBinsPt+1-cutBin[f][p][1]) continue;
                if(hInvMass[f][0][p][1][an]){
                    cInv[f][p][1]->cd(an+1);

                    matrixPhiPram[f][p][an] = new double[(nPol+1)*(nPol+1)];
                    TF1 *fULMPhi = (TF1*)fitHisto(hInvMass[f][0][p][1][an], hInvMass[f][1][p][1][an], matrixPhiPram[f][p][an], rebin, nPol, p, 1, an, fitPhiRangeLow[f][p][an]+shift[sLR], fitPhiRangeUp[f][p][an]+shift[sLR], fix);

                    double backBoundary = hInvMass[f][0][p][1][an]->GetBinCenter( hInvMass[f][0][p][1][an]->FindFirstBinAbove(0) );
                    //double mLowLimit = massLowPt[p];
                    double mLowLimit = fULMPhi->GetParameter(1) - 3*fULMPhi->GetParameter(2);
                    double mHighLimit = fULMPhi->GetParameter(1) + 3*fULMPhi->GetParameter(2);
                    if(backBoundary>massLowPt[p]) {
                        if(backBoundary<2.9 && massLowPt[p]<2.9) mLowLimit = 2.9;
                        else if(backBoundary<3.0 && massLowPt[p]<3.0) mLowLimit = 3.0;
                        else if(backBoundary>3.0 && backBoundary<3.1 && massLowPt[p]<3.0) mLowLimit = 3.0;
                    }
                    double tTotal_count = hInvMass[f][0][p][1][an]->Integral(hInvMass[f][0][p][1][an]->GetXaxis()->FindBin(mLowLimit+1e-6), hInvMass[f][0][p][1][an]->GetXaxis()->FindBin(mHighLimit-1e-6));

                    hInvMass[f][0][p][1][an]->Rebin(rebin);
                    hInvMass[f][0][p][1][an]->Draw();
                    hInvMass[f][1][p][1][an]->Rebin(rebin);
                    hInvMass[f][1][p][1][an]->Draw("samehist");
                    fULMPhi->Draw("same");
                    TF1 *fSigPhi = (TF1*)getFsig(fULMPhi, fitPhiRangeLow[f][p][an]+shift[sLR], fitPhiRangeUp[f][p][an]+shift[sLR]);
                    fSigPhi->Draw("same");
                    TF1 *fBkgPhi = (TF1*)getFbkg(fULMPhi, nPol, fitPhiRangeLow[f][p][an]+shift[sLR], fitPhiRangeUp[f][p][an]+shift[sLR]);
                    fBkgPhi->Draw("same");

                    TLegend *legFit = new TLegend(0.65,0.7,0.89,0.89);
                    legFit->SetBorderSize(0);
                    legFit->SetTextSize(0.04);
                    legFit->AddEntry(hInvMass[f][0][p][1][an],"unlike-sign pairs","lp");
                    legFit->AddEntry(hInvMass[f][1][p][1][an],"like-sign pairs","l");
                    legFit->AddEntry(fSigPhi,"signal","f");
                    legFit->AddEntry(fBkgPhi,"background","l");
                    legFit->Draw("same");

                    double *par = fULMPhi->GetParameters();
                    double *parErr = fULMPhi->GetParErrors();
                    Double_t tChi = fULMPhi->GetChisquare();
                    Double_t tNdf = fULMPhi->GetNDF();
                    double binWidth = hInvMass[f][0][p][1][an]->GetBinWidth(1);
                    Double_t tTotal = (fULMPhi->Integral(par[1]-5*par[2], par[1]+5*par[2]))/binWidth;
                    Double_t tSig = (fSigPhi->Integral(2,4))/binWidth;
                    Double_t tBkg = tTotal - tSig;
                    Double_t tSignificance = tSig/sqrt(tSig+tBkg);

                    Double_t tBkg_count = fBkgPhi -> Integral(mLowLimit,mHighLimit)/binWidth;
                    Double_t *tBkg_count_par = fBkgPhi->GetParameters();
                    Double_t tBkg_count_err = (fBkgPhi->IntegralError(mLowLimit, mHighLimit, tBkg_count_par, matrixPhiPram[f][p][an]))/binWidth;
                    Double_t tSig_count = tTotal_count - tBkg_count;
                    Double_t tSig_count_err = sqrt(tTotal_count + pow(tBkg_count_err,2));

                    drawLatex(0.15,0.85,Form("#chi^{2}/ndf= %5.2f/%2.0f",tChi,tNdf),22,0.045,1);
                    drawLatex(0.15,0.8,Form("N_{J/#psi}= %5.0f#pm%3.0f",par[0]/binWidth,parErr[0]/binWidth),22,0.045,1);
                    drawLatex(0.15,0.75,Form("Mean= %3.2f#pm%1.4f",par[1],parErr[1]),22,0.045,1);
                    drawLatex(0.15,0.7,Form("#sigma= %1.2f#pm%1.2f",par[2],parErr[2]),22,0.045,1);
                    drawLatex(0.15,0.65,Form("S/B= %2.2f",tSig/tBkg),22,0.045,1);
                    drawLatex(0.15,0.6,Form("S/#sqrt{S+B}= %2.2f",tSignificance),22,0.045,1);
                    drawLatex(0.15,0.55,Form("N_{J/#psi counting}= %5.0f#pm%3.0f",tSig_count,tSig_count_err),22,0.045,1);
                    drawLatex(0.15,0.5,Form("bin width: %1.2f",binWidth),22,0.045,1);
                    
                    Double_t tLineHigh = hInvMass[f][0][p][1][an]->GetBinContent(hInvMass[f][0][p][1][an]->GetMaximumBin());
                    drawLine(mLowLimit,0,mLowLimit,tLineHigh,2,2,8);
                    drawLine(mHighLimit, 0, mHighLimit, tLineHigh, 2, 2, 8);
                    
                    //-- reject data points with significance < 3 and bad fitting (nan)
                    if(tSignificance!=tSignificance) continue;
                    if(tSignificance<3) continue;

                    hPhiFit[f][p]->SetBinContent(an+1, par[0]/binWidth);
                    hPhiFit[f][p]->SetBinError(an+1, parErr[0]/binWidth);
                    hPhiCount[f][p]->SetBinContent(an+1, tSig_count);
                    hPhiCount[f][p]->SetBinError(an+1, tSig_count_err);
                }
            }
        }
    }

    // save pdf //
    if(savePDF){
        char pdfName[1024] = Form("Data%s/Rebin%dPol%dShift%dfix%d", outFile, rebin, nPol, sLR, fix);
        for(int f=0; f<2; f++){
            for(int p=0; p<nPtBins; p++){
                for(int an=0; an<2; an++){
                    if(f==0&&p==0&&an==0) cInv[f][p][an]->Print(Form("%s.pdf(",pdfName),Form("Title: %sPt%d%s", frameName[f], p, angleName[an]));
                    else if(f==1&&p==3&&an==1) cInv[f][p][an]->Print(Form("%s.pdf)",pdfName),Form("Title: %sPt%d%s", frameName[f], p, angleName[an]));
                    else cInv[f][p][an]->Print(Form("%s.pdf",pdfName),Form("Title: %sPt%d%s", frameName[f], p, angleName[an]));
                }
            }
        }
    }

    // save final histo //
    fout = new TFile(Form("%s.root", pdfName),"recreate");
    fout->cd();
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            hCosFit[f][p]->Write();
            hCosCount[f][p]->Write();
            hPhiFit[f][p]->Write();
            hPhiCount[f][p]->Write();
        }
    }

    return ;
}

//-----------------------------------------
TF1* getFsig(TF1* fUL,double rangeLow, double rangeUp)
{
    double *par = fUL->GetParameters();
    TF1 *fSig = new TF1("fSig","[0]/sqrt(2*TMath::Pi())/[2]*exp(-pow((x-[1])/sqrt(2)/[2],2))",rangeLow, rangeUp);
    fSig->SetParameters(par[0],par[1],par[2]);
    fSig->SetLineColor(4);
    fSig->SetFillStyle(3004);
    fSig->SetFillColor(4);
    return fSig;
}

//-----------------------------------------
TF1* getFbkg(TF1* fUL, int nPol, double rangeLow, double rangeUp)
{
    double *par = fUL->GetParameters();
    TF1 *fBkg;
    if(nPol==2){
        fBkg = new TF1("fBkg", "pol2(0)", rangeLow, rangeUp);
        fBkg->SetParameters(par[3],par[4],par[5]);
    }
    else if(nPol==3){
        fBkg = new TF1("fBkg", "pol3(0)", rangeLow, rangeUp);
        fBkg->SetParameters(par[3],par[4],par[5], par[6]);
    }
    else if(nPol==4){
        fBkg = new TF1("fBkg", "pol4(0)", rangeLow, rangeUp);
        fBkg->SetParameters(par[3],par[4],par[5], par[6], par[7]);
    }
    else if(nPol==5){
        fBkg = new TF1("fBkg", "pol5(0)", rangeLow, rangeUp);
        fBkg->SetParameters(par[3],par[4],par[5], par[6], par[7], par[8]);
    }
    else{
        cout<<"** getFbkg ERR: nPol=[2,5]**"<<endl;
    }
    fBkg->SetLineColor(8);
    return fBkg;
}

//----------------------------------------
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

//----------------------------------------
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
    TLine *l1 = new TLine(xlow,ylow,xup,yup);
    l1->SetLineWidth(lineWidth);
    l1->SetLineColor(lineColor);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
    return l1;
}



//-----------------------------------------//
void closeFile()
{
    fin->Close();
    fout->Close();
    return ;
}




