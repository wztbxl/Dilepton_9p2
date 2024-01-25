#include "ana_header.h"
#include "ana_Function.h"

int ana_ite(int nIter=1)
{
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->Divide(4,4);
    
    int nPad = 1;
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            correctEff(f,p,nIter);
            myFit(f, p);
            for(int i=0; i<2; i++){
                c1->cd(nPad);
                hCorr[f][p][i]->Draw("same");
                fLambda[f][p][i]->Draw("same");
                if(nPad%2==0){
                    drawLatex(0.6,0.3,Form("#lambda_{#theta} = %3.3f #pm %2.2f",parLambda[f][p][0],parLambdaErr[f][p][0]),22,0.045,1);
                    drawLatex(0.6,0.23,Form("#lambda_{#phi} = %3.3f #pm %2.2f",parLambda[f][p][1],parLambdaErr[f][1][0]),22,0.045,1);
                    drawLatex(0.6,0.16,Form("#chi^{2}/ndf = %4.1f / %d",parChi2[f][p],parNDF[f][p]),22,0.045,1);
                }
                if(nPad<=16) nPad++;
            }
        }
    }
    
    c1->SaveAs(Form("Iter%d.pdf",nIter));
    
    TFile *fout[2][nPtBins];
    for(int f=0; f<2; f++){
        for(int p=0; p<nPtBins; p++){
            fout[f][p] = new TFile(Form("%sPt%dIter%d.root", frameName[f], p, nIter),"recreate");
            fout[f][p]->cd();
            hIterChi[f][p] = new TH1D(Form("hIterChi%sPt%d", frameName[f], p),Form(""),10, 0, 10);
            hIterChi[f][p]->SetBinContent(nIter, parChi2[f][p]/parNDF[f][p]);
            hIterChi[f][p]->Write();
            hIterCovar[f][p] = new TH1D(Form("hIterCovar%sPt%d", frameName[f], p),Form(""),10, 0, 10);
            hIterCovar[f][p]->SetBinContent(nIter, parCovar[f][p]);
            hIterCovar[f][p]->Write();
            for(int an=0; an<3; an++){
                hIter[f][p][an] = new TH1D(Form("hIter%s%sPt%d", frameName[f], angleName[an], p),Form(""),10, 0, 10);
                hIter[f][p][an]->SetBinContent(nIter, parLambda[f][p][an]);
                hIter[f][p][an]->SetBinError(nIter, parLambdaErr[f][p][an]);
                hIter[f][p][an]->Write();
            }
            fout[f][p]->Close();
        }
    }

    return 0;
}

//---------------------------------
void correctEff(int Frame=0, int nPt=0, int nIterate=1)
{
    gStyle->SetOptStat("ne");
    gStyle->SetOptFit(111);
    
    // get histogram from real data //
    TFile *fData = TFile::Open(Form("./systematic_Signal.root",nPt));
    TH1D *hDataCos = (TH1D*)fData->Get(Form("h%s%sPt%d", frameName[Frame], angleName[0], nPt));
    TH1D *hDataPhi = (TH1D*)fData->Get(Form("h%s%sPt%d", frameName[Frame], angleName[1], nPt));
    
    // get histogram from embedding data //
    TFile *fEff;
    if(nIterate==1) fEff = TFile::Open(Form("efficiency%d.histo.root",nIterate));
    else{
        fEff = TFile::Open(Form("efficiencyFrame%dPt%dIter%d.histo.root", Frame, nPt, nIterate));
    }
    THnSparseF* hMcMPtCosPhi;
    THnSparseF* hRcMPtCosPhi;
    if(Frame==1){
        hMcMPtCosPhi = (THnSparseF*) fEff->Get("mhMcMPtCosPhiCS");
        hRcMPtCosPhi = (THnSparseF*) fEff->Get("mhRcMPtCosPhiCS");
    }
    else{
        hMcMPtCosPhi = (THnSparseF*) fEff->Get("mhMcMPtCosPhi");
        hRcMPtCosPhi = (THnSparseF*) fEff->Get("mhRcMPtCosPhi");
    }
    // for cos_theta
    TH1D *hMcCos = (TH1D*)getEffHist(hMcMPtCosPhi, nPt, 2, "hMcCos");
    TH1D *hRcCos = (TH1D*)getEffHist(hRcMPtCosPhi, nPt, 2, "hRcCos");
    TH1D *hRebinMcCos = (TH1D*)rebHisto(hMcCos, "hRebinMcCos", nCosBinsPt, xCosBinsPt, "Y");
    TH1D *hRebinRcCos = (TH1D*)rebHisto(hRcCos, "hRebinRcCos", nCosBinsPt, xCosBinsPt, "Y");
    // for phi //
    TH1D *hMcPhi = (TH1D*)getEffHist(hMcMPtCosPhi, nPt, 3, "hMcPhi");
    TH1D *hRcPhi = (TH1D*)getEffHist(hRcMPtCosPhi, nPt, 3, "hRcPhi");
    TH1D *hRebinMcPhi = (TH1D*)rebHisto(hMcPhi, "hRebinMcPhi", nPhiBinsPt, xPhiBinsPt, "Y");
    TH1D *hRebinRcPhi = (TH1D*)rebHisto(hRcPhi, "hRebinRcPhi", nPhiBinsPt, xPhiBinsPt, "Y");
    
    TH1D *hEffCos = (TH1D*)divideHisto(hRebinRcCos,hRebinMcCos);
    hEffCos->SetNameTitle(Form("hEffCos"), Form("hEffCos"));
    TH1D *hEffPhi = (TH1D*)divideHisto(hRebinRcPhi,hRebinMcPhi);
    hEffPhi->SetNameTitle(Form("hEffPhi"), Form("hEffPhi"));
    
    //--- normalize both data and efficience ---//
    hDataCos = (TH1D*)normHisto(hDataCos);
    hDataPhi = (TH1D*)normHisto(hDataPhi);
    hEffCos = (TH1D*)normHisto(hEffCos);
    hEffPhi = (TH1D*)normHisto(hEffPhi);
    
    hCorr[Frame][nPt][0] = (TH1D*) divideHisto(hDataCos,hEffCos,1,1,"");
    hCorr[Frame][nPt][0]->SetName(Form("hCorr%s%sPt%d", frameName[Frame], angleName[0], nPt));
    hCorr[Frame][nPt][1] = (TH1D*) divideHisto(hDataPhi,hEffPhi,1,1,"");
    hCorr[Frame][nPt][1]->SetName(Form("hCorr%s%sPt%d", frameName[Frame], angleName[1], nPt));
}

//---------------------------------
TH1D* getEffHist(THnSparseF *hn, int nPt, int axis, const char *name)
{
    THnSparseF *hClone = (THnSparseF*)hn->Clone("hClone");
    hClone->GetAxis(1)->SetRangeUser(xPtBins[nPt], xPtBins[nPt+1]);
    TH1D *h = (TH1D*)hClone->Projection(axis);
    h->SetName(name);
    return h;
}

//---------------------------------
void myFit( int chooseFrame=0, int nPt = 0) {
    //gStyle->SetOptFit();
    gStyle->SetStatY(0.6);
    TH1::SetDefaultSumw2();
    
    // read data
    
    h1 = hCorr[chooseFrame][nPt][0];
    h1->Scale(1./h1->Integral());
    h1->SetMinimum(0);
    h2 = hCorr[chooseFrame][nPt][1];
    h2->Scale(1./h2->Integral());
    h2->SetMinimum(0);
    
    //The default minimizer is Minuit
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter * minuit = TVirtualFitter::Fitter(0,4);
    for (int i = 0; i < 4; ++i) {
        minuit->SetParameter(i, Form(parName[i]), iniPar[i], stp[i], iniVlow[i], iniVhigh[i]);
    }
    //minuit->FixParameter(0);
    //minuit->FixParameter(1);
    //minuit->FixParameter(2);
    //minuit->FixParameter(3);
    minuit->SetFCN(myFcn);
    
    double arglist[10];
    arglist[0] = 0;
    // set print level
    minuit->ExecuteCommand("SET PRINT",arglist,2);
    
    // minimize
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance
    minuit->ExecuteCommand("MIGRAD",arglist,2);
    
    //get result
    double minParams[4];
    double parErrors[4];
    for (int i = 0; i < 4; ++i) {
        minParams[i] = minuit->GetParameter(i);
        parErrors[i] = minuit->GetParError(i);
    }
    double chi2, edm, errdef;
    int nvpar, nparx;
    minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
    int ndf = npfits-nvpar;
    
    fLambda[chooseFrame][nPt][0] = new TF1(Form("f%s%sPt%d", frameName[chooseFrame], angleName[0], nPt),"[0]*(3*(1+[1]*x*x))/(2*([1]+3))",-1,1);
    fLambda[chooseFrame][nPt][0]->SetParameter(0,minParams[0]);fLambda[chooseFrame][nPt][0]->SetParError(0,parErrors[0]);
    fLambda[chooseFrame][nPt][0]->SetParameter(1,minParams[1]);fLambda[chooseFrame][nPt][0]->SetParError(1,parErrors[1]);
    fLambda[chooseFrame][nPt][0]->SetParameter(2,minParams[2]);fLambda[chooseFrame][nPt][0]->SetParError(2,parErrors[2]);
    
    fLambda[chooseFrame][nPt][1] = new TF1(Form("f%s%sPt%d", frameName[chooseFrame], angleName[0], nPt),"[0]*((1+2*[2]*cos(2*x)/(3+[1]))/(pi))",0, pi);
    fLambda[chooseFrame][nPt][1]->SetParameter(0,minParams[3]);fLambda[chooseFrame][nPt][1]->SetParError(0,parErrors[3]);
    fLambda[chooseFrame][nPt][1]->SetParameter(1,minParams[1]);fLambda[chooseFrame][nPt][1]->SetParError(1,parErrors[1]);
    fLambda[chooseFrame][nPt][1]->SetParameter(2,minParams[2]);fLambda[chooseFrame][nPt][1]->SetParError(2,parErrors[2]);
    
    
    //--- save result to rootfile
    TString title;
    if(chooseFrame==1) title = "Collins Soper Frame";
    else if(chooseFrame==0) title = "Helicity Frame";
    TH1D *hFinalTheta = new TH1D("hFinalTheta",Form("%s: #lambda_{#theta} distribution; p_{T} (GeV/c); #lambda_{#theta}", title.Data()), nPtBins, xPtBins);
    TH1D *hFinalPhi = new TH1D("hFinalPhi",Form("%s: #lambda_{#varphi} distribution; p_{T} (GeV/c); #lambda_{#varphi}", title.Data()), nPtBins, xPtBins);
    
    double theta = minParams[1];
    double thetaErr = parErrors[1];
    
    double phi = minParams[2];
    double phiErr = parErrors[2];
    
    TH1D *hLambdaTheta = new TH1D("hLambdaTheta","; p_{T} GeV/c; #lambda_{#theta}", nPtBins, xPtBins);
    hLambdaTheta->SetBinContent(nPt+1, theta);
    hLambdaTheta->SetBinError(nPt+1, thetaErr);
    hLambdaTheta->SetMarkerStyle(21);
    hLambdaTheta->SetMarkerColor(2);
    hLambdaTheta->SetLineColor(2);
    TH1D *hLambdaPhi = new TH1D("hLambdaPhi","; p_{T} GeV/c; #lambda_{#varphi}", nPtBins, xPtBins);
    hLambdaPhi->SetBinContent(nPt+1, phi);
    hLambdaPhi->SetBinError(nPt+1, phiErr);
    hLambdaPhi->SetMarkerStyle(21);
    hLambdaPhi->SetMarkerColor(2);
    hLambdaPhi->SetLineColor(2);
    
    cout<<"pTheta: "<< theta <<" ; pPhi: "<< phi <<endl;
    
    double inv = 0;
    double invErr = 0;
    
    //--- get error matrix ---//
    Double_t matrix[4][4];
    gMinuit->mnemat(&matrix[0][0],4);
    
    cout<<"********** error matrix ***********"<<endl;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            cout<<matrix[i][j]<<";  ";
        }
        cout<<endl;
    }
    cout<<"***********************************"<<endl;
    
    cout<<matrix[1][2]<<endl;
    cout<<matrix[2][1]<<endl;
    
    calInvariant( theta, thetaErr, phi, phiErr, matrix[1][2], inv, invErr);
    
    TH1D *hLambdaInv = new TH1D("hLambdaInv","; p_{T} GeV/c; #lambda_{#varphi}", nPtBins, xPtBins);
    hLambdaInv->SetBinContent(nPt+1, inv);
    hLambdaInv->SetBinError(nPt+1, invErr);
    hLambdaInv->SetMarkerStyle(21);
    hLambdaInv->SetMarkerColor(2);
    hLambdaInv->SetLineColor(2);
    
    cout<<"***********************************"<<endl;
    cout<<"Theta: "<< theta <<" ; "<< thetaErr <<endl;
    cout<<"Phi:   "<< phi <<" ; "<< phiErr <<endl;
    
    parLambda[chooseFrame][nPt][0] = theta;
    parLambdaErr[chooseFrame][nPt][0] = thetaErr;
    parLambda[chooseFrame][nPt][1] = phi;
    parLambdaErr[chooseFrame][nPt][1] = phiErr;
    
    cout<<"***********************************"<<endl;
    cout<<"inv:     "<< inv <<endl;
    cout<<"inv err: "<< invErr <<endl;
    
    parLambda[chooseFrame][nPt][2] = inv;
    parLambdaErr[chooseFrame][nPt][2] = invErr;
    
    parChi2[chooseFrame][nPt] = chi2;
    parNDF[chooseFrame][nPt] = ndf;
    parCovar[chooseFrame][nPt] = matrix[1][2];
    
    return ;
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
