const char *frameName[2] = {"HX", "CS"};
const char* typeName[3] = {"Theta", "Phi", "inv"};
const char* polTitle[3] = {"#theta", "#varphi", "inv"};

//-----------------------------------------
void SetHistoMarker(TH2 *h, int marker, int color)
{
    if(!h) return;
    h->SetMarkerStyle(marker);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
}

//-----------------------------------------
void SetGraphMarker(TGraph *g, int marker, int color)
{
    if(!g) return;
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
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
void getCorrFactor(TF1 *fTheta, TF1 *fPhi, double p1, double p2, double &fac1, double &fac2)
{
    //-- theta
    double a1, b1;
    a1 = fTheta->GetParameter(0);
    b1 = fTheta->GetParameter(1);
    //cout<<a1<<", "<<b1<<endl;
    //cout<<(p1-b1)/a1<<endl;
    fac1  = (p1-b1)/a1 - p1;
    //-- phi
    double a2, b2;
    a2 = fPhi->GetParameter(0);
    b2 = fPhi->GetParameter(1);
    //cout<<a2<<", "<<b2<<endl;
    //cout<<(p2-b2)/a2<<endl;
    fac2  = (p2-b2)/a2 - p2;
}

//__________________________________________________
void getErrCorrFactor(TF1 *fTheta, TF1 *fPhi, double p1, double p2, double &fac1, double &fac2)
{
    //-- theta
    double a1;
    a1 = fTheta->GetParameter(0);
    fac1  = p1/a1 - p1;
    //-- phi
    double a2;
    a2 = fPhi->GetParameter(0);
    fac2  = p2/a2 - p2;
}

//-------------------------------------------------
double calInvariant(double theta, double phi)
{
    double par = (theta + 3.* phi)/(1 - phi);
    return par;
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
