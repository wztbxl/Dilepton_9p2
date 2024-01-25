const char* run_type = "Run15_pp200";
const bool gUseMean = 0;
const int gIter = 10;
const int nExpr = 500;

// Collins-Soper frame //
const int nthetaphi = 4;
const double ptheta[nthetaphi] = {-0.26, 0.48, 0.63, -0.047};
const double pphi[nthetaphi] = {-0.034, 0.008, -0.124, -0.3};
const double gxbins[nthetaphi] = {0.5, 1.5, 2.5, 3.5};

const int color[20] = {1, 2, 4, 6, 8, 1, 2, 4, 6, 8, 1, 2, 4, 6, 8, 1, 2, 4, 6, 8};
const double pthephi = 0.0;
const int nPtBins = 4;
const double lowPtBins[nPtBins] = {0,1,2,4};
const double highPtBins[nPtBins] = {1,2,4,10};
const char* hname[2] = {"Mc", "Rc"};
const char* typeName[3] = {"Theta", "Phi", "inv"};
const char* typeTitle[3] = {"cos#theta","#varphi", "inv"};
const char* polTitle[3] = {"#theta","#varphi", "inv"};
const double PI = TMath::Pi();

const int infColor[nPtBins] = {1, 2, 4, 6};
const int infMarker[nPtBins] = {21, 22, 23, 24};

//=================================================
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
