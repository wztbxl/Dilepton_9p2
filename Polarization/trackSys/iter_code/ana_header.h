const double pi = TMath::Pi();

const int nPtBins = 4;
double xPtBins[nPtBins+1] ={0, 1, 2, 4, 10};

const int nCosBinsPt = 10;
double xCosBinsPt[nCosBinsPt+1] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0 };

const int nPhiBinsPt = 15;
double xPhiBinsPt[nPhiBinsPt+1] = {
    0,1.*pi/15, 2.*pi/15, 3.*pi/15, 4.*pi/15, pi/3., pi*6./15., pi*7./15., pi*8./15., pi*9./15., pi*2./3., pi*11./15., pi*12./15., 13.*pi/15, 14.*pi/15, pi };

const char *frameName[2] = {"HX", "CS"};
const char *angleName[3] = {"Cos", "Phi", "inv"};
const char *labelName[3] = {"#lambda_{#theta}", "#lambda_{#varphi}", "#lambda_{inv}"};

const int dca[2] = {3, 1};
const int nhitsfit[2] = {15, 25};
const int nhitsdedx[2] = {10, 15};
const double nsigmapi[3] = {-1.5, -2.0, -2.5};

// histogram //
TH1D *hCorr[2][nPtBins][2];// frame & pt & angle
TF1 * fLambda[2][nPtBins][2];// frame & pt & angle
TH1D *hIter[2][nPtBins][3];
TH1D *hIterChi[2][nPtBins];
TH1D *hIterCovar[2][nPtBins];

// parameters //
double parLambda[2][nPtBins][3];
double parLambdaErr[2][nPtBins][3];// theta phi inv
double parChi2[2][nPtBins];
int parNDF[2][nPtBins];
double parCovar[2][nPtBins];

// data need to be globals to be visible by fcn
TH1D *h1, *h2;
Int_t npfits;//used in myFcn
const string parName[4] = {"norm1","lambdaTheta","lambdaPhi","norm2"};
double iniPar[4] = { 1, 0.15, 0.3, 1 };
//double iniVlow[4] = {0, -1, -1, 0};
//double iniVhigh[4] = {0, 1, 1, 0};
double iniVlow[4] = {0, 0, 0, 0};
double iniVhigh[4] = {0, 0, 0, 0};
double stp[4] = {0.01, 0.01, 0.01, 0.01};

void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
    TAxis *xaxis1  = h1->GetXaxis();
    TAxis *xaxis2  = h2->GetXaxis();
    
    int nbinX1 = h1->GetNbinsX();
    int nbinX2 = h2->GetNbinsX();
    
    double chi2 = 0;
    double x;
    double tmp;
    npfits = 0;
    for (int ix = 1; ix <= nbinX1; ++ix) {
        x = xaxis1->GetBinCenter(ix);
        if ( h1->GetBinError(ix) > 0 ) {
            tmp = ( h1->GetBinContent(ix) - p[0]*(3*(1+p[1]*x*x))/(2*(p[1]+3)))/h1->GetBinError(ix);
            chi2 += tmp*tmp;
            npfits++;
        }
    }
    for (int ix = 1; ix <= nbinX2; ++ix) {
        x = xaxis2->GetBinCenter(ix);
        if ( h2->GetBinError(ix) > 0 ) {
            tmp = (h2->GetBinContent(ix) - p[3]*((1+2*p[2]*cos(2*x)/(3+p[1]))/(pi)))/h2->GetBinError(ix);
            chi2 += tmp*tmp;
            npfits++;
        }
    }
    fval = chi2;
}

