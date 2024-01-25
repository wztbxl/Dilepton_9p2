const double varDca[2] = {3, 1};
const double varNhitsFit[2] = {15, 25};
const double varNhitsDedx[2] = {10, 15};
const double varNsigmaPi[3] = {-1.5, -2, -2.5};

const char *signName[2] = {"UL", "LS"};
const char *frameName[2] = {"HX", "CS"};
const char *angleName[2] = {"Cos", "Phi"};

const int nPtBins = 4;
double xPtBins[nPtBins+1] ={0, 1, 2, 4, 10};

const int nCosBinsPt = 10;
double xCosBinsPt[nCosBinsPt+1] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0 };

const double PI = TMath::Pi();
const int nPhiBinsPt = 15;
double xPhiBinsPt[nPhiBinsPt+1] = {
    0, 1.*PI/15., 2.*PI/15., 3.*PI/15., 4.*PI/15, PI/3., PI*6./15., PI*7./15., PI*8./15., PI*9./15., PI*2./3., PI*11./15., PI*12./15., 13.*PI/15, 14.*PI/15, PI };

int cutBin[2][nPtBins][2] = {
    {{0, 0}, {2, 5}, {2, 4}, {1, 3}},//HX
    {{2, 0}, {2, 4}, {2, 4}, {1, 4}}//CS
};// frame & pt & angle

double sigmaPt[nPtBins] = {0.04, 0.04, 0.05, 0.06};// frame & pt

double fitCosRangeLow[2][4][nCosBinsPt] = {
    {{2.8, 2.8, 2.8, 2.8, 2.8,  2.8, 2.8, 2.8, 2.8, 2.8},
    {3.2, 3.2, 2.8, 2.6, 2.6,  2.6, 2.6, 2.8, 3.2, 3.2},
    {3.2, 3.2, 2.6, 2.6, 2.6,  2.6, 2.4, 2.6, 3.2, 3.2},
        {3.2, 2.6, 2.6, 2.6, 2.4,  2.4, 2.4, 2.6, 2.6, 3.2}},//HX & pt
    {{3.2, 3.2, 2.8, 2.7, 2.6,  2.6, 2.7, 2.8, 3.2, 3.2},
        {3.2, 3.2, 2.4, 2.4, 2.4,  2.4, 2.4, 2.4, 3.2, 3.2},
        {3.2, 2.2, 2.4, 2.4, 2.4,  2.4, 2.4, 2.4, 2.2, 3.2},
        {3.2, 2.2, 2.4, 2.4, 2.4,  2.4, 2.4, 2.4, 2.2, 3.2} } //CS & pt
};
double fitCosRangeUp[2][4][nCosBinsPt] = {
    {{4, 4, 4, 4, 4,  4, 4, 4, 4, 4},
    {4, 4, 4, 4, 4,  4, 4, 4, 4, 4},
    {4, 4, 4, 4, 4,  4, 4, 4, 4, 4},
        {4, 4, 4, 4, 4,  4, 4, 4, 4, 4}},//HX & pt
    {{4, 4, 4, 4, 4,  4, 4, 4, 4, 4},
        {4, 4, 4, 4, 4,  4, 4, 4, 4, 4},
        {4, 4, 4, 4, 4,  4, 4, 4, 4, 4},
        {4, 4, 4, 4, 4,  4, 4, 4, 4, 4}}// HX & pt
};

double fitPhiRangeLow[2][4][nPhiBinsPt] = {
    {{2.8, 2.8, 2.8, 2.8,  2.8, 2.8, 2.6, 2.6,  2.6, 2.8, 2.8, 2.8,  2.8, 2.8, 2.8},
    {3.4, 3.4, 3.4, 3.4,  3.2, 2.6, 2.6, 2.6,  2.6, 2.6, 3.2, 3.2,  3.2, 3.2, 3.2},
    {3.2, 3.2, 3.2, 3.2,  2.6, 2.6, 2.6, 2.6,  2.6, 2.6, 2.7, 3.2,  3.2, 3.2, 3.2},
        {3.2, 3.2, 3.2, 2.6,  2.6, 2.6, 2.6, 2.6,  2.6, 2.6, 2.8, 2.6,  3.2, 3.2, 3.2}},
    {{2.8, 2.8, 2.8, 2.8,  2.8, 2.8, 2.6, 2.5,  2.6, 2.7, 2.8, 2.8,  2.8, 2.8, 2.8},
        {3.4, 3.4, 3.4, 3.4,  3., 2.7, 2.4, 2.6,  2.6, 2.7, 3., 3.2,  3.2, 3.2, 3.2 },
        {3.2, 3.2, 3.2, 3.2,  2.8, 2.6, 2.6, 2.6,  2.6, 2.4, 2.8, 3.2,  3.2, 3.2, 3.2},
        {3.2, 3.2, 3.2, 3.2,  2.6, 2.6, 2.6, 2.6,  2.6, 2.6, 2.6, 3.2,  3.2, 3.2, 3.2}}
};
double fitPhiRangeUp[2][4][nPhiBinsPt] ={
    {{4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4},
    {4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4},
    {4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4},
        {4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4}},
    {{4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4},
        {4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4},
        {4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4},
        {4, 4, 4, 3.4,  4, 4, 4, 4,  4, 4, 4, 3.4,  4, 4, 4}}
};

double massLowPt[nPtBins] = {2.8, 2.8, 2.8, 2.8 };
double massHighPt[nPtBins] = {3.4, 3.4, 3.4, 3.4};

double shift[5] = {0, 0.1, -0.1, 0.05, -0.05};

// File //
TFile *fin;
TFile *fout;

// histograme //
THnSparseF *hnMPtCosPhi[2][2];//frame & sign
TH1D **hInvMass[2][2][nPtBins][2];//frame & sign & pt & angle

TH1D *hCosFit[2][nPtBins];
TH1D *hCosCount[2][nPtBins];
TH1D *hPhiFit[2][nPtBins];
TH1D *hPhiCount[2][nPtBins];

double *matrixCosPram[2][nPtBins][nCosBinsPt];// frame & pt & angle
double *matrixPhiPram[2][nPtBins][nPhiBinsPt];

// functions //
bool reject=true;
double flinePol2( double *x, double *par);
double flinePol3( double *x, double *par);
double flinePol4( double *x, double *par);
double flinePol5( double *x, double *par);



