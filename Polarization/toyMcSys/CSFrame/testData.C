//--- including muon distribution histogram

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cmath>

#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFitResult.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "THnSparse.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "TTimer.h"
using namespace std;

const Double_t jpsiMass = 3.096;
const Double_t muMass = 0.105658367;
const Double_t PI = TMath::Pi();
const double ptLow = 0;
const double ptHigh = 10;
double eff1, eff2;
const int nTry = 500;
const int nPolBins = 150;

//-- histogram
TH2F  *mhMcThetaVsPt;
TH2F  *mhMcPhiVsPt;
TH2F  *mhRcThetaVsPt;
TH2F  *mhRcPhiVsPt;
TH2F  *mhMcThetaVsPtExpr[nTry];
TH2F  *mhMcPhiVsPtExpr[nTry];
TH2F  *mhRcThetaVsPtExpr[nTry];
TH2F  *mhRcPhiVsPtExpr[nTry];
THnSparse *mhnMcThetaVsPt;
THnSparse *mhnMcPhiVsPt;
THnSparse *mhnRcThetaVsPt;
THnSparse *mhnRcPhiVsPt;

//-- TRandom
TRandom3 myRandom(0);

//-- function
bool passJpsiPt(TF1 *f, double pt);
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass);
bool passPolar(TF2 *f2, double ptheta, double pphi, double pThetaPhi, double theta, double phi );
Double_t weightPola( double ptheta, double pphi, double pThetaPhi, double theta, double phi);
Double_t weightPola( double ptheta, double pphi, double pThetaPhi, double theta, double phi, TF2 *f2);

void calPolarizationHX(TLorentzVector iVec, TLorentzVector vec, double *ptheta, double *pphi, double pThetaPhi, TF2 *f2, bool doWeight, int type, int iter);
void calPolarizationCS(TLorentzVector iVec, TLorentzVector vec, double *ptheta, double *pphi, double pThetaPhi, TF2 *f2, bool doWeight, int type, int iter);

double calMuonEff(TH3D *hRc, double pt, double eta, double phi);
void bookHistograms(int iter, int loop);
void writeHistograms(const char* outFile, int n, int m);
bool passTrack(double pt, double eta, double phi);

bool debug = false;
const char *myDir = "/star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/toyMcSys/CSFrame";

int testData( const int nExpr=1e5, double LamThe= 0.5, double LamPhi= 0.5, double LamThePhi=0, int nLoop=1, int iter=0)
{
  myRandom.SetSeed(1000);//For efficiency set seed 1000; For generate pseudo data set seed 0
  TFile *fEff = new TFile(Form("%s/Rootfiles/muonEffNew.histo.root", myDir),"read");
  TH3D *mhMuonPtEtaPhi_MC_p = (TH3D*) fEff->Get("mhMuonPtEtaPhi_MC_p");
  TH3D *mhMuonPtEtaPhi_MtdTrig_p = (TH3D*)fEff->Get("mhMuonPtEtaPhi_MtdTrig_p");
  mhMuonPtEtaPhi_MtdTrig_p->Divide(mhMuonPtEtaPhi_MC_p);
  TH3D *mhMuonPtEtaPhi_MC_m = (TH3D*) fEff->Get("mhMuonPtEtaPhi_MC_m");
  TH3D *mhMuonPtEtaPhi_MtdTrig_m = (TH3D*)fEff->Get("mhMuonPtEtaPhi_MtdTrig_m");
  mhMuonPtEtaPhi_MtdTrig_m->Divide(mhMuonPtEtaPhi_MC_m);

  // CS frame //
  const int nthetaphi = 4;
  double lamTheArray[nthetaphi] = {-0.26, 0.48, 0.63, -0.047};
  double lamPhiArray[nthetaphi] = {-0.034, 0.008, -0.124, -0.3};

  int polIndex = -1;
  for(int it=0; it<nthetaphi; it++)
    {
      if( (lamTheArray[it]==LamThe) && (lamPhiArray[it]==LamPhi)) 
	{ polIndex = it; break; }
    }
  printf("[i] Input theta = %4.2f, phi = %4.2f, index = %d\n",LamThe,LamPhi,polIndex);

  double inLamThe[4][nTry];
  double inLamPhi[4][nTry];
  if(iter<=0 || iter==10)
    {
      for(int p=0; p<4; p++)
	{
	  for(int ibin=0; ibin<nTry; ibin++)
	    {
	      inLamThe[p][ibin] = LamThe;
	      inLamPhi[p][ibin] = LamPhi;
	    }
	}
    }
  else if(iter==1)
    {
      for(int p=0; p<4; p++)
	{
	  for(int ibin=0; ibin<nTry; ibin++)
	    {
	      inLamThe[p][ibin] = 0.0;
	      inLamPhi[p][ibin] = 0.0;
	    }
	}		  
    }
  else if(iter>=2 && iter<10)
    {
      if(polIndex ==-1)
	{
	  printf("[e] Wrong input values (%1.1f, %1.1f)\n",LamThe,LamPhi);
	  return -1;
	}
      TFile *fPar = new TFile(Form("%s/Resultfiles/Run15_pp200.JpsiPolPar.root", myDir),"read");
      for(int p=0; p<4; p++)
	{
	  TH1F *hPar = (TH1F*)fPar->Get(Form("Iter%d_hFitJpsiPol_Theta_Pt%d_%3.3f_%3.3f",iter-1,p,lamTheArray[polIndex],lamPhiArray[polIndex]));
	  for(int ibin=0; ibin<nTry; ibin++)
	    {
	      inLamThe[p][ibin] = hPar->GetBinContent(ibin+1);
	      if(inLamThe[p][ibin]>1) inLamThe[p][ibin] = 1;
	      if(inLamThe[p][ibin]<-1) inLamThe[p][ibin] = -1;
	    }
	  hPar = (TH1F*)fPar->Get(Form("Iter%d_hFitJpsiPol_Phi_Pt%d_%3.3f_%3.3f",iter-1,p,lamTheArray[polIndex],lamPhiArray[polIndex]));
	  for(int ibin=0; ibin<nTry; ibin++)
	    {
	      inLamPhi[p][ibin] = hPar->GetBinContent(ibin+1);
	      if(inLamPhi[p][ibin]>1) inLamPhi[p][ibin] = 1;
	      if(inLamPhi[p][ibin]<-1) inLamPhi[p][ibin] = -1;
	    }
	}
    }

  TF1* fPt = new TF1("f","([0]*TMath::Power(TMath::Exp([1]*x)+x/[2],[3]))*4*TMath::Pi()*x",0,20);//published
  fPt->SetParameters(7.287, -0.235, 3.412, -10.533);
  TF2 *fPolar = new TF2("fPolar", "1 + [0]*TMath::Cos(x)*TMath::Cos(x) + [1]*(1-TMath::Cos(x)*TMath::Cos(x))*TMath::Cos(2*y) + [2]*TMath::Sin(2.*x)*TMath::Cos(y)",-PI,PI,-PI,PI);

  gStyle->SetOptStat("ne");
  bookHistograms(iter, nLoop);
  for(int i=0; i<nExpr; i++)
    {
      if(i%(nExpr/10)==0) cout<<((double)i/nExpr*100)<<"% has compelete"<<endl;
      if(debug) cout<<i<<endl;

      double mc_pt  = fPt->GetRandom(ptLow, ptHigh);
      double mc_phi = myRandom.Uniform(-1*PI, PI);
      double mc_y   = myRandom.Uniform(-0.5, 0.5);
      double mc_px = mc_pt * TMath::Cos(mc_phi);
      double mc_py = mc_pt * TMath::Sin(mc_phi);
      double mc_pz = sqrt(mc_pt*mc_pt+jpsiMass*jpsiMass) * TMath::SinH(mc_y);
      TLorentzVector parent;
      parent.SetXYZM(mc_px,mc_py,mc_pz,jpsiMass);
      if(debug) cout<<" jpsi info; pt:"<<mc_pt<<"; phi:"<<mc_phi<<endl;

      //-- daughter
      TLorentzVector daughter1 = twoBodyDecay(parent,muMass);
      TLorentzVector daughter2 = parent - daughter1;
      int pt_index = -1;
      if(mc_pt < 1)      pt_index = 0;
      else if(mc_pt < 2) pt_index = 1;
      else if(mc_pt < 4) pt_index = 2;
      else               pt_index = 3;

      calPolarizationCS(daughter1, parent, inLamThe[pt_index], inLamPhi[pt_index], LamThePhi, fPolar, true, 0, iter);

      if(!passTrack(daughter1.Pt(),daughter1.Eta(),daughter1.Phi())) continue;
      if(!passTrack(daughter2.Pt(),daughter2.Eta(),daughter2.Phi())) continue;

      eff1 = calMuonEff(mhMuonPtEtaPhi_MtdTrig_p,daughter1.Pt(),daughter1.Eta(),daughter1.Phi());
      eff2 = calMuonEff(mhMuonPtEtaPhi_MtdTrig_m,daughter2.Pt(),daughter2.Eta(),daughter2.Phi());

      calPolarizationCS(daughter1, parent, inLamThe[pt_index], inLamPhi[pt_index], LamThePhi, fPolar, true, 1, iter);
    }
  char *outName;
  if(iter<=0) outName = Form("Data_%3.3f_%3.3f_%3.1f",LamThe,LamPhi,LamThePhi);
  else        outName = Form("Eff_%3.3f_%3.3f_%3.1f",LamThe,LamPhi,LamThePhi);
  writeHistograms(outName,nLoop,iter);
  fEff->Close();

  return 0;
}
//-------------------------------------------------------
bool passJpsiPt(TF1 *f, double pt)
{
  double fmax = f->GetMaximum();
  double eff = f->Eval(pt);
  double filter = myRandom.Uniform(0,fmax);
  if(filter<=eff) return true;
  else return false;
}

//-------------------------------------------------------
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter)
{
  float betax = parent.Px()/parent.E();
  float betay = parent.Py()/parent.E();
  float betaz = parent.Pz()/parent.E();
  daughter.Boost(betax,betay,betaz);
  return daughter;
}
//-------------------------------------------------------
TLorentzVector twoBodyDecay(TLorentzVector parent, double dmass)
{
  double e = parent.M()/2.;
  double p = sqrt(e*e-dmass*dmass);
  double costheta = myRandom.Uniform(-1.0,1.0);
  double phi = myRandom.Uniform(-TMath::Pi(),TMath::Pi());
  double pz = p*costheta;
  double px = p*sqrt(1.-costheta*costheta)*cos(phi);
  double py = p*sqrt(1.-costheta*costheta)*sin(phi);
  TLorentzVector daughter(px,py,pz,e);
  return myBoost(parent,daughter);
}

//-------------------------------------------------------
bool passPolar(TF2 *f2, double ptheta, double pphi, double pThetaPhi, double theta, double phi )
{
  //x_theta & y_phi
  f2->SetParameters(ptheta, pphi, pThetaPhi);
  double fmax = f2->GetMaximum();
  double eff = f2->Eval(theta, phi);
  double filter = myRandom.Uniform(0,fmax);
  if( filter <= eff ) return true;
  else return false;
}

//-------------------------------------------------------
Double_t weightPola( double ptheta, double pphi, double pThetaPhi, double theta, double phi )
{
  double costheta = TMath::Cos(theta);
  double w = 1 + ptheta*costheta*costheta + pphi*(1-costheta*costheta)*TMath::Cos(2*phi) + pThetaPhi*TMath::Sin(2.*theta)*TMath::Cos(phi);
  return w;
}

//-------------------------------------------------------
Double_t weightPola( double ptheta, double pphi, double pThetaPhi, double theta, double phi, TF2 *f2 )
{
  f2->SetParameters(ptheta,pphi,pThetaPhi);
  double w = f2->Eval(theta, phi);
  return w;
}

//-------------------------------------------------------
void calPolarizationHX(TLorentzVector iVec, TLorentzVector vec, double *ptheta, double *pphi, double pThetaPhi, TF2 *f2, bool doWeight, int type, int iter)
{
  // type - 0: Mc; type - 1: Rc;

  //iVec:positron TLorentzVector  vec:JPSI LorentzVector
  TLorentzVector mu(iVec);//positron
  TLorentzVector Proton1(0.,0.,100.,100.),Proton2(0.,0.,-100.,100.);//pp200GeV

  TVector3 XXHX,YYHX,ZZHX;//hilicity frame XYZaxis
  ZZHX = vec.Vect();
  YYHX = vec.Vect().Cross(Proton1.Vect());
  XXHX = YYHX.Cross(ZZHX);

  TVector3 jpsi = vec.BoostVector();

  mu.Boost(-1*jpsi);

  Float_t theta = mu.Angle(ZZHX);
  Float_t phi = TMath::ATan2((mu.Vect().Dot(YYHX.Unit())),(mu.Vect().Dot(XXHX.Unit())));
  Float_t absPhi = fabs(phi);
  Float_t cosTheta = TMath::Cos(theta);
  Float_t pt = vec.Pt();

  int pt_index = -1;
  if(pt < 1)      pt_index = 0;
  else if(pt < 2) pt_index = 1;
  else if(pt < 4) pt_index = 2;
  else            pt_index = 3;
    
  if(doWeight)
    {
      if(iter<=1 || iter==10)
	{
	  double wCS = weightPola(ptheta[0], pphi[0], pThetaPhi, theta, phi, f2);
	  if(type==0) 
	    {
	      mhMcThetaVsPt->Fill(pt, cosTheta, wCS);
	      mhMcPhiVsPt->Fill(pt, absPhi, wCS);
	    }
	  else if(type==1)
	    {
	      mhRcThetaVsPt->Fill(pt, cosTheta, wCS*eff1*eff2);
	      mhRcPhiVsPt->Fill(pt, absPhi, wCS*eff1*eff2); 
	    }
	}
      else if(iter==20)
	{
	  double step = 3./nPolBins;
	  for(int it=0; it<nPolBins; it++)
	    {
	      for(int ip=0; ip<nPolBins; ip++)
		{
		  double ltheta = -1.5 + it * step;
		  double lphi   = -1.5 + ip * step;
		  double wCS = weightPola(ltheta, lphi, pThetaPhi, theta, phi, f2);
		  double fill_theta[4] = {ltheta+step/2, lphi+step/2, cosTheta, pt_index*1.};
		  if(type==0) mhnMcThetaVsPt->Fill(fill_theta, wCS);
		  else        mhnRcThetaVsPt->Fill(fill_theta, wCS*eff1*eff2);
		      
		  double fill_phi[4] = {ltheta+step/2, lphi+step/2, absPhi, pt_index*1.};
		  if(type==0) mhnMcPhiVsPt->Fill(fill_phi, wCS);
		  else        mhnRcPhiVsPt->Fill(fill_phi, wCS*eff1*eff2);
		}
	    }
	}
      else
	{
	  for(int j=0; j<nTry; j++)
	    {
	      double wCS = weightPola(ptheta[j], pphi[j], pThetaPhi, theta, phi, f2);
	      if(type==0) 
		{
		  mhMcThetaVsPtExpr[j]->Fill(pt, cosTheta, wCS);
		  mhMcPhiVsPtExpr[j]->Fill(pt, absPhi, wCS);
		}
	      else if(type==1)
		{
		  mhRcThetaVsPtExpr[j]->Fill(pt, cosTheta, wCS*eff1*eff2);
		  mhRcPhiVsPtExpr[j]->Fill(pt, absPhi, wCS*eff1*eff2);
		}
	    }
	}
    }

}

//-------------------------------------------------------
void calPolarizationCS(TLorentzVector iVec, TLorentzVector vec, double *ptheta, double *pphi, double pThetaPhi, TF2 *f2, bool doWeight, int type, int iter)
{
    // type - 0: Mc; type - 1: Rc;
    
    //iVec:positron TLorentzVector  vec:JPSI LorentzVector
    TLorentzVector mu(iVec);//positron
    TLorentzVector Proton1(0.,0.,100.,100.),Proton2(0.,0.,-100.,100.);//pp200GeV
    
    TVector3 jpsi = vec.BoostVector();
    mu.Boost(-1*jpsi);
    
    Proton1.Boost(-1*jpsi);
    Proton2.Boost(-1*jpsi);
    TVector3 XX,YY,ZZ;//Collins-Soper frame
    ZZ = Proton1.Vect()*(1/Proton1.Vect().Mag())-Proton2.Vect()*(1/Proton2.Vect().Mag());
    YY = Proton1.Vect().Cross(Proton2.Vect());
    XX = YY.Cross(ZZ);
    
    Float_t theta = mu.Angle(ZZ);
    Float_t phi = TMath::ATan2((mu.Vect().Dot(YY.Unit())),(mu.Vect().Dot(XX.Unit())));
    Float_t absPhi = fabs(phi);
    Float_t cosTheta = TMath::Cos(theta);
    Float_t pt = vec.Pt();
    
    int pt_index = -1;
    if(pt < 1)      pt_index = 0;
    else if(pt < 2) pt_index = 1;
    else if(pt < 4) pt_index = 2;
    else            pt_index = 3;
    
    if(doWeight)
    {
        if(iter<=1 || iter==10)
        {
            double wCS = weightPola(ptheta[0], pphi[0], pThetaPhi, theta, phi, f2);
            if(type==0)
            {
                mhMcThetaVsPt->Fill(pt, cosTheta, wCS);
                mhMcPhiVsPt->Fill(pt, absPhi, wCS);
            }
            else if(type==1)
            {
                mhRcThetaVsPt->Fill(pt, cosTheta, wCS*eff1*eff2);
                mhRcPhiVsPt->Fill(pt, absPhi, wCS*eff1*eff2);
            }
        }
        else if(iter==20)
        {
            double step = 3./nPolBins;
            for(int it=0; it<nPolBins; it++)
            {
                for(int ip=0; ip<nPolBins; ip++)
                {
                    double ltheta = -1.5 + it * step;
                    double lphi   = -1.5 + ip * step;
                    double wCS = weightPola(ltheta, lphi, pThetaPhi, theta, phi, f2);
                    double fill_theta[4] = {ltheta+step/2, lphi+step/2, cosTheta, pt_index*1.};
                    if(type==0) mhnMcThetaVsPt->Fill(fill_theta, wCS);
                    else        mhnRcThetaVsPt->Fill(fill_theta, wCS*eff1*eff2);
                    
                    double fill_phi[4] = {ltheta+step/2, lphi+step/2, absPhi, pt_index*1.};
                    if(type==0) mhnMcPhiVsPt->Fill(fill_phi, wCS);
                    else        mhnRcPhiVsPt->Fill(fill_phi, wCS*eff1*eff2);
                }
            }
        }
        else
        {
            for(int j=0; j<nTry; j++)
            {
                double wCS = weightPola(ptheta[j], pphi[j], pThetaPhi, theta, phi, f2);
                if(type==0)
                {
                    mhMcThetaVsPtExpr[j]->Fill(pt, cosTheta, wCS);
                    mhMcPhiVsPtExpr[j]->Fill(pt, absPhi, wCS);
                }
                else if(type==1)
                {
                    mhRcThetaVsPtExpr[j]->Fill(pt, cosTheta, wCS*eff1*eff2);
                    mhRcPhiVsPtExpr[j]->Fill(pt, absPhi, wCS*eff1*eff2);
                }
            }
        }
    }
    
}

//-------------------------------------------------------
bool passTrack(double pt, double eta, double phi)
{
  if(pt>1.3 && eta<0.5 && eta>-0.5) return true;
  else return false;
}

//-------------------------------------------------------
double calMuonEff(TH3D *hRc, double pt, double eta, double phi)
{
  double ptnew = pt;
  if(pt<1.5) ptnew = 1.6;
  if(pt>10) ptnew = 9.5;
  int bin = hRc->FindBin(ptnew, eta, phi);
  double eff = hRc->GetBinContent(bin);
  return eff;
}

//-------------------------------------------------------
void bookHistograms(int iter, int loop)
{
  const int nPtBins = 4;
  const double xPtBins[nPtBins+1] = {0, 1, 2, 4, 10};
  const int nThetaBins = 10;
  const double minTheta = -1.0, maxTheta = 1.0;
  const int nPhiBins = 15;
  const double minPhi = 0, maxPhi = TMath::Pi();

  if(iter<=1 || iter==10)
    {
      char *suffix = Form("_Loop%d",loop);
      if(iter<0) suffix = Form("_Loop%d",0);
      else if(iter==1 || iter==10) suffix = Form("");
      mhMcThetaVsPt = new TH2F(Form("mhMcThetaVsPt%s",suffix), ";p_{T} (GeV/c);cos(#theta)", nPtBins, xPtBins, nThetaBins, minTheta, maxTheta);
      mhMcPhiVsPt   = new TH2F(Form("mhMcPhiVsPt%s",suffix), ";p_{T} (GeV/c);#varphi", nPtBins, xPtBins, nPhiBins, minPhi, maxPhi);
      mhRcThetaVsPt = new TH2F(Form("mhRcThetaVsPt%s",suffix), ";p_{T} (GeV/c);cos(#theta)", nPtBins, xPtBins, nThetaBins, minTheta, maxTheta);
      mhRcPhiVsPt   = new TH2F(Form("mhRcPhiVsPt%s",suffix), ";p_{T} (GeV/c);#varphi", nPtBins, xPtBins, nPhiBins, minPhi, maxPhi);
      mhMcThetaVsPt->Sumw2();
      mhMcPhiVsPt->Sumw2();
      mhRcThetaVsPt->Sumw2();
      mhRcPhiVsPt->Sumw2();
    }
  else if(iter==20)
    {
      const int dim = 4;
      const int nBinsTheta[dim] = {nPolBins, nPolBins, 10, 4};
      const double minBinTheta[dim] = {-1.5, -1.5, -1, 0};
      const double maxBinTheta[dim] = {1.5, 1.5, 1, 4};
      const int nBinsPhi[dim] = {nPolBins, nPolBins, 15, 4};
      const double minBinPhi[dim] = {-1.5, -1.5, 0, 0};
      const double maxBinPhi[dim] = {1.5, 1.5, TMath::Pi(), 4};
      char *suffix = Form("");

      mhnMcThetaVsPt = new THnSparseF(Form("mhnMcThetaVsPt%s",suffix), ";#lambda_{#theta};#lambda_{#varphi};cos(#theta);p_{T} (GeV/c)", dim, nBinsTheta, minBinTheta, maxBinTheta);
      mhnMcThetaVsPt->Sumw2();

      mhnMcPhiVsPt = new THnSparseF(Form("mhnMcPhiVsPt%s",suffix), ";#lambda_{#theta};#lambda_{#varphi};#varphi;p_{T} (GeV/c)", dim, nBinsPhi, minBinPhi, maxBinPhi);
      mhnMcPhiVsPt->Sumw2();

      mhnRcThetaVsPt = new THnSparseF(Form("mhnRcThetaVsPt%s",suffix), ";#lambda_{#theta};#lambda_{#varphi};cos(#theta);p_{T} (GeV/c)", dim, nBinsTheta, minBinTheta, maxBinTheta);
      mhnRcThetaVsPt->Sumw2();

      mhnRcPhiVsPt = new THnSparseF(Form("mhnRcPhiVsPt%s",suffix), ";#lambda_{#theta};#lambda_{#varphi};#varphi;p_{T} (GeV/c)", dim, nBinsPhi, minBinPhi, maxBinPhi);
      mhnRcPhiVsPt->Sumw2();
    }
  else
    {
     for(int i=0; i<nTry; i++)
	{
	  mhMcThetaVsPtExpr[i] = new TH2F(Form("mhMcThetaVsPt_Expr%d",i), ";p_{T} (GeV/c);cos(#theta)", nPtBins, xPtBins, nThetaBins, minTheta, maxTheta);
	  mhMcPhiVsPtExpr[i]   = new TH2F(Form("mhMcPhiVsPt_Expr%d",i), ";p_{T} (GeV/c);#varphi", nPtBins, xPtBins, nPhiBins, minPhi, maxPhi);
	  mhRcThetaVsPtExpr[i] = new TH2F(Form("mhRcThetaVsPt_Expr%d",i), ";p_{T} (GeV/c);cos(#theta)", nPtBins, xPtBins, nThetaBins, minTheta, maxTheta);
	  mhRcPhiVsPtExpr[i]   = new TH2F(Form("mhRcPhiVsPt_Expr%d",i), ";p_{T} (GeV/c);#varphi", nPtBins, xPtBins, nPhiBins, minPhi, maxPhi);
	  mhMcThetaVsPtExpr[i]->Sumw2();
	  mhMcPhiVsPtExpr[i]->Sumw2();
	  mhRcThetaVsPtExpr[i]->Sumw2();
	  mhRcPhiVsPtExpr[i]->Sumw2();
	}
    }

}
//-------------------------------------------------------
void writeHistograms(const char* outFile, int n, int m)
{
  char buf[1024];
  TFile *mFile = new TFile(Form("%s/output/Iter%d.Loop%d.%s.root", myDir, m,n,outFile),"recreate");
  mFile->cd();

  if(m<=1 || m==10)
    {
      mhMcThetaVsPt->Write();
      mhMcPhiVsPt->Write();
      mhRcThetaVsPt->Write();
      mhRcPhiVsPt->Write();
    }
  else if(m==20)
    {
      mhnMcThetaVsPt->Write();
      mhnMcPhiVsPt->Write();
      mhnRcThetaVsPt->Write();
      mhnRcPhiVsPt->Write();
    }
  else
    {
      for(int i=0; i<nTry; i++)
	{
	  mhMcThetaVsPtExpr[i]->Write();
	  mhMcPhiVsPtExpr[i]->Write();
	  mhRcThetaVsPtExpr[i]->Write();
	  mhRcPhiVsPtExpr[i]->Write();
	}
    }

  mFile->Close();
}
