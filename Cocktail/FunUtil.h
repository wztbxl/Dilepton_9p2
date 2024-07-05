#ifndef FUNUTIL_H 
#define FUNUTIL_H 

#include "MesonConstant.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"

//mass function for rho
//Double_t 	Gamma(Double_t mass); // p-wave Gamma(m) for rho->pipi
//Double_t 	GammaE(Double_t mass);// Gamma(m) for rho->ee
//Double_t 	PS(Double_t mass, Double_t pT, Double_t T);//phase-space factor
//Double_t 	massDistF0(Double_t *x, Double_t *par); // rho mass function without phase-space
//Double_t 	massDistF(Double_t *x, Double_t *par);  // rho mass function with phase-space
//Double_t tsallis(double* x, double* par)
//{
//	int   ptbin  = (int)(x[0]*100+0.5);
//	float ptlw   = histMeson->GetBinCenter(ptbin);
//	float ptup   = histMeson->GetBinCenter(ptbin+1);
//
//	float a1 = histMeson->GetBinContent(ptbin);
//	float a2 = histMeson->GetBinContent(ptbin+1);
//	float p0 = (a1-a2)/(ptlw-ptup);
//	float p1 = (a2*ptlw-a1*ptup)/(ptlw-ptup);
//
//	return par[0]*(p0*x[0]+p1)*x[0];
//}

const Double_t alpha = 1./137.036;
const Double_t Pi	 = 3.14159265;

/*
//------------------------------------------------------ //for new pT smear method,do not use crystallball function, assume pT resolution is decided by detector
TH1D *pTRes1D[19];
TFile *pTResInput = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/pTresFromEmbedding.root");
TH2D *PtRes2D = (TH2D*)pTResInput->Get("hPtResvsPtQ_zy");

void InitpTRes()
{
	for(int i = 0;i<19;i++)
	{
		int PtBinLow = PtRes2D->FindBin((i+1)*0.2+1.e-4);
		int PtBinHigh = PtRes2D->FindBin((i+2)*0.2-1.e-4); 
		pTRes1D[i] = pTRes2D->ProjectionY(Form("pTRes1D_%d",i));
	}//0.2 step is beacuse the embedding data is a small sample;
}
Double_t GetSmear(Double_t pT)
{
	double smear;
	int pTBin = int(pT/0.2);
	cout <<"pTBin = "<<pTBin<<endl;
	pTBin = pTBin-1;
	if (pTBin>=0 && pTBin <19)
		smear = pTRes1D[pTBin]->GetRandom();
	else if (pTBin < 0) smear = pTRes1D[0]->GetRandom();
	else if (pTBin >=19) smear = pTRes1D[18]->GetRandom();

	return smear;
}*/

//------------------------------------------------------ //bremsstrahlung tail
Double_t CrystalBall(Double_t *x, Double_t *par)
{
	Double_t N = par[0];
	Double_t mu = par[1];
	Double_t s = par[2];
	Double_t n = par[3];
	Double_t alpha1 = par[4];

	Double_t A = TMath::Power(n/fabs(alpha1), n) * TMath::Exp(-alpha1*alpha1/2.);
	Double_t B = n/fabs(alpha1) - fabs(alpha1);
	Double_t norm = (x[0]-mu)/s;

	if(norm > -alpha1) {
		return N * TMath::Exp(-0.5*norm*norm);
	} else {
		return N * A * TMath::Power(B-norm, -n);
	}
}

//------------------------------------------------------
Double_t CrystalBall2(Double_t *x, Double_t *par)
{
	Double_t N = par[0];
	Double_t mu = par[1];
	Double_t s = par[2];
	Double_t n = par[3];
	Double_t alpha2 = par[4];
	Double_t m = par[5];
	Double_t beta = par[6];

	Double_t A = TMath::Power(n/fabs(alpha2), n) * TMath::Exp(-alpha2*alpha2/2.);
	Double_t B = n/fabs(alpha2) - fabs(alpha2);

	Double_t C = TMath::Power(m/fabs(beta), m) * TMath::Exp(-beta*beta/2.);
	Double_t D = m/fabs(beta) - fabs(beta);

	Double_t norm = (x[0]-mu)/s;

	if(norm < -alpha2) {
		return N * A * TMath::Power(B-norm, -n);
	} else if(norm < beta) {
		return N * TMath::Exp(-0.5*norm*norm);
	} else {
		return N * C * TMath::Power(D+norm, -m);
	}
}
Double_t DoubleCrystalBall(double *x, double *p){//you can see twice as far into the future

  Double_t xx = x[0];
  Double_t alpha1 = p[0];
  Double_t alpha2 = p[1];
  Double_t n1 = p[2];
  Double_t n2 = p[3];
  Double_t xxmean = p[4];
  Double_t sigma = p[5];
  Double_t N = p[6];

  Double_t jmb = (xx - xxmean)/sigma;

  Double_t A = TMath::Power(n1/TMath::Abs(alpha1),n1)*TMath::Exp(-alpha1*alpha1/2.);
  Double_t B = n1/TMath::Abs(alpha1)-TMath::Abs(alpha1);
  Double_t C = TMath::Power(n2/TMath::Abs(alpha2),n2)*TMath::Exp(-alpha2*alpha2/2.);
  Double_t D = n2/TMath::Abs(alpha2)-TMath::Abs(alpha2);

  if(jmb < -alpha1){
    return N*A*TMath::Power(B-jmb,-n1);
  }//alpha1
  else if(jmb <= alpha2){//to match huck
    return N*TMath::Exp(-jmb*jmb/2.);
  }//alpha2
  else{
    return N*C*TMath::Power((D+jmb),-n2);
  }//else

}//doublecrystalball

//------------------------------------------------------
Double_t momRes(Double_t *x, Double_t *par)
{
	double pt = x[0];
	//double m = par[0];
	/*
	   const Double_t a = 8.e-3; //3.2e-3; // AuAu - see from embedding
	   const Double_t b = 7.6e-3;
	//  const Double_t a = 8.e-3; // pp - see from embedding
	//  const Double_t b = 8.e-3;
	return TMath::Sqrt(a*a*pt*pt+b*b*m*m/pt/pt+b*b); //pi mom-res
	*/
	const Double_t a = 6.0e-3;
	const Double_t b = 8.3e-3;
	return TMath::Sqrt(a*a*pt*pt+b*b);
}

//-------------------------------------------------------
Double_t Gamma(Double_t mass)
{
	const Double_t Gamma0 = 0.1491;
	if(mass>2.*Masspi0) 
		return Gamma0*Massrho/mass*TMath::Power((mass*mass-4.*Masspi0*Masspi0)/(Massrho*Massrho-4.*Masspi0*Masspi0),1.5);
	else
		return 0;
}

//-------------------------------------------------------
Double_t GammaE(Double_t mass)
{
	const Double_t Gamma0 = 0.1491;
	if(mass>2.*Masselectron) 
		return Gamma0*Massrho/mass*TMath::Power((mass*mass-4.*Masselectron*Masselectron)/(Massrho*Massrho-4.*Masselectron*Masselectron),0.5);
	else
		return 0;
}
//-------------------------------------------------------
Double_t PS(Double_t mass, Double_t pT, Double_t T)
{
	if(mass>0&&pT>0){
		return mass/TMath::Sqrt(mass*mass+pT*pT)*TMath::Exp(-TMath::Sqrt(mass*mass+pT*pT)/T);
	}else 
		return 0;
}
//-------------------------------------------------------
Double_t massDistF0(Double_t *x, Double_t *par)
{
	const Double_t Gamma2 = 4.72e-5;
	Double_t m = x[0];
	return m*Massrho*GammaE(m)/(TMath::Power((Massrho*Massrho-m*m),2.)+Massrho*Massrho*(Gamma(m)+GammaE(m)*Gamma2)*(Gamma(m)+GammaE(m)*Gamma2));
}

//-------------------------------------------------------
Double_t massDistRho(Double_t *x, Double_t *par)
{
	const Double_t Gamma2 = 4.72e-5;
	double pT = par[0];
	double T = par[1];
	double m = x[0];
	return m*Massrho*GammaE(m)/(TMath::Power((Massrho*Massrho-m*m),2.)+Massrho*Massrho*(Gamma(m)+GammaE(m)*Gamma2)*(Gamma(m)+GammaE(m)*Gamma2))*PS(m,pT,T)*1e2;
}
// end of rho->ee functions old formula from Zhangbu.


//-------------------------------------------------------
Double_t massDistRhoNA60(Double_t *x, Double_t *par)
{

	const double mRho = 0.770;//GeV
	const double GammaRhotot = 0.150;//GeV
	const double Masspi = 0.13957;
	double m  = x[0];
	double pT = par[0];
	double T  = par[1];

	double a0 = alpha*alpha*TMath::Power(mRho,4.)/3./TMath::Power(2.*Pi,4.);
	if(m<2.*Masspi) return 0.;
	double QEDterm = TMath::Power(1.-4.*Masspi*Masspi/m/m,1.5)*TMath::Power(1.-4.*Masselectron*Masselectron/m/m,0.5)*(1.+2.*Masselectron*Masselectron/m/m);

	double a1 = TMath::Power(m*m-mRho*mRho,2.)+m*m*GammaRhotot*GammaRhotot;

	//double psfac_na60 = TMath::Power(2.*Pi*m*T,1.5)*exp(-m/T);

	return a0*QEDterm/a1*PS(m,pT,T);
	//return a0*QEDterm/a1*psfac_na60; // NA60 formula
}


//-------------------------------------------------------
Double_t QEDee(Double_t mass){

	if(mass>0){
		Double_t a1 = (1.+2.*Masselectron*Masselectron/mass/mass);
		Double_t a2 = TMath::Power(1.-4.*Masselectron*Masselectron/mass/mass,0.5);
		return alpha*2./3./Pi/mass*a1*a2;
	}else return 0;
}

//-------------------------------------------------------
Double_t QEDeeomega(Double_t mass){
	const Double_t GammaToPi0gamma = 1.;
	//const Double_t GammaToPi0gamma = 0.0828*8.49e-3;

	if(mass>0&&mass<Massomega-Masspi0){
		//Double_t k = Massomega*Massomega-Masspi0*Masspi0;
		Double_t a1 = (1.+2.*Masselectron*Masselectron/mass/mass);
		Double_t a2 = TMath::Power(1.-4.*Masselectron*Masselectron/mass/mass,0.5);

		return 4.*alpha/3./Pi/mass*GammaToPi0gamma*a1*a2;
	}else 
		return 0;
}
//-------------------------------------------------------
Double_t QEDeePhi(Double_t mass){
	const Double_t GammaToEtagamma = 1.;
	//const Double_t GammaToEtagamma = 0.0131*4.26e-3;

	if(mass>0&&mass<Massphi-Masseta){
		//Double_t k = Massphi*Massomega-Masseta*Masseta;
		Double_t a1 = (1.+2.*Masselectron*Masselectron/mass/mass);
		Double_t a2 = TMath::Power(1.-4.*Masselectron*Masselectron/mass/mass,0.5);

		return 4.*alpha/3./Pi/mass*GammaToEtagamma*a1*a2;
	}else 
		return 0;
}
//-------------------------------------------------------
Double_t PSomega(Double_t mass){
	Double_t k = Massomega*Massomega-Masspi0*Masspi0;

	if(mass>0&&mass<Massomega-Masspi0){
		return TMath::Power((1.+mass*mass/k)*(1.+mass*mass/k)-4.*Massomega*Massomega*mass*mass/k/k,1.5); //PS
	}else 
		return 0;
}

//-------------------------------------------------------
Double_t PSPi0(Double_t mass){
	if(mass>0&&mass<Masspi0){
		return TMath::Power(1.-mass*mass/Masspi0/Masspi0,3);//PS
	}else 
		return 0;

}

//-------------------------------------------------------
Double_t PSEta(Double_t mass){

	if(mass>0&&mass<Masseta){
		return TMath::Power(1.-mass*mass/Masseta/Masseta,3);//PS
	}else 
		return 0;
}

//-------------------------------------------------------
Double_t PSEtaprim(Double_t mass){

	if(mass>0&&mass<Massetaprim){
		return TMath::Power(1.-mass*mass/Massetaprim/Massetaprim,3);//PS
	}else 
		return 0;
}

//-------------------------------------------------------
Double_t PSPhi(Double_t mass){
	Double_t k = Massphi*Massphi-Masseta*Masseta;

	if(mass>0&&mass<Massphi-Masseta){
		return TMath::Power((1.+mass*mass/k)*(1.+mass*mass/k)-4.*Massphi*Massphi*mass*mass/k/k,1.5); //PS
	}else return 0;
}

//-------------------------------------------------------
//Form Factor
Double_t FFSquare(Double_t *x, Double_t *par){
	Double_t m = x[0];
	Double_t LambdaInvSquare = par[0];
	return TMath::Power(1.-m*m*LambdaInvSquare,-2);
}

//-------------------------------------------------------
//Form Factor pi0
Double_t FFSquarePi0(Double_t *x, Double_t *par){
	Double_t m = x[0];
	Double_t LambdaInvSquare = par[0];
	return TMath::Power(1.+m*m*LambdaInvSquare,2);
}

//-------------------------------------------------------
//Form Factor
Double_t FFSquareEtaPrim(Double_t *x){
	TF1 *FEtaprim = new TF1("FEtaprim","1./(pow(1.-[0]*x*x,2)+[1]*[0])",2.*Masselectron,Massphi);//VDM eta' Breit-Wigner fit
	//TF1 *FEtaprim = new TF1("FEtaprim","1./(pow(1.-[0]*x*x,2)+[1]*[0])",2*Masselectron,Massetaprim);//VDM eta' Breit-Wigner fit
	FEtaprim->SetParameters(1.8396,0.01989); //from fit the scanned points.

	Double_t m = x[0];
	return FEtaprim->Eval(m);	
}
//-------------------------------------------------------
Double_t massDistomega(Double_t *x, Double_t *par){
	return FFSquare(x,par)*QEDeeomega(x[0])*PSomega(x[0]);
}

//-------------------------------------------------------
Double_t massDistEtaprim(Double_t *x, Double_t *par){
	return FFSquareEtaPrim(x)*QEDee(x[0])*PSEtaprim(x[0]); //from scanned fit.
	//return FFSquare(x,par)*QEDee(x[0])*PSEtaprim(x[0]);
}

//-------------------------------------------------------
Double_t massDistEta(Double_t *x, Double_t *par){
	return FFSquare(x,par)*QEDee(x[0])*PSEta(x[0]);
}

//-------------------------------------------------------
Double_t massDistPi0(Double_t *x, Double_t *par){
	return FFSquarePi0(x,par)*QEDee(x[0])*PSPi0(x[0]);
}

//-------------------------------------------------------
Double_t massDistPhi(Double_t *x, Double_t *par){
	return FFSquare(x,par)*QEDeePhi(x[0])*PSPhi(x[0]);
}


#endif
