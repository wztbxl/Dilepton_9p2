//////////////////////////////////////////////////////////
// This class is used for meson decay
// ver 1.0 2011/04/27 huangbc
//
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
#include "VecMeson.h"
#include "FunUtil.h"
#include "EffUtil.h"
#include "TTimer.h"

//-------------------------------------------------------
VecMeson::VecMeson(ParticleTypes particle, DecayMode dmode)
{
	if(particle==111) {sprintf(MesonType,"pi0"); 	mParIndex = 0;	mMass=Masspi0;		mWidth=Widthpi0;		}
	if(particle==221) {sprintf(MesonType,"eta");	mParIndex = 1;	mMass=Masseta;		mWidth=Widtheta;		}
	if(particle==113) {sprintf(MesonType,"rho");	mParIndex = 2;	mMass=Massrho;		mWidth=Widthrho;		}
	if(particle==223) {sprintf(MesonType,"omega");	mParIndex = 3;	mMass=Massomega;	mWidth=Widthomega;		}
	if(particle==333) {sprintf(MesonType,"phi");	mParIndex = 4;	mMass=Massphi;		mWidth=Widthphi;		}
	if(particle==331) {sprintf(MesonType,"etaprim");mParIndex = 5;	mMass=Massetaprim;	mWidth=Widthetaprim;	}
	if(particle==443) {sprintf(MesonType,"jpsi");	mParIndex = 6;	mMass=Massjpsi;		mWidth=Widthjpsi;		}
	if(particle==100443) {sprintf(MesonType,"psi");	mParIndex = 7;	mMass=Masspsi;		mWidth=Widthpsi;		}
	if(particle==0)   {sprintf(MesonType,"virtualphoton"); 	mParIndex = 8;	mMass=0.;	mWidth=0.;		}

	//mUseCocktailInput = kTRUE;
	mUseCocktailInput = kFALSE;
	mUseScaleEff = kFALSE;
	mUseTsaPtSpectra = 1;//0 is mT scaling, 1 is Tsallis fit.
	mDmode = dmode;
	mNTrks=100;

	/*
	 * Centrality Sequence:
	 * 0 - 0-80%
	 * 1 - 0-10%
	 * 2 - 10-40%
	 * 3 - 40-80%
	 * 4 - 40-60%
	 * 5 - 60-80%
	 * 6 - 60-70%
	 * 7 - 70-80%
	 */
	mCenIdx = 0;

	mMinRap = -1;
	mMaxRap = 1;

	for(int i=0;i<nSmearFac;i++){
		mPtSmearPar[i]=0.;
	}
	mPtSmearPar[0] = 0.009440;mPtSmearPar[1] = 0.01;mPtSmearPar[2]=0.000511; //from 54 GeV cocktail 
}

VecMeson::~VecMeson(){}

//-------------------------------------------------------
Int_t VecMeson::Init()
{
	cout<<"Initializing meson setting for meson: "<<MesonType<<endl;
	myRandom = new TRandom3();
	// 4 moment for meson, e+, e-, ee pair
	char name[256];

	hSampledPt = new TH1D("hSampledPt","hSampledPt",1000,0,10); 
	hEPSingleTrkEffvsPt = new TH2D("hEPSingleTrkEffvsPt","hEPSingleTrkEffvsPt",300,0,15,200,0,1);
	hEMSingleTrkEffvsPt = new TH2D("hEMSingleTrkEffvsPt","hEMSingleTrkEffvsPt",300,0,15,200,0,1);

	if(mDmode==2) sprintf(name,"hRCPairRapidityvsParentRapidity%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hRCPairRapidityvsParentRapidity%sdalitz",MesonType);
	hRCPairRapidityvsParentRapidity = new TH2D(name,name,400,-2,2,400,-2,2);
	if(mDmode==2) sprintf(name,"hMCPairPtvsParentPt%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hMCPairPtvsParentPt%sdalitz",MesonType);
	hMCPairPtvsParentPt = new TH2D(name,name,300,0,15,300,0,15);
	if(mDmode==2) sprintf(name,"hMCAcc0PairPtvsParentPt%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hMCAcc0PairPtvsParentPt%sdalitz",MesonType);
	hMCAcc0PairPtvsParentPt = new TH2D(name,name,300,0,15,300,0,15);
	if(mDmode==2) sprintf(name,"hMCAcc1PairPtvsParentPt%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hMCAcc1PairPtvsParentPt%sdalitz",MesonType);
	hMCAcc1PairPtvsParentPt = new TH2D(name,name,300,0,15,300,0,15);

	int npTBins = 200;
	int pTMin = 0;
	int pTMax = 5;
	int nMassBins = 1800;
	double MassMin = 0;
	double MassMax = 3.6;
	if(mDmode==2) sprintf(name,"hMCMvsPt%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hMCMvsPt%sdalitz",MesonType);
	hMCMvsPt = new TH2D(name,name,1500,0,15,3000,0,6);
	if(mDmode==2) sprintf(name,"hMCAcc0MvsPt%s2ee",MesonType); // Acc0 for the cut
	if(mDmode==3) sprintf(name,"hMCAcc0MvsPt%sdalitz",MesonType);
	hMCAcc0MvsPt = new TH2D(name,name,1500,0,15,3000,0,6);
	if(mDmode==2) sprintf(name,"hMCAcc1MvsPt%s2ee",MesonType); // Acc1 for the cut
	if(mDmode==3) sprintf(name,"hMCAcc1MvsPt%sdalitz",MesonType);
	hMCAcc1MvsPt = new TH2D(name,name,1500,0,15,3000,0,6);
	if(mDmode==2) sprintf(name,"hMCAcc2MvsPt%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hMCAcc2MvsPt%sdalitz",MesonType);
	hMCAcc2MvsPt = new TH2D(name,name,1500,0,15,3000,0,6);
	if(mDmode==2) sprintf(name,"hRCAcc1MvsPt3D%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hRCAcc1MvsPt3D%sdalitz",MesonType);
	hRCAcc1MvsPt3D = new TH2D(name,name,1500,0,15,3000,0,6);
	hRCAcc1MvsPt3D->Sumw2();
	if(mDmode==2) sprintf(name,"hRCAcc2MvsPt3D%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hRCAcc2MvsPt3D%sdalitz",MesonType);
	hRCAcc2MvsPt3D = new TH2D(name,name,1500,0,15,3000,0,6);
	hRCAcc2MvsPt3D->Sumw2();
	if(mDmode==2) sprintf(name,"hMCAcc1MvsPtwoSmear%s2ee",MesonType);
	if(mDmode==3) sprintf(name,"hMCAcc1MvsPtwoSmear%sdalitz",MesonType);
	hMCAcc1MvsPtwoSmear = new TH2D(name,name,1500,0,15,3000,0,6);
	hMCAcc1MvsPtwoSmear->Sumw2();
	// if(mDmode==2) sprintf(name,"hMCMvsPt%s2ee",MesonType);
	// if(mDmode==3) sprintf(name,"hMCMvsPt%sdalitz",MesonType);
	// hMCMvsPt = new TH2D(name,name,nMassBins,MassMin,MassMax,npTBins,pTMin,pTMax);
	// if(mDmode==2) sprintf(name,"hMCAcc0MvsPt%s2ee",MesonType);
	// if(mDmode==3) sprintf(name,"hMCAcc0MvsPt%sdalitz",MesonType);
	// hMCAcc0MvsPt = new TH2D(name,name,nMassBins,MassMin,MassMax,npTBins,pTMin,pTMax);
	// if(mDmode==2) sprintf(name,"hMCAcc1MvsPt%s2ee",MesonType);
	// if(mDmode==3) sprintf(name,"hMCAcc1MvsPt%sdalitz",MesonType);
	// hMCAcc1MvsPt = new TH2D(name,name,nMassBins,MassMin,MassMax,npTBins,pTMin,pTMax);
	// if(mDmode==2) sprintf(name,"hMCAcc2MvsPt%s2ee",MesonType);
	// if(mDmode==3) sprintf(name,"hMCAcc2MvsPt%sdalitz",MesonType);
	// hMCAcc2MvsPt = new TH2D(name,name,nMassBins,MassMin,MassMax,npTBins,pTMin,pTMax);
	// if(mDmode==2) sprintf(name,"hRCAcc1MvsPt3D%s2ee",MesonType);
	// if(mDmode==3) sprintf(name,"hRCAcc1MvsPt3D%sdalitz",MesonType);
	// hRCAcc1MvsPt3D = new TH2D(name,name,nMassBins,MassMin,MassMax,npTBins,pTMin,pTMax);
	// hRCAcc1MvsPt3D->Sumw2();
	// if(mDmode==2) sprintf(name,"hRCAcc2MvsPt3D%s2ee",MesonType);
	// if(mDmode==3) sprintf(name,"hRCAcc2MvsPt3D%sdalitz",MesonType);
	// hRCAcc2MvsPt3D = new TH2D(name,name,nMassBins,MassMin,MassMax,npTBins,pTMin,pTMax);
	// hRCAcc2MvsPt3D->Sumw2();
	// if(mDmode==2) sprintf(name,"hMCAcc1MvsPtwoSmear%s2ee",MesonType);
	// if(mDmode==3) sprintf(name,"hMCAcc1MvsPtwoSmear%sdalitz",MesonType);
	// hMCAcc1MvsPtwoSmear = new TH2D(name,name,nMassBins,MassMin,MassMax,npTBins,pTMin,pTMax);
	// hMCAcc1MvsPtwoSmear->Sumw2();

	hMCAcc1PairRapidity = new TH1D("hMCAcc1PairRapidity","hMCAcc1PairRapidity",400,-2,2);
	hMCAcc1PairEMPtvsEPPt = new TH2D("hMCAcc1PairEMPtvsEPPt","hMCAcc1PairEMPtvsEPPt",1000,0,10,1000,0,10);
	hPlusSmearvsPt = new TH2D("hPlusSmearvsPt","hPlusSmearvsPt",1500,0,15,600,-3,3);
	hMinusSmearvsPt = new TH2D("hMniusSmearvsPt","hMniusSmearvsPt",1500,0,15,600,-3,3);
	RapiditySTARAcc = new TH1D("RapiditySTARAcc","RapiditySTARAcc",240,-1.2,1.2);
	RapiditywoSTARAcc = new TH1D("RapiditywoSTARAcc",";Y;Counts",240,-1.2,1.2);
	RapidityOutSTARAcc = new TH1D("RapidityOutSTARAcc",";Y;Counts",240,-1.2,1.2);
	MeePtRapidityOver1 = new TH2D("MeePtRapidityOver1",";M_{ee};p_{T}",1000,0,10,3000,0,6);
	MeeFullRapidity = new TH1D("MeeFullRapidity",";M_{ee};Counts",3000,0,6);
	MeeWopTCut = new TH1D("MeeWopTCut",";M_{ee};Counts",3000,0,6);
	MeeWopTEtaCut = new TH1D("MeeWopTEtaCut",";M_{ee};Counts",3000,0,6);
	MeeWoEtaCut = new TH1D("MeeWoEtaCut",";M_{ee};Counts",3000,0,6);
	MeeWpTEtaWoRapidity = new TH1D("MeeWpTEtaWoRapidity",";M_{ee};Counts",3000,0,6);
	MeeWpTEtaWRapidity = new TH1D("MeeWpTEtaWRapidity",";M_{ee};Counts",3000,0,6);
	pTWithLargeMee = new TH2D("pTWithLargeMee",";p_{T};M_{ee}",1000,0,10,3000,0,6);
	PtMVsPtPLargeMass = new TH2D("PtMVsPtPLargeMass",";p_{T}^{e+};p_{T}^{e-}",1000,0,10,1000,0,10);
	EtaMVsEtaPLargeMass = new TH2D("EtaMVsEtaPLargeMass",";#eta_{e+};#eta_{e-}",240,-1.2,1.2,240,-1.2,1.2);
	CutRecorder = new TH1D("CutRecorder",";;Counts",15,0,15);
	CutRecorder->GetXaxis()->SetBinLabel(1,"N_{M > M_{#pi}}");
	CutRecorder->GetXaxis()->SetBinLabel(2,"Out of STAR Acc.");
	CutRecorder->GetXaxis()->SetBinLabel(4,"Out of p_{T}^{e+}");
	CutRecorder->GetXaxis()->SetBinLabel(5,"Out of p_{T}^{e-}");
	CutRecorder->GetXaxis()->SetBinLabel(6,"Out of #eta_{e+}");
	CutRecorder->GetXaxis()->SetBinLabel(7,"Out of #eta_{e-}");


	//****** pt spectra input ******
	//mT scaling for PHENIX
	TF1 *pi0 = new TF1("pi0","x*[0]*pow((exp(-([1]*x+[2]*pow(x,2)))+x/[3]),-[4])",0,10);
	pi0->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *eta = new TF1("eta","x*[0]*pow((exp(-([1]*sqrt(pow(0.5478,2)-pow(0.1350,2)+x*x)+[2]*pow(sqrt(pow(0.5478,2)-pow(0.1350,2)+x*x),2)))+sqrt(pow(0.5478,2)-pow(0.1350,2)+x*x)/[3]),-[4])",0,10);
	eta->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *etaprim = new TF1("etaprim","x*[0]*pow((exp(-([1]*sqrt(pow(0.9578,2)-pow(0.1350,2)+x*x)+[2]*pow(sqrt(pow(0.9578,2)-pow(0.1350,2)+x*x),2)))+sqrt(pow(0.9578,2)-pow(0.1350,2)+x*x)/[3]),-[4])",0,10);
	etaprim->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *omega = new TF1("omega","x*[0]*pow((exp(-([1]*sqrt(pow(0.7826,2)-pow(0.1350,2)+x*x)+[2]*pow(sqrt(pow(0.7826,2)-pow(0.1350,2)+x*x),2)))+sqrt(pow(0.7826,2)-pow(0.1350,2)+x*x)/[3]),-[4])",0,10);
	omega->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *phi = new TF1("phi","x*[0]*pow((exp(-([1]*sqrt(pow(1.019,2)-pow(0.1350,2)+x*x)+[2]*pow(sqrt(pow(1.019,2)-pow(0.1350,2)+x*x),2)))+sqrt(pow(1.019,2)-pow(0.1350,2)+x*x)/[3]),-[4])",0,10);
	phi->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *jpsi = new TF1("jpsi","x*[0]*pow((exp(-([1]*sqrt(pow(3.097,2)-pow(0.1350,2)+x*x)+[2]*pow(sqrt(pow(3.097,2)-pow(0.1350,2)+x*x),2)))+sqrt(pow(3.097,2)-pow(0.1350,2)+x*x)/[3]),-[4])",0,10);
	jpsi->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *psi = new TF1("psi","x*[0]*pow((exp(-([1]*sqrt(pow(3.686,2)-pow(0.1350,2)+x*x)+[2]*pow(sqrt(pow(3.686,2)-pow(0.1350,2)+x*x),2)))+sqrt(pow(3.686,2)-pow(0.1350,2)+x*x)/[3]),-[4])",0,10);
	psi->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *rho = new TF1("rho","x*[0]*pow((exp(-([1]*sqrt(pow(0.7755,2)-pow(0.1350,2)+x*x)+[2]*pow(sqrt(pow(0.7755,2)-pow(0.1350,2)+x*x),2)))+sqrt(pow(0.7755,2)-pow(0.1350,2)+x*x)/[3]),-[4])",0,10);
	rho->SetParameters(504.5,0.52,0.16,0.7,8.27);

	TF1 *virtualphotonpt = new TF1("virtualphotonpt","[0]",0,10.); //inv.yield*pt
	virtualphotonpt->SetParameter(0,1.);
	virtualphotonpt->SetNpx(10000);

	if(mParIndex==0) funMeson = (TF1 *)pi0; //pi0
	if(mParIndex==1) funMeson = (TF1 *)eta; //eta
	if(mParIndex==2) funMeson = (TF1 *)rho; //rho
	if(mParIndex==3) funMeson = (TF1 *)omega; //omega
	if(mParIndex==4) funMeson = (TF1 *)phi; //phi
	if(mParIndex==5) funMeson = (TF1 *)etaprim; //etaprim
	if(mParIndex==6) funMeson = (TF1 *)jpsi; //jpsi
	if(mParIndex==7) funMeson = (TF1 *)psi; //psi
	if(mParIndex==8) funMeson = (TF1 *)virtualphotonpt; //virtual photon, flat pt

	//-------------------------- now work at here, you need to think about how to change it for different energy -----------------------------//
	// load inv.yields pt spectra fit parameters
	//TFile *FTS  = new TFile(Form("inputFile/AuAu200_inputpT_Cen%d_%d.root",CentralityLow[mCenIdx],CentralityHi[mCenIdx])); 
	// 2024.03.12 now using the some of 54GeV files as input
	//inputFile is my path to cocktail
	TFile *FTS;
	if(mCenIdx>=6) {
		FTS = TFile::Open("./inputFile/pTTBW4Cocktail_54.root");  //temporary use 60-80% pt shape for 60-70%, 70-80%
																//because 62 GeV data is not divide the Centrality so use minbias data  temporary
	}
	else {
		FTS = TFile::Open("./inputFile/pTTBW4Cocktail_54.root"); 
	}
	FTS->Print();
	if (!FTS) cout << "Fail to open pT file" << endl;

	TFile *ftsa = TFile::Open("./inputFile/mesons_baryons_noOmega_080.root");
	TFile *phoInput  = TFile::Open("./inputFile/AuAu200_inputpT_Cen60_80.root");
	TFile *JpsiInput = TFile::Open("./inputFile/JPsiPtSpectra_62.root");

	//tsallis input
	//use 62.4 GeV input now and there only have the minbias data
	if(mParIndex==0) {//use same distribution in different centrality ,how to get the different distrinbution in different centriality use Qian`s method in different centrality 
		if(mCenIdx>=4) histMeson = (TH1D *)FTS->Get("pi0_to_gamma_eepT"); //40-60%, 60-80% pi0
		else 	       histMeson = (TH1D *)FTS->Get("His54GeV_pi0_to_gamma_eepT_Func"); //0-80%, 0-10%, 10-40%, and 40-80% pi0
	}
	if(mParIndex==1)   histMeson = (TH1D *)FTS->Get("His54GeV_eta_to_gamma_eepT_Func"); //eta
	if(mParIndex==2) {
		if(mCenIdx>=4) histMeson = (TH1D *)phoInput->Get("pT_rho"); //40-60%, 60-80% rho
		else           histMeson = (TH1D *)ftsa->Get("hFit22"); // 0-80% rho; For other centrality which has no rho spectrum, using the MB rho pT distribution
	}
	if(mParIndex==3)   histMeson = (TH1D *)FTS->Get("His54GeV_omega_to_eepT_Func"); //omega
	if(mParIndex==4)   histMeson = (TH1D *)FTS->Get("His54GeV_phi_to_eepT_Func"); //phi
	if(mParIndex==5)   histMeson = (TH1D *)FTS->Get("His54GeV_etaprime_to_gamma_eepT_Func"); //etaprim
	if(mParIndex==6)   histMeson = (TH1D *)JpsiInput->Get("hJPsidNdpT"); //jpsi//in Yi`s code what`s different in two file Jpsi spectra ?
	if(mParIndex==7){
		if(mCenIdx==0) histMeson = (TH1D *)phoInput->Get("pT_Psi"); //0-80% psi //use 200 GeV data
		else           histMeson = (TH1D *)JpsiInput->Get("hJPsidNdpT"); //psi - do not have psi pt shape for different centralities, using jpsi pt shape for psi
	}//for all the particles who decide the pT shape and what`s the different between the shape and yelid?//zhen`s question
	//What`s differnt about Jpsi and Psi?
	if(mParIndex==8)   histMeson = (TH1D *)virtualphotonpt->GetHistogram(); //virtual photon, flat pt

	if (!histMeson) cout << " Fail to get histogram " << endl;
	//histMeson->Print();
	cout << "mCenIdx = " << mCenIdx << " CentralityLow[mCenIdx] = " <<CentralityLow[mCenIdx] << " CentralityHi[mCenIdx] = " << CentralityHi[mCenIdx] << endl;
	TFile *fcocktail = TFile::Open(Form("/star/u/wangzhen/QA/wangzhen/Cocktail/genCocktail/output_old/cen%d%d_cocktail_withoutRho.root",CentralityLow[mCenIdx],CentralityHi[mCenIdx]));
	hCocktail = (TH2D *)fcocktail->Get("hMCAcc0MvsPt");// |Y_{ee}|<1

	//****** Smear for MC ******
	funSmearPt = new TF1("funSmearPt","sqrt([0]*[0]*x*x+[1]*[1])",0.,10.);
	funSmearPt->SetParameters(mPtSmearPar);
	funSmearPtEmb = new TF1("funSmearPtEmb","sqrt([0]*[0]*x*x+[1]*[1])",0,10);
	// funSmearPtEmb->SetParameters(0.005690,0.007473);
	funSmearPtEmb->SetParameters(0.00700,0.007473);
	//how to do the smear? and why we should do the smear?
	//****** pT Res From embedding ********
	TFile *pTResInput = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/pTresFromEmbedding.root");
	TFile *pTResFunInput = new TFile("/star/u/wangzhen/QA/wangzhen/Cocktail/DCBFunction.root");
	PtRes2D = (TH2D*)pTResInput->Get("hPtResvsPtQ_zy");
	PtRes2D->Print();
	for(int i = 0;i<19;i++)
        {
                int PtBinLow = PtRes2D->GetXaxis()->FindBin((i+1)*0.2+1.e-4);
                int PtBinHigh = PtRes2D->GetXaxis()->FindBin((i+2)*0.2-1.e-4);
                pTRes1D[i] = PtRes2D->ProjectionY(Form("pTRes1D_%d",i),PtBinLow,PtBinHigh);
				//pTResFun[i] = (TF1*)pTResFunInput->Get(Form("DCBpTBin%d",i));
				pTResFun[i] = new TF1(Form("DCBpTBin%d",i),DoubleCrystalBall,-2,2,7);
				for (int j = 0 ; j<7;j++)
				{
					pTResFun[i]->SetParameter(j,pTResPar[i][j]);
				}
		
		//cout<<1<<endl;
		pTRes1D[i]->Print();
		pTResFun[i]->Print();
        }//0.2 step is beacuse the embedding data is a small sample;

	momShape = new TF1("momShape",CrystalBall2,-1.,1.,7);
	//momShape->SetParameters(1., -1e-3, 0.01, 1.16, 1.84, 2.18, 1.92);//-3.9e-4
	momShape->SetParameters(1., -1e-3, 0.01, 1.29, 1.75, 2.92, 1.84);
	momShape->SetNpx(10000);

	//****** Efficiency for single track ******
	InitialzeEffHist(mCenIdx);

	//NA49
	//fRapidity = new TF1("fRapidity","[0]*(exp(-pow(x-[1],2)/2./[2]/[2])+exp(-pow(x+[1],2)/2./[2]/[2]))",-5,5);
	//fRapidity->SetParameters(107.6,0.72,1.18); //pi-

	//CERES Pb+Au
	fRapidity = new TF1("fRapidity","pow(cosh(3*x/4./sqrt(log([0]/2./[1]))/(1-pow(x,2)/2./[0]*[2])),-2.)",-5,5);
	fRapidity->SetParameter(0,9.2);//sqrt(s)//set as my energy
	fRapidity->SetParameter(1,0.938); // nucleon mass
	fRapidity->SetParameter(2,mMass); // meson mass
	if(mParIndex==8) fRapidity->SetParameter(2,Masspi0);
	fRapidity->SetRange(-1.,1.);

	//****** soure mass shape ******
	massfunMeson = new TF1("massfunMeson","2.*[0]*[1]/(pow((x-[2]),2)+pow(([1]/2),2))",2.*Masselectron,6);
	massfunMeson->SetParameters(1.,mWidth,mMass);
	massfunMeson->SetNpx(10000);

	if(mParIndex==2){
		massfunrho = new TF1("massfunrho",massDistRho,2.*Masselectron,6,2);
		massfunrho->SetParameter(0,0.5);
		massfunrho->SetParameter(1,0.16);
		massfunrho->SetNpx(10000);
	}

	if(mDmode==3){
		if(mParIndex==0){
			massfunee = new TF1("massfunee",massDistPi0,2.*Masselectron,Masspi0,1);
			// massfunee->SetParameter(0,1.756);
			massfunee->SetParameter(0,1.838);
		}
		if(mParIndex==1){
			massfunee = new TF1("massfunee",massDistEta,2.*Masselectron,Masseta,1);
			massfunee->SetParameter(0,1.95);
		}
		if(mParIndex==3){
			massfunee = new TF1("massfunee",massDistomega,2.*Masselectron,Massomega,1);
			massfunee->SetParameter(0,2.24);
		}
		if(mParIndex==4){
			massfunee = new TF1("massfunee",massDistPhi,2.*Masselectron,Massphi,1);
			massfunee->SetParameter(0,3.8);//3.8+-1.8 Physics Letters B 504 (2001) 275-281
			//massfunee->SetParameter(0,3.23);//Yifei
		}
		if(mParIndex==5){
			massfunee = new TF1("massfunee",massDistEtaprim,2.*Masselectron,Massetaprim,1);
			massfunee->SetParameter(0,1.7);
		}

		if(!massfunee){
			cout<<"ERROR: ee mass function can't be initialized, please check particle decay mode!"<<endl;
		}else{
			massfunee->SetNpx(10000);
			hmassfunee = (TH1D*)massfunee->GetHistogram();
		}
	}	

	cout<<"Initilized!"<<endl;
	return 1;
}

//-------------------------------------------------------
Double_t VecMeson::GetSmear(Double_t pT)
{
        double smear;
        int pTBin = int(pT/0.2);
        //cout <<"pT = "<<pT<<" pTBin = "<<pTBin<<endl;
        pTBin = pTBin-1;
        if (pTBin>0 && pTBin <=19)
                smear = pTRes1D[pTBin-1]->GetRandom();
        else if (pTBin <= 0) smear = pTRes1D[0]->GetRandom();
        else if (pTBin >=20) smear = pTRes1D[18]->GetRandom();

        return smear;
}
Double_t VecMeson::GetSmear2(Double_t pT)
{
        double smear;
        int pTBin = int(pT/0.2);
        //cout <<"pT = "<<pT<<" pTBin = "<<pTBin<<endl;
        pTBin = pTBin-1;
        if (pTBin>0 && pTBin <=19)
                smear = pTResFun[pTBin-1]->GetRandom();
        else if (pTBin <= 0) smear = pTResFun[0]->GetRandom();
        else if (pTBin >=20) smear = pTResFun[18]->GetRandom();

        return smear;
}


//-------------------------------------------------------
TLorentzVector VecMeson::myBoost(TLorentzVector parent, TLorentzVector daughter)
{
	float betax = parent.Px()/parent.E();
	float betay = parent.Py()/parent.E();
	float betaz = parent.Pz()/parent.E();
	daughter.Boost(betax,betay,betaz);
	return daughter;
}

//-------------------------------------------------------
TLorentzVector VecMeson::twoBodyDecay(TLorentzVector parent, Double_t dmass) {

	Double_t e = parent.M()/2.;
	Double_t p = sqrt(e*e-dmass*dmass);
	Double_t costheta = myRandom->Uniform(-1.,1.);
	Double_t phi = myRandom->Uniform(0,TMath::Pi()*2);
	Double_t pz = p*costheta;
	Double_t px = p*sqrt(1.-costheta*costheta)*cos(phi);
	Double_t py = p*sqrt(1.-costheta*costheta)*sin(phi);
	TLorentzVector daughter(px,py,pz,e);
	//cout<<px<<" "<<py<<"  "<<pz<<endl;
	return myBoost(parent,daughter);
}

//-------------------------------------------------------
Double_t VecMeson::SampleMesonMass()
{
	Double_t mass = 0.;
	if(mDmode==2){
		if(mParIndex!=2 && mParIndex!=8){ 
			//mass = mMass + myRandom->Gaus(0.,mWidth);
			mass = massfunMeson->GetRandom();
		}else if(mParIndex==2){
			mass = massfunrho->GetRandom();
		}else if(mParIndex==8){
			mass = myRandom->Uniform(2.*Masselectron,6.);
		}
	}

	if(mDmode==3){
		mass = mMass;
	}

	return mass;
}

//-------------------------------------------------------
void VecMeson::Rot(Double_t pin[3], Double_t pout[3], Double_t costheta, Double_t sintheta, Double_t cosphi, Double_t sinphi)
{
	pout[0] = pin[0]*costheta*cosphi-pin[1]*sinphi+pin[2]*sintheta*cosphi;
	pout[1] = pin[0]*costheta*sinphi+pin[1]*cosphi+pin[2]*sintheta*sinphi;
	pout[2] = -1.0  * pin[0] * sintheta + pin[2] * costheta;
	return;
}

//-------------------------------------------------------
void VecMeson::DalitzDecay(TLorentzVector parent)
{
	Double_t eeMass;
	Double_t e1, p1, e3, p3;
	Double_t betaSquare, lambda;
	Double_t costheta, sintheta, cosphi, sinphi, phi;

	Double_t NeMass = 0.;
	if(mParIndex==0||mParIndex==1||mParIndex==5) NeMass=0.;
	else if(mParIndex==3) NeMass=Masspi0;
	else if(mParIndex==4) NeMass=Masseta;
	else{ cout<<" Wrong input of dalitz decay particle"<<endl; return;}

	// Sample the lepton pair mass from function 
	//massfunee->SetRange(2.*Masselectron,mMass-NeMass);
	hmassfunee->SetAxisRange(2.*Masselectron,mMass-NeMass,"X");
	for( ;; ) {

		//eeMass = massfunee->GetRandom();
		eeMass = hmassfunee->GetRandom();
		//cout<<" random mass:"<<eeMass<<endl;

		if ( mMass - NeMass > eeMass && eeMass / 2. > Masselectron ) break; //energy conservation
	}

	// lepton pair kinematics in virtual photon rest frame
	e1 = eeMass / 2.;
	p1 = TMath::Sqrt((e1 + Masselectron) * (e1 - Masselectron));
	betaSquare = 1.0 - 4.0 * (Masselectron * Masselectron) / (eeMass * eeMass);
	lambda      = betaSquare / (2.0 - betaSquare);//what`s this? lorenzt transform?

	if ( NeMass<0.01 ){
		do{
			costheta = (2.0*myRandom->Rndm())-1.;
		}while ( (1.0+lambda*costheta*costheta)<(2.0*gRandom->Rndm()) );
	} else{
		costheta = (2.0*myRandom->Rndm())-1.;
	} //for dalitz decay (gamma+e+e) polarization
	// caculate by yourself

	//costheta = (2.0*myRandom->Rndm())-1.;

	sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
	phi      = 2.0 * TMath::ACos(-1.) * myRandom->Rndm();
	sinphi   = TMath::Sin(phi);
	cosphi   = TMath::Cos(phi);

	// momentum vectors of leptons in virtual photon rest frame
	Double_t pProd1[3] = {p1 * sintheta * cosphi, 
		p1 * sintheta * sinphi, 
		p1 * costheta};
	Double_t pProd2[3] = {-1.0 * p1 * sintheta * cosphi, 
		-1.0 * p1 * sintheta * sinphi, 
		-1.0 * p1 * costheta};

	// netural particle kinematics in meson rest frame
	e3       = (mMass * mMass + NeMass * NeMass - eeMass * eeMass)/(2. * mMass);
	p3       = TMath::Sqrt((e3 + NeMass)  * (e3 - NeMass));
	costheta = (2.0 * myRandom->Rndm()) - 1.;
	sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
	phi      = 2.0 * TMath::ACos(-1.) * myRandom->Rndm();
	sinphi   = TMath::Sin(phi);
	cosphi   = TMath::Cos(phi);       

	// netural particle 4-vector in meson rest frame
	fProducts[2].SetPx(p3 * sintheta * cosphi);
	fProducts[2].SetPy(p3 * sintheta * sinphi);
	fProducts[2].SetPz(p3 * costheta);
	fProducts[2].SetE(e3);

	// lepton 4-vectors in properly rotated virtual photon rest frame//rotated in the VP driection in the meson rest frame
	Double_t pRot1[3] = {0.};
	Rot(pProd1, pRot1, costheta, -sintheta, -cosphi, -sinphi);
	Double_t pRot2[3] = {0.};
	Rot(pProd2, pRot2, costheta, -sintheta, -cosphi, -sinphi); 
	fProducts[0].SetPx(pRot1[0]);
	fProducts[0].SetPy(pRot1[1]);
	fProducts[0].SetPz(pRot1[2]);
	fProducts[0].SetE(e1);
	fProducts[1].SetPx(pRot2[0]);
	fProducts[1].SetPy(pRot2[1]);
	fProducts[1].SetPz(pRot2[2]);
	fProducts[1].SetE(e1);

	// boost the dilepton into the meson's rest frame 
	Double_t eLPparent = TMath::Sqrt(p3 * p3 + eeMass * eeMass);
	TVector3 boostPair( -1.0 * fProducts[2].Px() / eLPparent, 
			-1.0 * fProducts[2].Py() / eLPparent,
			-1.0 * fProducts[2].Pz() / eLPparent);
	fProducts[0].Boost(boostPair);
	fProducts[1].Boost(boostPair);

	// boost all decay products into the lab frame 
	TVector3 boostLab( parent.Px() / parent.E(), 
			parent.Py() / parent.E(),
			parent.Pz() / parent.E());
	fProducts[0].Boost(boostLab);
	fProducts[1].Boost(boostLab);
	fProducts[2].Boost(boostLab);

	return;
}

//-------------------------------------------------------
Double_t VecMeson::EvalEff3D(TLorentzVector electron,int ipttpc,int ietatpc,int iphitpc,int ipttof,int ietatof,int iphitof,int charge)
{
	Double_t eff = 0.;
	Double_t pt = electron.Pt();
	Double_t p = electron.P();
	/*if (charge == 1)
		 eff = hEff_Tpc_Pos[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*hEff_Tof_Pos[ietatof][iphitof]->GetBinContent(ipttof+1)*fEPiEff_Tof->Eval(pt)*funnSigEffP->Eval(p)*funBetaEffP->Eval(p)*funndEdxEffPt->Eval(pt);
	else eff = hEff_Tpc_Neg[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*hEff_Tof_Neg[ietatof][iphitof]->GetBinContent(ipttof+1)*fEPiEff_Tof->Eval(pt)*funnSigEffP->Eval(p)*funBetaEffP->Eval(p)*funndEdxEffPt->Eval(pt);
	*/
	if (charge == 1)
	{
	eff = getEff(pt,hMBEff_Tof_Pos[ietatof][iphitof],ptl_Tof,pth_Tof)*f_ElecPoionRatio->Eval(pt);
	if(mDebug) cout << "Tof eff = " << eff << endl;
	eff = eff*getEff(pt,hMBEff_Tpc_Pos[ietatpc][iphitpc],ptl_Tpc,pth_Tpc);//cout<<"TPC"<<endl;
	if(mDebug) cout << "Tof*TPC eff = " << eff << endl;
	eff = eff*->GetParameter(0);//cout<<"beta"<<endl;
	if(mDebug) cout << "Tof*TPC*beta eff = " << eff << endl;
	if (pt > 0.8)
	{
		eff = eff*f_nSigmaEEff_HigpT->Eval(pt);
	if(mDebug) cout << "Tof*TPC*beta*nSigmaE eff = " << eff << endl;
	} else eff = eff*f_nSigmaEEff_lowpT->Eval(pt);
	if(mDebug) cout << "Tof*TPC*beta*nSigmaE eff = " << eff << endl;
	//cout<<"sigma"<<endl;//TOF*TPC*Beta*nSigmaE
	}
	else {
		eff = getEff(pt,hMBEff_Tof_Neg[ietatof][iphitof],ptl_Tof,pth_Tof)*f_ElecPoionRatio->Eval(pt);
		if(mDebug) cout << "Tof eff = " << eff << endl;
		eff = eff*getEff(pt,hMBEff_Tpc_Neg[ietatpc][iphitpc],ptl_Tpc,pth_Tpc);//cout<<"TPC"<<endl;
		if(mDebug) cout << "Tof*TPC eff = " << eff << endl;
		eff = eff*f_betaCutEff->GetParameter(0);//cout<<"beta"<<endl;
		if(mDebug) cout << "Tof*TPC*beta eff = " << eff << endl;
		if (pt > 0.8)
		{
			eff = eff*f_nSigmaEEff_HigpT->Eval(pt);
			if(mDebug) cout << "Tof*TPC*beta*nSigmaE eff = " << eff << endl;
		} else eff = eff*f_nSigmaEEff_lowpT->Eval(pt);
		if(mDebug) cout << "Tof*TPC*beta*nSigmaE eff = " << eff << endl;
		//cout<<"sigma"<<endl;//TOF*TPC*Beta*nSigmaE

	}
	return eff;
}

//-------------------------------------------------------
Double_t VecMeson::EvalScaleEff3D(TLorentzVector electron,int ipttpc,int ietatpc,int iphitpc,int ipttof,int ietatof,int iphitof,int charge)
{
	Double_t eff = 0.;
	Double_t pt = electron.Pt();
	Double_t p = electron.P();
	if (charge == 1)
		 eff = hMBEff_Tpc_Pos[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*fTpcEffRatioToMB->Eval(pt)*hMBEff_Tof_Pos[ietatof][iphitof]->GetBinContent(ipttof+1)*fMBEPiEff_Tof->Eval(pt)*fTofEffRatioToMB->Eval(pt)*funnSigEffP->Eval(p)*funBetaEffP->Eval(p)*funndEdxEffPt->Eval(pt);
	else eff = hMBEff_Tpc_Neg[ietatpc][iphitpc]->GetBinContent(ipttpc+1)*fTpcEffRatioToMB->Eval(pt)*hMBEff_Tof_Neg[ietatof][iphitof]->GetBinContent(ipttof+1)*fMBEPiEff_Tof->Eval(pt)*fTofEffRatioToMB->Eval(pt)*funnSigEffP->Eval(p)*funBetaEffP->Eval(p)*funndEdxEffPt->Eval(pt);
	return eff;
}

//-------------------------------------------------------
void VecMeson::GenerateDecay()
{
	TTimer *timer = new TTimer();

	for(int i=0;i<mNTrks;i++)
	{
		if(i%(mNTrks/10)==0) cout<<"processing "<<i<<" tracks..."<<endl;
		if(mDebug) cout << i << endl;
		// cout << i << endl;

		if((i%10000)==0){
			long long tmp = (long long)timer->GetAbsTime();
			unsigned int stime = abs(tmp)/myRandom->Rndm() ;
			//cout << " time is " << stime <<endl;
			gRandom->SetSeed(stime);
		}

		if(mDebug) cout << "before pass sample rho phi" << endl;
		//Double_t rap = fRapidity->GetRandom();
		Double_t rap = myRandom->Uniform(mMinRap, mMaxRap);
		Double_t phi = myRandom->Uniform(0.,2.*TMath::Pi());//for the mother particle?

		if(mDebug) cout << "pass sample rho phi" << endl;
		Double_t pt;
		if(mUseTsaPtSpectra == 1)      pt = histMeson->GetRandom();
		else if(mUseTsaPtSpectra == 0) pt = funMeson->GetRandom();
		else{
			cout<<"The input pT spectrum type is not right !"<<endl;
			return;
		}

		//Double_t smass = mMass + myRandom->Gaus(0.,mWidth);
		if(mParIndex==2) massfunrho->SetParameter(0,pt);
		Double_t smass = SampleMesonMass();//smear the width//???

		if(mParIndex==8 && mUseCocktailInput){ //use cocktail to be the virtual photon pt&mass input
			hCocktail->GetRandom2(pt, smass);
			if(smass<=2.*Masselectron) continue;
		}

		hSampledPt->Fill(pt);

		Double_t mt = sqrt(pt*pt+smass*smass);
		Double_t pz = mt*TMath::SinH(rap);
		Double_t ptot = sqrt(pt*pt+pz*pz);
		if(ptot==0){
			cout<<"p==0 continue"<<endl;
			continue;
		}
		Double_t eta = 0.5*log((ptot+pz)/(ptot-pz+1e-20));

		TLorentzVector parent;
		parent.SetPtEtaPhiM(pt,eta,phi,smass);

		TLorentzVector daughterP(0,0,0,0);
		TLorentzVector daughterN(0,0,0,0);
		if(mDmode==2){
			// two body for omega, phi, jpsi
			daughterP = twoBodyDecay(parent,Masselectron);
			daughterN = parent-daughterP;
		}
		if(mDmode==3){
			// dalitz   for pi0, eta, etaprim, omega, phi
			DalitzDecay(parent);
			daughterP = fProducts[0];
			daughterN = fProducts[1];
		}
		if(mDebug) cout<<"pt: "<<pt<<"   \tdaughter_pt:"<<daughterP.X()<<"  "<<daughterP.Y()<<"  "<<daughterP.Z()<<endl;

		if(mDebug) cout<<1<<endl;
		Double_t eppt  = daughterP.Pt();
		Double_t epeta = daughterP.Eta();
		Double_t epphi = daughterP.Phi();
		Double_t epsig = funSmearPt->Eval(eppt);
		Double_t epsigEmb = funSmearPtEmb->Eval(eppt);
		TVector3 eP(0,0,0);
		eP.SetPtEtaPhi(eppt,epeta,epphi);
		Double_t smeppt1  = eP.Pt()*(1. + momShape->GetRandom()*epsig/0.01);
		
 		// double smear = GetSmear(eP.Pt());
		double smear = GetSmear2(eP.Pt()); // new Smear for Function
		//cout<<epsig/epsigEmb<<endl; 
		// Double_t smeppt = eP.Pt();// try to remove pt smearing
		Double_t smeppt = eP.Pt()*(1.+smear*epsig/epsigEmb);
		double delatpT = smeppt-eP.Pt();
		hPlusSmearvsPt->Fill(smeppt, delatpT/smeppt);
		
		if(mDebug) cout<<smeppt1<<"/"<<smeppt<<endl;
		//eP.SetPtEtaPhi(smeppt,epeta,epphi);
		Double_t smepeta = eP.Eta();
		Double_t smepphi = eP.Phi();

		Double_t empt  = daughterN.Pt();
		Double_t emeta = daughterN.Eta();
		Double_t emphi = daughterN.Phi();
		Double_t emsig = funSmearPt->Eval(empt);
		Double_t emsigEmb = funSmearPtEmb->Eval(empt);
		TVector3 eN(0,0,0);
		eN.SetPtEtaPhi(empt,emeta,emphi);
		Double_t smempt1  = eN.Pt()*(1. + momShape->GetRandom()*emsig/0.01);
		
 		// smear = GetSmear(eN.Pt()); 
		smear = GetSmear2(eN.Pt()); // new Smear for Function
		if(mDebug) cout << smear <<endl;
		//cout<<emsig/emsigEmb<<endl;
		Double_t smempt = eN.Pt()*(1.+smear*emsig/emsigEmb);
		// Double_t smempt = eN.Pt();// try to remove the pT smear
		delatpT = smempt-eN.Pt();
    // cout << delatpT <<endl;
		hMinusSmearvsPt->Fill(smempt,delatpT/smempt);
		if(mDebug) cout<<smempt1<<"/"<<smempt<<endl;
		//eN.SetPtEtaPhi(smempt,emeta,emphi);
		Double_t smemeta = eN.Eta();
		Double_t smemphi = eN.Phi();

		TLorentzVector smdaughterP(0,0,0,0);
		TLorentzVector smdaughterN(0,0,0,0);
		TLorentzVector smdaughterPwoSM(0,0,0,0);
		TLorentzVector smdaughterNwoSM(0,0,0,0);
		
		smdaughterP.SetPtEtaPhiM(smeppt,smepeta,smepphi,Masselectron);
		smdaughterN.SetPtEtaPhiM(smempt,smemeta,smemphi,Masselectron);
		smdaughterPwoSM.SetPtEtaPhiM(eppt,smepeta,smepphi,Masselectron);
		smdaughterNwoSM.SetPtEtaPhiM(empt,smemeta,smemphi,Masselectron);

		if(mDebug) cout << "before eff" << endl;
		int ietatpc,iphitpc,ipttpc; 
		int ietatof,iphitof,ipttof; 
		Double_t epeff3d, emeff3d;
		if(mUseScaleEff){
			// tpcPtEtaPhi2Bin(0, smdaughterP.Pt(),smdaughterP.Eta(),smdaughterP.Phi(), &ipttpc,&ietatpc,&iphitpc); // MB tpc binning
			// tofPtEtaPhi2Bin(0, smdaughterP.Pt(),smdaughterP.Eta(),smdaughterP.Phi(), &ipttof,&ietatof,&iphitof); // MB tof binning
			// epeff3d = EvalScaleEff3D(smdaughterP,ipttpc,ietatpc,iphitpc,ipttof,ietatof,iphitof,1);

			// tpcPtEtaPhi2Bin(0, smdaughterN.Pt(),smdaughterN.Eta(),smdaughterN.Phi(), &ipttpc,&ietatpc,&iphitpc); // MB tpc binning
			// tofPtEtaPhi2Bin(0, smdaughterN.Pt(),smdaughterN.Eta(),smdaughterN.Phi(), &ipttof,&ietatof,&iphitof); // MB tof binning
			// emeff3d = EvalScaleEff3D(smdaughterN,ipttpc,ietatpc,iphitpc,ipttof,ietatof,iphitof,-1);
		}
		else{
			//tpcPtEtaPhi2Bin(mCenIdx, smdaughterP.Pt(),smdaughterP.Eta(),smdaughterP.Phi(), &ipttpc,&ietatpc,&iphitpc);
			//tofPtEtaPhi2Bin(mCenIdx, smdaughterP.Pt(),smdaughterP.Eta(),smdaughterP.Phi(), &ipttof,&ietatof,&iphitof);
			getEtaPhiBin_TPC(smdaughterP.Eta(),smdaughterP.Phi(),&ietatpc,&iphitpc);
			getEtaPhiBin_TOF(smdaughterP.Eta(),smdaughterP.Phi(),&ietatof,&iphitof);
			if(mDebug) cout << "getEtaPhiBin posi"<<endl;
			if(mDebug) cout << "eta = "<<smdaughterP.Eta()<<endl;
			if(mDebug) cout<< "ieta = "<< ietatpc<<"\n iphi = "<<iphitpc<<endl;
			epeff3d = EvalEff3D(smdaughterP,ipttpc,ietatpc,iphitpc,ipttof,ietatof,iphitof,1);
			if(mDebug) cout<<"get Eff posi, eff = "<<epeff3d<<endl;

			//tpcPtEtaPhi2Bin(mCenIdx, smdaughterN.Pt(),smdaughterN.Eta(),smdaughterN.Phi(), &ipttpc,&ietatpc,&iphitpc);
			//tofPtEtaPhi2Bin(mCenIdx, smdaughterN.Pt(),smdaughterN.Eta(),smdaughterN.Phi(), &ipttof,&ietatof,&iphitof);
			//cout<<"get Eff posi"<<endl;
			getEtaPhiBin_TPC(smdaughterN.Eta(),smdaughterN.Phi(),&ietatpc,&iphitpc);
			getEtaPhiBin_TOF(smdaughterN.Eta(),smdaughterN.Phi(),&ietatof,&iphitof);
			if(mDebug) cout << "getEtaPhiBin elec"<<endl;
			if(mDebug) cout << "eta = "<<smdaughterN.Eta()<<endl;
			if(mDebug) cout<< "ieta = "<< ietatpc<<"\niphi = "<<iphitpc<<endl;
			emeff3d = EvalEff3D(smdaughterN,ipttpc,ietatpc,iphitpc,ipttof,ietatof,iphitof,-1);
			if(mDebug) cout<<"get Eff elec, eff = "<<emeff3d<<endl;
		}
		if(mDebug) cout<<"after Eff"<<endl;

		TLorentzVector eepair(0,0,0,0);
		eepair = daughterP + daughterN;

		TLorentzVector rceepair(0,0,0,0);
		rceepair = smdaughterP + smdaughterN;
		hEPSingleTrkEffvsPt->Fill(smdaughterP.Pt(),epeff3d);
		hEMSingleTrkEffvsPt->Fill(smdaughterN.Pt(),emeff3d);
		TLorentzVector rceepair2(0,0,0,0);
		rceepair2 = smdaughterNwoSM+smdaughterPwoSM;

		hRCPairRapidityvsParentRapidity->Fill(parent.Rapidity(), rceepair.Rapidity());

		//hMCPairPtvsParentPt->Fill(parent.Pt(),rceepair.Pt());
		hMCPairPtvsParentPt->Fill(parent.Pt(),eepair.Pt());

		// hMCMvsPt->Fill(rceepair.M(),rceepair.Pt());
		hMCMvsPt->Fill(rceepair.Pt(),rceepair.M());

		if(fabs(rceepair.Rapidity())<=1.){ 
			// hMCAcc0MvsPt->Fill(rceepair.M(),rceepair.Pt());
			hMCAcc0MvsPt->Fill(rceepair.Pt(),rceepair.M());
			hMCAcc0PairPtvsParentPt->Fill(parent.Pt(),rceepair.Pt());
			if(smeppt>=0.2 && smempt>=0.2
					&&fabs(smepeta)<=1. && fabs(smemeta)<=1.){
				hMCAcc1PairEMPtvsEPPt->Fill(smeppt,smempt);
				hMCAcc1PairPtvsParentPt->Fill(parent.Pt(),rceepair.Pt());
				hMCAcc1PairRapidity->Fill(rceepair.Rapidity());

				// hMCAcc1MvsPt->Fill(rceepair.M(),rceepair.Pt());
				// hRCAcc1MvsPt3D->Fill(rceepair.M(),rceepair.Pt(),epeff3d*emeff3d);
				hMCAcc1MvsPt->Fill(rceepair.Pt(),rceepair.M());
				hRCAcc1MvsPt3D->Fill(rceepair.Pt(),rceepair.M(),epeff3d*emeff3d);
        		/*if(mParIndex==6) 
        		{
        		  hMCAcc1PairEMPtvsEPPt->Fill(smeppt,smempt,1./1.19);
        		  hMCAcc1PairPtvsParentPt->Fill(parent.Pt(),rceepair.Pt(),1./1.19);
        		  hMCAcc1PairRapidity->Fill(rceepair.Rapidity(),1./1.19);

        		  hMCAcc1MvsPt->Fill(rceepair.Pt(),rceepair.M(),1./1.19);
        		  hRCAcc1MvsPt3D->Fill(rceepair.Pt(),rceepair.M(),epeff3d*emeff3d*1./1.19);
        		}*/
			}
			// else()
		}
		if(fabs(rceepair2.Rapidity())<=1.){ 
			if(eppt>=0.2 && empt>=0.2
					&&fabs(smepeta)<=1. && fabs(smemeta)<=1.){

				hMCAcc1MvsPtwoSmear->Fill(rceepair2.Pt(),rceepair2.M());
				// hMCAcc1MvsPtwoSmear->Fill(rceepair2.M(),rceepair2.Pt());
        /*if(mParIndex==6) 
        {
          hMCAcc1PairEMPtvsEPPt->Fill(smeppt,smempt,1./1.19);
          hMCAcc1PairPtvsParentPt->Fill(parent.Pt(),rceepair.Pt(),1./1.19);
          hMCAcc1PairRapidity->Fill(rceepair.Rapidity(),1./1.19);

          hMCAcc1MvsPt->Fill(rceepair.Pt(),rceepair.M(),1./1.19);
          hRCAcc1MvsPt3D->Fill(rceepair.Pt(),rceepair.M(),epeff3d*emeff3d*1./1.19);
        }*/
			}
		}
		if(mDebug) cout << "after sample eff" << endl;

    	RapiditywoSTARAcc->Fill(rceepair.Rapidity());
		if(rceepair.M()>mMass) 
		{
			CutRecorder->Fill(0.5);
			PtMVsPtPLargeMass->Fill(smempt,smeppt);
			EtaMVsEtaPLargeMass->Fill(smemeta,smepeta);
			if( smeppt<0.2 || smempt<0.2 || fabs(smepeta)>1. || fabs(smemeta)>1. ) 	CutRecorder->Fill(1.5);
			if( smeppt<0.2 ) CutRecorder->Fill(3.5);
			if( smempt<0.2 ) CutRecorder->Fill(4.5); 
			if( smepeta > 1. ) CutRecorder->Fill(5.5);
			if( smemeta > 1. ) CutRecorder->Fill(6.5); 
		}
		if(fabs(rceepair.Rapidity())>1)
    	      {
    	        MeePtRapidityOver1->Fill(rceepair.Pt(),rceepair.M());
    	      }
    	// cout << "Mee = " <<  rceepair.M() << endl;
    	// cout << "Rapidity = " << rceepair.Rapidity() <<endl;

    	if (smeppt>=0.2 && smempt>=0.2 && fabs(smepeta)<=1. && fabs(smemeta)<=1.)
    	{  
			MeeFullRapidity->Fill(rceepair.M());// w.o any rapidity cut
			if (rceepair.Rapidity()<=1)
			{
				RapiditySTARAcc->Fill(rceepair.Rapidity());
			}//STAR rapidity cut

		  if (fabs(rceepair.Rapidity()) <= 1) MeeWpTEtaWRapidity->Fill(rceepair.M());
		  else MeeWpTEtaWoRapidity->Fill(rceepair.M());
    	}
		else
		{
			RapidityOutSTARAcc->Fill(rceepair.Rapidity());
		}
	
		if (fabs(rceepair.Rapidity()) <= 1)
		{
			if(smeppt>=0.2 && smempt>=0.2) MeeWoEtaCut->Fill(rceepair.M());
			if(fabs(smepeta)<=1. && fabs(smemeta)<=1.) MeeWopTCut->Fill(rceepair.M());
		}
		// if(eppt>=0.2 && empt>=0.2) MeeWoEtaCut->Fill(rceepair.M());
		// if(fabs(smepeta)<=1. && fabs(smemeta)<=1.) MeeWopTCut->Fill(rceepair.M());
		if(fabs(rceepair.Rapidity()) <= 1) MeeWopTEtaCut->Fill(rceepair.M());

	}
	cout << "after All"<<endl;

	return;
}
