#include "headers.h"
#include "StMiniTreeMaker.h"
#include "/star/u/wangzhen/run20/Dielectron/Efficiecny/StRoot/StMiniTreeMaker/RefMfun.h"
#include "/star/u/wangzhen/run20/Dielectron/Efficiecny/pileup.h"
//For Embedding QA, this ValidTrack and isElectron need change cut selection
//!!!!!!!!!!!!!!!!!!!!!!!!!
//change log 
//19.1.25 when i select the prue electron sample to TOF match Eff,i have nsigmaEcut only,
//for EID cut 1/beta i have nsigmaE cut and when I get nsigmaE Eff i will try have 1/beta cut and no 1/beta cut
//19.1.29 I had changed the PE Selection as mass<0.005
//19.2.25 it is helpful to add a fuction to select the mode or make the code separately
//Do it!!!!!
//
/////////////////////////////////////////////////////////////////
ClassImp(StMiniTreeMaker)

	//_____________________________________________________________________________
	StMiniTreeMaker::StMiniTreeMaker(const Char_t *name) : StMaker(name), mFillHisto(1), mPrintConfig(1), mPrintMemory(0), mPrintCpu(0), mStreamName(""), fOutFile(0), mOutFileName(""), mEvtTree(0), mDefaultVtx(1), mSelectVtxRank(1), mMaxVtxR(1.e4), mMaxVtxZ(1.e4), mMaxVzDiff(1.e4), mMinTrkPt(0.2), mMaxTrkEta(1.), mMinNHitsFit(15), mMinNHitsFitRatio(0.52), mMinNHitsDedx(10), mMaxDca(10.), mMaxnSigmaE(2.5), mMinnSigmaE(-0.75), mMaxBeta2TOF(0.03), mElectronCutMode(0), mPhotonicEMassCut(0.015)
{
	// default constructor

	// run15 st_mtd
	mTriggerIDs.clear();
	mTriggerIDs.push_back(780010); // AuAu@54  minbias
	mTriggerIDs.push_back(780020); // AuAu@54  minbias
	// mTriggerIDs.push_back(580021); // AuAu@54  minbias
}
//_____________________________________________________________________________
StMiniTreeMaker::~StMiniTreeMaker()
{
	// default destructor
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Init()
{
	//	refMultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMBnoVtx()
	//	refMultCorr = new StRefMultCorr("grefmult_VpdMBnoVtx");

	cout << "StMiniTreeMaker::Init()" << endl;
	if (!mOutFileName.Length())
	{
		LOG_ERROR << "StMiniTreeMaker:: no output file specified for tree and histograms." << endm;
		return kStERR;
	}
	fOutFile = new TFile(mOutFileName.Data(), "recreate");
	LOG_INFO << "StMiniTreeMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;

	if (mFillHisto)
		bookHistos();
	LOG_INFO << "BookHistos Finish" << endm;
	PileupUplimit = new TF1("PileupUplimit","pol5",0,400);
	PileupUplimit->SetParameters(-15.7887025834219, 0.789786364309292, -0.000637115144252616, 1.00019972792727e-05, -2.45208851616324e-08);
	PileupLowlimit = new TF1("PileupLowlimit","pol5",0,400);
	PileupLowlimit->SetParameters(16.4277056306649, 1.71652229539398, -0.00406847684302521, 1.65203560938885e-05, -2.96250329214512e-08);
	PileupLimit = new TF1("PileupLimit","[0]*x-[1]",0,1000);
	PileupLimit->SetParameters(0.7,10);

	indata_001.open("/star/u/wangzhen/run20/Dielectron/BadRunList/BadRunList.dat");
	mBadRunId_001.clear();
	if(indata_001.is_open()){
		cout<<"read in bad run list for 9.2 GeV Au+Au GeV ";
		Int_t oldId;
		Int_t newId=0;
		while(indata_001>>oldId){
			mBadRunId_001[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total bad run  list !!!"<<endl;
		return kFALSE;
	}
	indata_001.close();

	return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::InitRun(const Int_t runnumber)
{
	LOG_INFO << "Grab runnumber: " << runnumber << endm;
	return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Finish()
{
	if (fOutFile)
	{
		fOutFile->cd();
		fOutFile->Write();
		fOutFile->Close();
		LOG_INFO << "StMiniTreeMaker::Finish() -> write out tree in " << mOutFileName.Data() << endm;
	}
	if (mPrintConfig)
		printConfig();
	return kStOK;
}
//_____________________________________________________________________________
Int_t StMiniTreeMaker::Make()
{

	StTimer timer;
	if (mPrintMemory)
		StMemoryInfo::instance()->snapshot();
	if (mPrintCpu)
		timer.start();

	mPicoDstMaker = (StPicoDstMaker *)GetMaker("picoDst");
	if (Debug())
	{
		LOG_INFO << "PicoDstMaker pointer: " << mPicoDstMaker << endm;
	}

	if (mPicoDstMaker)
	{
		if (Debug())
			LOG_INFO << "Use Pico file as input" << endm;
		mPicoDst = mPicoDstMaker->picoDst();
		if (!mPicoDst)
		{
			LOG_WARN << "No PicoDst !" << endm;
			return kStOK;
		}
	}
	else
	{
		LOG_WARN << "No StPicoDstMaker !" << endm;
		return kStOK;
	}

	if (!processPicoEvent())
		return kStOK;

	if (mPrintMemory)
	{
		StMemoryInfo::instance()->snapshot();
		StMemoryInfo::instance()->print();
	}

	if (mPrintCpu)
	{
		timer.stop();
		LOG_INFO << "CPU time for StMiniTreeMaker::Make(): "
				 << timer.elapsedTime() << "sec " << endm;
	}

	return kStOK;
}
//_____________________________________________________________________________
Bool_t StMiniTreeMaker::processPicoEvent()
{
	int debugflag = 1;//0 is not print debug information
	if (mFillHisto)
		hEvent->Fill(0.5);

	StPicoEvent *picoEvent = mPicoDst->event();
	if (!picoEvent)
	{
		LOG_WARN << "No event level information !" << endm;
		return kFALSE;
	}

	Bool_t validTrigger = kFALSE;
	Bool_t minbias = kFALSE;

	Int_t nTrigs = 0;
	for (Int_t i = 0; i < mTriggerIDs.size(); i++)
	{
		if (picoEvent->isTrigger(mTriggerIDs[i]))
		{
			minbias = kTRUE;
			validTrigger = kTRUE;
			nTrigs++;
		}
	}
	//Zhen add it for the bad run removal
	int runID =  picoEvent->runId();
	map<Int_t, Int_t>::iterator iter_001 = mBadRunId_001.find(runID);
	if(iter_001 != mBadRunId_001.end()){
		//cout<<"bad run, continue"<<endl;
		return kFALSE;
	}

	if (!validTrigger)
	{
		return kFALSE;
	}

	if (mFillHisto)
	{
		if (minbias)
			hEvent->Fill(2.5);
	}

	TVector3 vtxPos = picoEvent->primaryVertex();
	double vpdvz = picoEvent->vzVpd();
	Int_t refMult = picoEvent->refMult();
	Int_t mnTOFMatch = (int)picoEvent->nBTOFMatch();
	Double_t refMultCorr = refMult;
	double reweight = 1.;
	hRefMultvsnTOFMatch->Fill(refMult,mnTOFMatch);
	hRefMultvsnTOFMatchvsVz->Fill(refMult,mnTOFMatch,vtxPos.Z());

	// get the centrality
	Int_t mCentrality = GetCentrality(refMultCorr);// using 9.2 temp Centrality defination
	hRefMult->Fill(refMult);
	hCentrality9->Fill(mCentrality); // 9 for 0-5%, 1 f0r 70-80%


	mCentrality = mCentrality-1;
	if (mCentrality<0 || mCentrality>8) return kFALSE;
	if (debugflag == 1) cout<<"event centrality = "<<mCentrality<<endl;	
	Float_t mTpceNSigmaECutLow;
	Float_t mTpceNSigmaECutHi;


	if (mFillHisto)
	{
		hVtxYvsVtxX->Fill(vtxPos.x(), vtxPos.y());
		hVPDVzvsTPCVz->Fill(vtxPos.z(), vpdvz);
		hVzDiff->Fill(vtxPos.z() - vpdvz);
	}
	if (debugflag == 1) cout<<"after fill histo "<<endl;	

	if (TMath::Abs(vtxPos.x()) < 1.e-5 && TMath::Abs(vtxPos.y()) < 1.e-5 && TMath::Abs(vtxPos.z()) < 1.e-5) // non zero vertem
		return kFALSE;
	if (mFillHisto)
		hEvent->Fill(7.5);
	if (sqrt(vtxPos.x() * vtxPos.x() + vtxPos.y() * vtxPos.y()) >= mMaxVtxR) // vr cut
		return kFALSE;
	if (mFillHisto)
		hEvent->Fill(8.5);
	if (TMath::Abs(vtxPos.z()) >= mMaxVtxZ)// Vz Cut
		return kFALSE;
	if (mFillHisto)
		hEvent->Fill(9.5);
	// if (mnTOFMatch < PileupLimit->Eval(refMult)) return kFALSE;//pile up cut
	if (debugflag == 1) cout<<"before pile up "<<endl;	

	if (!pileupRejection(vtxPos.z(), refMult, mnTOFMatch)) return kFALSE;//using vz dependence vz cut
	if (mFillHisto)
		hEvent->Fill(10.5);

	Int_t nNodes = mPicoDst->numberOfTracks();
	if (Debug())
	{
		LOG_INFO << "# of global tracks in picoDst: " << nNodes << endm;
	}
	if (debugflag == 1 ) cout<< "# of global tracks in picoDst: " << nNodes <<endl;

	TLorentzVector pair;
	TLorentzVector *Electron = new TLorentzVector[MaxNElectron];
	TLorentzVector *Positron = new TLorentzVector[MaxNElectron];
	double *ElectronID = new double[MaxNElectron];
	double *PositronID = new double[MaxNElectron];
	double *ElectronTag = new double[MaxNElectron];
	double *PositronTag = new double[MaxNElectron];
	int nElectron = 0;
	int nPositron = 0;
	// reset the arrays
	for (int i = 0; i < MaxNElectron; i++)
	{
		Electron[i].SetPtEtaPhiE(0, 0, 0, 0);
		Positron[i].SetPtEtaPhiE(0, 0, 0, 0);
		ElectronID[i] = 0;
		PositronID[i] = 0;
		ElectronTag[i] = 0;
		PositronTag[i] = 0;
	}
	if (debugflag == 1 ) cout<< "before track loop " << nNodes <<endl;

	
	int nChargeParticle = 0;
	for (Int_t i = 0; i < nNodes; i++)
	{
		StPicoTrack *pTrack = mPicoDst->track(i);
		if (!pTrack)
			continue;

		if (!isValidTrack(pTrack, vtxPos))
			continue;

		if (debugflag == 1 ) cout<<"select the track information OK"<<endl;
		
		Int_t charge = pTrack->charge();
		TVector3 pMom = pTrack->pMom();
		Float_t p = pMom.Mag();
		Float_t pt = pMom.Perp();
		Float_t eta = pMom.PseudoRapidity();
		Float_t phi = pMom.Phi();
		Float_t dedx = pTrack->dEdx();
		Float_t nSigmaE = pTrack->nSigmaElectron();
		Float_t nSigmaP = pTrack->nSigmaProton();
		Float_t nSigmaPi = pTrack->nSigmaPion();
		Float_t nSigmaK = pTrack->nSigmaKaon();
		double dca =  pTrack->gDCA(vtxPos).Mag();
		Float_t beta = -999;
		Float_t tofLocalY = -999;
		Float_t TOFMatchFlag = -1;
		Int_t bTofPidTraitsIndex = pTrack->bTofPidTraitsIndex();
		double msquare = -0.5;
		if (bTofPidTraitsIndex >= 0)
		{
			StPicoBTofPidTraits *btofPidTraitis = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
			beta = btofPidTraitis->btofBeta();
      		double CellId = btofPidTraitis->btofCellId();
			hTofCellID->Fill(CellId);
      		// if(CellId == 8994 || CellId == 8998 || CellId == 8999) continue;//new for 54.4 GeV hot TOF Cell ID
			hBeta->Fill(beta);
			tofLocalY = btofPidTraitis->btofYLocal();
			hTofLocalY->Fill(tofLocalY);
			TOFMatchFlag = btofPidTraitis->btofMatchFlag();
			hTofMatchFlag->Fill(TOFMatchFlag);
			msquare = pow(p, 2) * (1 - pow(beta, 2)) / pow(beta, 2);
		}

		hMsquare->Fill(msquare);
		hMsquraevsRefMult->Fill(msquare,refMult);
		hMsquarevsP->Fill(p*charge,msquare);
		if (isPiKP_masscut(msquare)) nChargeParticle++;

		if (TMath::Abs(nSigmaP) < 4 && TMath::Abs(msquare-0.879)<0.020)
		{
			hPureProtonNSigmaEvsP->Fill(p,nSigmaE);
			hPureProtonNSigmaEvsPt->Fill(pt,nSigmaE);
			hPureProtonPhivsEta->Fill(phi,eta);
			hPureProtonNSigmaEvsPCen->Fill(p,nSigmaE,mCentrality,reweight);
		}
		if (TMath::Abs(nSigmaK) < 4 && TMath::Abs(msquare-0.243)<0.005)
		{
			hPureKaonNSigmaEvsP->Fill(p,nSigmaE);
			hPureKaonNSigmaEvsPt->Fill(pt,nSigmaE);
			hPureKaonNSigmaEvsPCen->Fill(p,nSigmaE,mCentrality,reweight);
		}
		if (TMath::Abs(nSigmaPi) < 4 && TMath::Abs(msquare-0.019)<0.003)
		{
			hPurePionSigmaPionvsP->Fill(p,nSigmaPi);
			hPurePionNSigmaEvsP->Fill(p,nSigmaE);
			hPurePionNSigmaEvsPt->Fill(pt,nSigmaE);
			hPurePionNSigmaEvsPCen->Fill(p,nSigmaE,mCentrality,reweight);
		}
		if (nSigmaPi > 6 && TMath::Abs(msquare-0.019)<0.003)
		{
			hMergePionNSigmaEvsP->Fill(p,nSigmaE);
			hMergePionNSigmaEvsPt->Fill(pt,nSigmaE);
			hMergePionNSigmaEvsPCen->Fill(p,nSigmaE,mCentrality,reweight);
		}

		if (debugflag == 1 ) cout<<"select the pure track information OK"<<endl;
		if (charge > 0)
		{
			/*if (TMath::Abs(nSigmaP) < 0.5)
			{
				hDenPPlusTofEff->Fill(pt, eta, phi);
				if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
				{
					hNumPPlusTofEff->Fill(pt, eta, phi);
				}
			}
			if (TMath::Abs(nSigmaK) < 0.5)
			{
				hDenKPlusTofEff->Fill(pt, eta, phi);
				if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
				{
					hNumKPlusTofEff->Fill(pt, eta, phi);
				}
			}*/
			if (TMath::Abs(nSigmaPi) < 0.5)
			{
				hDenPiPlusTofEff->Fill(pt, eta, phi);
				hDenPiPlusTofEffCen[mCentrality]->Fill(pt, eta, phi,reweight);
				if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
				{
					hNumPiPlusTofEff->Fill(pt, eta, phi);
					hNumPiPlusTofEffCen[mCentrality]->Fill(pt, eta, phi,reweight);
				}
			}
		}
		if (charge < 0)
		{
			/*if (TMath::Abs(nSigmaP) < 0.5)
			{
				hDenPMinusTofEff->Fill(pt, eta, phi);
				if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
				{
					hNumPMinusTofEff->Fill(pt, eta, phi);
				}
			}
			if (TMath::Abs(nSigmaK) < 0.5)
			{
				hDenKMinusTofEff->Fill(pt, eta, phi);
				if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
				{
					hNumKMinusTofEff->Fill(pt, eta, phi);
				}
			}*/
			if (TMath::Abs(nSigmaPi) < 0.5)
			{
				hDenPiMinusTofEff->Fill(pt, eta, phi);
				hDenPiMinusTofEffCen[mCentrality]->Fill(pt, eta, phi,reweight);
				if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
				{
					hNumPiMinusTofEff->Fill(pt, eta, phi);
					hNumPiMinusTofEffCen[mCentrality]->Fill(pt, eta, phi,reweight);
				}
			}
		}
		hSignaEvsP->Fill(p,nSigmaE); // nSigmaE vs p no cut
		hBetavsPwonSigmaE->Fill(p,1./beta);
		if (TMath::Abs(nSigmaE) <2.0) 
		if (isElectronBetaCut(pTrack)) // no beta cut
    		hBetavsP->Fill(p,1./beta);
		if (TMath::Abs(1.0-1./beta)<0.025)
			hSigmaEvsPwithNSigEandBeta->Fill(p,nSigmaE);
		switch (mElectronCutMode) 
		{
		  case 0:
		    if (!isElectron(pTrack)) { continue;} //used for pure electron selection
		    break;
		  case 1:
		    if (!isElectronSigmaECut(pTrack)) {continue;} // used for the nSigmaE cut Study
		    break;
		  case 2:
		    if (!isElectronBetaCut(pTrack)) {continue;} // used for the 1/beta cut study
		    break;
		  default:
		    LOG_WARN << "Wrong Electron cut mode !" << endm;
		    break;
		}
		\
        	/*if(mCentrality<=3) { //60-80%
                	mTpceNSigmaECutLow = -0.75;  // old: -0.75
                	mTpceNSigmaECutHi  = 2.0;  // old: 1.25
       		 }
        	else { //0-60%
                mTpceNSigmaECutLow = -1*0.6;
                mTpceNSigmaECutHi  = 0.6;
        	}
        	if(nSigmaE<mTpceNSigmaECutLow || nSigmaE>mTpceNSigmaECutHi) continue;
		*/
		
		/////////////////////////////////////////////////////////////////////////
		if (debugflag == 1 ) cout<<"select the electron track information OK"<<endl;
		if (charge < 0)
		{
			Electron[nElectron].SetPtEtaPhiM(pt, eta, phi, Melectron);
			ElectronID[nElectron] = i;
			nElectron++;
		}
		if (charge > 0)
		{
			Positron[nPositron].SetPtEtaPhiM(pt, eta, phi, Melectron);
			PositronID[nPositron] = i;
			nPositron++;
		}
	}

	for (int i = 0; i < nElectron; i++)
	{
		for (int j = 0; j < nPositron; j++)
		{
			pair = Electron[i] + Positron[j];
			hULMvsPt->Fill(pair.Pt(), pair.M());
			hULMvsPtLarge->Fill(pair.Pt(), pair.M());
			if (pair.M() < mPhotonicEMassCut)// select photonic electrons
			{
				ElectronTag[i] = 1;
				PositronTag[j] = 1;
			}
		}
		for (int k = i + 1; k < nElectron; k++)
		{
			pair = Electron[i] + Electron[k];
			hLNegMvsPt->Fill(pair.Pt(), pair.M());
			hLMNegvsPt->Fill(pair.Pt(), pair.M());
			if (pair.M() < mPhotonicEMassCut)
			{
				ElectronTag[i] = 0;
				ElectronTag[k] = 0;
			}
		}
	}
	for (int i = 0; i < nPositron; i++)
	{
		for (int j = i + 1; j < nPositron; j++)
		{
			pair = Positron[i] + Positron[j];
			hLPosMvsPt->Fill(pair.Pt(), pair.M());
			hLMPosvsPt->Fill(pair.Pt(), pair.M());
			if (pair.M() < mPhotonicEMassCut)
			{
				PositronTag[i] = 0;
				PositronTag[j] = 0;
			}
		}
	}
	if (debugflag  == 1 ) cout<<"number of elcetron = "<<nElectron<<endl<<"number of positron = "<<nPositron<<endl<<"ratio of Elec and Posi = "<<(double)nElectron/(double)nPositron<<endl;
	for (int i = 0; i < nElectron; i++)
	{
		if (ElectronTag[i] == 0)
			continue;
		if(debugflag == 1) cout<<"photonic electron information OK"<<endl;
		StPicoTrack *ElecTrack = mPicoDst->track(ElectronID[i]);
		TVector3 pMom = ElecTrack->pMom();
		Float_t p = pMom.Mag();
		Float_t pt = pMom.Perp();
		Float_t eta = pMom.PseudoRapidity();
		Float_t phi = pMom.Phi();
		Float_t dedx = ElecTrack->dEdx();
		Int_t nHitsFit = ElecTrack->nHitsFit();
		Int_t nHitsDedx = ElecTrack->nHitsDedx();
		Int_t nHitsPoss = ElecTrack->nHitsMax();
		Float_t nSigmaE = ElecTrack->nSigmaElectron();
		Float_t dca = ElecTrack->gDCA(vtxPos).Mag();

		double beta = -999.0;
		double tofLocalY = -999.0;
		int TOFMatchFlag = -1;

		Int_t bTofPidTraitsIndex = ElecTrack->bTofPidTraitsIndex();
		if (bTofPidTraitsIndex >= 0)
		{
			StPicoBTofPidTraits *btofPidTraitis = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
			beta = btofPidTraitis->btofBeta();
			hBeta->Fill(beta);
			tofLocalY = btofPidTraitis->btofYLocal();
			hTofLocalY->Fill(tofLocalY);
			TOFMatchFlag = btofPidTraitis->btofMatchFlag();
			hTofMatchFlag->Fill(TOFMatchFlag);
		}
		//hDenPEMinusTofEff->Fill(pt, eta, phi);
		hDenPEMinusTofEffCen[mCentrality]->Fill(pt, eta, phi);
		if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
		{
			//hNumPEMinusTofEff->Fill(pt, eta, phi);
			hNumPEMinusTofEffCen[mCentrality]->Fill(pt, eta, phi,reweight);
			hPureElectronNSigmaEvsP->Fill(p,nSigmaE);
			hPureElectronNSigmaEvsPCen->Fill(p,nSigmaE,mCentrality,reweight);
      		hPEElectronnSigmaEvsP->Fill(p,nSigmaE);
			hPureElectronNSigmaEvsPt->Fill(pt,nSigmaE);
			hPureElectronNSigmaEvsPtvsCen->Fill(pt,nSigmaE,mCentrality,reweight);
			hPureElectronNSigmaEvsPvsCen->Fill(p,nSigmaE,mCentrality,reweight);
			hPureElectronNSigmaEvsPhivsCen->Fill(p,nSigmaE,mCentrality,reweight);
		}
		// 1/beta cut eff
		hPEMinusBetavsP->Fill(p,1.0/beta);
		hPEMinusBetavsPt->Fill(pt,1.0/beta);
		hPEMinusBetavsPCen->Fill(p,1.0/beta,mCentrality,reweight);
		hPEMinusBetavsPtCen->Fill(pt,1.0/beta,mCentrality,reweight);
	}
	for (int i = 0; i < nPositron; i++)
	{
		if (PositronTag[i] == 0)
			continue;
		if (debugflag == 1) cout<<"photonic positron information OK"<<endl;
		StPicoTrack *PosiTrack = mPicoDst->track(PositronID[i]);
		TVector3 pMom = PosiTrack->pMom();
		Float_t p = pMom.Mag();
		Float_t pt = pMom.Perp();
		Float_t eta = pMom.PseudoRapidity();
		Float_t phi = pMom.Phi();
		Float_t dedx = PosiTrack->dEdx();
		Int_t nHitsFit = PosiTrack->nHitsFit();
		Int_t nHitsDedx = PosiTrack->nHitsDedx();
		Int_t nHitsPoss = PosiTrack->nHitsMax();
		Float_t nSigmaE = PosiTrack->nSigmaElectron();
		Float_t dca = PosiTrack->gDCA(vtxPos).Mag();

		double beta = -999.0;
		double tofLocalY = -999.0;
		int TOFMatchFlag = -1;

		Int_t bTofPidTraitsIndex = PosiTrack->bTofPidTraitsIndex();
		if (bTofPidTraitsIndex >= 0)
		{
			StPicoBTofPidTraits *btofPidTraitis = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
			beta = btofPidTraitis->btofBeta();
			hBeta->Fill(beta);
			tofLocalY = btofPidTraitis->btofYLocal();
			hTofLocalY->Fill(tofLocalY);
			TOFMatchFlag = btofPidTraitis->btofMatchFlag();
			hTofMatchFlag->Fill(TOFMatchFlag);
		}
		//hDenPEPlusTofEff->Fill(pt, eta, phi);
		hDenPEPlusTofEffCen[mCentrality]->Fill(pt, eta, phi,reweight);
		//nSigmaE counting method to get nSigmaEeff
		//if (nSigmaE)
		//hDenNSigmaEPEPlusPtEtaPhi->Fill(pt, eta, phi);

		// TOF match eff
		if ((beta > 0.) && (TOFMatchFlag > 0) && (TMath::Abs(tofLocalY) < 1.8))
		{
			//hNumPEPlusTofEff->Fill(pt, eta, phi);
			hNumPEPlusTofEffCen[mCentrality]->Fill(pt, eta, phi,reweight);
			hPureElectronNSigmaEvsP->Fill(p,nSigmaE);
			hPureElectronNSigmaEvsPCen->Fill(p,nSigmaE,mCentrality,reweight);
			hPEPositronnSigmaEvsP->Fill(p,nSigmaE);
			hPurePositronNSigmaEvsPt->Fill(pt,nSigmaE);
			hPurePositronNSigmaEvsPtvsCen->Fill(pt,nSigmaE,mCentrality,reweight);
			hPurePositronNSigmaEvsPvsCen->Fill(p,nSigmaE,mCentrality,reweight);
			hPurePositronNSigmaEvsPhivsCen->Fill(p,nSigmaE,mCentrality,reweight);
		}
		// 1/beta cut eff
		hPEPlusBetavsP->Fill(p,1.0/beta);
		hPEPlusBetavsPCen->Fill(p,1.0/beta,mCentrality,reweight);
		hPEPlusBetavsPt->Fill(pt,1.0/beta);
		hPEPlusBetavsPtCen->Fill(pt,1.0/beta,mCentrality,reweight);

		/*hPEPlusBetavsP->Fill(p,beta);
		hPEPlusDcavsPt->Fill(pt, dca);
		hPEPlusDcavsEta->Fill(eta, dca);
		hPEPlusPhivsPt->Fill(pt, phi);
		hPEPlusEtavsPt->Fill(pt, eta);
		//hPEPlusYvsPt->Fill(pt,y);
		hPEPlusNHitsFitvsPt->Fill(pt, nHitsFit);
		hPEPlusNHitsFitvsPhi->Fill(phi, nHitsFit);
		hPEPlusNHitsFitvsEta->Fill(eta, nHitsFit);
		hPEPlusNHitsPossvsPt->Fill(pt, nHitsPoss);
		hPEPlusNHitsPossvsPhi->Fill(phi, nHitsPoss);
		hPEPlusNHitsPossvsEta->Fill(eta, nHitsPoss);
		hPEPlusNHitsDedxvsPt->Fill(pt, nHitsDedx);
		hPEPlusDedxvsP->Fill(p, dedx);
		hPEPlusDedxvsPhi->Fill(phi, dedx);
		hPEPlusDedxvsEta->Fill(eta, dedx);
		hPEPlusNSigmaEvsP->Fill(p, nSigmaE);
		hPEPlusNSigmaEvsPhi->Fill(phi, nSigmaE);
		hPEPlusNSigmaEvsEta->Fill(eta, nSigmaE);*/

	}
	hRefMultvsnChargeParticle->Fill(refMult,nChargeParticle);
	hnTOFMatchvsnChargePartile->Fill(mnTOFMatch,nChargeParticle);

	delete[] Electron;
	delete[] Positron;
	delete[] ElectronID;
	delete[] PositronID;
	delete[] ElectronTag;
	delete[] PositronTag;
	return kTRUE;
}
//_____________________________________________________________________________
//select the Pi/K/P/E by the mass square cut
bool StMiniTreeMaker::isPiKP_masscut(double msquare)
{
	if(TMath::Abs(msquare-0.879)<0.020 || TMath::Abs(msquare-0.243)<0.005 || TMath::Abs(msquare-0.019)<0.003) return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMiniTreeMaker::isValidTrack(StPicoTrack *pTrack, TVector3 vtxPos) const
{
	Float_t pt = pTrack->pMom().Perp();
	Float_t eta = pTrack->pMom().PseudoRapidity();
	//Float_t dca = pTrack->dca();
	//Float_t dca = (pTrack->dca()-vtxPos).mag();
	Float_t dca = pTrack->gDCA(vtxPos).Mag();

	if (pt < 0.2 || pt > 30.0)
		return kFALSE;
	if (TMath::Abs(eta) > 1.0)
		return kFALSE;
	if (pTrack->nHitsFit() < 20)
		return kFALSE;
	if (pTrack->nHitsFit() * 1. / pTrack->nHitsMax() < 0.52)
		return kFALSE;
	if (pTrack->nHitsDedx() < 15)
		return kFALSE;
	if (dca > 1)
		return kFALSE;

	return kTRUE;
}
//_____________________________________________________________________________
void StMiniTreeMaker::bookHistos()
{
	hEvent = new TH1D("hEvent", "Event statistics", 25, 0, 25);
	hEvent->GetXaxis()->SetBinLabel(1, "All events");
	hEvent->GetXaxis()->SetBinLabel(3, "minbias");
	hEvent->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
	hEvent->GetXaxis()->SetBinLabel(9, Form("|V_{r}|<%1.2f cm", mMaxVtxR));
	hEvent->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm", mMaxVtxZ));
	hEvent->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm", mMaxVzDiff));

	hVtxYvsVtxX = new TH2D("hVtxYvsVtxX", "hVtxYvsVtxX; V_{x} (cm); V_{y} (cm)", 120, -3, 3, 120, -3, 3);
	hVPDVzvsTPCVz = new TH2D("hVPDVzvsTPCVz", "hVPDVzvsTPCVz; TPC V_{z} (cm); VPD V_{z} (cm)", 800, -200, 200, 800, -200, 200);
	hVzDiff = new TH1D("hVzDiff", "hVzDiff; TPC V_{z} - VPD V_{z} (cm)", 400, -20, 20);
	hRefMult = new TH1D("hRefMult","hRefMult;RefMult;Counts",800,0,800);
	hCentrality9 = new TH1D("hCentrality9","hCentrality9;Centrality",11,0,11);

	hdEdxvsP = new TH2D("hdEdxvsP", "hdEdxvsP; p (GeV/c); dE/dx (KeV/cm)", 500, 0, 5, 400, 0, 20);
	hnSigEvsP = new TH2D("hnSigEvsP", "hnSigEvsP; p (GeV/c); n#sigma_{e}", 3000, 0, 3, 1200, -6, 6);
	hnSigEvsPWTOF = new TH2D("hnSigEvsPWTOF", "hnSigEvsP; p (GeV/c); n#sigma_{e}", 3000, 0, 3, 1200, -6, 6);
	hBetavsP = new TH2D("hBetavsP", "hBetavsP; p (GeV/c); 1/#beta", 3000, 0, 3, 800, 0.5, 1.3);
	hBetavsPwonSigmaE = new TH2D("hBetavsPwonSigmaE", "hBetavsPwonSigmaE; p (GeV/c); 1/#beta", 3000, 0, 3, 800, 0.5, 1.3);
	hULMvsphiV = new TH2D("hULMvsphiV", "UnLike-sign mass vs phiV;M_{primary} (GeV/c^{2});#phi_{V}", 2000, 0, 0.2, 180, 0, TMath::Pi());
	hLPosMvsphiV = new TH2D("hLPosMvsphiV", "Like-sign mass vs phiV;M_{primary} (GeV/c^{2});#phi_{V}", 2000, 0, 0.2, 180, 0, TMath::Pi());
	hLNegMvsphiV = new TH2D("hLNegMvsphiV", "Like-sign mass vs phiV;M_{primary} (GeV/c^{2});#phi_{V}", 2000, 0, 0.2, 180, 0, TMath::Pi());

	//histo for cal TOF match Eff
	hDenPiPlusTofEff = new TH3F("hDenPiPlusTofEff", "hDenPiPlusTofEff;p_{T} (GeV/c); #eta; #phi", 500, 0, 5, 100, -1, 1, 120, -PI - PI / 12, PI + PI / 12);
	hNumPiPlusTofEff = new TH3F("hNumPiPlusTofEff", "hNumPiPlusTofEff;p_{T} (GeV/c); #eta; #phi", 500, 0, 5, 100, -1, 1, 120, -PI - PI / 12, PI + PI / 12);
	hDenPiMinusTofEff = new TH3F("hDenPiMinusTofEff", "hDenPiMinusTofEff;p_{T} (GeV/c); #eta; #phi", 500, 0, 5, 100, -1, 1, 120, -PI - PI / 12, PI + PI / 12);
	hNumPiMinusTofEff = new TH3F("hNumPiMinusTofEff", "hNumPiMinusTofEff;p_{T} (GeV/c); #eta; #phi", 500, 0, 5, 100, -1, 1, 120, -PI - PI / 12, PI + PI / 12);
	hULMvsPt = new TH2F("hULMvsPt", "hULMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})", 500, 0, 5, 500, 0, 0.2);
	hLPosMvsPt = new TH2F("hLPosMvsPt", "hLPosMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})", 500, 0, 5, 500, 0, 0.2);
	hLNegMvsPt = new TH2F("hLNegMvsPt", "hLNegMvsPt;p_{T} (GeV/c); M_{ee} (GeV/c^{2})", 500, 0, 5, 500, 0, 0.2);
	hPEPlusBetavsP = new TH2F("hPEPlusBetavsP", "hPEPlusBetavsP;p (GeV/c); 1/#beta", 500, 0, 5, 200, 0.9, 1.1);
	hPEMinusBetavsP = new TH2F("hPEMinusBetavsP", "hPEMinusBetavsP;p (GeV/c); 1/#beta", 500, 0, 5, 200, 0.9, 1.1);
	hPEPlusBetavsPt = new TH2F("hPEPlusBetavsPt", "hPEPlusBetavsPt;p_{T} (GeV/c); 1/#beta", 500, 0, 5, 200, 0.9, 1.1);
	hPEMinusBetavsPt = new TH2F("hPEMinusBetavsPt", "hPEMinusBetavsPt;p_{T} (GeV/c); 1/#beta", 500, 0, 5, 200, 0.9, 1.1);
	
	hBeta = new TH1D("hBeta", "hbeta;#beta;counts", 200, 0, 2);
	hTofLocalY = new TH1D("hTofLocalY", "hTofLocalY;localY (cm);counts", 270, -0.2, 2.5);
	hTofMatchFlag = new TH1D("TofMatchFlag", "TofMatchFlag;TofMatchFlag;counts", 13, -3, 10);
	hTofCellID = new TH1D("hTofCellID","hTofCellID;CellID;nEvts",70000,0,70000);

	// define histogram for check untag and likesign method
	htrkMassdistrbution = new TH1D("htrkMassdistrbution", "htrkMassdistrbution;mass", 200, -0.4, 1.6);
	hULMvsPtLarge = new TH2D("hULMvsPtLarge", "hULMvsPtLarge;p_{T} (GeV/c);M_{ee}", 500, 0, 5, 500, 0, 5);
	hLMPosvsPt = new TH2D("hLMPosvsPt", "hLMPosvsPt;p_{T} (GeV/c);M_{ee}", 500, 0, 5, 500, 0, 5);
	hLMNegvsPt = new TH2D("hLMNegvsPt", "hLMNegvsPt;p_{T} (GeV/c);M_{ee}", 500, 0, 5, 500, 0, 5);
	hULDcavsPt = new TH2D("hULDcavsPt", "hULDcavsPt;p_{T} (GeV/c);dca (cm)", 500, 0, 5, 300, 0, 3);
	hLPosDcavsPt = new TH2D("hLPDcavsPt", "hLPDcavsPt;p_{T} (GeV/c);dca (cm)", 500, 0, 5, 300, 0, 3);
	hLNegDcavsPt = new TH2D("hLNDcavsPt", "hLNDcavsPt;p_{T} (GeV/c);dca (cm)", 500, 0, 5, 300, 0, 3);
	hRefMultvsnTOFMatch = new TH2D("hRefMultvsnTOFMatch","hRefMultvsnTOFMatch; RefMult; nTOFMatch",500,0,500,500,0,500);
	hMsquraevsRefMult = new TH2D("hMsquraevsRefMult","hMsquraevsRefMult; m^{2} (GeV/c)^2; RefMult",1000,-0.1,1.2,500,0,500);
	hRefMultvsnChargeParticle = new TH2D("hRefMultvsnChargeParticle","hRefMultvsnChargeParticle; RefMult; nPi/K/P",500,0,500,500,0,500);
	hnTOFMatchvsnChargePartile = new TH2D("hnTOFMatchvsnChargePartile","hnTOFMatchvsnChargePartile; nTOFMatch;nChargePartile",500,0,500,500,0,500);
	hRefMultvsnTOFMatchvsVz = new TH3D("hRefMultvsnTOFMatchvsVz","hRefMultvsnTOFMatchvsVz;RefMult;nTOFMatch;Vz (cm)",500,0,500,500,0,500,5,-145,145);

	//histo for Different Centrality efficiency
	for(int i  = 0;i<9;i++)
	{
		int nPtBins = 400;
		int nEtaBins = 120;
		int nPhiBins = 120;
		hDenPEPlusTofEffCen[i] = new TH3F(Form("hDenPEPlusTofEffCen%d",i), Form("hDenPEPlusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		hNumPEPlusTofEffCen[i] = new TH3F(Form("hNumPEPlusTofEffCen%d",i), Form("hNumPEPlusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		hDenPEMinusTofEffCen[i] = new TH3F(Form("hDenPEMinusTofEffCen%d",i), Form("hDenPEMinusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		hNumPEMinusTofEffCen[i] = new TH3F(Form("hNumPEMinusTofEffCen%d",i), Form("hNumPEMinusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		hDenPiPlusTofEffCen[i] = new TH3F(Form("hDenPiPlusTofEffCen%d",i), Form("hDenPiPlusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		hNumPiPlusTofEffCen[i] = new TH3F(Form("hNumPiPlusTofEffCen%d",i), Form("hNumPiPlusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		hDenPiMinusTofEffCen[i] = new TH3F(Form("hDenPiMinusTofEffCen%d",i), Form("hDenPiMinusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		hNumPiMinusTofEffCen[i] = new TH3F(Form("hNumPiMinusTofEffCen%d",i), Form("hNumPiMinusTofEffCen%d;p_{T} (GeV/c); #eta; #phi",i), nPtBins, 0, 5, nEtaBins, -1.2, 1.2, nPhiBins, -PI - PI / 12, PI + PI / 12);
		cout << i << endl;
  }
	hPureElectronNSigmaEvsPtvsCen = new TH3D(Form("hPureElectronNSigmaEvsPtvsCen"),Form("hPureElectronNSigmaEvsPtvsCen;p_{T} (GeV/C);n#sigma_{e}^{e}"),800,0,4,400,-19.005,20.995,16,0,16);	
	hPurePositronNSigmaEvsPtvsCen = new TH3D(Form("hPurePositronNSigmaEvsPtvsCen"),Form("hPurePositronNSigmaEvsPtvsCen;p_{T} (GeV/C);n#sigma_{e}^{e}"),800,0,4,400,-19.005,20.995,16,0,16);
	hPureElectronNSigmaEvsPvsCen = new TH3D(Form("hPureElectronNSigmaEvsPvsCen"),Form("hPureElectronNSigmaEvsPvsCen;p (GeV/C);n#sigma_{e}^{e}"),800,0,4,400,-19.005,20.995,16,0,16);	
	hPurePositronNSigmaEvsPvsCen = new TH3D(Form("hPurePositronNSigmaEvsPvsCen"),Form("hPurePositronNSigmaEvsPvsCen;p (GeV/C);n#sigma_{e}^{e}"),800,0,4,400,-19.005,20.995,16,0,16);
	hPureElectronNSigmaEvsPhivsCen = new TH3D(Form("hPureElectronNSigmaEvsPhivsCen"),Form("hPureElectronNSigmaEvsPhivsCen;phi ;n#sigma_{e}^{e}"),600,-PI - PI / 12, PI + PI / 12,400,-19.005,20.995,16,0,16);	
	hPurePositronNSigmaEvsPhivsCen = new TH3D(Form("hPurePositronNSigmaEvsPhivsCen"),Form("hPurePositronNSigmaEvsPhivsCen;phi ;n#sigma_{e}^{e}"),600,-PI - PI / 12, PI + PI / 12,400,-19.005,20.995,16,0,16);
	hPEPlusBetavsPCen = new TH3F(Form("hPEPlusBetavsPCen"), Form("hPEPlusBetavsPCen;p (GeV/c); 1/#beta"), 500, 0, 5, 200, 0.9, 1.1,16,0,16);
	hPEMinusBetavsPCen = new TH3F(Form("hPEMinusBetavsPCen"), Form("hPEMinusBetavsPCen;p (GeV/c); 1/#beta"), 500, 0, 5, 200, 0.9, 1.1,16,0,16);
	hPEPlusBetavsPtCen = new TH3F(Form("hPEPlusBetavsPtCen"), Form("hPEPlusBetavsPtCen;p (GeV/c); 1/#beta"), 500, 0, 5, 200, 0.9, 1.1,16,0,16);
	hPEMinusBetavsPtCen = new TH3F(Form("hPEMinusBetavsPtCen"), Form("hPEMinusBetavsPtCen;p (GeV/c); 1/#beta"), 500, 0, 5, 200, 0.9, 1.1,16,0,16);

	// for check something that in the start of analysis
	hMsquare = new TH1D("hMsquare","hMsquare;mass^(2) (GeV/c^{2})^{2};counts",510,-0.5,1.2);
	hMsquarevsP = new TH2D("hMsquarevsP","hMsquarevsP; p*q (GeV/C); mass^(2) (GeV/c^{2})^{2}",3000,-3,3,1800,-0.4,1.4);
	hPurePionNSigmaEvsP = new TH2D("hPurePionNSigmaEvsP","hPurePionNSigmaEvsP;p (GeV/C);n#sigma_{e}^{#pi}",500,0,4,400,-19.005,20.995);
	hMergePionNSigmaEvsP = new TH2D("hMergePionNSigmaEvsP","hMergePionNSigmaEvsP;p (GeV/C);n#sigma_{e}^{#pi}",500,0,4,400,-20.005,19.995);
	hPureKaonNSigmaEvsP = new TH2D("hPureKaonNSigmaEvsP","hPureKaonNSigmaEvsP;p (GeV/C);n#sigma_{e}^{K}",500,0,4,400,-19.005,20.995);
	hPureProtonNSigmaEvsP = new TH2D("hPureProtonNSigmaEvsP","hPureProtonNSigmaEvsP;p (GeV/C);n#sigma_{e}^{P}",500,0,4,400,-19.005,20.995);
	hPurePionNSigmaEvsPCen = new TH3D("hPurePionNSigmaEvsPCen",";p (GeV/C);n#sigma_{e}^{#pi};Centrality",500,0,4,400,-19.005,20.995,16,0,16);
	hMergePionNSigmaEvsPCen = new TH3D("hMergePionNSigmaEvsPCen",";p (GeV/C);n#sigma_{e}^{#Merge#pi};Centrality",500,0,4,400,-20.005,19.995,16,0,16);
	hPureKaonNSigmaEvsPCen = new TH3D("hPureKaonNSigmaEvsPCen",";p (GeV/C);n#sigma_{e}^{K};Centrality",500,0,4,400,-19.005,20.995,16,0,16);
	hPureProtonNSigmaEvsPCen = new TH3D("hPureProtonNSigmaEvsPCen",";p (GeV/C);n#sigma_{e}^{P};Centrality",500,0,4,400,-19.005,20.995,16,0,16);
	hPureElectronNSigmaEvsPCen = new TH3D("hPureElectronNSigmaEvsPCen",";p (GeV/C);n#sigma_{e}^{P};Centrality",500,0,4,400,-19.005,20.995,16,0,16);
	hPEElectronnSigmaEvsP = new TH2D("hPEElectronnSigmaEvsP","hPEElectronnSigmaEvsP;p (GeV/C);n#sigma_{e}^{P}",500,0,4,400,-19.005,20.995);
	hPEPositronnSigmaEvsP = new TH2D("hPEPositronnSigmaEvsP","hPEPositronnSigmaEvsP;p (GeV/C);n#sigma_{e}^{P}",500,0,4,400,-19.005,20.995);
	//Pt
	hPurePionNSigmaEvsPt = new TH2D("hPurePionNSigmaEvsPt","hPurePionNSigmaEvsPt;p_{T} (GeV/C);n#sigma_{e}^{#pi}",800,0,4,400,-19.005,20.995);
  hMergePionNSigmaEvsPt = new TH2D("hMergePionNSigmaEvsPt","hMergePionNSigmaEvsPt;p_{T} (GeV/C);n#sigma_{e}^{#pi}",800,0,4,400,-20.995,19);
  hPureKaonNSigmaEvsPt = new TH2D("hPureKaonNSigmaEvsPt","hPureKaonNSigmaEvsPt;p_{T} (GeV/C);n#sigma_{e}^{K})",800,0,4,400,-19.005,20.995);
  hPureProtonNSigmaEvsPt = new TH2D("hPureProtonNSigmaEvsPt","hPureProtonNSigmaEvsPt;p_{T} (GeV/C);n#sigma_{e}^{P})",800,0,4,400,-20.995,29);
  hPurePionSigmaPionvsP = new TH2D("hPurePionSigmaPionvsP","p (GeV/c);n#sigma_{#pi}^{pi}",800,0,4,400,-2,2);

	hPureElectronNSigmaEvsPt = new TH2D("hPureElectronNSigmaEvsPt","hPureElectronNSigmaEvsPt;p_{T} (GeV/C);n#sigma_{e}^{e})",800,0,4,400,-19.005,20.995);
	hPurePositronNSigmaEvsPt = new TH2D("hPurePositronNSigmaEvsPt","hPurePositronNSigmaEvsPt;p_{T} (GeV/C);n#sigma_{e}^{e})",800,0,4,400,-19.005,20.995);
	//
	hPureProtonPhivsEta = new TH2D("hPureProtonPhivsEta","hPureProtonPhivsEta;#phi;#eta",640, -3.2, 3.2,260,-1.3,1.3);
	hPureElectronNSigmaEvsP = new TH2D("hPureElectronNSigmaEvsP","hPureElectronNSigmaEvsP;p (GeV/C);n#sigma_{e}^{e})",500,0,4,400,-19.005,20.995);
	hSignaEvsP = new TH2D("hSigmaEvsP","hSignaEvsP;p (GeV/C);n#sigma_{e}",500,0,5,400,-19.005,20.995);
	hSigmaEvsPwithNSigmaE = new TH2D("hSigmaEvsPwithNSigmaE","hSigmaEvsPwithNSigmaE;p (GeV/C);n#sigma_{e}",500,0,5,400,-19.005,20.995);
	hSigmaEvsPwithNSigEandBeta = new TH2D("hSigmaEvsPwithBeta","hSigmaEvsPwithBeta;p (GeV/C);n#sigma_{e}",500,0,5,400,-19.005,20.995);
}
//_____________________________________________________________________________
void StMiniTreeMaker::printConfig()
{
	const char *decision[2] = {"no", "yes"};
	printf("=== Configuration for StMiniTreeMaker ===\n");
	printf("Fill the QA histo: %s\n", decision[mFillHisto]);
	printf("Use default vertex: %s\n", decision[mDefaultVtx]);
	printf("Select positive vertex ranking: %s\n", decision[mSelectVtxRank]);
	printf("Maximum |Vr|: %1.2f\n", mMaxVtxR);
	printf("Maximum |Vz|: %1.2f\n", mMaxVtxZ);
	printf("Maximum |VzDiff|: %1.2f\n", mMaxVzDiff);
	printf("Minimum track pt: %1.2f\n", mMinTrkPt);
	printf("Maximum track |eta| : %1.2f\n", mMaxTrkEta);
	printf("Minimum number of fit hits: %d\n", mMinNHitsFit);
	printf("Minimum ratio of fit hits: %1.2f\n", mMinNHitsFitRatio);
	printf("Minimum number of dedx hits: %d\n", mMinNHitsDedx);
	printf("Maximum dca: %1.2f\n", mMaxDca);
	printf("Maximum min nSigmaE for TPCe: %1.2f\n", mMinnSigmaE);
	printf("Maximum max nSigmaE  for TPCe: %1.2f\n", mMaxnSigmaE);
	printf("Maximum |1-1/beta| for TPCe: %1.3f\n", mMaxBeta2TOF);
	printf("=======================================\n");
}
double StMiniTreeMaker::phiVangle(TLorentzVector e1, TLorentzVector e2, int q1, int q2)
{
	Double_t pt1 = e1.Pt();
	Double_t eta1 = e1.Eta();
	Double_t phi1 = e1.Phi();

	Double_t pt2 = e2.Pt();
	Double_t eta2 = e2.Eta();
	Double_t phi2 = e2.Phi();

	TRandom3 *myRandom = new TRandom3();
	;
	TVector3 e1Mom, e2Mom;
	if (q1 > 0 && q2 < 0)
	{
		e2Mom.SetPtEtaPhi(pt1, eta1, phi1); //e+
		e1Mom.SetPtEtaPhi(pt2, eta2, phi2); //e-
	}
	else if (q1 < 0 && q2 > 0)
	{
		e2Mom.SetPtEtaPhi(pt2, eta2, phi2); //e+
		e1Mom.SetPtEtaPhi(pt1, eta1, phi1); //e-
	}
	else if (q1 == q2 && TMath::Abs(q1) == 1)
	{
		Double_t ran = myRandom->Uniform(-1, 1);
		if (ran > 0)
		{
			e2Mom.SetPtEtaPhi(pt1, eta1, phi1);
			e1Mom.SetPtEtaPhi(pt2, eta2, phi2);
		}
		else
		{
			e2Mom.SetPtEtaPhi(pt2, eta2, phi2);
			e1Mom.SetPtEtaPhi(pt1, eta1, phi1);
		}
	}
	else
		return -1;

	TVector3 pu = e1Mom + e2Mom;
	TVector3 pv = e1Mom.Cross(e2Mom);
	TVector3 pw = pu.Cross(pv);
	TVector3 pnz(0., 0., -1);
	TVector3 pwc = pu.Cross(pnz);

	Double_t angleV = pw.Angle(pwc);

	return angleV;
}
bool StMiniTreeMaker::isElectron(StPicoTrack *pTrack) // full cut for electrons
{

	Int_t bTofPidTraitsIndex = pTrack->bTofPidTraitsIndex();
	double beta = -999;
	if (bTofPidTraitsIndex >= 0)
	{
		StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
		beta = btofPidTraits->btofBeta();
	}
	if(beta<=0 || TMath::Abs(1.- 1./beta) > mMaxBeta2TOF) return kFALSE;
	TVector3 pMom = pTrack->pMom();
	double p = pMom.Mag();
	Float_t mTpceNSigmaECutLow;
	if (p < 0.8)
	{
		mTpceNSigmaECutLow = 3.0*p - 3.15;
	}
	else
	{
		mTpceNSigmaECutLow = -0.75;
	}
	double nsigmaE = pTrack->nSigmaElectron();
	// if (nsigmaE < (mTpceNSigmaECutLow) || nsigmaE >= 2.0) return kFALSE;
	if (TMath::Abs(nsigmaE)>2.0) return kFALSE;//loose cut to pruity
	return kTRUE;
}

bool StMiniTreeMaker::isElectronSigmaECut(StPicoTrack *pTrack) // wider nSigmaE cut
{

	Int_t bTofPidTraitsIndex = pTrack->bTofPidTraitsIndex();
	double beta = -999;
	if (bTofPidTraitsIndex >= 0)
	{
		StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
		beta = btofPidTraits->btofBeta();
		double CellId = btofPidTraits->btofCellId();
		if(CellId == 8994 || CellId == 8998 || CellId == 8999) return kFALSE;
	}
	if(beta<=0 || TMath::Abs(1.- 1./beta) > mMaxBeta2TOF) return kFALSE;
	double nsigmaE = pTrack->nSigmaElectron();
	if(TMath::Abs(nsigmaE)>2.) return kFALSE;
	return kTRUE;
}

bool StMiniTreeMaker::isElectronBetaCut(StPicoTrack *pTrack) // no beta cut
{
	Int_t bTofPidTraitsIndex = pTrack->bTofPidTraitsIndex();
	TVector3 pMom = pTrack->pMom();
	double p = pMom.Mag();
	Float_t mTpceNSigmaECutLow;
	if (p < 0.8)
	{
		mTpceNSigmaECutLow = 3.0*p - 3.15;
	}
	else
	{
		mTpceNSigmaECutLow = -0.75;
	}
	double nsigmaE = pTrack->nSigmaElectron();
	if (nsigmaE < (mTpceNSigmaECutLow) || nsigmaE >= 1.66) return kFALSE;
	return kTRUE;
}
