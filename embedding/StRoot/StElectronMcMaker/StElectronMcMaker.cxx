#include "StElectronMcMaker.h"

#include "TROOT.h"
#include "TSystem.h"
#include "StChain.h"
#include "TF1.h"
#include "TRandom.h"
#include "TChain.h"
#include "TTreeHelper.h"
#include "TDatime.h"
#include "StarRoot/TUnixTime.h"
#include "StMessMgr.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"
#include "StThreeVector.hh"
#include "StThreeVectorF.hh"
#include "StThreeVectorD.hh"
#include "StLorentzVectorD.hh"

#include "StTriggerIdCollection.h"
#include "StTriggerId.h"

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "StEvent/StBTofHeader.h"
//StEmc
#include "StEvent/StEmcCollection.h"
#include "StEmcCluster.h"
#include "StEmcDetector.h"
#include "StEmcModule.h"
#include "StEmcPoint.h"
#include "StEmcRawHit.h"
#include "StEvent/StEmcCollection.h"
#include "StEmcClusterCollection.h"
#include "StEvent/StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/StEmcRawMaker.h"
//#include "StEmcTriggerMaker/StEmcTriggerMaker.h"
//#include "StEmcTriggerMaker/StBemcTrigger.h"
#include "StTriggerUtilities/StTriggerSimuMaker.h"
#include "StTriggerUtilities/StTriggerSimuResult.h"
#include "StTriggerUtilities/Bemc/StBemcTriggerSimu.h"

#include "StEvent.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTriggerIdCollection.h"
#include "StEvent/StTriggerId.h"
#include "StEvent/StPrimaryVertex.h"

#include "StMcEvent/StMcTpcHit.hh"

#include "StEvent.h"
#include "StEvent/StTrack.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHit.h"
#include "StEvent/StBTofPidTraits.h"

#include "StMcEvent/StMcEvent.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcTrack.hh"
#include "StEventUtilities/StuRefMult.hh"

#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"
#include "StAssociationMaker/StMcParameterDB.h"

// #include "StMuDSTMaker/COMMON/StMuDstMaker.h"
// #include "StMuDSTMaker/COMMON/StMuEvent.h"
// #include "StMuDSTMaker/COMMON/StMuDst.h"
// #include "StMuDSTMaker/COMMON/StMuTrack.h"

#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDebug.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include "StBTofHeader.h"
#include "StPhysicalHelixD.hh"

#include <math.h>
#include "TRandom.h"
#include "Random.h"
#include "RanluxEngine.h"
#include "RandGauss.h"

static RanluxEngine engine;
static RandGauss ranGauss(engine);

StRefMultCorr * refmultCorrUtil;

//_________________________________________________
StGlobalTrack*  partnerTrack(mcTrackMapType* map, StMcTrack* mT, int *nCommon) {
	mcTrackMapIter i = map->find(mT);
	StGlobalTrack* rT = 0;
	*nCommon = 0;
	if (i != map->end()) {
		if( (int) ((*i).second->commonTpcHits()) > (*nCommon) ) {
			rT = const_cast<StGlobalTrack*>((*i).second->partnerTrack());
			*nCommon = (*i).second->commonTpcHits();
		}
	}
	return rT; 
} 
//_________________________________________________
StMcTrack*  partnerMcTrack(rcTrackMapType* map, StGlobalTrack* rT, int *nCommon) {
	rcTrackMapIter i = map->find(rT);
	StMcTrack* mT = 0;
	*nCommon = 0;
	if (i != map->end()) {
		if( (int) ((*i).second->commonTpcHits()) > (*nCommon) ) {
			mT = const_cast<StMcTrack*>((*i).second->partnerMcTrack());
			*nCommon = (*i).second->commonTpcHits();
		}
	}
	return mT;
} 

ClassImp(StElectronMcMaker)

	//_________________________________________________
StElectronMcMaker::StElectronMcMaker(const char *name, const char *file):StMaker(name)
{
	//book tree and auto-save to TFile
	mOutputFile = new TFile(file,"RECREATE");
	mOutputFile->SetCompressionLevel(1);
	int BufSize=(int)pow(2.,16.);
	mTree = new TTree("mcT","mcT",BufSize);
	mTree->SetAutoSave(1000000); // autosave every 1 Mbytes
	mTree->Branch("runId",&mElectron.runId,"runId/I");
	mTree->Branch("triggerId",&mElectron.triggerId,"triggerId[4]/I");//MB for isobar 2018
	mTree->Branch("muEvtId",&mElectron.muEvtId,"muEvtId/I");
	mTree->Branch("muPriVertexX",&mElectron.muPriVertexX,"muPriVertexX/F");
	mTree->Branch("muPriVertexY",&mElectron.muPriVertexY,"muPriVertexY/F");
	mTree->Branch("muPriVertexZ",&mElectron.muPriVertexZ,"muPriVertexZ/F");
	//mTree->Branch("rcEvtId",&mElectron.rcEvtId,"rcEvtId/I");
	//mTree->Branch("rcVertexRanking",&mElectron.rcVertexRanking,"rcVertexRanking/F");
	//mTree->Branch("rcVertexX",&mElectron.rcVertexX,"rcVertexX/F");
	//mTree->Branch("rcVertexY",&mElectron.rcVertexY,"rcVertexY/F");
	//mTree->Branch("rcVertexZ",&mElectron.rcVertexZ,"rcVertexZ/F");
	mTree->Branch("muVpdVz",&mElectron.muVpdVz,"muVpdVz/F");
	mTree->Branch("muRefMult",&mElectron.muRefMult,"muRefMult/I");
	mTree->Branch("mugRefMult",&mElectron.mugRefMult,"mugRefMult/I");
	mTree->Branch("munBtofMatch",&mElectron.munBtofMatch,"munBtofMatch/I");
	//mTree->Branch("rcRefMultCorr",&mElectron.rcRefMultCorr,"rcRefMultCorr/F");
	//mTree->Branch("gWeight",&mElectron.gWeight,"gWeigth/F");
	mTree->Branch("zdcX",&mElectron.zdcX,"zdcX/F");
	mTree->Branch("Centrality_16",&mElectron.Centrality_16,"Centrality_16/I");
        mTree->Branch("Centrality_9",&mElectron.Centrality_9,"Centrality_9/I");
        mTree->Branch("refmult_corr",&mElectron.refmult_corr,"refmult_corr/F");
        mTree->Branch("reweight",&mElectron.reweight,"reweight/F");

	mTree->Branch("nMcTrks",&mElectron.nMcTrks,"nMcTrks/I");
	mTree->Branch("nMcVertices",&mElectron.nMcVertices,"nMcVertices/I");
	mTree->Branch("nMcE",&mElectron.nMcE,"nMcE/I");
	mTree->Branch("nRcE",&mElectron.nRcE,"nRcE/I");
	mTree->Branch("nMcTrack",&mElectron.nMcTrack,"nMcTrack/I");
	mTree->Branch("geantId",mElectron.geantId,"geantId[nMcE]/I");
	mTree->Branch("mcTrackId",mElectron.mcTrackId,"mcTrackId[nMcE]/I");
	mTree->Branch("mcCharge",mElectron.mcCharge,"mcCharge[nMcE]/S");
	//mTree->Branch("mcTrkKey",mElectron.mcTrkKey,"mcTrkKey[nMcE]/I");
	//mTree->Branch("parentGeantId",mElectron.parentGeantId,"parentGeantId[nMcE]/I");
	//mTree->Branch("mcParentTrkKey",mElectron.mcParentTrkKey,"mcParentTrkKey[nMcE]/I");
	mTree->Branch("mcTrack_Vr",mElectron.mcTrack_Vr,"mcTrack_Vr[nMcE]/F");
	mTree->Branch("mcTrack_par_Geantid",mElectron.mcTrack_par_Geantid,"mcTrack_par_Geantid[nMcE]/I");
	mTree->Branch("mcPt",mElectron.mcPt,"mcPt[nMcE]/F");
	mTree->Branch("mcEta",mElectron.mcEta,"mcEta[nMcE]/F");
	mTree->Branch("mcPhi",mElectron.mcPhi,"mcPhi[nMcE]/F");
	//mTree->Branch("mcPtFirst",mElectron.mcPtFirst,"mcPtFirst[nMcE]/F");
	//mTree->Branch("mcEtaFirst",mElectron.mcEtaFirst,"mcEtaFirst[nMcE]/F");
	//mTree->Branch("mcPhiFirst",mElectron.mcPhiFirst,"mcPhiFirst[nMcE]/F");
	//mTree->Branch("mcPtLast",mElectron.mcPtLast,"mcPtLast[nMcE]/F");
	//mTree->Branch("mcEtaLast",mElectron.mcEtaLast,"mcEtaLast[nMcE]/F");
	//mTree->Branch("mcPhiLast",mElectron.mcPhiLast,"mcPhiLast[nMcE]/F");
	//mTree->Branch("rcTrkKey",mElectron.rcTrkKey,"rcTrkKey[nMcE]/I");
	mTree->Branch("rcidtruth",mElectron.rcidtruth,"rcidtruth[nRcE]/I");
	mTree->Branch("rcQaTruth",mElectron.rcQaTruth,"rcQaTruth[nRcE]/I");
	mTree->Branch("rcflag",mElectron.rcflag,"rcflag[nRcE]/I");
	mTree->Branch("rcCharge",mElectron.rcCharge,"rcCharge[nRcE]/S");
	mTree->Branch("rcPt",mElectron.rcPt,"rcPt[nRcE]/F");
	mTree->Branch("rcEta",mElectron.rcEta,"rcEta[nRcE]/F");
	mTree->Branch("rcPhi",mElectron.rcPhi,"rcPhi[nRcE]/F");
	// mTree->Branch("rcPtFirst",mElectron.rcPtFirst,"rcPtFirst[nMcE]/F");
	// mTree->Branch("rcEtaFirst",mElectron.rcEtaFirst,"rcEtaFirst[nMcE]/F");
	// mTree->Branch("rcPhiFirst",mElectron.rcPhiFirst,"rcPhiFirst[nMcE]/F");
	// mTree->Branch("rcPtLast",mElectron.rcPtLast,"rcPtLast[nMcE]/F");
	// mTree->Branch("rcEtaLast",mElectron.rcEtaLast,"rcEtaLast[nMcE]/F");
	// mTree->Branch("rcPhiLast",mElectron.rcPhiLast,"rcPhiLast[nMcE]/F");
	mTree->Branch("rcNHitsFit",mElectron.rcNHitsFit,"rcNHitsFit[nRcE]/I");
	mTree->Branch("rcNHitsPoss",mElectron.rcNHitsPoss,"rcNHitsPoss[nRcE]/I");
	mTree->Branch("rcNHitsDedx",mElectron.rcNHitsDedx,"rcNHitsDedx[nRcE]/I");
	//mTree->Branch("rcNHitsCommon",mElectron.rcNHitsCommon,"rcNHitsCommon[nMcE]/I");
	mTree->Branch("rcDedx",mElectron.rcDedx,"rcDedx[nRcE]/F");
	mTree->Branch("rcNSigmaE",mElectron.rcNSigmaE,"rcNSigmaE[nRcE]/F");
	mTree->Branch("rcNSigmaPi",mElectron.rcNSigmaPi,"rcNSigmaPi[nRcE]/F");
	mTree->Branch("rcNSigmaK",mElectron.rcNSigmaK,"rcNSigmaK[nRcE]/F");
	mTree->Branch("rcNSigmaP",mElectron.rcNSigmaP,"rcNSigmaP[nRcE]/F");
	mTree->Branch("rcDca",mElectron.rcDca,"rcDca[nRcE]/F");
	// mTree->Branch("rcDca1",mElectron.rcDca1,"rcDca1[nMcE]/F");
	// mTree->Branch("rcDca2",mElectron.rcDca2,"rcDca2[nMcE]/F");
	// mTree->Branch("rcDca3",mElectron.rcDca3,"rcDca3[nMcE]/F");
	// mTree->Branch("rcDca4",mElectron.rcDca4,"rcDca4[nMcE]/F");
	// mTree->Branch("rcDca5",mElectron.rcDca5,"rcDca5[nMcE]/F");
	// mTree->Branch("rcDca6",mElectron.rcDca6,"rcDca6[nMcE]/F");


	// - zero all pointers defined in the header file
	mMuDst=0;
	//mMuEvent=0;
	//mMcVertex.setX(0.); mMcVertex.setY(0.); mMcVertex.setZ(0.);
	//mAssocMaker=0;
	//mRcTrackMap=0;
	//mMcTrackMap=0;
}
//_________________________________________________
StElectronMcMaker::~StElectronMcMaker()
{
}
//_________________________________________________
void StElectronMcMaker::Clear(const char*)
{
	StMaker::Clear();
}
//_________________________________________________
Int_t StElectronMcMaker::Finish()
{
	mOutputFile->Write();
	mOutputFile->Close();
	return StMaker::Finish();
}
//_________________________________________________
Int_t StElectronMcMaker::Init()
{
    cout<<"ustc Init"<<endl;
	Clear("C");

	//Zhen change it ,no centrality definition
	// for 7.7 using the temporary centrality defination
	//intialize the StRefMultCorr
	//mRefMultCorr = new StRefMultCorr();

	return kStOK;
}
//_________________________________________________
Int_t StElectronMcMaker::Make()
{
   LOG_INFO << "StElectronMcMaker Make!!!!!!" << endm;
    memset(&mElectron,0,sizeof(mElectron)); 
     StMuDstMaker * mMuDstMaker = (StMuDstMaker* )GetMaker("MuDst");
    if(!mMuDstMaker) {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    StMuEvent * mMuEvent = (StMuEvent*)mMuDst->event();
    if( !mMuEvent ){
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }
	
	//vertex selection
     int const originalVertexId = mMuDst->currentVertexIndex();
     StMuPrimaryVertex* selectedVertex = nullptr;
     // choose the default vertex, i.e. the first vertex
     mMuDst->setVertexIndex(0);
     selectedVertex = mMuDst->primaryVertex();
     // fall back to default vertex if no vertex is selected in the algorithm above.
     // should skip this event in the event cuts below.
     if ( ! selectedVertex ){
	  LOG_INFO << "Vertex is not valid" << endm;
	  //cout<<originalVertexId<<endl;
	  mMuDst->setVertexIndex(originalVertexId);
     }
     //end of vertex selection
	 
	 
	 //vertex is not selected
     if ( ! selectedVertex ) return kStOk;
     //trigger
     if ( ! mMuEvent->triggerIdCollection().nominal().isTrigger(780010) && ! mMuEvent->triggerIdCollection().nominal().isTrigger(780020) ) return kStOK ; //7.7 GeV trigger
    //  if ( ! mMuEvent->triggerIdCollection().nominal().isTrigger(810010) && ! mMuEvent->triggerIdCollection().nominal().isTrigger(810020) && ! mMuEvent->triggerIdCollection().nominal().isTrigger(810030) && ! mMuEvent->triggerIdCollection().nominal().isTrigger(810040) ) return kStOK ; //7.7 GeV trigger
     //Vz
     if ( fabs(mMuEvent->primaryVertexPosition().z()) > 70.0 ) return kStOk ;
     //Vr
     if ( mMuEvent->primaryVertexPosition().perp() > 2.0 ) return kStOk ;
	 
	 
	 for(int i = 0 ;i<4;i++){
		 mElectron.triggerId[i] = -999;
	 }
	 
	 
	 if(mMuEvent->triggerIdCollection().nominal().isTrigger(780010)) mElectron.triggerId[0] =780010;
	 if(mMuEvent->triggerIdCollection().nominal().isTrigger(780020)) mElectron.triggerId[1] =780020;
	//  if(mMuEvent->triggerIdCollection().nominal().isTrigger(810030)) mElectron.triggerId[2] =810030;
	//  if(mMuEvent->triggerIdCollection().nominal().isTrigger(810040)) mElectron.triggerId[3] =810040;
	 
	 
	 mElectron.runId = -999;
     mElectron.muEvtId = -999;
	 
	 mElectron.runId = mMuEvent->runId();
     mElectron.muEvtId = mMuEvent->eventId();

	 mElectron.muPriVertexX = mMuEvent->primaryVertexPosition().x();
	 mElectron.muPriVertexY = mMuEvent->primaryVertexPosition().y();
	 mElectron.muPriVertexZ = mMuEvent->primaryVertexPosition().z();
	 
	 mElectron.muVpdVz = -999;
         if(StBTofHeader * header = mMuDst->btofHeader()){
            mElectron.muVpdVz = header->vpdVz();
            cout<<"vpdVz: "<<header->vpdVz()<<endl;
         }
//	 mElectron.muVpdVz = mMuEvent->vpdVz();
	 
	 mElectron.muRefMult = -999;
	 mElectron.mugRefMult = -999;
	 mElectron.zdcX = -999;
	 mElectron.muRefMult = mMuEvent->refMult();
	 mElectron.mugRefMult = mMuEvent->grefmult();
	 mElectron.zdcX = mMuEvent->runInfo().zdcCoincidenceRate();//hz
     mElectron.munBtofMatch = mMuDst->primaryVertex()->nBTOFMatch();

   

    //    StTriggerSimuMaker* mTriggerSimuMaker = (StTriggerSimuMaker*) GetMaker("StarTrigSimu");
    //    assert(mTriggerSimuMaker);
    //    mBemcTriggerSimu = (StBemcTriggerSimu*) mTriggerSimuMaker->bemc;

       for(int i=1;i<=4800;i++)
       {
          //if(mBemcTriggerSimu->getHT6bitAdc(i)>11) cout<<i<<" "<<mBemcTriggerSimu->getHT6bitAdc(i)<<endl;
       }


        refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
        refmultCorrUtil->init(mMuEvent->runId());
        Bool_t isPileUpEvt_Cen = !refmultCorrUtil->passnTofMatchRefmultCut(mMuEvent->refMult()*1.0,mMuDst->primaryVertex()->nBTOFMatch()*1.0);
        if(isPileUpEvt_Cen) return kStOk;

        refmultCorrUtil->initEvent(mMuEvent->refMult()*1.0,mMuEvent->primaryVertexPosition().z(),mMuEvent->runInfo().zdcCoincidenceRate());
        mElectron.Centrality_16 = refmultCorrUtil->getCentralityBin16();
        mElectron.Centrality_9 = refmultCorrUtil->getCentralityBin9();
        mElectron.refmult_corr = refmultCorrUtil->getRefMultCorr();
        mElectron.reweight = refmultCorrUtil->getWeight();

	 
	 //fill MC histograms
     //The MC arrays in MuDST
     TClonesArray *MuMcVertices   = mMuDst->mcArray(0);
     Int_t NoMuMcVertices = MuMcVertices->GetEntriesFast();
     TClonesArray *MuMcTracks     = mMuDst->mcArray(1); 
     Int_t NoMuMcTracks = MuMcTracks->GetEntriesFast();
     LOG_INFO <<"# of MC tracks = "<< NoMuMcTracks <<" # of MC vertices = "<< NoMuMcVertices << endm;
     if (! NoMuMcVertices || ! NoMuMcTracks) {
	  LOG_WARN << "Ev.  has no MC information ==> skip it" << endm;
	  return kStWarn;
     }
     
	 
	 Int_t nMc = 0; //# of Mc track
	 Int_t nMcE = 0;//# of Mc electron
	//looking for electon
	mElectron.nMcTrks = NoMuMcTracks;
	mElectron.nMcVertices = NoMuMcVertices;
	
    // Loop for MC tracks
     for(Int_t itrk=0; itrk<NoMuMcTracks; itrk++){
	  StMuMcTrack *mcTrack = (StMuMcTrack *) MuMcTracks->UncheckedAt(itrk);
	  if (! mcTrack) continue;

	  // Select only Triggered Mc Vertex, i.e. the MC track should originate from PV (IdVx=1)
	  Int_t IdVx = mcTrack->IdVx();
	  if (IdVx != 1) continue;
	   StMuMcVertex *mcVertex = (StMuMcVertex *) MuMcVertices->UncheckedAt(IdVx-1);//test
          //if(pow((mcVertex->XyzV().x()*mcVertex->XyzV().x()+mcVertex->XyzV().y()*mcVertex->XyzV().y()),0.5)<1) continue;//test
      if (! mcTrack->Charge()) continue;
   
	  const int Gid = mcTrack->GePid();

	  nMc++;  // # of MC tracks
	  
	  
	  mElectron.geantId[nMcE] = Gid;
	  mElectron.mcCharge[nMcE] = mcTrack->Charge();
	  if(Gid==2 || Gid==3){//e+ or e-
	
                mElectron.mcTrackId[nMcE] = mcTrack->Id();    
		StThreeVectorF mcMom = mcTrack->Pxyz();
 		
                StThreeVectorF mcTrack_vertex_pos = mcVertex->XyzV();
                Float_t mcTrack_Vx = mcTrack_vertex_pos.x();               
                Float_t mcTrack_Vy = mcTrack_vertex_pos.y(); 
                mElectron.mcTrack_Vr[nMcE] = sqrt(pow(mcTrack_Vx,2)+pow(mcTrack_Vy,2)); 

                Int_t idMcTrack = mcVertex->IdParTrk();
                if(!idMcTrack){
                    mElectron.mcTrack_par_Geantid[nMcE] = -999;
                }else{
                    StMuMcTrack * mcTrack_par = (StMuMcTrack*)MuMcTracks->UncheckedAt(idMcTrack-1);
                    const int mcTrack_par_Gid = mcTrack_par->GePid();
                    mElectron.mcTrack_par_Geantid[nMcE] = mcTrack_par_Gid;
                }
                             
       
		mElectron.mcPt[nMcE] = mcMom.perp();
		mElectron.mcEta[nMcE] = mcMom.pseudoRapidity();
		mElectron.mcPhi[nMcE] = mcMom.phi();
		
		
		nMcE++;
	  }
	  else {
	     LOG_WARN << "Gid: "<<Gid<<" in Ev. "<<endm;
	  }
     }
	 
	 mElectron.nMcE = nMcE;
	 mElectron.nMcTrack = nMc;
  
	 
	 
	 Int_t nRcE = 0;//# of Rc electron
	 
	 //fill Event QA histograms
     TObjArray* tracks = mMuDstMaker->muDst()->primaryTracks() ;
     TObjArrayIter GetTracks(tracks) ;
     StMuTrack* ptrack ;
	 StMuTrack* ptrack_primary ;
     while ( ( ptrack = (StMuTrack*)GetTracks.Next() ) )
     {
		 
	  
	  if (ptrack->idTruth() <= 0 || ptrack->idTruth() > NoMuMcTracks) {
	     //cout << "Illegal idTruth " << ptrack->idTruth() << " The track is ignored" << endl;
	     continue;
	  }
	  StMuMcTrack *mcTrack = (StMuMcTrack *) MuMcTracks->UncheckedAt(ptrack->idTruth()-1);
	  if (!mcTrack) {
	     LOG_WARN << "Inconsistency in mcArray(1), ignored" << endm;
	     continue;
	  }
	  if (mcTrack->Id() != ptrack->idTruth()) {
	     LOG_WARN << "Mismatched idTruth " << ptrack->idTruth() << " and mcTrack Id " <<  mcTrack->Id()
		  << " this track is ignored" <<  endm;
	  }
	  Int_t idMcVx = mcTrack->IdVx();
	  //while (idMcVx != 1) {
	     //StMuMcVertex *mcVertex = (StMuMcVertex *) MuMcVertices->UncheckedAt(idMcVx-1);
	     //Int_t idMcTrack = mcVertex->IdParTrk();
	     //if (! idMcTrack) break;
	     //StMuMcTrack *mcTrackP = (StMuMcTrack *) MuMcTracks->UncheckedAt(idMcTrack-1);
	     //idMcVx = mcTrackP->IdVx();
	     //if (! idMcVx) break;
	 // }
	  if (idMcVx != 1) continue; //this MC track is not eventually originated from PV

	  if(mcTrack->GePid() != 2 && mcTrack->GePid() != 3) continue;//
	  if(mcTrack->IdVx() != 1) {
	     LOG_WARN<<"mc track may not directly originate from PV!"<<endm;
	  }
	  if (! ptrack->charge())  continue;
	  
	  //ptrack_primary = (StMcTrack*)ptrack->primaryTrack();
	 // if(!ptrack_primary){
		///  LOG_WARN << "this is not primary track" << endm;
	   //  continue;
	//  }
	  
	  mElectron.rcidtruth[nRcE] = ptrack->idTruth();
	  mElectron.rcQaTruth[nRcE] = ptrack->qaTruth();
	  mElectron.rcflag[nRcE] = ptrack->flag();
	  mElectron.rcCharge[nRcE] = ptrack->charge();
	  
	  mElectron.rcPt[nRcE] = ptrack->momentum().perp();
	  mElectron.rcEta[nRcE] = ptrack->momentum().pseudoRapidity();
	  mElectron.rcPhi[nRcE] = ptrack->momentum().phi();
	  
	  mElectron.rcNHitsFit[nRcE] = ptrack->nHitsFit(kTpcId);
	  mElectron.rcNHitsPoss[nRcE] = ptrack->nHitsPoss(kTpcId);
	  mElectron.rcNHitsDedx[nRcE] = ptrack->nHitsDedx();
	  
	  mElectron.rcDedx[nRcE] = ptrack->dEdx();//?
	  mElectron.rcNSigmaE[nRcE] = ptrack->nSigmaElectron();
	  mElectron.rcNSigmaPi[nRcE] = ptrack->nSigmaPion();
	  mElectron.rcNSigmaK[nRcE] = ptrack->nSigmaKaon();
	  mElectron.rcNSigmaP[nRcE] = ptrack->nSigmaProton();
	  
	  mElectron.rcDca[nRcE] = (Float_t)ptrack->dcaGlobal().mag();//?
	  
	  nRcE++;
	  
	  //end of the filling
     }
	 
     mElectron.nRcE = nRcE;


	mTree->Fill();

	return kStOK;
}
//_________________________________________________
StTrackPairInfo* StElectronMcMaker::findBestMatchedGlobal(StMcTrack *mcTrack)
{
	pair<mcTrackMapIter,mcTrackMapIter> mcBounds
		= mMcTrackMap->equal_range(mcTrack);
	StTrackPairInfo* candTrackPair = 0;  //used for finding the best matched track
	StGlobalTrack* candTrack = 0;
	for(mcTrackMapIter mcMapIter = mcBounds.first ; mcMapIter != mcBounds.second; ++mcMapIter){
		StTrackPairInfo* assocPair = (*mcMapIter).second;
		StGlobalTrack* globTrack =(StGlobalTrack*)assocPair->partnerTrack();
		if(!globTrack || globTrack->flag()<=0) continue;
		if(globTrack->fitTraits().numberOfFitPoints(kTpcId)>=10){
			if(!candTrackPair){
				candTrackPair = assocPair;
				candTrack = globTrack;
			}
			else if(globTrack->fitTraits().numberOfFitPoints(kTpcId) > candTrack->fitTraits().numberOfFitPoints(kTpcId)){
				candTrackPair = assocPair;
				candTrack = globTrack;
			}
		} //nHitsFit requirement
	} //bounds loop
	return candTrackPair; //Note that candTrack might be zero, for example if only one track is matched and has 9 tpc fit pts.
}


void StElectronMcMaker::matchEmc(StMuTrack* track, StMuDst *mudst, Float_t* emcInfo)
{
  //cout<<"star matchEmc"<<endl;
  for(int i=0;i<11;i++) emcInfo[i] = -999.;
  if(!track||track->flag()<=0) return;
  //if(!isGoodTrack(track)) return;
  StEmcCollection* pEmcCol = mudst->emcCollection();
  //if(track->geometry()->momentum().perp()<1.5) return;
  //cout<<" pass momentum cut"<<endl;
  Int_t mod, eta, sub;
  StThreeVectorD position, momentum;
  StThreeVectorD positionBSMDE, momentumBSMDE;
  StThreeVectorD positionBSMDP, momentumBSMDP;
  Double_t bFld;
  bFld = mudst->event()->magneticField()/10.; // bFld in Tesla
  bool ok = false;
  bool okBSMDE = false;
  bool okBSMDP = false;
 
 
  //cout<<"matchEmc1"<<endl;
  StEmcPosition* pEmcPosition = new StEmcPosition();
  StEmcGeom* pGeo[4];
  for(int i=0;i<4;i++) {
     if(i!=1) pGeo[i] = StEmcGeom::getEmcGeom(detname[i].Data());
  }
  if(pEmcPosition){
    ok      = pEmcPosition->projTrack(&position, &momentum, track, bFld, pGeo[0]->Radius());
    okBSMDE = pEmcPosition->projTrack(&positionBSMDE, &momentumBSMDE, track, bFld, pGeo[2]->Radius());
    okBSMDP = pEmcPosition->projTrack(&positionBSMDP, &momentumBSMDP, track, bFld, pGeo[3]->Radius());  
  }

 // cout<<"matchEmc2"<<endl;
  if(ok && okBSMDE && okBSMDP) {
    //cout<<" projection ok"<<endl;
    StSPtrVecEmcPoint& bEmcPoints = pEmcCol->barrelPoints();
    float EmcE = 0.;
    float maxtowerE = -999.;
    float phiDist = -999.;
    float zDist = -999.;
    int nEta = 0;
    int nPhi = 0;
    float minDist = 999.;
    unsigned int maxadc = 0;
    unsigned int maxdsmadc = 0;
    unsigned int softId = -1;
    //cout<<"matchEmc skf"<<endl;

    if(!pGeo[0]) pGeo[0]=StEmcGeom::getEmcGeom(detname[0].Data());
    pGeo[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub);
    //cout<<" (pt,p,eta) "<<track->geometry()->momentum().perp()<<"  "<<track->geometry()->momentum().mag()<<"  "<<track->geometry()->momentum().pseudoRapidity()<<endl;
    for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it!= bEmcPoints.end(); it++) {
     // cout<<"matchEmc skf1"<<endl;
      bool associated = false;
      StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
      StPtrVecEmcCluster& smdeClusters = (*it)->cluster(kBarrelSmdEtaStripId);
      StPtrVecEmcCluster& smdpClusters = (*it)->cluster(kBarrelSmdPhiStripId);
      //if(smdeClusters.size()==0 || smdpClusters.size()==0 || bEmcClusters.size()==0) continue;
      //if(smdeClusters[0]==NULL || smdpClusters[0]==NULL || bEmcClusters[0]==NULL) continue;
      if( bEmcClusters.size()==0 || bEmcClusters[0]==NULL) continue;
      for(StPtrVecEmcClusterIterator cIter = bEmcClusters.begin(); cIter != bEmcClusters.end(); cIter++) {
       // cout<<"matchEmc skf2"<<endl;
	StPtrVecEmcRawHit& bEmcHits = (*cIter)->hit();
	for(StPtrVecEmcRawHitIterator hIter = bEmcHits.begin(); hIter != bEmcHits.end(); hIter++) {
         // cout<<"matchEmc skf3"<<endl;
	  if(mod == (Int_t) (*hIter)->module() && eta == (Int_t) (*hIter)->eta() && sub == (Int_t) (*hIter)->sub()) {
	    //cout<<"(mod,eta,sub,adc,energy)=("<<mod<<", "<<eta<<", "<<sub<<", "<<(*hIter)->adc()<<", "<<(*hIter)->energy()<<")"<<endl;
	    associated = true;
	    break;
	  }
	}//end of raw hits loop
        //cout<<"matchEmc skf"<<endl;
	if(associated) {
	  for(StPtrVecEmcRawHitIterator hitit = bEmcHits.begin(); hitit != bEmcHits.end(); hitit++) {
           // cout<<"matchEmc ustc1"<<endl;
	    if((*hitit)->energy()>maxtowerE) maxtowerE = (*hitit)->energy();
	    if((*hitit)->adc()>maxadc) maxadc = (*hitit)->adc();
        //     cout<<"maxtowerE: "<<maxtowerE<<endl;
          //   cout<<"current towerE: "<<(*hitit)->energy()<<endl;
           //  cout<<"maxadc: "<<maxadc<<endl;
           //  cout<<"current adc: "<<(*hitit)->adc()<<endl;
	    softId = (*hitit)->softId(1);
	  //  cout<<" softId = "<<softId<<endl;
	    if(mBemcTriggerSimu->barrelHighTowerAdc(softId)>maxdsmadc) maxdsmadc = mBemcTriggerSimu->barrelHighTowerAdc(softId);
           // cout<<"matchEmc ustc3"<<endl;
	  }
	}
      }//end of cluster loop
      //cout<<"matchEmc skf4"<<endl;
      if(associated) {
	EmcE += (*it)->energy();
	float deltaphi = (*it)->position().phi()-positionBSMDP.phi();
	if(deltaphi>=TMath::Pi()) deltaphi -= TMath::TwoPi();
	if(deltaphi<-TMath::Pi()) deltaphi += TMath::TwoPi();
	float rsmdp = pGeo[3]->Radius();
	float pointz = (*it)->position().z();
	float deltaz = pointz - positionBSMDE.z();
        //cout<<"matchEmc skf"<<endl;
	if(sqrt(deltaphi*deltaphi*rsmdp*rsmdp + deltaz*deltaz)<minDist) {
	  phiDist = deltaphi;
	  zDist = deltaz;
	  if(smdeClusters.size()>=1) nEta = smdeClusters[0]->nHits();
	  if(smdpClusters.size()>=1) nPhi = smdpClusters[0]->nHits();
	  minDist = sqrt(deltaphi*deltaphi*rsmdp*rsmdp + deltaz*deltaz);
	}
      }//end of if
      //cout<<"matchEmc skf5"<<endl;
    }//end of bEmcPoints loop

    //cout<<"matchEmc3"<<endl;
    emcInfo[0] = EmcE;
    emcInfo[1] = nEta;
    emcInfo[2] = nPhi;
    emcInfo[3] = zDist;
    emcInfo[4] = phiDist;
    emcInfo[5] = minDist;
    emcInfo[6] = maxadc;
    emcInfo[7] = maxdsmadc;
    emcInfo[8] = positionBSMDE.pseudoRapidity();
    emcInfo[9] = positionBSMDP.phi();
    emcInfo[10] = maxtowerE;
  }

  
  if(emcInfo[0]>0) {
    //cout<<" emc: "<<flush;
    for(int i=0;i<11;i++){
      //cout<<"  "<<emcInfo[i]<<flush;
    }
    cout<<endl;
  }
 
}


/*void StElectronMcMaker::matchEmc(StTrack* track, StEvent *event, Float_t* emcInfo)
{
  cout<<"star matchEmc"<<endl;
  for(int i=0;i<11;i++) emcInfo[i] = -999.;
  if(!track||track->flag()<=0) return;
  //if(!isGoodTrack(track)) return;
  StEmcCollection* pEmcCol = event->emcCollection();
  //if(track->geometry()->momentum().perp()<1.5) return;
  //cout<<" pass momentum cut"<<endl;
  Int_t mod, eta, sub;
  StThreeVectorD position, momentum;
  StThreeVectorD positionBSMDE, momentumBSMDE;
  StThreeVectorD positionBSMDP, momentumBSMDP;
  Double_t bFld;
  bFld = event->summary()->magneticField()/10.; // bFld in Tesla
  bool ok = false;
  bool okBSMDE = false;
  bool okBSMDP = false;
 
 
  //cout<<"matchEmc1"<<endl;
  StEmcPosition* pEmcPosition = new StEmcPosition();
  StEmcGeom* pGeo[4];
  for(int i=0;i<4;i++) {
     if(i!=1) pGeo[i] = StEmcGeom::getEmcGeom(detname[i].Data());
  }
  if(pEmcPosition){
    ok      = pEmcPosition->projTrack(&position, &momentum, track, bFld, pGeo[0]->Radius());
    okBSMDE = pEmcPosition->projTrack(&positionBSMDE, &momentumBSMDE, track, bFld, pGeo[2]->Radius());
    okBSMDP = pEmcPosition->projTrack(&positionBSMDP, &momentumBSMDP, track, bFld, pGeo[3]->Radius());  
  }

  //cout<<"matchEmc2"<<endl;
  if(ok && okBSMDE && okBSMDP) {
    //cout<<" projection ok"<<endl;
    StSPtrVecEmcPoint& bEmcPoints = pEmcCol->barrelPoints();
    float EmcE = 0.;
    float maxtowerE = -999.;
    float phiDist = -999.;
    float zDist = -999.;
    int nEta = 0;
    int nPhi = 0;
    float minDist = 999.;
    unsigned int maxadc = 0;
    unsigned int maxdsmadc = 0;
    unsigned int softId = -1;
    //cout<<"matchEmc skf"<<endl;

    if(!pGeo[0]) pGeo[0]=StEmcGeom::getEmcGeom(detname[0].Data());
    pGeo[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub);
    //cout<<" (pt,p,eta) "<<track->geometry()->momentum().perp()<<"  "<<track->geometry()->momentum().mag()<<"  "<<track->geometry()->momentum().pseudoRapidity()<<endl;
    for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it!= bEmcPoints.end(); it++) {
      //cout<<"matchEmc skf"<<endl;
      bool associated = false;
      StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
      StPtrVecEmcCluster& smdeClusters = (*it)->cluster(kBarrelSmdEtaStripId);
      StPtrVecEmcCluster& smdpClusters = (*it)->cluster(kBarrelSmdPhiStripId);
      //if(smdeClusters.size()==0 || smdpClusters.size()==0 || bEmcClusters.size()==0) continue;
      //if(smdeClusters[0]==NULL || smdpClusters[0]==NULL || bEmcClusters[0]==NULL) continue;
      if( bEmcClusters.size()==0 || bEmcClusters[0]==NULL) continue;
      for(StPtrVecEmcClusterIterator cIter = bEmcClusters.begin(); cIter != bEmcClusters.end(); cIter++) {
        //cout<<"matchEmc skf"<<endl;
	StPtrVecEmcRawHit& bEmcHits = (*cIter)->hit();
	for(StPtrVecEmcRawHitIterator hIter = bEmcHits.begin(); hIter != bEmcHits.end(); hIter++) {
          //cout<<"matchEmc skf"<<endl;
	  if(mod == (Int_t) (*hIter)->module() && eta == (Int_t) (*hIter)->eta() && sub == (Int_t) (*hIter)->sub()) {
	    //cout<<"(mod,eta,sub,adc,energy)=("<<mod<<", "<<eta<<", "<<sub<<", "<<(*hIter)->adc()<<", "<<(*hIter)->energy()<<")"<<endl;
	    associated = true;
	    break;
	  }
	}//end of raw hits loop
        //cout<<"matchEmc skf"<<endl;
	if(associated) {
	  for(StPtrVecEmcRawHitIterator hitit = bEmcHits.begin(); hitit != bEmcHits.end(); hitit++) {
            //cout<<"matchEmc skf"<<endl;
	    if((*hitit)->energy()>maxtowerE) maxtowerE = (*hitit)->energy();
	    if((*hitit)->adc()>maxadc) maxadc = (*hitit)->adc();
            //cout<<"matchEmc skf"<<endl;
	    softId = (*hitit)->softId(1);
	    //cout<<" softId = "<<softId<<endl;
	    if(mBemcTriggerSimu->barrelHighTowerAdc(softId)>maxdsmadc) maxdsmadc = mBemcTriggerSimu->barrelHighTowerAdc(softId);
            //cout<<"matchEmc skf"<<endl;
	  }
	}
      }//end of cluster loop
      //cout<<"matchEmc skf"<<endl;
      if(associated) {
	EmcE += (*it)->energy();
	float deltaphi = (*it)->position().phi()-positionBSMDP.phi();
	if(deltaphi>=TMath::Pi()) deltaphi -= TMath::TwoPi();
	if(deltaphi<-TMath::Pi()) deltaphi += TMath::TwoPi();
	float rsmdp = pGeo[3]->Radius();
	float pointz = (*it)->position().z();
	float deltaz = pointz - positionBSMDE.z();
        //cout<<"matchEmc skf"<<endl;
	if(sqrt(deltaphi*deltaphi*rsmdp*rsmdp + deltaz*deltaz)<minDist) {
	  phiDist = deltaphi;
	  zDist = deltaz;
	  if(smdeClusters.size()>=1) nEta = smdeClusters[0]->nHits();
	  if(smdpClusters.size()>=1) nPhi = smdpClusters[0]->nHits();
	  minDist = sqrt(deltaphi*deltaphi*rsmdp*rsmdp + deltaz*deltaz);
	}
      }//end of if
    }//end of bEmcPoints loop

    //cout<<"matchEmc3"<<endl;
    emcInfo[0] = EmcE;
    emcInfo[1] = nEta;
    emcInfo[2] = nPhi;
    emcInfo[3] = zDist;
    emcInfo[4] = phiDist;
    emcInfo[5] = minDist;
    emcInfo[6] = maxadc;
    emcInfo[7] = maxdsmadc;
    emcInfo[8] = positionBSMDE.pseudoRapidity();
    emcInfo[9] = positionBSMDP.phi();
    emcInfo[10] = maxtowerE;
  }

  
  if(emcInfo[0]>0) {
    cout<<" emc: "<<flush;
    for(int i=0;i<11;i++){
      cout<<"  "<<emcInfo[i]<<flush;
    }
    cout<<endl;
  }
 
}*/



