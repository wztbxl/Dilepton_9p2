#include "StElectronMcMaker.h"

#include "StChain.h"
#include "TF1.h"
#include "TRandom.h"

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

#include <math.h>
#include "TRandom.h"
#include "Random.h"
#include "RanluxEngine.h"
#include "RandGauss.h"

static RanluxEngine engine;
static RandGauss ranGauss(engine);

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
	mTree->Branch("triggerId",&mElectron.triggerId,"triggerId/I");
	mTree->Branch("mcEvtId",&mElectron.mcEvtId,"mcEvtId/I");
	mTree->Branch("mcVertexX",&mElectron.mcVertexX,"mcVertexX/F");
	mTree->Branch("mcVertexY",&mElectron.mcVertexY,"mcVertexY/F");
	mTree->Branch("mcVertexZ",&mElectron.mcVertexZ,"mcVertexZ/F");
	mTree->Branch("rcEvtId",&mElectron.rcEvtId,"rcEvtId/I");
	mTree->Branch("rcVertexRanking",&mElectron.rcVertexRanking,"rcVertexRanking/F");
	mTree->Branch("rcVertexX",&mElectron.rcVertexX,"rcVertexX/F");
	mTree->Branch("rcVertexY",&mElectron.rcVertexY,"rcVertexY/F");
	mTree->Branch("rcVertexZ",&mElectron.rcVertexZ,"rcVertexZ/F");
	mTree->Branch("rcVpdVz",&mElectron.rcVpdVz,"rcVpdVz/F");
	mTree->Branch("rcRefMult",&mElectron.rcRefMult,"rcRefMult/I");
	mTree->Branch("rcRefMultCorr",&mElectron.rcRefMultCorr,"rcRefMultCorr/F");
	mTree->Branch("gWeight",&mElectron.gWeight,"gWeigth/F");
	mTree->Branch("zdcX",&mElectron.zdcX,"zdcX/F");

	mTree->Branch("nMcTrks",&mElectron.nMcTrks,"nMcTrks/I");
	mTree->Branch("nMcE",&mElectron.nMcE,"nMcE/I");
	mTree->Branch("geantId",mElectron.geantId,"geantId[nMcE]/I");
	mTree->Branch("mcTrkKey",mElectron.mcTrkKey,"mcTrkKey[nMcE]/I");
	mTree->Branch("parentGeantId",mElectron.parentGeantId,"parentGeantId[nMcE]/I");
	mTree->Branch("mcParentTrkKey",mElectron.mcParentTrkKey,"mcParentTrkKey[nMcE]/I");
	mTree->Branch("mcPt",mElectron.mcPt,"mcPt[nMcE]/F");
	mTree->Branch("mcEta",mElectron.mcEta,"mcEta[nMcE]/F");
	mTree->Branch("mcPhi",mElectron.mcPhi,"mcPhi[nMcE]/F");
	mTree->Branch("mcPtFirst",mElectron.mcPtFirst,"mcPtFirst[nMcE]/F");
	mTree->Branch("mcEtaFirst",mElectron.mcEtaFirst,"mcEtaFirst[nMcE]/F");
	mTree->Branch("mcPhiFirst",mElectron.mcPhiFirst,"mcPhiFirst[nMcE]/F");
	mTree->Branch("mcPtLast",mElectron.mcPtLast,"mcPtLast[nMcE]/F");
	mTree->Branch("mcEtaLast",mElectron.mcEtaLast,"mcEtaLast[nMcE]/F");
	mTree->Branch("mcPhiLast",mElectron.mcPhiLast,"mcPhiLast[nMcE]/F");
	mTree->Branch("rcTrkKey",mElectron.rcTrkKey,"rcTrkKey[nMcE]/I");
	mTree->Branch("rcPtFirst",mElectron.rcPtFirst,"rcPtFirst[nMcE]/F");
	mTree->Branch("rcEtaFirst",mElectron.rcEtaFirst,"rcEtaFirst[nMcE]/F");
	mTree->Branch("rcPhiFirst",mElectron.rcPhiFirst,"rcPhiFirst[nMcE]/F");
	mTree->Branch("rcPtLast",mElectron.rcPtLast,"rcPtLast[nMcE]/F");
	mTree->Branch("rcEtaLast",mElectron.rcEtaLast,"rcEtaLast[nMcE]/F");
	mTree->Branch("rcPhiLast",mElectron.rcPhiLast,"rcPhiLast[nMcE]/F");
	mTree->Branch("rcNHitsFit",mElectron.rcNHitsFit,"rcNHitsFit[nMcE]/I");
	mTree->Branch("rcNHitsPoss",mElectron.rcNHitsPoss,"rcNHitsPoss[nMcE]/I");
	mTree->Branch("rcNHitsDedx",mElectron.rcNHitsDedx,"rcNHitsDedx[nMcE]/I");
	mTree->Branch("rcNHitsCommon",mElectron.rcNHitsCommon,"rcNHitsCommon[nMcE]/I");
	mTree->Branch("rcDedx",mElectron.rcDedx,"rcDedx[nMcE]/F");
	mTree->Branch("rcNSigmaE",mElectron.rcNSigmaE,"rcNSigmaE[nMcE]/F");
	mTree->Branch("rcNSigmaPi",mElectron.rcNSigmaPi,"rcNSigmaPi[nMcE]/F");
	mTree->Branch("rcNSigmaK",mElectron.rcNSigmaK,"rcNSigmaK[nMcE]/F");
	mTree->Branch("rcNSigmaP",mElectron.rcNSigmaP,"rcNSigmaP[nMcE]/F");
	mTree->Branch("rcDca",mElectron.rcDca,"rcDca[nMcE]/F");
	mTree->Branch("rcDca1",mElectron.rcDca1,"rcDca1[nMcE]/F");
	mTree->Branch("rcDca2",mElectron.rcDca2,"rcDca2[nMcE]/F");
	mTree->Branch("rcDca3",mElectron.rcDca3,"rcDca3[nMcE]/F");
	mTree->Branch("rcDca4",mElectron.rcDca4,"rcDca4[nMcE]/F");
	mTree->Branch("rcDca5",mElectron.rcDca5,"rcDca5[nMcE]/F");
	mTree->Branch("rcDca6",mElectron.rcDca6,"rcDca6[nMcE]/F");

	// - zero all pointers defined in the header file
	mRcEvent=0;
	mMcEvent=0;
	mMcVertex.setX(0.); mMcVertex.setY(0.); mMcVertex.setZ(0.);
	mAssocMaker=0;
	mRcTrackMap=0;
	mMcTrackMap=0;
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
	Clear("C");

	//Zhen change it ,no centrality definition
	//intialize the StRefMultCorr
	//mRefMultCorr = new StRefMultCorr();

	return kStOK;
}
//_________________________________________________
Int_t StElectronMcMaker::Make()
{
	Clear("C");
	memset(&mElectron,0,sizeof(mElectron)); //initialize the mElectron struct

	mMcEvent = (StMcEvent*)GetDataSet("StMcEvent");
	if(!mMcEvent) return kStErr;

	mRcEvent=(StEvent *)GetInputDS("StEvent");
	if(!mRcEvent){
		LOG_WARN << "No StEvent!" << endm;
		return kStOK;
	}

	mAssocMaker = (StAssociationMaker*)GetMaker("StAssociationMaker");
	if(mAssocMaker){
		mRcTrackMap = mAssocMaker->rcTrackMap();
		mMcTrackMap = mAssocMaker->mcTrackMap();
	}

	StTriggerIdCollection *trgIdColl = mRcEvent->triggerIdCollection();
	if(trgIdColl && trgIdColl->nominal() && //zaochen's trigger
			(trgIdColl->nominal()->isTrigger(610001) || 
			 trgIdColl->nominal()->isTrigger(610011) ||
			 trgIdColl->nominal()->isTrigger(610031) ||
			 trgIdColl->nominal()->isTrigger(610041) ||
			 trgIdColl->nominal()->isTrigger(610051) ||
			 trgIdColl->nominal()->isTrigger(610021)) 
	  ){
	}else{
		cout << "Not a minbias trigger event!" << endl;
		return kStOK;
	}
	if(trgIdColl->nominal()->isTrigger(610001)) mElectron.triggerId =610001;
	if(trgIdColl->nominal()->isTrigger(610011)) mElectron.triggerId =610011;
	if(trgIdColl->nominal()->isTrigger(610021)) mElectron.triggerId =610021;
	if(trgIdColl->nominal()->isTrigger(610031)) mElectron.triggerId =610031;
	if(trgIdColl->nominal()->isTrigger(610041)) mElectron.triggerId =610041;
	if(trgIdColl->nominal()->isTrigger(610051)) mElectron.triggerId =610051;

	mElectron.runId = mRcEvent->info()->runId();
	mElectron.rcEvtId = mRcEvent->info()->id();

	StPrimaryVertex  *pVertex = mRcEvent->primaryVertex();
	StThreeVectorF rcVertex(-999.,-999.,-999.);
	if(pVertex){
		mElectron.rcVertexRanking = pVertex->ranking();
		rcVertex = pVertex->position();
		if(fabs(rcVertex.x())<1.e-5 
				&& fabs(rcVertex.y())<1.e-5
				&& fabs(rcVertex.z())<1.e-5){
			return kStOK;//reject event with bad vertex
		}else{
			mElectron.rcVertexX = rcVertex.x();
			mElectron.rcVertexY = rcVertex.y();
			mElectron.rcVertexZ = rcVertex.z();
		}
	}else{
		return kStOK;
	}

	mElectron.rcVpdVz = -999.;
	if(mRcEvent->btofCollection()){
		StBTofHeader *btofHeader = mRcEvent->btofCollection()->tofHeader();
		if(btofHeader){
			mElectron.rcVpdVz = btofHeader->vpdVz();
		}
	}

	mElectron.mcEvtId = mMcEvent->eventNumber();
	if(mMcEvent->primaryVertex()){
		mMcVertex = mMcEvent->primaryVertex()->position();
		mElectron.mcVertexX = mMcVertex.x();
		mElectron.mcVertexY = mMcVertex.y();
		mElectron.mcVertexZ = mMcVertex.z();
	}else{
		return kStOK;
	}

	mElectron.rcRefMult = uncorrectedNumberOfPositivePrimaries(*mRcEvent) + uncorrectedNumberOfNegativePrimaries(*mRcEvent);
	mElectron.zdcX = mRcEvent->runInfo()->zdcCoincidenceRate();

	//wangzhen change it for no RefMultCorr
	//beacuse there is no centrality definition for 54GeV
	/*
	if(fabs(mElectron.mcVertexZ)<30.){ //refMultCorr requirment
		mRefMultCorr->init(mElectron.runId);
		mRefMultCorr->initEvent(mElectron.rcRefMult,mElectron.mcVertexZ,mElectron.zdcX);//use mc vertexz to do the refMult correction
		mElectron.rcRefMultCorr = mRefMultCorr->getRefMultCorr() ;
		//mElectron.gWeight = mRefMultCorr->getWeight() ;
		mElectron.gWeight = 1; //now, the correction parameters are unvaliable.
	}
	//cout<<"refMult:"<<mElectron.rcRefMult<<"  refMultCorr:"<<mElectron.rcRefMultCorr<<"  gWeight:"<<mElectron.gWeight<<endl;
	*/

	//smear the vertex
	StThreeVectorF rcVertexSmear[6];
	const float sigma[6] = {0.02, 0.05, 0.10, 0.15, 0.20, 0.25};
	for(int i=0;i<6;i++) {
		rcVertexSmear[i].set(rcVertex.x()+ranGauss.shoot()*sigma[i], rcVertex.y()+ranGauss.shoot()*sigma[i], rcVertex.z());
	}

	//looking for electon
	StSPtrVecMcTrack mcTracks = mMcEvent->tracks();
	mElectron.nMcTrks = mcTracks.size();
	//cout<<"# of mcTracks:"<<mcTracks.size()<<endl;
	
	int nMcE = 0;
	for(int i=0;i<(int)mcTracks.size();i++){
		StMcTrack *mcTrack = dynamic_cast<StMcTrack *>(mcTracks[i]);
		if(!mcTrack) continue;
		if(mcTrack->key()==0 && mcTrack->geantId()==0) continue; //not geant tracks

		mElectron.geantId[nMcE] = mcTrack->geantId();
		if(mElectron.geantId[nMcE]!=2 
				&& mElectron.geantId[nMcE]!=3) continue; //not a electron;
		mElectron.mcTrkKey[nMcE] = mcTrack->key();

		//cout<<"mTrackId:"<<i<<endl;

		const StMcTrack *mcParentTrack = mcTrack->parent();
		if(mcParentTrack){
			if(mcParentTrack->geantId()!=0) continue; 
			mElectron.parentGeantId[nMcE] =  mcParentTrack->geantId();
			mElectron.mcParentTrkKey[nMcE] = mcParentTrack->key();
		}

		//if(!mcTrack->startVertex() ||
		//		(mcTrack->startVertex()->position()-mMcVertex).mag()>1.e-5 ) continue;

		StThreeVectorF mcMom = mcTrack->momentum();
		//if(fabs(mcMom.pseudoRapidity())>1.5) continue;

		mElectron.mcPt[nMcE] = mcMom.perp();
		mElectron.mcEta[nMcE] = mcMom.pseudoRapidity();
		mElectron.mcPhi[nMcE] = mcMom.phi();

		StPtrVecMcTpcHit& mcTpcHit = mcTrack->tpcHits();
		if(mcTpcHit.size()>0){
			mElectron.mcPtFirst[nMcE] = mcTpcHit[0]->localMomentum().perp();
			mElectron.mcEtaFirst[nMcE] = mcTpcHit[0]->localMomentum().pseudoRapidity();
			mElectron.mcPhiFirst[nMcE] = mcTpcHit[0]->localMomentum().phi();
			mElectron.mcPtLast[nMcE] = mcTpcHit[mcTpcHit.size()-1]->localMomentum().perp();
			mElectron.mcEtaLast[nMcE] = mcTpcHit[mcTpcHit.size()-1]->localMomentum().pseudoRapidity();
			mElectron.mcPhiLast[nMcE] = mcTpcHit[mcTpcHit.size()-1]->localMomentum().phi();
		} 

		StGlobalTrack *rcGTrack = 0;
		StPrimaryTrack *rcPTrack = 0;

		// find the associate track
		if(mAssocMaker){
			if(mMcTrackMap){
				StTrackPairInfo* mTrackPairInfo = findBestMatchedGlobal(mcTrack);
				//cout<<"mTrackPairInfo Pointer:"<<mTrackPairInfo<<endl;
				if(mTrackPairInfo){
					rcGTrack = const_cast<StGlobalTrack*>(mTrackPairInfo->partnerTrack());
					if(rcGTrack) rcPTrack = dynamic_cast<StPrimaryTrack *>(rcGTrack->node()->track(primary));
					mElectron.rcNHitsCommon[nMcE] = mTrackPairInfo->commonTpcHits();
				}
			}
		}
		cout << "nCommonHits = " <<mElectron.rcNHitsCommon[nMcE]<< endl;
		//if(mElectron.rcNHitsCommon[nMcE]>10) cout<<"gTrkPointer:"<<rcGTrack<<"  pTrkPointer:"<<rcPTrack<<endl;

		if(rcPTrack){
			mElectron.rcTrkKey[nMcE] = rcPTrack->key();
			mElectron.rcPtFirst[nMcE] = rcPTrack->geometry()->momentum().perp();
			mElectron.rcEtaFirst[nMcE] = rcPTrack->geometry()->momentum().pseudoRapidity();
			mElectron.rcPhiFirst[nMcE] = rcPTrack->geometry()->momentum().phi();
			mElectron.rcPtLast[nMcE] = rcPTrack->outerGeometry()->momentum().perp();
			mElectron.rcEtaLast[nMcE] = rcPTrack->outerGeometry()->momentum().pseudoRapidity();
			mElectron.rcPhiLast[nMcE] = rcPTrack->outerGeometry()->momentum().phi();
			mElectron.rcNHitsFit[nMcE] = rcPTrack->fitTraits().numberOfFitPoints(kTpcId);
			mElectron.rcNHitsPoss[nMcE] = rcPTrack->numberOfPossiblePoints(kTpcId);

			static StTpcDedxPidAlgorithm pidAlgorithm;
			static StPionPlus* Pion = StPionPlus::instance();
			static StKaonPlus* Kaon = StKaonPlus::instance();
			static StProton* Proton = StProton::instance();
			static StElectron* Electron = StElectron::instance();
			const StParticleDefinition* p1 = rcPTrack->pidTraits(pidAlgorithm);
			if(p1 && pidAlgorithm.traits()){
				mElectron.rcNHitsDedx[nMcE] = pidAlgorithm.traits()->numberOfPoints();
				mElectron.rcDedx[nMcE] = pidAlgorithm.traits()->mean()*1.e6;
				mElectron.rcNSigmaE[nMcE] = pidAlgorithm.numberOfSigma(Electron);
				mElectron.rcNSigmaPi[nMcE] = pidAlgorithm.numberOfSigma(Pion);
				mElectron.rcNSigmaK[nMcE] = pidAlgorithm.numberOfSigma(Kaon);
				mElectron.rcNSigmaP[nMcE] = pidAlgorithm.numberOfSigma(Proton);
			}
			if(rcGTrack){
				mElectron.rcDca[nMcE] = rcGTrack->geometry()->helix().distance(rcVertex);
				mElectron.rcDca1[nMcE] = rcGTrack->geometry()->helix().distance(rcVertexSmear[0]);
				mElectron.rcDca2[nMcE] = rcGTrack->geometry()->helix().distance(rcVertexSmear[1]);
				mElectron.rcDca3[nMcE] = rcGTrack->geometry()->helix().distance(rcVertexSmear[2]);
				mElectron.rcDca4[nMcE] = rcGTrack->geometry()->helix().distance(rcVertexSmear[3]);
				mElectron.rcDca5[nMcE] = rcGTrack->geometry()->helix().distance(rcVertexSmear[4]);
				mElectron.rcDca6[nMcE] = rcGTrack->geometry()->helix().distance(rcVertexSmear[5]);
			}
		}
		nMcE++;
	}
	mElectron.nMcE = nMcE;

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
