#ifndef StElectronMcMaker_HH
#define StElectronMcMaker_HH

#ifndef StMaker_H
#include "StMaker.h"
#endif


class StTrack;
class StMcTrack;
class StMcEvent;
class TRandom;
class StFileI;
class TChain;
class TClonesArray;

//-- add class
class StEvent;
class StMcEvent;
class StTrack;
class StGlobalTrack;
class StMcTrack;
class StAssociationMaker;
class StRefMultCorr;
#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"
#include "StThreeVectorF.hh"
#include "StRefMultCorr/StRefMultCorr.h"

#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"

#include <vector>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

const int mMax = 500;
class StElectronMcMaker : public StMaker {

	protected:

	public:

		StMaker* currentChain;
		StElectronMcMaker(const char* name = "StElectronMcMaker",
				const char* file = "undefined");
		~StElectronMcMaker();
		void  Clear(const char* opt="");
		Int_t Init();
		Int_t Make();
		Int_t Finish();

	private:

		struct ElectronEvent {
			int   runId;
			int   triggerId;
			int   mcEvtId;
			float mcVertexX;
			float mcVertexY;
			float mcVertexZ;
			int   rcEvtId;
			float rcVertexRanking;
			float rcVertexX;
			float rcVertexY;
			float rcVertexZ;
			float rcVpdVz;
			int   rcRefMult;
			float rcRefMultCorr;
			float gWeight;
			float zdcX;

			int   nMcTrks;
			int   nMcE;
			int   geantId[mMax];
			int   mcTrkKey[mMax];
			int   parentGeantId[mMax];
			int   mcParentTrkKey[mMax];
			float mcPt[mMax];
			float mcEta[mMax];
			float mcPhi[mMax];
			float mcPtFirst[mMax];
			float mcEtaFirst[mMax];
			float mcPhiFirst[mMax];
			float mcPtLast[mMax];
			float mcEtaLast[mMax];
			float mcPhiLast[mMax];
			int   rcTrkKey[mMax];
			float rcPtFirst[mMax];
			float rcEtaFirst[mMax];
			float rcPhiFirst[mMax];
			float rcPtLast[mMax];
			float rcEtaLast[mMax];
			float rcPhiLast[mMax];
			int   rcNHitsFit[mMax];            
			int   rcNHitsPoss[mMax];
			int   rcNHitsDedx[mMax];
			int   rcNHitsCommon[mMax];
			float rcDedx[mMax];
			float rcNSigmaE[mMax];
			float rcNSigmaPi[mMax];
			float rcNSigmaK[mMax];
			float rcNSigmaP[mMax];
			float rcDca[mMax];
			float rcDca1[mMax]; // vertex XY smear by 0.02 cm
			float rcDca2[mMax]; // vertex XY smear by 0.05 cm
			float rcDca3[mMax]; // vertex XY smear by 0.10 cm
			float rcDca4[mMax]; // vertex XY smear by 0.15 cm
			float rcDca5[mMax]; // vertex XY smear by 0.20 cm
			float rcDca6[mMax]; // vertex XY smear by 0.25 cm
		};
		ElectronEvent mElectron;

		TFile* mOutputFile;
		TTree* mTree;

		StEvent*       mRcEvent;
		StMcEvent*     mMcEvent;

		StThreeVectorF       mMcVertex;
		StAssociationMaker*  mAssocMaker;
		rcTrackMapType*      mRcTrackMap;
		mcTrackMapType*      mMcTrackMap;
		StRefMultCorr*       mRefMultCorr;

		StTrackPairInfo*     findBestMatchedGlobal(StMcTrack*);

#ifndef ST_NO_TEMPLATE_DEF_ARGS
		typedef vector<Int_t> idVector;
#else
		typedef vector<Int_t,allocator<Int_t>> idVector;
#endif
		typedef idVector::iterator idVectorIter;

		//    idVector Kaons;

		// the following is a ROOT macro  that is needed in all ROOT accessible code
		ClassDef(StElectronMcMaker, 1)

};

#endif

