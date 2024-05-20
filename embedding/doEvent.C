#include "iostream.h"

class StChain;
class StRefMultCorr;
StChain *chain = 0;

Int_t iEvt = 0, istat = 0, nEvents = 0;
// ProtoTypes

// ------------------ Here is the actual method -----------------------------------------
void doEvent(const Int_t nEvents = 2000,
			 const char *fMcFile = "/star/data105/embedding/production_7p7GeV_2021/Electron_100_20224105/P22ib.SL22b/2021/037/22037024/st_physi     cs_adc_22037024_raw_0000000.MuDst.root",
			//  const char *fMcFile = "/star/data105/embedding/production_9p2GeV_2020/Electron_100_20233801/P23ia.SL23a/2020/219/21219011/st_physics_adc_21219011_raw_7000003.MuDst.root",
			 // const char *fMcFile = "/star/embed/embedding/production_isobar_2018/ElectronHighPt_200_20214218/P20ic.SL20c/2018/129/19129013/st_physics_adc_19129013_raw_2500008.MuDst.root",
			 // const char *fMcFile = "/star/embed/embedding/production_isobar_2018/JpsiMB_200_20215101/P20ic.SL20c/2018/129/19129013/st_physics_adc_19129013_raw_2500008.MuDst.root",
			 // const char *fMcFile = "/star/data09/reco/production_pp500_2022/ReversedFullField/dev/2022/063/23063039/st_fwd_adc_23063039_raw_1000008.MuDst.root",

			 const Char_t *outDir = "./test")
{
	// First load some shared libraries we need
	if (gClassTable->GetID("TTable") < 0)
	{
		gSystem->Load("libStar");
		gSystem->Load("libPhysics");
	}
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	cout << "skf" << endl;
	gROOT->Macro("loadMuDst.C");
	gSystem->Load("StarMagField");
	gSystem->Load("StMagF");
	gSystem->Load("StTpcDb");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("StEEmcUtil");
	gSystem->Load("StEmcUtil");
	gSystem->Load("StEEmcDbMaker");
	gSystem->Load("StMcEvent");
	gSystem->Load("StMcEventMaker");
	gSystem->Load("StDaqLib");
	gSystem->Load("St_g2t.so");
	gSystem->Load("St_geant_Maker.so");
	gSystem->Load("StAssociationMaker");
	gSystem->Load("StMcAnalysisMaker");
	gSystem->Load("StDbLib");
	gSystem->Load("StDbUtilities");
	gSystem->Load("StDbBroker");
	gSystem->Load("St_db_Maker");
	gSystem->Load("libgeometry_Tables");
	gSystem->Load("libgen_Tables");
	gSystem->Load("libsim_Tables");
	gSystem->Load("libglobal_Tables");
	gSystem->Load("StEmcRawMaker");
	gSystem->Load("StEmcADCtoEMaker");
	gSystem->Load("StPreEclMaker");
	gSystem->Load("StEpcMaker");
	gSystem->Load("StEmcSimulatorMaker");
	gSystem->Load("StTriggerUtilities");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("StBTofUtil");
	gSystem->Load("StVpdCalibMaker");
	// gSystem->Load("RTS");

	// gSystem->Load("StDetectorDbMaker");
	// gSystem->Load("StDbUtilities");

	gSystem->Load("StMuDstMaker");
	gSystem->Load("StMuEvent");
	gSystem->Load("StMuTrack");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("StEvent");
	gSystem->Load("StEventMaker");
	gSystem->Load("StElectronMcMaker");

	cout << "skf" << endl;

	// chain  = new StChain("StChain");
	// chain  = new StChain("bfc");
	chain = new StChain();

	StMuDstMaker *muDstMaker = new StMuDstMaker(0, 1, "", fMcFile, "", 100000, "MuDst");
	muDstMaker->SetStatus("*", 0);
	muDstMaker->SetStatus("MuEvent", 1);
	muDstMaker->SetStatus("PrimaryVertices", 1);
	muDstMaker->SetStatus("PrimaryTracks", 1);
	muDstMaker->SetStatus("GlobalTracks", 1);
	muDstMaker->SetStatus("CovGlobTrack", 1);
	muDstMaker->SetStatus("MCAll", 1);
	muDstMaker->SetStatus("BTofHeader", 1);
	muDstMaker->SetDebug(0);

	StIOMaker *ioMaker = new StIOMaker("IO");
	ioMaker->SetFile(fMcFile);
	ioMaker->SetIOMode("r");
	ioMaker->SetBranch("*", 0, "0");		   // deactivate all branches
	ioMaker->SetBranch("geantBranch", 0, "r"); // activate geant Branch
	ioMaker->SetBranch("eventBranch", 0, "r"); // activate event Branch

	// StMcEventMaker *mcEventReader = new StMcEventMaker(); // Make an instance..
	// mcEventReader->doPrintEventInfo = false;
	// mcEventReader->doPrintMemoryInfo = false;
	/*StMuDstMaker* muDstMaker = new StMuDstMaker(0,1,"",fMcFile,"",100000,"MuDst");
	muDstMaker->SetStatus("*",0);
	muDstMaker->SetStatus("MuEvent",1);
	muDstMaker->SetStatus("PrimaryVertices",1);
	muDstMaker->SetStatus("PrimaryTracks",1);
	muDstMaker->SetStatus("GlobalTracks",1);
	muDstMaker->SetStatus("CovGlobTrack",1);
	muDstMaker->SetStatus("MCAll",1);
	muDstMaker->SetStatus("BTofHeader",1);
	muDstMaker->SetDebug(0);*/

	// cout<<"skf"<<endl;
	StAssociationMaker *associator = new StAssociationMaker();
	cout << "created new StAssociationMaker" << endl;
	associator->useIdAssoc();
	associator->useInTracker();

	StMcParameterDB *parameterDB = StMcParameterDB::instance();
	parameterDB->setReqCommonHitsTpc(10); // Require 10 hits in common for tracks


	///star/data105/embedding/production_7p7GeV_2021/Electron_100_20224105/P22ib.SL22b/2021/056/22056005/st_physics_adc_22056005_raw_0000000.MuDst.root
	TString outputname = fMcFile;
	TString embedrun = fMcFile;
	int embedRunIndex = embedrun.Index("production_9p2GeV_2020", 0);
	embedrun.Remove(0, embedRunIndex);
	embedRunIndex = embedrun.Index("P23ia", 0);
	embedrun.Remove(embedRunIndex);
	cout << "embedrun:   " << embedrun << endl;
	int fileBeginIndex = outputname.Index("st_physics_adc", 0);
	outputname.Remove(0, fileBeginIndex + 15);
	cout << "output: " << outputname << endl;
	outputname.Prepend("_");
	outputname.Prepend(embedrun);
	outputname.Prepend("emb");
	outputname.ReplaceAll("/", "x");
	outputname.ReplaceAll("MuDst", "myminimc");
	outputname.Prepend("/");
	outputname.Prepend(outDir);
	cout << "Output: " << outputname << endl;

	StElectronMcMaker *electronMcMaker = new StElectronMcMaker("electron", outputname);

	// Initialize chain
	Int_t iInit = chain->Init();
	if (iInit)
		chain->Fatal(iInit, "on init");
	chain->PrintInfo();

	// Event loop
	int istat = 0, i = 1;
EventLoop:
	if (i <= nEvents && istat != 2)
	{

		cout << endl
			 << "============================ Event " << i
			 << " start ============================" << endl;

		chain->Clear();
		cout << "ustc" << endl;
		istat = chain->Make(i);
		cout << "ustc" << endl;
		if (istat == 2)
		{
			cout << "Last  event processed. Status = " << istat << endl;
		}
		if (istat == 3)
		{
			cout << "Error event processed. Status = " << istat << endl;
		}

		// gObjectTable->Print();
		i++;
		goto EventLoop;
	}

	i--;
	cout << endl
		 << "============================ Event " << i << " finish ============================" << endl;

	// Chain Finish
	chain->Finish();
	delete chain;
}
