#include "iostream.h"

class     StChain;
class StRefMultCorr;
StChain  *chain=0;

Int_t iEvt=0,istat=0,nEvents=0;
//ProtoTypes

// ------------------ Here is the actual method -----------------------------------------
void doEvent(const Int_t nEvents = 1000, 
		const char *fMcFile = "/star/embed/embedding/27GeV_production_2018/Electron_200_20193201/P19ib.SL19b/2018/150/19150006/st_physics_adc_19150006_raw_4000016.event.root",

		const Char_t *outDir="./test")
{
	// First load some shared libraries we need
	if (gClassTable->GetID("TTable") < 0) {
		gSystem->Load("libStar");
		gSystem->Load("libPhysics");
	}  
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StTpcDb");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("StEmcUtil");
	gSystem->Load("StEEmcUtil");
	gSystem->Load("StMcEvent");
	gSystem->Load("StMcEventMaker");
	gSystem->Load("StDaqLib");

	gSystem->Load("StDbBroker");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("St_db_Maker");

	gSystem->Load("StAssociationMaker");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("StElectronMcMaker");

	chain  = new StChain("StChain");

	StIOMaker* ioMaker = new StIOMaker();
	ioMaker->SetFile(fMcFile);
	ioMaker->SetIOMode("r");
	ioMaker->SetBranch("*",0,"0");             //deactivate all branches
	ioMaker->SetBranch("geantBranch",0,"r");   //activate geant Branch
	ioMaker->SetBranch("eventBranch",0,"r");   //activate event Branch

	StMcEventMaker*     mcEventReader = new StMcEventMaker; // Make an instance..
	mcEventReader->doPrintEventInfo = false;
	mcEventReader->doPrintMemoryInfo = false;

	StAssociationMaker* associator    = new StAssociationMaker;
	cout<<"created new StAssociationMaker"<<endl;
	associator->useIdAssoc();
	associator->useInTracker();

	StMcParameterDB* parameterDB = StMcParameterDB::instance();
	parameterDB->setReqCommonHitsTpc(10); // Require 10 hits in common for tracks

    TString outputname = fMcFile;
	TString embedrun = fMcFile;
	int embedRunIndex = embedrun.Index("27GeV_production_2018",0);
	embedrun.Remove(0,embedRunIndex);
	embedRunIndex = embedrun.Index("P19ib",0);
	embedrun.Remove(embedRunIndex);
	cout<<"embedrun:   "<<embedrun<<endl;
	int fileBeginIndex = outputname.Index("st_physics_adc",0);
	outputname.Remove(0,fileBeginIndex+15);
	cout<<"output: "<< outputname<<endl;
	outputname.Prepend("_");
	outputname.Prepend(embedrun);
	outputname.Prepend("emb");
	outputname.ReplaceAll("/","x");
	outputname.ReplaceAll("geant","myminimc");
	outputname.Prepend("/");
	outputname.Prepend(outDir);
	cout<<"Output: "<<outputname<<endl;

	StElectronMcMaker *electronMcMaker = new StElectronMcMaker("electron",outputname);

	// Initialize chain
	Int_t iInit = chain->Init();
	if (iInit) chain->Fatal(iInit,"on init");
	chain->PrintInfo();

	// Event loop
	int istat = 0, i = 1;
EventLoop: if(i <= nEvents && istat != 2) {

		cout << endl << "============================ Event " << i
			<< " start ============================" << endl;

		chain->Clear();
		istat = chain->Make(i);
		if (istat == 2) 
		{cout << "Last  event processed. Status = " << istat << endl;}
		if (istat == 3) 
		{cout << "Error event processed. Status = " << istat << endl;}

		//gObjectTable->Print();
		i++;
		goto EventLoop;
	}

	i--;
	cout<<endl<<"============================ Event "<<i<<" finish ============================"<<endl;

	// Chain Finish
	chain->Finish();
	delete chain;
}





