#include <TSystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void load(){
	//Load all the System libraries
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StTpcDb");
	gSystem->Load("StEvent");
	gSystem->Load("StMcEvent");
	gSystem->Load("StMcEventMaker");
	gSystem->Load("StDaqLib");
	gSystem->Load("libgen_Tables");
	gSystem->Load("libsim_Tables");
	gSystem->Load("libglobal_Tables");
	gSystem->Load("StMagF");

	gSystem->Load("St_g2t.so");
	gSystem->Load("St_geant_Maker.so");
	gSystem->Load("StAssociationMaker");
	gSystem->Load("StMcAnalysisMaker");
	gSystem->Load("libgeometry_Tables");   
	gSystem->Load("StTriggerUtilities");

	gSystem->Load("StEmcUtil");
	gSystem->Load("StEmcRawMaker");
	gSystem->Load("StEmcADCtoEMaker");
	gSystem->Load("StPreEclMaker");
	gSystem->Load("StEpcMaker");
	gSystem->Load("StEmcSimulatorMaker");

	gSystem->Load("StDbLib");
	gSystem->Load("StDbUtilities");
	gSystem->Load("StDbBroker");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("St_db_Maker");

	gSystem->Load("StMtdHitMaker");
	gSystem->Load("StMtdUtil");
	gSystem->Load("StMtdMatchMaker");
	gSystem->Load("StMtdCalibMaker");
	gSystem->Load("StBTofUtil");
	gSystem->Load("StVpdCalibMaker");

	gSystem->Load("StPicoEvent");
	gSystem->Load("StPicoDstMaker");
	// gSystem->Load("StPicoEvent");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("StMiniTreeMaker");
}

void doEvent(Int_t nEvents=-1, const Char_t *inputFile="test.list", const TString outputFile="test/test.root", const Bool_t debug = kFALSE)
{
	load();

	StChain *chain = new StChain();
	chain->SetDebug(0);
	StMuDebug::setLevel(0); // switch of some debug output

	StMuTimer timer;
	timer.start();

	StPicoDstMaker  *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst"); 

	St_db_Maker *dbMk = new St_db_Maker("StarDb", "MySQL:StarDb", "$STAR/StarDb","StarDb");
	// dbMk->SetDateTime(20170101,0); //for run17 picoDst

	cout<<"before StMiniTreeMaker;"<<endl;
	StMiniTreeMaker *miniTreeMaker = new StMiniTreeMaker();
	miniTreeMaker->setOutFileName(outputFile);
	miniTreeMaker->setUseDefaultVtx(1);
	// miniTreeMaker->setMaxVtxR(2.5);//for QA
	miniTreeMaker->setMaxVtxR(3.);//for QA
	// miniTreeMaker->setMaxVtxR(3.0);//for pure sample
	miniTreeMaker->setMaxVtxZ(70); // for QA
	// miniTreeMaker->setMaxVtxZ(60.);//for pure sample
	miniTreeMaker->setMaxVzDiff(13.);
	miniTreeMaker->setMinNHitsFitRatio(0.5);//for QA
	// miniTreeMaker->setMinNHitsFitRatio(0.52);//for pure sample
	miniTreeMaker->setMaxTrackEta(2);//for QA
	// miniTreeMaker->setMaxTrackEta(1.);//for pure sample
	miniTreeMaker->setMinTrackPt(0.01); // for QA
	// miniTreeMaker->setMinTrackPt(0.2);
	miniTreeMaker->setMinNHitsFit(5); // for QA
	// miniTreeMaker->setMinNHitsFit(20);//for pure sample
	miniTreeMaker->setMinNHitsDedx(5);
  // miniTreeMaker->setParticle(1);// 1 is pion,2 is kaon,3 is proton
	miniTreeMaker->setMaxBeta2TOF(0.05);
	//miniTreeMaker->setPrintMemory(1);
	//miniTreeMaker->setPrintCpu(1);
	miniTreeMaker->setPrintConfig(1);
	//if(name.find("st_hlt")!=std::string::npos)
	//	miniTreeMaker->setStreamName("st_hlt");
	//if(name.find("st_mtd")!=std::string::npos)
	//	miniTreeMaker->setStreamName("st_mtd");
	//if(name.find("st_physics")!=std::string::npos)
	//	miniTreeMaker->setStreamName("st_physics");
	if(debug)
		miniTreeMaker->SetDebug(1);

	cout<<"before chain->Init();"<<endl;
	if(chain->Init()==kStERR) return;
	cout<<"chain->Init();"<<endl;

	if(nEvents<0){
		nEvents = picoMaker->chain()->GetEntries();
    cout << nEvents << endl;
	}

	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;

	for(Int_t i=0; i<nEvents; i++) {
		if(debug) {
			cout<<endl;
			cout<<"Working on eventNumber "<< i <<endl;
		} else {
			if(i%1000==0)
				cout << "Working on eventNumber " << i << endl;
		}

    // cout << "Working on eventNumber " << i << endl;
		chain->Clear();
		int iret = chain->Make(i);

		if(iret) { cout << "Bad return code!" << iret << endl; break;}
	}

	chain->Finish();
	delete chain;

	timer.stop();
	cout << "Total time = " << timer.elapsedTime() << " s" << endl;

	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
}
