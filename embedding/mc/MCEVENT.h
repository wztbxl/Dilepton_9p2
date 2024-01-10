//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  9 20:40:36 2024 by ROOT version 5.34/30
// from TTree mcT/mcT
// found on file: embproduction_7p7GeV_2021xxElectron_100_20224105xP22ib.SL22bx2021x121x22121011xst_physics_adc_22121011_raw_7500003.myminimc.root_22121011_raw_7500003.myminimc.root
//////////////////////////////////////////////////////////

#ifndef MCEVENT_h
#define MCEVENT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MCEVENT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runId;
   Int_t           triggerId[4];
   Int_t           muEvtId;
   Float_t         muPriVertexX;
   Float_t         muPriVertexY;
   Float_t         muPriVertexZ;
   Float_t         muVpdVz;
   Int_t           muRefMult;
   Int_t           mugRefMult;
   Int_t           munBtofMatch;
   Float_t         zdcX;
   Int_t           Centrality_16;
   Int_t           Centrality_9;
   Float_t         refmult_corr;
   Float_t         reweight;
   Int_t           nMcTrks;
   Int_t           nMcVertices;
   Int_t           nMcE;
   Int_t           nRcE;
   Int_t           nMcTrack;
   Int_t           geantId[10];   //[nMcE]
   Int_t           mcTrackId[10];   //[nMcE]
   Short_t         mcCharge[10];   //[nMcE]
   Float_t         mcTrack_Vr[10];   //[nMcE]
   Int_t           mcTrack_par_Geantid[10];   //[nMcE]
   Float_t         mcPt[10];   //[nMcE]
   Float_t         mcEta[10];   //[nMcE]
   Float_t         mcPhi[10];   //[nMcE]
   Int_t           rcidtruth[10];   //[nRcE]
   Int_t           rcQaTruth[10];   //[nRcE]
   Int_t           rcflag[10];   //[nRcE]
   Short_t         rcCharge[10];   //[nRcE]
   Float_t         rcPt[10];   //[nRcE]
   Float_t         rcEta[10];   //[nRcE]
   Float_t         rcPhi[10];   //[nRcE]
   Int_t           rcNHitsFit[10];   //[nRcE]
   Int_t           rcNHitsPoss[10];   //[nRcE]
   Int_t           rcNHitsDedx[10];   //[nRcE]
   Float_t         rcDedx[10];   //[nRcE]
   Float_t         rcNSigmaE[10];   //[nRcE]
   Float_t         rcNSigmaPi[10];   //[nRcE]
   Float_t         rcNSigmaK[10];   //[nRcE]
   Float_t         rcNSigmaP[10];   //[nRcE]
   Float_t         rcDca[10];   //[nRcE]

   // List of branches
   TBranch        *b_runId;   //!
   TBranch        *b_triggerId;   //!
   TBranch        *b_muEvtId;   //!
   TBranch        *b_muPriVertexX;   //!
   TBranch        *b_muPriVertexY;   //!
   TBranch        *b_muPriVertexZ;   //!
   TBranch        *b_muVpdVz;   //!
   TBranch        *b_muRefMult;   //!
   TBranch        *b_mugRefMult;   //!
   TBranch        *b_munBtofMatch;   //!
   TBranch        *b_zdcX;   //!
   TBranch        *b_Centrality_16;   //!
   TBranch        *b_Centrality_9;   //!
   TBranch        *b_refmult_corr;   //!
   TBranch        *b_reweight;   //!
   TBranch        *b_nMcTrks;   //!
   TBranch        *b_nMcVertices;   //!
   TBranch        *b_nMcE;   //!
   TBranch        *b_nRcE;   //!
   TBranch        *b_nMcTrack;   //!
   TBranch        *b_geantId;   //!
   TBranch        *b_mcTrackId;   //!
   TBranch        *b_mcCharge;   //!
   TBranch        *b_mcTrack_Vr;   //!
   TBranch        *b_mcTrack_par_Geantid;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_rcidtruth;   //!
   TBranch        *b_rcQaTruth;   //!
   TBranch        *b_rcflag;   //!
   TBranch        *b_rcCharge;   //!
   TBranch        *b_rcPt;   //!
   TBranch        *b_rcEta;   //!
   TBranch        *b_rcPhi;   //!
   TBranch        *b_rcNHitsFit;   //!
   TBranch        *b_rcNHitsPoss;   //!
   TBranch        *b_rcNHitsDedx;   //!
   TBranch        *b_rcDedx;   //!
   TBranch        *b_rcNSigmaE;   //!
   TBranch        *b_rcNSigmaPi;   //!
   TBranch        *b_rcNSigmaK;   //!
   TBranch        *b_rcNSigmaP;   //!
   TBranch        *b_rcDca;   //!

   MCEVENT(TTree *tree=0);
   virtual ~MCEVENT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MCEVENT_cxx
MCEVENT::MCEVENT(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("embproduction_7p7GeV_2021xxElectron_100_20224105xP22ib.SL22bx2021x121x22121011xst_physics_adc_22121011_raw_7500003.myminimc.root_22121011_raw_7500003.myminimc.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("embproduction_7p7GeV_2021xxElectron_100_20224105xP22ib.SL22bx2021x121x22121011xst_physics_adc_22121011_raw_7500003.myminimc.root_22121011_raw_7500003.myminimc.root");
      }
      f->GetObject("mcT",tree);

   }
   Init(tree);
}

MCEVENT::~MCEVENT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MCEVENT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MCEVENT::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MCEVENT::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runId", &runId, &b_runId);
   fChain->SetBranchAddress("triggerId", triggerId, &b_triggerId);
   fChain->SetBranchAddress("muEvtId", &muEvtId, &b_muEvtId);
   fChain->SetBranchAddress("muPriVertexX", &muPriVertexX, &b_muPriVertexX);
   fChain->SetBranchAddress("muPriVertexY", &muPriVertexY, &b_muPriVertexY);
   fChain->SetBranchAddress("muPriVertexZ", &muPriVertexZ, &b_muPriVertexZ);
   fChain->SetBranchAddress("muVpdVz", &muVpdVz, &b_muVpdVz);
   fChain->SetBranchAddress("muRefMult", &muRefMult, &b_muRefMult);
   fChain->SetBranchAddress("mugRefMult", &mugRefMult, &b_mugRefMult);
   fChain->SetBranchAddress("munBtofMatch", &munBtofMatch, &b_munBtofMatch);
   fChain->SetBranchAddress("zdcX", &zdcX, &b_zdcX);
   fChain->SetBranchAddress("Centrality_16", &Centrality_16, &b_Centrality_16);
   fChain->SetBranchAddress("Centrality_9", &Centrality_9, &b_Centrality_9);
   fChain->SetBranchAddress("refmult_corr", &refmult_corr, &b_refmult_corr);
   fChain->SetBranchAddress("reweight", &reweight, &b_reweight);
   fChain->SetBranchAddress("nMcTrks", &nMcTrks, &b_nMcTrks);
   fChain->SetBranchAddress("nMcVertices", &nMcVertices, &b_nMcVertices);
   fChain->SetBranchAddress("nMcE", &nMcE, &b_nMcE);
   fChain->SetBranchAddress("nRcE", &nRcE, &b_nRcE);
   fChain->SetBranchAddress("nMcTrack", &nMcTrack, &b_nMcTrack);
   fChain->SetBranchAddress("geantId", geantId, &b_geantId);
   fChain->SetBranchAddress("mcTrackId", mcTrackId, &b_mcTrackId);
   fChain->SetBranchAddress("mcCharge", mcCharge, &b_mcCharge);
   fChain->SetBranchAddress("mcTrack_Vr", mcTrack_Vr, &b_mcTrack_Vr);
   fChain->SetBranchAddress("mcTrack_par_Geantid", mcTrack_par_Geantid, &b_mcTrack_par_Geantid);
   fChain->SetBranchAddress("mcPt", mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcEta", mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("rcidtruth", rcidtruth, &b_rcidtruth);
   fChain->SetBranchAddress("rcQaTruth", rcQaTruth, &b_rcQaTruth);
   fChain->SetBranchAddress("rcflag", rcflag, &b_rcflag);
   fChain->SetBranchAddress("rcCharge", rcCharge, &b_rcCharge);
   fChain->SetBranchAddress("rcPt", rcPt, &b_rcPt);
   fChain->SetBranchAddress("rcEta", rcEta, &b_rcEta);
   fChain->SetBranchAddress("rcPhi", rcPhi, &b_rcPhi);
   fChain->SetBranchAddress("rcNHitsFit", rcNHitsFit, &b_rcNHitsFit);
   fChain->SetBranchAddress("rcNHitsPoss", rcNHitsPoss, &b_rcNHitsPoss);
   fChain->SetBranchAddress("rcNHitsDedx", rcNHitsDedx, &b_rcNHitsDedx);
   fChain->SetBranchAddress("rcDedx", rcDedx, &b_rcDedx);
   fChain->SetBranchAddress("rcNSigmaE", rcNSigmaE, &b_rcNSigmaE);
   fChain->SetBranchAddress("rcNSigmaPi", rcNSigmaPi, &b_rcNSigmaPi);
   fChain->SetBranchAddress("rcNSigmaK", rcNSigmaK, &b_rcNSigmaK);
   fChain->SetBranchAddress("rcNSigmaP", rcNSigmaP, &b_rcNSigmaP);
   fChain->SetBranchAddress("rcDca", rcDca, &b_rcDca);
   Notify();
}

Bool_t MCEVENT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MCEVENT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MCEVENT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MCEVENT_cxx
