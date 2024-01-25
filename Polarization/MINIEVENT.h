//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 15 15:42:24 2017 by ROOT version 5.34/30
// from TTree qaTree/data for QA
// found on file: BBE00A021D5852556F14571E7EBCCF4A_9.root
//////////////////////////////////////////////////////////

#ifndef MINIEVENT_h
#define MINIEVENT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MINIEVENT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runId;
   Int_t           eventId;
   Float_t         bField;
   Double_t        zdcRate;
   Double_t        bbcRate;
   Float_t         tpcVx;
   Float_t         tpcVy;
   Float_t         tpcVz;
   Float_t         vpdVz;
   Int_t           npTracks;
   Int_t           nHitsFit[50];   //[npTracks]
   Int_t           nHitsPoss[50];   //[npTracks]
   Int_t           nHitsDedx[50];   //[npTracks]
   Float_t         nSigmaPion[50];   //[npTracks]
   Float_t         pt[50];   //[npTracks]
   Float_t         eta[50];   //[npTracks]
   Float_t         phi[50];   //[npTracks]
   Float_t         dca[50];   //[npTracks]
   Short_t         charge[50];   //[npTracks]
   Int_t           backleg[50];   //[npTracks]
   Int_t           module[50];   //[npTracks]
   Int_t           matchFlag[50];   //[npTracks]
   Float_t         dy[50];   //[npTracks]
   Float_t         dz[50];   //[npTracks]
   Float_t         dtof[50];   //[npTracks]

   // List of branches
   TBranch        *b_runId;   //!
   TBranch        *b_eventId;   //!
   TBranch        *b_bField;   //!
   TBranch        *b_zdcRate;   //!
   TBranch        *b_bbcRate;   //!
   TBranch        *b_tpcVx;   //!
   TBranch        *b_tpcVy;   //!
   TBranch        *b_tpcVz;   //!
   TBranch        *b_vpdVz;   //!
   TBranch        *b_npTracks;   //!
   TBranch        *b_nHitsFit;   //!
   TBranch        *b_nHitsPoss;   //!
   TBranch        *b_nHitsDedx;   //!
   TBranch        *b_nSigmaPion;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_dca;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_backleg;   //!
   TBranch        *b_module;   //!
   TBranch        *b_matchFlag;   //!
   TBranch        *b_dy;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_dtof;   //!

   MINIEVENT(TTree *tree=0);
   virtual ~MINIEVENT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MINIEVENT_cxx
MINIEVENT::MINIEVENT(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("BBE00A021D5852556F14571E7EBCCF4A_9.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("BBE00A021D5852556F14571E7EBCCF4A_9.root");
      }
      f->GetObject("qaTree",tree);

   }
   Init(tree);
}

MINIEVENT::~MINIEVENT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MINIEVENT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MINIEVENT::LoadTree(Long64_t entry)
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

void MINIEVENT::Init(TTree *tree)
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
   fChain->SetBranchAddress("eventId", &eventId, &b_eventId);
   fChain->SetBranchAddress("bField", &bField, &b_bField);
   fChain->SetBranchAddress("zdcRate", &zdcRate, &b_zdcRate);
   fChain->SetBranchAddress("bbcRate", &bbcRate, &b_bbcRate);
   fChain->SetBranchAddress("tpcVx", &tpcVx, &b_tpcVx);
   fChain->SetBranchAddress("tpcVy", &tpcVy, &b_tpcVy);
   fChain->SetBranchAddress("tpcVz", &tpcVz, &b_tpcVz);
   fChain->SetBranchAddress("vpdVz", &vpdVz, &b_vpdVz);
   fChain->SetBranchAddress("npTracks", &npTracks, &b_npTracks);
   fChain->SetBranchAddress("nHitsFit", nHitsFit, &b_nHitsFit);
   fChain->SetBranchAddress("nHitsPoss", nHitsPoss, &b_nHitsPoss);
   fChain->SetBranchAddress("nHitsDedx", nHitsDedx, &b_nHitsDedx);
   fChain->SetBranchAddress("nSigmaPion", nSigmaPion, &b_nSigmaPion);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("dca", dca, &b_dca);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("backleg", backleg, &b_backleg);
   fChain->SetBranchAddress("module", module, &b_module);
   fChain->SetBranchAddress("matchFlag", matchFlag, &b_matchFlag);
   fChain->SetBranchAddress("dy", dy, &b_dy);
   fChain->SetBranchAddress("dz", dz, &b_dz);
   fChain->SetBranchAddress("dtof", dtof, &b_dtof);
   Notify();
}

Bool_t MINIEVENT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MINIEVENT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MINIEVENT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MINIEVENT_cxx
