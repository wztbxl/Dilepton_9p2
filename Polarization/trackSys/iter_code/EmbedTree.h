//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May  1 05:07:22 2017 by ROOT version 5.34/30
// from TTree EmbedTree/Embedding tree
// found on file: ../Embedding15/submitDir/rootfile_rank/Jpsi_1D412F9AB2E90093B3BB5BA2413B25DF_9.root
//////////////////////////////////////////////////////////

#ifndef EmbedTree_h
#define EmbedTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class EmbedTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runID;
   Int_t           tofMult;
   Float_t         mBField;
   Float_t         bbcRate;
   Float_t         zdcRate;
   Float_t         tpcVx;
   Float_t         tpcVy;
   Float_t         tpcVz;
   Float_t         vpdVz;
   Int_t           rank;
   Int_t           nMcTracks;
   Int_t           nMcMuon;
   Int_t           nMcJpsi;
   Int_t           nRcJpsi;
   Int_t           mcPkey[500];   //[nMcTracks]
   Int_t           mcGeantId[500];   //[nMcTracks]
   Double_t        mcpt[500];   //[nMcTracks]
   Double_t        mcphi[500];   //[nMcTracks]
   Double_t        mceta[500];   //[nMcTracks]
   Int_t           mccharge[500];   //[nMcTracks]
   Int_t           rcPkey[500];   //[nMcTracks]
   Int_t           rcNHitsFit[500];   //[nMcTracks]
   Int_t           rcNHitsPoss[500];   //[nMcTracks]
   Int_t           rcNHitsDedx[500];   //[nMcTracks]
   Double_t        rcDca[500];   //[nMcTracks]
   Double_t        rcpt[500];   //[nMcTracks]
   Double_t        rcphi[500];   //[nMcTracks]
   Double_t        rceta[500];   //[nMcTracks]
   Double_t        rcNSigmaPi[500];   //[nMcTracks]
   Int_t           rcCharge[500];   //[nMcTracks]
   Int_t           rcBackleg[500];   //[nMcTracks]
   Int_t           rcModule[500];   //[nMcTracks]
   Double_t        rcDz[500];   //[nMcTracks]
   Double_t        rcDy[500];   //[nMcTracks]
   Double_t        rcDtof[500];   //[nMcTracks]
   Int_t           passTrkCut[500];   //[nMcTracks]
   Int_t           passMuonCut[500];   //[nMcTracks]

   // List of branches
   TBranch        *b_runID;   //!
   TBranch        *b_tofMult;   //!
   TBranch        *b_mBField;   //!
   TBranch        *b_bbcRate;   //!
   TBranch        *b_zdcRate;   //!
   TBranch        *b_tpcVx;   //!
   TBranch        *b_tpcVy;   //!
   TBranch        *b_tpcVz;   //!
   TBranch        *b_vpdVz;   //!
   TBranch        *b_rank;   //!
   TBranch        *b_nMcTracks;   //!
   TBranch        *b_nMcMuon;   //!
   TBranch        *b_nMcJpsi;   //!
   TBranch        *b_nRcJpsi;   //!
   TBranch        *b_mcPkey;   //!
   TBranch        *b_mcGeantId;   //!
   TBranch        *b_mcpt;   //!
   TBranch        *b_mcphi;   //!
   TBranch        *b_mceta;   //!
   TBranch        *b_mccharge;   //!
   TBranch        *b_rcPkey;   //!
   TBranch        *b_rcNHitsFit;   //!
   TBranch        *b_rcNHitsPoss;   //!
   TBranch        *b_rcNHitsDedx;   //!
   TBranch        *b_rcDca;   //!
   TBranch        *b_rcpt;   //!
   TBranch        *b_rcphi;   //!
   TBranch        *b_rceta;   //!
   TBranch        *b_rcNSigmaPi;   //!
   TBranch        *b_rcCharge;   //!
   TBranch        *b_rcBackleg;   //!
   TBranch        *b_rcModule;   //!
   TBranch        *b_rcDz;   //!
   TBranch        *b_rcDy;   //!
   TBranch        *b_rcDtof;   //!
   TBranch        *b_passTrkCut;   //!
   TBranch        *b_passMuonCut;   //!

   EmbedTree(TTree *tree=0);
   virtual ~EmbedTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EmbedTree_cxx
EmbedTree::EmbedTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Embedding15/submitDir/rootfile_rank/Jpsi_1D412F9AB2E90093B3BB5BA2413B25DF_9.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../Embedding15/submitDir/rootfile_rank/Jpsi_1D412F9AB2E90093B3BB5BA2413B25DF_9.root");
      }
      f->GetObject("EmbedTree",tree);

   }
   Init(tree);
}

EmbedTree::~EmbedTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EmbedTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EmbedTree::LoadTree(Long64_t entry)
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

void EmbedTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("tofMult", &tofMult, &b_tofMult);
   fChain->SetBranchAddress("mBField", &mBField, &b_mBField);
   fChain->SetBranchAddress("bbcRate", &bbcRate, &b_bbcRate);
   fChain->SetBranchAddress("zdcRate", &zdcRate, &b_zdcRate);
   fChain->SetBranchAddress("tpcVx", &tpcVx, &b_tpcVx);
   fChain->SetBranchAddress("tpcVy", &tpcVy, &b_tpcVy);
   fChain->SetBranchAddress("tpcVz", &tpcVz, &b_tpcVz);
   fChain->SetBranchAddress("vpdVz", &vpdVz, &b_vpdVz);
   fChain->SetBranchAddress("rank", &rank, &b_rank);
   fChain->SetBranchAddress("nMcTracks", &nMcTracks, &b_nMcTracks);
   fChain->SetBranchAddress("nMcMuon", &nMcMuon, &b_nMcMuon);
   fChain->SetBranchAddress("nMcJpsi", &nMcJpsi, &b_nMcJpsi);
   fChain->SetBranchAddress("nRcJpsi", &nRcJpsi, &b_nRcJpsi);
   fChain->SetBranchAddress("mcPkey", mcPkey, &b_mcPkey);
   fChain->SetBranchAddress("mcGeantId", mcGeantId, &b_mcGeantId);
   fChain->SetBranchAddress("mcpt", mcpt, &b_mcpt);
   fChain->SetBranchAddress("mcphi", mcphi, &b_mcphi);
   fChain->SetBranchAddress("mceta", mceta, &b_mceta);
   fChain->SetBranchAddress("mccharge", mccharge, &b_mccharge);
   fChain->SetBranchAddress("rcPkey", rcPkey, &b_rcPkey);
   fChain->SetBranchAddress("rcNHitsFit", rcNHitsFit, &b_rcNHitsFit);
   fChain->SetBranchAddress("rcNHitsPoss", rcNHitsPoss, &b_rcNHitsPoss);
   fChain->SetBranchAddress("rcNHitsDedx", rcNHitsDedx, &b_rcNHitsDedx);
   fChain->SetBranchAddress("rcDca", rcDca, &b_rcDca);
   fChain->SetBranchAddress("rcpt", rcpt, &b_rcpt);
   fChain->SetBranchAddress("rcphi", rcphi, &b_rcphi);
   fChain->SetBranchAddress("rceta", rceta, &b_rceta);
   fChain->SetBranchAddress("rcNSigmaPi", rcNSigmaPi, &b_rcNSigmaPi);
   fChain->SetBranchAddress("rcCharge", rcCharge, &b_rcCharge);
   fChain->SetBranchAddress("rcBackleg", rcBackleg, &b_rcBackleg);
   fChain->SetBranchAddress("rcModule", rcModule, &b_rcModule);
   fChain->SetBranchAddress("rcDz", rcDz, &b_rcDz);
   fChain->SetBranchAddress("rcDy", rcDy, &b_rcDy);
   fChain->SetBranchAddress("rcDtof", rcDtof, &b_rcDtof);
   fChain->SetBranchAddress("passTrkCut", passTrkCut, &b_passTrkCut);
   fChain->SetBranchAddress("passMuonCut", passMuonCut, &b_passMuonCut);
   Notify();
}

Bool_t EmbedTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EmbedTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EmbedTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EmbedTree_cxx
