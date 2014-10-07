//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  7 16:11:52 2014 by ROOT version 5.34/18
// from TTree deuterons/deuteron candidates
// found on file: mpuccio_Flowdnt.root
//////////////////////////////////////////////////////////

#ifndef DeutSelector_h
#define DeutSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>              

class TH2F;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class DeutSelector : public TSelector {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types
  Float_t         centrality;
  Float_t         eta;
  Float_t         TPCnClust;
  Float_t         TPCsignal;
  Float_t         TPCnSignal;
  Float_t         TPCchi2;
  Float_t         TPCshared;
  Float_t         ITSsignal;
  Float_t         ITSnClust;
  Float_t         ITSnClustPID;
  Float_t         ITSchi2;
  Float_t         TOFtime;
  Float_t         TOFsignalDz;
  Float_t         TOFsignalDx;
  Float_t         DCAxy;
  Float_t         DCAz;
  Float_t         p;
  Float_t         pTPC;
  Float_t         pT;
  Float_t         length;
  Float_t         sigmaQP;
  
  // List of branches
  TBranch        *b_centrality;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_TPCnClust;   //!
  TBranch        *b_TPCsignal;   //!
  TBranch        *b_TPCnSignal;   //!
  TBranch        *b_TPCchi2;   //!
  TBranch        *b_TPCshared;   //!
  TBranch        *b_ITSsignal;   //!
  TBranch        *b_ITSnClust;   //!
  TBranch        *b_ITSnClustPID;   //!
  TBranch        *b_ITSchi2;   //!
  TBranch        *b_TOFtime;   //!
  TBranch        *b_TOFsignalDz;   //!
  TBranch        *b_TOFsignalDx;   //!
  TBranch        *b_DCAxy;   //!
  TBranch        *b_DCAz;   //!
  TBranch        *b_p;   //!
  TBranch        *b_pTPC;   //!
  TBranch        *b_pT;   //!
  TBranch        *b_length;   //!
  TBranch        *b_sigmaQP;   //!
  
  DeutSelector(TTree * /*tree*/ =0) : fChain(0) { }
  virtual ~DeutSelector() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
private:
  TH2F*            fdEdxTPC;
  ClassDef(DeutSelector,0);
};

#endif

#ifdef DeutSelector_cxx
void DeutSelector::Init(TTree *tree)
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
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("TPCnClust", &TPCnClust, &b_TPCnClust);
  fChain->SetBranchAddress("TPCsignal", &TPCsignal, &b_TPCsignal);
  fChain->SetBranchAddress("TPCnSignal", &TPCnSignal, &b_TPCnSignal);
  fChain->SetBranchAddress("TPCchi2", &TPCchi2, &b_TPCchi2);
  fChain->SetBranchAddress("TPCshared", &TPCshared, &b_TPCshared);
  fChain->SetBranchAddress("ITSsignal", &ITSsignal, &b_ITSsignal);
  fChain->SetBranchAddress("ITSnClust", &ITSnClust, &b_ITSnClust);
  fChain->SetBranchAddress("ITSnClustPID", &ITSnClustPID, &b_ITSnClustPID);
  fChain->SetBranchAddress("ITSchi2", &ITSchi2, &b_ITSchi2);
  fChain->SetBranchAddress("TOFtime", &TOFtime, &b_TOFtime);
  fChain->SetBranchAddress("TOFsignalDz", &TOFsignalDz, &b_TOFsignalDz);
  fChain->SetBranchAddress("TOFsignalDx", &TOFsignalDx, &b_TOFsignalDx);
  fChain->SetBranchAddress("DCAxy", &DCAxy, &b_DCAxy);
  fChain->SetBranchAddress("DCAz", &DCAz, &b_DCAz);
  fChain->SetBranchAddress("p", &p, &b_p);
  fChain->SetBranchAddress("pTPC", &pTPC, &b_pTPC);
  fChain->SetBranchAddress("pT", &pT, &b_pT);
  fChain->SetBranchAddress("length", &length, &b_length);
  fChain->SetBranchAddress("sigmaQP", &sigmaQP, &b_sigmaQP);
}

Bool_t DeutSelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

#endif // #ifdef DeutSelector_cxx
