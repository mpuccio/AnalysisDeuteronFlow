//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 26 09:36:41 2014 by ROOT version 5.34/22
// from TTree Deuterons/Deuterons MC
// found on file: FileMerger.root
//////////////////////////////////////////////////////////

#ifndef EfficiencySelector_h
#define EfficiencySelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class TH1F;
class TH2F;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class EfficiencySelector : public TSelector {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types
  Float_t         pMC;
  Float_t         pTMC;
  Float_t         etaMC;
  Float_t         phiMC;
  Float_t         yMC;
  Float_t         p;
  Float_t         pT;
  Float_t         eta;
  Float_t         phi;
  Float_t         beta;
  Float_t         DCAxy;
  Float_t         DCAz;
  Float_t         chi2;
  Float_t         centrality;
  UShort_t        TPCnClusters;
  UShort_t        TPCnSignal;
  Char_t          ITSnClusters;
  Char_t          ITSnSignal;
  Bool_t          IsPrimary;
  Bool_t          IsSecondaryFromMaterial;
  
  // List of branches
  TBranch        *b_pMC;   //!
  TBranch        *b_pTMC;   //!
  TBranch        *b_etaMC;   //!
  TBranch        *b_phiMC;   //!
  TBranch        *b_yMC;   //!
  TBranch        *b_p;   //!
  TBranch        *b_pT;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_beta;   //!
  TBranch        *b_DCAxy;   //!
  TBranch        *b_DCAz;   //!
  TBranch        *b_chi2;   //!
  TBranch        *b_centrality;   //
  TBranch        *b_TPCnClusters;   //!
  TBranch        *b_TPCnSignal;   //!
  TBranch        *b_ITSnClusters;   //!
  TBranch        *b_ITSnSignal;   //!
  TBranch        *b_IsPrimary;   //!
  TBranch        *b_IsSecondaryFromMaterial;   //!
  
  EfficiencySelector(TTree * /*tree*/ =0)
  : fChain(0)
  , fAntiDMCYield(0x0)
  , fAntiDYield(0x0)
  , fAntiDYieldTOF(0x0)
  , fDMCYield(0x0)
  , fDYield(0x0)
  , fDYieldTOF(0x0)
  , fDdcaXYprimaries()
  , fDdcaZprimaries()
  , fDdcaXYsecondaries()
  , fDdcaZsecondaries() { }
  
  virtual ~EfficiencySelector() { }
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
  TH1F    *fAntiDMCYield;
  TH1F    *fAntiDYield;
  TH1F    *fAntiDYieldTOF;
  TH1F    *fDMCYield;
  TH1F    *fDYield;
  TH1F    *fDYieldTOF;
  
  TH2F    *fDdcaXYprimaries[5];
  TH2F    *fDdcaZprimaries[5];
  TH2F    *fDdcaXYsecondaries[5];
  TH2F    *fDdcaZsecondaries[5];

  ClassDef(EfficiencySelector,0);
};

#endif

#ifdef EfficiencySelector_cxx
void EfficiencySelector::Init(TTree *tree)
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
  
  fChain->SetBranchAddress("pMC", &pMC, &b_pMC);
  fChain->SetBranchAddress("pTMC", &pTMC, &b_pTMC);
  fChain->SetBranchAddress("etaMC", &etaMC, &b_etaMC);
  fChain->SetBranchAddress("phiMC", &phiMC, &b_phiMC);
  fChain->SetBranchAddress("yMC", &yMC, &b_yMC);
  fChain->SetBranchAddress("p", &p, &b_p);
  fChain->SetBranchAddress("pT", &pT, &b_pT);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("phi", &phi, &b_phi);
  fChain->SetBranchAddress("beta", &beta, &b_beta);
  fChain->SetBranchAddress("DCAxy", &DCAxy, &b_DCAxy);
  fChain->SetBranchAddress("DCAz", &DCAz, &b_DCAz);
  fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
  fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
  fChain->SetBranchAddress("TPCnClusters", &TPCnClusters, &b_TPCnClusters);
  fChain->SetBranchAddress("TPCnSignal", &TPCnSignal, &b_TPCnSignal);
  fChain->SetBranchAddress("ITSnClusters", &ITSnClusters, &b_ITSnClusters);
  fChain->SetBranchAddress("ITSnSignal", &ITSnSignal, &b_ITSnSignal);
  fChain->SetBranchAddress("IsPrimary", &IsPrimary, &b_IsPrimary);
  fChain->SetBranchAddress("IsSecondaryFromMaterial", &IsSecondaryFromMaterial, &b_IsSecondaryFromMaterial);
}

Bool_t EfficiencySelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

#endif // #ifdef EfficiencySelector_cxx
