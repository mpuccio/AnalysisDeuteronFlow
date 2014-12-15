//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 31 10:18:44 2014 by ROOT version 5.34/18
// from TTree deuterons/deuteron candidates
// found on file: FileMerger0.root
//////////////////////////////////////////////////////////

#ifndef AODSelector_h
#define AODSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TF1.h>

class TH1F;
class TH2F;
class TF1;


#define EPSILON 1E-16
#define MD 1.875612859f
#define M2D MD*MD
#define kNCent 5
#define kNBins 28
#define kNBinsTPC 3
#define kNDCABins 9
#define kNBinsTOF (kNBins - kNBinsTPC)
#define kNDCAbinsTOF (kNDCABins - kNBinsTPC)

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AODSelector : public TSelector {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types
  Float_t         centrality;
  Float_t         eta;
  Float_t         TPCnClust;
  Float_t         TPCsignal;
  Float_t         TPCnSignal;
  Float_t         TPCchi2;
  Float_t         ITSsignal;
  Float_t         ITSnClust;
  Float_t         ITSnClustPID;
  Float_t         TOFtime;
  Float_t         DCAxy;
  Float_t         DCAz;
  Float_t         p;
  Float_t         pTPC;
  Float_t         pT;
  Float_t         length;
  
  // List of branches
  TBranch        *b_centrality;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_TPCnClust;   //!
  TBranch        *b_TPCsignal;   //!
  TBranch        *b_TPCnSignal;   //!
  TBranch        *b_TPCchi2;   //!
  TBranch        *b_ITSsignal;   //!
  TBranch        *b_ITSnClust;   //!
  TBranch        *b_ITSnClustPID;   //!
  TBranch        *b_TOFtime;   //!
  TBranch        *b_DCAxy;   //!
  TBranch        *b_DCAz;   //!
  TBranch        *b_p;   //!
  TBranch        *b_pTPC;   //!
  TBranch        *b_pT;   //!
  TBranch        *b_length;   //!
  
  
  
  AODSelector(TTree * /*tree*/ =0) : fChain(0),
  fCorrectionAD("fCorrectionAD","[0]+[1]*exp([2]*x)",0,10),
  fCorrectionD("fCorrectionD","[0]+[1]*exp([2]*x)",0,10) {
    float bins[kNBins + 1] = {
      0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
      1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
      3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f,8.0f,10.f
    };
    for (int i = 0; i < kNBins + 1; ++i) fBins[i] = bins[i];
    
    double cent[kNCent + 1] = {0.,10.,20.,40.,60.,80.};
    for (int i = 0; i < kNCent + 1; ++i) fCentralityClasses[i] = cent[i];
    
    fCorrectionAD.SetParameters(-2.10154e-03,-4.53472e-01,-3.01246e+00);
    fCorrectionD.SetParameters(-2.00277e-03,-4.93461e-01,-3.05463e+00);
  }
  virtual ~AODSelector() { }
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
  Int_t            GetCentBin(float cent);
  Int_t            GetPtBin(float pt);
  TH1F*            fBeta;
  TH2F*            fBeta2D;
  TH2F*            fBeta2DPt;
  Double_t         fBins[kNBins + 1];      // length = fkNBins + 1
  Double_t         fCentralityClasses[kNCent + 1];
  TH1F*            fDCAxy;
  TH1F*            fDCAz;
  TH2F*            fDCA2D;
  TH1F*            fSignalAD[kNBinsTOF * kNCent]; // length = nbinstof x 5
  TH1F*            fSignalD[kNBinsTOF * kNCent];  // length = nbinstof x 5
  TF1*             fDeutBB;
  TH1F*            fTPCSignalN;
  TH2F*            fdEdxTPC;
  TH2F*            fdEdxTPCpT;
  TH2F*            fdEdxTPCproj;
  TH2F*            fdEdxTPCSignal;
  TH1F*            fdEdxTPCSignalCounts[kNCent];
  TH1F*            fdEdxTPCSignalCountsAD[kNCent];
  TH1F*            fCentrality;
  TH1F*            fCentralityClass;
  TH1F*            fGamma;
  
  TH2F            *fDdcaXY[kNCent];
  TH2F            *fDdcaZ[kNCent];
  TH2F            *fDCASignal[kNDCAbinsTOF * kNCent];  // length = fkNDCAbinsTOF x 5
  
  TF1              fCorrectionAD;
  TF1              fCorrectionD;
  
  Float_t          fPrevious;
  ClassDef(AODSelector,0);
};

#endif

#ifdef AODSelector_cxx
void AODSelector::Init(TTree *tree)
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
  fChain->SetBranchAddress("ITSsignal", &ITSsignal, &b_ITSsignal);
  fChain->SetBranchAddress("ITSnClust", &ITSnClust, &b_ITSnClust);
  fChain->SetBranchAddress("ITSnClustPID", &ITSnClustPID, &b_ITSnClustPID);
  fChain->SetBranchAddress("TOFtime", &TOFtime, &b_TOFtime);
  fChain->SetBranchAddress("DCAxy", &DCAxy, &b_DCAxy);
  fChain->SetBranchAddress("DCAz", &DCAz, &b_DCAz);
  fChain->SetBranchAddress("p", &p, &b_p);
  fChain->SetBranchAddress("pTPC", &pTPC, &b_pTPC);
  fChain->SetBranchAddress("pT", &pT, &b_pT);
  fChain->SetBranchAddress("length", &length, &b_length);
}

Bool_t AODSelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

#endif // #ifdef AODSelector_cxx
