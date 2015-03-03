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
#include <TString.h>
#include <TArrayD.h>

class TH1F;
class TH2F;
class TH3F;


#define EPSILON 1E-16
#define MD 1.875612859f
#define M2D MD*MD




// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AODSelector : public TSelector {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  void SetPtBins(Int_t n, Double_t *arr) { fBins.Set(n, arr); }
  void SetCentBins(Int_t n, Double_t *arr) { fCentralityBins.Set(n, arr); }

  void SetRequireITSrecPoints (int rec = 4) { fRequireITSrecPoints = rec; }
  void SetRequireITSsignal (int sig = 3) { fRequireITSsignal = sig; }
  void SetRequireTPCrecPoints (int rec = 70) { fRequireITSrecPoints = rec; }
  void SetRequireTPCsignal (int sig = 70) { fRequireITSsignal = sig; }
  void SetRequireSPDrecPoints (int rec = 1) { fRequireSPDrecPoints = rec; }
  void SetEtaRange (float emin, float emax) { fRequireEtaMin = emin; fRequireEtaMax = emax; }
  void SetYRange (float ymin, float ymax) { fRequireYmin = ymin; fRequireYmax = ymax; }
  void SetRequireMaxChi2 (float maxChi2 = 4.f) { fRequireMaxChi2 = maxChi2; }
  void SetRequireMaxDCAxy (float maxDCA) { fRequireMaxDCAxy = maxDCA; }
  void SetRequireMaxDCAz (float maxDCA) { fRequireMaxDCAz = maxDCA; }
  void SetRequireTPCpidSigmas (float sig) { fRequireTPCpidSigmas = (sig > 0) ? sig : -sig; }
  void SetRequireITSpidSigmas (float sig) { fRequireITSpidSigmas = sig; }
  
  void SetOutputOption(TString name, Bool_t rec = kFALSE) { fTaskName = name; fRecreate = rec; }
  // Declaration of leaf types
  Float_t         centrality;
  Float_t         eta;
  Float_t         TPCsignal;
  Float_t         chi2NDF;
  Float_t         TOFtime;
  Float_t         DCAxy;
  Float_t         DCAz;
  Float_t         p;
  Float_t         pTPC;
  Float_t         pT;
  Float_t         length;
  Float_t         ITSsigmad;
  Float_t         ITSsigmat;
  Float_t         TPCsigmad;
  Float_t         TPCsigmat;
  Int_t           FilterMap;
  UShort_t        TPCnClust;
  UShort_t        TPCnClustShared;
  UShort_t        TPCnSignal;
  UShort_t        trigger;
  UChar_t         ITSnClust;
  UChar_t         ITSnSignal;
  
  // List of branches
  TBranch        *b_centrality;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_TPCsignal;   //!
  TBranch        *b_chi2NDF;   //!
  TBranch        *b_TOFtime;   //!
  TBranch        *b_DCAxy;   //!
  TBranch        *b_DCAz;   //!
  TBranch        *b_p;   //!
  TBranch        *b_pTPC;   //!
  TBranch        *b_pT;   //!
  TBranch        *b_length;   //!
  TBranch        *b_ITSsigmad;   //!
  TBranch        *b_ITSsigmat;   //!
  TBranch        *b_TPCsigmad;   //!
  TBranch        *b_TPCsigmat;   //!
  TBranch        *b_FilterMap;   //!
  TBranch        *b_TPCnClust;   //!
  TBranch        *b_TPCnClustShared;   //!
  TBranch        *b_TPCnSignal;   //!
  TBranch        *b_trigger;   //!
  TBranch        *b_ITSnClust;   //!
  TBranch        *b_ITSnSignal;   //!
  
  
  AODSelector(TTree * /*tree*/ =0) : fChain(0)
  ,fTaskName("standardname")
  ,fRecreate(kTRUE)
  ,fCorrectionAD("fCorrectionAD","[0]+[1]*exp([2]*x)",0,10)
  ,fCorrectionD("fCorrectionD","[0]+[1]*exp([2]*x)",0,10)
  ,fRequireITSrecPoints(2u)
  ,fRequireITSsignal(0u)
  ,fRequireSPDrecPoints(1u)
  ,fRequireTPCrecPoints(70u)
  ,fRequireTPCsignal(70u)
  ,fRequireEtaMin(-0.8f)
  ,fRequireEtaMax(0.8f)
  ,fRequireYmin(-0.5f)
  ,fRequireYmax(0.5f)
  ,fRequireMaxChi2(4.f)
  ,fRequireMaxDCAxy(0.5f)
  ,fRequireMaxDCAz(1.f)
  ,fRequireTPCpidSigmas(3.f)
  ,fRequireITSpidSigmas(-1.f)
  {
    for (int i = 0; i < kNBins + 1; ++i) fBins[i] = kBins[i];
    for (int i = 0; i < kNCent + 1; ++i) fCentralityBins[i] = kCent[i];
    
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
  Bool_t           Flatten(float cent);
  
  TString          fTaskName;
  Bool_t           fRecreate;
  TArrayD          fBins;
  TArrayD          fCentralityBins;
  TF1*             fDeutBB;
  TH2F*            fTPCSignal;
  
  Bool_t           fRequireITSrefit;       ///< Cut on tracks: set true to require ITS refit
  Bool_t           fRequireTPCrefit;       ///< Cut on tracks: set true to require TPC refit
  Bool_t           fRequireNoKinks;        ///< Cut on tracks: set true to exclude tracks from kink vertices
  UShort_t         fRequireITSrecPoints;   ///< Cut on tracks: minimum number of required ITS recpoints
  UShort_t         fRequireITSsignal;      ///< Cut on tracks: minimum number of required ITS PID recpoints
  UShort_t         fRequireSPDrecPoints;   ///< Cut on tracks: minimum number of required SPD recpoints
  UShort_t         fRequireTPCrecPoints;   ///< Cut on tracks: minimum number of required TPC recpoints
  UShort_t         fRequireTPCsignal;      ///< Cut on tracks: minimum number of required TPC PID recpoints
  Float_t          fRequireEtaMin;         ///< Cut on tracks: minimum eta for the track
  Float_t          fRequireEtaMax;         ///< Cut on tracks: maximum eta for the track
  Float_t          fRequireYmin;           ///< Cut on tracks: mimimum y for the track (using PDG mass)
  Float_t          fRequireYmax;           ///< Cut on tracks: maximum y for the track (using PDG mass)
  Float_t          fRequireMaxChi2;        ///< Cut on tracks: maximum TPC \f$\chi^{2}/NDF\f$
  Float_t          fRequireMaxDCAxy;       ///< Cut on tracks: maximum \f$DCA_{xy}\f$ for the track
  Float_t          fRequireMaxDCAz;        ///< Cut on tracks: maximum \f$DCA_{z}\f$ for the track
  Float_t          fRequireTPCpidSigmas;   ///< Cut on TPC PID number of sigmas
  Float_t          fRequireITSpidSigmas;

  
  TH1F            *fCentrality;            //!< Events centrality distribution
  TH1F            *fFlattenedCentrality;   //!< Events centrality distribution after the flattening
  TH1F            *fCentralityClasses;     //!< Events statistics per centrality classes
  
  TH3F            *fATOFsignal;            //!<
  TH3F            *fATPCcounts;            //!<
  TH3F            *fMDCAxy;                //!<
  TH3F            *fMDCAz;                 //!<
  TH3F            *fMTOFsignal;            //!<
  TH3F            *fMTPCcounts;            //!<
  
  TF1              fCorrectionAD;
  TF1              fCorrectionD;
  
  Bool_t           fSkipEvent;
  
  enum Trigger_m {
    kMB = 1,
    kCentral = 2,
    kSemiCentral = 4,
    kEMCEJE = 8,
    kEMCEGA = 16
  };

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
  fChain->SetBranchAddress("TPCsignal", &TPCsignal, &b_TPCsignal);
  fChain->SetBranchAddress("chi2NDF", &chi2NDF, &b_chi2NDF);
  fChain->SetBranchAddress("TOFtime", &TOFtime, &b_TOFtime);
  fChain->SetBranchAddress("DCAxy", &DCAxy, &b_DCAxy);
  fChain->SetBranchAddress("DCAz", &DCAz, &b_DCAz);
  fChain->SetBranchAddress("p", &p, &b_p);
  fChain->SetBranchAddress("pTPC", &pTPC, &b_pTPC);
  fChain->SetBranchAddress("pT", &pT, &b_pT);
  fChain->SetBranchAddress("length", &length, &b_length);
  fChain->SetBranchAddress("ITSsigmad", &ITSsigmad, &b_ITSsigmad);
  fChain->SetBranchAddress("ITSsigmat", &ITSsigmat, &b_ITSsigmat);
  fChain->SetBranchAddress("TPCsigmad", &TPCsigmad, &b_TPCsigmad);
  fChain->SetBranchAddress("TPCsigmat", &TPCsigmat, &b_TPCsigmat);
  fChain->SetBranchAddress("FilterMap", &FilterMap, &b_FilterMap);
  fChain->SetBranchAddress("TPCnClust", &TPCnClust, &b_TPCnClust);
  fChain->SetBranchAddress("TPCnClustShared", &TPCnClustShared, &b_TPCnClustShared);
  fChain->SetBranchAddress("TPCnSignal", &TPCnSignal, &b_TPCnSignal);
  fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
  fChain->SetBranchAddress("ITSnClust", &ITSnClust, &b_ITSnClust);
  fChain->SetBranchAddress("ITSnSignal", &ITSnSignal, &b_ITSnSignal);
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
