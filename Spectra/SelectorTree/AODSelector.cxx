#define AODSelector_cxx
// The class definition in AODSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("AODSelector.cxx")
// Root > T->Process("AODSelector.cxx","some options")
// Root > T->Process("AODSelector.cxx+")
//

#include "AODSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>

static void BinLogAxis(const TH1 *h)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = const_cast<TAxis*>(h->GetXaxis());
  const Int_t bins = axis->GetNbins();
  
  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
  
  newBins[0] = from;
  Double_t factor = pow(to / from, 1. / bins);
  
  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i - 1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
}

Double_t BetheBlochAleph(Double_t bg,
                         Double_t kp1,
                         Double_t kp2,
                         Double_t kp3,
                         Double_t kp4,
                         Double_t kp5) {
  Double_t beta = bg / TMath::Sqrt(1. + bg * bg);
  Double_t aa = TMath::Power(beta,kp4);
  Double_t bb = TMath::Power(1. / bg,kp5);
  bb = TMath::Log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa;
}

Double_t DeuteronTPC(Double_t *x, Double_t *) {
  // Deuteron expected signal in TPC, taken from AntiHe4 task
  return BetheBlochAleph(x[0] / MD ,4.69637,7.51827,0.0183746,2.60,2.7);
}

unsigned int Log2Int(const unsigned int v) {
  // Taken from https://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious
  static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000};
  unsigned int r = (v & b[0]) != 0;
  r |= ((v & b[4]) != 0) << 4;
  r |= ((v & b[3]) != 0) << 3;
  r |= ((v & b[2]) != 0) << 2;
  r |= ((v & b[1]) != 0) << 1;
  return r;
}

void AODSelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
}

Int_t AODSelector::GetCentBin(float cent) {
  if (cent < fCentralityBins[0]) return -1;

  for (int i = 0; i < kNCent; ++i)
    if (cent <= fCentralityBins[i + 1])
      return i;
  
  return -1;
}

Int_t AODSelector::GetPtBin(float pt) {
  for (int i = 0; i < kNBins; ++i) {
    if (pt < fBins[i+1] && pt >= fBins[i]) {
      return i;
    }
  }
  return -1;
}

Bool_t AODSelector::Flatten(float cent) {
  float prob[13] = {
    0.855566,0.846964,0.829618,0.829259,0.830984,
    0.85094,0.844346,0.851818,0.874758,1,
    0.374767,0.650491,0.946963
  };
  return gRandom->Rndm() > prob[int(cent)];
}

void AODSelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  const Int_t nPtBins = kNBins;
  const Int_t nCentBins = kNCent;
  const Double_t *pTbins = fBins;
  const Double_t *centBins = fCentralityBins;
  
  const Int_t nDCAbins = 80;
  const Float_t range = 0.5;
  Double_t dcaBins[nDCAbins+1];
  for (int i = 0; i <= nDCAbins; ++i)
    dcaBins[i] = i * range * 2 / nDCAbins - range;
  
  fCentrality = new TH1F("fCentrality",";Centrality (%);Events / 1%;",100,0.,100.);
  fCentralityClasses = new TH1F("fCentralityClasses",";Centrality classes(%);Events / Class;",nCentBins,centBins);
  fFlattenedCentrality = new TH1F("fFlattenCentrality","Centrality distribution after the flattening;Centrality (%);Events / 1%;",100,0.,100.);
  GetOutputList()->Add(fCentrality);
  GetOutputList()->Add(fCentralityClasses);
  GetOutputList()->Add(fFlattenedCentrality);
  
  const int tofNbins = 75;
  const float tofLowBoundary = -2.4,tofHighBoundary = 3.6;
  Double_t tofBins[tofNbins + 1];
  const float deltaTOF = (tofHighBoundary - tofLowBoundary) / tofNbins;
  for (int i = 0; i <= tofNbins; ++i)
    tofBins[i] = i * deltaTOF + tofLowBoundary;
  const int dcaZnBins = 400;
  const float dcaZlimit = 10.;
  Double_t dcazBins[dcaZnBins + 1];
  const float deltaDCAz = 2.f * dcaZlimit / dcaZnBins;
  for (int i = 0; i <= dcaZnBins; ++i)
    dcazBins[i] = i * deltaDCAz - dcaZlimit;
  
  
  fATOFsignal = new TH3F("fATOFsignal",
                         ";Centrality (%);p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                         nCentBins,centBins,nPtBins,pTbins,tofNbins,tofBins);
  fATPCcounts = new TH2F("fATPCcounts",";Centrality (%);p_{T} (GeV/c); TPC counts",
                         nCentBins,centBins,nPtBins,pTbins);
  fMDCAxy = new TH3F("fMDCAxy",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                     nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
  fMDCAz = new TH3F("fMDCAz",";Centrality (%);p_{T} (GeV/c); DCA_{z} (cm)",
                    nCentBins,centBins,nPtBins,pTbins,dcaZnBins,dcazBins);
  fMTOFsignal = new TH3F("fMTOFsignal",
                         ";Centrality (%);p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                         nCentBins,centBins,nPtBins,pTbins,tofNbins,tofBins);
  fMTPCcounts = new TH2F("fMTPCcounts",";Centrality (%);p_{T} (GeV/c); TPC counts",
                         nCentBins,centBins,nPtBins,pTbins);
  GetOutputList()->Add(fATOFsignal);
  GetOutputList()->Add(fATPCcounts);
  GetOutputList()->Add(fMDCAxy);
  GetOutputList()->Add(fMDCAz);
  GetOutputList()->Add(fMTOFsignal);
  GetOutputList()->Add(fMTPCcounts);
  
  fDeutBB = new TF1("fDeutBB",DeuteronTPC,0.3,6,0);
}

Bool_t AODSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either AODSelector::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  GetEntry(entry);
  
  if (centrality < 0) {
    fSkipEvent = kFALSE;
    if (-centrality < 13.f)
      fSkipEvent = Flatten(-centrality);
    
    if (!fSkipEvent) {
      fFlattenedCentrality->Fill(-centrality);
      fCentralityClasses->Fill(-centrality);
    }
    fCentrality->Fill(-centrality);
  }
  if (fSkipEvent) return kTRUE;
  
  const int cent = GetCentBin(TMath::Abs(centrality));
  if (cent < 0)
    return kTRUE;

  
  if (TMath::Abs(eta) > 0.8)
    return kTRUE;
  
  if (TMath::Abs(chi2NDF) > kChi2Cut)
    return kTRUE;
  
  Double_t pz = TMath::Sqrt(p * p - pT * pT);
  Double_t e = TMath::Sqrt(p * p + M2D);
  Double_t y = 0.5 * TMath::Log((e + pz) / (e - pz));
  if (TMath::Abs(y) > 0.5) {
    return kTRUE;
  }
  
  if (TPCnSignal < kTPCsig)
    return kTRUE;
  
  
  if (ITSnClust - ITSnSignal <= 0) return kTRUE;
  if (TMath::Abs(DCAz) > kDCAz) return kTRUE;
  
  if (TMath::Abs(DCAxy) > 0.5f) return kTRUE;
  
  Float_t c_pT = pT;
  if (c_pT > 0) {
    c_pT -= fCorrectionD(c_pT);
  } else {
    c_pT += fCorrectionAD(-c_pT);
  }
  
  if (TPCsignal > 0.7f * fDeutBB->Eval(pTPC) && TPCsignal < 1.3f * fDeutBB->Eval(pTPC)) {
    if(c_pT > 0.) {
      fMTPCcounts->Fill(centrality,c_pT);
      fMDCAxy->Fill(centrality,c_pT,DCAxy);
      fMDCAz->Fill(centrality,c_pT,DCAz);
    } else {
      fATPCcounts->Fill(centrality,-c_pT);
    }
    if (TOFtime > 0.f && length > 350.f) {
      Float_t beta = length / (2.99792457999999984e-02f * TOFtime);
      if (beta < (1.f - EPSILON)) {
        Float_t gamma = 1. / TMath::Sqrt(1.f - (beta * beta));
        const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
        if(c_pT > 0.) {
          fMTOFsignal->Fill(centrality,c_pT,dm);
        } else {
          fATOFsignal->Fill(centrality,-c_pT,dm);
        }
      }
    }
  }
  
  return kTRUE;
}

void AODSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
}

void AODSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  TFile f("nuclei.root","update");
  f.cd();
  TList l;
  l.SetOwner(kTRUE);
  l.SetName(kName);
  l.Add((TH1F*)((TH1F*)GetOutputList()->FindObject("fCentrality"))->Clone());
  l.Add((TH1F*)((TH1F*)GetOutputList()->FindObject("fFlattenedCentrality"))->Clone());
  l.Add((TH1F*)((TH1F*)GetOutputList()->FindObject("fCentralityClasses"))->Clone());
  l.Add((TH3F*)((TH3F*)GetOutputList()->FindObject("fATOFsignal"))->Clone());
  l.Add((TH2F*)((TH2F*)GetOutputList()->FindObject("fATPCcounts"))->Clone());
  l.Add((TH3F*)((TH3F*)GetOutputList()->FindObject("fMDCAxy"))->Clone());
  l.Add((TH3F*)((TH3F*)GetOutputList()->FindObject("fMDCAz"))->Clone());
  l.Add((TH3F*)((TH3F*)GetOutputList()->FindObject("fMTOFsignal"))->Clone());
  l.Add((TH2F*)((TH2F*)GetOutputList()->FindObject("fMTPCcounts"))->Clone());
  l.Write();
  f.Close();
}
