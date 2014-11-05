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
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>

#define EPSILON 1E-5
#define MD 1.875612859f
#define M2D MD*MD


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

void AODSelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
}

Int_t AODSelector::GetCentBin(float cent) {
  Float_t centBin[] = {0.f,5.f,20.f,40.f,60.f};
  for (int i = 0; i < 4; ++i) {
    if (cent < centBin[i+1] && cent >= centBin[i]) {
      switch (i) {
        case 0:
          return 0;
          break;
        case 1:
          return -1;
          break;
        default:
          return i - 1;
          break;
      }
    }
  }
  return -1;
}

Int_t AODSelector::GetPtBin(float pt) {
  for (int i = 0; i < 13; ++i) {
    if (pt < fBins[i+1] && pt >= fBins[i]) {
      return i;
    }
  }
  return -1;
}

void AODSelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  fBeta = new TH1F("fBeta",";#beta;Entries",100,0.f,1.f);
  fBeta2D = new TH2F("fBeta2D","TOF; #frac{p}{z} (GeV/c); #beta", 500,0.01,5.,2250,0.1,1.);
  fCentrality = new TH1F("fCentrality",";Centrality;Events / 1%",100,0,100);
  fGamma = new TH1F("fGamma",";#gamma;Entries",1000,1.f,1000.f);
  fdEdxTPC = new TH2F("fdEdxTPC",";p (GeV/c);dE/dx (a.u.)",500,0.1,5.,500,0,2000);
  fdEdxTPCSignal = new TH2F("fdEdxTPCSignal",";p (GeV/c);dE/dx (a.u.)",500,0.1,5.,2000,0,2000);
  Double_t d[7] = {0.35,0.5,0.6,0.7,0.8,1.};
  fdEdxTPCproj = new TH2F("fdEdxTPCproj",";p (GeV/c);dE/dx (a.u.)",93,0.4,3.2,2000,0,2000);
  BinLogAxis(fdEdxTPC);
  float bins[14] = {0.8f,1.0f,1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,3.0f,3.5f,4.0f,5.0f};
  fBins[0] = bins[0];
  for (int cent = 0; cent < 3; ++cent) {
    for (int i = 0; i < 13; ++i) {
      fBins[i+1] = bins[i+1];
      fSignal[cent*13+i] = new TH1F(Form("fSignal%i_%i",cent,i),Form("%4.1f #leq p_{T} < %4.1f ; m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",fBins[i],fBins[i+1]),60,-2.4,2.4);
      GetOutputList()->Add(fSignal[cent*13+i]);
    }
  }
  GetOutputList()->Add(fdEdxTPC);
  GetOutputList()->Add(fdEdxTPCSignal);
  for (int i = 0; i < 3; ++i) {
    fdEdxTPCSignalCounts[i] = new TH1F(Form("fdEdxTPCSignalCounts%i",i),";p_{T};Entries",5,d);
    GetOutputList()->Add(fdEdxTPCSignalCounts[i]);
    fdEdxTPCSignalCountsAD[i] = new TH1F(Form("fdEdxTPCSignalCountsAD%i",i),";p_{T};Entries",5,d);
    GetOutputList()->Add(fdEdxTPCSignalCountsAD[i]);
  }
  GetOutputList()->Add(fdEdxTPCproj);
  GetOutputList()->Add(fCentrality);
  GetOutputList()->Add(fBeta);
  GetOutputList()->Add(fBeta2D);
  GetOutputList()->Add(fGamma);
  
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
  if (TMath::Abs(eta) > 0.8)
    return kTRUE;
  if (pTPC == p) {
    return kTRUE;
  }

  Double_t pz = TMath::Sqrt(p * p - pT * pT);
  Double_t e = TMath::Sqrt(p * p + M2D);
  Double_t y = 0.5 * TMath::Log((e + pz) / (e - pz));
  if (TMath::Abs(y) > 0.5) {
    return kTRUE;
  }
  
  const int cent = GetCentBin(centrality);
  if (cent < 0) return kTRUE;
  
  fdEdxTPC->Fill(pTPC,TPCsignal);
  fdEdxTPCproj->Fill(pTPC,TPCsignal);
  
  if (centrality != fPrevious) {
    fPrevious = centrality;
    fCentrality->Fill(centrality);
  }
  
  if (TPCsignal > 0.7f * fDeutBB->Eval(pTPC) && TPCsignal < 1.3f * fDeutBB->Eval(pTPC)) {
    fdEdxTPCSignal->Fill(pTPC,TPCsignal);
    if (pTPC > 5.f) return kTRUE;
    if (pTPC < 1.f) {
      fdEdxTPCSignalCounts[cent]->Fill(pT);
      fdEdxTPCSignalCountsAD[cent]->Fill(-pT);
    }
    if (TOFtime > 0.f && length > 0.f) {
      Float_t beta = length / (2.99792457999999984e-02 * TOFtime);
      fBeta->Fill(beta);
      fBeta2D->Fill(pTPC,beta);
      if (beta < (1.f - EPSILON)) {
        Float_t gamma = 1 / TMath::Sqrt(1 - (beta * beta));
        fGamma->Fill(gamma);
        const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
        const int j = GetPtBin(-pT);
        if (j < 0) {
          return kTRUE;
        }
        fSignal[13 * cent + j]->Fill(dm);
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
  TFile f("AODSelector.root","recreate");
  TF1 fbb("BBAleph",DeuteronTPC,0.4f,5.f,0);
  TCanvas c;
  fdEdxTPC = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPC"));
  if (fdEdxTPC) {
    c.cd();
    fdEdxTPC->DrawClone();
    fbb.DrawClone("same");
    f.cd();
    c.Write();
    fdEdxTPC->Write();
  }
  
  float bins[14] = {0.8f,1.0f,1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,3.0f,3.5f,4.0f,5.0f};
  fBins[0] = bins[0];
  for (int i = 0; i < 13; ++i)
    fBins[i+1] = bins[i+1];
  
  fdEdxTPCSignal = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPCSignal"));
  if (fdEdxTPCSignal) {
    fdEdxTPCSignal->Write();
  }
  
  for (int i = 0; i < 3; ++i) {
    fdEdxTPCSignalCounts[i] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fdEdxTPCSignalCounts%i",i)));
    fdEdxTPCSignalCountsAD[i] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fdEdxTPCSignalCountsAD%i",i)));
    if (fdEdxTPCSignalCounts[i])
      fdEdxTPCSignalCounts[i]->Write();
    if (fdEdxTPCSignalCountsAD[i])
      fdEdxTPCSignalCountsAD[i]->Write();
  }
  
  fCentrality = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fCentrality"));
  if (fCentrality) {
    fCentrality->Write();
  }
  
  fdEdxTPCproj = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPCproj"));
  if (fdEdxTPCproj) {
    fdEdxTPCproj->Write();
  }
  
  fBeta = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fBeta"));
  if (fBeta) {
    fBeta->Write();
  }
  
  fGamma = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fGamma"));
  if (fGamma) {
    fGamma->Write();
  }
  
  for (int cent=0; cent < 3; ++cent) {
    for (int i = 0; i < 13; ++i) {
      fSignal[cent*13+i] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fSignal%i_%i",cent,i)));
      if (fSignal[cent*13+i]) {
        fSignal[cent*13+i]->Write();
      }
    }
  }
  
  fBeta2D = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fBeta2D"));
  if (fBeta2D) {
    fBeta2D->Write();
  }
  f.Close();

  
}
