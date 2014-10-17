#define DeutSelector_cxx
// The class definition in DeutSelector.h has been generated automatically
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
// Root > T->Process("DeutSelector.C")
// Root > T->Process("DeutSelector.C","some options")
// Root > T->Process("DeutSelector.C+")
//

#include "DeutSelector.h"
#include <TH2.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TStyle.h>
#include <Riostream.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>

#define EPSILON 1E-5

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
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t aa = TMath::Power(beta,kp4);
  Double_t bb = TMath::Power(1./bg,kp5);
  bb=TMath::Log(kp3+bb);
  return (kp2-aa-bb)*kp1/aa;
}

Double_t DeuteronTPC(Double_t *x, Double_t *par) {
  // Deuteron expected signal in TPC, taken from AntiHe4 task
  return BetheBlochAleph(x[0]/1.875612,4.69637,7.51827,0.0183746,2.60,2.7);
}

void DeutSelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();
}

void DeutSelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  fBeta = new TH1D("fBeta",";#beta;Entries",100,0.f,1.f);
  fGamma = new TH1D("fGamma",";#gamma;Entries",1000,1.f,1000.f);
  fdEdxTPC = new TH2F("fdEdxTPC",";p (GeV/c);dE/dx (a.u.)",500,0.1,5.,500,0,2000);
  fdEdxTPCSignal = new TH2F("fdEdxTPCSignal",";p (GeV/c);dE/dx (a.u.)",500,0.1,5.,2000,0,2000);
  fdEdxTPCproj = new TH2F("fdEdxTPCproj",";p (GeV/c);dE/dx (a.u.)",93,0.4,3.2,2000,0,2000);
  BinLogAxis(fdEdxTPC);
  fSignal = new TH1D("fSignal",";",600,0.1,6.1);
  GetOutputList()->Add(fdEdxTPC);
  GetOutputList()->Add(fdEdxTPCSignal);
  GetOutputList()->Add(fdEdxTPCproj);
  GetOutputList()->Add(fSignal);
  GetOutputList()->Add(fBeta);
  GetOutputList()->Add(fGamma);
  
  fDeutBB = new TF1("fDeutBB",DeuteronTPC,0.3,6,0);

}

Bool_t DeutSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either DeutSelector::GetEntry() or TBranch::GetEntry()
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
  if (pTPC == p) {
    return kTRUE;
  }
  fdEdxTPC->Fill(pTPC,TPCsignal);
  fdEdxTPCproj->Fill(pTPC,TPCsignal);
  
  if (TOFtime > 0.f && length > 0.f) {
    Float_t beta = length / (2.99792457999999984e-02 * TOFtime);
    fBeta->Fill(beta);
    if (beta < (1.f - EPSILON)) {
      Float_t gamma = 1 / TMath::Sqrt(1 - (beta * beta));
      fGamma->Fill(gamma);
      fSignal->Fill(p / (beta * gamma));
    }
  }
  
  if (TPCsignal > 0.7f * fDeutBB->Eval(pTPC) && TPCsignal < 1.3f * fDeutBB->Eval(pTPC)) {
    fdEdxTPCSignal->Fill(pTPC,TPCsignal);
  }
  return kTRUE;
}

void DeutSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  delete fDeutBB;
}

void DeutSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  TFile f("final.root","recreate");
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

  fdEdxTPCSignal = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPCSignal"));
  if (fdEdxTPCSignal) {
    fdEdxTPCSignal->Write();
  }

  fdEdxTPCproj = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPCproj"));
  if (fdEdxTPCproj) {
    fdEdxTPCproj->Write();
  }
  
  fBeta = dynamic_cast<TH1D*>(GetOutputList()->FindObject("fBeta"));
  if (fBeta) {
    fBeta->Write();
  }

  fGamma = dynamic_cast<TH1D*>(GetOutputList()->FindObject("fGamma"));
  if (fGamma) {
    fGamma->Write();
  }

  fSignal = dynamic_cast<TH1D*>(GetOutputList()->FindObject("fSignal"));
  if (fSignal) {
    fSignal->Write();
  }
  f.Close();

}
