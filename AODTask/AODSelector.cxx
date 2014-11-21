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
  float bins[17] = {0.8f,1.0f,1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,3.0f,3.5f,4.0f,5.0f,6.f,8.f,10.f};
  for (int i = 0; i < 16; ++i) {
    if (pt < bins[i+1] && pt >= bins[i]) {
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
  fBeta2D = new TH2F("fBeta2D","TOF; #frac{p}{z} (GeV/c); #beta", 1000,0.01,10.,2250,0.1,1.);
  fCentrality = new TH1F("fCentrality",";Centrality;Events / 1%",100,0,100);
  fGamma = new TH1F("fGamma",";#gamma;Entries",1000,1.f,1000.f);
  fdEdxTPC = new TH2F("fdEdxTPC",";p (GeV/c);dE/dx (a.u.)",1000,0.1,10.,500,0,2000);
  fdEdxTPCpT = new TH2F("fdEdxTPCpT",";p_{T} (GeV/c);dE/dx (a.u.)",1000,0.1,10.,500,0,2000);
  fdEdxTPCSignal = new TH2F("fdEdxTPCSignal",";p (GeV/c);dE/dx (a.u.)",1000,0.1,10.,2000,0,2000);
  fTPCSignalN = new TH1F("fTPCSignalN",";#clusters for PID;Entries",161,-0.5f,160.5f);
  Double_t d[7] = {0.35,0.5,0.6,0.7,0.8,1.};
  fdEdxTPCproj = new TH2F("fdEdxTPCproj",";p (GeV/c);dE/dx (a.u.)",93,0.4,3.2,2000,0,2000);
  BinLogAxis(fdEdxTPC);
  float bins[17] = {0.8f,1.0f,1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,3.0f,3.5f,4.0f,5.0f,6.f,8.f,10.f};
  fBins[0] = bins[0];
  for (int cent = 0; cent < 3; ++cent) {
    for (int i = 0; i < 16; ++i) {
      fBins[i+1] = bins[i+1];
      float hlim = (i > 12) ? 2.4f : 2.4f;
      float llim = (i > 12) ? -3.4f : -2.4f;
      int nbin = (i > 12) ? 85 : 60;
      fSignal[cent*16+i] = new TH1F(Form("fSignal%i_%i",cent,i),Form("%4.1f #leq p_{T} < %4.1f ; m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",fBins[i],fBins[i+1]),nbin,llim,hlim);
      GetOutputList()->Add(fSignal[cent*16+i]);
      fSignalAD[cent*16+i] = new TH1F(Form("fSignalAD%i_%i",cent,i),Form("%4.1f #leq p_{T} < %4.1f ; m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",fBins[i],fBins[i+1]),nbin,llim,hlim);
      GetOutputList()->Add(fSignalAD[cent*16+i]);
      fMassSpectra[cent*16+i] = new TH1F(Form("fMassSpectra%i_%i",cent,i),Form("%4.1f #leq p_{T} < %4.1f ; m (GeV/c);Entries",fBins[i],fBins[i+1]),100,0,3.);
      GetOutputList()->Add(fMassSpectra[cent*16+i]);
      fMassSpectraAD[cent*16+i] = new TH1F(Form("fMassSpectraAD%i_%i",cent,i),Form("%4.1f #leq p_{T} < %4.1f ; m (GeV/c);Entries",fBins[i],fBins[i+1]),100,0,3.);
      GetOutputList()->Add(fMassSpectraAD[cent*16+i]);
      fMassdEdxD[cent*16+i] = new TH2F(Form("fMassdEdxD%i_%i",cent,i),Form("%4.1f #leq p_{T} < %4.1f ; m (GeV/c); #beta; Entries",fBins[i],fBins[i+1]),600,0,3.,500,0,500);
      GetOutputList()->Add(fMassdEdxD[cent*16+i]);
      fMassdEdxAD[cent*16+i] = new TH2F(Form("fMassdEdxAD%i_%i",cent,i),Form("%4.1f #leq p_{T} < %4.1f ; m (GeV/c); #beta; Entries",fBins[i],fBins[i+1]),600,0,3.,500,0,500);
      GetOutputList()->Add(fMassdEdxAD[cent*16+i]);
    }
  }
  fCompleteSignalAD = new TH2F("fCompleteSignalAD",";p_{T};m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",112,0.8,12,85,-3.4,3.4);
  fCompleteSignalD = new TH2F("fCompleteSignalD",";p_{T};m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",112,0.8,12,85,-3.4,3.4);
  GetOutputList()->Add(fCompleteSignalD);
  GetOutputList()->Add(fCompleteSignalAD);
  
  GetOutputList()->Add(fdEdxTPC);
  GetOutputList()->Add(fdEdxTPCpT);
  GetOutputList()->Add(fdEdxTPCSignal);
  for (int i = 0; i < 3; ++i) {
    fdEdxTPCSignalCounts[i] = new TH1F(Form("fdEdxTPCSignalCounts%i",i),";p_{T};Entries",5,d);
    GetOutputList()->Add(fdEdxTPCSignalCounts[i]);
    fdEdxTPCSignalCountsAD[i] = new TH1F(Form("fdEdxTPCSignalCountsAD%i",i),";p_{T};Entries",5,d);
    GetOutputList()->Add(fdEdxTPCSignalCountsAD[i]);
  }
  GetOutputList()->Add(fdEdxTPCproj);
  GetOutputList()->Add(fTPCSignalN);
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
  
  if (TPCnSignal < 70)
    return kTRUE;
  
  fTPCSignalN->Fill(TPCnSignal);
  
  const int cent = GetCentBin(centrality);
  if (cent < 0) return kTRUE;
  
  
  fdEdxTPC->Fill(pTPC,TPCsignal);
  fdEdxTPCpT->Fill(TMath::Abs(pT),TPCsignal);
  fdEdxTPCproj->Fill(pTPC,TPCsignal);
  
  if (centrality != fPrevious) {
    fPrevious = centrality;
    fCentrality->Fill(centrality);
  }
  
  if (pTPC < 1.3) {
    if (TPCsignal > 0.7f * fDeutBB->Eval(pTPC) && TPCsignal < 1.3f * fDeutBB->Eval(pTPC)) {
      fdEdxTPCSignal->Fill(pTPC,TPCsignal);
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
          const int j = GetPtBin(TMath::Abs(pT));
          if (j < 0) {
            return kTRUE;
          }
          
          if(pT > 0) {
            fSignal[16 * cent + j]->Fill(dm);
            fMassSpectra[16 * cent + j]->Fill(p/(beta*gamma));
            fMassdEdxD[16 * cent + j]->Fill(p/(beta*gamma),TPCsignal);
            fCompleteSignalD->Fill(pT,dm);
          } else {
            fSignalAD[16 * cent + j]->Fill(dm);
            fMassSpectraAD[16 * cent + j]->Fill(p/(beta*gamma));
            fMassdEdxAD[16 * cent + j]->Fill(p/(beta*gamma),TPCsignal);
            fCompleteSignalAD->Fill(-pT,dm);
          }     
        }
      }
    }
  } else {
    fdEdxTPCSignal->Fill(pTPC,TPCsignal);
    if (TOFtime > 0.f && length > 0.f) {
      Float_t beta = length / (2.99792457999999984e-02 * TOFtime);
      fBeta->Fill(beta);
      fBeta2D->Fill(pTPC,beta);
      if (beta < (1.f - EPSILON)) {
        Float_t gamma = 1 / TMath::Sqrt(1 - (beta * beta));
        fGamma->Fill(gamma);
        const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
        const int j = GetPtBin(TMath::Abs(pT));
        if (j < 0) {
          return kTRUE;
        }
        
        if(pT > 0) {
          fSignal[16 * cent + j]->Fill(dm);
          fMassSpectra[16 * cent + j]->Fill(p/(beta*gamma));
          fMassdEdxD[16 * cent + j]->Fill(p/(beta*gamma),TPCsignal);
          fCompleteSignalD->Fill(pT,dm);
        } else {
          fSignalAD[16 * cent + j]->Fill(dm);
          fMassSpectraAD[16 * cent + j]->Fill(p/(beta*gamma));
          fMassdEdxAD[16 * cent + j]->Fill(p/(beta*gamma),TPCsignal);
          fCompleteSignalAD->Fill(-pT,dm);
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
  TFile f("AODSelector.root","recreate");
  fdEdxTPC = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPC"));
  if (fdEdxTPC) {
    fdEdxTPC->Write();
  }
  
  fdEdxTPCpT = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPCpT"));
  if (fdEdxTPCpT) {
    fdEdxTPCpT->Write();
  }
  
  float bins[17] = {0.8f,1.0f,1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,3.0f,3.5f,4.0f,5.0f,6.f,8.f,10.f};
  fBins[0] = bins[0];
  for (int i = 0; i < 16; ++i)
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
  
  fCompleteSignalAD = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fCompleteSignalAD"));
  if (fCompleteSignalAD) {
    fCompleteSignalAD->Write();
  }
  
  fCompleteSignalD = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fCompleteSignalD"));
  if (fCompleteSignalD) {
    fCompleteSignalD->Write();
  }
  
  fTPCSignalN = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fTPCSignalN"));
  if (fTPCSignalN) {
    fTPCSignalN->Write();
  }
  
  
  fBeta = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fBeta"));
  if (fBeta) {
    fBeta->Write();
  }
  
  fGamma = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fGamma"));
  if (fGamma) {
    fGamma->Write();
  }
  
  fBeta2D = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fBeta2D"));
  if (fBeta2D) {
    fBeta2D->Write();
  }
  
  f.mkdir("MassSpectra");
  f.mkdir("MassdEdx");
  f.mkdir("Signal");
  for (int cent=0; cent < 3; ++cent) {
    for (int i = 0; i < 16; ++i) {
      f.cd("Signal");
      fSignal[cent*16+i] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fSignal%i_%i",cent,i)));
      if (fSignal[cent*16+i])
        fSignal[cent*16+i]->Write();

      fSignalAD[cent*16+i] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fSignalAD%i_%i",cent,i)));
      if (fSignalAD[cent*16+i])
        fSignalAD[cent*16+i]->Write();
    }
    for (int i = 0; i < 16; ++i) {
      f.cd("MassSpectra");
      fMassSpectra[cent*16+i] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fMassSpectra%i_%i",cent,i)));
      if (fMassSpectra[cent*16+i])
        fMassSpectra[cent*16+i]->Write();
      
      fMassSpectraAD[cent*16+i] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fMassSpectraAD%i_%i",cent,i)));
      if (fMassSpectraAD[cent*16+i])
        fMassSpectraAD[cent*16+i]->Write();
    }
    for (int i = 0; i < 16; ++i) {
      f.cd("MassdEdx");
      fMassdEdxD[cent*16+i] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fMassdEdxD%i_%i",cent,i)));
      if (fMassdEdxD[cent*16+i])
        fMassdEdxD[cent*16+i]->Write();
      
      fMassdEdxAD[cent*16+i] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fMassdEdxAD%i_%i",cent,i)));
      if (fMassdEdxAD[cent*16+i])
        fMassdEdxAD[cent*16+i]->Write();
    }
  }


  f.Close();

  
}
