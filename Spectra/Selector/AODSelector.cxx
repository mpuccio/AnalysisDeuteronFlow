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
#include <TFractionFitter.h>

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

Double_t GausPol0(Double_t *x, Double_t *p) {
  return p[0] * TMath::Gaus(x[0],p[1],p[2]) + p[3];
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
  if (cent < 0) return -1;
  if (cent <= 10) return 0;
  else if (cent <= 20) return 1;
  else if (cent <= 40) return 2;
  else if (cent <= 60) return 3;
  else if (cent <= 80) return 4;
  return -1;
}

Int_t AODSelector::GetPtBin(float pt) {
  for (int i = 0; i < fkNBins; ++i) {
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
  fBeta2D = new TH2F("fBeta2D","TOF; #frac{p}{z} (GeV/c); #beta", 1000,0.01,10.,2250,0.1,1.);
  fBeta2DPt = new TH2F("fBeta2DPt","TOF; p_{T} (GeV/c); #beta", 1000,0.01,10.,2250,0.1,1.);
  fCentrality = new TH1F("fCentrality",";Centrality;Events / 1%",100,0,100);
  fGamma = new TH1F("fGamma",";#gamma;Entries",1000,1.f,1000.f);
  fDCAxy = new TH1F("fDCAxy",";DCA_{xy} (cm);Entries",4000,-4.f,4.f);
  fDCA2D = new TH2F("fDCA2D",";DCA_{xy} (cm);DCA_{z} (cm);Entries",500,-0.5f,0.5f,600,-0.6f,0.6f);
  fDCAz = new TH1F("fDCAz",";DCA_{z} (cm);Entries",4000,-4.f,4.f);
  fdEdxTPC = new TH2F("fdEdxTPC",";p (GeV/c);dE/dx (a.u.)",1000,0.1,10.,500,0,2000);
  fdEdxTPCpT = new TH2F("fdEdxTPCpT",";p_{T} (GeV/c);dE/dx (a.u.)",1000,0.1,10.,500,0,2000);
  fdEdxTPCSignal = new TH2F("fdEdxTPCSignal",";p (GeV/c);dE/dx (a.u.)",1000,0.1,10.,2000,0,2000);
  fTPCSignalN = new TH1F("fTPCSignalN",";#clusters for PID;Entries",161,-0.5f,160.5f);
  fdEdxTPCproj = new TH2F("fdEdxTPCproj",";p (GeV/c);dE/dx (a.u.)",93,0.4,3.2,2000,0,2000);
  BinLogAxis(fdEdxTPC);

  float binM[86];
  binM[0] = -3.4f;
  const float step = 6.8f / 85.f;
  for (int i = 1; i < 85; ++i) binM[i] = -3.4f + i * step;
  for (int cent = 0; cent < 5; ++cent) {
    for (int i = 0; i < fkNBins; ++i) {
      fSignalAD[cent * fkNBins + i] = new TH1F(Form("fSignalAD%i_%i",cent,i),Form("%4.2f<p_{T}#leq%4.2f;m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",fBins[i],fBins[i+1]),50,-2.0,2.0);
      fSignalD[cent * fkNBins + i] = new TH1F(Form("fSignalD%i_%i",cent,i),Form("%4.2f<p_{T}#leq%4.2f;m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",fBins[i],fBins[i+1]),50,-2.0,2.0);
      GetOutputList()->Add(fSignalD[cent * fkNBins + i]);
      GetOutputList()->Add(fSignalAD[cent * fkNBins + i]);
    }
    Double_t binDCA[] = {
      0.35, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
      2.0
    };
    fDdcaXY[cent] = new TH2F(Form("fDdcaXY%i",cent),";p_{T} (GeV/c); DCA_{xy} (cm)",10,binDCA,40,-0.5f,0.5f);
    fDdcaZ[cent] = new TH2F(Form("fDdcaZ%i",cent),";p_{T} (GeV/c); DCA_{z} (cm)",10,binDCA,40,-0.5f,0.5f);
    GetOutputList()->Add(fDdcaXY[cent]);
    GetOutputList()->Add(fDdcaZ[cent]);
  }
  
  GetOutputList()->Add(fdEdxTPC);
  GetOutputList()->Add(fdEdxTPCpT);
  GetOutputList()->Add(fdEdxTPCSignal);
  for (int i = 0; i < 5; ++i) {
    fdEdxTPCSignalCounts[i] = new TH1F(Form("fdEdxTPCSignalCounts%i",i),";p_{T};Entries",5,binDCA);
    GetOutputList()->Add(fdEdxTPCSignalCounts[i]);
    fdEdxTPCSignalCountsAD[i] = new TH1F(Form("fdEdxTPCSignalCountsAD%i",i),";p_{T};Entries",5,binDCA);
    GetOutputList()->Add(fdEdxTPCSignalCountsAD[i]);
  }
  
  for (int k = 0; k < 5; ++k) {
    for (int i = 0; i < 2; ++i) {
      fDCASignal[k * 2 + i] = new TH2F(Form("fDCASignal%i_%i",k,i),";m^{2} - m^{2}_{PDG} (GeV/c)^{2};DCA_{z} (cm);Entries",50,-2.0,2.0,40,-0.5f,0.5f);
      GetOutputList()->Add(fDCASignal[k * 2 + i]);
    }
  }
  
  GetOutputList()->Add(fdEdxTPCproj);
  GetOutputList()->Add(fTPCSignalN);
  GetOutputList()->Add(fCentrality);
  GetOutputList()->Add(fBeta);
  GetOutputList()->Add(fBeta2D);
  GetOutputList()->Add(fBeta2DPt);
  GetOutputList()->Add(fGamma);
  GetOutputList()->Add(fDCAxy);
  GetOutputList()->Add(fDCAz);
  GetOutputList()->Add(fDCA2D);
  
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
  
  if (ITSnClust - ITSnClustPID <= 0) return kTRUE;
  if (TMath::Abs(DCAz) > 2.f) return kTRUE;
  fDCAxy->Fill(DCAxy);
  fDCAz->Fill(DCAz);
  fDCA2D->Fill(DCAxy,DCAz);
  
  if (TMath::Abs(DCAxy) > 0.5f) return kTRUE;
  
  if (pTPC < 3.5) {
    if (TPCsignal > 0.7f * fDeutBB->Eval(pTPC) && TPCsignal < 1.3f * fDeutBB->Eval(pTPC)) {
      fdEdxTPCSignal->Fill(pTPC,TPCsignal);
      if (TMath::Abs(pT) < 0.8f) {
        if (pT >= 0.35f) {
          fdEdxTPCSignalCounts[cent]->Fill(pT);
          fDdcaXY[cent]->Fill(pT,DCAxy);
          fDdcaZ[cent]->Fill(pT,DCAz);
        } else if (pT <= -0.35f) {
          fdEdxTPCSignalCountsAD[cent]->Fill(TMath::Abs(pT));
        }
      } else if (TOFtime > 0.f && length > 0.f) {
        Float_t beta = length / (2.99792457999999984e-02 * TOFtime);
        fBeta->Fill(beta);
        fBeta2D->Fill(pTPC,beta);
        fBeta2DPt->Fill(TMath::Abs(pT),beta);
        if (beta < (1.f - EPSILON)) {
          Float_t gamma = 1 / TMath::Sqrt(1 - (beta * beta));
          fGamma->Fill(gamma);
          const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
          const int j = GetPtBin(TMath::Abs(pT));
          if (j < 0) {
            return kTRUE;
          }
          
          if(pT > 0) {
            fSignalD[fkNBins * cent + j]->Fill(dm);
            fDdcaXY[cent]->Fill(pT,DCAxy);
            fDdcaZ[cent]->Fill(pT,DCAz);
            if (j < 2) {
              fDCASignal[2 * cent + j]->Fill(dm, DCAxy);
            }
          } else {
            fSignalAD[fkNBins * cent + j]->Fill(dm);
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
      fBeta2DPt->Fill(TMath::Abs(pT),beta);
      if (beta < (1.f - EPSILON)) {
        Float_t gamma = 1 / TMath::Sqrt(1 - (beta * beta));
        fGamma->Fill(gamma);
        const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
        const int j = GetPtBin(TMath::Abs(pT));
        if (j < 0) {
          return kTRUE;
        }
        
        if(pT > 0) {
          fSignalD[cent * fkNBins + j]->Fill(dm);
        } else {
          fSignalAD[cent * fkNBins + j]->Fill(dm);
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
  
  fdEdxTPCSignal = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPCSignal"));
  if (fdEdxTPCSignal) {
    fdEdxTPCSignal->Write();
  }
  
  for (int i = 0; i < 5; ++i) {
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
  
  fBeta2DPt = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fBeta2DPt"));
  if (fBeta2DPt) {
    fBeta2DPt->Write();
  }
  
  fDCAz = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fDCAz"));
  if (fDCAz) {
    fDCAz->Write();
  }
  
  fDCAxy = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fDCAxy"));
  if (fDCAxy) {
    fDCAxy->Write();
  }
  
  fDCA2D = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fDCA2D"));
  if (fDCA2D) {
    fDCA2D->Write();
  }
  f.cd();
  f.mkdir("MassSpectra");
  f.mkdir("MassdEdx");
  for (int cent=0; cent < 5; ++cent) {
    f.cd();
    f.mkdir(Form("Signal%i",cent));
    f.cd(Form("Signal%i",cent));
    for (int j = 0; j < fkNBins; ++j) {
      fSignalAD[cent * fkNBins + j] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fSignalAD%i_%i",cent,j)));
      if (fSignalAD[cent * fkNBins + j]) {
        fSignalAD[cent * fkNBins + j]->Write();
      }
      
      fSignalD[cent * fkNBins + j] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fSignalD%i_%i",cent,j)));
      if (fSignalD[cent * fkNBins + j]) {
        fSignalD[cent * fkNBins + j]->Write();
      }
    }
    fDdcaXY[cent] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaXY%i",cent)));
    fDdcaZ[cent] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaZ%i",cent)));
    if (!fDdcaXY[cent] || !fDdcaZ[cent]) {f.Close(); return;}
    f.mkdir(Form("DCA%i",cent));
    f.cd(Form("DCA%i",cent));
    Double_t binDCA[] = {
      0.35, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
      2.0
    };
    
    for (int i = 0; i < 10; ++i) {
      TH1D *hprim = fDdcaXY[cent]->ProjectionY(Form("dcaxy_%i",i),i + 1, i + 2);
      TH1D *hseco = fDdcaZ[cent]->ProjectionY(Form("dcaz_%i",i),i + 1, i + 2);
      hprim->SetTitle(Form("%4.2f < p_{T} #leq %4.2f;DCA_{xy} (cm);Entries",binDCA[i],binDCA[i+1]));
      hseco->SetTitle(Form("%4.2f < p_{T} #leq %4.2f;DCA_{z} (cm);Entries",binDCA[i],binDCA[i+1]));
      hprim->Write();
      hseco->Write();
    }

  }
  
  f.cd();
  f.mkdir("DCASignal");
  f.cd("DCASignal");
  for (int k = 0; k < 2; ++k) {
    for (int i = 0; i < 5; ++i) {
      fDCASignal[k * 5 + i] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDCASignal%i_%i",i,k)));
      if (fDCASignal[k * 5 + i]) {
        fDCASignal[k * 5 + i]->Write();
      }
    }
  }

  f.Close();
}
