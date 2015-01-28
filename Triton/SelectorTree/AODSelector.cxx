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
#include <TRandom3.h>

Double_t SigmaITS(Double_t sig, Double_t p, Int_t nClust, Bool_t isDeuteron) {
  // isDeuteron == kFALSE -> triton
  // taken from AliITSPIDResponse
  
  Double_t bg = p;
  const Double_t bbParamDeu[5] = {76.43,-34.21,113.2,-18.12,0.6019};
  const Double_t bbParamTri[5] = {13.34,55.17,66.41,-6.601,-0.4134};
  Double_t parResolDeu3[3]={0.06918,0.02498,1.1};
  Double_t parResolDeu4[3]={0.06756,0.02078,1.05};
  Double_t parResolTri3[3]={0.07239,0.0192,1.1};
  Double_t parResolTri4[3]={0.06083,0.02579,1.15};
  const Double_t *par,*resPar;
  if (isDeuteron) {
    bg /= 1.875612;
    resPar = parResolDeu3;
    par = bbParamDeu;
    if (nClust == 4) {
      resPar = parResolDeu4;
    }
  } else {
    bg /= 2.808921;
    resPar = parResolTri3;
    par = bbParamTri;
    if (nClust == 4) {
      resPar = parResolTri4;
    }
  }
  const Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  const Double_t gamma=bg / beta;
  Double_t bb = 1.;
  
  if(gamma>=0. && beta>0.)
    bb = par[0] + par[1]/bg + par[2]/(bg*bg) + par[3]/(bg*bg*bg) + par[4]/(bg*bg*bg*bg);
  
  const Double_t &c = resPar[2];
  const Double_t &r = resPar[0] + resPar[1] * p;
  
  return (sig - bb) / (r * bb * c);
}


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

unsigned int Log2Int(const unsigned int v) {
  // Taken from https://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious
  static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000};
  register unsigned int r = (v & b[0]) != 0;
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
  if (cent < fCentralityClasses[0]) return -1;

  for (int i = 0; i < kNCent; ++i)
    if (cent <= fCentralityClasses[i + 1])
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
  fBeta = new TH1F("fBeta",";#beta;Entries",100,0.f,1.f);
  fBeta2D = new TH2F("fBeta2D","TOF; #frac{p}{z} (GeV/c); #beta", 1000,0.01,10.,2250,0.1,1.);
  fBeta2DPt = new TH2F("fBeta2DPt","TOF; p_{T} (GeV/c); #beta", 1000,0.01,10.,2250,0.1,1.);
  fCentrality = new TH1F("fCentrality",";Centrality;Events / 1%",100,0,100);
  fCentralityClass = new TH1F("fCentralityClass",";Centrality Class; Events / Class",kNCent,fCentralityClasses);
  fGamma = new TH1F("fGamma",";#gamma;Entries",1000,1.f,1000.f);
  fDCAxy = new TH1F("fDCAxy",";DCA_{xy} (cm);Entries",4000,-4.f,4.f);
  fDCA2D = new TH2F("fDCA2D",";DCA_{xy} (cm);DCA_{z} (cm);Entries",500,-0.5f,0.5f,600,-0.6f,0.6f);
  fDCAz = new TH1F("fDCAz",";DCA_{z} (cm);Entries",4000,-4.f,4.f);
  fdEdxTPC = new TH2F("fdEdxTPC",";p (GeV/c);dE/dx (a.u.)",1000,0.1,10.,500,0,2000);
  fdEdxTPCpT = new TH2F("fdEdxTPCpT",";p_{T} (GeV/c);dE/dx (a.u.)",1000,0.1,10.,500,0,2000);
  fdEdxTPCSignal = new TH2F("fdEdxTPCSignal",";p (GeV/c);dE/dx (a.u.)",1000,0.1,10.,2000,0,2000);
  fTPCSignalN = new TH1F("fTPCSignalN",";#clusters for PID;Entries",161,-0.5f,160.5f);
  fTriggerHist = new TH1F("fTriggerHist",";Trigger;Entries",5,-0.5,4.5);
  fdEdxTPCproj = new TH2F("fdEdxTPCproj",";p (GeV/c);dE/dx (a.u.)",93,0.4,3.2,2000,0,2000);
  fdEdxTriton = new TH2F("fdEdxTriton",";p (GeV/c);dE/dx (a.u.)",560,0.4,3.2,2000,0,2000);
  fdEdxDeuteron = new TH2F("fdEdxDeuteron",";p (GeV/c);dE/dx (a.u.)",560,0.4,3.2,2000,0,2000);

  BinLogAxis(fdEdxTPC);

  const Char_t* aTriggerNames[] = { "kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" };
  for (Int_t ii = 1; ii <= 5; ++ii)
    fTriggerHist->GetXaxis()->SetBinLabel(ii, aTriggerNames[ii - 1]);
  
  for (int i = 0; i < kNBinsTOF; ++i) {
    fTOFSignal[i] = new TH1F(Form("fTOFSignal%i",i),
                             Form("%4.2f<p_{T}#leq%4.2f;m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",
                                  fBins[i+kNBinsTPC],fBins[i+1+kNBinsTPC]),50,-2.0,2.0);
    GetOutputList()->Add(fTOFSignal[i]);
  }
  
  for (int cent = 0; cent < kNCent; ++cent) {
    for (int iB = 0; iB < kNBins; ++iB) {
      fdEdxTPCSignalSlicesD[cent * kNBins + iB] = new TH1F(Form("fTPCSignalD%i_%i",cent,iB),Form("%4.2f<p_{T}#leq%4.2f;p (GeV/c);dE/dx (a.u.)",fBins[iB],fBins[iB+1]),2500,0,2500);
      fdEdxTPCSignalSlicesAD[cent * kNBins + iB] = new TH1F(Form("fTPCSignalAD%i_%i",cent,iB),Form("%4.2f<p_{T}#leq%4.2f;p (GeV/c);dE/dx (a.u.)",fBins[iB],fBins[iB+1]),2500,0,2500);
      fdEdxITSSignalSlicesD[cent * kNBins + iB] = new TH1F(Form("fITSSignalD%i_%i",cent,iB),Form("%4.2f<p_{T}#leq%4.2f;p (GeV/c);dE/dx (a.u.)",fBins[iB],fBins[iB+1]),2500,0,2500);
      fdEdxITSSignalSlicesAD[cent * kNBins + iB] = new TH1F(Form("fITSSignalAD%i_%i",cent,iB),Form("%4.2f<p_{T}#leq%4.2f;p (GeV/c);dE/dx (a.u.)",fBins[iB],fBins[iB+1]),2500,0,2500);
      GetOutputList()->Add(fdEdxTPCSignalSlicesAD[cent * kNBins + iB]);
      GetOutputList()->Add(fdEdxTPCSignalSlicesD[cent * kNBins + iB]);
      GetOutputList()->Add(fdEdxITSSignalSlicesAD[cent * kNBins + iB]);
      GetOutputList()->Add(fdEdxITSSignalSlicesD[cent * kNBins + iB]);
    }
    for (int i = 0; i < kNBinsTOF; ++i) {
      int j = i + kNBinsTPC;
      fSignalAD[cent * kNBinsTOF + i] = new TH1F(Form("fSignalAD%i_%i",cent,i),
                                               Form("%4.2f<p_{T}#leq%4.2f;m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",fBins[j],fBins[j+1]),75,-2.0,4.0);
      fSignalD[cent * kNBinsTOF + i] = new TH1F(Form("fSignalD%i_%i",cent,i),
                                              Form("%4.2f<p_{T}#leq%4.2f;m^{2} - m^{2}_{PDG} (GeV/c)^{2};Entries",fBins[j],fBins[j+1]),75,-2.0,4.0);
      GetOutputList()->Add(fSignalD[cent * kNBinsTOF + i]);
      GetOutputList()->Add(fSignalAD[cent * kNBinsTOF + i]);
    }
    fDdcaXY[cent] = new TH2F(Form("fDdcaXY%i",cent),";p_{T} (GeV/c); DCA_{xy} (cm)",
                             10,fBins,160,-0.5f,0.5f);
    fDdcaZ[cent] = new TH2F(Form("fDdcaZ%i",cent),";p_{T} (GeV/c); DCA_{z} (cm)",
                            10,fBins,160,-0.5f,0.5f);
    GetOutputList()->Add(fDdcaXY[cent]);
    GetOutputList()->Add(fDdcaZ[cent]);
  }
  
  GetOutputList()->Add(fdEdxTPC);
  GetOutputList()->Add(fdEdxTPCpT);
  GetOutputList()->Add(fdEdxTPCSignal);
  for (int i = 0; i < kNCent; ++i) {
    fdEdxTPCSignalCounts[i] = new TH1F(Form("fdEdxTPCSignalCounts%i",i),";p_{T};Entries",kNBinsTPC,fBins);
    GetOutputList()->Add(fdEdxTPCSignalCounts[i]);
    fdEdxTPCSignalCountsAD[i] = new TH1F(Form("fdEdxTPCSignalCountsAD%i",i),";p_{T};Entries",kNBinsTPC,fBins);
    GetOutputList()->Add(fdEdxTPCSignalCountsAD[i]);
  }
  
  for (int k = 0; k < kNCent; ++k) {
    for (int i = 0; i < kNDCAbinsTOF; ++i) {
      fDCASignal[k * kNDCAbinsTOF + i] = new TH2F(Form("fDCASignal%i_%i",k,i),
                                       ";m^{2} - m^{2}_{PDG} (GeV/c)^{2};DCA_{z} (cm);Entries",
                                       75,-2.0,4.0,160,-0.5f,0.5f);
      GetOutputList()->Add(fDCASignal[k * kNDCAbinsTOF + i]);
    }
  }
  
  GetOutputList()->Add(fdEdxTPCproj);
  GetOutputList()->Add(fdEdxTriton);
  GetOutputList()->Add(fdEdxDeuteron);
  GetOutputList()->Add(fTPCSignalN);
  GetOutputList()->Add(fTriggerHist);
  GetOutputList()->Add(fCentrality);
  GetOutputList()->Add(fCentralityClass);
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
  
  if (centrality < 0) {
    fSkipEvent = kFALSE;
    if (-centrality < 13.f)
      fSkipEvent = Flatten(-centrality);
    
    if (!fSkipEvent) {
      fCentrality->Fill(-centrality);
      fCentralityClass->Fill(-centrality);
      fTriggerHist->Fill(Log2Int(trigger));
    }
  }
  if (fSkipEvent) return kTRUE;
  
  const int cent = GetCentBin(TMath::Abs(centrality));
  if (cent < 0)
    return kTRUE;

  
  if (TMath::Abs(eta) > 0.8)
    return kTRUE;
  
  if (TMath::Abs(chi2NDF) > 6)
    return kTRUE;
  
  Double_t pz = TMath::Sqrt(p * p - pT * pT);
  Double_t e = TMath::Sqrt(p * p + M2D);
  Double_t y = 0.5 * TMath::Log((e + pz) / (e - pz));
  if (TMath::Abs(y) > 0.5) {
    return kTRUE;
  }
  
  if (TPCnSignal < 70)
    return kTRUE;
  
  fTPCSignalN->Fill(TPCnSignal);
  
  
  fdEdxTPC->Fill(pTPC,TPCsignal);
  fdEdxTPCpT->Fill(TMath::Abs(pT),TPCsignal);
  fdEdxTPCproj->Fill(pTPC,TPCsignal);
  if (TMath::Abs(TPCsigmad) < 3.) {
    fdEdxDeuteron->Fill(pTPC,TPCsignal);
  }
  if (TMath::Abs(TPCsigmat) < 3.) {
    fdEdxTriton->Fill(pTPC,TPCsignal);
  }
  
  if (ITSnClust - ITSnSignal <= 0) return kTRUE;
  if (TMath::Abs(DCAz) > 2.f) return kTRUE;
  fDCAxy->Fill(DCAxy);
  fDCAz->Fill(DCAz);
  fDCA2D->Fill(DCAxy,DCAz);
  
  if (TMath::Abs(DCAxy) > 0.5f) return kTRUE;
  
  Float_t c_pT = pT;
  if (c_pT > 0) {
    c_pT -= fCorrectionD(c_pT);
  } else {
    c_pT += fCorrectionAD(-c_pT);
  }
  
  const int j = GetPtBin(TMath::Abs(c_pT));
  if (TMath::Abs(TPCsigmat) < 3.) {
    if (c_pT > 0) {
      fdEdxITSSignalSlicesD[cent * kNBins + j]->Fill(ITSsignal);
      fdEdxTPCSignalSlicesD[cent * kNBins + j]->Fill(ITSsignal);
    } else {
      fdEdxITSSignalSlicesAD[cent * kNBins + j]->Fill(ITSsignal);
      fdEdxTPCSignalSlicesAD[cent * kNBins + j]->Fill(ITSsignal);
    }
  }
  
  if (pTPC < 3.5f) {
    if (TMath::Abs(TPCsigmat) < 3.) {
      fdEdxTPCSignal->Fill(pTPC,TPCsignal);
      if (TMath::Abs(c_pT) < fBins[kNBinsTPC]) {
        if (c_pT >= fBins[0]) {
          fdEdxTPCSignalCounts[cent]->Fill(c_pT);
          fDdcaXY[cent]->Fill(c_pT,DCAxy);
          fDdcaZ[cent]->Fill(c_pT,DCAz);
        } else if (c_pT <= -fBins[0]) {
          fdEdxTPCSignalCountsAD[cent]->Fill(TMath::Abs(c_pT));
        }
      }
      if (TOFtime > 0.f && length > 0.f) {
        Float_t beta = length / (2.99792457999999984e-02f * TOFtime);
        fBeta->Fill(beta);
        fBeta2D->Fill(pTPC,beta);
        fBeta2DPt->Fill(TMath::Abs(c_pT),beta);
        if (beta < (1.f - EPSILON)) {
          Float_t gamma = 1. / TMath::Sqrt(1.f - (beta * beta));
          fGamma->Fill(gamma);
          const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
          if (j < kNBinsTPC) {
            return kTRUE;
          }
          
          if(c_pT > 0.) {
            fTOFSignal[j - kNBinsTPC]->Fill(dm);
            fSignalD[cent * kNBinsTOF + j - kNBinsTPC]->Fill(dm);
            fDdcaXY[cent]->Fill(c_pT, DCAxy);
            fDdcaZ[cent]->Fill(c_pT, DCAz);
            if (j - kNBinsTPC < kNDCAbinsTOF) {
              fDCASignal[cent * kNDCAbinsTOF + j - kNBinsTPC]->Fill(dm, DCAxy);
            }
          } else {
            fSignalAD[cent * kNBinsTOF + j - kNBinsTPC]->Fill(dm);
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
      fBeta2DPt->Fill(TMath::Abs(c_pT),beta);
      if (beta < (1.f - EPSILON)) {
        Float_t gamma = 1 / TMath::Sqrt(1 - (beta * beta));
        fGamma->Fill(gamma);
        const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
        if (j < kNBinsTPC) {
          return kTRUE;
        }
        
        if(c_pT > 0) {
          fSignalD[cent * kNBinsTOF + j - kNBinsTPC]->Fill(dm);
          fTOFSignal[j - kNBinsTPC]->Fill(dm);
        } else {
          fSignalAD[cent * kNBinsTOF + j - kNBinsTPC]->Fill(dm);
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
  
  for (int i = 0; i < kNCent; ++i) {
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
  
  fCentralityClass = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fCentralityClass"));
  if (fCentralityClass) {
    fCentralityClass->Write();
  }
  
  fdEdxTPCproj = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTPCproj"));
  if (fdEdxTPCproj) {
    fdEdxTPCproj->Write();
  }
  
  fdEdxTriton = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxTriton"));
  if (fdEdxTriton) {
    fdEdxTriton->Write();
  }
  
  fdEdxDeuteron = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fdEdxDeuteron"));
  if (fdEdxDeuteron) {
    fdEdxDeuteron->Write();
  }
  
  fTPCSignalN = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fTPCSignalN"));
  if (fTPCSignalN) {
    fTPCSignalN->Write();
  }
  
  fTriggerHist = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fTriggerHist"));
  if (fTriggerHist) {
    fTriggerHist->Write();
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
  for (int cent = 0; cent < kNCent; ++cent) {
    f.cd();
    f.mkdir(Form("Signal%i",cent));
    f.cd(Form("Signal%i",cent));
    for (int j = 0; j < kNBinsTOF; ++j) {
      fSignalAD[cent * kNBinsTOF + j] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fSignalAD%i_%i",cent,j)));
      if (fSignalAD[cent * kNBinsTOF + j]) {
        fSignalAD[cent * kNBinsTOF + j]->Write();
      }
      
      fSignalD[cent * kNBinsTOF + j] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fSignalD%i_%i",cent,j)));
      if (fSignalD[cent * kNBinsTOF + j]) {
        fSignalD[cent * kNBinsTOF + j]->Write();
      }
    }
    fDdcaXY[cent] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaXY%i",cent)));
    fDdcaZ[cent] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaZ%i",cent)));
    if (!fDdcaXY[cent] || !fDdcaZ[cent]) {f.Close(); return;}
    f.mkdir(Form("DCA%i",cent));
    f.cd(Form("DCA%i",cent));
    
    for (int i = 0; i < kNDCABins; ++i) {
      TH1D *hprim = fDdcaXY[cent]->ProjectionY(Form("dcaxy_%i",i),i + 1, i + 1);
      TH1D *hseco = fDdcaZ[cent]->ProjectionY(Form("dcaz_%i",i),i + 1, i + 1);
      hprim->SetTitle(Form("%4.2f < p_{T} #leq %4.2f;DCA_{xy} (cm);Entries",fBins[i],fBins[i + 1]));
      hseco->SetTitle(Form("%4.2f < p_{T} #leq %4.2f;DCA_{z} (cm);Entries",fBins[i],fBins[i + 1]));
      hprim->Write();
      hseco->Write();
    }
    f.cd();
    f.mkdir(Form("TPCs%i",cent));
    f.mkdir(Form("ITSs%i",cent));
    for (int iB = 0; iB < kNBins; ++iB) {
      fdEdxITSSignalSlicesD[cent * kNBins + iB] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fITSSignalD%i_%i",cent,iB)));
      fdEdxITSSignalSlicesAD[cent * kNBins + iB] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fITSSignalAD%i_%i",cent,iB)));
      fdEdxTPCSignalSlicesD[cent * kNBins + iB] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fTPCSignalD%i_%i",cent,iB)));
      fdEdxTPCSignalSlicesAD[cent * kNBins + iB] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fTPCSignalAD%i_%i",cent,iB)));
      if (!fdEdxITSSignalSlicesD[cent * kNBins + iB] || !fdEdxITSSignalSlicesAD[cent * kNBins + iB] ||
          !fdEdxTPCSignalSlicesD[cent * kNBins + iB] || !fdEdxTPCSignalSlicesAD[cent * kNBins + iB] ) {
        cout << "Missing TPC/ITS signal " << cent << "\t" << iB << endl;
        continue;
      }
      f.cd(Form("ITSs%i",cent));
      fdEdxITSSignalSlicesD[cent * kNBins + iB]->Write();
      fdEdxITSSignalSlicesAD[cent * kNBins + iB]->Write();
      f.cd();
      f.mkdir(Form("TPCs%i",cent));
      fdEdxTPCSignalSlicesD[cent * kNBins + iB]->Write();
      fdEdxTPCSignalSlicesAD[cent * kNBins + iB]->Write();
      f.cd();
    }
  }
  
  f.mkdir("TOFSignal");
  f.cd("TOFSignal");
  for (int j = 0; j < kNBinsTOF; ++j) {
    fTOFSignal[j] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fTOFSignal%i",j)));
    if (fTOFSignal[j]) {
      fTOFSignal[j]->Write();
    }
  }
    
  f.cd();
  f.mkdir("DCASignal");
  f.cd("DCASignal");
  for (int k = 0; k < kNCent; ++k) {
    for (int i = 0; i < kNDCAbinsTOF; ++i) {
      fDCASignal[k * kNDCAbinsTOF + i] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDCASignal%i_%i",i,k)));
      if (fDCASignal[k * kNDCAbinsTOF + i]) {
        fDCASignal[k * kNDCAbinsTOF + i]->Write();
      }
    }
  }

  f.Close();
}
