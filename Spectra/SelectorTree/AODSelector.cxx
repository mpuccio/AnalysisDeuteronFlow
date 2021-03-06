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
#include <Riostream.h>

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

void AODSelector::BinLogAxis(const TH1 *h) {
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


void AODSelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
}

Int_t AODSelector::GetCentBin(float cent) {
  if (cent < fCentralityBins[0]) return -1;

  for (int i = 0; i < fCentralityBins.GetSize() - 1; ++i)
    if (cent <= fCentralityBins[i + 1])
      return i;
  
  return -1;
}

Int_t AODSelector::GetPtBin(float pt) {
  for (int i = 0; i < fBins.GetSize() - 1; ++i) {
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
  enum {kEtaMin=1,kEtaMax,kYMin,kYMax,kTPCsig,kTPCchi2,kSPDrec,kDCAxy,kDCAz,kRecreate};
  TList *l = GetInputList();
  TH2F *h2 = (TH2F*)l->FindObject("hBins");
  fBins = *(h2->GetYaxis()->GetXbins());
  fCentralityBins = *(h2->GetXaxis()->GetXbins());
  TH1D *cuts = (TH1D*)l->FindObject("hCuts");
  fRequireEtaMin = cuts->GetBinContent(kEtaMin);
  fRequireEtaMax = cuts->GetBinContent(kEtaMax);
  fRequireYmin = cuts->GetBinContent(kYMin);
  fRequireYmax = cuts->GetBinContent(kYMax);
  fRequireTPCsignal = cuts->GetBinContent(kTPCsig);
  fRequireMaxChi2 = cuts->GetBinContent(kTPCchi2);
  fRequireSPDrecPoints = cuts->GetBinContent(kSPDrec);
  fRequireMaxDCAxy = cuts->GetBinContent(kDCAxy);
  fRequireMaxDCAz = cuts->GetBinContent(kDCAz);

  
  cout << "======= CUTS SUMMARY =======" << endl;
  cout << "Eta range " << fRequireEtaMin << " " << fRequireEtaMax << endl;
  cout << "Y range " << fRequireYmin << " " << fRequireYmax << endl;
  cout << "Minimum number of TPC clusters " << fRequireTPCsignal << endl;
  cout << "Maximum chi2 " << fRequireMaxChi2 << endl;
  cout << "Maximum DCAs (xy and z) " << fRequireMaxDCAxy << " " << fRequireMaxDCAz << endl;
  cout << "============================" << endl;
  cout << "\n======= BINS SUMMARY =======" << endl;
  h2->Print("base");
  cout << "============================" << endl;

  
  
  const Int_t nPtBins = fBins.GetSize() - 1;
  const Int_t nCentBins = fCentralityBins.GetSize() - 1;
  const Double_t *pTbins = fBins.GetArray();
  const Double_t *centBins = fCentralityBins.GetArray();
  
  const Int_t nDCAbins = 80;
  const Float_t range = 0.5;
  Double_t dcaBins[nDCAbins+1];
  for (int i = 0; i <= nDCAbins; ++i)
    dcaBins[i] = i * range * 2 / nDCAbins - range;
  
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
  const int nTpcBins = 40;
  Double_t tpcBins[nTpcBins + 1];
  for (int i = 0; i <= nTpcBins; ++i) {
    tpcBins[i] = -4.f + i * 0.2f;
  }
  
  Double_t dcaBB[6] = {0.,0.5,0.75,1.,1.5,2};
  Double_t tpcBB[6] = {0,60.5,65.5,70.5,75.5,80.5};
  Double_t chiBB[7] = {0,3.5,4.,4.5,5.,5.5,6.};
  fTPCstudyClusters = new TH3F("fTPCstudyClusters",";Centrality;TPC clusters;p_{T}",nCentBins,centBins,5,tpcBB,8,pTbins);
  fTPCstudyChi2 = new TH3F("fTPCstudyChi2",";Centrality;TPC #chi^{2}/NDF;p_{T}",nCentBins,centBins,6,chiBB,8,pTbins);
  fTPCstudyDCAz = new TH3F("fTPCstudyDCAz",";Centrality;DCA_{z} (cm);p_{T}",nCentBins,centBins,5,dcaBB,8,pTbins);
  GetOutputList()->Add(fTPCstudyChi2);
  GetOutputList()->Add(fTPCstudyClusters);
  GetOutputList()->Add(fTPCstudyDCAz);
  
  fCentrality = new TH1F("fCentrality",";Centrality (%);Events / 1%;",100,0.,100.);
  fCentralityClasses = new TH1F("fCentralityClasses",";Centrality classes(%);Events / Class;",nCentBins,centBins);
  fFlattenedCentrality = new TH1F("fFlattenCentrality","Centrality distribution after the flattening;Centrality (%);Events / 1%;",100,0.,100.);
  GetOutputList()->Add(fCentrality);
  GetOutputList()->Add(fCentralityClasses);
  GetOutputList()->Add(fFlattenedCentrality);
  
  fTPCSignal = new TH2F("fTPCSignal","p_{T} (GeV/c);dE/dx (a.u.)",290,0.1,3.,500,0,2500);
  BinLogAxis(fTPCSignal);
  
  fATOFsignal = new TH3F("fATOFsignal",
                         ";Centrality (%);p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                         nCentBins,centBins,nPtBins,pTbins,tofNbins,tofBins);
  fATPCcounts = new TH2F("fATPCcounts",";Centrality (%);p_{T} (GeV/c); Entries",
                         nCentBins,centBins,nPtBins,pTbins);
  fMDCAxyTPC = new TH3F("fMDCAxyTPC",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                        nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
  fMDCAzTPC = new TH3F("fMDCAzTPC",";Centrality (%);p_{T} (GeV/c); DCA_{z} (cm)",
                       nCentBins,centBins,nPtBins,pTbins,dcaZnBins,dcazBins);
  fMDCAxyTOF = new TH3F("fMDCAxyTOF",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                        nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
  fMDCAzTOF = new TH3F("fMDCAzTOF",";Centrality (%);p_{T} (GeV/c); DCA_{z} (cm)",
                       nCentBins,centBins,nPtBins,pTbins,dcaZnBins,dcazBins);
  fMTOFsignal = new TH3F("fMTOFsignal",
                         ";Centrality (%);p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                         nCentBins,centBins,nPtBins,pTbins,tofNbins,tofBins);
  fMTPCcounts = new TH2F("fMTPCcounts",";Centrality (%);p_{T} (GeV/c); Entries",
                         nCentBins,centBins,nPtBins,pTbins);
  GetOutputList()->Add(fATOFsignal);
  GetOutputList()->Add(fATPCcounts);
  GetOutputList()->Add(fMDCAxyTPC);
  GetOutputList()->Add(fMDCAzTPC);
  GetOutputList()->Add(fMDCAxyTOF);
  GetOutputList()->Add(fMDCAzTOF);
  GetOutputList()->Add(fMTOFsignal);
  GetOutputList()->Add(fMTPCcounts);
  GetOutputList()->Add(fTPCSignal);
  GetOutputList()->Add(cuts);
  
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

  
  if (eta < fRequireEtaMin || eta > fRequireEtaMax)
    return kTRUE;
  
  Double_t pz = TMath::Sqrt(p * p - pT * pT);
  Double_t e = TMath::Sqrt(p * p + M2D);
  Double_t y = 0.5 * TMath::Log((e + pz) / (e - pz));
  if (y < fRequireYmin || y > fRequireYmax) {
    return kTRUE;
  }
  
  if (TMath::Abs(DCAxy) > fRequireMaxDCAxy) return kTRUE;
  if (ITSnClust - ITSnSignal < fRequireSPDrecPoints) return kTRUE;
  
  if (TPCsignal > 0.7f * fDeutBB->Eval(pTPC) && TPCsignal < 1.3f * fDeutBB->Eval(pTPC)) {
    if (TMath::Abs(DCAz) <= fRequireMaxDCAz && TPCnSignal >= fRequireTPCsignal)
      fTPCstudyChi2->Fill(centrality, TMath::Abs(chi2NDF), pT);
    if (TMath::Abs(chi2NDF) <= fRequireMaxChi2 && TPCnSignal >= fRequireTPCsignal)
      fTPCstudyDCAz->Fill(centrality, TMath::Abs(DCAz), pT);
    if (TMath::Abs(chi2NDF) <= fRequireMaxChi2 && TMath::Abs(DCAz) <= fRequireMaxDCAz)
      fTPCstudyClusters->Fill(centrality, TMath::Abs(TPCnSignal), pT);
  }
  
  if (TPCnSignal < fRequireTPCsignal)
    return kTRUE;
  
  if (TMath::Abs(DCAz) > fRequireMaxDCAz)
    return kTRUE;
  
  if (TMath::Abs(chi2NDF) > fRequireMaxChi2)
    return kTRUE;
  
  Float_t c_pT = pT;
  if (c_pT > 0) {
    c_pT -= fCorrectionD(c_pT);
  } else {
    c_pT += fCorrectionAD(-c_pT);
  }

  fTPCSignal->Fill(pTPC, TPCsignal);
  
  if (TPCsignal > 0.7f * fDeutBB->Eval(pTPC) && TPCsignal < 1.3f * fDeutBB->Eval(pTPC)) {
    if(c_pT > 0.) {
      fMTPCcounts->Fill(centrality,c_pT);
      fMDCAxyTPC->Fill(centrality,c_pT,DCAxy);
      fMDCAzTPC->Fill(centrality,c_pT,DCAz);
    } else {
      fATPCcounts->Fill(centrality,-c_pT);
    }
    if (TOFtime > 0.f && length > 350.f) {
      Float_t beta = length / (2.99792457999999984e-02f * TOFtime);
      if (beta < (1.f - EPSILON)) {
        Float_t gamma = 1. / TMath::Sqrt(1.f - (beta * beta));
        const float dm = p * p / (beta * beta * gamma * gamma) - M2D;
        if(c_pT > 0.) {
          if (TMath::Abs(dm) < 2.) {
            fMDCAxyTOF->Fill(centrality,c_pT,DCAxy);
            fMDCAzTOF->Fill(centrality,c_pT,DCAz);
          }
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
  TList *l = GetInputList();
  TNamed *name = (TNamed*)l->FindObject("taskName");
  TH1D *cuts = (TH1D*)l->FindObject("hCuts");
  fRecreate = cuts->GetBinContent(10) > 0.5;
  fTaskName = name->GetTitle();
  
  TFile f("nuclei.root",(fRecreate) ? "recreate" : "update");
  f.mkdir(fTaskName.Data());
  f.cd(fTaskName.Data());
  if (GetOutputList()->FindObject("fCentrality") &&
      GetOutputList()->FindObject("fFlattenCentrality") &&
      GetOutputList()->FindObject("fCentralityClasses") &&
      GetOutputList()->FindObject("fATOFsignal") &&
      GetOutputList()->FindObject("fATPCcounts") &&
      GetOutputList()->FindObject("fMDCAxyTPC") &&
      GetOutputList()->FindObject("fMDCAzTPC") &&
      GetOutputList()->FindObject("fMDCAxyTOF") &&
      GetOutputList()->FindObject("fMDCAzTOF") &&
      GetOutputList()->FindObject("fMTOFsignal") &&
      GetOutputList()->FindObject("fMTPCcounts") &&
      GetOutputList()->FindObject("fTPCstudyClusters") &&
      GetOutputList()->FindObject("fTPCstudyDCAz") &&
      GetOutputList()->FindObject("fTPCstudyChi2")) {
    cuts->Write();
    ((TH1F*)GetOutputList()->FindObject("fCentrality"))->Write();
    ((TH1F*)GetOutputList()->FindObject("fFlattenCentrality"))->Write();
    ((TH1F*)GetOutputList()->FindObject("fCentralityClasses"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fATOFsignal"))->Write();
    ((TH2F*)GetOutputList()->FindObject("fATPCcounts"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fMDCAxyTPC"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fMDCAzTPC"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fMDCAxyTOF"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fMDCAzTOF"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fMTOFsignal"))->Write();
    ((TH2F*)GetOutputList()->FindObject("fMTPCcounts"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fTPCstudyClusters"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fTPCstudyDCAz"))->Write();
    ((TH3F*)GetOutputList()->FindObject("fTPCstudyChi2"))->Write();
  }
  
  if (GetOutputList()->FindObject("fTPCSignal")) {
    fTPCSignal = (TH2F*)GetOutputList()->FindObject("fTPCSignal");
    fTPCSignal->Write();
  }
  f.Close();

}
