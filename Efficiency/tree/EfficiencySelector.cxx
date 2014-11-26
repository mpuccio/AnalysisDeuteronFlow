#define EfficiencySelector_cxx
// The class definition in EfficiencySelector.h has been generated automatically
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
// Root > T->Process("EfficiencySelector.cxx")
// Root > T->Process("EfficiencySelector.cxx","some options")
// Root > T->Process("EfficiencySelector.cxx+")
//

#include "EfficiencySelector.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <Riostream.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>


void EfficiencySelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
}

void EfficiencySelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  fDYield = new TH1F("fDYield",";p_{T} (GeV/c);Number of d tracks",120,0.1,6.1);
  fAntiDYield = new TH1F("fAntiDYield",";p_{T} (GeV/c);Number of #bar{d} tracks",120,0.1,6.1);
  fDYieldTOF = new TH1F("fDYieldTOF",";p_{T} (GeV/c);Number of d tracks",120,0.1,6.1);
  fAntiDYieldTOF = new TH1F("fAntiDYieldTOF",";p_{T} (GeV/c);Number of #bar{d} tracks",120,0.1,6.1);
  fAntiDMCYield = new TH1F("fAntiDMCYield",";p_{T} (GeV/c);Number of MC #bar{d}",120,0.1,6.1);
  fDMCYield = new TH1F("fDMCYield",";p_{T} (GeV/c);Number of MC d",120,0.1,6.1);
  
  GetOutputList()->Add(fAntiDMCYield);
  GetOutputList()->Add(fAntiDYield);
  GetOutputList()->Add(fAntiDYieldTOF);
  GetOutputList()->Add(fDMCYield);
  GetOutputList()->Add(fDYield);
  GetOutputList()->Add(fDYieldTOF);
}

Bool_t EfficiencySelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either EfficiencySelector::GetEntry() or TBranch::GetEntry()
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
  
  if (TMath::Abs(eta) < 0.8 && IsPrimary && TMath::Abs(yMC) < 0.5) {
    if (ITSnClusters > 1 && TPCnClusters > 70u && TPCnSignal > 70u) { // reconstructed (anti) deuterons
      if (pTMC > 0.f) { // deuteron
        fDYield->Fill(TMath::Abs(pTMC));
        if (beta > 0.f && beta < 1.f) fDYieldTOF->Fill(TMath::Abs(pTMC));
      } else {        // anti-deuteron
        fAntiDYield->Fill(TMath::Abs(pTMC));
        if (beta > 0.f && beta < 1.f) fAntiDYieldTOF->Fill(TMath::Abs(pTMC));
      }
    }
    if (pTMC > 0.f) fDMCYield->Fill(TMath::Abs(pTMC));
    else            fAntiDMCYield->Fill(TMath::Abs(pTMC));
  }
  
  return kTRUE;
}

void EfficiencySelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
}

void EfficiencySelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  fAntiDMCYield = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fAntiDMCYield"));
  fAntiDYield = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fAntiDYield"));
  fAntiDYieldTOF = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fAntiDYieldTOF"));
  fDMCYield = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fDMCYield"));
  fDYield = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fDYield"));
  fDYieldTOF = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fDYieldTOF"));
  if (!fAntiDYieldTOF || !fAntiDYield || !fAntiDMCYield || !fDMCYield || !fDYieldTOF || !fDYield) {
    return;
  }
  
  fDMCYield->Rebin(4);
  fAntiDMCYield->Rebin(4);
  fDYield->Rebin(4);
  fAntiDYield->Rebin(4);
  fAntiDYieldTOF->Rebin(4);
  fDYieldTOF->Rebin(4);
  TCanvas *effCv = new TCanvas();
  effCv->cd();
  
  TEfficiency *effAcc = new TEfficiency(*fDYield,*fDMCYield);
  effAcc->Paint("");
  TGraphAsymmErrors *gr = effAcc->GetPaintedGraph();
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetRangeUser(0.f,1.1f);
  effAcc->SetName("d");
  effAcc->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
  effAcc->SetMarkerStyle(20);
  effAcc->SetMarkerSize(0.7);
  effAcc->SetLineColor(kRed);
  effAcc->SetMarkerColor(kRed);
  effAcc->Draw("AP");
  
  
  TEfficiency *effAccTof = new TEfficiency(*fDYieldTOF,*fDMCYield);
  effAccTof->SetName("dTOF");
  effAccTof->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
  effAccTof->SetMarkerStyle(33);
  effAccTof->SetMarkerSize(0.7);
  effAccTof->SetLineColor(kRed);
  effAccTof->SetMarkerColor(kRed);
  effAccTof->Draw("same");
  
  TEfficiency *effAccAD = new TEfficiency(*fAntiDYield,*fAntiDMCYield);
  effAccAD->SetName("ad");
  effAccAD->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
  effAccAD->SetMarkerStyle(20);
  effAccAD->SetMarkerSize(0.7);
  effAccAD->SetLineColor(kBlue);
  effAccAD->SetMarkerColor(kBlue);
  effAccAD->Draw("same");
  
  TEfficiency *effAccTofAD = new TEfficiency(*fAntiDYieldTOF,*fAntiDMCYield);
  effAccTofAD->SetName("adTOF");
  effAccTofAD->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
  effAccTofAD->SetMarkerStyle(33);
  effAccTofAD->SetMarkerSize(0.7);
  effAccTofAD->SetLineColor(kBlue);
  effAccTofAD->SetMarkerColor(kBlue);
  effAccTofAD->Draw("same");
  
  TLegend *ll = new TLegend(0.6,0.125,0.85,0.31);
  ll->SetFillColor(kWhite);
  ll->SetLineColor(kWhite);
  ll->AddEntry(effAcc,"d tracking","lp");
  ll->AddEntry(effAccTof,"d tracking + TOF","lp");
  ll->AddEntry(effAccAD,"#bar{d} tracking","lp");
  ll->AddEntry(effAccTofAD,"#bar{d} tracking + TOF","lp");
  ll->Draw();
}
