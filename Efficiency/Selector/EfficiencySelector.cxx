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
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TString.h>
#include <TObjArray.h>
#include <TF1.h>
#include <TRandom3.h>

static int GetCentrality(float cent) {
  if (cent < 0) return -1;
  if (cent <= 10) return 0;
  else if (cent <= 20) return 1;
  else if (cent <= 40) return 2;
  else if (cent <= 60) return 3;
  else if (cent <= 80) return 4;
  return -1;
}

Bool_t Skip(float cent) {
  if (cent > 0.6f) return kFALSE;
  if (cent < 0.f)  return kTRUE;
  float prob[15] = {
    0.058812,0.510204,0.606061,0.831025,0.884956,1.22951,1.53061,1.66667,1.74419,1.98675,
    2.0979  ,2.5641  ,2.52101 ,2.54237 ,3.09278
  };
  float rnd = gRandom->Rndm();
  return (rnd > prob[int(cent * 25.)]); // * 25 == / 0.04
}

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
  
  const double bins[] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f,8.0f,10.0f
  };
  
  fDYield = new TH1F("fDYield",";p_{T} (GeV/c);Number of d tracks",kNBins,bins);
  fAntiDYield = new TH1F("fAntiDYield",";p_{T} (GeV/c);Number of #bar{d} tracks",kNBins,bins);
  fDYieldTOF = new TH1F("fDYieldTOF",";p_{T} (GeV/c);Number of d tracks",kNBins,bins);
  fAntiDYieldTOF = new TH1F("fAntiDYieldTOF",";p_{T} (GeV/c);Number of #bar{d} tracks",kNBins,bins);
  fAntiDMCYield = new TH1F("fAntiDMCYield",";p_{T} (GeV/c);Number of MC #bar{d}",kNBins,bins);
  fDMCYield = new TH1F("fDMCYield",";p_{T} (GeV/c);Number of MC d",kNBins,bins);
  fPtCorrectionD = new TH2F("fPtCorrectionD","d;p_{T} (GeV/c);p_{T} - p_{T_{rec}};entries",57,0.3,6,200,-1,1);
  fPtCorrectionAD = new TH2F("fPtCorrectionAD","#bar{d};p_{T} (GeV/c);p_{T} - p_{T_{rec}};entries",57,0.3,6,200,-1,1);
  fCentrality = new TH1F("fCentrality",";Centrality;Events/1%",400,0,4);

  for (int i = 0; i < kNCent; ++i) {
    fAntiDMCYieldCent[i] = new TH1F(Form("fAntiDMCYieldCent%i",i),";p_{T} (GeV/c);Number of MC #bar{d}",kNBins,bins);
    fAntiDYieldCent[i] = new TH1F(Form("fAntiDYieldCent%i",i),";p_{T} (GeV/c);Number of #bar{d} tracks",kNBins,bins);
    fAntiDYieldTOFCent[i] = new TH1F(Form("fAntiDYieldTOFCent%i",i),";p_{T} (GeV/c);Number of #bar{d} tracks",kNBins,bins);
    fDMCYieldCent[i] = new TH1F(Form("fDMCYieldCent%i",i),";p_{T} (GeV/c);Number of MC d",kNBins,bins);
    fDYieldCent[i] = new TH1F(Form("fDYieldCent%i",i),";p_{T} (GeV/c);Number of d tracks",kNBins,bins);
    fDYieldTOFCent[i] = new TH1F(Form("fDYieldTOFCent%i",i),";p_{T} (GeV/c);Number of d tracks",kNBins,bins);
    fDdcaXYprimaries[i] = new TH2F(Form("fDdcaXYprimaries_%i",i),";p_{T} (GeV/c); DCA_{xy} (cm)",kNBins,bins,160,-0.5f,0.5f);
    fDdcaZprimaries[i] = new TH2F(Form("fDdcaZprimaries_%i",i),";p_{T} (GeV/c); DCA_{z} (cm)",kNBins,bins,160,-0.5f,0.5f);
    fDdcaXYsecondaries[i] = new TH2F(Form("fDdcaXYsecondaries_%i",i),";p_{T} (GeV/c); DCA_{xy} (cm)",kNBins,bins,160,-0.5f,0.5f);
    fDdcaZsecondaries[i] = new TH2F(Form("fDdcaZsecondaries_%i",i),";p_{T} (GeV/c); DCA_{z} (cm)",kNBins,bins,160,-0.5f,0.5f);
    GetOutputList()->Add(fAntiDMCYieldCent[i]);
    GetOutputList()->Add(fAntiDYieldCent[i]);
    GetOutputList()->Add(fAntiDYieldTOFCent[i]);
    GetOutputList()->Add(fDMCYieldCent[i]);
    GetOutputList()->Add(fDYieldCent[i]);
    GetOutputList()->Add(fDYieldTOFCent[i]);
    GetOutputList()->Add(fDdcaXYprimaries[i]);
    GetOutputList()->Add(fDdcaZprimaries[i]);
    GetOutputList()->Add(fDdcaXYsecondaries[i]);
    GetOutputList()->Add(fDdcaZsecondaries[i]);
  }
  
  GetOutputList()->Add(fCentrality);
  GetOutputList()->Add(fAntiDMCYield);
  GetOutputList()->Add(fAntiDYield);
  GetOutputList()->Add(fAntiDYieldTOF);
  GetOutputList()->Add(fDMCYield);
  GetOutputList()->Add(fDYield);
  GetOutputList()->Add(fDYieldTOF);
  GetOutputList()->Add(fPtCorrectionD);
  GetOutputList()->Add(fPtCorrectionAD);

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
  
  if (centrality != fPrevious) {
    fPrevious = centrality;
    if (Skip(centrality)) return kTRUE;
    fCentrality->Fill(centrality);
  }
  int cent = GetCentrality(centrality);
  if (cent < 0) return kTRUE;
  
  if (TMath::Abs(etaMC) < 0.8 && TMath::Abs(yMC) < 0.5) {
    if ((ITSnClusters - ITSnSignal) > 0 && TPCnClusters > 70u && TPCnSignal > 70 && chi2 >= 0.f &&
        chi2 < 6) { // reconstructed (anti) deuterons
      if (IsPrimary) {
        if (pTMC > 0.f) { // deuteron
          fDYield->Fill(TMath::Abs(pTMC));
          fDYieldCent[cent]->Fill(TMath::Abs(pTMC));
          fDdcaXYprimaries[cent]->Fill(TMath::Abs(pT), DCAxy);
          fDdcaZprimaries[cent]->Fill(TMath::Abs(pT), DCAz);
          if (beta > 0.f && beta < 1.f) {
            fDYieldTOF->Fill(TMath::Abs(pTMC));
            fDYieldTOFCent[cent]->Fill(TMath::Abs(pTMC));
          }
        } else {        // anti-deuteron
          fAntiDYield->Fill(TMath::Abs(pTMC));
          fAntiDYieldCent[cent]->Fill(TMath::Abs(pTMC));
          if (beta > 0.f && beta < 1.f) {
            fAntiDYieldTOF->Fill(TMath::Abs(pTMC));
            fAntiDYieldTOFCent[cent]->Fill(TMath::Abs(pTMC));
          }
        }
      } else if(IsSecondaryFromMaterial && pTMC > 0.f) {
        if (pT < 0.8 || beta > 0) {
          fDdcaXYsecondaries[cent]->Fill(TMath::Abs(pT), DCAxy);
          fDdcaZsecondaries[cent]->Fill(TMath::Abs(pT), DCAz);
        }
      }
      if (pTMC > 0) {
        fPtCorrectionD->Fill(TMath::Abs(pT), TMath::Abs(pT) - TMath::Abs(pTMC));
      } else {
        fPtCorrectionAD->Fill(TMath::Abs(pT), TMath::Abs(pT) - TMath::Abs(pTMC));
      }
    }
    if (IsPrimary) {
      if (pTMC > 0.f) {
        fDMCYield->Fill(TMath::Abs(pTMC));
        fDMCYieldCent[cent]->Fill(TMath::Abs(pTMC));
      }
      else {
        fAntiDMCYield->Fill(TMath::Abs(pTMC));
        fAntiDMCYieldCent[cent]->Fill(TMath::Abs(pTMC));
      }
    }
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
  
  
  TCanvas *effCv = new TCanvas();
  effCv->cd();
  
  TEfficiency *effAcc = new TEfficiency(*fDYield,*fDMCYield);
  effAcc->Paint("");
  TGraphAsymmErrors *gr = effAcc->GetPaintedGraph();
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetRangeUser(0.f,1.1f);
  gr->SetName("d");
  gr->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.7);
  gr->SetLineColor(kRed);
  gr->SetMarkerColor(kRed);
  gr->Draw("AP");
  effAcc->SetName("d");
  effAcc->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
  effAcc->SetMarkerStyle(20);
  effAcc->SetMarkerSize(0.7);
  effAcc->SetLineColor(kRed);
  effAcc->SetMarkerColor(kRed);
  effAcc->Draw("same");
  
  
  TEfficiency *effAccTof = new TEfficiency(*fDYieldTOF,*fDMCYield);
  effAccTof->SetName("dTOF");
  effAccTof->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
  effAccTof->SetMarkerStyle(24);
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
  effAccTofAD->SetMarkerStyle(24);
  effAccTofAD->SetMarkerSize(0.7);
  effAccTofAD->SetLineColor(kBlue);
  effAccTofAD->SetMarkerColor(kBlue);
  effAccTofAD->Draw("same");
  
  TLegend *ll = new TLegend(0.61,0.705,0.86,0.89);
  ll->SetFillColor(kWhite);
  ll->SetLineColor(kWhite);
  ll->AddEntry(gr,"d tracking","lp");
  ll->AddEntry(effAccTof,"d tracking + TOF","lp");
  ll->AddEntry(effAccAD,"#bar{d} tracking","lp");
  ll->AddEntry(effAccTofAD,"#bar{d} tracking + TOF","lp");
  ll->Draw();
  
  TFile f("EfficiencyOutput.root","recreate");
  f.cd();
  effCv->Write();
  f.mkdir("Efficiencies");
  f.cd("Efficiencies");
  effAcc->Write();
  effAccTof->Write();
  effAccAD->Write();
  effAccTofAD->Write();
  f.cd();
  const double bins[] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f,8.0f,10.0f
  };
  TString centString[kNCent] = {"0-10%","10-20%","20-40%","40-60%","60-80%"};
  for (int k = 0; k < kNCent; ++k) {
    f.mkdir(Form("Cent%i",k));
    fDdcaXYprimaries[k] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaXYprimaries_%i",k)));
    fDdcaXYsecondaries[k] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaXYsecondaries_%i",k)));
    fDdcaZprimaries[k] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaZprimaries_%i",k)));
    fDdcaZsecondaries[k] = dynamic_cast<TH2F*>(GetOutputList()->FindObject(Form("fDdcaZsecondaries_%i",k)));
    
    f.cd(Form("Cent%i",k));
    fDdcaXYprimaries[k]->Write();
    fDdcaZprimaries[k]->Write();
    fDdcaXYsecondaries[k]->Write();
    fDdcaZsecondaries[k]->Write();
    f.cd();
    f.mkdir(Form("Cent%i/dca_xy",k));
    f.cd(Form("Cent%i/dca_xy",k));
    
    TLatex latex;
    latex.SetTextSize(0.042);
    latex.SetTextAlign(13);
    latex.SetTextFont(12);
    for (int i = 0; i < 10; ++i) {
      TCanvas c(Form("DCAxy_%i",i),Form("DCAxy %4.2f < pT < %4.2f",bins[i],bins[i + 1]));
      c.cd();
      TH1D *hprim = fDdcaXYprimaries[k]->ProjectionY(Form("primxy_%i",i),i + 1, i + 2);
      TH1D *hseco = fDdcaXYsecondaries[k]->ProjectionY(Form("secoxy_%i",i),i + 1, i + 2);
      TH1D sum(*hprim);
      sum.GetXaxis()->CenterTitle();
      sum.GetYaxis()->CenterTitle();
      sum.Add(hseco);
      sum.SetLineWidth(3);
      sum.SetLineColor(kBlue - 2);
      sum.SetName(Form("sumxy_%i",i));
      sum.SetTitle(";DCA_{xy} (cm); Entries");
      sum.Draw();
      hprim->SetLineColor(kBlue);
      hprim->Draw("same");
      hseco->SetLineColor(kRed);
      hseco->Draw("same");
      latex.DrawLatexNDC(0.135,0.85,Form("#splitline{%4.2f GeV/c < p_{T} #leq %4.2f GeV/c}"
                                         "{Centrality %s}",bins[i],bins[i+1],centString[k].Data()));
      c.Write();
    }
    f.cd();
    f.mkdir(Form("Cent%i/Efficiency",k));
    f.cd(Form("Cent%i/Efficiency",k));

    fAntiDMCYieldCent[k] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fAntiDMCYieldCent%i",k)));
    fAntiDYieldCent[k] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fAntiDYieldCent%i",k)));
    fAntiDYieldTOFCent[k] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fAntiDYieldTOFCent%i",k)));
    fDMCYieldCent[k] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fDMCYieldCent%i",k)));
    fDYieldCent[k] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fDYieldCent%i",k)));
    fDYieldTOFCent[k] = dynamic_cast<TH1F*>(GetOutputList()->FindObject(Form("fDYieldTOFCent%i",k)));
    
    TCanvas cv(Form("Cent%i",k),Form("Cent%i",k));
    cv.cd();
    TEfficiency effAcc2(*fDYieldCent[k],*fDMCYieldCent[k]);
    effAcc2.Paint("");
    gr = effAcc2.GetPaintedGraph();
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->CenterTitle();
    gr->GetYaxis()->SetRangeUser(0.f,1.1f);
    gr->SetName("d");
    gr->SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.7);
    gr->SetLineColor(kRed);
    gr->SetMarkerColor(kRed);
    gr->Draw("AP");
    effAcc2.SetName(Form("d%i",k));
    effAcc2.SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
    effAcc2.SetMarkerStyle(20);
    effAcc2.SetMarkerSize(0.7);
    effAcc2.SetLineColor(kRed);
    effAcc2.SetMarkerColor(kRed);
    effAcc2.Draw("same");
    
    TEfficiency effAccTof2(*fDYieldTOFCent[k],*fDMCYieldCent[k]);
    effAccTof2.SetName(Form("dTOF%i",k));
    effAccTof2.SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
    effAccTof2.SetMarkerStyle(24);
    effAccTof2.SetMarkerSize(0.7);
    effAccTof2.SetLineColor(kRed);
    effAccTof2.SetMarkerColor(kRed);
    effAccTof2.Draw("same");
    
    TEfficiency effAccAD2(*fAntiDYieldCent[k],*fAntiDMCYieldCent[k]);
    effAccAD2.SetName(Form("ad%i",k));
    effAccAD2.SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
    effAccAD2.SetMarkerStyle(20);
    effAccAD2.SetMarkerSize(0.7);
    effAccAD2.SetLineColor(kBlue);
    effAccAD2.SetMarkerColor(kBlue);
    effAccAD2.Draw("same");
    
    TEfficiency effAccTofAD2(*fAntiDYieldTOFCent[k],*fAntiDMCYieldCent[k]);
    effAccTofAD2.SetName(Form("adTOF%i",k));
    effAccTofAD2.SetTitle(";p_{T} (GeV/c);Efficiency x Acceptance");
    effAccTofAD2.SetMarkerStyle(24);
    effAccTofAD2.SetMarkerSize(0.7);
    effAccTofAD2.SetLineColor(kBlue);
    effAccTofAD2.SetMarkerColor(kBlue);
    effAccTofAD2.Draw("same");
    
    ll->Draw();
    cv.Write();
    effAccTofAD2.Write();
    effAccAD2.Write();
    effAccTof2.Write();
    effAcc2.Write();
  }
  delete ll;
  
  f.mkdir("PtCorrection");
  fPtCorrectionD = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fPtCorrectionD"));
  fPtCorrectionAD = dynamic_cast<TH2F*>(GetOutputList()->FindObject("fPtCorrectionAD"));
  if (!fPtCorrectionD || !fPtCorrectionAD) {
    f.Close();
    return;
  }
  TF1 *fitF = new TF1("fitF","[0]+[1]*exp([2]*x)",0,6);
  fitF->SetParameters(-0.002,-0.45,-3.);
  fitF->SetParNames("a","b","c");
  f.cd("PtCorrection");
  TObjArray *objAD = new TObjArray();
  fPtCorrectionAD->FitSlicesY(0,0,-1,0,"QNR",objAD);
  objAD->SetName("CorrectionAD");
  ((TH1D*)objAD->At(1))->Fit(fitF);
  objAD->Write();
  fPtCorrectionAD->Write();
  
  TObjArray *objD = new TObjArray();
  fPtCorrectionD->FitSlicesY(0,0,-1,0,"QNR",objD);
  objD->SetName("CorrectionD");
  ((TH1D*)objD->At(1))->Fit(fitF);
  objD->Write();
  fPtCorrectionD->Write();
  f.cd();
  fCentrality = dynamic_cast<TH1F*>(GetOutputList()->FindObject("fCentrality"));
  if (fCentrality) fCentrality->Write();
  f.Close();

}
