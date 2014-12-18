//
//  MakeFinalResults.C
//
//  Created by Maximiliano Puccio on 10/12/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#include <Riostream.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TList.h>
#include <TFractionFitter.h>
#include <TVirtualFitter.h>

#pragma mark Helper functions declaration
enum { kDeuteron, kAntideuteron };
Float_t  CorrectForMaterial(TH1F* data, TObjArray &obj, TFile* output); //!< Function used to return the primary fraction in a pT bin
void     CorrectForEfficiency(TH1F* rawcounts, TEfficiency* eff, TEfficiency* effTOF);
Double_t CrystalBall(Double_t *x_, Double_t *p_);
Double_t ExpBkg(Double_t *x, Double_t *p);
Double_t ExpCB(Double_t *x_, Double_t *p);
Double_t ExpExpBkg(Double_t *x, Double_t *p);
Double_t ExpExpGaus(Double_t *x, Double_t *p);
Double_t GausSgl(Double_t *x, Double_t *p);
TH1F*    HistoFromFunction(TString title,Double_t (*func)(Double_t*,Double_t*),Double_t *par);
TH1F*    HistoFromFunction(TString title,TF1 *func);
TF1*     TOFfitter(TH1F* data, TH1F* &back, TH1F* &sig, int iCent, int iPt, int k = kDeuteron);

#pragma mark Global information
const TString kMCFile = "EfficiencyOutput.root";  
const TString kDataFile = "AODSelector.root";     
const TString kOutputFile = "FinalResults.root";
const double kBins[] = {
  0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
  1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
  3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f
};
const int kTPCsize = 3;
const int kTPCTOFsize = 23;
const int kNBins = kTPCsize + kTPCTOFsize;
const int kNBinsSecondaries = 8; // TODO: automatic counting of the bins below 1.2 GeV/c
const int kNCentralities = 4;
const float kEfficiencyTolerance = 0.05f;
const int kRebin = 4;
const TString kTitles[kNCentralities] = {"Centrality 0-10%","Centrality 10-20%","Centrality 20-40%","Centrality 40-60%"};//,"Centrality 60-80%"};
const Color_t kColors[kNCentralities] = {kBlue-4,kGreen-3,kOrange,kOrange+8};//,kRed};
const int kNTPCSecondaries = ((kNBinsSecondaries - kTPCsize >= 0)
                              ? kTPCsize : kNBinsSecondaries);
const int kNTOFSecondaries = (kNBinsSecondaries - kNTPCSecondaries);

//__________________________________________________________________________________________________
void MakeFinalResults() {
  
#pragma mark Opening files
  TFile *fDATA = TFile::Open(kDataFile.Data());
  TFile *fDATAH = TFile::Open("mpuccio_Flowd.root");
  TFile *fMC = TFile::Open(kMCFile.Data());
  TFile *fOUTPUT = new TFile(kOutputFile.Data(),"recreate");
  if (!fDATA || !fMC || !fOUTPUT) {
    cout << "Missing files." << endl;
    return;
  }
  fOUTPUT->cd();
  fOUTPUT->mkdir("Deuteron");
  fOUTPUT->mkdir("Deuteron/Spectra");
  fOUTPUT->mkdir("Deuteron/Counts");
  fOUTPUT->mkdir("Deuteron/Efficiency");
  fOUTPUT->mkdir("Deuteron/Summary");
  fOUTPUT->mkdir("Deuteron/Summary/Fit");
  fOUTPUT->mkdir("Deuteron/Summary/SB");
  fOUTPUT->mkdir("Deuteron/SB");
  fOUTPUT->mkdir("Deuteron/Fit");
  fOUTPUT->mkdir("Deuteron/Fractions");
  fOUTPUT->mkdir("Deuteron/Fractions/Fit");
  fOUTPUT->mkdir("Deuteron/Fractions/debug");
  fOUTPUT->mkdir("AntiDeuteron");
  fOUTPUT->mkdir("AntiDeuteron/Spectra");
  fOUTPUT->mkdir("AntiDeuteron/Counts");
  fOUTPUT->mkdir("AntiDeuteron/Efficiency");
  fOUTPUT->mkdir("AntiDeuteron/Summary");
  fOUTPUT->mkdir("AntiDeuteron/Summary/Fit");
  fOUTPUT->mkdir("AntiDeuteron/Summary/SB");
  fOUTPUT->mkdir("AntiDeuteron/SB");
  fOUTPUT->mkdir("AntiDeuteron/Fit");


#pragma mark Taking all Monte Carlo information
  TEfficiency *dEff     = (TEfficiency*)fMC->Get("Efficiencies/d");
  TEfficiency *dEffTOF  = (TEfficiency*)fMC->Get("Efficiencies/dTOF");
  TEfficiency *adEff    = (TEfficiency*)fMC->Get("Efficiencies/ad");
  TEfficiency *adEffTOF = (TEfficiency*)fMC->Get("Efficiencies/adTOF");
  if (!dEff || !dEffTOF || !adEff || !adEffTOF) {
    cout << "Missing efficiencies." << endl;
    return;
  }

  TH1F hMCPrimaries[kNBinsSecondaries * kNCentralities];
  TH1F hMCSecondaries[kNBinsSecondaries * kNCentralities];
  for (int iC = 0; iC < kNCentralities; ++iC) {
    TH2F *hMC2Dp = (TH2F*)fMC->Get(Form("Cent%i/fDdcaXYprimaries_%i",iC,iC));
    TH2F *hMC2Ds = (TH2F*)fMC->Get(Form("Cent%i/fDdcaXYsecondaries_%i",iC,iC));
    if (!hMC2Dp || !hMC2Ds) {
      cout << "Missing MC template." << endl;
      return;
    }
    
    for (int iB = 0; iB < kNBinsSecondaries; ++iB) {
      TH1D* hProjP = hMC2Dp->ProjectionY(Form("hMCProjP%i_%i",iC,iB), iB + 1, iB + 1);
      TH1D* hProjS = hMC2Ds->ProjectionY(Form("hMCProjS%i_%i",iC,iB), iB + 1, iB + 1);
      TH1F &hPrim = hMCPrimaries[iC * kNBinsSecondaries + iB];
      TH1F &hSec = hMCSecondaries[iC * kNBinsSecondaries + iB];
      hPrim.SetNameTitle(Form("hMCPrimaries%i_%i",iC,iB),
                         Form("%4.2f < p_{T} #leq %4.2f; DCA_{xy} (cm); Entries",
                              kBins[iB],kBins[iB + 1]));
      hSec.SetNameTitle(Form("hMCSecondaries%i_%i",iC,iB),
                        Form("%4.2f < p_{T} #leq %4.2f; DCA_{xy} (cm); Entries",
                             kBins[iB],kBins[iB + 1]));
      hPrim.SetBins(hProjP->GetNbinsX(), hProjP->GetXaxis()->GetXmin(),
                    hProjP->GetXaxis()->GetXmax());
      hSec.SetBins(hProjS->GetNbinsX(), hProjS->GetXaxis()->GetXmin(),
                   hProjS->GetXaxis()->GetXmax());
      hPrim.Add(hProjP);
      hSec.Add(hProjS);
      hPrim.Rebin(kRebin);
      hSec.Rebin(kRebin);
      delete hProjP;
      delete hProjS;
    }
  }
  
#pragma mark Common variables
  TH1F hRawCounts[kNCentralities];
  TH1F hSpectra[kNCentralities];
  for (int iC = 0; iC < kNCentralities; ++iC) {
    hRawCounts[iC].SetBins(kNBins, kBins);
    hSpectra[iC].SetBins(kNBins, kBins);
  }
  TList * l = (TList*)fDATAH->Get("mpuccio_Flowd");
  TH1F* hEventCount = (TH1F*)l->FindObject("fHistCentralityClass10");
  if (!hEventCount) {
    cout << "Missing event counter" << endl;
    return;
  }
  double kNEvents[kNCentralities];
  double kNorm[kNCentralities];

  kNEvents[0] = hEventCount->GetBinContent(1);
  kNEvents[1] = hEventCount->GetBinContent(2);
  kNEvents[2] = hEventCount->GetBinContent(3) + hEventCount->GetBinContent(4);
  kNEvents[3] = hEventCount->GetBinContent(5) + hEventCount->GetBinContent(6);
//  kNEvents[4] = hEventCount->GetBinContent(6) + hEventCount->GetBinContent(7);

  for (int i = 0; i < kNCentralities; ++i)
    kNorm[i] = 1. / (kNEvents[i] * TMath::TwoPi());
  
  TH1F hFractions[kNCentralities];
  for (int iC = 0; iC < kNCentralities; ++iC) {
    hFractions[iC].SetName(Form("hFractions%i",iC));
    hFractions[iC].SetTitle(Form("%s;p_{T} (GeV/c); Primary fractions",kTitles[iC].Data()));
    hFractions[iC].SetBins(kNBins, kBins);
  }
  
#pragma mark Analysing deuterons
  for (int iC = 0; iC < kNCentralities; ++iC) {
    hRawCounts[iC].Reset();
    hSpectra[iC].Reset();
    hRawCounts[iC].SetName(Form("hRawCountsD%i",iC));
    hRawCounts[iC].SetTitle(Form("%s;p_{T} (GeV/c); d raw counts",kTitles[iC].Data()));
    hSpectra[iC].SetName(Form("hSpectraD%i",iC));
    hSpectra[iC].SetTitle(";p_{T} (GeV/c);1/N_{ev}1/(2#pip_{T})d^{2}N/(dp_{T}dy) (GeV/c)^{-2}");
    TH1F* hTPCcounts = (TH1F*)fDATA->Get(Form("fdEdxTPCSignalCounts%i",iC));
    if (!hTPCcounts && kTPCsize > 0) {
      cout << "Missing " << Form("fdEdxTPCSignalCounts%i",iC) << endl;
      return;
    }
    TCanvas *summaryCv[3], *sbCv[3];
    for (int i = 0; i < 3; ++i) {
      summaryCv[i] = new TCanvas(Form("keynoted%i_%i",iC,i));
      summaryCv[i]->Divide(3,3);
      sbCv[i] = new TCanvas(Form("keynotedSB%i_%i",iC,i));
      sbCv[i]->Divide(3,3);
    }
    
    // DCA preparations
    TObjArray obj(2);
    obj.SetOwner(kFALSE);
    
    //
    for (int iB = 0; iB < kNBins; ++iB) {
      if (iB < kTPCsize) {
        if (iB < kNBinsSecondaries) {
          TH1F *hd = (TH1F*)fDATA->Get(Form("DCA%i/dcaxy_%i",iC,iB));
          if (!hd) {
            cout << "Missing " << Form("DCA%i/dcaxy_%i",iC,iB) << endl;
            return;
          }
          Float_t prim = 1.f;
          hd->Rebin(kRebin);
          obj.Add(&hMCPrimaries[kNBinsSecondaries * iC + iB]);
          obj.Add(&hMCSecondaries[kNBinsSecondaries * iC + iB]);
          prim = CorrectForMaterial(hd, obj, fOUTPUT);
          fOUTPUT->mkdir(Form("Deuteron/Fractions/debug/%i_%i",iC,iB));
          fOUTPUT->cd(Form("Deuteron/Fractions/debug/%i_%i",iC,iB));
          hd->Write();
          hMCPrimaries[kNBinsSecondaries * iC + iB].Write();
          hMCSecondaries[kNBinsSecondaries * iC + iB].Write();
          obj.Remove(&hMCPrimaries[kNBinsSecondaries * iC + iB]);
          obj.Remove(&hMCSecondaries[kNBinsSecondaries * iC + iB]);
          cout << "CAZZO " << prim << endl;
          hFractions[iC].SetBinContent(iB + 1,prim);
          cout<< "CAZZO " << prim << endl;
          hRawCounts[iC].SetBinContent(iB + 1, prim * hTPCcounts->GetBinContent(iB + 1));
        } else
          hRawCounts[iC].SetBinContent(iB + 1, hTPCcounts->GetBinContent(iB + 1));
      } else {
        TH1F *data = (TH1F*)fDATA->Get(Form("Signal%i/fSignalD%i_%i",iC,iC,iB - kTPCsize));
        if (!data) {
          cout << "Missing " << Form("Signal%i/fSignalD%i_%i",iC,iC,iB - kTPCsize);
          cout << " from DATA file\n";
          return;
        }
        const int n = (iB - kTPCsize) / 9;
        int shift = -1;
        if (n == 1) shift = 9 - 1;
        else if (n == 2) shift = 18 -1;
        TH1F *back, *sigl;
        summaryCv[n]->cd(iB - kTPCsize - shift);
        TF1 *tpl = TOFfitter(data, back, sigl, iC, iB);
        data->DrawCopy("e");
        back->DrawCopy("lsame");
        Double_t *pars = tpl->GetParameters();
        TH1F *hTPL = HistoFromFunction("tpl", tpl);
        hTPL->Add(back,-1.f);
        Float_t prim = 1.f;
        if (iB < kNBinsSecondaries) {
          TH2F *hDCASig = (TH2F*)fDATA->Get(Form("DCASignal/fDCASignal%i_%i",iC,iB - kTPCsize));
          if (!hDCASig) {
            cout << "Missing " << Form("DCASignal/fDCASignal%i_%i",iC,iB) << endl;
            return;
          }
          TH1D *hP = hDCASig->ProjectionY(Form("hDCATOFp%i_%i",iC,iB));/*,
                                          hTPL->FindBin(pars[1] - 3.f * pars[2]),
                                          hTPL->FindBin(pars[1] + 3.f * pars[2]));*/
          TH1F hd(Form("hDCATOF%i_%i",iC,iB),";DCA_{xy} (cm);Entries",hP->GetNbinsX(),
                  hP->GetXaxis()->GetXmin(),hP->GetXaxis()->GetXmax());
          hd.Add(hP);
          hd.Rebin(4);
          obj.Add(&hMCPrimaries[kNBinsSecondaries * iC + iB]);
          obj.Add(&hMCSecondaries[kNBinsSecondaries * iC + iB]);
          prim = CorrectForMaterial(&hd, obj, fOUTPUT);
          fOUTPUT->mkdir(Form("Deuteron/Fractions/debug/%i_%i",iC,iB));
          fOUTPUT->cd(Form("Deuteron/Fractions/debug/%i_%i",iC,iB));
          hd.Write();
          hMCPrimaries[kNBinsSecondaries * iC + iB].Write();
          hMCSecondaries[kNBinsSecondaries * iC + iB].Write();
          obj.Remove(&hMCPrimaries[kNBinsSecondaries * iC + iB]);
          obj.Remove(&hMCSecondaries[kNBinsSecondaries * iC + iB]);

          delete hP;
        }
        hFractions[iC].SetBinContent(iB + 1, prim);
        hRawCounts[iC].SetBinContent(iB + 1, prim * hTPL->Integral(hTPL->FindBin(pars[1] - 3.f * pars[2]),
                                                                   hTPL->FindBin(pars[1] + 3.f * pars[2])));
        sbCv[n]->cd(iB - kTPCsize - shift);
        sigl->Divide(tpl);
        back->Divide(tpl);
        sigl->DrawCopy("l");
        back->DrawCopy("lsame");
        delete hTPL;
        delete tpl;
        delete back;
        delete sigl;
      }
    }
  
    for (int i = 0; i < 3; ++i) {
      fOUTPUT->cd("Deuteron/Summary/Fit");
      summaryCv[i]->Write();
      delete summaryCv[i];
      fOUTPUT->cd("Deuteron/Summary/SB");
      sbCv[i]->Write();
      delete sbCv[i];
      fOUTPUT->cd();
    }
    fOUTPUT->cd("Deuteron/Counts");
    hSpectra[iC].Add(hRawCounts);
    hRawCounts[iC].Write();
    
    fOUTPUT->cd("Deuteron/Fractions");
    hFractions[iC].Write();
    
    fOUTPUT->cd("Deuteron/Spectra");
    TEfficiency *eff = (TEfficiency*)fMC->Get(Form("Cent%i/Efficiency/d%i",iC,iC));
    TEfficiency *effTOF = (TEfficiency*)fMC->Get(Form("Cent%i/Efficiency/dTOF%i",iC,iC));
    if (!eff || !effTOF) {
      cout << "Missing efficiencies!" << endl;
      return;
    }
    CorrectForEfficiency(&hSpectra[iC], eff, effTOF);
    hSpectra[iC].Scale(kNorm[iC]);
    for (int i = 1; i <= hSpectra[iC].GetNbinsX(); ++i) {
      hSpectra[iC].SetBinContent(i, hSpectra[iC].GetBinContent(i) /
                                 (hSpectra[iC].GetBinCenter(i) * hSpectra[iC].GetBinWidth(i)));
    }
    hSpectra[iC].Write();
  }

  fOUTPUT->cd("Deuteron/Efficiency");
  TGraphAsymmErrors* grDeff = adEff->CreateGraph();
  grDeff->SetName("dEffTPC");
  grDeff->Write();
  delete grDeff;
  grDeff = dEffTOF->CreateGraph();
  grDeff->SetName("dEffTPCTOF");
  grDeff->Write();
  delete grDeff;

#pragma mark Analysing anti deuterons
  for (int iC = 0; iC < kNCentralities; ++iC) {
    hRawCounts[iC].Reset();
    hSpectra[iC].Reset();
    hRawCounts[iC].SetName(Form("hRawCountsAD%i",iC));
    hRawCounts[iC].SetTitle(Form("%s;p_{T} (GeV/c); #bar{d} raw counts",kTitles[iC].Data()));
    hSpectra[iC].SetName(Form("hSpectraAD%i",iC));
    hSpectra[iC].SetTitle(";p_{T} (GeV/c);1/N_{ev}1/(2#pip_{T})d^{2}N/(dp_{T}dy) (GeV/c)^{-2}");
    TH1F* hTPCcounts = (TH1F*)fDATA->Get(Form("fdEdxTPCSignalCountsAD%i",iC));
    if (!hTPCcounts && kTPCsize > 0) {
      cout << "Missing " << Form("fdEdxTPCSignalCountsAD%i",iC) << endl;
      return;
    }
    TCanvas *summaryCv[3], *sbCv[3];
    for (int i = 0; i < 3; ++i) {
      summaryCv[i] = new TCanvas(Form("keynotead%i_%i",iC,i));
      summaryCv[i]->Divide(3,3);
      sbCv[i] = new TCanvas(Form("keynoteadSB%i_%i",iC,i));
      sbCv[i]->Divide(3,3);
    }
    for (int iB = 0; iB < kNBins; ++iB) {
      if (iB < kTPCsize) {
        hRawCounts[iC].SetBinContent(iB + 1, hTPCcounts->GetBinContent(iB + 1));
      } else {
        if (iB < kTPCsize) {
          hRawCounts[iC].SetBinContent(iB + 1, hTPCcounts->GetBinContent(iB + 1));
        } else {
          TH1F *data = (TH1F*)fDATA->Get(Form("Signal%i/fSignalAD%i_%i",iC,iC,iB - kTPCsize));
          if (!data) {
            cout << "Missing " << Form("Signal%i/fSignalAD%i_%i",iC,iC,iB - kTPCsize);
            cout << " from DATA file\n";
            return;
          }
          const int n = (iB - kTPCsize) / 9;
          int shift = -1;
          if (n == 1) shift = 9 - 1;
          else if (n == 2) shift = 18 -1;
          TH1F *back, *sigl;
          summaryCv[n]->cd(iB - kTPCsize - shift);
          TF1 *tpl = TOFfitter(data, back, sigl, iC, iB, kAntideuteron);
          TH1F *hTPL = HistoFromFunction("tpl", tpl);
          data->DrawCopy("e");
          back->DrawCopy("lsame");
          hTPL->Add(back,-1.f);
          Double_t *pars = tpl->GetParameters();
          hRawCounts[iC].SetBinContent(iB, hTPL->Integral(hTPL->FindBin(pars[1] - 3.f * pars[2]),
                                                          hTPL->FindBin(pars[1] + 3.f * pars[2])));
          sbCv[n]->cd(iB - kTPCsize - shift);
          sigl->Divide(tpl);
          back->Divide(tpl);
          sigl->DrawCopy("l");
          back->DrawCopy("lsame");
          delete hTPL;
          delete tpl;
          delete back;
          delete sigl;
        }

      }
    }
    for (int i = 0; i < 3; ++i) {
      fOUTPUT->cd("AntiDeuteron/Summary/Fit");
      summaryCv[i]->Write();
      delete summaryCv[i];
      fOUTPUT->cd("AntiDeuteron/Summary/SB");
      sbCv[i]->Write();
      delete sbCv[i];
      fOUTPUT->cd();
    }
    fOUTPUT->cd("AntiDeuteron/Counts");
    hSpectra[iC].Add(hRawCounts);
    hRawCounts[iC].Write();
    
    fOUTPUT->cd("AntiDeuteron/Spectra");
    TEfficiency *eff = (TEfficiency*)fMC->Get(Form("Cent%i/Efficiency/ad%i",iC,iC));
    TEfficiency *effTOF = (TEfficiency*)fMC->Get(Form("Cent%i/Efficiency/adTOF%i",iC,iC));
    if (!eff || !effTOF) {
      cout << "Missing efficiencies!" << endl;
      return;
    }
    CorrectForEfficiency(&hSpectra[iC], eff, effTOF);
    hSpectra[iC].Scale(kNorm[iC]);
    for (int i = 1; i <= hSpectra[iC].GetNbinsX(); ++i) {
      hSpectra[iC].SetBinContent(i, hSpectra[iC].GetBinContent(i) /
                                 (hSpectra[iC].GetBinCenter(i) * hSpectra[iC].GetBinWidth(i)));
    }
    hSpectra[iC].Write();
  }
  
  fOUTPUT->cd("AntiDeuteron/Efficiency");
  TGraphAsymmErrors* grADeff = adEff->CreateGraph();
  grADeff->SetName("adEffTPC");
  grADeff->Write();
  delete grADeff;
  grADeff = adEffTOF->CreateGraph();
  grADeff->SetName("adEffTPCTOF");
  grADeff->Write();
  delete grADeff;
  fOUTPUT->cd();

#pragma mark Closing files
  fDATA->Close();
  fMC->Close();
  fOUTPUT->Close();
  delete fOUTPUT;
}

#pragma mark Helper function definition
//__________________________________________________________________________________________________
Float_t CorrectForMaterial(TH1F* hd, TObjArray &obj, TFile* output) {
  TVirtualFitter::SetMaxIterations(10000000);

  Int_t binLow = ((TH1F*)obj[0])->FindBin(-0.5);
  Int_t binUp = ((TH1F*)obj[0])->FindBin(0.5);
  TFractionFitter fitter(hd,&obj);
  fitter.SetRangeX(binLow,binUp);
  fitter.Constrain(0,0.,1.);
  fitter.Constrain(1,0.,1.);
  Int_t result = fitter.Fit();
  Double_t yieldSec = 0., yieldPri = 1., error = 0.;
  if (result == 0) {
    TH1F* hp = (TH1F*)fitter.GetMCPrediction(0);
    TH1F* hs = (TH1F*)fitter.GetMCPrediction(1);
    fitter.GetResult(0, yieldPri, error);
    fitter.GetResult(1, yieldSec, error);
    TH1F* hfit = (TH1F*)fitter.GetPlot();
    Float_t dataIntegral = hfit->Integral();
    hfit->SetLineColor(kGreen + 1);
    hfit->SetLineWidth(3);
    hs->Scale(1. / hs->Integral());
    hp->Scale(1. / hp->Integral());
    hs->Scale(yieldSec * dataIntegral);
    hp->Scale(yieldPri * dataIntegral);
    TCanvas cv(hd->GetName());
    cv.cd();
    hd->SetMarkerStyle(20);
    hd->SetMarkerSize(0.5);
    hd->SetMarkerColor(kBlack);
    hd->Draw("e");
    
    hfit->DrawCopy("same");
    hs->SetLineColor(kRed);
    hp->SetLineColor(kBlue);
    hs->DrawCopy("same");
    hp->DrawCopy("same");
    output->cd("Deuteron/Fractions/Fit");
    cv.Write();
  }
//  else {
//    hd->Rebin(2);
//    ((TH1F*)obj[0])->Rebin(2);
//    ((TH1F*)obj[1])->Rebin(2);
//    if (hd->GetNbinsX() < 4) return 1.f;
//    return CorrectForMaterial(hd, obj, output);
//  }
    cout << "LOREM IPSUM " << yieldPri << " " << yieldSec << " " << obj.GetEntries() << " " << hd->GetName() << endl;
  return yieldPri;
  
}

//__________________________________________________________________________________________________
void CorrectForEfficiency(TH1F* rawcounts, TEfficiency* eff, TEfficiency* effTOF) {
  for (int iB = 1; iB < kNBins; ++iB) {
    float e = 0.f;
    if (iB <= kTPCsize) e = eff->GetEfficiency(iB);
    else e = effTOF->GetEfficiency(iB);
    if (e < kEfficiencyTolerance) rawcounts->SetBinContent(iB, 0.f);
    else rawcounts->SetBinContent(iB, rawcounts->GetBinContent(iB) / e);
  }
}

//__________________________________________________________________________________________________
Double_t CrystalBall(Double_t *x_, Double_t *p_)
{
  Double_t &m = x_[0];
  Double_t &norm = p_[0];
  Double_t &m0 = p_[1];
  Double_t &sigma = p_[2];
  Double_t &alpha = p_[3];
  Double_t &n = p_[4];
  
  Double_t t = (m - m0)/sigma;
  if (alpha > 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)alpha);
  
  if (t >= -absAlpha) {
    return norm * exp(-0.5 * t * t);
  }
  else {
    Double_t a =  TMath::Power(n / absAlpha,n) * exp(-0.5 * absAlpha * absAlpha);
    Double_t b = n / absAlpha - absAlpha;
    return norm * a / TMath::Power(b - t, n);
  }
}

//__________________________________________________________________________________________________
Double_t ExpBkg(Double_t *x, Double_t *p) {
  return p[0] * TMath::Exp(p[1] * x[0]);
}

//__________________________________________________________________________________________________
Double_t ExpCB(Double_t *x_, Double_t *p) {
  Double_t &x = x_[0];
  return CrystalBall(x_,p) + p[5] * TMath::Exp(p[6] * x);
}

//__________________________________________________________________________________________________
Double_t ExpExpBkg(Double_t *x, Double_t *p) {
  return p[0] * TMath::Exp(p[1] * x[0]) + p[2] * TMath::Exp(p[3] * x[0]);
}

//__________________________________________________________________________________________________
Double_t ExpExpGaus(Double_t *x, Double_t *p) {
  return GausSgl(x,p) + ExpBkg(x,&p[3]) + ExpBkg(x,&p[5]);
}

//__________________________________________________________________________________________________
Double_t GausSgl(Double_t *x, Double_t *p) {
  return p[0] * TMath::Gaus(x[0],p[1],p[2]);
}

//__________________________________________________________________________________________________
TH1F * HistoFromFunction(TString title,Double_t (*func)(Double_t*,Double_t*),Double_t *par) {
  TH1F * a = new TH1F(title.Data(),title.Data(),100,-2.4,2.4);
  Double_t x[1];
  for (int i = 1; i <= 100; ++i) {
    x[0] = a->GetBinCenter(i);
    a->SetBinContent(i,func(x,par));
  }
  return a;
}

//__________________________________________________________________________________________________
TH1F * HistoFromFunction(TString title,TF1 *func) {
  TH1F * a = new TH1F(title.Data(),title.Data(),100,-2.4,2.4);
  Double_t x;
  for (int i = 1; i <= 100; ++i) {
    x = a->GetBinCenter(i);
    a->SetBinContent(i,func->Eval(x));
  }
  return a;
}

//__________________________________________________________________________________________________
TF1* TOFfitter(TH1F* data, TH1F* &back, TH1F* &sigl, int iCent, int iPt, int k) {
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.5);
  data->SetMarkerColor(kBlue);
  TF1 *tpl;
  if (kBins[iPt] < 1.6f) {
    tpl = new TF1(Form("tpl%i_%i",iCent,iPt),ExpCB,-2.4,2.4,7);
    tpl->SetParameters((iCent == 0) ? 20000 : 2500, 0.f, 2.86394e-01, 1, 2, 42,-0.6);
    if (iCent > 2 && k == kDeuteron && kBins[iPt] < 1.) {
      tpl->SetParameters(500, 0.f, 1.86394e-01, 1, 20, 0.84,-0.6);
    }
    if (k == kAntideuteron) {
      if (kBins[iPt] < 1.) {
        if (iCent == 0) {
          tpl->SetParameters((kBins[iPt] < 0.9) ? 300 : 800, 0.f, 1.86394e-01, 1, 20, 0.84,-0.6);
        } else if (iCent == 2) {
          tpl->SetParameters((kBins[iPt] < 0.9) ? 140 : 800, 0.09f, 1.86394e-01, 1, 20, 0.84,-0.6);
        } else {
          if (kBins[iPt] >= 0.9)
            tpl->SetParameters(800, 0.09f, 1.86394e-01, 1, 20, 0.84,-0.6);
          else
            tpl->SetParameters((iCent == 1) ? 120 : 80, 0.09f, 0.095, 0.7, 12, 0.84,-0.3);
        }
      }
    }
    tpl->SetParNames("N","#mu","#sigma","#alpha","n","A","#tau");
  } else {
    tpl = new TF1(Form("tpl%i_%i",iCent,iPt),ExpExpGaus,-2.4,2.4,7);
    tpl->SetParameters((kBins[iPt] > 4.f) ? 500 : 5000, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);
    if (iCent >= 1) {
      tpl->SetParameters((kBins[iPt] > 4.f) ? 40 : 500, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);
    }
    if (iCent > 2) {
      int height = 600;
      if (kBins[iPt] > 3.2f) height = 20;
      else if (kBins[iPt] > 1.8f) height = 300;
      tpl->SetParameters(height, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);
    }
    if (k == kAntideuteron) {
      if (iCent == 3) {
        tpl->SetParameters(700, 5.66003e-02, 0.12, 4.28095e-01,-4.04126e-01, 6.51728e+02,-3.00000e-01);
        if (kBins[iPt] > 3.2f)
          tpl->SetParameters(20, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);
        else if (kBins[iPt] > 1.8f)
          tpl->SetParameters(300, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);

      }
    }
    tpl->SetParNames("N","#mu","#sigma","A_{1}","#tau_{1}","A_{2}","#tau_{2}");
  }
  tpl->SetParLimits(1,-0.2,0.2);
  tpl->SetParLimits(2,0.,0.6);
  tpl->SetNpx(300);
  tpl->SetDrawOption("e");
  data->Fit(tpl,"LRQ","",-2,2);
  data->Fit(tpl,"LRQ","",-2,2);
  Double_t *parameters = tpl->GetParameters();

  if (kBins[iPt] < 1.6f) {
    back = HistoFromFunction(Form("background%i_%i",iCent,iPt),ExpBkg,&parameters[5]);
    sigl = HistoFromFunction(Form("signl%i_%i",iCent,iPt),CrystalBall,parameters);
  } else {
    back = HistoFromFunction(Form("background%i_%i",iCent,iPt),ExpExpBkg,&parameters[3]);
    sigl = HistoFromFunction(Form("signl%i_%i",iCent,iPt),GausSgl,parameters);
  }
  back->SetDirectory(0);
  sigl->SetDirectory(0);
  back->SetLineColor(kBlue);
  back->SetNameTitle(Form("back%i_%i",iCent,iPt),
                            Form("%4.1f #leq p_{T} < %4.1f ;m_{TOF}^{2}-m_{PDG}^{2};S/B",
                                 kBins[iPt],kBins[iPt+1]));
  back->GetYaxis()->SetRangeUser(0,1.05);
  sigl->SetNameTitle(Form("sigl%i_%i",iCent,iPt),
                            Form("%4.1f #leq p_{T} < %4.1f ;m_{TOF}^{2}-m_{PDG}^{2};S/B",
                                 kBins[iPt],kBins[iPt+1]));
  sigl->SetLineColor(kRed);
  return tpl;
}
