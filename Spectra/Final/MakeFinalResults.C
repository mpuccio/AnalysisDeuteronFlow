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
#include "RooGlobalFunc.h"
#include <RooDataHist.h>
#include <RooGExpModel.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooExponential.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooDecay.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooCBShape.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooArgSet.h>
#include <TLatex.h>

using TMath::Sqrt2;
using TMath::Sqrt;
using TMath::Exp;
using namespace RooFit;
#define sq(x) x * x

#pragma mark Helper functions declaration
enum { kDeuteron, kAntideuteron };
Float_t  CorrectForMaterial(TH1F* data, TObjArray &obj, TFile* output); //!< Function used to return the primary fraction in a pT bin
void     CorrectForEfficiency(TH1F* rawcounts, TF1* eff, TF1* effTOF);
void     CorrectForEfficiency(TH1F* rawcounts, TEfficiency* eff, TEfficiency* effTOF);
Double_t FitTemplateLowPt(Double_t *x_, Double_t *p);
Double_t FitTemplateHighPt(Double_t *x_, Double_t *p);
Double_t ExpoPdf(Double_t *x_, Double_t *p);
Double_t DoubleExpoPdf(Double_t *x_, Double_t *p);
Double_t TailedGausPdf(Double_t *x_, Double_t *p);

TH1F*    HistoFromFunction(TString title,Double_t (*func)(Double_t*,Double_t*),Double_t *par);
TH1F*    HistoFromFunction(TString title,TF1 *func);

#pragma mark Global information
const TString kMCFile = "EfficiencyOutput.root";
const TString kDataFile = "AODSelector.root";
const TString kOutputFile = "FinalResults.root";
const double kBins[] = {
  0.4f,0.5f,0.6f,
  0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
  1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
  3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f
};
const int kTPCsize = 3;
const int kTPCTOFsize = 23;
const int kNBins = kTPCsize + kTPCTOFsize;
const int kNBinsSecondaries = 8; // TODO: automatic counting of the bins below 1.2 GeV/c
const int kNCentralities = 4;
const float kEfficiencyTolerance = 0.05f;
const float kExpoPdfShift = -2.f;
const int kNBinsLowPt = 4;
const int kRebin = 8;
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
  TFile *fEff = TFile::Open("FittedEfficiencies.root");
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
  

  
  
  //TH1F *dat = (TH1F*)fDATA->Get("Signal0/fSignalAD0_2");

#pragma mark Setup fit framework
  RooRealVar m("dm2","m^{2} - m^{2}_{PDG} (GeV/c^{2})^{2}",-2.,2.5);
  m.setBins(10000,"cache");
  RooRealVar s1("sigma1","#sigma1",0.08,0.2);
  RooRealVar s2("sigma2","#sigma2",0.1,1.);
  RooRealVar mu("mu","#mu",-0.5,0.5);
  RooRealVar tau1("tau1","#tau_{1}",-6.5,-1.5);
  RooRealVar tau2("tau2","#tau_{2}",-1.5,0.);
  RooRealVar cbkg("cbkg","Background coefficient",0.,1.);
  RooRealVar csig("csig","Signal coefficient",0.,1.);
  RooRealVar fbkg("fbkg","Background fraction",5000.,0.,500000.);
  RooRealVar fsig("fsig","Signal Fraction",5000.,0.,500000.);
//  RooCBShape sig("signal", "signal", m, mu, s, alpha, n);
  RooGaussian sig1("sig1", "sig1", m, mu, s1);
  RooGaussian sig2("sig2", "sig2", m, mu, s2);
  RooAddPdf sig0("signal","signal",RooArgList(sig1,sig2),csig);

  RooExponential bkg1("background1","background1",m,tau1);
  RooExponential bkg2("background2","background2",m,tau2);
  RooAddPdf bkg0("background","background",RooArgList(bkg1,bkg2),cbkg);
  RooAddPdf model1("model1","model1",RooArgList(sig0,bkg2),RooArgList(fsig,fbkg));
  RooAddPdf model2("model2","model2",RooArgList(sig2,bkg0),RooArgList(fsig,fbkg));
  
#pragma mark Taking all Monte Carlo information
  TF1 *dEff[4],*adEff[4],*dEffTOF[4],*adEffTOF[4];
  for (int iC = 0; iC < 4; ++iC) {
    TGraphAsymmErrors* grd = (TGraphAsymmErrors*)fEff->Get(Form("d%i",iC));
    TGraphAsymmErrors* grad = (TGraphAsymmErrors*)fEff->Get(Form("ad%i",iC));
    TGraphAsymmErrors* grdTOF = (TGraphAsymmErrors*)fEff->Get(Form("dTOF%i",iC));
    TGraphAsymmErrors* gradTOF = (TGraphAsymmErrors*)fEff->Get(Form("adTOF%i",iC));
    if (!grd || !grad || !grdTOF || !gradTOF) {
      cout << "Missing efficiencies." << endl;
      return;
    }
    dEff[iC] = static_cast<TF1*>(grd->GetListOfFunctions()->FindObject("pol1"));
    adEff[iC] = static_cast<TF1*>(grad->GetListOfFunctions()->FindObject("pol1"));
    dEffTOF[iC] = static_cast<TF1*>(grdTOF->GetListOfFunctions()->FindObject("fTOF"));
    adEffTOF[iC] = static_cast<TF1*>(gradTOF->GetListOfFunctions()->FindObject("f"));
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
      delete hProjP;
      delete hProjS;
    }
  }
  //
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
  Color_t colors[4] = {kRed,kOrange,kYellow,kGreen};
  
  TH1F* centralityClass = (TH1F*)fDATA->Get("fCentralityClass");
  kNEvents[0] = centralityClass->GetBinContent(1);
  kNEvents[1] = centralityClass->GetBinContent(2);
  kNEvents[2] = centralityClass->GetBinContent(3);
  kNEvents[3] = centralityClass->GetBinContent(4);
  
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
//      sbCv[i] = new TCanvas(Form("keynotedSB%i_%i",iC,i));
//      sbCv[i]->Divide(3,3);
    }
    
    for (int iB = 0; iB < kNBins; ++iB) {
      if (iB < kTPCsize) {
        hRawCounts[iC].SetBinContent(iB + 1, hTPCcounts->GetBinContent(iB + 1));
      } else {
        TH1F *dat = (TH1F*)fDATA->Get(Form("Signal%i/fSignalD%i_%i",iC,iC,iB - kTPCsize));
        dat->Rebin(2);
        if (!dat) {
          cout << "Missing " << Form("Signal%i/fSignalD%i_%i",iC,iC,iB - kTPCsize);
          cout << " from DATA file\n";
          return;
        }
        const int n = (iB - kTPCsize) / 9;
        int shift = -1;
        if (n == 1) shift = 9 - 1;
        else if (n == 2) shift = 18 -1;
        summaryCv[n]->cd(iB - kTPCsize - shift);
        Float_t prim = 1.f;
       
        RooAddPdf &model = (kBins[iB + 1] > 1.4) ? model2 : model1;
        RooAbsPdf &bkg = (kBins[iB + 1] > 1.4) ? static_cast<RooAbsPdf&>(bkg0) : static_cast<RooAbsPdf&>(bkg2);
        RooAbsPdf &sig = (kBins[iB + 1] > 1.4) ? static_cast<RooAbsPdf&>(sig2) : static_cast<RooAbsPdf&>(sig0);
        RooDataHist data("data","data",RooArgList(m),Import(*dat));
        RooPlot *plot = m.frame();
        plot->SetTitle(Form("%1.1f #leq p_{T} < %1.1f",kBins[iB],kBins[iB + 1]));
        plot->SetName(Form("d%i_%i",iC,iB - kTPCsize));
        RooFitResult *res = model.fitTo(data,Extended());
        data.plotOn(plot);
        model.plotOn(plot);
        model.plotOn(plot,Components(bkg),LineStyle(kDashed));
        model.plotOn(plot,Components(sig),LineStyle(kDotted));
        model.paramOn(plot,Label(Form("#chi^{2}/NDF = %2.4f",plot->chiSquare())));
        plot->Draw();
        
        hRawCounts[iC].SetBinContent(iB + 1, fsig.getVal());
        hRawCounts[iC].SetBinError(iB + 1, fsig.getError());
//        sbCv[n]->cd(iB - kTPCsize - shift);

      }
    }
    //
    for (int i = 0; i < 3; ++i) {
      fOUTPUT->cd("Deuteron/Summary/Fit");
      summaryCv[i]->Write();
      delete summaryCv[i];
//      fOUTPUT->cd("Deuteron/Summary/SB");
//      sbCv[i]->Write();
//      delete sbCv[i];
      fOUTPUT->cd();
    }
    
    fOUTPUT->cd("Deuteron/Counts");
    hRawCounts[iC].Write();
    
    TObjArray obj(2);
    obj.SetOwner(kFALSE);
    for (int iB = 0; iB < kNBinsSecondaries; ++iB) {
      TH1F *hd;
      if (iB < kNTPCSecondaries) {
        hd = (TH1F*)fDATA->Get(Form("DCA%i/dcaxy_%i",iC,iB));
        if (!hd) {
          cout << "Missing " << Form("DCA%i/dcaxy_%i",iC,iB) << endl;
          return;
        }
      } else {
        TH2F *hDCASig = (TH2F*)fDATA->Get(Form("DCASignal/fDCASignal%i_%i",iC,iB - kTPCsize));
        if (!hDCASig) {
          cout << "Missing " << Form("DCASignal/fDCASignal%i_%i",iC,iB) << endl;
          return;
        }
        TH1D *hP = hDCASig->ProjectionY(Form("hDCATOFp%i_%i",iC,iB));
        hd = new TH1F(Form("hDCATOF%i_%i",iC,iB),";DCA_{xy} (cm);Entries",hP->GetNbinsX(),
                      hP->GetXaxis()->GetXmin(),hP->GetXaxis()->GetXmax());
        hd->Add(hP);
        delete hP;
      }
      
      Float_t prim = 1.f;
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
      hFractions[iC].SetBinContent(iB + 1,prim);
      hRawCounts[iC].SetBinContent(iB + 1, prim * hRawCounts[iC].GetBinContent(iB + 1));
      if (iB >= kNTPCSecondaries) delete hd;
    }
    
    fOUTPUT->cd("Deuteron/Fractions");
    hFractions[iC].Write();
    
    fOUTPUT->cd("Deuteron/Spectra");
    hSpectra[iC].SetMarkerStyle(26);
    hSpectra[iC].SetMarkerColor(colors[iC]);
    hSpectra[iC].SetLineColor(colors[iC]);
    for (int i = 1; i <= hSpectra[iC].GetNbinsX(); ++i) {
      float x = hSpectra[iC].GetBinCenter(i);
      float dx = hSpectra[iC].GetBinWidth(i);
      float eff = (i <= kTPCsize) ? dEff[iC]->Eval(x) : dEffTOF[iC]->Eval(x);
      if (eff > 0.05) {
        hSpectra[iC].SetBinContent(i, kNorm[iC] * hRawCounts[iC].GetBinContent(i) / (x * dx * eff));
        if (hRawCounts[iC].GetBinContent(i) > 0.f) {

          hSpectra[iC].SetBinError(i,
                                   hSpectra[iC].GetBinContent(i) *
                                   sqrt(sq(hRawCounts[iC].GetBinError(i) /
                                           hRawCounts[iC].GetBinContent(i)) + 1.f / kNEvents[iC]
                                        ));
        }
      }
    }
    hSpectra[iC].Write();

  }

  
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
          TH1F *dat = (TH1F*)fDATA->Get(Form("Signal%i/fSignalAD%i_%i",iC,iC,iB - kTPCsize));
          if (!dat) {
            cout << "Missing " << Form("Signal%i/fSignalAD%i_%i",iC,iC,iB - kTPCsize);
            cout << " from DATA file\n";
            return;
          }
          const int n = (iB - kTPCsize) / 9;
          int shift = -1;
          if (n == 1) shift = 9 - 1;
          else if (n == 2) shift = 18 -1;
          summaryCv[n]->cd(iB - kTPCsize - shift);
          
          RooAddPdf &model = (kBins[iB + 1] > 1.4) ? model2 : model1;
          RooAbsPdf &bkg = (kBins[iB + 1] > 1.4) ? static_cast<RooAbsPdf&>(bkg0) : static_cast<RooAbsPdf&>(bkg2);
          RooAbsPdf &sig = (kBins[iB + 1] > 1.4) ? static_cast<RooAbsPdf&>(sig2) : static_cast<RooAbsPdf&>(sig0);
          RooDataHist data("data","data",RooArgList(m),Import(*dat));
          RooPlot *plot = m.frame();
          plot->SetTitle(Form("%1.1f #leq p_{T} < %1.1f",kBins[iB],kBins[iB + 1]));
          plot->SetName(Form("d%i_%i",iC,iB - kTPCsize));
          RooFitResult *res = model.fitTo(data,Extended());
          data.plotOn(plot);
          model.plotOn(plot);
          model.plotOn(plot,Components(bkg),LineStyle(kDashed));
          model.plotOn(plot,Components(sig),LineStyle(kDotted));
          model.paramOn(plot,Label(Form("#chi^{2}/NDF = %2.4f",plot->chiSquare())));
          plot->Draw();
          
          hRawCounts[iC].SetBinContent(iB + 1, fsig.getVal());
          hRawCounts[iC].SetBinError(iB + 1, fsig.getError());
          //        sbCv[n]->cd(iB - kTPCsize - shift);
        }
      }
      for (int i = 0; i < 3; ++i) {
        fOUTPUT->cd("AntiDeuteron/Summary/Fit");
        summaryCv[i]->Write();
        delete summaryCv[i];
//        fOUTPUT->cd("AntiDeuteron/Summary/SB");
//        sbCv[i]->Write();
//        delete sbCv[i];
//        fOUTPUT->cd();
      }
      fOUTPUT->cd("AntiDeuteron/Counts");
      hSpectra[iC].Add(&hRawCounts[iC]);
      hRawCounts[iC].Write();
  
      fOUTPUT->cd("AntiDeuteron/Spectra");
      hSpectra[iC].SetMarkerStyle(26);
      hSpectra[iC].SetMarkerColor(colors[iC]);
      hSpectra[iC].SetLineColor(colors[iC]);
      for (int i = 1; i <= hSpectra[iC].GetNbinsX(); ++i) {
        float x = hSpectra[iC].GetBinCenter(i);
        float dx = hSpectra[iC].GetBinWidth(i);
        float eff = (i <= kTPCsize) ? adEff[iC]->Eval(x) : adEffTOF[iC]->Eval(x);
        if (eff > 0.05) {
          hSpectra[iC].SetBinContent(i, kNorm[iC] * hRawCounts[iC].GetBinContent(i) / (x * dx * eff));
          if (hRawCounts[iC].GetBinContent(i) > 0.f) {
            
            hSpectra[iC].SetBinError(i,
                                     hSpectra[iC].GetBinContent(i) *
                                     sqrt(sq(hRawCounts[iC].GetBinError(i) /
                                             hRawCounts[iC].GetBinContent(i)) + 1.f / kNEvents[iC]
                                          ));
          }
        }
      }
      hSpectra[iC].Write();
    }
  
  #pragma mark Closing files
    fDATA->Close();
    fMC->Close();
    fOUTPUT->Close();
    delete fOUTPUT;
}

#pragma mark Helper function definition
//__________________________________________________________________________________________________
Float_t CorrectForMaterial(TH1F* hd, TObjArray &obj, TFile* output) {
  //  TVirtualFitter::SetMaxIterations(10000000);
  hd->Rebin(kRebin);
  ((TH1F*)obj[0])->Rebin(kRebin);
  ((TH1F*)obj[1])->Rebin(kRebin);
  Int_t binLow = ((TH1F*)obj[0])->FindBin(-0.5);
  Int_t binUp = ((TH1F*)obj[0])->FindBin(0.5);
  TFractionFitter fitter(hd,&obj);
  fitter.SetRangeX(binLow,binUp);
  fitter.Constrain(0,0.,1.);
  fitter.Constrain(1,0.,0.8);
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
  return yieldPri;
}

//__________________________________________________________________________________________________
void CorrectForEfficiency(TH1F* rawcounts, TF1* eff, TF1* effTOF) {
  for (int iB = 1; iB < kNBins; ++iB) {
    float e = 0.f;
    if (iB <= kTPCsize) e = eff->Eval(rawcounts->GetBinCenter(iB));
    else e = effTOF->Eval(rawcounts->GetBinCenter(iB));
    if (e < kEfficiencyTolerance) rawcounts->SetBinContent(iB, 0.f);
    else rawcounts->SetBinContent(iB, rawcounts->GetBinContent(iB) / e);
  }
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
Double_t FitTemplateLowPt(Double_t *x_, Double_t *p) { // 5 parameters
  return p[0] * TailedGausPdf(x_, &p[1]) + p[3] * ExpoPdf(x_, &p[4]);
}

//__________________________________________________________________________________________________
Double_t FitTemplateHighPt(Double_t *x_, Double_t *p) { // 8 parameters
  return p[0] * TailedGausPdf(x_, &p[1]) + p[3] * DoubleExpoPdf(x_, &p[4]);
}

//__________________________________________________________________________________________________
Double_t ExpoPdf(Double_t *x_, Double_t *p) {
  Double_t x = x_[0] + kExpoPdfShift;
  return p[0] * TMath::Exp(-p[0] * x);
}

//__________________________________________________________________________________________________
Double_t DoubleExpoPdf(Double_t *x_, Double_t *p) {
  return (p[2] * ExpoPdf(x_, p) + p[3] * ExpoPdf(x_, &p[1])) / (p[2] + p[3]);
}

//__________________________________________________________________________________________________
Double_t TailedGausPdf(Double_t *x_, Double_t *p) {
  Double_t &x = x_[0];
  Double_t norm = 1.f / (Sqrt(TMath::PiOver2()) * p[0] * (1. + TMath::Erf(p[1] / Sqrt2())) +
                         p[0] * exp(-0.5 * p[1] * p[1]) / p[1]);
  if (x >= p[0] * p[1])
    return norm * exp(-p[1] * (x - 0.5 * p[0] * p[1])/p[0]);
  else
    return norm * exp(- x * x / (2 * p[0] * p[0]));
}
