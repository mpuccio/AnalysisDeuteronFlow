#include <TMath.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TFolder.h>
#include <TList.h>
#include <TH2F.h>
#include <TObjArray.h>
#include <TFractionFitter.h>

TH1F * HistoFromFunction(char *title,Double_t (*func)(Double_t*,Double_t*),Double_t *par) {
	TH1F * a =new TH1F(title,title,100,-2.4,2.4);
	Double_t x[1];
	for (int i = 1; i <= 100; ++i)
	{
		x[0] = a->GetBinCenter(i);
		a->SetBinContent(i,func(x,par));
	}
	return a;
}

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

Double_t ExpCB(Double_t *x_, Double_t *p) {
	Double_t &x = x_[0];
	return CrystalBall(x_,p) + p[5] * TMath::Exp(p[6] * x); //+ p[7];// + TMath::Gaus(x,p[3],p[4]) * p[5];	
}

Double_t GausSgl(Double_t *x, Double_t *p) {
	return p[0] * TMath::Gaus(x[0],p[1],p[2]);
}

Double_t ExpBkg(Double_t *x, Double_t *p) {
	return p[0] * TMath::Exp(p[1] * x[0]);
}

Double_t ExpExpBkg(Double_t *x, Double_t *p) {
	return p[0] * TMath::Exp(p[1] * x[0]) + p[2] * TMath::Exp(p[3] * x[0]);
}

Double_t ExpGaus(Double_t *x, Double_t *p) {
	return GausSgl(x,p) + ExpBkg(x,&p[3]);
}

Double_t ExpExpGaus(Double_t *x, Double_t *p) {
	return GausSgl(x,p) + ExpBkg(x,&p[3]) + ExpBkg(x,&p[5]);
}

void DeuteronSpectra() {
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	Double_t bins[] = {0.35,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.2,3.4,3.6,4.,4.5,5.,6.};
	Float_t binsPtTOF[19] = {0.8f,1.0f,1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,3.6f,4.0f,4.5f,5.0f,6.0f};
	TString titles[5] = {"Centrality 0-10%","Centrality 10-20%","Centrality 20-40%","Centrality 40-60%","Centrality 60-80%"};
	Color_t colors[5] = {kBlue-4,kGreen-3,kOrange,kOrange+8,kRed};
  
	TFile ResultFile("Results.root","recreate");
	ResultFile.cd();
	ResultFile.mkdir("Deuteron");
	ResultFile.mkdir("Deuteron/Spectra");
	ResultFile.mkdir("Deuteron/Efficiency");
	ResultFile.mkdir("Deuteron/Summary");
	ResultFile.mkdir("Deuteron/SB");
	ResultFile.mkdir("Deuteron/Fit");
  ResultFile.mkdir("Deuteron/Fractions");
  ResultFile.mkdir("Deuteron/Fractions/Fit");
  ResultFile.mkdir("AntiDeuteron");
	ResultFile.mkdir("AntiDeuteron/Spectra");
	ResultFile.mkdir("AntiDeuteron/Efficiency");
	ResultFile.mkdir("AntiDeuteron/Summary");
	ResultFile.mkdir("AntiDeuteron/SB");
	ResultFile.mkdir("AntiDeuteron/Fit");

  // Files necessary to finalize the analysis
  TFile *f = TFile::Open("AODSelector.root");
  TFile *fMC = TFile::Open("EfficiencyOutput.root");
  if (!f) {
    cout << "Missing data file" << endl;
    return;
  }
  if (!fMC) {
    cout << "Missing MC file" << endl;
    return;
  }
  
  // Efficiencies
  TEfficiency *dEff = (TEfficiency*)fMC->Get("Efficiencies/d");
  TEfficiency *dEffTOF = (TEfficiency*)fMC->Get("Efficiencies/dTOF");
  TEfficiency *adEff = (TEfficiency*)fMC->Get("Efficiencies/ad");
  TEfficiency *adEffTOF = (TEfficiency*)fMC->Get("Efficiencies/adTOF");
  
	// DEUTERON
  cout << "DEUTERON IN THE BOTTLE" << endl;
  cout << "===================================================" << endl;
  Double_t Nev = ((TH1F*)f->Get("fCentrality"))->GetEntries();
  
  TH1F *countsD[5];
	TH1F *rawSpectraD[5];
	TCanvas *spectraCanvas = new TCanvas("FullSpectra");
	spectraCanvas->SetLogy();
	for (int i = 0; i < 5; ++i)	{
		countsD[i] = new TH1F(Form("counts%i",i),";p_{T} (GeV/c);Raw counts",22,bins);
		rawSpectraD[i] = new TH1F(Form("spectra%i",i),";p_{T} (GeV/c);1/N_{ev}1/(2#pip_{T})d^{2}N/(dp_{T}dy) (GeV/c)^{-2}",22,bins);
	}

  cout << "PERFORMING THE QUANTUM LOOP" << endl;
	for (int k = 0; k < 5; ++k)
	{	
		TH1F* tpcY = (TH1F*)f->Get(Form("fdEdxTPCSignalCounts%i",k));
		for (int i = 1; i <= 4; ++i)
		{
      countsD[k]->SetBinContent(i,tpcY->GetBinContent(i));
			rawSpectraD[k]->SetBinContent(i,tpcY->GetBinContent(i)/(TMath::TwoPi()*rawSpectraD[k]->GetBinCenter(i)*Nev*rawSpectraD[k]->GetBinWidth(i)*dEff->GetEfficiency(i)));
		}
	}
  
  cout << "OPTIMIZING FOR INEFFICIENT CONVERSION" << endl;
  TH1D *hDCATOF[10];
	for (int iCent = 0; iCent < 5; ++iCent)
	{
		TCanvas *cKeyNote = new TCanvas(Form("keynoted%i",iCent));
		cKeyNote->Divide(4,3);
		TCanvas *cKeyNoteSB = new TCanvas(Form("keynoteSBd%i",iCent));
		cKeyNoteSB->Divide(4,3);
		for (int iPt = 0; iPt < 17; ++iPt)
		{
			TCanvas *cv = new TCanvas(Form("Cent%i_Pt_%i",iCent,iPt),Form("Cent%i_Pt_%i",iCent,iPt));
			TH1F *data = (TH1F*)f->Get(Form("Signal%i/fSignalD%i_%i",iCent,iCent,iPt));
      if (!data) {
        cout << "Missing " << Form("Signal%i/fSignalD%i_%i",iCent,iCent,iPt) << " from DATA file\n";
        return;
      }
			data->SetMarkerStyle(20);
			data->SetMarkerSize(0.5);
			data->SetMarkerColor(kBlue);
			TF1 *tpl;
			if (iPt < 4) {
				tpl = new TF1(Form("tpl%i_%i",iCent,iPt),ExpCB,-2.4,2.4,7);
				tpl->SetParameters((iCent == 0) ? 20000 : 4000, 0.f, 2.86394e-01, 1, 2, 42,-0.6);
				tpl->SetParNames("N","#mu","#sigma","#alpha","n","A","#tau");
			} else {
				tpl = new TF1(Form("tpl%i_%i",iCent,iPt),ExpExpGaus,-2.4,2.4,7);
				tpl->SetParameters((iPt > 10) ? 500 : 5000, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);
        if (iCent >= 1) {
          tpl->SetParameters((iPt > 14) ? 40 : 500, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);
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
			TH1F* background;
			TH1F* signl;
			TH1F* templ;
			if (iPt < 4) {
				background = HistoFromFunction(Form("background%i_%i",iCent,iPt),ExpBkg,&parameters[5]);
				signl = HistoFromFunction(Form("signl%i_%i",iCent,iPt),CrystalBall,parameters);
				templ = HistoFromFunction(Form("templ%i_%i",iCent,iPt),ExpCB,parameters);
			} else {
				background = HistoFromFunction(Form("background%i_%i",iCent,iPt),ExpExpBkg,&parameters[3]);
				signl = HistoFromFunction(Form("signl%i_%i",iCent,iPt),GausSgl,parameters);
				templ = HistoFromFunction(Form("templ%i_%i",iCent,iPt),ExpExpGaus,parameters);
			}
			TCanvas *sb = new TCanvas(Form("SB%i_%i",iCent,iPt));
			sb->cd();
			TH1F backFraction(*background);
			TH1F siglFraction(*signl);
			backFraction.SetLineColor(kBlack);
			backFraction.SetNameTitle(Form("backFraction%i_%i",iCent,iPt),Form("%4.1f #leq p_{T} < %4.1f ;m_{TOF}^{2}-m_{PDG}^{2};S/B",binsPtTOF[iPt],binsPtTOF[iPt+1]));
			backFraction.Divide(templ);
			backFraction.GetYaxis()->SetRangeUser(0,1.05);
			siglFraction.SetNameTitle(Form("siglFraction%i_%i",iCent,iPt),Form("%4.1f #leq p_{T} < %4.1f ;m_{TOF}^{2}-m_{PDG}^{2};S/B",binsPtTOF[iPt],binsPtTOF[iPt+1])); 
			siglFraction.SetLineColor(kRed);
			siglFraction.Divide(templ);
			backFraction.Draw("l");
			siglFraction.Draw("lsame");
			TLegend ll1(0.72,0.425,0.97,0.61);
			ll1.SetFillColor(kWhite);
			ll1.SetLineColor(kBlack);
			ll1.AddEntry(&backFraction,"Background fraction","l");
			ll1.AddEntry(&siglFraction,"Signal fraction","l");
			ll1.Draw();
			cKeyNoteSB->cd(iPt+1);
			backFraction.Draw("l");
			siglFraction.Draw("lsame");
			ll1.Draw();
			background->SetLineColor(kRed);
			background->SetLineWidth(1);
			background->SetLineStyle(2);
			data->GetXaxis()->SetRangeUser(-2.,2.);
			cv->cd();
			data->Draw("e");
			background->Draw("lsame");
			if (iPt < 12) {
				cKeyNote->cd(iPt+1);
				data->Draw("e");
				background->Draw("lsame");
			}
			templ->Add(background,-1.f);
      if (iCent < 4 && iPt < 16) {
        countsD[iCent]->SetBinContent(iPt+5,(Int_t)templ->Integral(
                                                                   templ->FindBin(parameters[1] - 3.f * parameters[2]),
                                                                   templ->FindBin(parameters[1] + 3.f * parameters[2])
                                                                   ));
        if (iPt < 2) {
          TH2F *hhh = (TH2F*)f->Get(Form("DCASignal/fDCASignal%i_%i",iCent,iPt));
          hDCATOF[iCent * 2 + iPt] = hhh->ProjectionY(Form("hDCATOF%i_%i",iCent,iPt),templ->FindBin(parameters[1] - 3.f * parameters[2]),
                                                      templ->FindBin(parameters[1] + 3.f * parameters[2]));
        }
//        rawSpectraD[iCent]->SetBinContent(iPt+5,countsD[iCent]->GetBinContent(iPt+5) /
//                                          (TMath::TwoPi()*rawSpectraD[iCent]->GetBinCenter(iPt+5)*
//                                           Nev*rawSpectraD[iCent]->GetBinWidth(iPt+5)*
//                                           dEffTOF->GetEfficiency(iPt+5)));
      }
			ResultFile.cd("Deuteron/Fit");
			cv->Write();
			delete cv;

			ResultFile.cd("Deuteron/SB");
			sb->Write();
			delete sb;
		}
		ResultFile.cd("Deuteron/Summary");
		cKeyNote->Write();
		cKeyNoteSB->Write();
		delete cKeyNote;
		delete cKeyNoteSB;
	}
  
  // Primary / secondary fractions
  TH1F *hPrimaryFraction[5];
  TH1F *hSecondaryFraction[5];
  for (int i = 0; i < 4; ++i) {
    hPrimaryFraction[i] = new TH1F(Form("hPrimaryFraction%i",i),";p_{T} (GeV/c);Primary fraction",22,bins);
    hSecondaryFraction[i] = new TH1F(Form("hPSecondaryFraction%i",i),";p_{T} (GeV/c);Secondary fraction",22,bins);
    TH2F *hDCAp = (TH2F*)fMC->Get(Form("Cent%i/fDdcaXYprimaries_%i",i,i));
    TH2F *hDCAs = (TH2F*)fMC->Get(Form("Cent%i/fDdcaXYsecondaries_%i",i,i));
    if (!hDCAp || !hDCAs) {
      cout << "Missing MC DCAs" << endl;
      return;
    }
    TH1F* hp = new TH1F("hp",";DCA_{xy} (cm); Entries",hDCAp->GetNbinsY(),
                        hDCAp->GetYaxis()->GetXmin(), hDCAp->GetYaxis()->GetXmax());
    TH1F* hs = new TH1F("hs",";DCA_{xy} (cm); Entries",hDCAs->GetNbinsY(),
                        hDCAs->GetYaxis()->GetXmin(), hDCAs->GetYaxis()->GetXmax());
    TObjArray obj(2);
    obj.SetOwner(kTRUE);
    obj.Add(hp);
    obj.Add(hs);
    for (int k = 0; k < 6; ++k) {
      TH1F* hd;
      if (k < 4) hd = (TH1F*)f->Get(Form("DCA%i/dcaxy_%i",i,k));
      else hd = (TH1F*)hDCATOF[i * 2 + k - 4];
      if (!hd) {
        cout << "Missing data " << Form("DCA%i/dcaxy_%i",i,k) << endl;
        return;
      }
      hd->Rebin(4);
      hp->Add(hDCAp->ProjectionY("py",k + 1));
      hs->Add(hDCAs->ProjectionY("sy",k + 1));
      
      Int_t binLow = hs->FindBin(-0.5);
      Int_t binUp = hs->FindBin(0.5);
      TFractionFitter fitter(hd,&obj);
      fitter.SetRangeX(binLow,binUp);
      fitter.Constrain(0,0.,1.);
      fitter.Constrain(1,0.,1.);
      Int_t result = fitter.Fit();
      if (result == 0) {
        Double_t yieldSec,yieldPri,error;
        fitter.GetResult(0, yieldPri, error);
        fitter.GetResult(1,yieldSec,error);
        Float_t dataIntegral = hd->Integral();
        cout << "YIELD: " << yieldPri << " " << yieldSec << " " << dataIntegral <<  endl;
        TH1F* hfit = (TH1F*)fitter.GetPlot();
        hfit->SetLineColor(kGreen + 1);
        hfit->SetLineWidth(3);
        hs->Scale(1. / hs->Integral());
        hp->Scale(1. / hp->Integral());
        hs->Scale(yieldSec * dataIntegral);
        hp->Scale(yieldPri * dataIntegral);
        hPrimaryFraction[i]->SetBinContent(k + 1, yieldPri * dataIntegral);
        hSecondaryFraction[i]->SetBinContent(k + 1, yieldSec * dataIntegral);
        TCanvas cv(Form("DCA%i_%i",i,k));
        cv.cd();
        hd->SetMarkerStyle(20);
        hd->SetMarkerSize(0.5);
        hd->SetMarkerColor(kBlack);
        hd->Draw("e");
        
        hfit->Draw("same");
        hs->SetLineColor(kRed);
        hp->SetLineColor(kBlue);
        hs->Draw("same");
        hp->Draw("same");
        ResultFile.cd("Deuteron/Fractions/Fit");
        cv.Write();
      }
      hp->Reset();
      hs->Reset();
    }
    ResultFile.cd("Deuteron/Fractions/");
    hSecondaryFraction[i]->Write();
  }
  
  // Finalising spectra
  for (int iCent = 0; iCent < 5; ++iCent) {
    countsD[iCent]->Add(hSecondaryFraction[iCent],-1.f);
    for (int iBin = 1; iBin <= countsD[iCent]->GetNbinsX(); ++iBin) {
      Float_t efficiency = (iBin < 5) ? dEff->GetEfficiency(iBin) : dEffTOF->GetEfficiency(iBin);
      rawSpectraD[iCent]->SetBinContent(iBin,countsD[iCent]->GetBinContent(iBin) /
                                        (TMath::TwoPi() * rawSpectraD[iCent]->GetBinCenter(iBin) *
                                         Nev * rawSpectraD[iCent]->GetBinWidth(iBin) * efficiency));
    }
  }
  
  ResultFile.cd();
  
  cout << "WRITING DOWN OBSCURE FORMULAS" << endl;
	TLegend ll2(0.62,0.625,0.87,0.81);
	ll2.SetFillColor(kWhite);
	ll2.SetLineColor(kWhite);
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 1; j <= rawSpectraD[i]->GetEntries(); ++j)
		{
			rawSpectraD[i]->SetBinError(j,1e-12);
		}
		rawSpectraD[i]->SetMarkerStyle(8);
		rawSpectraD[i]->SetMarkerColor(colors[i]);
		rawSpectraD[i]->GetYaxis()->SetRangeUser(1e-7,20);
		rawSpectraD[i]->GetXaxis()->SetRangeUser(1e-2,10);
    if(i!=4) ll2.AddEntry(rawSpectraD[i],titles[i].Data(),"lp");
    spectraCanvas->cd();
    switch (i) {
      case 0:
        rawSpectraD[i]->Draw("e");
        break;
      case 4:
        break;
      default:
        rawSpectraD[i]->Draw("esame");
        break;
    }
		ResultFile.cd("Deuteron/Spectra");
		rawSpectraD[i]->Write();
		countsD[i]->Write();
  }
  spectraCanvas->cd();
  ll2.Draw();
  ResultFile.cd("Deuteron/Spectra");
  spectraCanvas->Write();
  delete spectraCanvas;
	for (int i = 0; i < 5; ++i)
	{
		delete rawSpectraD[i];
		delete countsD[i];
	}

	ResultFile.cd("Deuteron/Efficiency");
	TGraphAsymmErrors* grdeff = dEff->CreateGraph();
	grdeff->SetName("dEffTPC");
	grdeff->Write();
	TGraphAsymmErrors* grdeffTOF = dEffTOF->CreateGraph();
	grdeffTOF->SetName("dEffTPCTOF");
	grdeffTOF->Write();
	cout << "===================================================" << endl;
  
	// ANTIDEUTERON
  cout << "GOING TO THE OTHER SIDE OF C" << endl;
  cout << "===================================================" << endl;
	TH1F *countsAD[5];
	TH1F *rawSpectraAD[5];
	TCanvas *spectraCanvasAD = new TCanvas("FullSpectra");
	spectraCanvasAD->SetLogy();
	for (int i = 0; i < 5; ++i)	{
		countsAD[i] = new TH1F(Form("countsAD%i",i),";p_{T} (GeV/c);Raw counts",20,bins);
		rawSpectraAD[i] = new TH1F(Form("spectraAD%i",i),";p_{T} (GeV/c);1/N_{ev}1/(2#pip_{T})d^{2}N/(dp_{T}dy) (GeV/c)^{-2}",20,bins);
	}

  cout << "PROJECTION IN A DARK CHAMBER" << endl;
	for (int k = 0; k < 5; ++k)
	{	
		TH1F* tpcY = (TH1F*)f->Get(Form("fdEdxTPCSignalCountsAD%i",k));
    if (!tpcY) {
      cout << "Missing " << Form("fdEdxTPCSignalCountsAD%i",k) << " from DATA file\n";
      return;
    }
		for (int i = 1; i <= 4; ++i)
		{
			countsAD[k]->SetBinContent(i,tpcY->GetBinContent(i));
			if(adEff->GetEfficiency(i) == 0) rawSpectraAD[k]->SetBinContent(i,0);
			else rawSpectraAD[k]->SetBinContent(i,tpcY->GetBinContent(i)/(TMath::TwoPi()*rawSpectraAD[k]->GetBinCenter(i)*Nev*rawSpectraAD[k]->GetBinWidth(i)*adEff->GetEfficiency(i)));
		}
	}
	
  cout << "WRITING THE HISTORY WITH THE FUTURE SIMPLE" << endl;
	for (int iCent = 0; iCent < 5; ++iCent)
	{
		TCanvas *cKeyNote = new TCanvas(Form("keynoted%i",iCent));
		cKeyNote->Divide(4,3);
		TCanvas *cKeyNoteSB = new TCanvas(Form("keynoteSBd%i",iCent));
		cKeyNoteSB->Divide(4,3);
		for (int iPt = 0; iPt < 17; ++iPt)
		{
			TCanvas *cv = new TCanvas(Form("Cent%i_Pt_%i",iCent,iPt),Form("Cent%i_Pt_%i",iCent,iPt));
			TH1F *data = (TH1F*)f->Get(Form("Signal%i/fSignalAD%i_%i",iCent,iCent,iPt));
      if (!data) {
        cout << "Missing " << Form("Signal%i/fSignalAD%i_%i",iCent,iCent,iPt) << " from DATA file\n";
        return;
      }
      data->SetMarkerStyle(20);
			data->SetMarkerSize(0.5);
			data->SetMarkerColor(kBlue);
			TF1 *tpl;
			if (iPt < 4) {
				tpl = new TF1(Form("tpl%i_%i",iCent,iPt),ExpCB,-2.4,2.4,7);
				tpl->SetParameters(700, 0.f, 2.86394e-01, 0.4, 3, (iCent==0&&iPt==3) ? 1000 : 40,-0.6);
				tpl->SetParLimits(4,1,20);
				tpl->SetParNames("N","#mu","#sigma","#alpha","n","A","#tau");
			} else {
				tpl = new TF1(Form("tpl%i_%i",iCent,iPt),ExpExpGaus,-2.4,2.4,7);
				Float_t high = 400.f, low = 100.f;
        tpl->SetParameters((iPt > 9) ? low : high, 5.66003e-02, 2.86394e-01, 1.28095e+00,-2.04126e+00, 3.51728e+02,-2.00000e-01);
				if (iCent >= 2)
				{
					high = 400.f;
					low = 50.f;
          tpl->SetParameters((iPt > 9) ? low : high, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-1.00000e-01);
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
			TH1F* background;
			TH1F* signl;
			TH1F* templ;
			if (iPt < 4) {
				background = HistoFromFunction(Form("background%i_%i",iCent,iPt),ExpBkg,&parameters[5]);
				signl = HistoFromFunction(Form("signl%i_%i",iCent,iPt),CrystalBall,parameters);
				templ = HistoFromFunction(Form("templ%i_%i",iCent,iPt),ExpCB,parameters);
			} else {
				background = HistoFromFunction(Form("background%i_%i",iCent,iPt),ExpExpBkg,&parameters[3]);
				signl = HistoFromFunction(Form("signl%i_%i",iCent,iPt),GausSgl,parameters);
				templ = HistoFromFunction(Form("templ%i_%i",iCent,iPt),ExpExpGaus,parameters);
			}
			TCanvas *sb = new TCanvas(Form("SB%i_%i",iCent,iPt));
			sb->cd();
			TH1F backFraction(*background);
			TH1F siglFraction(*signl);
			backFraction.SetLineColor(kBlack);
			backFraction.SetNameTitle(Form("backFraction%i_%i",iCent,iPt),
                                Form("%4.1f #leq p_{T} < %4.1f ;m_{TOF}^{2}-m_{PDG}^{2};S/B",
                                     binsPtTOF[iPt],binsPtTOF[iPt+1]));
			backFraction.Divide(templ);
			backFraction.GetYaxis()->SetRangeUser(0,1.05);
			siglFraction.SetNameTitle(Form("siglFraction%i_%i",iCent,iPt),
                                Form("%4.1f #leq p_{T} < %4.1f ;m_{TOF}^{2}-m_{PDG}^{2};S/B",
                                     binsPtTOF[iPt],binsPtTOF[iPt+1]));
			siglFraction.SetLineColor(kRed);
			siglFraction.Divide(templ);
			backFraction.Draw("l");
			siglFraction.Draw("lsame");
			TLegend ll1(0.72,0.425,0.97,0.61);
			ll1.SetFillColor(kWhite);
			ll1.SetLineColor(kBlack);
			ll1.AddEntry(&backFraction,"Background fraction","l");
			ll1.AddEntry(&siglFraction,"Signal fraction","l");
			ll1.Draw();
			cKeyNoteSB->cd(iPt+1);
			backFraction.Draw("l");
			siglFraction.Draw("lsame");
			ll1.Draw();
			background->SetLineColor(kRed);
			background->SetLineWidth(1);
			background->SetLineStyle(2);
			data->GetXaxis()->SetRangeUser(-2.,2.);
			cv->cd();
			data->Draw("e");
			background->Draw("lsame");
			if (iPt < 12) {
				cKeyNote->cd(iPt+1);
				data->Draw("e");
				background->Draw("lsame");
			}
			templ->Add(background,-1.f);
      if (iCent < 4 && iPt < 16) {
        countsAD[iCent]->SetBinContent(iPt+5,(Int_t)templ->Integral(
                                                                    templ->FindBin(parameters[1] - 3.f * parameters[2]),
                                                                    templ->FindBin(parameters[1] + 3.f * parameters[2])
                                                                    ));
        if (iPt < 2) {
          
        }
        rawSpectraAD[iCent]->SetBinContent(iPt+5,countsAD[iCent]->GetBinContent(iPt+5) /
                                           (TMath::TwoPi()*rawSpectraAD[iCent]->GetBinCenter(iPt+5)*
                                            Nev*rawSpectraAD[iCent]->GetBinWidth(iPt+5)*
                                            adEffTOF->GetEfficiency(iPt+5)));
      }
			ResultFile.cd("AntiDeuteron/Fit");
			cv->Write();
			delete cv;

			ResultFile.cd("AntiDeuteron/SB");
			sb->Write();
			delete sb;
		}
		ResultFile.cd("AntiDeuteron/Summary");
		cKeyNote->Write();
		cKeyNoteSB->Write();
	}

  cout << "FIXING THE UNKNOWN" << endl;
	TLegend llAD(0.62,0.625,0.87,0.81);
	llAD.SetFillColor(kWhite);
	llAD.SetLineColor(kWhite);
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 1; j <= rawSpectraAD[i]->GetEntries(); ++j)
		{
			rawSpectraAD[i]->SetBinError(j,1e-12);
		}
		rawSpectraAD[i]->SetMarkerStyle(8);
		rawSpectraAD[i]->SetMarkerColor(colors[i]);
		rawSpectraAD[i]->GetYaxis()->SetRangeUser(1e-7,20);
		rawSpectraAD[i]->GetXaxis()->SetRangeUser(1e-2,10);
    if(i!=4)llAD.AddEntry(rawSpectraAD[i],titles[i].Data(),"lp");
		spectraCanvasAD->cd();
    switch (i) {
      case 0:
        rawSpectraAD[i]->Draw("e");
        break;
      case 4:
        break;
      default:
        rawSpectraAD[i]->Draw("esame");
        break;
    }
		ResultFile.cd("AntiDeuteron/Spectra");
		rawSpectraAD[i]->Write();
		countsAD[i]->Write();
  }
  spectraCanvasAD->cd();
  llAD.Draw();

	ResultFile.cd("AntiDeuteron/Spectra");
	spectraCanvasAD->Write();
	delete spectraCanvasAD;
	for (int i = 0; i < 5; ++i)
	{
		delete countsAD[i];
		delete rawSpectraAD[i];
	}

	ResultFile.cd("AntiDeuteron/Efficiency");
	// adEff->Paint();
	TGraphAsymmErrors* gradeff = adEff->CreateGraph();
	gradeff->SetName("adEffTPC");
	gradeff->Write();
	// adEffTOF->Paint();
	TGraphAsymmErrors* gradeffTOF = adEffTOF->CreateGraph();
	gradeffTOF->SetName("adEffTPCTOF");
	gradeffTOF->Write();
	
	f->Close();
  fMC->Close();

	ResultFile.Close();
}