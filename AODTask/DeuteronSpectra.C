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
	Double_t bins[] = {0.35,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.};
	Float_t binsPtTOF[14] = {0.8f,1.0f,1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,3.0f,3.5f,4.0f,5.0f};
	TString titles[3] = {"Centrality 0-5%","Centrality 20-40%","Centrality 40-60%"};
	Color_t colors[] = {kRed,kOrange,kBlue};

	TFile fEff("finalEff.root");
	TEfficiency *dEff = (TEfficiency*)fEff.Get("db");
	TEfficiency *dEffTOF = (TEfficiency*)fEff.Get("dTOFb");
	TEfficiency *adEff = (TEfficiency*)fEff.Get("adb");
	TEfficiency *adEffTOF = (TEfficiency*)fEff.Get("adTOFb");
	TFile ResultFile("Results.root","recreate");
	ResultFile.cd();
	ResultFile.mkdir("Deuteron");
	ResultFile.mkdir("Deuteron/Spectra");
	ResultFile.mkdir("Deuteron/Efficiency");
	ResultFile.mkdir("Deuteron/Summary");
	ResultFile.mkdir("Deuteron/SB");
	ResultFile.mkdir("Deuteron/Fit");
	ResultFile.mkdir("AntiDeuteron");
	ResultFile.mkdir("AntiDeuteron/Spectra");
	ResultFile.mkdir("AntiDeuteron/Efficiency");
	ResultFile.mkdir("AntiDeuteron/Summary");
	ResultFile.mkdir("AntiDeuteron/SB");
	ResultFile.mkdir("AntiDeuteron/Fit");

	// DEUTERON
	TFile *f = TFile::Open("AODSelector.root");
	Double_t Nev = ((TH1F*)f->Get("fCentrality"))->GetEntries();

	TH1F *countsD[3];
	TH1F *rawSpectraD[3];
	TCanvas *spectraCanvas = new TCanvas("FullSpectra");
	spectraCanvas->SetLogy();
	for (int i = 0; i < 3; ++i)	{
		countsD[i] = new TH1F(Form("counts%i",i),";p_{T} (GeV/c);Raw counts",17,bins);
		rawSpectraD[i] = new TH1F(Form("spectra%i",i),";p_{T} (GeV/c);1/N_{ev}1/(2#pip_{T})d^{2}N/(dp_{T}dy) (GeV/c)^{-2}",17,bins);
	}

	for (int k = 0; k < 3; ++k)
	{	
		TH1F* tpcY = (TH1F*)f->Get(Form("fdEdxTPCSignalCounts%i",k));
		for (int i = 1; i <= 4; ++i)
		{
			countsD[k]->SetBinContent(i,tpcY->GetBinContent(i));
			rawSpectraD[k]->SetBinContent(i,tpcY->GetBinContent(i)/(TMath::TwoPi()*rawSpectraD[k]->GetBinCenter(i)*Nev*rawSpectraD[k]->GetBinWidth(i)*dEff->GetEfficiency(i)));
		}
	}
	
	for (int iCent = 0; iCent < 3; ++iCent)
	{
		TCanvas *cKeyNote = new TCanvas(Form("keynoted%i",iCent));
		cKeyNote->Divide(4,3);
		TCanvas *cKeyNoteSB = new TCanvas(Form("keynoteSBd%i",iCent));
		cKeyNoteSB->Divide(4,3);
		for (int iPt = 0; iPt < 13; ++iPt)
		{
			TCanvas *cv = new TCanvas(Form("Cent%i_Pt_%i",iCent,iPt),Form("Cent%i_Pt_%i",iCent,iPt));
			TH1F *data = (TH1F*)f->Get(Form("fSignal%i_%i",iCent,iPt));
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
			backFraction.DrawClone("l");
			siglFraction.DrawClone("lsame");
			TLegend ll1(0.72,0.425,0.97,0.61);
			ll1.SetFillColor(kWhite);
			ll1.SetLineColor(kBlack);
			ll1.AddEntry(&backFraction,"Background fraction","l");
			ll1.AddEntry(&siglFraction,"Signal fraction","l");
			ll1.DrawClone();
			cKeyNoteSB->cd(iPt+1);
			backFraction.DrawClone("l");
			siglFraction.DrawClone("lsame");
			ll1.DrawClone();
			background->SetLineColor(kRed);
			background->SetLineWidth(1);
			background->SetLineStyle(2);
			data->GetXaxis()->SetRangeUser(-2.,2.);
			cv->cd();
			data->DrawClone("e");
			background->DrawClone("lsame");
			if (iPt < 12) {
				cKeyNote->cd(iPt+1);
				data->DrawClone("e");
				background->DrawClone("lsame");
			}
			templ->Add(background,-1.f);
			countsD[iCent]->SetBinContent(iPt+5,(Int_t)templ->Integral(
				templ->FindBin(parameters[1] - 3.f * parameters[2]),
				templ->FindBin(parameters[1] + 3.f * parameters[2])
				));
			rawSpectraD[iCent]->SetBinContent(iPt+5,countsD[iCent]->GetBinContent(iPt+5) / (TMath::TwoPi()*rawSpectraD[iCent]->GetBinCenter(iPt+5)*Nev*rawSpectraD[iCent]->GetBinWidth(iPt+5)*dEffTOF->GetEfficiency(iPt+5)));

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

	TLegend ll2(0.62,0.625,0.87,0.81);
	ll2.SetFillColor(kWhite);
	ll2.SetLineColor(kWhite);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 1; j <= rawSpectraD[i]->GetEntries(); ++j)
		{
			rawSpectraD[i]->SetBinError(j,1e-12);
		}
		rawSpectraD[i]->SetMarkerStyle(8);
		rawSpectraD[i]->SetMarkerColor(colors[i]);
		rawSpectraD[i]->GetYaxis()->SetRangeUser(1e-6,20);
		rawSpectraD[i]->GetXaxis()->SetRangeUser(1e-2,6);		
    ll2.AddEntry(rawSpectraD[i],titles[i].Data(),"lp");
		spectraCanvas->cd();
		rawSpectraD[i]->DrawClone((i==0) ? "e" : "esame");
		ResultFile.cd("Deuteron/Spectra");
		rawSpectraD[i]->Write();
		countsD[i]->Write();
  }
  spectraCanvas->cd();
  ll2.Draw();
  ResultFile.cd("Deuteron/Spectra");
  spectraCanvas->Write();
  delete spectraCanvas;
	for (int i = 0; i < 3; ++i)
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
	
	// ANTIDEUTERON

	TH1F *countsAD[3];
	TH1F *rawSpectraAD[3];
	TCanvas *spectraCanvasAD = new TCanvas("FullSpectra");
	spectraCanvasAD->SetLogy();
	for (int i = 0; i < 3; ++i)	{
		countsAD[i] = new TH1F(Form("countsAD%i",i),";p_{T} (GeV/c);Raw counts",17,bins);
		rawSpectraAD[i] = new TH1F(Form("spectraAD%i",i),";p_{T} (GeV/c);1/N_{ev}1/(2#pip_{T})d^{2}N/(dp_{T}dy) (GeV/c)^{-2}",17,bins);
	}

	for (int k = 0; k < 3; ++k)
	{	
		TH1F* tpcY = (TH1F*)f->Get(Form("fdEdxTPCSignalCountsAD%i",k));
		for (int i = 1; i <= 4; ++i)
		{
			countsAD[k]->SetBinContent(i,tpcY->GetBinContent(i));
			if(adEff->GetEfficiency(i) == 0) rawSpectraAD[k]->SetBinContent(i,0);
			else rawSpectraAD[k]->SetBinContent(i,tpcY->GetBinContent(i)/(TMath::TwoPi()*rawSpectraAD[k]->GetBinCenter(i)*Nev*rawSpectraAD[k]->GetBinWidth(i)*adEff->GetEfficiency(i)));
		}
	}
	
	for (int iCent = 0; iCent < 3; ++iCent)
	{
		TCanvas *cKeyNote = new TCanvas(Form("keynoted%i",iCent));
		cKeyNote->Divide(4,3);
		TCanvas *cKeyNoteSB = new TCanvas(Form("keynoteSBd%i",iCent));
		cKeyNoteSB->Divide(4,3);
		for (int iPt = 0; iPt < 13; ++iPt)
		{
			TCanvas *cv = new TCanvas(Form("Cent%i_Pt_%i",iCent,iPt),Form("Cent%i_Pt_%i",iCent,iPt));
			TH1F *data = (TH1F*)f->Get(Form("fSignalAD%i_%i",iCent,iPt));
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
				Float_t high = 500.f, low = 100.f;
				if (iCent != 2)
				{
					high = 2000.f;
					low = 500.f;
				}
				tpl->SetParameters((iPt > 9) ? low : high, 5.66003e-02, 2.86394e-01, 4.28095e+00,-4.04126e+00, 6.51728e+02,-3.00000e-01);
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
			backFraction.DrawClone("l");
			siglFraction.DrawClone("lsame");
			TLegend ll1(0.72,0.425,0.97,0.61);
			ll1.SetFillColor(kWhite);
			ll1.SetLineColor(kBlack);
			ll1.AddEntry(&backFraction,"Background fraction","l");
			ll1.AddEntry(&siglFraction,"Signal fraction","l");
			ll1.DrawClone();
			cKeyNoteSB->cd(iPt+1);
			backFraction.DrawClone("l");
			siglFraction.DrawClone("lsame");
			ll1.DrawClone();
			background->SetLineColor(kRed);
			background->SetLineWidth(1);
			background->SetLineStyle(2);
			data->GetXaxis()->SetRangeUser(-2.,2.);
			cv->cd();
			data->DrawClone("e");
			background->DrawClone("lsame");
			if (iPt < 12) {
				cKeyNote->cd(iPt+1);
				data->DrawClone("e");
				background->DrawClone("lsame");
			}
			templ->Add(background,-1.f);
			countsAD[iCent]->SetBinContent(iPt+5,(Int_t)templ->Integral(
				templ->FindBin(parameters[1] - 3.f * parameters[2]),
				templ->FindBin(parameters[1] + 3.f * parameters[2])
				));
			rawSpectraAD[iCent]->SetBinContent(iPt+5,countsAD[iCent]->GetBinContent(iPt+5) / (TMath::TwoPi()*rawSpectraAD[iCent]->GetBinCenter(iPt+5)*Nev*rawSpectraAD[iCent]->GetBinWidth(iPt+5)*adEffTOF->GetEfficiency(iPt+5)));

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

	TLegend llAD(0.62,0.625,0.87,0.81);
	llAD.SetFillColor(kWhite);
	llAD.SetLineColor(kWhite);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 1; j <= rawSpectraAD[i]->GetEntries(); ++j)
		{
			rawSpectraAD[i]->SetBinError(j,1e-12);
		}
		rawSpectraAD[i]->SetMarkerStyle(8);
		rawSpectraAD[i]->SetMarkerColor(colors[i]);
		rawSpectraAD[i]->GetYaxis()->SetRangeUser(1e-6,20);
		rawSpectraAD[i]->GetXaxis()->SetRangeUser(1e-2,6);		
    llAD.AddEntry(rawSpectraAD[i],titles[i].Data(),"lp");
		spectraCanvasAD->cd();
		rawSpectraAD[i]->DrawClone((i==0) ? "e" : "esame");
		ResultFile.cd("AntiDeuteron/Spectra");
		rawSpectraAD[i]->Write();
		countsAD[i]->Write();
  }
  spectraCanvasAD->cd();
  llAD.DrawClone();

	ResultFile.cd("AntiDeuteron/Spectra");
	spectraCanvasAD->Write();
	delete spectraCanvasAD;
	for (int i = 0; i < 3; ++i)
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

	ResultFile.Close();
}