//
//  AliAnalysisTaskFlowd.cxx
//
//  Created by Maximiliano Puccio on 14/09/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#include "AliAnalysisTaskFlowd.h"
#include <Riostream.h>
#include "TChain.h"
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliTOFPIDResponse.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskFlowd)

//__________________________________________________________________________________________________
inline static Double_t DeuteronTPC(Double_t x) {
  // Deuteron expected signal in TPC
  return AliExternalTrackParam::BetheBlochAleph(x/1.875612,4.69637,7.51827,0.0183746,2.60,2.7);
}

//__________________________________________________________________________________________________
AliAnalysisTaskFlowd::AliAnalysisTaskFlowd(const char* name)
:AliAnalysisTaskSE(name)
,fEvent(0x0)
,fESDtrackCuts(0x0)
,fPid(0x0)
,fEventHandler(0x0)
,fFillTree(kFALSE)
,fHistCentralityClass10(0x0)
,fHistCentralityPercentile(0x0)
,fHistDeDx(0x0)
,fHistDeDxITS(0x0)
,fHistDeuteron(0x0)
,fHistTOF2D(0x0)
,fHistTOFnuclei(0x0)
,fHistTriggerStat(0x0)
,fHistTriggerStatAfterEventSelection(0x0)
,fNCounter(0)
,fOutputContainer(0x0)
,fTrigger(0)
,fTree(0x0)
,fkNTriggers(5)
{
  // Standard c-tor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  
  // cuts for candidates
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,0);
  fESDtrackCuts->SetEtaRange(-0.9f,0.9f);
  fESDtrackCuts->SetMaxDCAToVertexZ(10.f);
  fESDtrackCuts->SetMaxNOfMissingITSPoints(4);
}

//__________________________________________________________________________________________________
AliAnalysisTaskFlowd::~AliAnalysisTaskFlowd()
{
  delete fTree;
  delete fESDtrackCuts;
}

//__________________________________________________________________________________________________
void AliAnalysisTaskFlowd::BinLogAxis(const TH1 *h)
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

//__________________________________________________________________________________________________
Bool_t AliAnalysisTaskFlowd::IsTriggered() {
  // Check if Event is triggered and fill Trigger Histogram
  
  if ((fEventHandler->IsEventSelected() & AliVEvent::kMB))          fTrigger |= kMB;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kCentral))     fTrigger |= kCentral;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) fTrigger |= kSemiCentral;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      fTrigger |= kEMCEJE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      fTrigger |= kEMCEGA;
  
  for (Int_t ii = 0; ii < fkNTriggers; ++ii)
  {
    if(fTrigger & (1 << ii))
    {
      fHistTriggerStat->Fill(ii);
    }
  }
  
  return (Bool_t)fTrigger;
}

//__________________________________________________________________________________________________
void AliAnalysisTaskFlowd::ResetEvent()
{
  // Reset event
  
  fTrigger = 0;
  return;
}

//__________________________________________________________________________________________________
Int_t AliAnalysisTaskFlowd::SetupEvent()
{
  // Reset and filling of the trigger information
  ResetEvent();
  IsTriggered();
  return 0;
}

//__________________________________________________________________________________________________
void AliAnalysisTaskFlowd::UserCreateOutputObjects()
{
  // Creation of the histograms, this is called once
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);
  TString s = GetName();
  s.Append("_results");
  fOutputContainer->SetName(s.Data());
  // Histograms to count number of events
  fHistCentralityClass10  = new TH1F("fHistCentralityClass10",
                                     "Centrality with class10;Centrality;Entries",11,-0.5f,10.5f);
  //
  fHistCentralityPercentile  = new TH1F("fHistCentralityPercentile",
                                        "Centrality with percentile;Centrality;Entries",
                                        101,-0.1f,100.1f);
  
  // Trigger statitics histogram
  const Char_t* aTriggerNames[] = { "kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" };
  fHistTriggerStat = new TH1F("fHistTriggerStat","Trigger statistics",
                              fkNTriggers,-0.5f,fkNTriggers - 0.5f);
  //
  fHistTriggerStatAfterEventSelection = new TH1F("fHistTriggerStatAfterEventSelection",
                                                 "Trigger statistics after event selection",
                                                 fkNTriggers,-0.5f,fkNTriggers-0.5f);
  //
  for (Int_t ii = 0; ii < fkNTriggers; ++ii)
  {
    fHistTriggerStat->GetXaxis()->SetBinLabel(ii + 1, aTriggerNames[ii]);
    fHistTriggerStatAfterEventSelection->GetXaxis()->SetBinLabel(ii + 1, aTriggerNames[ii]);
  }
  
  // dE/dx performance
  fHistDeDx = new TH2F("fHistDeDx", "dE/dx", 3000, 0.1, 18.0, 1000, 0.0, 1000);
  fHistDeDx->GetYaxis()->SetTitle("TPC dE/dx signal (a.u.)");
  fHistDeDx->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDx);
  
  fHistDeDxITS = new TH2F("fHistDeDxITS", "dE/dx ITS", 676, 0.1, 4.0, 1000, 0.0, 1000);
  fHistDeDxITS->GetYaxis()->SetTitle("ITS dE/dx signal (a.u.)");
  fHistDeDxITS->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDxITS);
  
  // TOF performance
  fHistTOF2D = new TH2F("fHistTOF2D", "TOF2D; #frac{p}{z} (GeV/c); #beta", 500, 0.01, 10., 2250, 0.2, 1.1);
  //
  fHistTOFnuclei = new TH2F("fHistTOFnuclei","TOF; #frac{p}{z} (GeV/c); #beta",
                            500,0.01,10.,2250,0.7,1.);
  
  //alphas
  fHistDeuteron = new TH1F("fHistDeuteron", "Deuteron", 38, 1.12, 6.63);
  fHistDeuteron->GetYaxis()->SetTitle("Counts");
  fHistDeuteron->GetXaxis()->SetTitle("#frac{m^{2}}{z^{2}} (GeV^{2}/c^{4})");
  
  //
  
  fOutputContainer->Add(fHistCentralityClass10);
  fOutputContainer->Add(fHistCentralityPercentile);
  fOutputContainer->Add(fHistTriggerStat);
  fOutputContainer->Add(fHistTriggerStatAfterEventSelection);
  fOutputContainer->Add(fHistDeDx);
  fOutputContainer->Add(fHistDeDxITS);
  fOutputContainer->Add(fHistDeuteron);
  fOutputContainer->Add(fHistTOF2D);
  fOutputContainer->Add(fHistTOFnuclei);
  
  PostData(1,fOutputContainer);
  if(fFillTree) {
    OpenFile(2);
    fTree = new TTree("deuterons","Deuterons candidates");
    fTree = new TTree("Deuterons","Deuterons MC");
    fTree->SetAutoSave(100000000);
    fTree->Branch("centrality", &fTCentrality,"centrality/F");
    fTree->Branch("eta", &fTEta,"eta/F");
    fTree->Branch("phi", &fTPhi,"phi/F");
    fTree->Branch("chi2", &fTTPCchi2,"chi2/F");
    fTree->Branch("TPCsignal", &fTTPCsignal,"TPCsignal/F");
    fTree->Branch("ITSsignal", &fTITSsignal,"ITSsignal/F");
    fTree->Branch("beta", &fTBeta,"beta/F");
    fTree->Branch("DCAxy", &fTDCAxy,"DCAxy/F");
    fTree->Branch("DCAz", &fTDCAz,"DCAz/F");
    fTree->Branch("p", &fTMomentum,"p/F");
    fTree->Branch("pTPC", &fTPtpc,"pTPC/F");
    fTree->Branch("pT", &fTSignedPt,"pT/F");
    fTree->Branch("sigmaPt", &fTSigmaQoverPt,"sigmaPt/F");
    fTree->Branch("pT", &fTSignedPt,"pT/F");
    fTree->Branch("TPCnClusters", &fTTPCnClust,"TPCnClusters/s");
    fTree->Branch("TPCnSignal", &fTTPCnSignal,"TPCnSignal/s");
    fTree->Branch("ITSnClusters", &fTITSnClusters,"ITSnClusters/B");
    fTree->Branch("ITSnSignal", &fTITSnSignal,"ITSnSignal/B");
    fTree->Branch("trigger", &fTTrigger,"trigger/B");
    fTree->Branch("EventCounter",&fTEventCounter,"EventCounter/O");
    PostData(2,fTree);
  } else {
    fTree = new TTree();
  }
}

//__________________________________________________________________________________________________
void AliAnalysisTaskFlowd::UserExec(Option_t *)
{
  // Called for each event
  
  //get Event-Handler for the trigger information
  fEventHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->
                                                      GetInputEventHandler());
  if (!fEventHandler)
  {
    AliError("Could not get InputHandler");
    return;
  }
  
  fEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fEvent) {
    return;
  }
  
  if (SetupEvent() < 0) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTree);
    return;
  }
  
  const AliESDVertex *vertex = fEvent->GetPrimaryVertex();
  if(vertex->GetNContributors() < 1) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTree);
    return;
  }
  
  // check if event is selected by physics selection class
  if (TMath::Abs(vertex->GetZ()) > 10) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTree);
    return;
  }
  
  // Centrality selection in PbPb
  Int_t centralityClass10 = -1;
  fTCentrality = -1;
  AliCentrality *aodCentrality = fEvent->GetCentrality();
  // centrality percentile determined with V0
  centralityClass10 = aodCentrality->GetCentralityClass10("V0M");
  fTCentrality = aodCentrality->GetCentralityPercentile("V0M");
  //select only events with centralities between 0 and 80 %
  if (fTCentrality < 0. || fTCentrality > 80.) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTree);
    return;
  }
  
  //
  // select only events which pass kMB, kCentral, kSemiCentral
  if (!(fTrigger & kMB) && !(fTrigger & kCentral) && !(fTrigger & kSemiCentral)) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTree);
    return;
  }
  //
  fHistCentralityClass10->Fill(centralityClass10);
  fHistCentralityPercentile->Fill(fTCentrality);
  //
  if(fTrigger & kMB)          fTTrigger = 0;
  if(fTrigger & kCentral)     fTTrigger = 1;
  if(fTrigger & kSemiCentral) fTTrigger = 2;
  if(fTrigger & kEMCEJE)      fTTrigger = 3;
  if(fTrigger & kEMCEGA)      fTTrigger = 4;
  fHistTriggerStatAfterEventSelection->Fill(fTTrigger);
  //
  if (!fPid) fPid = fEventHandler->GetPIDResponse();
  if (!fPid) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTree);
    return;
  }
  
  fTEventCounter = kTRUE;
  // Track loop to fill the spectra
  for (Int_t iTracks = 0; iTracks < fEvent->GetNumberOfTracks(); iTracks++)
  {
    AliESDtrack* track = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(iTracks));
    
    // Track cuts
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    if (!track->GetInnerParam()) continue;
    if (track->GetTPCSharedMap().CountBits() > 1) continue;
    if (track->GetID() < 0) continue;
    // End of track cuts
    
    unsigned int nSPD = 0, nITS = 0;
    for (int i = 0; i < 6; ++i) {
      if (track->HasPointOnITSLayer(i)) {
        if(i < 2) nSPD++;
        nITS++;
      }
    }

    Double_t sign = track->Charge() > 0 ? 1 : -1;

    fTEta = track->Eta();
    fTPhi = track->Phi();
    fTTPCchi2 = track->GetTPCchi2() / track->GetTPCNcls();
    fTTPCsignal = track->GetTPCsignal();
    fTITSsignal = track->GetITSsignal();
    track->GetImpactParameters(fTDCAxy, fTDCAz);
    fTMomentum = track->P();
    fTPtpc = track->GetTPCmomentum();
    fTSignedPt = sign * track->Pt();
    fTSigmaQoverPt = track->GetSigma1Pt2();
    fTTPCnClust = track->GetTPCNcls();
    fTTPCnSignal = track->GetTPCsignalN();
    fTITSnClusters = nITS;
    fTITSnSignal = nITS - nSPD;
    
    Double_t expSignalDeuteron = DeuteronTPC(fTPtpc);
    
    //
    if(sign > 0)
    {
      fHistDeDx->Fill(fTPtpc, fTTPCsignal);
      if ((nITS - nSPD) > 2)
        fHistDeDxITS->Fill(fTPtpc,track->GetITSsignal());
    }
    
    //
    ULong_t  status = track->GetStatus();
    Float_t mass = 0;
    Float_t time = -1;
    Float_t gamma = -1;
    Bool_t hasTOFout  = status & AliVTrack::kTOFout;
    Bool_t hasTOFtime = status & AliVTrack::kTIME;
    Float_t length = track->GetIntegratedLength();
    Bool_t hasTOF = (hasTOFout & hasTOFtime) && length > 350.f;
    if (hasTOF == kTRUE)
    {
      Float_t time0 = fPid->GetTOFResponse().GetStartTime(fTMomentum);
      time = track->GetTOFsignal() - time0;
      if (time > 0)
      {
        fTBeta = length / (2.99792457999999984e-02 * time);
        hasTOF &= (Bool_t)(fTBeta < 1.f);
        if (fTBeta < 0.999)
        {
          gamma = 1.f / TMath::Sqrt(1.f - fTBeta * fTBeta);
          mass = fTPtpc / TMath::Sqrt(gamma * gamma - 1); // using inner TPC mom. as approx.
          if((fTTPCsignal - 0.7f * expSignalDeuteron) > 0.f)
          {
            fHistDeuteron->Fill(mass * mass);
            fHistTOFnuclei->Fill(fTPtpc,fTBeta);
          }
          
        }
        fHistTOF2D->Fill(fTPtpc,fTBeta);
      } else hasTOF = kFALSE;
    }
    
    if (!fFillTree) continue;
    if (fTPtpc > 1.6f && !hasTOF) continue;
    
    
    const bool isOKwithTPC = (fTTPCsignal - 0.7f * expSignalDeuteron > 0.f);
    
    if((isOKwithTPC && fTPtpc < 1.2f) || fTPtpc >= 1.2f)
    {
      if (!hasTOF) fTBeta = -1.f;
      fTree->Fill();
      PostData(2, fTree);
      if (fTEventCounter) fTEventCounter = kFALSE;
    }
  }//end loop over tracks
  
  // Post output data.
  PostData(1, fOutputContainer);
  if(fFillTree) PostData(2, fTree);
}

//________________________________________________________________________
void AliAnalysisTaskFlowd::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
}
