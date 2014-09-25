//
//  AliAnalysisTaskFlowd.cxx
//
//  Created by Maximiliano Puccio on 14/09/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#include "AliAnalysisTaskFlowd.h"
#include <Riostream.h>
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TCanvas.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "TFile.h"
#include "TString.h"
#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskFlowd)

//__________________________________________________________________________________________________
inline static Double_t DeuteronTPC(Double_t x) {
  // Deuteron expected signal in TPC
  return AliExternalTrackParam::BetheBlochAleph(x / 1.876,1.45802,27.4992,4.00313e-15,
                                                2.48485,8.31768);
}

//__________________________________________________________________________________________________
inline static Int_t NumberOfPIDClustersITS(AliESDtrack* track) {
  Char_t cmap = track->GetITSClusterMap();
  Int_t nclPID = 0;
  for (int i = 2; i < 6; ++i) {
    if (cmap & (1 << i)) {
      ++nclPID;
    }
  }
  return nclPID;
}

//__________________________________________________________________________________________________
AliAnalysisTaskFlowd::AliAnalysisTaskFlowd(const char* name)
:AliAnalysisTaskSE(name)
,fESD(0x0)
,fESDpid(0x0)
,fESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts")
,fESDtrackCutsStrict("AliESDtrackCuts","AliESDtrackCuts")
,fEventHandler(0x0)
,fHistCentralityClass10(0x0)
,fHistCentralityPercentile(0x0)
,fHistDeDx(0x0)
,fHistDeDxITS(0x0)
,fHistDeDxITSsa(0x0)
,fHistDeuteron(0x0)
,fHistTOF2D(0x0)
,fHistTOFnuclei(0x0)
,fHistTriggerStat(0x0)
,fHistTriggerStatAfterEventSelection(0x0)
,fNCounter(0)
,fOutputContainer(0x0)
,fTrigger(0)
,fkNTriggers(5)
{
  // Standard c-tor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
//  DefineOutput(2, TTree::Class());
  
  // cuts for candidates
  //
  fESDtrackCuts.SetAcceptKinkDaughters(kFALSE);
//  fESDtrackCuts.SetMinNClustersTPC(70);
//  fESDtrackCuts.SetMaxChi2PerClusterTPC(6);
  fESDtrackCuts.SetMaxDCAToVertexXY(3);
  fESDtrackCuts.SetMaxDCAToVertexZ(2);
  //fESDtrackCuts.SetRequireTPCRefit(kTRUE);
  fESDtrackCuts.SetRequireITSRefit(kTRUE);
  fESDtrackCuts.SetMinNClustersITS(1);
  fESDtrackCuts.SetEtaRange(-1.0,1.0);
  
  fESDtrackCutsStrict.SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsStrict.SetMinNClustersTPC(70);
  fESDtrackCutsStrict.SetMaxChi2PerClusterTPC(6);
  fESDtrackCutsStrict.SetMaxDCAToVertexXY(3);
  fESDtrackCutsStrict.SetMaxDCAToVertexZ(2);
  fESDtrackCutsStrict.SetRequireTPCRefit(kTRUE);
  fESDtrackCutsStrict.SetRequireITSRefit(kTRUE);
  fESDtrackCutsStrict.SetMinNClustersITS(1);
  fESDtrackCutsStrict.SetEtaRange(-1.0,1.0);

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
void AliAnalysisTaskFlowd::BinLogAxis(const TH3 *h, Int_t axisNumber)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = const_cast<TAxis*>(h->GetXaxis());
  if (axisNumber == 1) axis = const_cast<TAxis*>(h->GetYaxis());
  if (axisNumber == 2) axis = const_cast<TAxis*>(h->GetZaxis());
  const Int_t bins = axis->GetNbins();
  
  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
  
  newBins[0] = from;
  const Double_t factor = pow(to / from, 1. / bins);
  
  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
}

//__________________________________________________________________________________________________
Int_t AliAnalysisTaskFlowd::Initialize()
{
  // Initialisation === Reset
  ResetEvent();
  return 0;
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
  fHistDeDx = new TH2F("fHistDeDx", "dE/dx", 1000, 0.1, 6.0, 1000, 0.0, 1000);
  fHistDeDx->GetYaxis()->SetTitle("TPC dE/dx signal (a.u.)");
  fHistDeDx->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDx);
  
  fHistDeDxITS = new TH2F("fHistDeDxITS", "dE/dx ITS", 677, 0.01, 4.0, 1000, 0.0, 1000);
  fHistDeDxITS->GetYaxis()->SetTitle("ITS dE/dx signal (a.u.)");
  fHistDeDxITS->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDxITS);

  
  fHistDeDxITSsa = new TH2F("fHistDeDxITSsa", "dE/dx ITSsa", 677, 0.01, 4.0, 1000, 0.0, 1000);
  fHistDeDxITSsa->GetYaxis()->SetTitle("ITS dE/dx signal (a.u.)");
  fHistDeDxITSsa->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDxITSsa);

  
  // TOF performance
  fHistTOF2D = new TH2F("fHistTOF2D", "TOF2D; #frac{p}{z} (GeV/c); #beta", 500, 0.0, 10., 2250, 0.2, 1.1);
  //
  fHistTOFnuclei = new TH2F("fHistTOFnuclei","TOF; #frac{p}{z} (GeV/c); #beta",
                            500,0.,10.,2250,0.7,1.);

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
  fOutputContainer->Add(fHistDeDxITSsa);
  fOutputContainer->Add(fHistDeuteron);
  fOutputContainer->Add(fHistTOF2D);
  fOutputContainer->Add(fHistTOFnuclei);
  PostData(1,fOutputContainer);
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
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    return;
  }
  
  if (SetupEvent() < 0) {
    PostData(1, fOutputContainer);
    return;
  }
  
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors() < 1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors() < 1) {
      PostData(1, fOutputContainer);
      return;
    }
  }
  
  // check if event is selected by physics selection class
  Bool_t isSelected = kFALSE;
  isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->
                                        GetInputEventHandler()))->IsEventSelected();
  if (!isSelected || TMath::Abs(vertex->GetZv()) > 10) {
    PostData(1, fOutputContainer);
    return;
  }
  
  // Centrality selection in PbPb
  Int_t centralityClass10 = -1;
  Float_t centralityPercentile = -1;
  //
  if (fESD->GetEventSpecie() == 4) { //species == 4 == PbPb
    // PbPb
    AliCentrality *esdCentrality = fESD->GetCentrality();
    // centrality percentile determined with V0
    centralityClass10 = esdCentrality->GetCentralityClass10("V0M");
    centralityPercentile = esdCentrality->GetCentralityPercentile("V0M");
    //select only events with centralities between 0 and 80 %
    if (centralityPercentile < 0. || centralityPercentile > 80. )
      return;
  }
  //
  // select only events which pass kMB, kCentral, kSemiCentral
  if (!(fTrigger & kMB) && !(fTrigger & kCentral) && !(fTrigger & kSemiCentral))
    return;
  //
  fHistCentralityClass10->Fill(centralityClass10);
  fHistCentralityPercentile->Fill(centralityPercentile);
  //
  if(fTrigger & kMB )          fHistTriggerStatAfterEventSelection->Fill(0);
  if(fTrigger & kCentral )     fHistTriggerStatAfterEventSelection->Fill(1);
  if(fTrigger & kSemiCentral ) fHistTriggerStatAfterEventSelection->Fill(2);
  if(fTrigger & kEMCEJE )      fHistTriggerStatAfterEventSelection->Fill(3);
  if(fTrigger & kEMCEGA )      fHistTriggerStatAfterEventSelection->Fill(4);
  //
  if (!fESDpid)
  {
    fESDpid = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->
                                     GetInputEventHandler()))->GetESDpid();
  }
  if (!fESDpid)
  {
    fESDpid = new AliESDpid(); // HACK FOR MC PBPB --> PLEASE REMOVE AS SOON AS POSSIBLE
    fESDpid->GetTPCResponse().SetBetheBlochParameters(1.28778e+00 / 50., 3.13539e+01,
                                                      TMath::Exp(-3.16327e+01), 1.87901e+00,
                                                      6.41583e+00);
  }
    //
//  Bool_t fillTree = kFALSE;
  // Track loop to fill the spectram
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++)
  {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!fESDtrackCuts.AcceptTrack(track)) continue;
    UInt_t  status = track->GetStatus();
    Double_t ptot = track->GetInnerParam()->GetP();
    Double_t ptotInc = track->GetP(); // total momentum of the incoming particle
    Double_t sign = track->GetSign();

    Double_t expSignalDeuteron = DeuteronTPC(ptot);

    // fill final histograms
    if (NumberOfPIDClustersITS(track) > 2 && !(status & AliVTrack::kTPCrefit)) {
      fHistDeDxITSsa->Fill(ptotInc,track->GetITSsignal());
    }
    if(!fESDtrackCutsStrict.AcceptTrack(track))
      continue;
    Double_t nClustersTPCPID = track->GetTPCsignalN();
    TBits shared = track->GetTPCSharedMap();
    if (shared.CountBits() > 1 ||
        nClustersTPCPID < 80 || track->GetKinkIndex(0) != 0 || track->GetTPCNclsIter1() < 80)
      continue;
    
    //
    if(sign > 0) {
      fHistDeDx->Fill(ptot, track->GetTPCsignal());
      if (NumberOfPIDClustersITS(track) > 2)
        fHistDeDxITS->Fill(ptot,track->GetITSsignal());
    }
    
    //
    Float_t mass = 0;
    Float_t time = -1;
    Float_t beta = 0;
    Float_t gamma = -1;
    Bool_t hasTOFout  = status & AliESDtrack::kTOFout;
    Bool_t hasTOFtime = status & AliESDtrack::kTIME;
    Float_t length = track->GetIntegratedLength();
    Bool_t hasTOF = (hasTOFout & hasTOFtime) && length > 350.f;
    if (hasTOF == kTRUE && ptot < 5) {
      Float_t time0 = fESDpid->GetTOFResponse().GetStartTime(track->P());//fESDpid->GetTOFResponse().GetTimeZero();
      time = track->GetTOFsignal() - time0;
      if (time > 0) {
        beta = length / (2.99792457999999984e-02 * time);
        if (beta < 0.975) {
          gamma = 1/TMath::Sqrt(1 - beta*beta);
          mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
          if (TMath::Sqrt(track->GetTOFsignalDz() * track->GetTOFsignalDz() +
                          track->GetTOFsignalDx()*track->GetTOFsignalDx()) < 5.) {
            if((track->GetTPCsignal()-expSignalDeuteron)/expSignalDeuteron > -0.3 &&
               (track->GetTPCsignal()-expSignalDeuteron)/expSignalDeuteron < 0.3) {
              fHistDeuteron->Fill(mass*mass);
              fHistTOFnuclei->Fill(ptot,beta);
            }
          }
        }
        fHistTOF2D->Fill(ptot,beta);
      }
    }
    
    
  }//end loop over tracks
  
  
  // Post output data.
  PostData(1, fOutputContainer);
}      

//________________________________________________________________________
void AliAnalysisTaskFlowd::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
  //if (!GetOutputData(0)) return;
}
