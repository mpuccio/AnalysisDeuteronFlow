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
#include "THn.h"
#include "TCanvas.h"
#include "AliAODpidUtil.h"
#include "AliAODVertex.h"
#include "TFile.h"
#include "TString.h"
#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliTOFPIDResponse.h"
#include "AliExternalTrackParam.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

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
,fAOD(0x0)
,fAODpid(0x0)
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
,fTTree(0x0)
,fTCentrality(0)
,fTEta(0)
,fTTPCsignal(0)
,fTchi2NDF(0)
,fTITSsignal(0)
,fTTOFtime(0)
,fTDCAxy(0)
,fTDCAz(0)
,fTp(0)
,fTpTPC(0)
,fTpT(0)
,fTLength(0)
,fTTPCsigmad(0)
,fTTPCsigmat(0)
,fTFilterMap(0)
,fTITSnClust(0)
,fTITSnSignal(0)
,fTTPCnClust(0)
,fTTPCnClustShared(0)
,fTTPCnSignal(0)
,fTTrigger(0)
,fkNTriggers(5)
{
  // Standard c-tor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());


}

//__________________________________________________________________________________________________
AliAnalysisTaskFlowd::~AliAnalysisTaskFlowd()
{
  if (fTTree) delete fTTree;
  
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
  
  if ((fEventHandler->IsEventSelected() & AliVEvent::kMB))          fTTrigger |= kMB;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kCentral))     fTTrigger |= kCentral;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) fTTrigger |= kSemiCentral;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      fTTrigger |= kEMCEJE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      fTTrigger |= kEMCEGA;
  
  for (Int_t ii = 0; ii < fkNTriggers; ++ii)
  {
    if(fTTrigger & (1 << ii))
    {
      fHistTriggerStat->Fill(ii);
    }
  }
  
  return (Bool_t)fTTrigger;
}

//__________________________________________________________________________________________________
void AliAnalysisTaskFlowd::ResetEvent()
{
  // Reset event
  
  fTTrigger = 0;
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
    fTTree = new TTree("deuterons","deuteron candidates");
    fTTree->Branch("centrality", &fTCentrality);
    fTTree->Branch("eta",&fTEta);
    fTTree->Branch("TPCsignal",&fTTPCsignal);
    fTTree->Branch("chi2NDF", &fTchi2NDF);
    fTTree->Branch("ITSsignal", &fTITSsignal);
    fTTree->Branch("TOFtime", &fTTOFtime);
    fTTree->Branch("DCAxy", &fTDCAxy);
    fTTree->Branch("DCAz", &fTDCAz);
    fTTree->Branch("p", &fTp);
    fTTree->Branch("pTPC", &fTpTPC);
    fTTree->Branch("pT", &fTpT);
    fTTree->Branch("length", &fTLength);
    fTTree->Branch("TPCsigmad", &fTTPCsigmad);
    fTTree->Branch("TPCsigmat", &fTTPCsigmat);
    fTTree->Branch("FilterMap", &fTFilterMap);
    fTTree->Branch("ITSnClust", &fTITSnClust);
    fTTree->Branch("ITSnSignal", &fTITSnSignal);
    fTTree->Branch("TPCnClust", &fTTPCnClust);
    fTTree->Branch("TPCnClustShared", &fTTPCnClustShared);
    fTTree->Branch("TPCnSignal", &fTTPCnSignal);
    fTTree->Branch("trigger", &fTTrigger);

    fTTree->SetAutoSave(100000000);
    PostData(2,fTTree);
  } else {
    fTTree = new TTree();
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
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    return;
  }
  
  if (SetupEvent() < 0) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTTree);
    return;
  }
  
  const AliAODVertex *vertex = fAOD->GetPrimaryVertex();
  if(vertex->GetNContributors() < 1) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTTree);
    return;
  }
  
  // check if event is selected by physics selection class
  Bool_t isSelected = kFALSE;
  isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->
                                        GetInputEventHandler()))->IsEventSelected();
  if (!isSelected || TMath::Abs(vertex->GetZ()) > 10) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTTree);
    return;
  }
  
  // Centrality selection in PbPb
  Int_t centralityClass10 = -1;
  AliCentrality *aodCentrality = fAOD->GetCentrality();
  // centrality percentile determined with V0
  centralityClass10 = aodCentrality->GetCentralityClass10("V0M");
  fTCentrality = aodCentrality->GetCentralityPercentile("V0M");
  //select only events with centralities between 0 and 80 %
  if (fTCentrality < 0. || fTCentrality > 80.) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTTree);
    return;
  }
  
  //
  // select only events which pass kMB, kCentral, kSemiCentral
  if (!(fTTrigger & kMB) && !(fTTrigger & kCentral) && !(fTTrigger & kSemiCentral)) {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTTree);
    return;
  }
  //
  fHistCentralityClass10->Fill(centralityClass10);
  fHistCentralityPercentile->Fill(fTCentrality);
  //
  if(fTTrigger & kMB)          fHistTriggerStatAfterEventSelection->Fill(0);
  if(fTTrigger & kCentral)     fHistTriggerStatAfterEventSelection->Fill(1);
  if(fTTrigger & kSemiCentral) fHistTriggerStatAfterEventSelection->Fill(2);
  if(fTTrigger & kEMCEJE)      fHistTriggerStatAfterEventSelection->Fill(3);
  if(fTTrigger & kEMCEGA)      fHistTriggerStatAfterEventSelection->Fill(4);
  //
  if (!fAODpid)
  {
    fAODpid = ((AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->
                                     GetInputEventHandler()))->GetAODpidUtil();
  }
  if (!fAODpid)
  {
    PostData(1, fOutputContainer);
    if(fFillTree) PostData(2, fTTree);
    return;
  }

  bool first = true;
  // Track loop to fill the spectra
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++)
  {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
    ULong_t  status = track->GetStatus();
    Double_t ptot = track->GetTPCmomentum();
    
    // Track cuts
    fTTPCnClustShared = track->GetTPCnclsS();
    if (!(status & AliVTrack::kITSrefit)) continue;
    if (!(status & AliVTrack::kTPCrefit)) continue;
    if (track->GetID() < 0) continue;
    if (track->GetTPCsignalN() < 70) continue;
    if (track->GetTPCNcls() < 70) continue;
    if (TMath::Abs(track->Eta()) > 0.8f) continue;
    if (track->Chi2perNDF() > 4.f) continue;
    AliAODVertex *vtx1 = (AliAODVertex*)track->GetProdVertex();
    if(Int_t(vtx1->GetType()) == AliAODVertex::kKink) continue;
    unsigned int nSPD = 0, nITS = 0;
    for (int i = 0; i < 6; ++i) {
      if (track->HasPointOnITSLayer(i)) {
        if(i < 2) nSPD++;
        nITS++;
      }
    }
    if (nITS < 2) continue;
    Double_t dca[2],cov[3];
    if (!track->PropagateToDCA(vertex, fAOD->GetMagneticField(), 100, dca, cov)) continue;
    if (TMath::Abs(dca[1]) > 2.f) continue;
    // End of track cuts
    
    Double_t sign = track->Charge() > 0 ? 1 : -1;
    Double_t expSignalDeuteron = DeuteronTPC(ptot);
    
    //
    if(sign > 0)
    {
      fHistDeDx->Fill(ptot, track->GetTPCsignal());
      if ((nITS - nSPD) > 2)
        fHistDeDxITS->Fill(ptot,track->GetITSsignal());
    }
    
    //
    Float_t mass = 0;
    Float_t time = -1;
    Float_t beta = 0;
    Float_t gamma = -1;
    Bool_t hasTOFout  = status & AliVTrack::kTOFout;
    Bool_t hasTOFtime = status & AliVTrack::kTIME;
    Float_t length = track->GetIntegratedLength();
    Bool_t hasTOF = (hasTOFout & hasTOFtime) && length > 350.f;
    if (hasTOF == kTRUE)
    {
      Float_t time0 = fAODpid->GetTOFResponse().GetStartTime(track->P());
      time = track->GetTOFsignal() - time0;
      if (time > 0)
      {
        beta = length / (2.99792457999999984e-02 * time);
        hasTOF &= (Bool_t)(beta < 1.f);
        if (beta < 0.975)
        {
          gamma = 1.f / TMath::Sqrt(1.f - beta * beta);
          mass = ptot / TMath::Sqrt(gamma * gamma - 1); // using inner TPC mom. as approx.
          if((track->GetTPCsignal() - expSignalDeuteron) / expSignalDeuteron > -0.3)
          {
            fHistDeuteron->Fill(mass * mass);
            fHistTOFnuclei->Fill(ptot,beta);
          }
          
        }
        fHistTOF2D->Fill(ptot,beta);
      } else hasTOF = kFALSE;
    }
    
    if (!fFillTree)
      continue;
    
    if (ptot > 1.2f && !hasTOF)
      continue;
    
    
    if((fAODpid->NumberOfSigmasTPC(track,AliPID::kProton) >= -3.f && ptot < 1.2f) || ptot >= 1.2f)
    {
      fTEta        = track->Eta();                                                  // eta
      fTTPCnClust  = track->GetTPCNcls();                                           // TPCnClust
      fTTPCsignal  = track->GetTPCsignal();                                         // TPCsignal
      fTTPCnSignal = track->GetTPCsignalN();                                        // TPCnSignal
      fTchi2NDF    = track->Chi2perNDF();                                           // TPCchi2
      fTITSsignal  = track->GetITSsignal();                                         // ITSsignal
      fTITSnClust  = nITS;                                                          // ITSnClust
      fTITSnSignal = nITS - nSPD;                                                   // ITSnClustPID
      fTTOFtime    = hasTOF ? time : -1.f;                                          // TOFtime
      fTDCAxy      = dca[0];                                                        // DCAxy
      fTDCAz       = dca[1];                                                        // DCAz
      fTp          = track->P();                                                    // p
      fTpTPC       = ptot;                                                          // pTPC
      fTpT         = sign * track->Pt();                                            // pT
      fTLength     = length;                                                        // length
      fTTPCsigmad  = fAODpid->NumberOfSigmasTPC(track, AliPID::kDeuteron);          // TPCsigmad
      fTTPCsigmat  = fAODpid->NumberOfSigmasTPC(track, AliPID::kTriton);            // TPCsigmat
      fTFilterMap  = track->GetFilterMap();
      if (first) {
        fTCentrality = -fTCentrality;
        fTTree->Fill();
        fTCentrality = -fTCentrality;
      } else
        fTTree->Fill();
      PostData(2, fTTree);
    }
  }//end loop over tracks
  
  // Post output data.
  PostData(1, fOutputContainer);
  if(fFillTree) PostData(2, fTTree);
}

//________________________________________________________________________
void AliAnalysisTaskFlowd::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
}
