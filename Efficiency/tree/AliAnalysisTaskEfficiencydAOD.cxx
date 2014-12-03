//
//  AliAnalysisTaskEfficiencydAOD.cxx
//
//
//  Created by Maximiliano Puccio on 28/10/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#include "AliAnalysisTaskEfficiencydAOD.h"

// ROOT includes
#include <Rtypes.h>
#include <TAxis.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TParticle.h>
#include <TTree.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliPhysicsSelection.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"
#include "AliStack.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

using TMath::TwoPi;

ClassImp(AliAnalysisTaskEfficiencydAOD)

//__________________________________________________________________________________________________
AliAnalysisTaskEfficiencydAOD::AliAnalysisTaskEfficiencydAOD():
AliAnalysisTaskSE("defficiencyAOD"),						//
fPIDResponse(),
fTree(),
fTpMC(0.f),
fTpTMC(0.f),
fTetaMC(0.f),
fTyMC(0.f),
fTphiMC(0.f),
fTp(0.f),
fTpT(0.f),
fTeta(0.f),
fTphi(0.f),
fTbeta(0.f),
fTDCAxy(0.f),
fTDCAz(0.f),
fTchi2(-1.f),
fTcentrality(-1.f),
fTTPCnClusters(0u),
fTTPCnSignal(0u),
fTITSnClusters(0),
fTITSnSignal(0),
fTIsPrimary(false),
fTIsSecondaryFromMaterial(false)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
}

//__________________________________________________________________________________________________
AliAnalysisTaskEfficiencydAOD::~AliAnalysisTaskEfficiencydAOD(){
  
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  // Destructor
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEfficiencydAOD::UserCreateOutputObjects(){
  OpenFile(1);
  fTree = new TTree("Deuterons","Deuterons MC");
  fTree->SetAutoSave(100000000);
  fTree->Branch("pMC", &fTpMC,"pMC/F");
  fTree->Branch("pTMC", &fTpTMC,"pTMC/F");
  fTree->Branch("etaMC", &fTetaMC,"etaMC/F");
  fTree->Branch("phiMC", &fTphiMC,"phiMC/F");
  fTree->Branch("yMC", &fTyMC,"yMC/F");
  fTree->Branch("p", &fTp,"p/F");
  fTree->Branch("pT", &fTpT,"pT/F");
  fTree->Branch("eta", &fTeta,"eta/F");
  fTree->Branch("phi", &fTphi,"phi/F");
  fTree->Branch("beta", &fTbeta,"beta/F");
  fTree->Branch("DCAxy", &fTDCAxy,"DCAxy/F");
  fTree->Branch("DCAz", &fTDCAz,"DCAz/F");
  fTree->Branch("chi2", &fTchi2,"chi2/F");
  fTree->Branch("centrality", &fTcentrality,"centrality/F");
  fTree->Branch("TPCnClusters", &fTTPCnClusters,"TPCnClusters/s");
  fTree->Branch("TPCnSignal", &fTTPCnSignal,"TPCnSignal/s");
  fTree->Branch("ITSnClusters", &fTITSnClusters,"ITSnClusters/B");
  fTree->Branch("ITSnSignal", &fTITSnSignal,"ITSnSignal/B");
  fTree->Branch("IsPrimary", &fTIsPrimary,"IsPrimary/O");
  fTree->Branch("IsSecondaryFromMaterial",&fTIsSecondaryFromMaterial,"IsSecondaryFromMaterial/O");
  PostData(1,fTree);
  
  
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEfficiencydAOD::UserExec(Option_t *){
  AliVEvent *ev = dynamic_cast<AliVEvent*> (InputEvent());
  
  // Check event selection mask
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t mask = handl->IsEventSelected();
  if (!(mask & 0xffffffff)) {
    PostData(1, fTree);
    return;
  }
  if (!((mask & AliVEvent::kMB) || (mask & AliVEvent::kCentral) || (mask & AliVEvent::kSemiCentral))) {
    PostData(1, fTree);
    return;
  }
  
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  AliStack* stack = mcEvent->Stack();
  
  AliCentrality *centr=ev->GetCentrality();
  
  if (!fPIDResponse)
  {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->
                                            GetInputEventHandler()))->GetPIDResponse();
  }
  if (!fPIDResponse)
  {
    PostData(1, fTree);
    return;
  }
  
  // Centrality selection in PbPb, percentile determined with V0
  fTcentrality = centr->GetCentralityPercentile("V0M");
  //select only events with centralities between 0 and 80 %
  if (fTcentrality < 0. || fTcentrality > 80.) {
    PostData(1, fTree);
    return;
  }
  // Primary vertex displacement cut
  const AliVVertex *primaryVtx = ev->GetPrimaryVertex();
  if(TMath::Abs(primaryVtx->GetZ()) > 10.) {
    PostData(1, fTree);
    return;
  }
  
  // Making the list of deuterons in acceptance
  TList mcD, mcDbar;
  mcD.SetOwner(kFALSE);
  mcDbar.SetOwner(kFALSE);
  for (int iMC = 0; iMC < stack->GetNtrack(); ++iMC) {
    TParticle *part = stack->Particle(iMC);
    const int pdg = part->GetPdgCode();
    if (stack->IsPhysicalPrimary(iMC))
      part->SetFirstMother(-1);
    else if (stack->IsSecondaryFromMaterial(iMC))
      part->SetFirstMother(-2);
    else
      part->SetFirstMother(0);
    if (pdg == 1000010020) {
      mcD.Add(part);
    } else if (pdg == -1000010020) {
      mcDbar.Add(part);
    }
  }
  
  // Checking how many deuterons in acceptance are reconstructed well
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(ev->GetTrack(iT));
    // Skip fake tracks
    if (track->GetLabel() < 0) continue;
    if (track->GetID() <= 0) continue;
    // Track cuts
    ULong_t status = track->GetStatus();
    if (!(status & AliVTrack::kITSrefit)) continue;
    if (!(status & AliVTrack::kTPCrefit)) continue;
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
    // Has the track a nice prolongation in TOF? Answer in hasTOF.
    Bool_t hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
    Bool_t hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
    Float_t length = track->GetIntegratedLength();
    Bool_t hasTOF = Bool_t(hasTOFout & hasTOFtime) && length > 350.f;
    
    // Getting the particle and checking if it is in one of the two lists of MC particles
    TParticle *part = stack->Particle(TMath::Abs(track->GetLabel()));
    if (!part) continue;
    
    int isDeuteron = 2;
    if (mcD.Contains(part)) {
      isDeuteron = 1;
    } else if (mcDbar.Contains(part)) {
      isDeuteron = -1;
    }
    if (isDeuteron == 2) continue;
    
    fTpMC = part->P();
    fTpTMC = isDeuteron * part->Pt();
    fTetaMC = part->Eta();
    fTphiMC = part->Phi();
    fTyMC = part->Y();
    fTp = track->P();
    fTpT = track->Charge() * track->Pt() / TMath::Abs(track->Charge());
    fTeta = track->Eta();
    fTphi = track->Phi();
    if (hasTOF) {
      const float len = track->GetIntegratedLength();
      const float tim = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(fTp);
      fTbeta = len / (2.99792457999999984e-02 * tim);
    } else
      fTbeta = -1.f;
    track->GetImpactParameters(fTDCAxy, fTDCAz);
    fTchi2 = track->Chi2perNDF();
    fTTPCnClusters = track->GetTPCNcls();
    fTTPCnSignal = track->GetTPCsignalN();
    fTITSnClusters = track->GetNcls(0);
    fTITSnSignal = nITS - nSPD;
    fTIsPrimary = Bool_t(part->GetFirstMother() == -1);
    fTIsSecondaryFromMaterial = Bool_t(part->GetFirstMother() == -2);
    
    fTree->Fill();
    if (isDeuteron == 1)
      mcD.Remove(part);
    else
      mcDbar.Remove(part);
    
  } // End AOD track loop
  
  fTp = -1.f;
  fTpT = 0.f;
  fTeta = 0.f;
  fTphi = 0.f;
  fTbeta = -1.f;
  fTDCAxy = -999.f;
  fTDCAz = -999.f;
  fTchi2 = -1.f;
  fTTPCnClusters = 0;
  fTTPCnSignal = 0;
  fTITSnClusters = 0;
  fTITSnSignal = 0;
  
  TListIter nextD(&mcD);
  TParticle *part;
  while ((part = (TParticle*)nextD())) {
    fTpMC = part->P();
    fTpTMC = part->Pt();
    fTetaMC = part->Eta();
    fTphiMC = part->Phi();
    fTyMC = part->Y();
    fTIsPrimary = Bool_t(part->GetFirstMother() == -1);
    fTIsSecondaryFromMaterial = Bool_t(part->GetFirstMother() == -2);
    fTree->Fill();
  }
  
  TListIter nextDbar(&mcDbar);
  while ((part = (TParticle*)nextDbar())) {
    fTpMC = part->P();
    fTpTMC = -part->Pt();
    fTetaMC = part->Eta();
    fTphiMC = part->Phi();
    fTyMC = part->Y();
    fTIsPrimary = Bool_t(part->GetFirstMother() == -1);
    fTIsSecondaryFromMaterial = Bool_t(part->GetFirstMother() == -2);
    fTree->Fill();
  }
  //  Post output data.
  PostData(1,fTree);
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEfficiencydAOD::Terminate(Option_t *) {
  // Merge output
  // Called once at the end of the query
  
  return;
}
