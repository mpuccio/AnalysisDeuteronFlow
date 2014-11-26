//
//  AliAnalysisTaskEfficiencyd.cxx
//
//
//  Created by Maximiliano Puccio on 28/10/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#include "AliAnalysisTaskEfficiencyd.h"

// ROOT includes
#include <Rtypes.h>
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TTree.h>

// ALIROOT includes
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDtrack.h"
#include "AliPhysicsSelection.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

using TMath::TwoPi;

ClassImp(AliAnalysisTaskEfficiencyd)

//__________________________________________________________________________________________________
AliAnalysisTaskEfficiencyd::AliAnalysisTaskEfficiencyd():
AliAnalysisTaskSE("Task_d_efficiency"),						//
fAOD(),
fAODpid(),
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
fTIsPrimary(false)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
}

//__________________________________________________________________________________________________
AliAnalysisTaskEfficiencyd::~AliAnalysisTaskEfficiencyd(){
  
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  // Destructor
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEfficiencyd::UserCreateOutputObjects(){
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
  PostData(1,fTree);
  

}

//__________________________________________________________________________________________________
void AliAnalysisTaskEfficiencyd::UserExec(Option_t *){
  AliAODEvent *ev = dynamic_cast<AliAODEvent*> (InputEvent());
  
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
  
  // get branch "mcparticles"
  TClonesArray *stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!stack) {
    PostData(1, fTree);
    return;
  }
  
  // Centrality quality (do we need this with AODs?)
  AliCentrality *centr=ev->GetCentrality();
  if (centr->GetQuality() != 0) {
    PostData(1, fTree);
    return;
  }
  
  if (!fAODpid)
  {
    fAODpid = ((AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->
                                     GetInputEventHandler()))->GetAODpidUtil();
  }
  if (!fAODpid)
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
  
  // Used later for the track propagation
  Double_t b = (Double_t)ev->GetMagneticField();
  Double_t dca_tr[2] = {0.,0.}, covar_tr[2] = {0.,0.};
  
  // Making the list of deuterons in acceptance
  TList mcD, mcDbar;
  mcD.SetOwner(kFALSE);
  mcDbar.SetOwner(kFALSE);
  for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {
    AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
    const int pdg = part->GetPdgCode();
    if (pdg == 1000010020) {
      mcD.Add(part);
    } else if (pdg == -1000010020) {
      mcDbar.Add(part);
    }
  }
  
  // Checking how many deuterons in acceptance are reconstructed well
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliAODTrack *aodtr = dynamic_cast<AliAODTrack*>(ev->GetTrack(iT));
    // Skip fake tracks
    if (aodtr->GetLabel() < 0) continue;
    if (aodtr->GetID() <= 0) continue;
    // Track cuts
    ULong_t status = aodtr->GetStatus();
    if (!(status & AliVTrack::kITSrefit)) continue;
    if (!(status & AliVTrack::kTPCrefit)) continue;
    AliAODVertex *vtx1 = (AliAODVertex*)aodtr->GetProdVertex();
    if(Int_t(vtx1->GetType()) == AliAODVertex::kKink) continue;
    unsigned int nSPD = 0, nITS = 0;
    for (int i = 0; i < 6; ++i) {
      if (aodtr->HasPointOnITSLayer(i)) {
        if(i < 2) nSPD++;
        nITS++;
      }
    }
    if (nITS < 2) continue;
    if (!aodtr->PropagateToDCA(primaryVtx, b, 100, dca_tr, covar_tr)) continue;
    // Has the track a nice prolongation in TOF? Answer in hasTOF.
    Bool_t hasTOFout  = aodtr->GetStatus() & AliESDtrack::kTOFout;
    Bool_t hasTOFtime = aodtr->GetStatus() & AliESDtrack::kTIME;
    Float_t length = aodtr->GetIntegratedLength();
    Bool_t hasTOF = Bool_t(hasTOFout & hasTOFtime) && length > 350.f;
    
    // Getting the particle and checking if it is in one of the two lists of MC particles
    AliAODMCParticle *part = (AliAODMCParticle*) stack->At(aodtr->GetLabel());
    if (!part) continue;
    int isDeuteron = -1;
    if (mcD.Contains(part)) {
      isDeuteron = 1;
    } else if (mcDbar.Contains(part)) {
      isDeuteron = 2;
    }
    if (isDeuteron == -1) continue;
    
    fTpMC = part->P();
    fTpTMC = part->Charge() * part->Pt() / TMath::Abs(part->Charge());
    fTetaMC = part->Eta();
    fTphiMC = part->Phi();
    fTyMC = part->Y();
    fTp = aodtr->P();
    fTpT = aodtr->Charge() * aodtr->Pt() / TMath::Abs(aodtr->Charge());
    fTeta = aodtr->Eta();
    fTphi = aodtr->Phi();
    if (hasTOF) {
      const float len = aodtr->GetIntegratedLength();
      const float tim = aodtr->GetTOFsignal() - fAODpid->GetTOFResponse().GetStartTime(fTp);
      fTbeta = len / (2.99792457999999984e-02 * tim);
    } else
      fTbeta = -1.f;
    fTDCAxy = dca_tr[0];
    fTDCAz = dca_tr[1];
    fTchi2 = aodtr->Chi2perNDF();
    fTTPCnClusters = aodtr->GetTPCNcls();
    fTTPCnSignal = aodtr->GetTPCsignalN();
    fTITSnClusters = aodtr->GetITSNcls();
    fTITSnSignal = nITS - nSPD;
    fTIsPrimary = part->IsPhysicalPrimary();
    
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
  AliAODMCParticle *part;
  while ((part = (AliAODMCParticle*)nextD())) {
    fTpMC = part->P();
    fTpTMC = part->Charge() * part->Pt() / TMath::Abs(part->Charge());
    fTetaMC = part->Eta();
    fTphiMC = part->Phi();
    fTyMC = part->Y();
    fTIsPrimary = part->IsPhysicalPrimary();
    fTree->Fill();
  }
  
  TListIter nextDbar(&mcDbar);
  while ((part = (AliAODMCParticle*)nextDbar())) {
    fTpMC = part->P();
    fTpTMC = part->Charge() * part->Pt() / TMath::Abs(part->Charge());
    fTetaMC = part->Eta();
    fTphiMC = part->Phi();
    fTyMC = part->Y();
    fTIsPrimary = part->IsPhysicalPrimary();
    fTree->Fill();
  }
  //  Post output data.
  PostData(1,fTree);
}

//__________________________________________________________________________________________________
void AliAnalysisTaskEfficiencyd::Terminate(Option_t *) {
  // Merge output
  // Called once at the end of the query
  
  return;
}
