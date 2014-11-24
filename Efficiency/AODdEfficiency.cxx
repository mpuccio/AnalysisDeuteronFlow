//
//  AODdEfficiency.cxx
//
//
//  Created by Maximiliano Puccio on 28/10/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#include "AODdEfficiency.h"

// ROOT includes
#include <Rtypes.h>
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

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
#include "AliVTrack.h"
#include "AliVVertex.h"

using TMath::TwoPi;

ClassImp(AODdEfficiency)

//__________________________________________________________________________________________________
AODdEfficiency::AODdEfficiency():
AliAnalysisTaskSE("Task d efficiency"),						//
fAOD(0),
fOutput(0),
fAntiDMCYield(0),
fAntiDYield(0),
fAntiDYieldTOF(0),
fDCAxyPrimariesAD(),
fDCAxySecondariesAD(),
fDCAzPrimariesAD(),
fDCAzSecondariesAD(),
fDCAxyPrimariesD(),
fDCAxySecondariesD(),
fDCAzPrimariesD(),
fDCAzSecondariesD(),
fDMCYield(0),
fDYield(0),
fDYieldTOF(0),
fEtaPhiCoverage(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//__________________________________________________________________________________________________
AODdEfficiency::~AODdEfficiency(){
  
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//__________________________________________________________________________________________________
void AODdEfficiency::UserCreateOutputObjects(){
  
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("dEfficiency");
  
  fDYield = new TH1F("fDYield",";p_{T} (GeV/c);Number of d tracks",120,0.1,6.1);
  fAntiDYield = new TH1F("fAntiDYield",";p_{T} (GeV/c);Number of #bar{d} tracks",120,0.1,6.1);
  fDYieldTOF = new TH1F("fDYieldTOF",";p_{T} (GeV/c);Number of d tracks",120,0.1,6.1);
  fAntiDYieldTOF = new TH1F("fAntiDYieldTOF",";p_{T} (GeV/c);Number of #bar{d} tracks",120,0.1,6.1);
  fAntiDMCYield = new TH1F("fAntiDMCYield",";p_{T} (GeV/c);Number of MC #bar{d}",120,0.1,6.1);
  fDMCYield = new TH1F("fDMCYield",";p_{T} (GeV/c);Number of MC d",120,0.1,6.1);
  fEtaPhiCoverage = new TH2F("fEtaPhiCoverage",";#eta;#phi;Entries",160,-0.8,0.8,628,0,TwoPi());
  
  fOutput->Add(fDYield);
  fOutput->Add(fAntiDYield);
  fOutput->Add(fDYieldTOF);
  fOutput->Add(fAntiDYieldTOF);
  fOutput->Add(fDMCYield);
  fOutput->Add(fAntiDMCYield);
  fOutput->Add(fEtaPhiCoverage);
  for (int iCent = 0; iCent < 3; ++iCent) {
    fDCAxyPrimariesD[iCent] = new TH1F(Form("fDCAxyPrimariesD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fDCAxySecondariesD[iCent] = new TH1F(Form("fDCAxySecondariesD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fDCAzPrimariesD[iCent] = new TH1F(Form("fDCAzPrimariesD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fDCAzSecondariesD[iCent] = new TH1F(Form("fDCAzSecondariesD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fOutput->Add(fDCAxyPrimariesD[iCent]);
    fOutput->Add(fDCAxySecondariesD[iCent]);
    fOutput->Add(fDCAzPrimariesD[iCent]);
    fOutput->Add(fDCAzSecondariesD[iCent]);
    fDCAxyPrimariesAD[iCent] = new TH1F(Form("fDCAxyPrimariesAD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fDCAxySecondariesAD[iCent] = new TH1F(Form("fDCAxySecondariesAD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fDCAzPrimariesAD[iCent] = new TH1F(Form("fDCAzPrimariesAD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fDCAzSecondariesAD[iCent] = new TH1F(Form("fDCAzSecondariesAD%i",iCent),";DCA_{xy} (cm);Entries",500,-5,5);
    fOutput->Add(fDCAxyPrimariesAD[iCent]);
    fOutput->Add(fDCAxySecondariesAD[iCent]);
    fOutput->Add(fDCAzPrimariesAD[iCent]);
    fOutput->Add(fDCAzSecondariesAD[iCent]);
  }
  PostData(1,fOutput);
}

//__________________________________________________________________________________________________
void AODdEfficiency::UserExec(Option_t *){
  AliAODEvent *ev = dynamic_cast<AliAODEvent*> (InputEvent());
  
  // Check event selection mask
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t mask = handl->IsEventSelected();
  if (!(mask & 0xffffffff)) {
    PostData(1, fOutput);
    return;
  }
  if (!((mask & AliVEvent::kMB) || (mask & AliVEvent::kCentral) || (mask & AliVEvent::kSemiCentral))) {
    PostData(1, fOutput);
    return;
  }
  
  // get branch "mcparticles"
  TClonesArray *stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!stack) {
    PostData(1, fOutput);
    return;
  }
  
  // Centrality quality (do we need this with AODs?)
  AliCentrality *centr=ev->GetCentrality();
  if (centr->GetQuality() != 0) {
    PostData(1, fOutput);
    return;
  }
  
  // Centrality selection in PbPb, percentile determined with V0
  Float_t centralityPercentile = centr->GetCentralityPercentile("V0M");
  //select only events with centralities between 0 and 80 %
  if (centralityPercentile < 0. || centralityPercentile > 80.) {
    PostData(1, fOutput);
    return;
  }
  
  // Primary vertex displacement cut
  const AliVVertex *primaryVtx = ev->GetPrimaryVertex();
  if(TMath::Abs(primaryVtx->GetZ()) > 10.) {
    PostData(1, fOutput);
    return;
  }
  
  Int_t cent = -1;
  if (centralityPercentile < 5.f)
    cent = 0;
  else if (centralityPercentile > 20.f && centralityPercentile < 40.f)
    cent = 1;
  else if (centralityPercentile > 40.f && centralityPercentile < 60.f)
    cent = 2;
  
  // Used later for the track propagation
  Double_t b = (Double_t)ev->GetMagneticField();
  Double_t dca_tr[2] = {0.,0.}, covar_tr[2] = {0.,0.};
  
  // Making the list of deuterons in acceptance
  TList mcD, mcDbar;
  mcD.SetOwner(kFALSE);
  mcDbar.SetOwner(kFALSE);
  for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {
    AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
    //if (!part->IsPhysicalPrimary()) continue;
    if (TMath::Abs(part->Eta()) > 0.8) continue;
    if (TMath::Abs(part->Y()) > 0.5) continue;
    const float pt = part->Pt();
    const int pdg = part->GetPdgCode();
    if (pdg == 1000010020) {
      fDMCYield->Fill(pt);
      mcD.Add(part);
    } else if (pdg == -1000010020) {
      fAntiDMCYield->Fill(pt);
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
    if (TMath::Abs(aodtr->Eta()) > 0.8f) continue;
    if (aodtr->Chi2perNDF() > 4.f) continue;
    AliAODVertex *vtx1 = (AliAODVertex*)aodtr->GetProdVertex();
    if(Int_t(vtx1->GetType()) == AliAODVertex::kKink) continue;
    if (aodtr->GetTPCClusterMap().CountBits() < 70) continue;
    unsigned int nSPD = 0, nITS = 0;
    for (int i = 0; i < 6; ++i) {
      if (aodtr->HasPointOnITSLayer(i)) {
        if(i < 2) nSPD++;
        nITS++;
      }
    }
    if (nSPD < 1 || nITS < 2) continue;
    if (!aodtr->PropagateToDCA(primaryVtx, b, 100, dca_tr, covar_tr)) continue;
    if (dca_tr[1] > 2.f) continue;
    // Has the track a nice prolongation in TOF? Answer in hasTOF.
    Bool_t hasTOFout  = aodtr->GetStatus() & AliESDtrack::kTOFout;
    Bool_t hasTOFtime = aodtr->GetStatus() & AliESDtrack::kTIME;
    Float_t length = aodtr->GetIntegratedLength();
    Bool_t hasTOF = Bool_t(hasTOFout & hasTOFtime) && length > 350.f;
    
    // Getting the particle and checking if it is in one of the two lists of MC particles
    AliAODMCParticle *part = (AliAODMCParticle*) stack->At(aodtr->GetLabel());
    if (!part) continue;
    if (mcD.Contains(part)) {
      if (part->IsPhysicalPrimary()) {
        fEtaPhiCoverage->Fill(part->Eta(),part->Phi());
        fDYield->Fill(part->Pt());
        if (hasTOF)
          fDYieldTOF->Fill(part->Pt());
        if (cent >= 0) {
          fDCAxyPrimariesD[cent]->Fill(dca_tr[0]);
          fDCAzPrimariesD[cent]->Fill(dca_tr[1]);
        }
      } else {
        if (cent >= 0) {
          fDCAxySecondariesD[cent]->Fill(dca_tr[0]);
          fDCAzSecondariesD[cent]->Fill(dca_tr[1]);
        }
      }
    } else if (mcDbar.Contains(part)) {
      if (part->IsPhysicalPrimary()) {
        fEtaPhiCoverage->Fill(part->Eta(),part->Phi());
        fAntiDYield->Fill(part->Pt());
        if (hasTOF)
          fAntiDYieldTOF->Fill(part->Pt());
        if (cent >= 0) {
          fDCAxyPrimariesAD[cent]->Fill(dca_tr[0]);
          fDCAzPrimariesAD[cent]->Fill(dca_tr[1]);
        }
      } else {
        if (cent >= 0) {
          fDCAxySecondariesAD[cent]->Fill(dca_tr[0]);
          fDCAzSecondariesAD[cent]->Fill(dca_tr[1]);
        }
      }
    }
  } // End AOD track loop
  
  //  Post output data.
  PostData(1,fOutput);
}

//__________________________________________________________________________________________________
void AODdEfficiency::Terminate(Option_t *) {
  // Merge output
  // Called once at the end of the query
  
  return;
}
