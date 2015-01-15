//
//  AliAnalysisTaskEfficiencydAOD.h
//
//  Created by Maximiliano Puccio on 28/10/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#ifndef __AliAnalysisTaskEfficiencydAOD__
#define __AliAnalysisTaskEfficiencydAOD__

class TTree;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEfficiencydAOD : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskEfficiencydAOD();
  virtual ~AliAnalysisTaskEfficiencydAOD();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
private:
  AliAnalysisTaskEfficiencydAOD(const AliAnalysisTaskEfficiencydAOD &source);
  AliAnalysisTaskEfficiencydAOD &operator=(const AliAnalysisTaskEfficiencydAOD &source);
  
  AliPIDResponse *fPIDResponse;
  TTree          *fTree;
  Float_t         fTpMC;
  Float_t         fTpTMC;
  Float_t         fTetaMC;
  Float_t         fTyMC;
  Float_t         fTphiMC;
  Float_t         fTp;
  Float_t         fTpT;
  Float_t         fTeta;
  Float_t         fTphi;
  Float_t         fTbeta;
  Float_t         fTDCAxy;
  Float_t         fTDCAz;
  Float_t         fTchi2;
  Float_t         fTcentrality;
  UShort_t        fTTPCnClusters;
  UShort_t        fTTPCnSharedClusters;
  UShort_t        fTTPCnSignal;
  Char_t          fTITSnClusters;
  Char_t          fTITSnSignal;
  Char_t          fTParticleSpecie;
  Bool_t          fTIsPrimary;
  Bool_t          fTIsSecondaryFromMaterial;
  
  ClassDef(AliAnalysisTaskEfficiencydAOD, 1);
};


#endif /* defined(__AliAnalysisTaskEfficiencydAOD__) */
