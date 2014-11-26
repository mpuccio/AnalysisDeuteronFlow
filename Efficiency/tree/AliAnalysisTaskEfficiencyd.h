//
//  AliAnalysisTaskEfficiencyd.h
//
//  Created by Maximiliano Puccio on 28/10/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#ifndef __AliAnalysisTaskEfficiencyd__
#define __AliAnalysisTaskEfficiencyd__

class AliAODEvent;
class AliAODTrack;
class AliAODpidUtil;
class TTree;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEfficiencyd : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskEfficiencyd();
  virtual ~AliAnalysisTaskEfficiencyd();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
private:
  AliAnalysisTaskEfficiencyd(const AliAnalysisTaskEfficiencyd &source);
  AliAnalysisTaskEfficiencyd &operator=(const AliAnalysisTaskEfficiencyd &source);
  
  AliAODEvent     *fAOD;              // AOD object
  AliAODpidUtil   *fAODpid;           //! basic TPC object for n-sigma cuts
  
  TTree   *fTree;
  Float_t  fTpMC;
  Float_t  fTpTMC;
  Float_t  fTetaMC;
  Float_t  fTyMC;
  Float_t  fTphiMC;
  Float_t  fTp;
  Float_t  fTpT;
  Float_t  fTeta;
  Float_t  fTphi;
  Float_t  fTbeta;
  Float_t  fTDCAxy;
  Float_t  fTDCAz;
  Float_t  fTchi2;
  Float_t  fTcentrality;
  UShort_t fTTPCnClusters;
  UShort_t fTTPCnSignal;
  Char_t   fTITSnClusters;
  Char_t   fTITSnSignal;
  Bool_t   fTIsPrimary;
  
  ClassDef(AliAnalysisTaskEfficiencyd, 1);
};


#endif /* defined(__AliAnalysisTaskEfficiencyd__) */
