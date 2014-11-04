//
//  AODdEfficiency.h
//
//  Created by Maximiliano Puccio on 28/10/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#ifndef __AODdEfficiency__
#define __AODdEfficiency__

class AliAODEvent;
class AliAODTrack;
class TH1F;
class TList;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AODdEfficiency: public AliAnalysisTaskSE {
public:
  
  AODdEfficiency();
  virtual ~AODdEfficiency();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
  // }
  
  
private:
  AODdEfficiency (const AODdEfficiency &source);
  AODdEfficiency &operator=(const AODdEfficiency &source);
  
  AliAODEvent     *fAOD;              // AOD object
  
  TList   *fOutput;                   //! tlist with output
  
  TH1F    *fAntiDMCYield;
  TH1F    *fAntiDYield;
  TH1F    *fAntiDYieldTOF;
  TH1F    *fDCAxyPrimariesAD[3];
  TH1F    *fDCAxySecondariesAD[3];
  TH1F    *fDCAzPrimariesAD[3];
  TH1F    *fDCAzSecondariesAD[3];
  TH1F    *fDCAxyPrimariesD[3];
  TH1F    *fDCAxySecondariesD[3];
  TH1F    *fDCAzPrimariesD[3];
  TH1F    *fDCAzSecondariesD[3];
  TH1F    *fDMCYield;
  TH1F    *fDYield;
  TH1F    *fDYieldTOF;
  TH2F    *fEtaPhiCoverage;
  
  ClassDef(AODdEfficiency, 1);
};


#endif /* defined(__AODdEfficiency__) */
