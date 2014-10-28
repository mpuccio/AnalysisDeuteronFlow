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
  
  AliAODEvent     *fAOD; //AOD object
  
  TList   *fOutput;                   //! tlist with output
  
  TH1F    *fDYield;
  TH1F    *fAntiDYield;
  TH1F    *fDYieldTOF;
  TH1F    *fAntiDYieldTOF;
  TH1F    *fAntiDMCYield;
  TH1F    *fDMCYield;
  TH1F    *fDClones;
  
  ClassDef(AODdEfficiency, 1);
};


#endif /* defined(__AODdEfficiency__) */
