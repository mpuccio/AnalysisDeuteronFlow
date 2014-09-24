//
//  AliAnalysisTaskFlowd.h
//
//  Created by Maximiliano Puccio on 14/09/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#ifndef ALIANALYSISTASKFLOWD_H
#define ALIANALYSISTASKFLOWD_H

class TF1;
class TH1F;
class TH2F;
class TH3F;
class AliESDEvent;
class AliESDVertex;
class AliESDpid;

#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "TH3F.h"
#include "TGraph.h"
#include "AliStack.h"

class AliAnalysisTaskFlowd : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskFlowd(const char *name = "DeuteronAnalysis");
  virtual ~AliAnalysisTaskFlowd() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  Int_t  Initialize();
  Int_t  SetupEvent();
  void   ResetEvent();
  
  
private:
  // Private methods
  AliAnalysisTaskFlowd(const AliAnalysisTaskFlowd&);            //! Not implemented
  AliAnalysisTaskFlowd& operator=(const AliAnalysisTaskFlowd&); //! Not implemented
  void BinLogAxis(const TH1 *h);
  void BinLogAxis(const TH3 *h, Int_t axisNumber);
  Bool_t IsTriggered();
  void SetBBParameters(Int_t runNumber);
  
  // Private variables
  Double_t              fBBParametersLightParticles[5];     //! Bethe Bloch params light particles
  AliESDEvent          *fESD;                               //! ESD object
  AliESDpid            *fESDpid;                            //! basic TPC object for n-sigma cuts
  AliESDtrackCuts       fESDtrackCuts;                      //  basic cut variables
  AliESDtrackCuts       fESDtrackCutsSharp;                 //  sharp cut variables -> final results
  AliInputEventHandler *fEventHandler;                      //! for ESDs or AODs
  TH1F                 *fHistCentralityClass10;             //! centrality distribution
  TH1F                 *fHistCentralityPercentile;          //! centrality distribution
  TH2F                 *fHistDeDx;                          //! Histo for a dE/dx
  TH3F                 *fHistDeDxRegion;                    //! Histo for a dE/dx per Region
  TH2F                 *fHistDeDxSharp;                     //! Histo for a dE/dx with sharp cuts
  TH1F                 *fHistDeuteron;                      //! d plot TOF mass
  TH1F                 *fHistDeuteronSignal;                //! d TOF mass for signal candidates
  TH2F                 *fHistTOF2D;                         //! Histo for a TOF
  TH2F                 *fHistTOFnuclei;                     //! Histo for a TOF nuclei
  TH1F                 *fHistTriggerStat;                   //! Trigger statistic
  TH1F                 *fHistTriggerStatAfterEventSelection;//! Trigger statistic after selection
  Int_t                 fNCounter;                          //  # points in the signal graph
  TList                *fOutputContainer;                   //! Output data container
  Int_t                 fTrigger;                           //  "Trigger mask"
  
  // Constants
  const Int_t           fkNTriggers;                        //  Number of used triggers
  enum Trigger_m {
    kMB = 1,
    kCentral = 2,
    kSemiCentral = 4,
    kEMCEJE = 8,
    kEMCEGA = 16
  };
ClassDef(AliAnalysisTaskFlowd, 1)
};
#endif /* defined(__AliAnalysisTaskFlowd__) */
