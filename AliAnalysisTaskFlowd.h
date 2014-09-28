//
//  AliAnalysisTaskFlowd.h
//
//  Created by Maximiliano Puccio on 14/09/14.
//  Copyright (c) 2014 ALICE. All rights reserved.
//

#ifndef ALIANALYSISTASKFLOWD_H
#define ALIANALYSISTASKFLOWD_H

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDVertex;
class AliESDpid;
class TNtuple;

#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowd : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskFlowd(const char *name = "DeuteronAnalysis");
  virtual ~AliAnalysisTaskFlowd() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  Int_t          Initialize();
  Int_t          SetupEvent();
  void           ResetEvent();
  
private:
  // Private methods
  AliAnalysisTaskFlowd(const AliAnalysisTaskFlowd&);            //! Not implemented
  AliAnalysisTaskFlowd& operator=(const AliAnalysisTaskFlowd&); //! Not implemented
  void BinLogAxis(const TH1 *h);
  Bool_t IsTriggered();
  
  // Private variables
  AliESDEvent          *fESD;                               //! ESD object
  AliESDpid            *fESDpid;                            //! basic TPC object for n-sigma cuts
  AliESDtrackCuts       fESDtrackCuts;                      //  basic cut variables
  AliESDtrackCuts       fESDtrackCutsStrict;                //  basic cut variables
  AliInputEventHandler *fEventHandler;                      //! for ESDs or AODs
  Bool_t                fFillTree;                          //  TTree switch
  TH1F                 *fHistCentralityClass10;             //! centrality distribution
  TH1F                 *fHistCentralityPercentile;          //! centrality distribution
  TH2F                 *fHistDeDx;                          //! Histo for a dE/dx
  TH2F                 *fHistDeDxITS;                       //! Histo for a dE/dx ITS
  TH2F                 *fHistDeDxITSsa;                     //! Histo for a dE/dx ITSsa
  TH1F                 *fHistDeuteron;                      //! d plot TOF mass
  TH2F                 *fHistTOF2D;                         //! Histo for a TOF
  TH2F                 *fHistTOFnuclei;                     //! Histo for a TOF nuclei
  TH1F                 *fHistTriggerStat;                   //! Trigger statistic
  TH1F                 *fHistTriggerStatAfterEventSelection;//! Trigger statistic after selection
  Int_t                 fNCounter;                          //  # points in the signal graph
  TNtuple              *fNtuple;                            //! Deuteron ntuple
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
