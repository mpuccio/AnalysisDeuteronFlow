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
class AliPIDResponse;
class TTree;
class AliESDtrackCuts;
class AliInputEventHandler;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowd : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskFlowd(const char *name = "DeuteronAnalysis");
  virtual ~AliAnalysisTaskFlowd();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  Int_t          SetupEvent();
  void           ResetEvent();
  void           SetFillTree(Bool_t io = kFALSE) { fFillTree = io; }
  
private:
  // Private methods
  AliAnalysisTaskFlowd(const AliAnalysisTaskFlowd&);            //! Not implemented
  AliAnalysisTaskFlowd& operator=(const AliAnalysisTaskFlowd&); //! Not implemented
  void BinLogAxis(const TH1 *h);
  Bool_t IsTriggered();
  
  // Private variables
  AliESDEvent          *fEvent;                             //! ESD object
  AliESDtrackCuts      *fESDtrackCuts;                      //  Cuts
  AliPIDResponse       *fPid;                               //! basic TPC object for n-sigma cuts
  AliInputEventHandler *fEventHandler;                      //! for ESDs or AODs
  Bool_t                fFillTree;                          //  TTree switch
  TH1F                 *fHistCentralityClass10;             //! centrality distribution
  TH1F                 *fHistCentralityPercentile;          //! centrality distribution
  TH2F                 *fHistDeDx;                          //! Histo for a dE/dx
  TH2F                 *fHistDeDxITS;                       //! Histo for a dE/dx ITS
  TH1F                 *fHistDeuteron;                      //! d plot TOF mass
  TH2F                 *fHistTOF2D;                         //! Histo for a TOF
  TH2F                 *fHistTOFnuclei;                     //! Histo for a TOF nuclei
  TH1F                 *fHistTriggerStat;                   //! Trigger statistic
  TH1F                 *fHistTriggerStatAfterEventSelection;//! Trigger statistic after selection
  Int_t                 fNCounter;                          //  # points in the signal graph
  TList                *fOutputContainer;                   //! Output data container
  Int_t                 fTrigger;                           //  "Trigger mask"
  
  TTree                *fTree;
  Float_t               fTCentrality;
  Float_t               fTEta;
  Float_t               fTPhi;
  Float_t               fTTPCchi2;
  Float_t               fTTPCsignal;
  Float_t               fTITSsignal;
  Float_t               fTBeta;
  Float_t               fTDCAxy;
  Float_t               fTDCAz;
  Float_t               fTMomentum;
  Float_t               fTPtpc;
  Float_t               fTSignedPt;
  Float_t               fTSigmaQoverPt;
  UShort_t              fTTPCnClust;
  UShort_t              fTTPCnSignal;
  Char_t                fTITSnClusters;
  Char_t                fTITSnSignal;
  Char_t                fTTrigger;
  Bool_t                fTEventCounter;
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
#endif /* defined(ALIANALYSISTASKFLOWD_H) */
