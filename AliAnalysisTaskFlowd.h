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
class TTree;

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
  void ResetTreeVariables();
  
  // Private variables
  AliESDEvent          *fESD;                               //! ESD object
  AliESDpid            *fESDpid;                            //! basic TPC object for n-sigma cuts
  AliESDtrackCuts       fESDtrackCuts;                      //  basic cut variables
  AliESDtrackCuts       fESDtrackCutsStrict;                //  basic cut variables
  AliInputEventHandler *fEventHandler;                      //! for ESDs or AODs
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
  TList                *fOutputContainer;                   //! Output data container
  Int_t                 fTrigger;                           //  "Trigger mask"
  TTree                *fTree;                              //! Deuteron tree
  
  // tree variables
  Char_t     fName[1000];
  Int_t      fEvnt;
  Char_t     fFileName[1000];
  Int_t      fEventNumber;
  Float_t    fCentrality;
  //
  Int_t      fItrk;
  //
  Double_t   fEta[1000];
  Int_t      fKinkIndex[1000];
  //
  UShort_t   fTPCNsignal[1000];
  UShort_t   fTPCnCluster[1000];
  Double_t   fChi2PerClusterTPC[1000];
  Bool_t     fTPCRefit[1000];
  Int_t      fTPCSharedClusters[1000];
  UShort_t   fTPCNclsIter1[1000];
  //
  Double_t   fITSsignal[1000];
  Int_t      fITSnCluster[1000];
  Int_t      fITSnClusterPID[1000];
  Double_t   fChi2PerClusterITS[1000];
  Bool_t     fITSRefit[1000];
  //
  Bool_t     fTOFRefit[1000];
  Bool_t     fTOFtime[1000];
  Bool_t     fTOFout[1000];
  Double_t   fTOFsignalDz[1000];
  Double_t   fTOFsignalDx[1000];
  //
  Float_t    fDCAZ[1000];
  Float_t    fDCAXY[1000];
  //
  Double_t   fTrkPtot[1000];
  Double_t   fTPCPtot[1000];
  Double_t   fTrackPt[1000];
  Double_t   fDeDx[1000];
  Double_t   fSign[1000];
  Float_t    fMass[1000];
  Float_t    fTime[1000];
  Float_t    fLength[1000];
  Double_t   fSigmaQP[1000];
  
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
