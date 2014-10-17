void RunPoD(
  TString dataset = "Data;Period=LHC11h_2;Run=170593;Variant=ESD;Pass=pass2",
  Bool_t usePhysicsSelection = kTRUE,
  Int_t numEvents = 999999999,
  Int_t firstEvent = 0
) {

  // Not needed on the VAF
  //gEnv->SetValue("XSec.GSI.DelegProxy","2");
   
  TString alirootMode = "ALIROOT";
  TString extraLibs = "ANALYSIS:ANALYSISalice"; // extraLibs = "ANALYSIS:OADB:ANALYSISalice:CORRFW:OADB:PWGmuon";

  TList *list = new TList(); 
  list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));  // important: creates token on every PROOF worker

  // Not needed on the VAF
  //TProof::Mgr("alice-caf.cern.ch")->SetROOTVersion("VO_ALICE@ROOT::v5-34-08");

  // Note the difference between CAF and VAF
  //TProof::Open("alice-caf.cern.ch");
  TProof::Open("pod://");

  // Check the dataset before running the analysis!
  gProof->ShowDataSet( dataset.Data() );
  //return;
  gProof->SetParameter("PROOF_UseMergers", (Int_t)0);
  // Not needed on the VAF
  //gProof->EnablePackage("VO_ALICE@AliRoot::v5-04-81-AN", list);

  // A single AliRoot package for *all* AliRoot versions: new on VAF
  TFile::Cp("http://personalpages.to.infn.it/~berzano/cloud/extras/AliRoot.par", "AliRoot.par");
  gProof->UploadPackage("AliRoot.par");
  gProof->EnablePackage("AliRoot.par", list);  // this "list" is the same as always

  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train");

  AliESDInputHandler *esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSel = AddTaskPhysicsSelection(kFALSE);
  // Centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentr = AddTaskCentrality();
  // PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(kFALSE); // useMC
  // PID QA
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQATask = AddTaskPIDqa();
  
  
  gProof->Load("AliAnalysisTaskFlowd.cxx+");  // DON'T use double '+' when running multiple times: it uselessly recompiles everything!
  gROOT->LoadMacro("AddTaskFlowd.C");
  AliAnalysisTaskFlowd *flowTask = AddTaskFlowd();
  //if (usePhysicsSelection) {
  //flowTask->SelectCollisionCandidates(AliVEvent::kAny);
    //}
 
  if (!mgr->InitAnalysis()) return;

  mgr->StartAnalysis("proof", dataset, numEvents, firstEvent);

}
