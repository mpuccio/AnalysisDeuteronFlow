// Note: it is sufficient to run this macro from ROOT (AliRoot is not needed!)
//   root -l runProof.C
{

  // Not needed on TAF
  //gEnv->SetValue("XSec.GSI.DelegProxy","2");

  TList *list = new TList(); 
  list->Add( new TNamed("ALIROOT_MODE", "AliRoot") );
  list->Add( new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:OADB:ANALYSISalice:CORRFW:OADB:PWGmuon") );
  list->Add( new TNamed("ALIROOT_EXTRA_INCLUDES", "PWGHF/vertexingHF") );

  // Different on TAF
  TProof::Open("pod://");

  // New on TAF: a single AliRoot package (correct version automatically OK)
  // Uncomment the following line to upgrade the PARfile
  // TFile::Cp("http://personalpages.to.infn.it/~berzano/cloud/extras/AliRoot.par", "AliRoot.par");
  gProof->UploadPackage("AliRoot.par");
  gProof->EnablePackage("AliRoot.par", list);  // this "list" is the same as always

  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train");

  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  if (gProof->Load("AODdEfficiency.cxx+g") != 0) return; // task
  if (gROOT->LoadMacro("AddTaskAODdEfficiency.C") != 0) return;   // addtask
  AODdEfficiency *taskITSPIdMC = AddTaskAODdEfficiency();

  // Calibration
  if (gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C")) return;
  // AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kFALSE, kTRUE); // dati reali
  AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kTRUE, kTRUE); // MC
 
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  //  mgr->StartAnalysis("proof", "Sim;Period=LHC11b9_1;Run=137161;Variant=AOD048|Sim;Period=LHC11b9_1;Run=137162;Variant=AOD048");
  mgr->StartAnalysis("proof", "Find;BasePath=/alice/sim/2014/LHC14a6;FileName=AOD/*/AliAOD.root");
}
