AODdEfficiency *AddTaskAODdEfficiency(Bool_t readMC = kTRUE,
    TString outputFileName = "dEfficiency.root") {

  // Creates, configures and attaches to the train a HyperTriton MC search task
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskITSPIdMC", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskITSPIdMC", "This task requires an input event handler");
    return NULL;
  }

  // TString type = "AOD";
  //mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if(type.Contains("ESD")) {
    ::Error("AddTaskITSPIdMC", "This task requires to run on AOD");
    return NULL;
  }

  // Create and configure the task
  AODdEfficiency *taskITSPIdMC = new AODdEfficiency();

  //taskITSPIdMC->SetReadMC(readMC);
  //taskITSPIdMC->SetCollidingSystems(lCollidingSystems);
  //taskITSPIdMC->SetAnalysisType(type);

  // TODO: inserire le funzioni chiamate e usate nella classe AODtaskLambda.cxx
  mgr->AddTask(taskITSPIdMC);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via manager as below
  if (outputFileName.IsNull()) {
    outputFileName = AliAnalysisManager::GetCommonFileName();
  }

  // if (lCollidingSystems) outputFileName += "_AA";
  //else outputFileName += "_PP";
  // if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("dEfficiency",
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    outputFileName);

  mgr->ConnectInput(taskITSPIdMC, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskITSPIdMC, 1, coutput1);

  return taskITSPIdMC;
}
