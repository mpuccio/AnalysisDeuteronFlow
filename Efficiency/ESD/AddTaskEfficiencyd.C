AliAnalysisTask *AddTaskEfficiencyd(){
  
  
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEfficiencyd", "No analysis manager found.");
    return 0;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEfficiencyd", "This task requires an input event handler");
    return NULL;
  }
  
  if (!mgr->GetMCtruthEventHandler()) {
    Error("AddAnalysisTaskParticleEfficiency", "cannot get MC truth event handler");
    return NULL;
  }
  
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskEfficiencyd* task = new AliAnalysisTaskEfficiencyd();
  
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *coutput1 =
  mgr->CreateContainer("mpuccio_Efficiency", TTree::Class(),
                       AliAnalysisManager::kOutputContainer,"mpuccio_Efficiency.root");
  coutput1->SetSpecialOutput();
  
  //connect containers
  mgr->ConnectInput  (task,  0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1,  coutput1);
  
  return task;
}
