AliAnalysisTask *AddTaskFlowd(Bool_t fillTree = kFALSE){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowd", "No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskFlowd", "This task requires an input event handler");
      return NULL;
   }  


  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskFlowd *task = new AliAnalysisTaskFlowd("Flowd");

  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("mpuccio_Flowd", TList::Class(),
                           AliAnalysisManager::kOutputContainer,"mpuccio_Flowd.root");

  //connect containers
  mgr->ConnectInput  (task,  0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput1);
  task->SetFillTree(fillTree);
  if (fillTree) {
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("mpuccio_FlowdTree", TTree::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              "mpuccio_Flowdnt.root");
    coutput2->SetSpecialOutput();
    mgr->ConnectOutput (task,  2, coutput2);
  }

  return task;
}
