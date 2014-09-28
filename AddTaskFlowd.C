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

  //Int_t iResult = task->Initialize();
  // if (!iResult)
  mgr->AddTask(task);
    // else {
    //AliError("NO pt ranges specfied, not adding the task !!!");
    //return -1;
    //}

  //mgr->AddTask(task);
  
  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  //AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  //dumm output container
  // AliAnalysisDataContainer *coutput0 =
  //     mgr->CreateContainer("mpuccio_treeFlowd",
  //                          TTree::Class(),
  //                          AliAnalysisManager::kExchangeContainer,
  //                          "mpuccio_default");

  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("mpuccio_Flowd", TList::Class(),AliAnalysisManager::kOutputContainer,"mpuccio_Flowd.root");

  //connect containers
  mgr->ConnectInput  (task,  0,  mgr->GetCommonInputContainer());
  //mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);
  task->SetFillTree(fillTree);
//  if (fillTree) {
//    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("mpuccio_FlowdTree", TNtuple::Class(),
//                                                              AliAnalysisManager::kOutputContainer,
//                                                              "mpuccio_Flowd_tree.root");
//    coutput2->SetSpecialOutput();
//    mgr->ConnectOutput (task,  2, coutput2);
//  }

  return task;
}
