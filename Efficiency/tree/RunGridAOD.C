class AliAnalysisGrid;
const char *dataset   = "";
TString anaLibs = "";
//Int_t iESDfilter       = 1;
Int_t iESDfilter       = 0;
Int_t iAODTagCreation  = 1;
Int_t iAODAddMCBranch  = 0;

void GetRunList(TString who = "all", Int_t &start, Int_t &nRuns) {
  if (who.Contains("maximiliano",TString::kIgnoreCase)) {
    start = 0;
    nRuns = 20;
  }
  else if (who.Contains("simone",TString::kIgnoreCase)) {
    start = 20;
    nRuns = 20;
  }
  else if (who.Contains("stefania",TString::kIgnoreCase)) {
    start = 40;
    nRuns = 20;
  }
  else if (who.Contains("stefano",TString::kIgnoreCase)) {
    start = 60;
    nRuns = 20;
  }
  else if (who.Contains("massimo",TString::kIgnoreCase)) {
    start = 80;
    nRuns = 20;
  }
  else if (who.Contains("elena",TString::kIgnoreCase)) {
    start = 100;
    nRuns = 20;
  }
  else {
    start = 0;
    nRuns = 45;
  }
  return;
}

//______________________________________________________________________________
void RunGridAOD(TString runtype = "grid", // local, proof or grid
                TString gridmode = "test", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
                const Long64_t nentries = 400, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
                const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
                const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
                const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
                const char *taskname = "Efficiencyd"
                )
{
  // check run type
  TString who = runtype;
  if (who.Contains("stefano",TString::kIgnoreCase) || who.Contains("simone",TString::kIgnoreCase) || who.Contains("maximiliano",TString::kIgnoreCase) || who.Contains("stefania",TString::kIgnoreCase) || who.Contains("elena",TString::kIgnoreCase) || who.Contains("massimo",TString::kIgnoreCase) || who.Contains("all",TString::kIgnoreCase)) {
    runtype = "grid";
    gridmode = "full";
  }
  
  if(runtype != "local" && runtype != "proof" && runtype != "grid"){
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tint runtype = local, proof or grid\n");
    return;
  }
  Printf("%s analysis chosen",runtype.Data());
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGHF/vertexingHF");
  // Load analysis specific libraries
  //=====================================================================
  
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -g");
  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGHFbase.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");
  gSystem->Load("libPWGHFvertexingHF.so");
  
  // Load analysis specific libraries
  //=====================================================================
  //------ Create AlienPlugin ---------------------
  AliAnalysisGrid *plugin = 0x0;
  TChain *chain = 0x0;
  if (runtype != "local") {
    plugin = CreateAlienHandler(who,taskname, gridmode.Data(), proofcluster, proofdataset);
    if(!plugin) return;
  } else {
    gROOT->LoadMacro("$ALICE_ROOT/PWGCF/Correlations/macros/dphicorrelations/CreateESDChain.C");
    chain = CreateESDChain("ESDs.txt");
  }
  
  //---- Create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
  if(plugin) mgr->SetGridHandler(plugin);
  
  //  Input
  
  AliAODInputHandler* iH = new AliAODInputHandler("handler","handler for my analisys");
  mgr->SetInputEventHandler(iH);
  
  //--------------------------------------------------------------
  // Other tasks

  // PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(kFALSE); // useMC
  // PID QA
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQATask = AddTaskPIDqa();
  
  gROOT->LoadMacro("./AliAnalysisTaskEfficiencydAOD.cxx+g");//$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/AliAnalysisTaskCheckCascadePbPb.cxx++g");
  gROOT->LoadMacro("./AddTaskEfficiencyd.C");//$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/macros/AddTaskCheckCascadePbPb.C");
  AliAnalysisTaskEfficiencydAOD *task = AddTaskEfficiencyd(kTRUE);//kTRUE);
  

  
  //__________________________________________________________________________
  // Disable debug printouts
  mgr->SetDebugLevel(3);
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  AliLog::SetGlobalDebugLevel(0);
  
  //__________________________________________________________________________
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  // Start analysis in grid.
  if (runtype == "local")
    mgr->StartAnalysis(runtype,chain,nentries,firstentry);
  else
    mgr->StartAnalysis(runtype,nentries,firstentry);
}


//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(TString &who,const char *taskname, const char *gridmode,
                                    const char *proofcluster, const char *proofdataset)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(gridmode);
  
  // Set versions of used packages
  
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-08");
  plugin->SetAliROOTVersion("vAN-20140915");
  plugin->SetExecutableCommand("aliroot -b -q");
  
  // Declare input data to be processed.
  
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  //    plugin->SetGridDataDir("/alice/data/2010/LHC10h");
  
  // plugin->SetGridDataDir(" /alice/data/2011/LHC11h_2/"); //sim
  // plugin->SetDataPattern("pass2/*AliAOD.root"); // sim
  
  plugin->SetGridDataDir("/alice/sim/2013/LHC13d15/"); //sim
//  plugin->SetGridDataDir("/alice/sim/2014/LHC14a6/"); //sim
  plugin->SetDataPattern("*/AliAOD.root"); // sim
  plugin->SetRunPrefix("");   // mc
  
  Int_t runlist[89] = {
//    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328,
//    139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029,
//    139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653,
//    138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469,
//    138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192,
//    138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
//    137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546,
//    137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434,
//    137432, 137431, 137366, 137243, 137236, 137235, 137232, 137231, 137162, 137161 // Counter
    137161, 137162, 137231, 137232, 137235, 137236, 137243, 137366, 137431, 137432,
    137434, 137439, 137440, 137441, 137443, 137530, 137531, 137539, 137541, 137544,
    137546, 137549, 137595, 137608, 137638, 137639, 137685, 137686, 137691, 137692,
    137693, 137704, 137718, 137722, 137724, 137751, 137752, 137844, 137848, 138190,
    138192, 138197, 138201, 138225, 138275, 138364, 138396, 138438, 138439, 138442,
    138469, 138534, 138578, 138579, 138582, 138583, 138621, 138638, 138652, 138653,
    138662, 138666, 138730, 138732, 138837, 138870, 138871, 138872, 139028, 139029,
    139036, 139037, 139038, 139105, 139107, 139173, 139309, 139310, 139314, 139328,
    139329, 139360, 139437, 139438, 139465, 139503, 139505, 139507, 139510
  };
  
  Int_t start = 0, nrun = 0;
  GetRunList(who,start,nrun);
  
  for(Int_t i = start; i < start + nrun; i++) {
    plugin->AddRunNumber(runlist[i]);
  }

  plugin->SetNrunsPerMaster(15);
  plugin->SetOutputToRunNo();
  
  
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("Efficiencyd2013");
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/taskname/out
  
//  plugin->SetAdditionalLibs("libTree.so libGeom.so libPhysics.so libVMC.so libMinuit.so libSTEERBase.so libESD.so libAOD.so  libANALYSIS.so libOADB.so libANALYSISalice.so libCORRFW.so libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so");
  
  // plugin->SetAdditionalLibs("libCORRFW.so libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so");
  
  plugin->SetAnalysisSource("AliAnalysisTaskEfficiencydAOD.cxx");
  //plugin->SetAdditionalLibs("AliAnalysisTaskEfficiencyd.h AliAnalysisTaskEfficiencyd.cxx ");
  cout<<"-->>>>>>>>>>>>>>>>>>>>>>>>> 1"<<endl;
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  //plugin->SetAdditionalLibs("libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so AliAODMuonReplicator0.so");
  
  cout<<"-->>>>>>>>>>>>>>>>>>>>>>>>> 2"<<endl;
  // plugin->SetAdditionalLibs("libTree.so libGeom.so libPhysics.so libVMC.so libMinuit.so libSTEERBase.so libESD.so libAOD.so libANALYSIS.so libOADB.so libANALYSISalice.so libCORRFW.so libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so");
  
  // plugin->SetAnalysisSource("AliAnalysisTaskESDMuonFilterO.cxx");
  //plugin->SetAnalysisSource("AliAODMuonReplicator0.cxx");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  // plugin->SetAdditionalLibs("AliAODMuonReplicator0_cxx.so");
  // plugin->SetAdditionalLibs("AliAnalysisTaskESDMuonFilterO_cxx.so");
  
  //questo
  plugin->SetAdditionalLibs("libSTEERBase.so libESD.so AliAnalysisTaskEfficiencydAOD.h AliAnalysisTaskEfficiencydAOD.cxx libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -g");
  
  
  // plugin->SetAdditionalLibs("AliAODMuonReplicator0.h AliAODMuonReplicator0.cxx");
  // plugin->SetAdditionalLibs("AliAnalysisTaskESDMuonFilterO.h AliAnalysisTaskESDMuonFilterO.cxx");
  
  //plugin->SetAdditionalLibs("AliAODMuonReplicator0.h AliAODMuonReplicator0.cxx AliAnalysisTaskESDMuonFilterO.h AliAnalysisTaskESDMuonFilterO.cxx");
  
  cout<<"-->>>>>>>>>>>>>>>>>>>>>>>>> 3"<<endl;
	
  // plugin->SetAdditionalLibs("AliAODMuonReplicator0.h AliAODMuonReplicator0.cxx");
  
  // plugin->SetAdditionalLibs("AliAnalysisTaskESDMuonFilterO.h AliAnalysisTaskESDMuonFilterO.cxx");
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs();
  //plugin->SetOutputFiles("list.root");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s.C",taskname));
  
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(35);
  
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable(Form("%s.sh",taskname));
  
  // set number of test files to use in "test" mode
  plugin->SetNtestFiles(2);
  
  // file containing a list of chuncks to be used for testin
  plugin->SetFileForTestMode("testdata");
  
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName(Form("%s.jdl",taskname));
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  //----------------------------------------------------------
  //---      PROOF MODE SPECIFIC SETTINGS         ------------
  //----------------------------------------------------------
  // Proof cluster
  plugin->SetProofCluster(proofcluster);
  // Dataset to be used
  plugin->SetProofDataSet(proofdataset);
  // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
  plugin->SetProofReset(0);
  // May limit number of workers
  plugin->SetNproofWorkers(0);
  // May limit the number of workers per slave
  plugin->SetNproofWorkersPerSlave(1);
  // May use a specific version of root installed in proof
  plugin->SetRootVersionForProof("current");
  // May set the aliroot mode. Check http://aaf.cern.ch/node/83
  plugin->SetAliRootMode("default"); // Loads AF libs by default
  // May request ClearPackages (individual ClearPackage not supported)
  plugin->SetClearPackages(kFALSE);
  // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
  plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root etc
  // Request connection to alien upon connection to grid
  plugin->SetProofConnectGrid(kFALSE);
  
  cout<<"-->>>>>>>>>>>>>>>>>>>>>>>>> 4"<<endl;
	
  
  return plugin;
  
  cout<<"-->>>>>>>>>>>>>>>>>>>>>>>>> 5"<<endl;
	
}
