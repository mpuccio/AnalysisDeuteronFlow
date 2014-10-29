class AliAnalysisGrid;
const char *dataset   = "";
TString anaLibs = "";
//Int_t iESDfilter       = 1;
Int_t iESDfilter       = 0;
Int_t iAODTagCreation  = 1;
Int_t iAODAddMCBranch  = 0;

//______________________________________________________________________________
void RunGrid(
         const char* runtype = "grid", // local, proof or grid
         const char *gridmode = "test", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
         const Long64_t nentries = 400, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
         const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
         const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
         const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
         const char *taskname = "Flowd_PbPb2011"
         )
{
  // check run type
  if(runtype != "local" && runtype != "proof" && runtype != "grid"){
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tint runtype = local, proof or grid\n");
    return;
  }
  Printf("%s analysis chosen",runtype);
  
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
    plugin = CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset);
    if(!plugin) return;
  } else {
    gROOT->LoadMacro("$ALICE_ROOT/PWGCF/Correlations/macros/dphicorrelations/CreateESDChain.C");
    chain = CreateESDChain("ESDs.txt");
  }
  
  //---- Create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
  if(plugin) mgr->SetGridHandler(plugin);
  
  //  Input
  
  AliESDInputHandler* iH = new AliESDInputHandler("handler","handler for my analisys");
  mgr->SetInputEventHandler(iH);
  
  //--------------------------------------------------------------
  // Other tasks

  // Physics selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSel = AddTaskPhysicsSelection(kFALSE); // useMC
  // Centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentr = AddTaskCentrality();
  // PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(kFALSE); // useMC
  // PID QA
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQATask = AddTaskPIDqa();

  gROOT->LoadMacro("./AliAnalysisTaskFlowd.cxx+g");//$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/AliAnalysisTaskCheckCascadePbPb.cxx++g");
  gROOT->LoadMacro("./AddTaskFlowd.C");//$ALICE_ROOT/PWGLF/STRANGENESS/Cascades/macros/AddTaskCheckCascadePbPb.C");
  AliAnalysisTaskFlowd *task = AddTaskFlowd(kFALSE);
  

  
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
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode,
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
  
  plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/"); //sim
  plugin->SetDataPattern("*/pass2/*AliESDs.root"); // sim
  plugin->SetRunPrefix("000");   // real data
  
  Int_t runlist[] = {                                                               // Counter
    /*170309, 170308, 170306, 170270, 170269, 170268, 170230,*/ 170228, 170204, 170203, // 10
    170193, 170163, 170159, 170155, 170081, 170027, 169859, 169858, 169855, 169846, // 20
    169838, 169837, 169835, 169417, 169415, 169411, 169238, 169167, 169160, 169156, // 30
    169148, 169145, 169144, 169138, 169094, 169091, 169035, 168992, 168988, 168826, // 40
    168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, // 50
    168342, 168341, 168325, 168322, 168311, 168310, 167988, 167987                  // 58
  };
  
  for(Int_t i = 0; i < 1; i++)
    plugin->AddRunNumber(runlist[i]);

  plugin->SetNrunsPerMaster(1);
  plugin->SetOutputToRunNo();
  
  
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(taskname);
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/taskname/out
  
//  plugin->SetAdditionalLibs("libTree.so libGeom.so libPhysics.so libVMC.so libMinuit.so libSTEERBase.so libESD.so libAOD.so  libANALYSIS.so libOADB.so libANALYSISalice.so libCORRFW.so libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so");
  
  // plugin->SetAdditionalLibs("libCORRFW.so libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so");
  
  plugin->SetAnalysisSource("AliAnalysisTaskFlowd.cxx");
  //plugin->SetAdditionalLibs("AliAnalysisTaskFlowd.h AliAnalysisTaskFlowd.cxx ");
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
  plugin->SetAdditionalLibs("libSTEERBase.so libESD.so AliAnalysisTaskFlowd.h AliAnalysisTaskFlowd.cxx libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so");
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
  plugin->SetSplitMaxInputFileNumber(10);
  
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable(Form("%s.sh",taskname));
  
  // set number of test files to use in "test" mode
  plugin->SetNtestFiles(10);
  
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
