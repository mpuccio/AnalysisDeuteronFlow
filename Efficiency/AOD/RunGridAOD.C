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
    nRuns = 36;
  }
  else if (who.Contains("stefania",TString::kIgnoreCase)) {
    start = 36;
    nRuns = 36;
  }
  else if (who.Contains("massimo",TString::kIgnoreCase)) {
    start = 72;
    nRuns = 36;
  }
  else {
    start = 0;
    nRuns = 108;
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
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  
  // Load analysis specific libraries
  //=====================================================================
  
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
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
    chain = CreateAODChain("AODs.txt");
  }
  
  //---- Create the analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
  if(plugin) mgr->SetGridHandler(plugin);
  
  //  Input
  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  if (plugin) mgr->SetMCtruthEventHandler(mcHandler);
  
  AliAODInputHandler* iH = new AliAODInputHandler("handler","handler for my analisys");
  mgr->SetInputEventHandler(iH);
  

  //--------------------------------------------------------------
  // Other tasks

  // PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(kTRUE); // useMC
  // PID QA
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQATask = AddTaskPIDqa();
  
  gROOT->LoadMacro("./AliAnalysisTaskEfficiencydAOD.cxx+");
  gROOT->LoadMacro("./AddTaskEfficiencydAOD.C");
  AliAnalysisTaskEfficiencydAOD *task = AddTaskEfficiencydAOD();
  if(!task) ::Error("Task uninitialized.");

  
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
  plugin->SetROOTVersion("v5-34-08-7");
  plugin->SetAliROOTVersion("v5-06-00");
  plugin->SetAliPhysicsVersion("vAN-20150120");
  plugin->SetExecutableCommand("aliroot -b -q");
  
  // Declare input data to be processed.
  
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  
  plugin->SetGridDataDir("/alice/sim/2014/LHC14a6/"); //sim
  //plugin->SetGridDataDir("/alice/sim/2012/LHC12d3/"); //sim
  plugin->SetDataPattern("*/AliAOD.root"); // sim
  plugin->SetRunPrefix("");   // mc
  
  Int_t runlist[108] = {
    167915, 169099, 169858, 167920, 169138, 169859, 167985, 169144, 169923, 167987,
    169145, 169965, 167988, 169148, 170027, 168069, 169156, 170040, 168076, 169160,
    170081, 168105, 169167, 170083, 168107, 169238, 170084, 168108, 169411, 170085,
    168115, 169415, 170088, 168310, 169417, 170089, 168311, 169418, 170091, 168322,
    169419, 170155, 168325, 169420, 170159, 168341, 169475, 170163, 168342, 169498,
    170193, 168361, 169504, 170203, 168362, 169506, 170204, 168458, 169512, 170207,
    168460, 169515, 170228, 168464, 169550, 170230, 168467, 169553, 170268, 168511,
    169554, 170269, 168512, 169555, 170270, 168514, 169557, 170306, 168777, 169586,
    170308, 168826, 169587, 170309, 168988, 169588, 170311, 168992, 169590, 170312,
    169035, 169591, 170313, 169040, 169835, 170315, 169044, 169837, 170387, 169045,
    169838, 170388, 169091, 169846, 170572, 169094, 169855, 170593                  // LHC14a6
  };
  
  Int_t start = 0, nrun = 0;
  GetRunList(who,start,nrun);
  
  for(Int_t i = start; i < start + nrun; i++) {
    plugin->AddRunNumber(runlist[i]);
  }

  plugin->SetNrunsPerMaster(12);
  plugin->SetOutputToRunNo();
  
  
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("Efficiencyd2011");
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/taskname/out
  
  cout<<"-->>>>>>>>>>>>>>>>>>>>>>>>> 1"<<endl;
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisTaskEfficiencydAOD.cxx");
  
  cout<<"-->>>>>>>>>>>>>>>>>>>>>>>>> 2"<<endl;
  // plugin->SetAdditionalLibs("libTree.so libGeom.so libPhysics.so libVMC.so libMinuit.so libSTEERBase.so libESD.so libAOD.so libANALYSIS.so libOADB.so libANALYSISalice.so libCORRFW.so libPWGHFbase.so libPWGflowBase.so libPWGflowTasks.so libPWGHFvertexingHF.so");
  
  // plugin->SetAnalysisSource("AliAnalysisTaskESDMuonFilterO.cxx");
  //plugin->SetAnalysisSource("AliAODMuonReplicator0.cxx");
  
  //questo
  plugin->SetAdditionalLibs("libSTEERBase.so libESD.so AliAnalysisTaskEfficiencydAOD.h AliAnalysisTaskEfficiencydAOD.cxx libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so libOADB.so");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_PHYSICS/include -g");
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs();
  //plugin->SetOutputFiles("list.root");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s.C",taskname));
  
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(25);
  
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
