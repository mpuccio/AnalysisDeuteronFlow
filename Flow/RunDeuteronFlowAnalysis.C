enum anaModes {mLocal,mLocalPAR,mPROOF,mGrid,mGridPAR};
//mLocal: Analyze locally files in your computer using aliroot
//mLocalPAR: Analyze locally files in your computer using root + PAR files
//mPROOF: Analyze CAF files with PROOF
//mGrid: Analyze files on Grid via AliEn plug-in and using precompiled FLOW libraries
//       (Remark: When using this mode set also Bool_t bUseParFiles = kFALSE; in CreateAlienHandler.C)
//mGrid + par files: Analyze files on Grid via AliEn plug-in and using par files for FLOW package.
//                   Simply set Int_t mode = mGrid and Bool_t useFlowParFiles = kTRUE as arguments.

// CENTRALITY DEFINITION
Bool_t kUseCentrality = kFALSE;
Int_t binfirst = -1; //if kUseCentrality then change accordingly
Int_t binlast = -1;  //if kUseCentrality then change accordingly
const Int_t numberOfCentralityBins = 9;
Float_t centralityArray[numberOfCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.}; // in centrality percentile
//Int_t centralityArray[numberOfCentralityBins+1] = {41,80,146,245,384,576,835,1203,1471,10000}; // in terms of TPC only reference multiplicity

TString commonOutputFileName = "AnalysisResults"; // e.g.: result for centrality bin 0 will be in the file "outputCentrality0.root", etc


//void runFlowTask(Int_t mode=mLocal, Int_t nRuns = 10,
//Bool_t DATA = kFALSE, const Char_t* dataDir="/Users/snelling/alice_data/Therminator_midcentral", Int_t offset = 0)

//void runFlowTask(Int_t mode = mGridPAR, Int_t nRuns = 50000000,
//		 Bool_t DATA = kTRUE, const Char_t* dataDir="/alice/data/LHC10h_000137161_p1_plusplusplus", Int_t offset=0)
//void runFlowTask(Int_t mode = mLocal, Int_t nRuns = 50000000,
//		 Bool_t DATA = kTRUE, const Char_t* dataDir="./data/", Int_t offset=0)
//void runFlowTask(Int_t mode = mGridPAR, Bool_t DATA = kTRUE)

AliAnalysisGrid* CreateAlienHandler(Bool_t bUseParFiles=kFALSE);

void RunDeuteronFlowAnalysis(Int_t mode = mGrid,
                             Bool_t useFlowParFiles = kFALSE,
                             Bool_t DATA = kTRUE,
                             Bool_t useTender = kFALSE)
{
  // Time:
  TStopwatch timer;
  timer.Start();
  // Cross-check user settings before starting:
  //  CrossCheckUserSettings(DATA);
  // Load needed libraries:
  LoadLibraries(mode,useFlowParFiles);
  // Create and configure the AliEn plug-in:
  if(mode == mGrid || mode == mGridPAR)
  {
    AliAnalysisGrid *alienHandler = CreateAlienHandler(useFlowParFiles);
    if(!alienHandler) return;
  }
  // Chains:
  if(mode == mLocal || mode == mLocalPAR) {
    TChain *chain = new TChain("esdTree");
    //chain->Add("/home/pchrist/ALICE/HeavyIons/Data/137161/Set1/AliESDs.root");
    TChain* chain = CreateESDChain(dataDir, nRuns, offset);
    //TChain* chain = CreateAODChain(dataDir, nRuns, offset);
  }
  
  // Create analysis manager:
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
  // Connect plug-in to the analysis manager:
  if(mode == mGrid || mode == mGridPAR)
  {
    mgr->SetGridHandler(alienHandler);
  }
  
  // Event handlers:
  AliVEventHandler* esdH = new AliAODInputHandler;
  mgr->SetInputEventHandler(esdH);
  if (!DATA) {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mc);
  }
  
  // Task to check the offline trigger:
//  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
//  AddTaskPhysicsSelection(!DATA);
  
  //Add the centrality determination task
  if(kUseCentrality) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AddTaskCentrality();
  }
  
  //Add the TOF tender
  //gROOT->LoadMacro("$ALICE_ROOT/PWG/FLOW/macros/AddTaskTenderFlow.C");
  //AddTaskTenderFlow();
  
  // Setup analysis and usage of centrality bins
  gROOT->LoadMacro("AddTaskFlowCentralityPIDdeuteron.C");
  Float_t kLowCentralityBin = -1.;
  Float_t kHighCentralityBin = -1;
  if(kUseCentrality) {
    kLowCentralityBin = centralityArray[binfirst];
    kHighCentralityBin = centralityArray[binlast];
  }
  AddTaskFlowCentralityPID(kLowCentralityBin,
                           kHighCentralityBin,
                           commonOutputFileName);
  
  // Enable debug printouts:
  mgr->SetDebugLevel(2);
  // Run the analysis:
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(mode == mLocal || mode == mLocalPAR) {
    mgr->StartAnalysis("local",chain);
  } else if(mode == mPROOF) {
    mgr->StartAnalysis("proof",dataDir,nRuns,offset);
  } else if(mode == mGrid || mode == mGridPAR) {
    mgr->StartAnalysis("grid");
  }
  
  // Print real and CPU time used for analysis:
  timer.Stop();
  timer.Print();
  
} // end of void runFlowTask(...)

//===============================================================================================
/*
 void CrossCheckUserSettings(Bool_t bData)
 {
 // Check in this method if the user settings make sense.
 if(LYZ1SUM && LYZ2SUM) {cout<<" WARNING: You cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1 !!!!"<<endl; exit(0); }
 if(LYZ1PROD && LYZ2PROD) {cout<<" WARNING: You cannot run LYZ1 and LYZ2 at the same time! LYZ2 needs the output from LYZ1 !!!!"<<endl; exit(0); }
 if(LYZ2SUM && LYZEP) {cout<<" WARNING: You cannot run LYZ2 and LYZEP at the same time! LYZEP needs the output from LYZ2 !!!!"<<endl; exit(0); }
 if(LYZ1SUM && LYZEP) {cout<<" WARNING: You cannot run LYZ1 and LYZEP at the same time! LYZEP needs the output from LYZ2 !!!!"<<endl; exit(0); }
 } // end of void CrossCheckUserSettings()
 */
//===============================================================================================

void LoadLibraries(const anaModes mode, Bool_t useFlowParFiles )
{
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");
  
  if (mode==mLocal || mode==mGrid)
  {
    gSystem->Load("libSTEERBase");
    gSystem->Load("libCDB");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libRAWDatarec");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libSTEER");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libTPCbase");
    gSystem->Load("libTOFbase");
    gSystem->Load("libTOFrec");
    gSystem->Load("libTRDbase");
    gSystem->Load("libVZERObase");
    gSystem->Load("libVZEROrec");
    gSystem->Load("libT0base");
    gSystem->Load("libT0rec");
    gSystem->Load("libTENDER");
    gSystem->Load("libTENDERSupplies");
    
    if (useFlowParFiles)
    {
      AliAnalysisAlien::SetupPar("PWGflowBase");
      AliAnalysisAlien::SetupPar("PWGflowTasks");
    }
    else
    {
      gSystem->Load("libPWGflowBase");
      gSystem->Load("libPWGflowTasks");
    }
  }
  else if (mode==mPROOF)
  {
    TList* list = new TList();
    list->Add(new TNamed("ALIROOT_MODE", "ALIROOT"));
    if (useFlowParFiles)
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:TENDER:TENDERSupplies"));
    else
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "ANALYSIS:ANALYSISalice:TENDER:TENDERSupplies:PWGflowBase:PWGflowTasks"));
    
    //list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES","PWG/FLOW/Base:PWG/FLOW/Tasks"));
    
    // Connect to proof
    printf("*** Connect to PROOF ***\n");
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    //TProof* proof = TProof::Open("alice-caf.cern.ch");
    TProof* proof = TProof::Open("skaf.saske.sk");
    
    // list the data available
    //gProof->ShowDataSets("/*/*");
    //gProof->ShowDataSets("/alice/sim/"); //for MC Data
    //gProof->ShowDataSets("/alice/data/"); //for REAL Data
    
    proof->ClearPackages();
    proof->EnablePackage("VO_ALICE@AliRoot::v4-21-14-AN",list);
    
    if (useFlowParFiles)
    {
      gProof->UploadPackage("PWGflowBase.par");
      gProof->UploadPackage("PWGflowTasks.par");
    }
    
    // Show enables Packages
    gProof->ShowEnabledPackages();
  }
} // end of void LoadLibraries(const anaModes mode)

// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...
  
  if (!aDataDir)
    return 0;
  
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    printf("%s not found.\n", aDataDir);
    return 0;
  }
  
  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;
  
  if (flags & 2)
  {
    TString execDir(gSystem->pwd());
    TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
    TList* dirList            = baseDir->GetListOfFiles();
    Int_t nDirs               = dirList->GetEntries();
    gSystem->cd(execDir);
    
    Int_t count = 0;
    
    for (Int_t iDir=0; iDir<nDirs; ++iDir)
    {
      TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
      if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
        continue;
      
      if (offset > 0)
      {
        --offset;
        continue;
      }
      
      if (count++ == aRuns)
        break;
      
      TString presentDirName(aDataDir);
      presentDirName += "/";
      presentDirName += presentDir->GetName();
      chain->Add(presentDirName + "/AliESDs.root/esdTree");
      //  cerr<<presentDirName<<endl;
    }
    
  }
  else
  {
    // Open the input stream
    ifstream in;
    in.open(aDataDir);
    
    Int_t count = 0;
    
    // Read the input list of files and add them to the chain
    TString esdfile;
    while(in.good()) {
      in >> esdfile;
      if (!esdfile.Contains("root")) continue; // protection
      
      if (offset > 0)
      {
        --offset;
        continue;
      }
      
      if (count++ == aRuns)
        break;
      
      // add esd file
      chain->Add(esdfile);
    }
    
    in.close();
  }
  
  return chain;
  
} // end of TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)

//===============================================================================================

TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliAOD.root
  // <aDataDir>/<dir1>/AliAOD.root
  // ...
  
  if (!aDataDir)
    return 0;
  
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    printf("%s not found.\n", aDataDir);
    return 0;
  }
  
  TChain* chain = new TChain("aodTree");
  TChain* chaingAlice = 0;
  
  if (flags & 2)
  {
    TString execDir(gSystem->pwd());
    TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
    TList* dirList            = baseDir->GetListOfFiles();
    Int_t nDirs               = dirList->GetEntries();
    gSystem->cd(execDir);
    
    Int_t count = 0;
    
    for (Int_t iDir=0; iDir<nDirs; ++iDir)
    {
      TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
      if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 ||
          strcmp(presentDir->GetName(), "..") == 0)
        continue;
      
      if (offset > 0)
      {
        --offset;
        continue;
      }
      
      if (count++ == aRuns)
        break;
      
      TString presentDirName(aDataDir);
      presentDirName += "/";
      presentDirName += presentDir->GetName();
      chain->Add(presentDirName + "/AliAOD.root/aodTree");
      // cerr<<presentDirName<<endl;
    }
    
  }
  else
  {
    // Open the input stream
    ifstream in;
    in.open(aDataDir);
    
    Int_t count = 0;
    
    // Read the input list of files and add them to the chain
    TString aodfile;
    while(in.good()) {
      in >> aodfile;
      if (!aodfile.Contains("root")) continue; // protection
      
      if (offset > 0)
      {
        --offset;
        continue;
      }
      
      if (count++ == aRuns)
        break;
      
      // add aod file
      chain->Add(aodfile);
    }
    
    in.close();
  }
  
  return chain;
  
} // end of TChain* CreateAODChain(const char* aDataDir, Int_t aRuns, Int_t offset)

AliAnalysisGrid* CreateAlienHandler(Bool_t bUseParFiles) {
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init
  // then source /tmp/gclient_env_$UID in the current shell.
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  //plugin->SetRunMode("test");
  //plugin->SetRunMode("offline");
  //plugin->SetRunMode("submit");
  plugin->SetRunMode("test");
  //plugin->SetRunMode("terminate");
  plugin->SetNtestFiles(3); // Relevant only for run mode "test"
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-08");
  plugin->SetAliROOTVersion("vAN-20140915");
  plugin->SetExecutableCommand("aliroot -b -q");
  
  // Declare input data to be processed - can be done in two ways:
  // METHOD 1: Create automatically XML collections using alien 'find' command.
  // ============================================================================
  //  Example 1: MC production (set in macro runFlowTask.C: DATA = kFALSE)
  //plugin->SetGridDataDir("/alice/sim/LHC10d4");
  //plugin->SetDataPattern("*AliESDs.root"); // The default data pattern, other may be "*tag.root", "*ESD.tag.root", etc
  //plugin->AddRunNumber(119844); // Alternatively use e.g. plugin->SetRunRange(105044,106044); to add more runs in one go
  //plugin->SetOutputToRunNo();
  // ============================================================================
  //  Example 2: Real data (set in macro runFlowTask.C: DATA = kTRUE, MCEP = kFALSE)
  plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/");
  plugin->SetDataPattern("*/pass2/AOD145/*AliAOD.root");
  plugin->SetRunPrefix("000"); // IMPORTANT!
  int nrun = 62, start = 0;
  Int_t runlist[120] = {                                                                 // Counter
    170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170204, 170203,
    170193, 170163, 170159, 170155, 170081, 170027, 169859, 169858, 169855, 169846,
    169838, 169837, 169835, 169417, 169415, 169411, 169238, 169167, 169160, 169156,
    169148, 169145, 169144, 169138, 169094, 169091, 169035, 168992, 168988, 168826,
    168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361,
    168342, 168341, 168325, 168322, 168311, 168310, 167988, 167987,                     // from CF
    167902, 167903, 167909, 167986, 168066, 168068, 168103, 168104, 168212, 168461,
    168984, 169040, 169044, 169045, 169099, 169143, 169418, 169419, 169420, 169475,
    169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557,
    169584, 169586, 169587, 169588, 169590, 169591, 169922, 169956, 169975, 169981,
    170036, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152,
    170195, 170207, 170208, 170311, 170312, 170313, 170315, 170387, 170388, 170556,
    170572, 170593                                                                      // from HF
  };
  
  for(Int_t i = start; i < 1; i++) {
    plugin->AddRunNumber(runlist[i]);
  }
  // plugin->AddRunNumber(119844); // Alternatively use e.g. plugin->SetRunRange(104044,106044); to add more runs in one go
  plugin->SetOutputToRunNo();
  // ============================================================================
  
  // METHOD 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //plugin->AddDataFile("hijingWithoutFlow10000Evts.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("data");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes:
  // ... (if this is needed see in official tutorial example how to do it!)
  plugin->SetAnalysisSource("AliCustomPIDResponse.cxx");
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("AliCustomPIDResponse.h AliCustomPIDResponse.cxx");
  //plugin->SetAdditionalLibs("libCORRFW.so libTOFbase.so libTOFrec.so");
  if(!bUseParFiles)
  {
    plugin->SetAdditionalLibs("libGui.so libProof.so libMinuit.so libXMLParser.so "
                              "libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so "
                              "libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so "
                              "libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so "
                              "libT0rec.so libTENDER.so libTENDERSupplies.so "
                              "libPWGflowBase.so libPWGflowTasks.so");
  }
  else // load libs via par files
  {
    plugin->SetAdditionalLibs("libGui.so libProof.so libMinuit.so libXMLParser.so "
                              "libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so "
                              "libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so "
                              "libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so "
                              "libT0rec.so libTENDER.so libTENDERSupplies.so");
    plugin->EnablePackage("PWGflowBase.par");
    plugin->EnablePackage("PWGflowTasks.par");
  }
  // Do not specify your outputs by hand anymore:
  plugin->SetDefaultOutputs(kTRUE);
  // To specify your outputs by hand set plugin->SetDefaultOutputs(kFALSE); and comment in line plugin->SetOutputFiles("...");
  // and plugin->SetOutputArchive("..."); bellow.
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // plugin->SetOutputFiles("AnalysisResults.root");
  // Optionally define the files to be archived.
  // plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  // plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // plugin->SetOutputArchive("log_archive.zip:");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("flowAnalysis.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);
  // Optionally set number of runs per masterjob:
  plugin->SetNrunsPerMaster(1);
  // Optionally set overwrite mode. Will trigger overwriting input data colections AND existing output files:
  plugin->SetOverwriteMode(kTRUE);
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(100000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("flowAnalysis.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  return plugin;
}

