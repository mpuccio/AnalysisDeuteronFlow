{
  TProof::Open("pod://");
  if (!gProof) return;

  TDSet *manual_dset = new TDSet("TTree", "deuterons");
  TString user[3] = {"m/mpuccio","m/masera","s/sbufalin"};
  for (int i = 0; i < 3; ++i)
  {
    TString ddset = Form("Find;BasePath=/alice/cern.ch/user/%s/NucleiPbPb2011test/output;FileName=*/root_archive.zip;Anchor=mpuccio_Flowdnt.root;Tree=/deuterons;Mode=local",user[i].Data());

    TFileCollection *fc = gProof->GetDataSet(ddset.Data());

  // Create TDSet manually (terrible hack): apparently PROOF is stupid
  // and cannot create it correctly
    TFileInfo *fi;
    TIter next( fc->GetList() );
    while (( fi = (TFileInfo *)next() )) {
      const char *fn = fi->GetCurrentUrl()->GetUrl();
      Printf("adding: %s", fn);
      manual_dset->Add( fn );
    }
  }

  const double kCent[4] = {0.,10.,30.,50.};
  const double kBins[29] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f,8.0f,10.f
  };
  TH2F *hBins = new TH2F("hBins","hBins",3,kCent,28,kBins);
  TH1D *hCuts = new TH1D("hCuts","hCuts",10,0,10);
  TNamed *taskTitle = new TNamed("taskName","deuteron3cent");
  enum cutsName {kEtaMin=1,kEtaMax,kYMin,kYMax,kTPCsig,kTPCchi2,kSPDrec,kDCAxy,kDCAz,kRecreate};
  double cutsA[10] = {-0.8,0.8,-0.5,0.5,70,4.,1,0.5,1,1};
  for (int i = 1; i < 11; ++i) {
    hCuts->SetBinContent(i,cutsA[i - 1]);
  }
  gProof->AddInput(hBins);
  gProof->AddInput(hCuts);
  gProof->AddInput(taskTitle);
  
  // Process the TDset
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  
  // Other cuts
  hCuts->SetBinContent(10,0.);
  // chi2
  hCuts->SetBinContent(5,3.5);
  taskTitle->SetTitle("deuteron3cent_chi0");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(5,4.5);
  taskTitle->SetTitle("deuteron3cent_chi1");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(5,5);
  taskTitle->SetTitle("deuteron3cent_chi2");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(5,5.5);
  taskTitle->SetTitle("deuteron3cent_chi3");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(5,6.);
  taskTitle->SetTitle("deuteron3cent_chi4");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(5,4.);
  // dca_z
  hCuts->SetBinContent(9,0.5);
  taskTitle->SetTitle("deuteron3cent_dcaz0");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(9,0.75);
  taskTitle->SetTitle("deuteron3cent_dcaz1");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(9,1.5);
  taskTitle->SetTitle("deuteron3cent_dcaz2");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(9,2.);
  taskTitle->SetTitle("deuteron3cent_dcaz3");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(9,1.);
  // tpc
  hCuts->SetBinContent(4,60);
  taskTitle->SetTitle("deuteron3cent_tpc0");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(4,65);
  taskTitle->SetTitle("deuteron3cent_tpc1");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(4,75);
  taskTitle->SetTitle("deuteron3cent_tpc2");
  gProof->Process(manual_dset, "AODSelector.cxx+g");
  hCuts->SetBinContent(4,80);
  taskTitle->SetTitle("deuteron3cent_tpc3");
  gProof->Process(manual_dset, "AODSelector.cxx+g");







}
