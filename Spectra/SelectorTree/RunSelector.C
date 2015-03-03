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

  gSystem->Load("AODSelector.cxx++g");

  const double kCent[5] = {0.,10.,20.,40.,60.};
  const double kBins[29] = {
    0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,1.2f,1.4f,
    1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,3.2f,3.4f,
    3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f,8.0f,10.f
  };
  AODSelector *sel = new AODSelector();
  sel->SetPtBins(29,kBins);
  sel->SetCentBins(5,kCent);
  sel->SetOutputOption("deuterons3cent",kTRUE);
  // Process the TDset
  gProof->Process(manual_dset, sel);

}
