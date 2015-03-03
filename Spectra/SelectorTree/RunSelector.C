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
  AODSelector *sel = new AODSelector();
  sel->SetOutputOption("deuterons3cent",kTRUE);
  // Process the TDset
  gProof->Process(manual_dset, sel);

}
