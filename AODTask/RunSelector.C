{
  TProof::Open("pod://");
  if (!gProof) return;

  TDSet *manual_dset = new TDSet("TTree", "deuterons");
  TString user[4] = {"m/mpuccio","s/scapodic","e/ebruna","s/strogolo"}; 
  for (int i = 0; i < 4; ++i)
  {
    TString ddset = Form("Find;BasePath=/alice/cern.ch/user/%s/NucleiPbPb2011/output;FileName=*/root_archive.zip;Anchor=mpuccio_Flowdnt.root;Tree=/deuterons;Mode=local",user[i].Data());

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

  // Process the TDset
  gProof->Process(manual_dset, "AODSelector.cxx++g");

}
