{
  TProof::Open("pod://");
  if (!gProof) return;

  // OCCHIO: non e' il tuo dataset, e' una versione "ridotta"
  TString ddset = "Find;BasePath=/alice/cern.ch/user/m/mpuccio/Flowd_PbPb2011/output;FileName=*/root_archive.zip;Regexp=.*000170(2[3-9][0-9]|30[0-9]).*;Anchor=mpuccio_Flowdnt.root;Tree=/deuterons;Mode=local";

  TFileCollection *fc = gProof->GetDataSet(ddset.Data());

  // Create TDSet manually (terrible hack): apparently PROOF is stupid
  // and cannot create it correctly
  TDSet *manual_dset = new TDSet("TTree", "deuterons");
  TFileInfo *fi;
  TIter next( fc->GetList() );
  while (( fi = (TFileInfo *)next() )) {
    const char *fn = fi->GetCurrentUrl()->GetUrl();
    Printf("adding: %s", fn);
    manual_dset->Add( fn );
  }

  // Process the TDset
  gProof->Process(manual_dset, "DeutSelector.C++g");

}
