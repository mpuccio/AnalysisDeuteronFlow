int fkNTriggers = 5;

void BinLogAxis(const TH1 *h)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = const_cast<TAxis*>(h->GetXaxis());
  const Int_t bins = axis->GetNbins();
  
  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
  
  newBins[0] = from;
  Double_t factor = pow(to / from, 1. / bins);
  
  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i - 1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
}


void MergeFlowd() 
{
	TGrid *alien = TGrid::Connect("alien");
	if (alien->IsZombie())
	{
		delete alien;
		cout << "Fatal: Alien is a zombie!" << endl;
		return;
	}

  Int_t runlist[] = {                                                               // Counter
    170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170204, 170203,
    170193, 170163, 170159, 170155, 170081, 170027, 169859, 169858, 169855, 169846,
    169838, 169837, 169835, 169417, 169415, 169411, 169238, 169167, 169160, 169156,
    169148, 169145, 169144, 169138, 169094, 169091, 169035, 168992, 168988, 168826,
    168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361,
    168342, 168341, 168325, 168322, 168311, 168310, 167988, 167987
  };	

  TString dataBaseDir = "/alice/cern.ch/user/m/mpuccio/Flowd_PbPb2011/outputAOD/000";
 
  for (int iRun = 0; iRun < 58; ++iRun)
  {
    TFileMerger merger;
    merger.OutputFile(Form("$HOME/data/Run_%i.root",runlist[iRun]));
  	TString runDir = Form("%s%i",dataBaseDir.Data(),runlist[iRun]);
  	TGridResult *res = alien->Command(Form("find %s */mpuccio_Flowdnt.root",runDir.Data()));
  	TIter iter(res);
  	TMap *map;
  	while ((map = (TMap*)iter())) 
  	{
  		TObjString *obj = dynamic_cast<TObjString*>(map->GetValue("turl"));
  		if (!obj || !obj->String().Length())
  		{
  			delete res;
  			break;
  		}

  		TFile *f = TFile::Open(TFile::AsyncOpen(obj->String().Data(),"TIMEOUT=5"));
  		if (!f->IsOpen())
  		{
  			cout << "File " << obj->String().Data() << " has some problems..." << endl;
  			continue;
  		}
  		merger.AddFile(f);
  		merger.PartialMerge();
  		f->Close();
      delete f;
    }
    merger.Merge();
  }
}
