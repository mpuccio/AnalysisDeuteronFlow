{
	TProof::Open("");
	TChain *chain = new TChain("Deuterons","Deuterons");
	chain->Add("$HOME/data/LHC14a6_0.root");
	chain->Add("$HOME/data/LHC14a6_1.root");
	chain->SetProof();
	chain->Process("EfficiencySelector.cxx+");
}