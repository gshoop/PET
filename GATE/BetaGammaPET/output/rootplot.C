// This file is run using root. In order to execute the commands contained enter the command
// "root -l rootplot.C"
{
	TFile *_file0 = TFile::Open("pet.root");

	Hits->Draw("edep");

}
