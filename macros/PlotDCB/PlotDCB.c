// Plot sensitivity of the image space of a fixed bin in the measurement space

void PlotDCB(){

	TString pathToSystemMatrix = "../../data/SystemMatrix/Original/SPCIBase441.root";
	TFile* input = new TFile(pathToSystemMatrix, "READ");

	if (!input->IsOpen()){
		cout << "File not found!\nShut down macro.\n";	
	}

	TString dc = "0510";
	Int_t b = 650;
	TString name;
	name.Form("Sensitivity of Image Space for %s_%03i", dc.Data(), b);

	TH2F* sensitivityOfDCB = new TH2F("sens", name,
					  21, -50, 50,
				          21, -50, 50);
	sensitivityOfDCB->GetXaxis()->SetTitle("#it{x} / mm");
	sensitivityOfDCB->GetYaxis()->SetTitle("#it{y} / mm");
	sensitivityOfDCB->SetStats(kFALSE);

	Int_t yVoxel = 0;
	TIter nextVoxel(input->GetListOfKeys());
	TKey* keyVoxel;
	while ((keyVoxel = (TKey*)nextVoxel())){

		TString nameOfVoxel = keyVoxel->GetName();
		TDirectory* directoryOfVoxel = (TDirectory*)input->Get(nameOfVoxel);

		Int_t globalVoxel = nameOfVoxel.Atoi();
		Int_t xVoxel = globalVoxel % 21 + 1;

		if ((globalVoxel % 21) == 0){
			++yVoxel;
		}

		TIter nextEDS(directoryOfVoxel->GetListOfKeys());
		TKey* keyEDS;
		while ((keyEDS = (TKey*)nextEDS())){

			TString nameOfEDS = keyEDS->GetName();
			TString nameOfEDSReduced = (TString)nameOfEDS(3, 4);

			if (nameOfEDSReduced == dc){

				TH1F* EDSInVoxel = (TH1F*)directoryOfVoxel->Get(nameOfEDS);
				sensitivityOfDCB->SetBinContent(xVoxel, yVoxel,
							        EDSInVoxel->GetBinContent(b));				
			}
		}		

		delete keyEDS;
	}
	delete keyVoxel;

	name.Form("%s_%03i.png", dc.Data(), b);
	TCanvas* c = new TCanvas("c", "Plot", 600, 600);
	c->SetRightMargin(0.15);
	sensitivityOfDCB->Draw("COLZ2");
	c->Update();
	c->SaveAs(name);

}
