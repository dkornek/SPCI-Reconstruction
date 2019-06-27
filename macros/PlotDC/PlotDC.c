// Plot sensitivity of the image space of a fixed detector pair in the measurement space

void PlotDC(){

	TString pathToSystemMatrix = "../../data/SystemMatrix/Original/SPCIBase441.root";
	TFile* input = new TFile(pathToSystemMatrix, "READ");

	if (!input->IsOpen()){
		cout << "File not found!\nShut down macro.\n";	
	}

	TString dc = "0010";
	TString name;
	name.Form("Sensitivity of Image Space for %s", dc.Data());

	TH2F* sensitivityOfDC = new TH2F("sens", name,
					  21, -50, 50,
				          21, -50, 50);
	sensitivityOfDC->GetXaxis()->SetTitle("#it{x} / mm");
	sensitivityOfDC->GetYaxis()->SetTitle("#it{y} / mm");
	sensitivityOfDC->SetStats(kFALSE);

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

				Double_t binContent = 0;
				for (Int_t bin = 1;; ++bin){

					if (EDSInVoxel->IsBinOverflow(bin)){
						break;
					}

					binContent += EDSInVoxel->GetBinContent(bin);
				}

				sensitivityOfDC->SetBinContent(xVoxel, yVoxel,
							       binContent);
								
			}
		}		

		delete keyEDS;
	}
	delete keyVoxel;

	name.Form("%s.png", dc.Data());
	TCanvas* c = new TCanvas("c", "Plot", 600, 600);
	c->SetRightMargin(0.15);
	sensitivityOfDC->Draw("COLZ2");
	c->Update();
	c->SaveAs(name);

}
