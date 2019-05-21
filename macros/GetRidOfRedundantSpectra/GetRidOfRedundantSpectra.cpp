// ROOT Macro to get rid of redundant spectra

void getStats(TFile* f, Int_t& nDet){
    // extract index of both the detectors from the name of the spectrum

    Int_t detectorD;
    Int_t detectorC;

    TString nameOfLastVoxel = f->GetListOfKeys()->Last()->GetName();
    TDirectory* dirOfVoxel = (TDirectory*)f->Get(nameOfLastVoxel);
    TString nameOfLastSpectrum = dirOfVoxel->GetListOfKeys()->Last()->GetName();

    detectorD = ((TString)nameOfLastSpectrum(3, 2)).Atoi();
    detectorC = ((TString)nameOfLastSpectrum(5, 2)).Atoi();
    nDet = std::max(detectorD, detectorC) + 1;
}

void GetRidOfRedundantSpectra(){

    // open system matrix file with redundant spectra
    TString pathToMeasurement = "../folder/subfolder/measurement.root";
    pathToMeasurement = "../../data/SystemMatrix/Original/SPCIBase49_240.root";
    TFile* input = new TFile(pathToMeasurement, "READ");
    if (!input->IsOpen()){
        std::cout << "Measurement file not found!\n";
        return;
    }

    // open file for saving system matrix
    TFile* output = new TFile("SPCIBase49_120.root", "RECREATE");
    if (!output->IsOpen()){
        std::cout << "Output file could not be created!\n";
        return;
    }

    // get number of detectors
    Int_t numberOfDetectors;
    getStats(input, numberOfDetectors);

    TIter nextVoxel(input->GetListOfKeys());
    TKey* keyVoxel;
    while ((keyVoxel = (TKey*)nextVoxel())){
        // iterate through every voxel
        TString nameOfVoxel = keyVoxel->GetName();
        TString titleOfVoxel = keyVoxel->GetTitle();
        TDirectory* subDir = output->mkdir(nameOfVoxel, titleOfVoxel);
        subDir->cd();

        TDirectory* directoryOfVoxel = (TDirectory*)input->Get(nameOfVoxel);
        TIter nextSpectrum(directoryOfVoxel->GetListOfKeys());
        TKey* keySpectrum;
        while ((keySpectrum = (TKey*)nextSpectrum())){

            TString nameOfSpectrum = keySpectrum->GetName();

            Int_t detectorD = ((TString)nameOfSpectrum(3, 2)).Atoi();
            Int_t detectorC = ((TString)nameOfSpectrum(5, 2)).Atoi();

            for (Int_t d = 0; d < numberOfDetectors; ++d){
                if (d == detectorD){
                    for (Int_t c = d + 1; c < numberOfDetectors; ++c){
                        if (c == detectorC){

                            TH1F* spectrum = (TH1F*)directoryOfVoxel->Get(nameOfSpectrum);
                            spectrum->Write();
                            delete spectrum;
                        }
                    }
                }
            }
        }
    }

    // quit
    output->Close();
    input->Close();
    std::cout << "File written successfully." << std::endl;
}
