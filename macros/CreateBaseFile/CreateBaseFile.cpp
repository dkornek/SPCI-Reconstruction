// ROOT Macro to create a base file out of the GEANT4 simulation

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-05-10

TString getNameOfSpectrum(TString nameOfDir, TString nameOfS){
    // form the name of the spectrum according to the needs of the reconstruction

    Int_t d = ((TString)nameOfS(0, 2)).Atoi();
    Int_t c = ((TString)nameOfS(4, 2)).Atoi();

    TString newName;
    newName.Form("%s%02i%02i",
                 nameOfDir.Data(), d, c);
    return newName;
}

TString getDescriptionOfSpectrum(TString nameOfDir, TString titleOfDir, TString nameOfS){
    // form the description of the spectrum according to already used convention

    Int_t d = ((TString)nameOfS(0, 2)).Atoi();
    Int_t c = ((TString)nameOfS(4, 2)).Atoi();

    TString newDescription;
    newDescription.Form("#Delta#it{E} of die %i and die %i in Voxel %s %s",
                        d, c, nameOfDir.Data(), titleOfDir.Data());
    return newDescription;
}

void getVoxelIndices(TString titleOfVoxel, Int_t& x, Int_t& y, Int_t& z){
    // find out the location in the image space of a given voxel

    // actual coordinates
    std::vector<Double_t> location = { 0, 0, 0 };

    TString title = titleOfVoxel.Copy();
    title.Remove(0, 1);
    title.Replace(title.Length() - 1, 1, ',');

    Int_t positionOfSeparator;
    Int_t lengthOfString;
    for (Int_t i = 0; i < 3; ++i){
        positionOfSeparator = title.First(',');
        lengthOfString = title.Last(',');

        location[i] =  ((TString)title(0, positionOfSeparator)).Atof();
        title = title(positionOfSeparator + 1, lengthOfString - positionOfSeparator + 1);
    }

    Double_t voxelLength = 5.0;
    x = location[0] / voxelLength + 10;
    y = abs(location[1] / voxelLength - 10);
    // z = location[2] / voxelLength;
    z = 0.0;
}

void createVoxelList(TFile* fIn, TFile* fOut){
    // create the subdirectories for the voxels and fill them with the y-projections

    // voxel data
    Int_t nameOfVoxel;
    TString nameOfDirectory;
    TString titleOfDirectory;
    TDirectory* voxelDirectory;

    Int_t indexX;
    Int_t indexY;
    Int_t indexZ;

    TH2F* data;

    // spectrum data
    TString nameOfSpectrum;
    TString descriptionOfSpectrum;
    TH1F* spectrum;

    TIter nextVoxel(fIn->GetListOfKeys());
    TKey* keyVoxel;
    while ((keyVoxel = (TKey*)nextVoxel())){
        // get data of voxel
        voxelDirectory = (TDirectory*)fIn->Get(keyVoxel->GetName());

        // create name of subdirectory
        nameOfVoxel = ((TString)keyVoxel->GetName()).Atoi();
        nameOfDirectory.Form("%03i", nameOfVoxel);

        // create title of subdirectory
        getVoxelIndices(keyVoxel->GetTitle(), indexX, indexY, indexZ);
        titleOfDirectory.Form("(%i,%i,%i)", indexX, indexY, indexZ);

        // create subdirectory
        TDirectory* subDir = fOut->mkdir(nameOfDirectory, titleOfDirectory);
        subDir->cd();

        TIter nextSpectrum(voxelDirectory->GetListOfKeys());
        TKey* keySpectrum;
        while ((keySpectrum = (TKey*)nextSpectrum())){
            nameOfSpectrum = keySpectrum->GetName();

            if (nameOfSpectrum.BeginsWith('G')){
                continue;
            }

            // get base data for spectrum
            data = (TH2F*)voxelDirectory->Get(nameOfSpectrum);

            // get projection
            descriptionOfSpectrum = getDescriptionOfSpectrum(nameOfDirectory, titleOfDirectory, nameOfSpectrum);
            nameOfSpectrum = getNameOfSpectrum(nameOfDirectory, nameOfSpectrum);
            spectrum = (TH1F*)data->ProjectionY(nameOfSpectrum.Data(), 661, 662);
            spectrum->SetTitle(descriptionOfSpectrum.Data());
            spectrum->Write();

            delete data;
            delete spectrum;
        }
    }
}

void CreateBaseFile(){

    // open simulation file
    TString pathToSimulation = "../folder/subfolder/simulation.root";
    TFile* input = new TFile(pathToSimulation, "READ");
    if (!input->IsOpen()){
        std::cout << "Simulation file not found!\n";
        return;
    }

    // create new system matrix file
    TString pathToProjection = "../folder/subfolder/basefile.root";
    TFile* output = new TFile(pathToProjection, "RECREATE");
    if (!output->IsOpen()){
        std::cout << "Output file could not be created!\n";
        return;
    }

    createVoxelList(input, output);

    // save data and close
    output->Close();
    input->Close();
    std::cout << "File written successfully." << std::endl;
}
