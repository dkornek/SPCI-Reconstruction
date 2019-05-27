// ROOT Macro to create simulated sources out of the BASE file

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-05-08

void getStats(TFile* f, Int_t& nDet, Int_t& nBin, Double_t& rMin, Double_t& rMax){
    // extract index of both the detectors from the name of the spectrum

    Int_t detectorD;
    Int_t detectorC;

    TString nameOfLastVoxel = f->GetListOfKeys()->Last()->GetName();
    TDirectory* dirOfVoxel = (TDirectory*)f->Get(nameOfLastVoxel);
    TString nameOfLastSpectrum = dirOfVoxel->GetListOfKeys()->Last()->GetName();
    nBin = ((TH1F*)dirOfVoxel->Get(nameOfLastSpectrum))->GetNbinsX();
    rMin = ((TH1F*)dirOfVoxel->Get(nameOfLastSpectrum))->GetXaxis()->GetXmin();
    rMax = ((TH1F*)dirOfVoxel->Get(nameOfLastSpectrum))->GetXaxis()->GetXmax();

    detectorD = ((TString)nameOfLastSpectrum(3, 2)).Atoi();
    detectorC = ((TString)nameOfLastSpectrum(5, 2)).Atoi();
    nDet = std::max(detectorD, detectorC) + 1;
}

TList* createEmptySpectraList(TFile* f){
    // create all spectra for the measurements file

    // get info of base file
    Int_t numberOfDetectors;
    Int_t numberOfBins;
    Double_t rangeMin;
    Double_t rangeMax;

    getStats(f,
             numberOfDetectors, numberOfBins,
             rangeMin, rangeMax);

    TList* list = new TList();
    for (Int_t d = 0; d < numberOfDetectors; ++d){
        for (Int_t c = 0; c < numberOfDetectors; ++c){

            if (d != c){
                TString nameOfSpectrum;
                nameOfSpectrum.Form("%02i%02i", d, c);

                TString description;
                description.Form("Sum energy spectrum Die %i with Die %i", d, c);
                TH1F* spectrum = new TH1F(nameOfSpectrum, description,
                                          numberOfBins, rangeMin, rangeMax);
                list->AddLast(spectrum);
            }
        }
    }

    return list;
}

void addSpectraToList(TFile* f, TList* l, std::vector<Int_t> loc, std::vector<Double_t> w){
    // Add all spectra of specified voxels to the measurements files

    TIter nextVoxel(f->GetListOfKeys());
    TKey* keyVoxel;
    while ((keyVoxel = (TKey*)nextVoxel())){
        // iterate through every voxel
        TString nameOfVoxel = keyVoxel->GetName();

        for (Int_t n = 0; n <= loc.size(); ++n){
            // iterate through all wanted sources

            if (nameOfVoxel.Atoi() == loc[n]){
                TDirectory* directoryOfVoxel = (TDirectory*)f->Get(nameOfVoxel);

                TIter nextSpectrum(directoryOfVoxel->GetListOfKeys());
                TKey* keySpectrum;
                while ((keySpectrum = (TKey*)nextSpectrum())){
                    // add all spectra of wanted source to the measurements spectra
                    TString nameOfSpectrum = keySpectrum->GetName();

                    // get spectrum of specified voxel
                    TH1F* spectrumInVoxel = (TH1F*)directoryOfVoxel->Get(nameOfSpectrum)->Clone("");
                    spectrumInVoxel->Scale(w[n]);

                    // get corresponding spectrum in measurements file
                    TH1F* spectrum = (TH1F*)l->FindObject((TString)nameOfSpectrum(3, 4));
                    spectrum->Add(spectrumInVoxel);

                    delete spectrumInVoxel;
                }
            }
        }
    }
}

void CreateSource(){

    // Open base file = system matrix
    TString pathToBase = "../folder/subfolder/*.root";
    pathToBase = "../../data/SystemMatrix/Original/SPCIBase441.root";
    TFile* input = new TFile(pathToBase, "READ");
    if (!input->IsOpen()){
        std::cout << "Input file not found!\n";
        return;
    }

    // Create source file = simulated measurements
    TFile* output = new TFile("Source_Difficult.root", "RECREATE");
    if (!output->IsOpen()){
        std::cout << "Output file could not be created!\n";
        return;
    }

    // Create all spectra
    TList* spectraList = createEmptySpectraList(input);

    // point source locations
    std::vector<Int_t> sourceLocations =    {  50,  51,  52,  53,  54,  69,  70,  71,  72,  73,
                                               74,  75,  76,  77,  90,  91,  92,  93,  94,  95,
                                               96,  97,  98, 109, 110, 111, 112, 113, 114, 115,
                                              116, 117, 118, 119, 120, 121, 130, 131, 132, 133,
                                              134, 135, 136, 137, 138, 139, 140, 141, 142, 149,
                                              150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                              160, 161, 162, 163, 164, 165, 170, 171, 172, 173,
                                              174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                                              184, 185, 186, 190, 191, 192, 193, 194, 195, 196,
                                              197, 198, 199, 200, 201, 202, 203, 204, 205, 206,
                                              207, 208, 211, 212, 213, 214, 215, 216, 217, 218,
                                              219, 220, 221, 222, 223, 224, 225, 226, 227, 228,
                                              229, 232, 233, 234, 235, 236, 237, 238, 239, 240,
                                              241, 242, 243, 244, 245, 246, 247, 248, 249, 250,
                                              254, 255, 256, 257, 258, 259, 260, 261, 262, 263,
                                              264, 265, 266, 267, 268, 269, 270, 275, 276, 277,
                                              278, 279, 280, 281, 282, 283, 284, 285, 286, 287,
                                              288, 289, 290, 291, 298, 299, 300, 301, 302, 303,
                                              304, 305, 306, 307, 308, 309, 310, 319, 320, 321,
                                              322, 323, 324, 325, 326, 327, 328, 329, 330, 331,
                                              342, 343, 344, 345, 346, 347, 348, 349, 350, 363,
                                              364, 365, 366, 367, 368, 369, 370, 371, 386, 387,
                                              388, 389, 390 };

    std::vector<Double_t> weights =         { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              1, 2, 2, 2, 2, 2, 2, 2, 5, 2,
                                              2, 2, 2, 2, 2, 2, 2, 1, 1, 1,
                                              2, 2, 2, 2, 2, 5, 5, 5, 2, 2,
                                              2, 2, 2, 2, 2, 2, 1, 1, 1, 2,
                                              2, 2, 2, 2, 5, 5, 5, 2, 2, 2,
                                              2, 2, 2, 2, 2, 1, 1, 1, 2, 2,
                                              2, 2, 2, 5, 5, 5, 2, 2, 2, 2,
                                              2, 2, 2, 2, 1, 2, 2, 2, 2, 2,
                                              2, 2, 5, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                              2, 2, 2 };

    if (sourceLocations.size() != weights.size()){
        std::cout << "Number of locations and weights do not match!\n";
        return;
    }

    // add specified sources
    addSpectraToList(input, spectraList, sourceLocations, weights);

    // save and close
    spectraList->Write();
    output->Close();
    input->Close();
    std::cout << "File written successfully." << std::endl;
}
