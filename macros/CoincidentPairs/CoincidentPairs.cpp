// ROOT Macro to define the coincident pairs given all pairs

void getStats(TFile* f, Int_t& nDet, Int_t& nBin, Double_t& rMin, Double_t& rMax){
    // extract index of both the detectors from the name of the spectrum

    Int_t detectorD;
    Int_t detectorC;

    TString nameOfLastSpectrum = f->GetListOfKeys()->Last()->GetName();
    TH1F* lastSpectrum = (TH1F*)f->Get(nameOfLastSpectrum);

    nBin = lastSpectrum->GetNbinsX();
    rMin = lastSpectrum->GetXaxis()->GetXmin();
    rMax = lastSpectrum->GetXaxis()->GetXmax();

    detectorD = ((TString)nameOfLastSpectrum(0, 2)).Atoi();
    detectorC = ((TString)nameOfLastSpectrum(2, 2)).Atoi();
    nDet = std::max(detectorD, detectorC) + 1;

    delete lastSpectrum;
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

void addSpectraToList(TFile* f, TList* l, std::vector<Int_t> detID){

    TIter nextSpectrum(f->GetListOfKeys());
    TKey* keySpectrum;
    while ((keySpectrum = (TKey*)nextSpectrum())){
        // iterate through every spectrum
        TString nameOfSpectrum = keySpectrum->GetName();
        Int_t detectorD = ((TString)nameOfSpectrum(0, 2)).Atoi();
        Int_t detectorC = ((TString)nameOfSpectrum(2, 2)).Atoi();

        for (UInt_t d = 0; d <= detID.size(); ++d){
            // iterate through all wanted detectors
            if (detectorD == detID[d]){
                for (UInt_t c = 0; c <= detID.size(); ++c){
                    // iterate through all wanted detectors
                    if (detectorC == detID[c]){

                        // get spectrum in file
                        TH1F* spectrumInFile = (TH1F*)f->Get(nameOfSpectrum);

                        // get spectrum in list
                        TH1F* spectrumInList = (TH1F*)l->FindObject(nameOfSpectrum);
                        spectrumInList->Add(spectrumInFile);

                        delete spectrumInFile;
                    }
                }
            }
        }
    }
}

void CoincidentPairs(){

    // open measurements file
    TString pathToMeasurement = "../folder/subfolder/measurement.root";
    pathToMeasurement = "../../data/Measurements/SixPoints/SixPoints_B50_N9E8_Noise.root";
    TFile* input = new TFile(pathToMeasurement, "READ");
    if (!input->IsOpen()){
        std::cout << "Measurement file not found!\n";
        return;
    }

    // open file for saving new detector arrangement
    TFile* output = new TFile("SixPoints_B50_N9E8_Noise_5_7_13_15.root", "RECREATE");
    if (!output->IsOpen()){
        std::cout << "Output file could not be created!\n";
        return;
    }

    // choose detectors to be considered
    TList* spectraList = createEmptySpectraList(input);
    std::vector<Int_t> detectorID = { 5, 7, 13, 15 };  // change here!

    // add spectra to list
    addSpectraToList(input, spectraList, detectorID);

    // save data
    spectraList->Write();
    output->Close();
    input->Close();
    std::cout << "File written successfully." << std::endl;
}


// ##### ORIGINAL ARRANGEMENT === 4 x 4 DETECTOR #####
//^ y
//|
//3  7 11 15
//2  6 10 14
//1  5  9 13
//0  4  8 12
//-> x

// ##### 2 x 2 DETECTOR #####
// 5, 6, 9, 10

// ##### 3 x 3 DETECTOR #####
// 0, 1, 2, 4, 5, 6, 8, 9, 10

// ##### NEW DESIGN BROAD #####
// 0, 3, 12, 15

// ##### NEW DESIGN NARROW #####
// 0, 2, 8, 10
