// ROOT Macro to generate data sets of fictive sources based on a given system matrix (inverse crime)
// A possibility to add pseudo Poisson noise is offered.

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
        for (Int_t c = d + 1; c < numberOfDetectors; ++c){

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

void addSpectraToList(TFile* f, TList* l, std::vector<Int_t> loc, std::vector<Double_t> w, Bool_t addNoise){
    // Add all spectra of specified voxels to the measurements files

    int numberOfEvents = 0;

    TIter nextVoxel(f->GetListOfKeys());
    TKey* keyVoxel;
    while ((keyVoxel = (TKey*)nextVoxel())){
        // iterate through every voxel
        TString nameOfVoxel = keyVoxel->GetName();

        for (Int_t n = 0; n < loc.size(); ++n){
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
                    Int_t integral = Int_t(spectrumInVoxel->Integral());

                    Int_t nBin = spectrumInVoxel->GetNbinsX();

                    if (addNoise){
                        Double_t rMin = spectrumInVoxel->GetXaxis()->GetXmin();
                        Double_t rMax = spectrumInVoxel->GetXaxis()->GetXmax();

                        for (Int_t i = 0; i < 10; ++i){
                            TH1F* randomSpectrum = new TH1F("randSpec", "Random Spectrum", nBin, rMin, rMax);
                            randomSpectrum->FillRandom(spectrumInVoxel, integral);

                            for (Int_t b = 1; b <= nBin; ++b){
                                spectrumInVoxel->SetBinContent(b, randomSpectrum->GetBinContent(b));
                            }

                            delete randomSpectrum;
                        }
                    }

                    for (Int_t b = 1; b <= nBin; ++b){
                        Int_t binContent = Int_t(w[n] * spectrumInVoxel->GetBinContent(b));
                        spectrumInVoxel->SetBinContent(b, binContent);
                    }

                    // get corresponding spectrum in measurements file
                    TH1F* spectrum = (TH1F*)l->FindObject((TString)nameOfSpectrum(3, 4));
                    spectrum->Add(spectrumInVoxel);

                    numberOfEvents += spectrumInVoxel->Integral();

                    delete spectrumInVoxel;
                }
            }
        }
    }

    std::cout << "Total number of events:\t" << numberOfEvents << "\n";
}


void CreateSource(){

    // Open base file = system matrix
    TString pathToBase = "../folder/subfolder/*.root";
    pathToBase = "../../data/SystemMatrix/SPCIBase441_B50.root";
    TFile* input = new TFile(pathToBase, "READ");
    if (!input->IsOpen()){
        std::cout << "Input file not found!\n";
        return;
    }

    // Create source file = simulated measurements
    TFile* output = new TFile("TwoPoints_B50_Noise_5.root", "RECREATE");
    if (!output->IsOpen()){
        std::cout << "Output file could not be created!\n";
        return;
    }

    // Create all spectra
    TList* spectraList = createEmptySpectraList(input);

    // point source locations
    std::vector<Int_t> sourceLocations =    { 283, 157 };

    std::vector<Double_t> weights =         {   1,   1 };

    if (sourceLocations.size() != weights.size()){
        std::cout << "Number of locations and weights do not match!\n";
        return;
    }

    // add noise
    int choice;
    std::cout << "Do you want to add pseudo-noise? (0 = NO; 1 = YES): ";
    std::cin >> choice;

    if (choice == 1){
        // add specified sources
        addSpectraToList(input, spectraList, sourceLocations, weights, kTRUE);
        std::cout << "Pseudo-noise has been added.\n";
    } else {
        // add specified sources
        addSpectraToList(input, spectraList, sourceLocations, weights, kFALSE);
        std::cout << "No noise was added.\n";
    }

    // save and close
    spectraList->Write();
    output->Close();
    input->Close();
    std::cout << "File written successfully." << std::endl;
}
