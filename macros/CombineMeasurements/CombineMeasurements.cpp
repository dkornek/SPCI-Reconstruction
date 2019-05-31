// ROOT macro to create a single TH2F Histogram out of all 1D-Measurement Histograms

TH2F* calculateMeasurementsHistogram(TFile* f){
    // Retrieve needed information
    // x = Detector pairs: 0 = 0001 -> 120 = 1415
    // y = Energy range

    TString nameOfLastSpectrum = f->GetListOfKeys()->Last()->GetName();
    TH1F* lastSpectrum = (TH1F*)f->Get(nameOfLastSpectrum);
    Int_t numberOfEnergyBins = lastSpectrum->GetNbinsX();
    Double_t eMin = lastSpectrum->GetXaxis()->GetXmin();
    Double_t eMax = lastSpectrum->GetXaxis()->GetXmax();
    delete lastSpectrum;

    Int_t numberOfPairs = 0;
    TIter next(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())){
        // iterate through every spectrum

        numberOfPairs += 1;
    }

    TH2F* measurements = new TH2F("measurements", "Measurements",
                                  numberOfPairs, 0, numberOfPairs,
                                  numberOfEnergyBins, eMin, eMax);

    numberOfPairs = 1;
    next.Reset();
    while ((key = (TKey*)next())){
        // iterate through every spectrum

        TString nameOfSpectrum = key->GetName();
        TH1F* spectrum = (TH1F*)f->Get(nameOfSpectrum);

        for (Int_t bin = 1; bin <= numberOfEnergyBins; ++bin){
            Double_t energyBinContent = spectrum->GetBinContent(bin);

            measurements->SetBinContent(numberOfPairs, bin,
                                        energyBinContent);
        }

        delete spectrum;
        ++numberOfPairs;
    }

    return measurements;
}

void CombineMeasurements(){

    // read input file = measurements
    TString pathToMeasurement = "../../data/Measurements/SPCIBase49/Bins50/SourceCross.root";
    TFile* input = new TFile(pathToMeasurement, "READ");

    if (input->IsOpen()){
        cout << "File opened successfully!\n";
    }
    else {
        cout << "File not found. Quit Application.\n";
    }

    TH2F* measurements = calculateMeasurementsHistogram(input);
    measurements->SetStats(kFALSE);

    TCanvas* mCanvas = new TCanvas("measurements", "Measurements", 700, 700);
    measurements->Draw("COLZ");
    mCanvas->Update();

//    input->Close();
}
