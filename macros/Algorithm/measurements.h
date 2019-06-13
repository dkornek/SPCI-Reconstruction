// measurements.h

#pragma once
#include <algorithm>
#include <TH3.h>
#include <TKey.h>
#include <TFile.h>

#include "utilities.h"

// ##### MEASUREMENTS #####
class Measurements{
public:
    Measurements(const TString pathToMeasurements);
    ~Measurements();

    void createN_dcb();
    void fillN_dcb(const Bool_t normalized);

    Int_t numberOfDetectors;
    Int_t numberOfBins;

    TH3F* N_dcb = nullptr;  // contains all detected events in a 3D histogram

private:
    void getNumbers();

    TFile* measurementsFile = nullptr;
    TIter* nextS_dcMeasurements = nullptr;
    TKey* keyS_dcMeasurements = nullptr;
};

Measurements::Measurements(const TString pathToMeasurements) : numberOfDetectors(0), numberOfBins(0){
    // open the measurement *.root file

    measurementsFile = new TFile(pathToMeasurements, "READ");
    nextS_dcMeasurements = new TIter(measurementsFile->GetListOfKeys());
    getNumbers();
}

Measurements::~Measurements(){
    delete N_dcb;

    delete keyS_dcMeasurements;
    delete nextS_dcMeasurements;
    delete measurementsFile;
}

void Measurements::createN_dcb(){
    // create empty 3d histogram of N_dcb
    // x = detector d, y = detector c, z = bins in spectrum
    // histogram is based on the last spectrum

    N_dcb = new TH3F("N_dcb", "Number of events in bin b in S_dc",
                     numberOfDetectors, 0, numberOfDetectors,
                     numberOfDetectors, 0, numberOfDetectors,
                     numberOfBins, 0, numberOfBins);
}

void Measurements::fillN_dcb(const Bool_t normalized){
    // fill N_dcb with the number of events extracted from the spectra

    nextS_dcMeasurements->Reset();

    if (normalized){
        while ((keyS_dcMeasurements = (TKey*)nextS_dcMeasurements->Next())){
            // iterate through every measurement spectrum

            Int_t detectorD;
            Int_t detectorC;
            TString nameOfS_dc = keyS_dcMeasurements->GetName();
            Utilities::getDetectorIndices(nameOfS_dc, detectorD, detectorC);

            // get spectrum
            TH1F* S_dc = (TH1F*)measurementsFile->Get(nameOfS_dc);

            Double_t integral = S_dc->Integral();
            if (integral != 0){

                S_dc->Scale(1.0 / integral);
                for (Int_t bin = 1; bin <= numberOfBins; ++bin){
                    N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         S_dc->GetBinContent(bin));
                }

            } else{

                for (Int_t bin = 1; bin <= numberOfBins; ++bin){
                    N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         0);
                }
            }

            delete S_dc;
        }
    }

    else {
        while ((keyS_dcMeasurements = (TKey*)nextS_dcMeasurements->Next())){
            // iterate through every measurement spectrum

            Int_t detectorD;
            Int_t detectorC;
            TString nameOfS_dc = keyS_dcMeasurements->GetName();
            Utilities::getDetectorIndices(nameOfS_dc, detectorD, detectorC);

            // get spectrum
            TH1F* S_dc = (TH1F*)measurementsFile->Get(nameOfS_dc);

            if (S_dc->Integral() != 0){
                for (Int_t bin = 1; bin <= numberOfBins; ++bin){
                    N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         S_dc->GetBinContent(bin));
                }
            }

            delete S_dc;
        }
    }
}

void Measurements::getNumbers(){
    // extract number of detectors and bins

    // title of spectra is expected to be of "ddcc" format
    Int_t detD;
    Int_t detC;
    TString nameOfLastS_dc;
    nameOfLastS_dc = measurementsFile->GetListOfKeys()->Last()->GetName();
    Utilities::getDetectorIndices(nameOfLastS_dc, detD, detC);

    // get number of detectors
    numberOfDetectors = std::max(detD, detC) + 1;

    // get number of bins
    TH1F* lastS_dc = (TH1F*)measurementsFile->Get(nameOfLastS_dc);
    numberOfBins = lastS_dc->GetNbinsX();
    delete lastS_dc;
}
