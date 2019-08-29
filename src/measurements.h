// measurements.h

#pragma once
#include <algorithm>
#include <TH3.h>
#include <TKey.h>
#include <TFile.h>

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
