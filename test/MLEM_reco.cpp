#include <iostream>
//#include <TFile.h>
//#include <TCollection.h>
//#include <TH3.h>
//#include <TString.h>

void getDetectorIndices(TString energySpectrum, Int_t& indexD, Int_t& indexC){
    TString stringD = energySpectrum(0, 2);
    TString stringC = energySpectrum(2, 4);

    indexD = stringD.Atoi();
    indexC = stringC.Atoi();
}

TH3I* getEmptyN_dcb(TFile& f, Int_t numberOfBins){
    // extract information needed for the matrix N_dcb

    // iterate through all energy spectra
    TIter next(f.GetListOfKeys());
    TKey* key;
    Int_t iIndexD, iIndexC;  // indices of detectors d and c

    while ((key = (TKey*)next())){
        getDetectorIndices(key->GetName(), iIndexD, iIndexC);
    }

    TH3I* emptyN_dcb = new TH3I("N_dcb", "Number of events in energy bin b of S_dc",
                                iIndexD + 1, 0, iIndexD,            // detector d in [0, 15]
                                iIndexD + 1, 0, iIndexD,            // detector c in [0, 15]
                                numberOfBins, 1, numberOfBins);     // number of bins in conditioned spectra S_dc

    delete key;
    return emptyN_dcb;
}

void fillN_dcb(TFile& f, TH3I& matrixN_dcb){
    // fill contents of cells in the matrix N_dcb

    // iterate through all energy spectra
    TIter next(f.GetListOfKeys());
    TKey* key;
    TString iNameOfS_dc;
    TH1* iS_dc = nullptr;
    Int_t iIndexD, iIndexC;  // indices of detectors d and c
    Int_t nBinsZ = matrixN_dcb.GetNbinsZ();

    while ((key = (TKey*)next())){
        iNameOfS_dc = key->GetName();
        getDetectorIndices(iNameOfS_dc, iIndexD, iIndexC);

        iS_dc = (TH1*)f.Get(iNameOfS_dc);  // extract conditioned spectrum
        for (Int_t bin = 1; bin <= nBinsZ; ++bin){
            matrixN_dcb.SetBinContent(iIndexD, iIndexC, bin, iS_dc->GetBinContent(bin));
        }
    }

    delete key;
    delete iS_dc;
}

void MLEM_reco(TString pathToMeasurements = "Position_x1.5_y-1.5.root",
               TString pathToProjections = "SPCIBase49.root",
               Int_t nBins = 1000){

    // prepare the matrix N_dcb according to the measurement data
    TFile measurementsFile(pathToMeasurements, "READ");
    TH3I* N_dcb = getEmptyN_dcb(measurementsFile, nBins);
    fillN_dcb(measurementsFile, *N_dcb);
    measurementsFile.Close();

    // prepare the matrices p_dcbv
    TFile projectionsFile(pathToProjections, "READ");
    projectionsFile.Close();
}
