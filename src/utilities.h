// utilities.h
// Some handy functions that are used by both ML-EM and OE
//
// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-08-27

// #pragma once
#include <iostream>
#include <TROOT.h>

namespace Utilities {
    Bool_t checkForSameNumberOfBins(const Int_t NbinsMeasurements, const Int_t NbinsSystemMatrix);
    void getDetectorIndices(const TString nameOfSpectrum, Int_t& d, Int_t& c);
    void getImageSpaceIndices(const TString titleOfVoxel, Int_t &x, Int_t &y, Int_t &z);
}
