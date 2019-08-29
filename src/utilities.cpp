// utilities.cpp

#include "utilities.h"

Bool_t Utilities::checkForSameNumberOfBins(const Int_t NbinsMeasurements, const Int_t NbinsSystemMatrix){
    // image reconstruction is not possible unless the number of bins are identical

    if (NbinsMeasurements == NbinsSystemMatrix){
        return kTRUE;
    } else {
        std::cout << "Binning pattern is not consistent\n";
        return kFALSE;
    }
}

void Utilities::getDetectorIndices(const TString nameOfSpectrum, Int_t& d, Int_t& c){
    // extract index of both the detectors from the name of the spectrum

    d = ((TString)nameOfSpectrum(0, 2)).Atoi();
    c = ((TString)nameOfSpectrum(2, 2)).Atoi();
}

void Utilities::getImageSpaceIndices(const TString titleOfVoxel, Int_t &x, Int_t &y, Int_t &z){
    // extract location in image space from the title (=name of directory) of the voxel v

    // prepare title
    TString title = titleOfVoxel.Copy();
    title.Remove(0, 1);
    title.Replace(title.Length() - 1, 1, ',');

    // actual indices / coordinates
    std::vector<Int_t> location = { 0, 0, 0 };
    for (Int_t i = 0; i < 3; ++i){
        Int_t positionOfSeparator = title.First(',');
        Int_t lengthOfString = title.Last(',');

        location[i] =  ((TString)title(0, positionOfSeparator)).Atoi();
        title = title(positionOfSeparator + 1, lengthOfString - positionOfSeparator + 1);
    }

    // assign indices / coordinates
    x = location[0];
    y = location[1];
    z = location[2];
}
