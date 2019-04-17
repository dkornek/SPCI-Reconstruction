/* General comments:
 *
 * Measurement data should not be redundant for performance reasons
 * bins in spectrum in measurement = bins in spectrum in projection
 */


#include <iostream>
#include <vector>
#include <TFile.h>
#include <TCollection.h>
#include <TH3.h>
#include <TString.h>
#include <TKey.h>

using namespace std;


void getDetectorIndices(TString energySpectrum, Int_t& indexD, Int_t& indexC){
    TString stringD = energySpectrum(0, 2);
    TString stringC = energySpectrum(2, 4);

    indexD = stringD.Atoi();
    indexC = stringC.Atoi();
}

void getImageSpaceIndices(TString position, Int_t& iX, Int_t& iY, Int_t& iZ){
    TString stringX = position(1);
    TString stringY = position(3);
    TString stringZ = position(5);

    iX = stringX.Atoi();
    iY = stringY.Atoi();
    iZ = stringZ.Atoi();
}

TH3I* createEmptyN_dcbArray (TFile& f){
    // extract information needed for the matrix N_dcb

    // analyse name of last histogram in the file
    // convention of spectra titles e.g. for d = 7 and c = 11: S_dc = 0711
    Int_t iIndexD, iIndexC;  // indices of detectors d and c
    TString iNameOfLastSpectrum = f.GetListOfKeys()->Last()->GetName();

    getDetectorIndices(iNameOfLastSpectrum, iIndexD, iIndexC);
    Int_t nDetectors = max(iIndexC, iIndexD);
//    Int_t nBins = ((TH1*)f.Get(iNameOfLastSpectrum))->GetNbinsX();
    Int_t nBins = 250;

    // construct the 3d histogram of given size
    // this array will contain "e = 0.5 * nDetectors * (nDetectors - 1) * numberOfBins" elements, rest will be zeros
    TH3I* emptyN_dcbArray = new TH3I("N_dcb", "Number of events in energy bin b of S_dc",
                                     nDetectors + 1, 0, nDetectors,      // detector d
                                     nDetectors + 1, 0, nDetectors,      // detector c
                                     nBins, 1, nBins);                   // number of bins in conditioned spectra S_dc
    return emptyN_dcbArray;
}

void fillN_dcbArray(TFile& f, TH3I& N_dcbArr){
    // fill contents of cells in the matrix N_dcb

    // iterate through all energy spectra
    TIter next(f.GetListOfKeys());
    TKey* key;
    Int_t iIndexD, iIndexC;  // indices of detectors d and c
    TString iNameOfS_dc;  // string representation for title of spectrum
    TH1* iS_dc = nullptr;  // spectrum
//    Int_t nBinsZ = N_dcbArr.GetNbinsZ();  // number of bins

    while ((key = (TKey*)next())){
        iNameOfS_dc = key->GetName();
        getDetectorIndices(iNameOfS_dc, iIndexD, iIndexC);

        iS_dc = (TH1*)f.Get(iNameOfS_dc);  // extract conditioned spectrum
        iS_dc->Rebin(4);

//        for (Int_t bin = 1; bin <= nBinsZ; ++bin){
        for (Int_t bin = 1; bin <= 250; ++bin){
            N_dcbArr.SetBinContent(iIndexD, iIndexC, bin, iS_dc->GetBinContent(bin));
        }
    }

    delete key;
//    delete iS_dc;
}

TH3F* createEmptyA_vArray(TFile& f){
    // extract information needed for the activity distribution A_v
    // the initial activity distribution will be homogeneous = every cell will be 1;

    Int_t iImageSpaceX, iImageSpaceY, iImageSpaceZ;  // indices of image space
    TString iLastPosition = f.GetListOfKeys()->Last()->GetTitle();  // last folder in file will determine image space volume

    getImageSpaceIndices(iLastPosition, iImageSpaceX, iImageSpaceY, iImageSpaceZ);
    TH3F* emptyA_vArray = new TH3F("A_v", "Activity distribution",
                                    iImageSpaceX + 1, 0, iImageSpaceX,    // x axis
                                    iImageSpaceY + 1, 0, iImageSpaceY,    // y axis
                                    iImageSpaceZ + 1, 0, iImageSpaceZ);   // z axis

    return emptyA_vArray;
}

void setHomogeneousActivity(TFile& f, TH3F& A_vArr, Int_t content){
    // make activity distribution homogeneous

    // fill every cell with content
    TIter next(f.GetListOfKeys());
    TKey* key;

    Int_t iImageSpaceX, iImageSpaceY, iImageSpaceZ;  // indices of image space
    TString iPosition;

    while((key = (TKey*)next())){
        iPosition = key->GetTitle();
        getImageSpaceIndices(iPosition, iImageSpaceX, iImageSpaceY, iImageSpaceZ);
        A_vArr.SetBinContent(iImageSpaceX, iImageSpaceY, iImageSpaceZ, content);
    }

    delete key;
}

/*
vector<Int_t> calculateEmissionsFromVoxel(TFile& f){
    // data for iterating all voxels v in the image space
    TIter iVoxelNext(f.GetListOfKeys());
    TKey* iVoxelKey = nullptr;
    TString iVoxelName;
    TDirectory* iVoxelData = nullptr;

    // data for iterating all spectra in voxel v
    TKey* iSpectrumKey = nullptr;
    TString iSpectrumName;
    TH1* iSpectrum_dc = nullptr;

    Int_t iEmissionsFromVoxel;
    Int_t iEmissionsInSpectrum;
    vector<Int_t> iEmissionsFromVoxelVector;

    while ((iVoxelKey = (TKey*)iVoxelNext())){
        // iterate through all voxels

        iVoxelName = iVoxelKey->GetName();
        iVoxelData = (TDirectory*)f.Get(iVoxelName);  // all spectra for voxel v
        TIter iSpectrumNext(iVoxelData->GetListOfKeys());
        iEmissionsFromVoxel = 0;

        while ((iSpectrumKey = (TKey*)iSpectrumNext())){
            // iterate through all voxels; calculate emissions from all voxels by cumulative sum

            iSpectrumName = iSpectrumKey->GetName();
            iSpectrum_dc = (TH1*)iVoxelData->Get(iSpectrumName);
            iEmissionsInSpectrum = iSpectrum_dc->Integral();
            iEmissionsFromVoxel += iEmissionsInSpectrum;
        }

        iEmissionsFromVoxelVector.push_back(iEmissionsFromVoxel);
    }

    delete iVoxelKey;
    delete iSpectrumKey;
    delete iSpectrum_dc;

    return iEmissionsFromVoxelVector;
}

TList* calculateP_dcbv(TFile& f, vector<Int_t> iEmissionVector){
    // calculate the probabilities p_dcbv

    // data for iterating all voxels v in the image space
    TIter iVoxelNext(f.GetListOfKeys());
    TKey* iVoxelKey = nullptr;
    TString iVoxelName;
    TDirectory* iVoxelData = nullptr;

    // data for iterating all spectra in voxel v
    TKey* iSpectrumKey = nullptr;
    TString iSpectrumName;
    TH1* iSpectrum_dc = nullptr;

    TList* iP_dcbvList = new TList();  // all p_dcbv = scaled energy spectra
    Int_t iEmissionVectorIndex = 0;

    while ((iVoxelKey = (TKey*)iVoxelNext())){
        // iterate through all voxels

        iVoxelName = iVoxelKey->GetName();
        iVoxelData = (TDirectory*)f.Get(iVoxelName);  // all spectra for voxel v
        TIter iSpectrumNext(iVoxelData->GetListOfKeys());

        while ((iSpectrumKey = (TKey*)iSpectrumNext())){
            // iterate through all voxels; calculate emissions from all voxels by cumulative sum

            iSpectrumName = iSpectrumKey->GetName();
            iSpectrum_dc = (TH1*)iVoxelData->Get(iSpectrumName);
            iSpectrum_dc->Scale(1.0 / iEmissionVector[iEmissionVectorIndex]);
            iP_dcbvList->AddLast(iSpectrum_dc);
        }

        iEmissionVectorIndex += 1;
    }

    delete iVoxelKey;
    delete iSpectrumKey;
//    delete iSpectrum_dc;

    return iP_dcbvList;
}
*/

TList* createProjectionMatrices(TFile& f){
    // create all the projection matrices p_dcbv

//    vector<Int_t> iEmissionsFromVoxelVector = calculateEmissionsFromVoxel(f);
//    TList* iP_dcbvList = calculateP_dcbv(f, iEmissionsFromVoxelVector);

    // data for iterating all voxels v in the image space
    TIter folderNext(f.GetListOfKeys());
    TKey* folderKey = nullptr;
    TString folderName;

    // data for iterating all spectra in voxel v
    TDirectory* voxelData = nullptr;
    TKey* spectrumKey = nullptr;
    TString spectrumName;
    TH1* iS_dc = nullptr;

    Int_t emissionsFromVoxel;
    Int_t emissionsInSpectrum;

    TList* p_dcbv = new TList();  // all p_dcbv = scaled energy spectra

    while ((folderKey = (TKey*)folderNext())){
        // iterate through all voxels v

        folderName = folderKey->GetName();  // voxel name
        voxelData = (TDirectory*)f.Get(folderName);  // energy spectra for voxel v
        TIter spectrumNext(voxelData->GetListOfKeys()), spectrumNext2(voxelData->GetListOfKeys());
        emissionsFromVoxel = 0;

        while ((spectrumKey = (TKey*)spectrumNext())){
            // iterate through all conditioned spectra
            // calculate all emissions from voxel v

            spectrumName = spectrumKey->GetName();
            iS_dc = (TH1*)voxelData->Get(spectrumName);
            emissionsInSpectrum = iS_dc->Integral();
            emissionsFromVoxel += emissionsInSpectrum;
        }

        while ((spectrumKey = (TKey*)spectrumNext2())){
            // iterate through all conditioned spectra again
            // calculate all p_dcbv

            spectrumName = spectrumKey->GetName();
            iS_dc = (TH1*)voxelData->Get(spectrumName);
            iS_dc->Scale(1.0 / emissionsFromVoxel);
            p_dcbv->AddLast(iS_dc);
        }
    }

//    TH1F* test = (TH1F*)p_dcbv->FindObject("0481314");
//    cout << test->GetBinContent(80) << "\n";

    delete folderKey;
    delete spectrumKey;

    return p_dcbv;
}

//void calculateActivity(TH3F* iActivityArray, TH3I* iMeasurementsArray, TList* iProbabilities){

//    TH1F* test = (TH1F*)iProbabilities->FindObject("0481314");
//    cout << test->GetBinContent(80) << "\n";

//    // get number of voxels in image space
////    Int_t iImageSpaceX, iImageSpaceY, iImageSpaceZ;
////    TString iLastPosition = f.GetListOfKeys()->Last()->GetTitle();
////    getImageSpaceIndices(iLastPosition, iImageSpaceX, iImageSpaceY, iImageSpaceZ);
////    Int_t nVoxels = (iImageSpaceX + 1) * (iImageSpaceY + 1) * (iImageSpaceZ + 1);
//    Int_t nVoxels = 49;

//    // get number of bins
////    TDirectoryFile* iLastFolder = (TDirectoryFile*)f.Get(f.GetListOfKeys()->Last()->GetName());
////    Int_t nBins = ((TH1*)iLastFolder->Get(iLastFolder->GetListOfKeys()->Last()->GetName()))->GetNbinsX();
//    Int_t nBins = 250;
////    cout << test->GetNbinsX() << "\n";

//    cout << "Bins in N_dcb: " << iMeasurementsArray->GetNbinsZ() << ", Bins in p_dcv: " << test->GetNbinsX() << "\n";



//}

void MLEM_reco(TString pathToMeasurements = "Position_x1.5_y-1.5.root",
               TString pathToProjections = "SPCIBase49.root"){

    auto t1 = std::chrono::steady_clock::now();

    // prepare the matrix N_dcb according to the measurement data
    TFile measurementsFile(pathToMeasurements, "READ");
    TH3I* N_dcbArray = createEmptyN_dcbArray(measurementsFile);
    fillN_dcbArray(measurementsFile, *N_dcbArray);
    cout << N_dcbArray->GetBinContent(6, 10, 125) << endl;

    // prepare the initial activity distribution
    TFile projectionsFile(pathToProjections, "READ");
    TH3F* A_vArray = createEmptyA_vArray(projectionsFile);
    setHomogeneousActivity(projectionsFile, *A_vArray, 1);  // fill every cell with 1
    cout << A_vArray->GetBinContent(5, 5, 0) << endl;

    // prepare the projection matrices p_dcbv
    TList* p_dcbvList = createProjectionMatrices(projectionsFile);

//    TH1F* test = (TH1F*)p_dcbvList->FindObject("0481314");
//    cout << test->GetBinContent(80) << "\n";

//    calculateActivity(A_vArray, N_dcbArray, p_dcbvList);


    TIter next(projectionsFile.GetListOfKeys());
    TKey* key;

    Int_t iImageSpaceX, iImageSpaceY, iImageSpaceZ;  // indices of image space
    TString iPosition;

    Int_t binContentInActivity;
    Int_t voxel = 0;
    TString nameP_dcbv;

    TH1F* p_dcbSpectrum = nullptr;
    Double_t p_dcbv;
    Double_t n_dcb;
    Double_t correctionFactor = 0;
    Double_t correctionFactorCumulative = 0;
    Double_t newBinContent;
//    p_dcbSpectrum = (TH1F*)p_dcbvList->FindObject("0481514");

    while((key = (TKey*)next())){
        // iterate through every voxel in the image space

        iPosition = key->GetTitle();
        getImageSpaceIndices(iPosition, iImageSpaceX, iImageSpaceY, iImageSpaceZ);
        binContentInActivity = A_vArray->GetBinContent(iImageSpaceX, iImageSpaceY, iImageSpaceZ);
        correctionFactorCumulative = 0;

//        if ((voxel % 7) == 0){
//            cout << "\n";
//        }

//        cout << binContentInActivity << " ";


        for (Int_t d = 0; d <= 15; ++d){
            for (Int_t c = 0; c <= 15; ++c){

                if (d != c){

                    nameP_dcbv.Form("%03i%02i%02i", voxel, d, c);
                    p_dcbSpectrum = (TH1F*)p_dcbvList->FindObject(nameP_dcbv);

                    for (Int_t bin = 1; bin <= 250; ++bin){
                        p_dcbv = p_dcbSpectrum->GetBinContent(bin);
                        n_dcb = N_dcbArray->GetBinContent(d, c, bin);
                        correctionFactor = p_dcbv * n_dcb;
                        correctionFactorCumulative += correctionFactor;

                        cout << nameP_dcbv << ", Bin: " << bin << " " << correctionFactorCumulative << "\n";
                    }
                }
            }
        }

        voxel += 1;
    }

    projectionsFile.Close();
    measurementsFile.Close();

    auto t2 = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    cout << "\n" << elapsedTime << endl;
}

#ifndef __CINT__
int main(){
    MLEM_reco();
    return 0;
}
#endif
