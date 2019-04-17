/* General comments:
 *
 * Measurement data should not be redundant for performance reasons
 * bins in spectrum in measurement = bins in spectrum in projection
 */


//#include <iostream>
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
    Int_t nBins = ((TH1*)f.Get(iNameOfLastSpectrum))->GetNbinsX();

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
    Int_t nBinsZ = N_dcbArr.GetNbinsZ();  // number of bins

    while ((key = (TKey*)next())){
        iNameOfS_dc = key->GetName();
        getDetectorIndices(iNameOfS_dc, iIndexD, iIndexC);

        iS_dc = (TH1*)f.Get(iNameOfS_dc);  // extract conditioned spectrum
        for (Int_t bin = 1; bin <= nBinsZ; ++bin){
            N_dcbArr.SetBinContent(iIndexD, iIndexC, bin, iS_dc->GetBinContent(bin));
        }
    }

    delete key;
    delete iS_dc;
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

void createProjectionMatrices(TFile& f){
    // create all the projection matrices p_dcbv
    // let p_dcb = (p_dcb1, p_dcb2, ... p_dcbv)
    // all p_dcb will be stored in a TCollection object, containing TCollection objects of p_dcbv

    // get number of bins
    TDirectoryFile* iLastFolder = (TDirectoryFile*)f.Get(f.GetListOfKeys()->Last()->GetName());
    Int_t nBins = ((TH1*)iLastFolder->Get(iLastFolder->GetListOfKeys()->Last()->GetName()))->GetNbinsX();

    // get number of voxels in image space
    Int_t iImageSpaceX, iImageSpaceY, iImageSpaceZ;
    TString iLastPosition = f.GetListOfKeys()->Last()->GetTitle();
    getImageSpaceIndices(iLastPosition, iImageSpaceX, iImageSpaceY, iImageSpaceZ);
    Int_t nVoxels = (iImageSpaceX + 1) * (iImageSpaceY + 1) * (iImageSpaceZ + 1);

    // calculate all p_dcbv
    TIter folderNext(f.GetListOfKeys());
    TKey* folderKey;
    TString folderName;

    TDirectory* voxelData;
    TKey* spectrumKey;
    TString spectrumName;
    TH1* iS_dc = nullptr;

    Int_t emissionsFromVoxel;
    Int_t emissionsInSpectrum;

    TH3F* p_dcbv = nullptr;
    TString iTitleSystem;
    TString iTitel;

    while ((folderKey = (TKey*)folderNext())){
        // iterate through all voxels v

        folderName = folderKey->GetName();  // voxel name
        voxelData = (TDirectory*)f.Get(folderName);  // energy spectra for voxel v
        TIter spectrumNext(voxelData->GetListOfKeys());
        emissionsFromVoxel = 0;

        while ((spectrumKey = (TKey*)spectrumNext())){
            // iterate through all conditioned spectra
            // calculate all emissions from voxel v

            spectrumName = spectrumKey->GetName();
            iS_dc = (TH1*)voxelData->Get(spectrumName);
            emissionsInSpectrum = iS_dc->Integral();
            emissionsFromVoxel += emissionsInSpectrum;

            cout << spectrumName << "\n";
        }

        while ((spectrumKey = (TKey*)spectrumNext())){
            // iterate through all conditioned spectra
            // calculate all emissions from voxel v

            spectrumName = spectrumKey->GetName();
            iS_dc = (TH1*)voxelData->Get(spectrumName);
            emissionsInSpectrum = iS_dc->Integral();
            emissionsFromVoxel += emissionsInSpectrum;

            cout << spectrumName << "h\n";
        }

//        while ((spectrumKey = (TKey*)spectrumNext())){
//            // iterate through all conditioned spectra
//            // calculate all emissions from voxel v

//            cout << "Hi \n";
//        }


//        while ((spectrumKey = (TKey*)spectrumNext())){
//            // iterate through all conditioned spectra again
//            // calculate all emissions from voxel v

//////            spectrumName = spectrumKey->GetName();
//////            iS_dc = (TH1*)voxelData->Get(spectrumName);
//////            iS_dc->Scale(1.0 / emissionsFromVoxel);

////            cout << spectrumName << "\n";
//            cout << "Hi";

    }

    delete p_dcbv;
    delete folderKey;
    delete spectrumKey;
    delete iS_dc;

//    TH3F* emptyA_vArray = new TH3F("A_v", "Activity distribution",
//                                    iImageSpaceX + 1, 0, iImageSpaceX,    // x axis
//                                    iImageSpaceY + 1, 0, iImageSpaceY,    // y axis
//                                    iImageSpaceZ + 1, 0, iImageSpaceZ);   // z axis

}

void MLEM_reco(TString pathToMeasurements = "Position_x1.5_y-1.5.root",
               TString pathToProjections = "SPCIBase49.root"){

    // prepare the matrix N_dcb according to the measurement data
    TFile measurementsFile(pathToMeasurements, "READ");
    TH3I* N_dcbArray = createEmptyN_dcbArray(measurementsFile);
    fillN_dcbArray(measurementsFile, *N_dcbArray);
    cout << N_dcbArray->GetBinContent(6, 10, 500) << endl;
    measurementsFile.Close();

    // prepare the initial activity distribution
    TFile projectionsFile(pathToProjections, "READ");
    TH3F* A_vArray = createEmptyA_vArray(projectionsFile);
    setHomogeneousActivity(projectionsFile, *A_vArray, 1);  // fill every cell with 1
    cout << A_vArray->GetBinContent(5, 5, 0) << endl;

    // prepare the projection matrices p_dcbv
    createProjectionMatrices(projectionsFile);

    projectionsFile.Close();
}
