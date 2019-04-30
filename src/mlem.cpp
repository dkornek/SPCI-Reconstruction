// mlem.cpp

#include "mlem.h"

#include <iostream>
#include <chrono>

// ##### PUBLIC FUNCTIONS #####
ReconstructionMLEM::ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections){
    // prepare all data to start the reconstruction

    // open files
    this->openFile("MEASUREMENT", pathToMeasurements);
    this->nextS_dcMeasurements = new TIter(this->measurementsFile->GetListOfKeys());

    this->openFile("PROJECTION", pathToProjections);
    this->nextVoxel = new TIter(this->projectionsFile->GetListOfKeys());

    // prepare the measurement data N_dcb
    // N_dcb = number of events in detector pair (d, c) separated by bins
    this->createN_dcb();

    // prepare the probabilities p_dcbv
    // p_dcbv = voxel v specific probabilitiy to contribute one count to N_dcb;
    // there are v 3d histograms p_dcb in the list p_dcbv
    this->createP_dcbv();

    // prepare the homogeneous activity distribution A_v
    // A_v = activity in voxel v of the image space
    this->createA_v();
}

ReconstructionMLEM::~ReconstructionMLEM(){
    delete this->A_v;

    delete this->nextS_dcMeasurements;
    delete this->N_dcb;
    delete this->measurementsFile;

    delete this->nextVoxel;
    p_dcbv->Delete();
    delete this->p_dcbv;
    delete this->projectionsFile;
}

void ReconstructionMLEM::start(Int_t maxNumberOfIterations, Double_t stopCriterion){
    // execute the calculation

    this->checkValidity();
    if (!this->isCalculationValid){
        return;
    }

    Int_t numberOfIterations = 0;
    auto t1 = std::chrono::steady_clock::now();

    for (;;++numberOfIterations){
        this->calculate();

        if (numberOfIterations == maxNumberOfIterations){
            break;
        }

        if ((this->deviation >= stopCriterion) && (this->deviation <= 1.0)){
            break;
        }
    }

    auto t2 = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

    std::cout << "Image reconstruction done. Steps: " << numberOfIterations << "\n";
    std::cout << "\nCalculation time: " << elapsedTime << " seconds\n";
}

// ##### PRIVATE FUNCTIONS #####
// ##### UTILITIES #####
void ReconstructionMLEM::openFile(TString fileType, TString pathToFile){
    // open file of specified type

    if (fileType == "MEASUREMENT"){
        this->measurementsFile = new TFile(pathToFile, "READ");
        return;
    }

    if (fileType == "PROJECTION"){
        this->projectionsFile = new TFile(pathToFile, "READ");
        return;
    }
}

void ReconstructionMLEM::getDetectorIndices(TString nameOfSpectrum, Int_t &d, Int_t &c){
    // extract index of both the detectors from the name of the spectrum

    d = ((TString)nameOfSpectrum(0, 2)).Atoi();
    c = ((TString)nameOfSpectrum(2, 4)).Atoi();
}

void ReconstructionMLEM::getImageSpaceIndices(TString titleOfVoxel, Int_t &x, Int_t &y, Int_t &z){
    // extraxct location in image space from the title (=name of directory) of the voxel v

    x = ((TString)titleOfVoxel(1)).Atoi();
    y = ((TString)titleOfVoxel(3)).Atoi();
    z = ((TString)titleOfVoxel(5)).Atoi();
}

void ReconstructionMLEM::fillN_dcb(){
    // fill N_dcb with number of events extracted from the spectra

    Int_t detectorD;
    Int_t detectorC;

    TString nameOfS_dc;
    TH1F* S_dc;

    this->nextS_dcMeasurements->Reset();
    while ((this->keyS_dcMeasurements = (TKey*)this->nextS_dcMeasurements->Next())){
        nameOfS_dc = this->keyS_dcMeasurements->GetName();
        // get detector indices
        this->getDetectorIndices(nameOfS_dc, detectorD, detectorC);

        // get spectrum
        S_dc = (TH1F*)this->measurementsFile->Get(nameOfS_dc);

        // fill N_dcb
        for (Int_t bin = 1; bin <= this->NbinsMeasurements; ++bin){
            this->N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                       S_dc->GetBinContent(bin));
        }
    }
}

void ReconstructionMLEM::makeA_vHomogeneous(){
    // set all bin contents in A_v = 1

    TString titleOfVoxel;
    Int_t indexX;
    Int_t indexY;
    Int_t indexZ;

    // iterate through all voxels
    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){
        titleOfVoxel = this->keyVoxel->GetTitle();
        this->getImageSpaceIndices(titleOfVoxel, indexX, indexY, indexZ);
        this->A_v->SetBinContent(indexX + 1, indexY + 1, indexZ + 1,
                                 1.0);
    }
}

void ReconstructionMLEM::checkValidity(){
    // calculation is only possible if binning pattern in the measurement and
    // projection spectra is identical

    if (this->NbinsMeasurements == this->NbinsProjections){
        this->isCalculationValid = kTRUE;
    } else{
        std::cout << "Calculation is not valid.\nBinning pattern is not consistent\n";
        this->isCalculationValid = kFALSE;
    }
}

// ##### PREPARATION FUNCTIONS #####
void ReconstructionMLEM::createN_dcb(){
    // create empty 3d histogram of N_dcb
    // x = detector d, y = detector c, z = bins in spectrum

    // histogram is based on the last spectrum
    // title of spectra is expected to be of "ddcc" format
    Int_t detectorD;
    Int_t detectorC;

    TString nameOfLastSpectrum;
    nameOfLastSpectrum = this->measurementsFile->GetListOfKeys()->Last()->GetName();
    this->getDetectorIndices(nameOfLastSpectrum, detectorD, detectorC);

    // get number of detectors
    this->numberOfDetectors = std::max(detectorD, detectorC) + 1;

    // get number of bins
    this->NbinsMeasurements = ((TH1F*)this->measurementsFile->Get(nameOfLastSpectrum))->GetNbinsX();

    // create the 3d histogram
    this->N_dcb = new TH3F("N_dcb", "Number of events in bin b in S_dc",
                           this->numberOfDetectors, 0, this->numberOfDetectors,
                           this->numberOfDetectors, 0, this->numberOfDetectors,
                           this->NbinsMeasurements, 0, this->NbinsMeasurements);
    this->fillN_dcb();
}

void ReconstructionMLEM::createP_dcbv(){
    // create list of 3d histograms containing the probabilities for each voxel

    this->p_dcbv = new TList();

    // information about voxel
    TString nameOfVoxel;
    TDirectory* dirOfVoxel;
    Double_t emissionsInVoxel;

    // information about p_dcb
    TString internalName;
    TString description;

    // information about energy spectrum
    TString nameOfS_dc;
    TH1F* S_dc;
    Int_t detectorD;
    Int_t detectorC;

    // get number of bins
    nameOfVoxel = this->projectionsFile->GetListOfKeys()->Last()->GetName();
    dirOfVoxel = (TDirectory*)this->projectionsFile->Get(nameOfVoxel);
    nameOfS_dc = dirOfVoxel->GetListOfKeys()->Last()->GetName();
    this->NbinsProjections = ((TH1F*)dirOfVoxel->Get(nameOfS_dc))->GetNbinsX();

    // iterate through all voxels v
    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){
        emissionsInVoxel = 0.0;

        nameOfVoxel = this->keyVoxel->GetName();
        internalName.Form("p_dcb_%s", nameOfVoxel.Data());
        description.Form("Probabilities to contribute one count to N_dcb from voxel %s", nameOfVoxel.Data());

        TH3F* p_dcb = new TH3F(internalName, description,
                               this->numberOfDetectors, 0, this->numberOfDetectors,
                               this->numberOfDetectors, 0, this->numberOfDetectors,
                               this->NbinsProjections, 0, this->NbinsProjections);

        dirOfVoxel = (TDirectory*)this->projectionsFile->Get(nameOfVoxel);

        // iterate through all spectra in voxel v
        TIter nextS_dc(dirOfVoxel->GetListOfKeys());
        while ((this->keyS_dcVoxel = (TKey*)nextS_dc())){
            nameOfS_dc = this->keyS_dcVoxel->GetName();
            S_dc = (TH1F*)dirOfVoxel->Get(nameOfS_dc);
            emissionsInVoxel += S_dc->Integral();

            this->getDetectorIndices(nameOfS_dc(3, 7), detectorD, detectorC);
            for (Int_t bin = 1; bin <= this->NbinsProjections; ++bin){
                p_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                     S_dc->GetBinContent(bin));
            }
        }

        p_dcb->Scale(1.0 / emissionsInVoxel);
        this->p_dcbv->AddLast(p_dcb);
    }
}

void ReconstructionMLEM::createA_v(){
    // create the activity distribution A (=A_v)
    // all voxels of the image space contain an activity A_v

    // the image space is based on the last folder title in the projections
    // title is expected to have the form "(X,Y,Z)"
    Int_t indexX;
    Int_t indexY;
    Int_t indexZ;

    // get volume
    TString titleOfLastVoxel;
    titleOfLastVoxel = this->projectionsFile->GetListOfKeys()->Last()->GetTitle();
    this->getImageSpaceIndices(titleOfLastVoxel, indexX, indexY, indexZ);

    // create the 3d histogram
    this->A_v = new TH3F("A_v", "Activity distribution",
                         indexX + 1, 0, indexX + 1,
                         indexY + 1, 0, indexY + 1,
                         indexZ + 1, 0, indexZ + 1);
    this->makeA_vHomogeneous();
}

// ##### CALCULATION FUNCTIONS #####
void ReconstructionMLEM::calculate(){
    // execute the maximum likelihood expectation maximization algorithm
    // to calculate the activity distribution A (=A_v)

    TH3F* currentActivity = (TH3F*)this->A_v->Clone("Current Activity");
    Double_t activityInVoxel;

    // information about voxel v
    TString titleOfVoxel;
    TString nameOfVoxel;
    Int_t indexX;
    Int_t indexY;
    Int_t indexZ;

    // information about p_dcbv
    TH3F* p_dcb;

    // correction factors
    Double_t correctionFactor;
    TH3F* quotient;
    TH3F* denominator = new TH3F("denom", "Denominator for correction factor",
                                 this->numberOfDetectors, 0, this->numberOfDetectors,
                                 this->numberOfDetectors, 0, this->numberOfDetectors,
                                 this->NbinsProjections, 0, this->NbinsProjections);

    // calculate the denominator
    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){

        // get activity A_v in voxel v
        nameOfVoxel = this->keyVoxel->GetName();
        titleOfVoxel = this->keyVoxel->GetTitle();

        this->getImageSpaceIndices(titleOfVoxel, indexX, indexY, indexZ);
        activityInVoxel = currentActivity->GetBinContent(indexX + 1, indexY + 1, indexZ + 1);

        // get p_dcb for voxel v
        p_dcb = (TH3F*)this->p_dcbv->At(nameOfVoxel.Atoi());

        // multiply them and add them to the denominator
        denominator->Add(p_dcb, activityInVoxel);
    }

    // prepare the fraction part of the correction factor
    this->deviation = this->N_dcb->Integral() / denominator->Integral();
    quotient = (TH3F*)this->N_dcb->Clone("Measurements");
    quotient->Divide(denominator);

    // calculate new activity distribution
    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){

        // get activity A_v in voxel v
        nameOfVoxel = this->keyVoxel->GetName();
        titleOfVoxel = this->keyVoxel->GetTitle();

        this->getImageSpaceIndices(titleOfVoxel, indexX, indexY, indexZ);
        activityInVoxel = currentActivity->GetBinContent(indexX + 1, indexY + 1, indexZ + 1);

        // calculate correction factor
        p_dcb = (TH3F*)this->p_dcbv->At(nameOfVoxel.Atoi());
        p_dcb->Multiply(quotient);
        correctionFactor = p_dcb->Integral();

        // calculate new activity A_v
        this->A_v->SetBinContent(indexX + 1, indexY + 1, indexZ + 1,
                                 activityInVoxel * correctionFactor);
    }

    delete quotient;
    delete denominator;
    delete currentActivity;
}
