/*
 * TO DO
 * Remove Root Part
*/

#include <iostream>
#include <chrono>

#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TH3.h>


class ReconstructionMLEM{
public:
    ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections);
    ~ReconstructionMLEM();
    void start(Int_t numberOfIterations);

    // Activity distribution
    TH3F* A_v;

private:
    // ##### UTILITIES #####
    void openFile(TString fileType, TString pathToFile);
    void getDetectorIndices(TString nameOfSpectrum, Int_t& d, Int_t& c);
    void getImageSpaceIndices(TString titleOfVoxel, Int_t& x, Int_t& y, Int_t& z);
    void fillN_dcb();
    void makeA_vHomogeneous();
    void checkValidity();

    // ##### PREPARATION FUNCTIONS #####
    void createN_dcb();
    void createP_dcbv();
    void createA_v();

    // ##### CALCULATION FUNCTIONS #####
    void calculate();

    // ##### MEMBERS #####

    Bool_t isCalculationValid;

    // Measurement Data
    TFile* measurementsFile;
    TH3F* N_dcb;

    Int_t numberOfDetectors;
    Int_t NbinsMeasurements;

    TIter* nextS_dcMeasurements;
    TKey* keyS_dcMeasurements;

    // Projection Data
    // Probabilities
    TFile* projectionsFile;
    TList* p_dcbv;

    Int_t NbinsProjections;

    TIter* nextVoxel;
    TKey* keyVoxel;
    TKey* keyS_dcVoxel;
};

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
    delete this->nextS_dcMeasurements;
    delete this->N_dcb;
    delete this->measurementsFile;

    delete this->nextVoxel;
    p_dcbv->Delete();
    delete this->p_dcbv;
    delete this->projectionsFile;
}

void ReconstructionMLEM::start(Int_t numberOfIterations){
    // execute the calculation

    this->checkValidity();

    if (this->isCalculationValid){
        for (Int_t i = 0; i < numberOfIterations; ++i){
            this->calculate();
        }
    }
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
        for (Int_t bin = 0; bin <= this->NbinsMeasurements; ++bin){
            this->N_dcb->SetBinContent(detectorD, detectorC, bin,
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
        this->A_v->SetBinContent(indexX, indexY, indexZ, 1.0);
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

    // iterate through every voxel v
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
            for (Int_t bin = 0; bin <= this->NbinsProjections; ++bin){
                p_dcb->SetBinContent(detectorD, detectorC, bin,
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
    // to calculate the activiy distribution A (=A_v)

    TH3F currentActivity;
    this->A_v->Copy(currentActivity);
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
    TH3F* numerator;
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
        activityInVoxel = currentActivity.GetBinContent(indexX, indexY, indexZ);

        // get p_dcb for voxel v
        p_dcb = (TH3F*)this->p_dcbv->At(nameOfVoxel.Atoi());

        // multiply them and add them to the denominator
        denominator->Add(p_dcb, activityInVoxel);
    }

    // calculate new activity distribution
    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){

        // get activity A_v in voxel v
        nameOfVoxel = this->keyVoxel->GetName();
        titleOfVoxel = this->keyVoxel->GetTitle();

        this->getImageSpaceIndices(titleOfVoxel, indexX, indexY, indexZ);
        activityInVoxel = currentActivity.GetBinContent(indexX, indexY, indexZ);

        // calculate correction factor
        numerator = (TH3F*)this->p_dcbv->At(nameOfVoxel.Atoi());
        numerator->Multiply(this->N_dcb);
        numerator->Divide(denominator);
        correctionFactor = numerator->Integral();

        // calculate new activity A_v
        this->A_v->SetBinContent(indexX, indexY, indexZ,
                                 activityInVoxel * correctionFactor);
    }

    delete denominator;
}

// ##### ROOT #####
void mlem(){
    auto t1 = std::chrono::steady_clock::now();

    // create an instance of the reconstruction class
    TString pathM = "../data/measurement_data/bins_250/20181115_137Cs_notLinear7_sim_1_250_bins.root";
    TString pathP = "../data/projection_data/bins_250/SPCIBase49_250_bins.root";
    ReconstructionMLEM* reco = new ReconstructionMLEM(pathM, pathP);
    reco->start(10);

    TCanvas* canvas = new TCanvas("c", "Activity distribution", 600, 400);
    //TH2D* image = (TH2D*)reco->A_v->Project3D("xy o");  // not suitable
    TH2F* image = new TH2F("adi", "Reconstruction",
                           7, 0, 7,
                           7, 0, 7);

    Double_t binContent;
    Int_t voxel = 0;
    for (Int_t d = 0; d < 7; ++d){
        for (Int_t c = 0; c < 7; ++c){
            binContent = reco->A_v->GetBinContent(d, c, 0);

            if (voxel == 48){
                std::cout << voxel << ": " << binContent << "\n";
            }
            image->SetBinContent(d+1, c+1, binContent);

            ++voxel;
        }
    }

    std::cout << reco->A_v->GetNbinsX() <<std::endl;
    std::cout << reco->A_v->GetNbinsY() <<std::endl;
    std::cout << reco->A_v->GetNbinsZ() <<std::endl;

    reco->A_v->Draw("BOX");

    std::cout << image->GetBinContent(6, 6, 0);

    //gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(kBird);
    //image->Draw("COLZ");
    canvas->Update();

//    delete reco;

    auto t2 = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    std::cout << "\nIt took: " << elapsedTime << " seconds\n";
}

#ifndef __CINT__
int main(){
    mlem();
    return 0;
}
#endif