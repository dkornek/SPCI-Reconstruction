// ROOT Macro using the Origin Ensemble algorithm for SPCI-Reconstruction

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-05-27

/* NOMENCLATURE
* PAPER (Sitek, 2011)        MACRO                                 NAME
* i, i = 1...I               v, v = 1...V                          Voxel
* k, k = 1...K               dcb                                   Bin
* g_k                        N_dcb                                 Number of Events / Counts in Bin
* a_ki                       p_dcbv                                System Matrix Element (Probability)
* epsilon_i                  p_v = SUM_{dcb} (p_dcbv)              Sensitivity of Voxel
* f_i                        A_v                                   Estimate / Activity in Voxel
* s                          s                                     State
* c_si                       C_sv                                  Number of Event Origins in Voxel for State
* f_si = c_si / epsilon_i    A_sv = C_sv / p_v                     Estimate / Activity in Voxel fpr State
* rho_ski                    rho_sdcbv                             Probability of Event in Bin to have the Origin in Voxel
* rho_ski ~ a_ki f_si
*/

#include <iostream>
#include <vector>
#include <random>

#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2.h>
#include <TH3.h>
#include <TBenchmark.h>
#include <TCanvas.h>


// ##### STATE CHARACTERIZATION #####
class State{
public:
    std::vector<Int_t> event;                   // events n (1 ... N)
    std::vector<Int_t> voxelOrigin;             // voxel v containing the origin (1 ... V) of event n
    std::vector<std::array<Int_t, 3> > bin;     // assigned bin (d/c/bin) of event n

    std::vector<Int_t> countsInVoxel;           // counts C_sv in voxel v for all events
};

// ##### RECONSTRUCTION #####
class ReconstructionOE{
public:
    ReconstructionOE(TString pathToMeasurements, TString pathToProjections);
    ~ReconstructionOE();
    void start();

    // Activity distribution
    TH3F* A_v;
    TH2F* A_vProject3D;

private:
    // ##### UTILITIES #####
    void openFile(TString fileType, TString pathToFile);
    void getDetectorIndices(TString nameOfSpectrum, Int_t& d, Int_t& c);
    void getImageSpaceIndices(TString titleOfVoxel, Int_t& x, Int_t& y, Int_t& z);
    void fillImageSpaceIndices(Int_t x, Int_t y, Int_t z);
    void generateInitialState();
    void generateNextState();
    void checkValidity();
    void plotActivity();

    // ##### PREPARATION FUNCTIONS #####
    void extractInfoFromMeasurements();
    void createP_dcbv();
    void createA_v();

    // ##### CALCULATION FUNCTIONS #####
    void calculateActivity();
    Double_t calculateTransitionProbability();

    // ##### MEMBERS #####
    Bool_t isCalculationValid;

    TCanvas* canvasResults;
    TPad* plot3D;
    TPad* plot2D;

    // Activity Data
    Int_t A_vBinsX;
    Int_t A_vBinsY;
    Int_t A_vBinsZ;

    std::vector<Double_t> imageVolume;  // in mm, xMin, xMax, yMin, yMax, zMin, zMax
    std::vector<std::array<Int_t, 3> > imageIndices;
    Int_t numberOfVoxels;

    std::vector<State> beforeEquilibrium;  // states before in equilibrium
    std::vector<State> originEnsembles;  // states in equilibrium

    // Measurement Data
    TFile* measurementsFile;

    Int_t numberOfDetectors;
    Int_t NbinsMeasurements;
    Int_t numberOfEvents;

    TIter* nextS_dcMeasurements;
    TKey* keyS_dcMeasurements;

    // Projection Data / System Matrix
    TFile* projectionsFile;

    TList* p_dcbv;
    std::vector<Double_t> p_v;  // p_v = SUM_{dcb} (p_dcbv)

    Int_t NbinsProjections;

    TIter* nextVoxel;
    TKey* keyVoxel;
    TKey* keyS_dcVoxel;
};

// ##### PUBLIC FUNCTIONS #####
ReconstructionOE::ReconstructionOE(TString pathToMeasurements, TString pathToProjections) :
    canvasResults(nullptr), plot3D(nullptr), plot2D(nullptr){
    // prepare all data to start the reconstruction

    // open files
    openFile("MEASUREMENT", pathToMeasurements);
    nextS_dcMeasurements = new TIter(measurementsFile->GetListOfKeys());

    openFile("PROJECTION", pathToProjections);
    nextVoxel = new TIter(projectionsFile->GetListOfKeys());

    TBenchmark b;

    // prepare the measurement data
    extractInfoFromMeasurements();

    // prepare the system matrix
    b.Start("p_dcb");
    createP_dcbv();
    b.Stop("p_dcb");

    std::cout << "\np_dcb Creation Time:\t" << b.GetRealTime("p_dcb") << " s\n";
}

ReconstructionOE::~ReconstructionOE(){

    delete plot2D;
    delete plot3D;
    delete canvasResults;

    delete A_v;
    delete A_vProject3D;

    delete keyS_dcMeasurements;
    delete nextS_dcMeasurements;
    delete measurementsFile;

    delete keyS_dcVoxel;
    delete keyVoxel;
    delete nextVoxel;
    p_dcbv->Delete();
    delete p_dcbv;
    delete projectionsFile;
}

void ReconstructionOE::start(){
    // Reconstruct the image using the Origin Ensemble algorithm

    TBenchmark b;

    // Step 0: prepare the image space containing the activity
    b.Start("A_v");
    createA_v();
    b.Stop("A_v");
    std::cout << "\nA_v Creation Time:\t" << b.GetRealTime("A_v") << " s\n";

    // Step 1: create initial state s_0 by randomly selecting possible origins of the detected events
    b.Start("S_0");
    generateInitialState();
    b.Stop("S_0");
    std::cout << "\nS_0 Creation Time:\t" << b.GetRealTime("S_0") << " s\n";

    // Step 2: generate new states to reach equilibrium
    b.Start("Prep");
    generateNextState();
    b.Stop("Prep");
    std::cout << "\nTime to reach Equilibrium:\t" << b.GetRealTime("Prep") << " s\n";

//    // calculate final result
//    calculateActivity();

//    // Show results
//    plotActivity();
}

// ##### PRIVATE FUNCTIONS #####
// ##### UTILITIES #####
void ReconstructionOE::openFile(TString fileType, TString pathToFile){
    // open file of specified type

    if (fileType == "MEASUREMENT"){
        measurementsFile = new TFile(pathToFile, "READ");
        return;
    }

    if (fileType == "PROJECTION"){
        projectionsFile = new TFile(pathToFile, "READ");
        return;
    }
}

void ReconstructionOE::getDetectorIndices(TString nameOfSpectrum, Int_t &d, Int_t &c){
    // extract index of both the detectors from the name of the spectrum

    d = ((TString)nameOfSpectrum(0, 2)).Atoi();
    c = ((TString)nameOfSpectrum(2, 2)).Atoi();
}

void ReconstructionOE::getImageSpaceIndices(TString titleOfVoxel, Int_t &x, Int_t &y, Int_t &z){
    // extract location in image space from the title (=name of directory) of the voxel v

    // prepare title
    TString title = titleOfVoxel.Copy();
    title.Remove(0, 1);
    title.Replace(title.Length() - 1, 1, ',');

    // actual coordinates
    std::vector<Int_t> location = { 0, 0, 0 };
    for (Int_t i = 0; i < 3; ++i){
        Int_t positionOfSeparator = title.First(',');
        Int_t lengthOfString = title.Last(',');

        location[i] =  ((TString)title(0, positionOfSeparator)).Atoi();
        title = title(positionOfSeparator + 1, lengthOfString - positionOfSeparator + 1);
    }

    x = location[0];
    y = location[1];
    z = location[2];
}

void ReconstructionOE::fillImageSpaceIndices(Int_t x, Int_t y, Int_t z){
    // get a vector of the image space indices

    for (Int_t zIndex = 0; zIndex <= z; ++zIndex){
        for (Int_t yIndex = 0; yIndex <= y; ++yIndex){
            for (Int_t xIndex = 0; xIndex <= x; ++xIndex){

                std::array<Int_t, 3> coordinate = { xIndex + 1,
                                                    yIndex + 1,
                                                    zIndex + 1};
                imageIndices.push_back(coordinate);
            }
        }
    }
}

void ReconstructionOE::generateInitialState(){
    // generate random inital state for the origin locations in the image space

    std::random_device rd;  // obtain a seed for the random number engine
    std::mt19937 generate(rd());  // Mersenne-Twister random number engine
    std::uniform_real_distribution<Double_t> distribution(0, numberOfVoxels);

    State initialState;
    Int_t events = 0;

    nextS_dcMeasurements->Reset();
    while ((keyS_dcMeasurements = (TKey*)nextS_dcMeasurements->Next())){
        // iterate through every spectrum, retrieve every event, assign random location

        TString nameOfS_dc = keyS_dcMeasurements->GetName();
        TH1F* S_dc = (TH1F*)measurementsFile->Get(nameOfS_dc);

        Int_t detectorD;
        Int_t detectorC;
        getDetectorIndices(nameOfS_dc, detectorD, detectorC);

        for (Int_t bin = 1; bin <= NbinsMeasurements; ++bin){
            Int_t EventsInBin = S_dc->GetBinContent(bin);

            for (Int_t e = 0; e < EventsInBin; ++e){
                // save number of event n
                initialState.event.push_back(events);

                // save bin of event n
                std::array<Int_t, 3> binOfEvent;
                binOfEvent[0] = detectorD + 1;
                binOfEvent[1] = detectorC + 1;
                binOfEvent[2] = bin;
                initialState.bin.push_back(binOfEvent);

                // generate random origin (voxel) for event n
                Int_t origin;
                while (true){
                    origin = Int_t(distribution(generate));

                    // check, if probability for that origin is > 0
                    TH3F* p_dcb = (TH3F*)p_dcbv->At(origin);
                    Double_t probability = p_dcb->GetBinContent(detectorD + 1, detectorC + 1, bin);

                    if (probability > 0){
                        break;
                    }
                }

                initialState.voxelOrigin.push_back(origin);
                ++events;
            }
        }

        delete S_dc;
    }
    beforeEquilibrium.push_back(initialState);

    // save number of total events
    numberOfEvents = initialState.event.size();
}

void ReconstructionOE::generateNextState(){
    // generate a new state

    std::random_device rd;  // obtain a seed for the random number engine
    std::mt19937 generate(rd());  // Mersenne-Twister random number engine
    std::uniform_real_distribution<Double_t> distributionEvents(0, numberOfEvents);
    std::uniform_real_distribution<Double_t> distributionVoxels(0, numberOfVoxels);
    std::uniform_real_distribution<Double_t> transition(0, 1);

    for (Int_t n = 0; n < numberOfEvents; ++n){

        // randomly select event n and voxel v of event n
        Int_t randomEvent = Int_t(distributionEvents(generate));
        Int_t voxelOfEventOld = beforeEquilibrium.at(0).voxelOrigin.at(randomEvent);

        // randomly select a new voxel v' for event n
        Int_t voxelOfEventNew = Int_t(distributionVoxels(generate));

        // if they are the same, then don't do anything
        if (voxelOfEventNew == voxelOfEventOld){
            continue;
        }

        // calculate the transition probability
        TH3F* p_dcbOld = (TH3F*)p_dcbv->At(voxelOfEventOld);
        Double_t systemMatrixElementOld;
        Double_t sensitivityOld = p_v.at(voxelOfEventOld);
        Double_t countsInVoxelOld;

        TH3F* p_dcbNew = (TH3F*)p_dcbv->At(voxelOfEventNew);
        Double_t systemMatrixElementNew;
        Double_t sensitivityNew = p_v.at(voxelOfEventNew);
        Double_t countsInVoxelNew;
    }
}

void ReconstructionOE::checkValidity(){
    // calculation is only possible if binning pattern in the measurement and
    // projection spectra is identical

    if (NbinsMeasurements == NbinsProjections){
        isCalculationValid = kTRUE;

    } else{
        std::cout << "Calculation is not valid.\nBinning pattern is not consistent\n";
        isCalculationValid = kFALSE;
    }
}

void ReconstructionOE::plotActivity(){
    // show results of the activity

    if (canvasResults == nullptr){
        canvasResults = new TCanvas("results", "Image Reconstruction", 10, 10, 900, 450);
    }
    canvasResults->cd();

    if (plot3D == nullptr){
        plot3D = new TPad("activity3D", "3D Activity Distribution", 0.01, 0.01, 0.50, 0.99);
        plot3D->Draw();
    }

    if (plot2D == nullptr){
        plot2D = new TPad("activity2d", "2D Projection", 0.51, 0.01, 0.99, 0.99);
        plot2D->SetRightMargin(0.2);
        plot2D->Draw();
    }

    plot3D->cd();
    A_v->Draw("BOX2");

    plot2D->cd();
    A_vProject3D = (TH2F*)A_v->Project3D("yx");
    A_vProject3D->SetStats(kFALSE);
    A_vProject3D->SetTitle("2D Projection of 3D Activity Distribution");
    A_vProject3D->GetXaxis()->SetTitleOffset(1);
    A_vProject3D->GetYaxis()->SetTitleOffset(1);

    A_vProject3D->GetZaxis()->SetTitle("Relative Activity");
    A_vProject3D->GetZaxis()->SetTitleOffset(1.8);
    A_vProject3D->Draw("COLZ");

    canvasResults->Update();
}

// ##### PREPARATION FUNCTIONS #####
void ReconstructionOE::extractInfoFromMeasurements(){
    // extract number of detectors, number of bins in spectra

    // get number of detectors
    TString nameOfLastSpectrum;
    nameOfLastSpectrum = measurementsFile->GetListOfKeys()->Last()->GetName();

    Int_t detectorD;
    Int_t detectorC;
    getDetectorIndices(nameOfLastSpectrum, detectorD, detectorC);

    numberOfDetectors = std::max(detectorD, detectorC) + 1;

    // get number of bins
    TH1F* lastSpectrum = (TH1F*)measurementsFile->Get(nameOfLastSpectrum);
    NbinsMeasurements = lastSpectrum->GetNbinsX();
    delete lastSpectrum;
}

void ReconstructionOE::createP_dcbv(){
    // create list of 3d histograms containing the probabilities for each voxel

    p_dcbv = new TList();

    // get number of bins
    TString nameOfLastVoxel = projectionsFile->GetListOfKeys()->Last()->GetName();
    TDirectory* dirOfLastVoxel = (TDirectory*)projectionsFile->Get(nameOfLastVoxel);
    TString nameOfLastS_dc = dirOfLastVoxel->GetListOfKeys()->Last()->GetName();
    TH1F* lastS_dc = (TH1F*)dirOfLastVoxel->Get(nameOfLastS_dc);
    NbinsProjections = lastS_dc->GetNbinsX();
    delete lastS_dc;
    delete dirOfLastVoxel;

    checkValidity();
    if (!isCalculationValid){
        return;
    }

    // iterate through all voxels v
    nextVoxel->Reset();
    while ((keyVoxel = (TKey*)nextVoxel->Next())){
        TString nameOfVoxel = keyVoxel->GetName();

        TString description;
        description.Form("Probabilities to contribute one count to N_dcb from voxel %s", nameOfVoxel.Data());

        TString internalName;
        internalName.Form("p_dcb_%s", nameOfVoxel.Data());
        TH3F* p_dcb = new TH3F(internalName, description,
                               numberOfDetectors, 0, numberOfDetectors,
                               numberOfDetectors, 0, numberOfDetectors,
                               NbinsProjections, 0, NbinsProjections);

        TDirectory* dirOfVoxel = (TDirectory*)projectionsFile->Get(nameOfVoxel);
        TIter nextS_dc(dirOfVoxel->GetListOfKeys());
        while ((keyS_dcVoxel = (TKey*)nextS_dc())){
            // iterate through all spectra in voxel v

            TString nameOfS_dc = keyS_dcVoxel->GetName();
            TH1F* S_dc = (TH1F*)dirOfVoxel->Get(nameOfS_dc);
            // S_dc->Scale(1.0 / S_dc->Integral());  // faster convergence

            Int_t detectorD;
            Int_t detectorC;
            getDetectorIndices(nameOfS_dc(3, 4), detectorD, detectorC);
            for (Int_t bin = 1; bin <= NbinsProjections; ++bin){
                // set probability

                p_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                     S_dc->GetBinContent(bin));
            }

            delete S_dc;
        }

        p_dcbv->AddLast(p_dcb);
        p_v.push_back(p_dcb->Integral());

        delete dirOfVoxel;
    }
}

void ReconstructionOE::createA_v(){
    // create the activity distribution A (=A_v)
    // all voxels of the image space contain an activity A_v

    // the image space is based on the last folder title in the projections
    // title is expected to have the form "(X,Y,Z)"

    // get volume
    TString titleOfLastVoxel;
    titleOfLastVoxel = projectionsFile->GetListOfKeys()->Last()->GetTitle();
    getImageSpaceIndices(titleOfLastVoxel, A_vBinsX, A_vBinsY, A_vBinsZ);

    fillImageSpaceIndices(A_vBinsX, A_vBinsY, A_vBinsZ);
    numberOfVoxels = imageIndices.size();

    if (imageVolume.size() != 6){
        imageVolume = {0, static_cast<Double_t>(A_vBinsX + 1),
                       0, static_cast<Double_t>(A_vBinsY + 1),
                       0, static_cast<Double_t>(A_vBinsZ + 1)};
    }

    // create the 3d histogram
    A_v = new TH3F("A_v", "3D Activity Distribution",
                   A_vBinsX + 1, imageVolume[0], imageVolume[1],
                   A_vBinsY + 1, imageVolume[2], imageVolume[3],
                   A_vBinsZ + 1, imageVolume[4], imageVolume[5]);

    // Plotting options
    A_v->GetXaxis()->SetTitle("#it{x} / mm");
    A_v->GetXaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * A_vBinsX), "X");

    A_v->GetYaxis()->SetTitle("#it{y} / mm");
    A_v->GetXaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * A_vBinsY), "Y");

    A_v->GetZaxis()->SetTitle("#it{z} / mm");
    A_v->GetZaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * A_vBinsZ), "Z");
}

// ##### CALCULATION FUNCTIONS #####
void ReconstructionOE::calculateActivity(){
    // calculate the mean of the origin ensemble to get the expectation value of the image

    State currentState = beforeEquilibrium.at(0);

    for (UInt_t n = 0; n < currentState.event.size(); ++n){
        Int_t voxel = currentState.voxelOrigin.at(n);

        std::array<Int_t, 3> coordinates = imageIndices.at(voxel);

        Int_t binContent = A_v->GetBinContent(coordinates[0], coordinates[1], coordinates[2]) + 1;
        A_v->SetBinContent(coordinates[0], coordinates[1], coordinates[2], binContent);
    }
}

Double_t ReconstructionOE::calculateTransitionProbability(){
    // calculate the transition probability to move the origin to another voxel

    Double_t probability;

    return probability;
}

// ##### ROOT #####
void OE(){
    TBenchmark b;
    b.Start("total");

    // Specify the location of the measurement file
    TString pathToMeasurements = "../folder/subfolder/*.root";
    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0.root";

    // Specify the location of the projections file
    TString pathToProjection = "../folder/subfolder/*.root";
    pathToProjection = "../../data/SystemMatrix/Bins50/SPCIBase49.root";

    // create an instance of the reconstruction class
    ReconstructionOE* reco = new ReconstructionOE(pathToMeasurements, pathToProjection);
    reco->start();

    b.Stop("total");
    std::cout << "\nTotal Time:\t\t" << b.GetRealTime("total") << " seconds\n";
}

// ##### COMPILE FILE #####
#ifndef __CINT__
int main(){
    OE();

    return 0;
}
#endif
