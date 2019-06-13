// originensemble.h

#pragma once
#include <vector>
#include <TH2.h>
#include <TH3.h>
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "imagespace.h"
#include "systemmatrix.h"
#include "measurements.h"
#include "state.h"
#include "utilities.h"


const Bool_t normalizeSpectra = kFALSE;

// ##### RESULTS #####
class ResultsOE{
public:
    ResultsOE();
    ~ResultsOE();

    void plot(TH3F* A_v);

private:
    TH2F* A_vProject3D = nullptr;

    TCanvas* canvas;
    TPad* plot3D;
    TPad* plot2D;
};

ResultsOE::ResultsOE(){
    canvas = new TCanvas("results", "Image Reconstruction", 10, 10, 900, 450);
    canvas->cd();

    plot3D = new TPad("activity3D", "3D Activity Distribution", 0.01, 0.01, 0.50, 0.99);
    plot3D->Draw();

    plot2D = new TPad("activity2d", "2D Projection", 0.51, 0.01, 0.99, 0.99);
    plot2D->SetRightMargin(0.2);
    plot2D->Draw();
}

ResultsOE::~ResultsOE(){
    delete A_vProject3D;
    delete plot2D;
    delete plot3D;
    delete canvas;
}

void ResultsOE::plot(TH3F *A_v){
    plot3D->cd();
    A_v->Draw("BOX2");

    plot2D->cd();
    A_vProject3D = (TH2F*)A_v->Project3D("yx");
    A_vProject3D->SetStats(kTRUE);
    A_vProject3D->SetTitle("2D Projection of 3D Activity Distribution");

    A_vProject3D->GetXaxis()->SetTitleOffset(1);
    A_vProject3D->GetYaxis()->SetTitleOffset(1);
    A_vProject3D->GetZaxis()->SetTitle("Relative Activity");
    A_vProject3D->GetZaxis()->SetTitleOffset(1.8);
    A_vProject3D->Draw("COLZ");
}

// ##### RECONSTRUCTION #####
class ReconstructionOE{
public:
    ReconstructionOE(const TString pathToMeasurements,
                     const TString pathToProjections,
                     const std::vector<Double_t> volume);
    ~ReconstructionOE();

    void start(const Int_t numberOfIterationsForEquilibirum, const Int_t numberOfIterationsInEquilibrium);

private:
    // ##### CALCULATION FUNCTIONS #####
    void reachEquilibrium(const Int_t numberOfIterations);
    void sampleEquilibriumStates(const Int_t numberOfIterations);
    void calculateActivity();

    // ##### MEMBERS #####
    Bool_t isCalculationValid;

    std::vector<State> statesBeforeEquilibrium;
    std::vector<State> statesInEquilibrium;

    Measurements* measurementData = nullptr;
    SystemMatrix* systemMatrixData = nullptr;
    ImageSpace* image = nullptr;
    ResultsOE* results = nullptr;
};

ReconstructionOE::ReconstructionOE(const TString pathToMeasurements,
                                   const TString pathToProjections,
                                   const std::vector<Double_t> volume){
    // prepare data for reconstruction using OE

    TBenchmark b;

    // prepare the measurement data
    b.Start("N_dcb");
    measurementData = new Measurements(pathToMeasurements);
    measurementData->createN_dcb();
    measurementData->fillN_dcb(normalizeSpectra);
    b.Stop("N_dcb");

    // prepare the image & system matrix
    b.Start("p_dcb");
    systemMatrixData = new SystemMatrix(pathToProjections);

    isCalculationValid = Utilities::checkForSameNumberOfBins(measurementData->numberOfBins,
                                                             systemMatrixData->numberOfBins);

    if (isCalculationValid){
        systemMatrixData->createSystemMatrix(normalizeSpectra,
                                             measurementData->numberOfDetectors);
    }
    b.Stop("p_dcb");

    b.Start("A_v");
    image = new ImageSpace(volume, systemMatrixData->size);
    b.Stop("A_v");

    std::cout << "\nN_dcb Creation Time:\t" << b.GetRealTime("N_dcb") << " s\n";
    std::cout << "\np_dcb Creation Time:\t" << b.GetRealTime("p_dcb") << " s\n";
    std::cout << "\nA_v Creation Time:\t" << b.GetRealTime("A_v") << " s\n";
}

ReconstructionOE::~ReconstructionOE(){
    delete results;
    delete image;
    delete systemMatrixData;
    delete measurementData;
}

void ReconstructionOE::start(const Int_t numberOfIterationsForEquilibirum, const Int_t numberOfIterationsInEquilibrium){
    // start the OE image reconstruction

    TBenchmark b;

    results = new ResultsOE();

    // step 1: create initial state s_0 by randomly selecting possible origins for the detected events
    b.Start("stats");
    b.Start("S_0");
    State initalState(measurementData->N_dcb);
    initalState.generateRandomOrigins(image->numberOfVoxels, systemMatrixData->systemMatrix);
    statesBeforeEquilibrium.push_back(initalState);
    b.Stop("S_0");
    std::cout << "\nS_0 Creation Time:\t" << b.GetRealTime("S_0") << " s\n";

    // step 2: generate new states until equilibrium is reached
    b.Start("UntilE");
    reachEquilibrium(numberOfIterationsForEquilibirum);
    b.Stop("UntilE");
    std::cout << "\nReach Equilibrium Time:\t" << b.GetRealTime("UntilE") << " s\n";

    // step 3: generate new states in equilibrium
    b.Start("InE");
    sampleEquilibriumStates(numberOfIterationsInEquilibrium);
    b.Stop("InE");
    std::cout << "\nSampling Time:\t" << b.GetRealTime("InE") << " s\n";

    // step 4: calculate "mean state" = mean of counts in each voxel of all sampled states
    calculateActivity();
    b.Stop("stats");

    std::cout << "\nImage reconstruction done.\n";
    std::cout << "\nCalculation Time:\t" << b.GetRealTime("stats") << " seconds\n";

    // step 5: show results
    results->plot(image->A_v);
}

// ##### CALCULATION FUNCTIONS #####
void ReconstructionOE::reachEquilibrium(const Int_t numberOfIterations){
    // generate random states until equilibrium is reached

    for (Int_t n = 0; n < numberOfIterations; ++n){
        State nextState = statesBeforeEquilibrium.back();
        nextState.MCMCNextState(systemMatrixData->systemMatrix, systemMatrixData->sensitivities);
        statesBeforeEquilibrium.push_back(nextState);
    }
}

void ReconstructionOE::sampleEquilibriumStates(const Int_t numberOfIterations){
    // generate states in equilibrium

    State equilibriumState = statesBeforeEquilibrium.back();
    statesInEquilibrium.push_back(equilibriumState);

    for (Int_t n = 0; n < numberOfIterations; ++n){
        State nextState = statesInEquilibrium.back();
        nextState.MCMCNextState(systemMatrixData->systemMatrix, systemMatrixData->sensitivities);

        if (n % 1 == 0){
            // save every state, can be changed according to needs
            statesInEquilibrium.push_back(nextState);
        }
    }
}

void ReconstructionOE::calculateActivity(){
    // calculate the voxel-specific means of the detected counts

    Int_t numberOfStates = statesInEquilibrium.size();
    for (Int_t v = 0; v < image->numberOfVoxels; ++v){

        // calculate the mean counts
        Double_t meanCountInVoxel = 0.0;
        for (Int_t c = 0; c < numberOfStates; ++c){
            meanCountInVoxel += statesInEquilibrium.at(c).countsInVoxel.at(v);
        }
        meanCountInVoxel = meanCountInVoxel / numberOfStates;

        // calculate the mean activity
        Double_t sensitivityOfVoxel = systemMatrixData->sensitivities.at(v);
        Double_t activityInVoxel = meanCountInVoxel / sensitivityOfVoxel;

        // assign mean activity to voxel in image space
        std::array<Int_t, 3> coordinate = image->imageIndices.at(v);
        image->A_v->SetBinContent(coordinate[0], coordinate[1], coordinate[2],
                                  activityInVoxel);
    }

    image->A_v->Scale(1.0);  // bug handling to make Project3D work
}
