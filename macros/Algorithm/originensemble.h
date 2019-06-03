// originensemble.h

#pragma once
#include "imagespace.h"
#include "systemmatrix.h"
#include "measurements.h"
#include "state.h"
#include "utilities.h"


const Bool_t normalizeSpectra = kFALSE;

// ##### RESULTS #####
class ResultsOE{

};

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

    // ##### MEMBERS #####
    Bool_t isCalculationValid;

    std::vector<State> statesBeforeEquilibrium;

    Measurements* measurementData = nullptr;
    SystemMatrix* systemMatrixData = nullptr;
    ImageSpace* image = nullptr;
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
    image = new ImageSpace(volume,
                           systemMatrixData->size);
    b.Stop("A_v");

    std::cout << "\nN_dcb Creation Time:\t" << b.GetRealTime("N_dcb") << " s\n";
    std::cout << "\np_dcb Creation Time:\t" << b.GetRealTime("p_dcb") << " s\n";
    std::cout << "\nA_v Creation Time:\t" << b.GetRealTime("A_v") << " s\n";
}

ReconstructionOE::~ReconstructionOE(){
    delete image;
    delete systemMatrixData;
    delete measurementData;
}

void ReconstructionOE::start(const Int_t numberOfIterationsForEquilibirum, const Int_t numberOfIterationsInEquilibrium){
    // start the OE image reconstruction

    TBenchmark b;

    // step 1: create initial state s_0 by randomly selecting possible origins for the detected events
    b.Start("S_0");
    State initalState(measurementData->N_dcb);
    initalState.generateRandomOrigins(image->numberOfVoxels, systemMatrixData->systemMatrix);
    statesBeforeEquilibrium.push_back(initalState);
    b.Stop("S_0");

    // step 2: generate new states until equilibrium is reached
    b.Start("UntilE");
    reachEquilibrium(numberOfIterationsForEquilibirum);
    b.Stop("UntilE");

    std::cout << "\nS_0 Creation Time:\t" << b.GetRealTime("S_0") << " s\n";
    std::cout << "\nReach Equilibrium Time:\t" << b.GetRealTime("UntilE") << " s\n";
}

// ##### CALCULATION FUNCTIONS #####
void ReconstructionOE::reachEquilibrium(const Int_t numberOfIterations){

    State nextState = statesBeforeEquilibrium.back();
    nextState.MCMCNextState(systemMatrixData->systemMatrix);
    statesBeforeEquilibrium.push_back(nextState);

//    std::cout << state0.events.at(10) << "\t" << state0.origins.at(10) << "\t " << &state0 << "\n";
//    std::cout << state1.events.at(10) << "\t" << state1.origins.at(10) << "\t " << &state1 << "\n";

}
