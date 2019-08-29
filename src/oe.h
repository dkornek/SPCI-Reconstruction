// oe.h

#pragma once
#include <vector>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>

#include "imagespace.h"
#include "systemmatrix.h"
#include "measurements.h"
#include "state.h"

// ##### RESULTS #####
class ResultsOE{
public:
    ResultsOE();
    ~ResultsOE();

    void plot(TH3F* A_v);
    void plotTransitions(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis);

    TH2F* A_vProject3D = nullptr;

    TCanvas* canvas2D;
    TPad* plot2D;

    TCanvas* canvasTrans;
    TPad* plotTrans;
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
    void sampleEquilibriumStates(const Int_t numberOfIterations);
    void calculateActivity();

    // ##### MEMBERS #####
    Bool_t isCalculationValid;
    std::vector<std::vector<Int_t> > countsInVoxelsOfStates;

    Measurements* measurementData = nullptr;
    SystemMatrix* systemMatrixData = nullptr;
    ImageSpace* image = nullptr;
    ResultsOE* results = nullptr;
    State* state = nullptr;
};
