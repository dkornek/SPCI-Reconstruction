// mlem.h

#pragma once
#include <TH2.h>
#include <TH3.h>
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "imagespace.h"
#include "systemmatrix.h"
#include "measurements.h"
#include "utilities.h"

// It should always be kFALSE
// kTRUE will probably never work with noise or real measurements
const Bool_t normalizeSpectra = kFALSE;

// ##### RESULTS #####
class ResultsMLEM{
public:
    ResultsMLEM();
    ~ResultsMLEM();

    void plotActivity(TH3F* A_v);
    void plotChiSquare(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis);
    void plotLogLikelihood(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis);

    TH2F* A_vProject3D = nullptr;

    TCanvas* canvas2D;
    TPad* plot2D;

    TCanvas* canvasChi;
    TPad* plotChi;

    TCanvas* canvasLogLike;
    TPad* plotLogLike;
};

// ##### RECONSTRUCTION #####
class ReconstructionMLEM{
public:
    ReconstructionMLEM(const TString pathToMeasurements,
                       const TString pathToProjections,
                       const std::vector<Double_t> volume);
    ~ReconstructionMLEM();

    void start(const Int_t maxNumberOfIterations);
    void setAccelerator(const Double_t a){accelerator = a;}

private:
    // ##### PREPARATION FUNCTIONS #####
    void createProjections();

    // ##### CALCULATION FUNCTIONS #####
    void calculate();
    void projection();
    void backprojection();

    Double_t calculateChiSquare();
    Double_t calculateLogLike();

    // ##### MEMBERS #####
    Bool_t isCalculationValid;
    Double_t accelerator;

    TH3F* projections = nullptr;

    Measurements* measurementData = nullptr;
    SystemMatrix* systemMatrixData = nullptr;
    ImageSpace* image = nullptr;
    ResultsMLEM* results = nullptr;
};
