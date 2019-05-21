// mlem.h
// image reconstruction using maximum likelihood expectation maximization algorithm

#pragma once
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TH3.h>
#include <TGraph.h>


class ReconstructionMLEM{
public:
    ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections);
    ~ReconstructionMLEM();
    void start(Int_t maxNumberOfIterations, Double_t stoppingCriterium);
    void setAccelerator(Double_t a){accelerator = a;}
    void setActivityThreshold(Double_t t){activityThreshold = t;}
    void setImageVolume(std::vector<Double_t> v){imageVolume = v;}

    // Activity distribution
    TH3F* A_v;
    TList* A_vProject3DSteps;
    TGraph* plotChi;

private:
    // ##### UTILITIES #####
    void openFile(TString fileType, TString pathToFile);
    void getDetectorIndices(TString nameOfSpectrum, Int_t& d, Int_t& c);
    void getImageSpaceIndices(TString titleOfVoxel, Int_t& x, Int_t& y, Int_t& z);
    void fillImageSpaceIndices(Int_t x, Int_t y, Int_t z);
    void fillN_dcb();
    void makeA_vHomogeneous();
    void checkValidity();
    Double_t calculateChiSquareStatistic();
    void plotActivity(Int_t step);
    void plotChiStatistics(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis);

    // ##### PREPARATION FUNCTIONS #####
    void createN_dcb();
    void createP_dcbv();
    void createA_v();

    // ##### CALCULATION FUNCTIONS #####
    void calculate();
    void projection();
    void backprojection();

    // ##### MEMBERS #####
    Double_t accelerator;
    Double_t activityThreshold;
    std::vector<Double_t> imageVolume;  // in mm, xMin, xMax, yMin, yMax, zMin, zMax
    std::vector<std::array<Int_t, 3> > imageIndices;

    Bool_t isCalculationValid;

    // Measurement Data
    TFile* measurementsFile;

    TH3F* N_dcb;

    Int_t numberOfDetectors;
    Int_t NbinsMeasurements;

    TIter* nextS_dcMeasurements;
    TKey* keyS_dcMeasurements;

    // Projection Data
    TFile* projectionsFile;

    TList* p_dcbv;
    TList* p_dcbvPrime;  // p_dcbv * N_dcb
    std::vector<Double_t> p_dcbvSum;

    TH3F* projections;

    Int_t NbinsProjections;

    TIter* nextVoxel;
    TKey* keyVoxel;
    TKey* keyS_dcVoxel;
};
