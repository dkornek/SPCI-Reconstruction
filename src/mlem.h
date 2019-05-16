// mlem.h
// image reconstruction using maximum likelihood expectation maximization algorithm

#pragma once
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TH3.h>


class ReconstructionMLEM{
public:
    ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections);
    ~ReconstructionMLEM();
    void start(Int_t maxNumberOfIterations);

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
    void projection(TH3F& p);
    void backprojection(TH3F& p, TH3F& bp);

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
    TList* p_dcbvPrime;  // p_dcbv * N_dcb
    std::vector<Double_t> p_dcbvSum;

    Int_t NbinsProjections;

    TIter* nextVoxel;
    TKey* keyVoxel;
    TKey* keyS_dcVoxel;
};
