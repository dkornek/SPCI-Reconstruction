// image_processing.h
// class for saving reconstructed images

#pragma once
#include <TROOT.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>

class ActivityDistribution{
public:
    ActivityDistribution(TH3F* a, TString fileName);
    ~ActivityDistribution();
    void save3DHistogram();
    void save2DProjection();
    void save2DSlices();

    TFile* saveData;
    TH3F* activity;
};
