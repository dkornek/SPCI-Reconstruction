// image_processing.h
// class for saving reconstructed images

#pragma once

#include "TFile.h"
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>

class ActivityDistribution{
public:
    ActivityDistribution(TString fileName);
    ~ActivityDistribution();
    void save3DHistogram(TH3F* activity3D);
    void saveChiStatistics(TGraph* chi);
    void save2DSlices(TH3F* activity3D);
    void save2DSteps(TList* steps);
    void saveAll(TH3F* a3D, TList* s, TGraph* chi);

    TFile* saveData;
};
