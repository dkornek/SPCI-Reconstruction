// ROOT Macro using the Origin Ensemble algorithm for SPCI-Reconstruction

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-06-01

#include <iostream>
#include <TBenchmark.h>

#include "originensemble.h"


// ##### ROOT #####
void OE_macro(){
    TBenchmark b;
    b.Start("total");

    // Specify the location of the measurement file
    TString pathToMeasurements = "../folder/subfolder/*.root";
    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0+24+26+36.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0+24+26+36_NewDesignBroad.root";

    // Specify the location of the projections file
    TString pathToProjection = "../folder/subfolder/*.root";
    pathToProjection = "../../data/SystemMatrix/Bins50/SPCIBase49.root";

    // create an instance of the reconstruction class
    ReconstructionOE* reco = new ReconstructionOE(pathToMeasurements,
                                                  pathToProjection,
                                                  {-30, 30, -30, 30, 5, 15});  // {-50, 50, -50, 50, 5, 10}
    reco->start(50, 50);

    b.Stop("total");
    std::cout << "\nTotal Time:\t\t" << b.GetRealTime("total") << " seconds\n";
}

// ##### COMPILE FILE #####
#ifndef __CINT__
int main(){
    OE_macro();

    return 0;
}
#endif
