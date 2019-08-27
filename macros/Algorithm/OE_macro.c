// ROOT Macro using the Origin Ensemble algorithm for SPCI-Reconstruction
//
// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-08-27


#include <iostream>
#include <TBenchmark.h>

#include "originensemble.h"


// ##### ROOT #####
void OE_macro(){
    TBenchmark b;
    b.Start("total");

    // Specify the location of the measurement file
    TString pathToMeasurements = "../folder/subfolder/*.root";

    // Specify the location of the projections file
    TString pathToProjection = "../folder/subfolder/*.root";

    // create an instance of the reconstruction class
    ReconstructionOE* reco = new ReconstructionOE(pathToMeasurements,
                                                  pathToProjection,
                                                  {-52.5, 52.5, -52.5, 52.5, 5, 10});
    reco->start(10, 1);

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
