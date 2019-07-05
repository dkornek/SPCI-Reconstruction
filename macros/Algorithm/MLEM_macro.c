// ROOT Macro using the Maximum Likelihood Expectation Maximization algorithm for SPCI-Reconstruction
// Standardization is used = faster convergence rate

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-06-01

#include <iostream>
#include <TBenchmark.h>

#include "mlem.h"


// ##### ROOT #####
void MLEM_macro(){
    TBenchmark b;
    b.Start("total");

    // Specify the location of the measurement file
    TString pathToMeasurements = "../folder/subfolder/*.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0+24+26+36_NewDesignNarrow.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0+24+26+36_NewDesignBroad.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0+24+26+36.root";
    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins50/Source0.root";

    // Specify the location of the projections file
    TString pathToProjection = "../folder/subfolder/*.root";
    pathToProjection = "../../data/SystemMatrix/Bins50/SPCIBase441.root";

    // create an instance of the reconstruction class
    ReconstructionMLEM* reco = new ReconstructionMLEM(pathToMeasurements,
                                                      pathToProjection,
                                                      // {-30, 30, -30, 30, 5, 15});  //
                                                      {-50, 50, -50, 50, 5, 10});
    reco->setAccelerator(1.5);
    reco->start(50, 0.9999999);

    b.Stop("total");
    std::cout << "\nTotal Time:\t\t" << b.GetRealTime("total") << " seconds\n";
}

// ##### COMPILE FILE #####
#ifndef __CINT__
int main(){
    MLEM_macro();

    return 0;
}
#endif
