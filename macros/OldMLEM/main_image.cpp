#include "ImageReconstruction.h"
#include "TBenchmark.h"

using namespace std;


void main_image() {
    TBenchmark bench;

    bench.Start("total");
    ImageReconstruction reco(2,     // dim
                             441,   // voxels
                             5.0,   // voxel_len / mm
                             4,     // det_n
                             4,     // det_m
                             50,    // bins
                             1);    // rebin factor

    //set base file, input file and output file
    reco.SetBaseFileName("../../data/SystemMatrix/Bins50/SPCIBase441.root");
//    reco.SetInputFileName("../../data/Measurements/*.root");
    reco.SetInputFileName("../../data/Measurements/SPCIBase441/Bins50/SourceHZDR.root");
    reco.SetOutputFileName("Result.root");

    //start reconstruction
    reco.CalculateA_v(5000);          // number of iterations

    bench.Stop("total");
    std::cout << "\nTotal Time: " << bench.GetRealTime("total") << " seconds\n";
}
