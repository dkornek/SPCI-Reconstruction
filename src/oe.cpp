// oe.cpp

#include "oe.h"
#include <TBenchmark.h>
#include <TGraph.h>
#include <TStyle.h>
#include "utilities.h"

// ##### RESULTS #####
ResultsOE::ResultsOE(){
    canvas2D = new TCanvas("a2d_c", "OE Reconstruction", 10, 100, 450, 410);
    canvas2D->cd();

    gStyle->SetPalette(53, 0);
    plot2D = new TPad("a2d_p", "", 0, 0, 1, 1);
    plot2D->SetLeftMargin(0.12);
    plot2D->SetRightMargin(0.22);
    plot2D->SetBottomMargin(0.12);
    plot2D->Draw();

    canvasTrans = new TCanvas("trans_c", "Transitions", 500, 100, 450, 450);
    canvasTrans->cd();

    plotTrans = new TPad("trans_p", "", 0, 0, 1, 1);
    plotTrans->SetLeftMargin(0.20);
    plotTrans->SetRightMargin(0.03);
    plotTrans->SetBottomMargin(0.12);
    plotTrans->Draw();
}

ResultsOE::~ResultsOE(){
    delete A_vProject3D;
    delete plotTrans;
    delete plot2D;
    canvasTrans->Close();
    delete canvasTrans;
    canvas2D->Close();
    delete canvas2D;
}

void ResultsOE::plot(TH3F *A_v){
    plot2D->cd();
    A_vProject3D = (TH2F*)A_v->Project3D("yx");
    A_vProject3D->Smooth(1, "k3a");
    A_vProject3D->SetContour(99);
    A_vProject3D->SetStats(kFALSE);
    A_vProject3D->SetTitle("");

    A_vProject3D->GetXaxis()->SetTitleOffset(1);
    A_vProject3D->GetXaxis()->SetNdivisions(5);
    A_vProject3D->GetXaxis()->SetLabelSize(0.06);
    A_vProject3D->GetXaxis()->SetTitleSize(0.06);
    A_vProject3D->GetYaxis()->SetTitleOffset(1);
    A_vProject3D->GetYaxis()->SetNdivisions(5);
    A_vProject3D->GetYaxis()->SetLabelSize(0.06);
    A_vProject3D->GetYaxis()->SetTitleSize(0.06);
    A_vProject3D->GetZaxis()->SetTitle("Probability");
    A_vProject3D->GetZaxis()->SetTitleOffset(1.5);
    A_vProject3D->GetZaxis()->SetNdivisions(3);
    A_vProject3D->GetZaxis()->SetLabelSize(0.06);
    A_vProject3D->GetZaxis()->SetTitleSize(0.06);
    A_vProject3D->GetZaxis()->SetMaxDigits(2);

    A_vProject3D->Draw("COLZ");
    canvas2D->Update();
}

void ResultsOE::plotTransitions(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis){
    Int_t n = xAxis.size();
    if (n != 0){
        Double_t* x = &xAxis[0];
        Double_t* y = &yAxis[0];

        plotTrans->cd();
        TGraph* graph = new TGraph(n, x, y);

        graph->SetTitle(";Number of iterations;Rel. Number of Transitions");
        graph->GetXaxis()->SetLabelSize(0.06);
        graph->GetXaxis()->SetTitleSize(0.06);
        graph->GetYaxis()->SetLabelSize(0.06);
        graph->GetYaxis()->SetTitleSize(0.06);

        graph->Draw("AL");
        canvasTrans->Update();
    }
}

// ##### RECONSTRUCTION #####
ReconstructionOE::ReconstructionOE(const TString pathToMeasurements,
                                   const TString pathToProjections,
                                   const std::vector<Double_t> volume){
    // prepare data for reconstruction using OE

    TBenchmark b;

    // prepare the measurement data
    b.Start("N_dcb");
    measurementData = new Measurements(pathToMeasurements);
    measurementData->createN_dcb();
    measurementData->fillN_dcb(kFALSE);
    b.Stop("N_dcb");

    // prepare the image & system matrix
    b.Start("p_dcb");
    systemMatrixData = new SystemMatrix(pathToProjections);

    isCalculationValid = Utilities::checkForSameNumberOfBins(measurementData->numberOfBins,
                                                             systemMatrixData->numberOfBins);

    if (isCalculationValid){
        systemMatrixData->createSystemMatrix(measurementData->numberOfDetectors);
    }
    b.Stop("p_dcb");

    b.Start("A_v");
    image = new ImageSpace(volume, systemMatrixData->size);
    b.Stop("A_v");

    countsInVoxelsOfStates.reserve(10000);  // reserve enough memory for at least 10000 states

    std::cout << "\nN_dcb Creation Time:\t" << b.GetRealTime("N_dcb") << " s\n";
    std::cout << "\np_dcb Creation Time:\t" << b.GetRealTime("p_dcb") << " s\n";
    std::cout << "\nA_v Creation Time:\t" << b.GetRealTime("A_v") << " s\n";
}

ReconstructionOE::~ReconstructionOE(){
    delete state;
    delete results;
    delete image;
    delete systemMatrixData;
    delete measurementData;
}

void ReconstructionOE::start(const Int_t numberOfIterationsForEquilibirum, const Int_t numberOfIterationsInEquilibrium){
    // start the OE image reconstruction

    TBenchmark b;

    results = new ResultsOE();

    // step 1: create initial state s_0 by randomly selecting possible origins for the detected events
    b.Start("stats");
    b.Start("S_0");
    state = new State(measurementData->N_dcb);
    std::vector<Int_t> countsInVoxel = state->generateRandomOrigins(image->numberOfVoxels,
                                                                    systemMatrixData->systemMatrixVector);
    b.Stop("S_0");
    std::cout << "\nS_0 Creation Time:\t" << b.GetRealTime("S_0") << " s\n";

    // step 2: generate new states until equilibrium is reached
    b.Start("UntilE");
    reachEquilibrium(numberOfIterationsForEquilibirum);
    b.Stop("UntilE");
    std::cout << "\nReach Equilibrium Time:\t" << b.GetRealTime("UntilE") << " s\n";

    // step 3: generate new states in equilibrium
    b.Start("InE");
    sampleEquilibriumStates(numberOfIterationsInEquilibrium);
    b.Stop("InE");
    std::cout << "\nSampling Time:\t" << b.GetRealTime("InE") << " s\n";

    // step 4: calculate "mean state" = mean of counts in each voxel of all sampled states
    calculateActivity();
    b.Stop("stats");

    std::cout << "\nImage reconstruction done.\n";
    std::cout << "\nCalculation Time:\t" << b.GetRealTime("stats") << " seconds\n";

    // step 5: show results
    image->A_v->Scale(1.0 / image->A_v->Integral());
    results->plot(image->A_v);

    results->canvas2D->SaveAs("OE_EmissionDensity.pdf");
    results->canvasTrans->SaveAs("Transitions.pdf");
}

// ##### CALCULATION FUNCTIONS #####
void ReconstructionOE::reachEquilibrium(const Int_t numberOfIterations){
    // generate random states until equilibrium is reached

    std::vector<Double_t> relTransitionsXaxis;
    std::vector<Double_t> relTransitionsYaxis;

    for (Int_t n = 0; n < numberOfIterations; ++n){
        Double_t relTransitions;
        std::vector<Int_t> countsInVoxel = state->MCMCNextState(systemMatrixData->systemMatrixVector,
                                                                systemMatrixData->sensitivities,
                                                                relTransitions);
        relTransitionsXaxis.push_back(n + 1);
        relTransitionsYaxis.push_back(relTransitions);
    }

    results->plotTransitions(relTransitionsXaxis, relTransitionsYaxis);
}

void ReconstructionOE::sampleEquilibriumStates(const Int_t numberOfIterations){
    // generate states in equilibrium

    for (Int_t n = 0; n < numberOfIterations; ++n){
        Double_t relTransitions;
        std::vector<Int_t> countsInVoxel = state->MCMCNextState(systemMatrixData->systemMatrixVector,
                                                                systemMatrixData->sensitivities,
                                                                relTransitions);

        if (n % 1 == 0){
            // save every state, can be changed according to needs
            countsInVoxelsOfStates.push_back(countsInVoxel);
        }
    }
}

void ReconstructionOE::calculateActivity(){
    // calculate the voxel-specific means of the detected counts

    Int_t numberOfStates = countsInVoxelsOfStates.size();
    for (Int_t v = 0; v < image->numberOfVoxels; ++v){

        // calculate the mean counts
        Double_t meanCountsInVoxel = 0.0;
        for (Int_t c = 0; c < numberOfStates; ++c){
            meanCountsInVoxel += countsInVoxelsOfStates.at(c).at(v);
        }
        meanCountsInVoxel /= numberOfStates;

        // calculate the activity
        Double_t sensitivityOfVoxel = systemMatrixData->sensitivities.at(v);
        Double_t activityInVoxel = meanCountsInVoxel / sensitivityOfVoxel;

        // assign mean activity to voxel in image space
        std::array<Int_t, 3> coordinate = image->imageIndices.at(v);
        image->A_v->SetBinContent(coordinate[0], coordinate[1], coordinate[2],
                                  activityInVoxel);
    }

    image->A_v->Scale(1.0);  // bug handling to make Project3D work
}
