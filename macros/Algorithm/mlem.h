// mlem.h

#pragma once
#include <TH2.h>
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "imagespace.h"
#include "systemmatrix.h"
#include "measurements.h"
#include "utilities.h"


const Bool_t normalizeSpectra = kTRUE;

// ##### RESULTS #####
class ResultsMLEM{
public:
    ResultsMLEM();
    ~ResultsMLEM();

    void plot(TH3F* A_v, std::vector<Double_t> xAxis, std::vector<Double_t> yAxis);

private:
    TH2F* A_vProject3D = nullptr;

    TCanvas* canvas;
    TPad* plot3D;
    TPad* plot2D;
    TPad* plotChi;
};

ResultsMLEM::ResultsMLEM(){
    canvas = new TCanvas("results", "Image Reconstruction", 10, 10, 1350, 450);
    canvas->cd();

    plot3D = new TPad("activity3D", "3D Activity Distribution", 0.01, 0.01, 0.33, 0.99);
    plot3D->Draw();

    plot2D = new TPad("activity2d", "2D Projection", 0.34, 0.01, 0.66, 0.99);
    plot2D->SetRightMargin(0.2);
    plot2D->Draw();

    plotChi = new TPad("chi", "ChiSquare", 0.67, 0.01, 0.99, 0.99);
    plotChi->SetLogy();
    plotChi->SetLogx();
    plotChi->Draw();
}

ResultsMLEM::~ResultsMLEM(){
    delete A_vProject3D;
    delete plotChi;
    delete plot2D;
    delete plot3D;
    delete canvas;
}

void ResultsMLEM::plot(TH3F *A_v, std::vector<Double_t> xAxis, std::vector<Double_t> yAxis){
    plot3D->cd();
    A_v->Draw("BOX2");

    plot2D->cd();
    A_vProject3D = (TH2F*)A_v->Project3D("yx");
    A_vProject3D->SetStats(kFALSE);
    A_vProject3D->SetTitle("2D Projection of 3D Activity Distribution");

    A_vProject3D->GetXaxis()->SetTitleOffset(1);
    A_vProject3D->GetYaxis()->SetTitleOffset(1);
    A_vProject3D->GetZaxis()->SetTitle("Relative Activity");
    A_vProject3D->GetZaxis()->SetTitleOffset(1.8);
    A_vProject3D->Draw("COLZ");

    Int_t n = xAxis.size();
    if (n != 0){
        Double_t* x = &xAxis[0];
        Double_t* y = &yAxis[0];

        plotChi->cd();
        TGraph* graph = new TGraph(n, x, y);
        graph->SetTitle("#chi^{2}-Statistics;Number of iterations;#chi^{2}");
        graph->Draw("AL*");
    }

    canvas->Update();
}

// ##### RECONSTRUCTION #####
class ReconstructionMLEM{
public:
    ReconstructionMLEM(const TString pathToMeasurements,
                       const TString pathToProjections,
                       const std::vector<Double_t> volume);
    ~ReconstructionMLEM();

    void start(const Int_t maxNumberOfIterations, const Double_t stoppingCriterion);
    void setAccelerator(const Double_t a){accelerator = a;}

private:
    // ##### PREPARATION FUNCTIONS #####
    void createProjections();

    // ##### CALCULATION FUNCTIONS #####
    void calculate();
    void projection();
    void backprojection();

    Double_t calculateChiSquare();

    // ##### MEMBERS #####
    Bool_t isCalculationValid;
    Double_t accelerator;

    TH3F* projections = nullptr;

    Measurements* measurementData = nullptr;
    SystemMatrix* systemMatrixData = nullptr;
    ImageSpace* image = nullptr;
    ResultsMLEM* results = nullptr;
};

// ##### PUBLIC FUNCTIONS #####
ReconstructionMLEM::ReconstructionMLEM(const TString pathToMeasurements,
                                       const TString pathToProjections,
                                       const std::vector<Double_t> volume) :
    accelerator(1) {
    // prepare data for reconstruction using ML-EM

    TBenchmark b;

    // prepare the measurement data
    b.Start("N_dcb");
    measurementData = new Measurements(pathToMeasurements);
    measurementData->createN_dcb();
    measurementData->fillN_dcb(normalizeSpectra);
    b.Stop("N_dcb");

    // prepare the image & system matrix
    b.Start("p_dcb");
    systemMatrixData = new SystemMatrix(pathToProjections);

    isCalculationValid = Utilities::checkForSameNumberOfBins(measurementData->numberOfBins,
                                                             systemMatrixData->numberOfBins);

    if (isCalculationValid){
        systemMatrixData->createSystemMatrix(normalizeSpectra,
                                             measurementData->N_dcb,
                                             measurementData->numberOfDetectors);

        systemMatrixData->createP_dcbvPrime(measurementData->N_dcb);
        createProjections();
    }
    b.Stop("p_dcb");

    b.Start("A_v");
    image = new ImageSpace(volume,
                           systemMatrixData->size);
    image->makeA_vHomogeneous();
    b.Stop("A_v");

    std::cout << "\nN_dcb Creation Time:\t" << b.GetRealTime("N_dcb") << " s\n";
    std::cout << "\np_dcb Creation Time:\t" << b.GetRealTime("p_dcb") << " s\n";
    std::cout << "\nA_v Creation Time:\t" << b.GetRealTime("A_v") << " s\n";
}

ReconstructionMLEM::~ReconstructionMLEM(){
    delete projections;
    delete measurementData;
    delete systemMatrixData;
    delete image;
    delete results;
}

void ReconstructionMLEM::start(const Int_t maxNumberOfIterations, const Double_t stoppingCriterion){
    // execute the calculation

    if (!isCalculationValid){
        return;
    }

    TBenchmark b;

    Double_t chiSquare = 0;
    std::vector<Double_t> chiSquareXaxis;
    std::vector<Double_t> chiSquareYaxis;

    results = new ResultsMLEM();

    b.Start("stats");
    Int_t numberOfIterations = 0;
    for (;;){

        if (numberOfIterations == maxNumberOfIterations){
            break;
        }

        calculate();

        if (numberOfIterations % 10 == 0){
            results->plot(image->A_v, chiSquareXaxis, chiSquareYaxis);
        }

        if (numberOfIterations % 5 == 0){
            chiSquare = calculateChiSquare();
            chiSquareXaxis.push_back(numberOfIterations + 1);
            chiSquareYaxis.push_back(chiSquare);
        }

        if (numberOfIterations % 5 == 1){
            if (calculateChiSquare() / chiSquare >= stoppingCriterion){
                break;
            }
        }

        ++numberOfIterations;
    }
    b.Stop("stats");

    // Inform user
    std::cout << "\nImage reconstruction done. Steps: " << numberOfIterations << "\n";
    std::cout << "\nCalculation Time:\t" << b.GetRealTime("stats") << " seconds\n";

    results->plot(image->A_v, chiSquareXaxis, chiSquareYaxis);
}

// ##### PREPARATION FUNCTIONS #####
void ReconstructionMLEM::createProjections(){
    // create the 3D histogram for the forward projection

    Int_t numberOfDetectors = measurementData->numberOfDetectors;
    Int_t numberOfBins = measurementData->numberOfBins;

    projections = new TH3F("projections", "Projections",
                           numberOfDetectors, 0, numberOfDetectors,
                           numberOfDetectors, 0, numberOfDetectors,
                           numberOfBins, 0, numberOfBins);
}

// ##### CALCULATION FUNCTIONS #####
void ReconstructionMLEM::calculate(){
    // execute the maximum likelihood expectation maximization algorithm
    // to calculate the activity distribution A (=A_v)

    // calculate the projection
    projection();

    // calculate the backprojection
    backprojection();
}

void ReconstructionMLEM::projection(){
    // calculates the forward projections

    projections->Reset();
    for (Int_t v = 0; v < image->numberOfVoxels; ++v){

        std::array<Int_t, 3> coordinate = image->imageIndices.at(v);
        Double_t activityInVoxel = image->A_v->GetBinContent(coordinate[0], coordinate[1], coordinate[2]);

        TH3F* p_dcb = (TH3F*)systemMatrixData->systemMatrix->At(v);
        projections->Add(p_dcb, activityInVoxel);
    }
}

void ReconstructionMLEM::backprojection(){
    // create TH3 for activity correction factors

    for (Int_t v = 0; v < image->numberOfVoxels; ++v){

        std::array<Int_t, 3> coordinate = image->imageIndices.at(v);
        Double_t activityInVoxel = image->A_v->GetBinContent(coordinate[0], coordinate[1], coordinate[2]);

        TH3F* p_dcbPrime = (TH3F*)systemMatrixData->p_dcbvPrime->At(v);

        Double_t correctionFactor = 0.0;
        for (Int_t d = 1; d <= measurementData->numberOfDetectors; ++d){
            for (Int_t c = 1; c <= measurementData->numberOfDetectors; ++c){
                for (Int_t bin = 1; bin <= measurementData->numberOfBins; ++bin){

                    Double_t N_dcbPrimeBinContent = projections->GetBinContent(d, c, bin);
                    if (N_dcbPrimeBinContent != 0){

                        Double_t p_dcbPrimeBinContent = p_dcbPrime->GetBinContent(d, c, bin);
                        if (p_dcbPrimeBinContent != 0){

                            correctionFactor += p_dcbPrimeBinContent / N_dcbPrimeBinContent;
                        }
                    }
                }
            }
        }

        correctionFactor = correctionFactor / systemMatrixData->sensitivities[v];
        correctionFactor = std::pow(correctionFactor, accelerator);

        image->A_v->SetBinContent(coordinate[0], coordinate[1], coordinate[2],
                                  correctionFactor * activityInVoxel);
    }
}

Double_t ReconstructionMLEM::calculateChiSquare(){
    // The Chi²-statistics is used to stop the MLEM algorithm before convergence
    // The curve is monotonically decreasing -> Change of rate is proportional to convergence

    Double_t ChiSquareTestVariable = 0.0;
    for (Int_t d = 1; d <= measurementData->numberOfDetectors; ++d){
        for (Int_t c = 1; c <= measurementData->numberOfDetectors; ++c){
            for (Int_t bin = 1; bin <= measurementData->numberOfBins; ++bin){

                Double_t projectionBinContent = projections->GetBinContent(d, c, bin);
                if (projectionBinContent != 0){

                    Double_t N_dcbBinContent = measurementData->N_dcb->GetBinContent(d, c, bin);
                    if (N_dcbBinContent != 0){

                        ChiSquareTestVariable += std::pow(N_dcbBinContent - projectionBinContent, 2) / projectionBinContent;
                    } else{

                        ChiSquareTestVariable += projectionBinContent;
                    }
                }
            }
        }
    }

    return ChiSquareTestVariable;
}
