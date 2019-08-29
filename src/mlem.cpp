// mlem.cpp

#include "mlem.h"
#include <TGraph.h>

// ##### RESULTS #####
ResultsMLEM::ResultsMLEM(){
    canvas2D = new TCanvas("a2d_c", "MLEM Reconstruction", 10, 10, 450, 410);
    canvas2D->cd();

    gStyle->SetPalette(53, 0);
    plot2D = new TPad("a2d_p", "", 0, 0, 1, 1);
    plot2D->SetLeftMargin(0.12);
    plot2D->SetRightMargin(0.22);
    plot2D->SetBottomMargin(0.12);
    plot2D->Draw();

    canvasChi = new TCanvas("chi_c", "Chi Square", 340, 40, 450, 450);
    canvasChi->cd();

    plotChi = new TPad("chi_p", "", 0, 0, 1, 1);
    plotChi->SetLeftMargin(0.20);
    plotChi->SetRightMargin(0.03);
    plotChi->SetBottomMargin(0.12);
    plotChi->SetLogy();
    plotChi->Draw();

    canvasLogLike = new TCanvas("log_c", "log-Likelihood", 670, 70, 450, 450);
    canvasLogLike->cd();

    plotLogLike = new TPad("log_p", "", 0, 0, 1, 1);
    plotLogLike->SetLeftMargin(0.20);
    plotLogLike->SetRightMargin(0.03);
    plotLogLike->SetBottomMargin(0.12);
    plotLogLike->Draw();
}

ResultsMLEM::~ResultsMLEM(){
    delete A_vProject3D;
    delete plot2D;
    canvas2D->Close();
    delete canvas2D;

    delete plotChi;
    canvasChi->Close();
    delete canvasChi;

    delete plotLogLike;
    canvasLogLike->Close();
    delete canvasLogLike;
}

void ResultsMLEM::plotActivity(TH3F *A_v){
    plot2D->cd();
    A_vProject3D = (TH2F*)A_v->Project3D("yx");
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
    A_vProject3D->GetZaxis()->SetLabelSize(0.06);
    A_vProject3D->GetZaxis()->SetTitleOffset(1.5);
    A_vProject3D->GetZaxis()->SetNdivisions(4);
    A_vProject3D->GetZaxis()->SetMaxDigits(2);
    A_vProject3D->GetZaxis()->SetTitle("Probability");
    A_vProject3D->GetZaxis()->SetTitleSize(0.06);
    A_vProject3D->Draw("COLZ");
    canvas2D->Update();
    canvas2D->WaitPrimitive();
}

void ResultsMLEM::plotChiSquare(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis){

    Int_t n = xAxis.size();
    if (n != 0){
        Double_t* x = &xAxis[0];
        Double_t* y = &yAxis[0];

        plotChi->cd();
        TGraph* graph = new TGraph(n, x, y);
        graph->GetXaxis()->SetLabelSize(0.06);
        graph->GetXaxis()->SetTitleSize(0.06);
        graph->GetYaxis()->SetLabelSize(0.06);
        graph->GetYaxis()->SetTitleSize(0.06);
        graph->GetYaxis()->SetMaxDigits(2);

        graph->SetMarkerStyle(5);
        graph->SetTitle(";Number of iterations;#chi^{2}");
        graph->Draw("AL");
        canvasChi->Update();
    }
}

void ResultsMLEM::plotLogLikelihood(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis){

    Int_t n = xAxis.size();
    if (n != 0){
        Double_t* x = &xAxis[0];
        Double_t* y = &yAxis[0];

        plotLogLike->cd();
        TGraph* graph = new TGraph(n, x, y);
        graph->GetXaxis()->SetLabelSize(0.06);
        graph->GetXaxis()->SetTitleSize(0.06);
        graph->GetYaxis()->SetLabelSize(0.06);
        graph->GetYaxis()->SetTitleSize(0.06);
        graph->GetYaxis()->SetMaxDigits(2);

        graph->SetMarkerStyle(5);
        graph->SetTitle(";Number of iterations;-log L(#lambda)");
        graph->Draw("AL");
        canvasLogLike->Update();
    }
}

// ##### RECONSTRUCTION #####
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
    delete results;
    delete image;
    delete projections;
    delete measurementData;
    delete systemMatrixData;
}

void ReconstructionMLEM::start(const Int_t maxNumberOfIterations){
    // execute the calculation

    if (!isCalculationValid){
        return;
    }

    TBenchmark b;

    Double_t chiSquare = 0;
    std::vector<Double_t> chiSquareXaxis;
    std::vector<Double_t> chiSquareYaxis;

    Double_t logLike = 0;
    std::vector<Double_t> logLikeXaxis;
    std::vector<Double_t> logLikeYaxis;

    results = new ResultsMLEM();

    b.Start("stats");
    Int_t numberOfIterations = 0;
    for (;;){

        if (numberOfIterations == maxNumberOfIterations){
            break;
        }

        calculate();

        if (numberOfIterations % 1 == 0){
            chiSquare = calculateChiSquare();
            chiSquareXaxis.push_back(numberOfIterations + 1);
            chiSquareYaxis.push_back(chiSquare);

            logLike = calculateLogLike();
            logLikeXaxis.push_back(numberOfIterations + 1);
            logLikeYaxis.push_back(logLike);
        }

        ++numberOfIterations;
    }
    b.Stop("stats");

    // Inform user
    std::cout << "\nImage reconstruction done. Steps: " << numberOfIterations << "\n";
    std::cout << "\nCalculation Time:\t" << b.GetRealTime("stats") << " seconds\n";

    image->A_v->Scale(1.0 / image->A_v->Integral());
    results->plotActivity(image->A_v);
    results->plotChiSquare(chiSquareXaxis, chiSquareYaxis);
    results->plotLogLikelihood(logLikeXaxis, logLikeYaxis);

    results->canvas2D->SaveAs("MLEM_EmissionDensitiy.pdf");
    results->canvasLogLike->SaveAs("Likelihood.pdf");
    results->canvasChi->SaveAs("ChiSquare.pdf");
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
    // Calculate the Chi Square Statistics

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

Double_t ReconstructionMLEM::calculateLogLike(){
    // Calculate the Log-Likelihood-Function

    Double_t LogLike = 0.0;
    for (Int_t d = 1; d <= measurementData->numberOfDetectors; ++d){
        for (Int_t c = 1; c <= measurementData->numberOfDetectors; ++c){
            for (Int_t bin = 1; bin <= measurementData->numberOfBins; ++bin){

                Double_t projectionBinContent = projections->GetBinContent(d, c, bin);
                if (projectionBinContent != 0){

                    Double_t N_dcbBinContent = measurementData->N_dcb->GetBinContent(d, c, bin);
                    Double_t summand = projectionBinContent - N_dcbBinContent * std::log(projectionBinContent);

                    LogLike += summand;
                }
            }
        }
    }

    return LogLike;
}
