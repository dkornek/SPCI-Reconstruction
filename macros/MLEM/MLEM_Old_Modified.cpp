// ROOT Macro using the Maximum Likelihood Expectation Maximization algorithm for SPCI-Reconstruction
// Standardization is used = faster convergence rate

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-05-27


#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2.h>
#include <TH3.h>
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TGraph.h>


// ##### RECONSTRUCTION #####
class ReconstructionMLEM{
public:
    ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections);
    ~ReconstructionMLEM();
    void start(Int_t maxNumberOfIterations, Double_t stoppingCriterion, Bool_t doPlotting);
    void setAccelerator(Double_t a){accelerator = a;}
    void setActivityThreshold(Double_t t){activityThreshold = t;}
    void setImageVolume(std::vector<Double_t> v){imageVolume = v;}

    // Activity distribution
    TH3F* A_v;
    TH2F* A_vProject3D;

private:
    // ##### UTILITIES #####
    void openFile(TString fileType, TString pathToFile);
    void getDetectorIndices(TString nameOfSpectrum, Int_t& d, Int_t& c);
    void getImageSpaceIndices(TString titleOfVoxel, Int_t& x, Int_t& y, Int_t& z);
    void fillImageSpaceIndices(Int_t x, Int_t y, Int_t z);
    void fillN_dcb();
    void makeA_vHomogeneous();
    void checkValidity();
    Double_t calculateChiSquareStatistic();
    void plotActivity();
    void plotChiStatistics(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis);

    // ##### PREPARATION FUNCTIONS #####
    void createN_dcb();
    void createP_dcbv();
    void createA_v();

    // ##### CALCULATION FUNCTIONS #####
    void calculate();
    void projection();
    void backprojection();

    // ##### MEMBERS #####
    Double_t accelerator;
    Double_t activityThreshold;

    Bool_t isCalculationValid;

    TCanvas* canvasResults;
    TPad* plot3D;
    TPad* plot2D;
    TPad* plotChi;

    // Activity Data
    std::vector<Double_t> imageVolume;  // in mm, xMin, xMax, yMin, yMax, zMin, zMax
    std::vector<std::array<Int_t, 3> > imageIndices;

    // Measurement Data
    TFile* measurementsFile;

    TH3F* N_dcb;

    Int_t numberOfDetectors;
    Int_t NbinsMeasurements;

    TIter* nextS_dcMeasurements;
    TKey* keyS_dcMeasurements;

    // Projection Data
    TFile* projectionsFile;

    TList* p_dcbv;
    TList* p_dcbvPrime;  // p_dcbv * N_dcb
    std::vector<Double_t> p_dcbvSum;  // p_v = SUM_{dcb} (p_dcbv)

    TH3F* projections;

    Int_t NbinsProjections;

    TIter* nextVoxel;
    TKey* keyVoxel;
    TKey* keyS_dcVoxel;
};

// ##### PUBLIC FUNCTIONS #####
ReconstructionMLEM::ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections) :
    accelerator(1), activityThreshold(0),
    canvasResults(nullptr), plot3D(nullptr), plot2D(nullptr), plotChi(nullptr){
    // prepare all data to start the reconstruction

    // open files
    openFile("MEASUREMENT", pathToMeasurements);
    nextS_dcMeasurements = new TIter(measurementsFile->GetListOfKeys());

    openFile("PROJECTION", pathToProjections);
    nextVoxel = new TIter(projectionsFile->GetListOfKeys());

    TBenchmark b;

    // prepare the measurement data N_dcb
    // N_dcb = number of events in detector pair (d, c) separated by bins
    b.Start("N_dcb");
    createN_dcb();
    b.Stop("N_dcb");

    // prepare the probabilities p_dcbv
    // p_dcbv = voxel v specific probabilitiy to contribute one count to N_dcb;
    // there are v 3d histograms p_dcb in the list p_dcbv
    b.Start("p_dcb");
    createP_dcbv();
    b.Stop("p_dcb");

    std::cout << "\nN_dcb Creation Time:\t" << b.GetRealTime("N_dcb") << " s\n";
    std::cout << "\np_dcb Creation Time:\t" << b.GetRealTime("p_dcb") << " s\n";
}

ReconstructionMLEM::~ReconstructionMLEM(){

    delete plot3D;
    delete plot2D;
    delete plotChi;
    delete canvasResults;

    delete A_v;
    delete A_vProject3D;

    delete keyS_dcMeasurements;
    delete nextS_dcMeasurements;
    delete N_dcb;
    delete measurementsFile;

    delete keyS_dcVoxel;
    delete keyVoxel;
    delete nextVoxel;
    delete projections;
    p_dcbv->Delete();
    delete p_dcbv;
    p_dcbvPrime->Delete();
    delete p_dcbvPrime;
    delete projectionsFile;
}

void ReconstructionMLEM::start(Int_t maxNumberOfIterations, Double_t stoppingCriterion, Bool_t doPlotting){
    // execute the calculation

    checkValidity();
    if (!isCalculationValid){
        return;
    }

    TBenchmark b;

    // prepare the homogeneous activity distribution A_v
    // A_v = activity in voxel v of the image space
    b.Start("A_v");
    createA_v();
    b.Stop("A_v");
    std::cout << "\nA_v Creation Time:\t" << b.GetRealTime("A_v") << " s\n";

    Double_t chiSquare = 0;
    std::vector<Double_t> chiSquareXaxis;
    std::vector<Double_t> chiSquareYaxis;

    b.Start("stats");
    Int_t numberOfIterations = 0;
    for (;;){

        if (numberOfIterations == maxNumberOfIterations){
            break;
        }

        calculate();

        if (doPlotting && (numberOfIterations % 10 == 0)){
            plotActivity();

            if (stoppingCriterion){
                plotChiStatistics(chiSquareXaxis, chiSquareYaxis);
            }
        }

        if (stoppingCriterion){

            if (numberOfIterations % 5 == 0){
                chiSquare = calculateChiSquareStatistic();
                chiSquareXaxis.push_back(numberOfIterations + 1);
                chiSquareYaxis.push_back(chiSquare);
            }

            if (numberOfIterations % 5 == 1){

                Double_t ratio = calculateChiSquareStatistic() / chiSquare;
                std::cout << ratio << "\n";
                if (ratio >= stoppingCriterion){
                    break;
                }
            }
        }

        ++numberOfIterations;
    }
    b.Stop("stats");

    // Inform user
    std::cout << "\nImage reconstruction done. Steps: " << numberOfIterations << "\n";
    std::cout << "\nCalculation Time:\t" << b.GetRealTime("stats") << " seconds\n";

    // show results
    plotActivity();
    if (stoppingCriterion){
        plotChiStatistics(chiSquareXaxis, chiSquareYaxis);
    }
}

// ##### PRIVATE FUNCTIONS #####
// ##### UTILITIES #####
void ReconstructionMLEM::openFile(TString fileType, TString pathToFile){
    // open file of specified type

    if (fileType == "MEASUREMENT"){
        measurementsFile = new TFile(pathToFile, "READ");
        return;
    }

    if (fileType == "PROJECTION"){
        projectionsFile = new TFile(pathToFile, "READ");
        return;
    }
}

void ReconstructionMLEM::getDetectorIndices(TString nameOfSpectrum, Int_t &d, Int_t &c){
    // extract index of both the detectors from the name of the spectrum

    d = ((TString)nameOfSpectrum(0, 2)).Atoi();
    c = ((TString)nameOfSpectrum(2, 2)).Atoi();
}

void ReconstructionMLEM::getImageSpaceIndices(TString titleOfVoxel, Int_t &x, Int_t &y, Int_t &z){
    // extract location in image space from the title (=name of directory) of the voxel v

    // prepare title
    TString title = titleOfVoxel.Copy();
    title.Remove(0, 1);
    title.Replace(title.Length() - 1, 1, ',');

    // actual coordinates
    std::vector<Int_t> location = { 0, 0, 0 };
    for (Int_t i = 0; i < 3; ++i){
        Int_t positionOfSeparator = title.First(',');
        Int_t lengthOfString = title.Last(',');

        location[i] =  ((TString)title(0, positionOfSeparator)).Atoi();
        title = title(positionOfSeparator + 1, lengthOfString - positionOfSeparator + 1);
    }

    x = location[0];
    y = location[1];
    z = location[2];
}

void ReconstructionMLEM::fillImageSpaceIndices(Int_t x, Int_t y, Int_t z){
    // get a vector of the image space indices

    for (Int_t zIndex = 0; zIndex <= z; ++zIndex){
        for (Int_t yIndex = 0; yIndex <= y; ++yIndex){
            for (Int_t xIndex = 0; xIndex <= x; ++xIndex){

                std::array<Int_t, 3> coordinate = { xIndex + 1,
                                                    yIndex + 1,
                                                    zIndex + 1};
                imageIndices.push_back(coordinate);
            }
        }
    }
}

void ReconstructionMLEM::fillN_dcb(){
    // fill N_dcb with number of events extracted from the spectra

    nextS_dcMeasurements->Reset();
    while ((keyS_dcMeasurements = (TKey*)nextS_dcMeasurements->Next())){
        // iterate through every measurement spectrum

        Int_t detectorD;
        Int_t detectorC;
        TString nameOfS_dc = keyS_dcMeasurements->GetName();
        getDetectorIndices(nameOfS_dc, detectorD, detectorC);

        // get spectrum
        TH1F* S_dc = (TH1F*)measurementsFile->Get(nameOfS_dc);

        Double_t integral = S_dc->Integral();
        if (integral != 0){
            S_dc->Scale(1.0 / integral);  // faster convergence

            // fill N_dcb
            for (Int_t bin = 1; bin <= NbinsMeasurements; ++bin){
                N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                     S_dc->GetBinContent(bin));
            }

        } else{

            // fill N_dcb
            for (Int_t bin = 1; bin <= NbinsMeasurements; ++bin){
                N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                     0);
            }
        }

        delete S_dc;
    }
}

void ReconstructionMLEM::makeA_vHomogeneous(){
    // set all bin contents in A_v = 1

    for (UInt_t v = 0; v < imageIndices.size(); ++v){
        // iterate through all voxels

        std::array<Int_t, 3> coordinate = imageIndices.at(v);
        A_v->SetBinContent(coordinate[0], coordinate[1], coordinate[2],
                           1.0);
    }
}

void ReconstructionMLEM::checkValidity(){
    // calculation is only possible if binning pattern in the measurement and
    // projection spectra is identical

    if (NbinsMeasurements == NbinsProjections){
        isCalculationValid = kTRUE;

    } else{
        std::cout << "Calculation is not valid.\nBinning pattern is not consistent\n";
        isCalculationValid = kFALSE;
    }
}

Double_t ReconstructionMLEM::calculateChiSquareStatistic(){
    // The Chi²-statistics is used to stop the MLEM algorithm before convergence
    // The curve is monotonically decreasing -> Change of rate is proportional to convergence

    Double_t ChiSquareTestVariable = 0.0;
    for (Int_t d = 1; d <= numberOfDetectors; ++d){
        for (Int_t c = 1; c <= numberOfDetectors; ++c){
            for (Int_t bin = 1; bin <= NbinsProjections; ++bin){

                Double_t projectionContent = projections->GetBinContent(d, c, bin);
                if (projectionContent != 0){

                    Double_t N_dcbContent = N_dcb->GetBinContent(d, c, bin);
                    if (N_dcbContent != 0){

                        ChiSquareTestVariable += std::pow(N_dcbContent - projectionContent, 2) / projectionContent;
                    } else{

                        ChiSquareTestVariable += projectionContent;
                    }
                }
            }
        }
    }

    return ChiSquareTestVariable;
}

void ReconstructionMLEM::plotActivity(){
    // show results of the activity

    if (canvasResults == nullptr){
        canvasResults = new TCanvas("results", "Image Reconstruction", 10, 10, 1350, 450);
    }
    canvasResults->cd();

    if (plot3D == nullptr){
        plot3D = new TPad("activity3D", "3D Activity Distribution", 0.01, 0.01, 0.33, 0.99);
        plot3D->Draw();
    }

    if (plot2D == nullptr){
        plot2D = new TPad("activity2d", "2D Projection", 0.34, 0.01, 0.66, 0.99);
        plot2D->SetRightMargin(0.2);
        plot2D->Draw();
    }

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

    canvasResults->Update();
}

void ReconstructionMLEM::plotChiStatistics(std::vector<Double_t> xAxis, std::vector<Double_t> yAxis){
    // show results of Chi test variable

    if (canvasResults == nullptr){
        canvasResults = new TCanvas("results", "Image Reconstruction", 10, 10, 1350, 450);
    }
    canvasResults->cd();

    if (plotChi == nullptr){
        plotChi = new TPad("chi", "ChiSquare", 0.67, 0.01, 0.99, 0.99);
        plotChi->SetLogy();
        plotChi->SetLogx();
        plotChi->Draw();
    }

    Int_t n = xAxis.size();
    if (n == 0){
        return;
    }

    Double_t* x = &xAxis[0];
    Double_t* y = &yAxis[0];

    plotChi->cd();
    TGraph* graph = new TGraph(n, x, y);
    graph->SetTitle("#chi^{2}-Statistics;Number of iterations;#chi^{2}");
    graph->Draw("AL*");

    canvasResults->Update();
}

// ##### PREPARATION FUNCTIONS #####
void ReconstructionMLEM::createN_dcb(){
    // create empty 3d histogram of N_dcb
    // x = detector d, y = detector c, z = bins in spectrum

    // histogram is based on the last spectrum
    // title of spectra is expected to be of "ddcc" format
    Int_t detectorD;
    Int_t detectorC;
    TString nameOfLastSpectrum;

    nameOfLastSpectrum = measurementsFile->GetListOfKeys()->Last()->GetName();
    getDetectorIndices(nameOfLastSpectrum, detectorD, detectorC);

    // get number of detectors
    numberOfDetectors = std::max(detectorD, detectorC) + 1;

    // get number of bins
    NbinsMeasurements = ((TH1F*)measurementsFile->Get(nameOfLastSpectrum))->GetNbinsX();

    // create the 3d histogram
    N_dcb = new TH3F("N_dcb", "Number of events in bin b in S_dc",
                     numberOfDetectors, 0, numberOfDetectors,
                     numberOfDetectors, 0, numberOfDetectors,
                     NbinsMeasurements, 0, NbinsMeasurements);
    fillN_dcb();

    // TCanvas* N_dcbCanvas = new TCanvas("n_dcbCanvas", "Measurements", 700, 700);
    // N_dcb->Draw("BOX2Z");
    // N_dcbCanvas->Update();
}

void ReconstructionMLEM::createP_dcbv(){
    // create list of 3d histograms containing the probabilities for each voxel

    p_dcbv = new TList();
    p_dcbvPrime = new TList();

    // get number of bins
    TString nameOfLastVoxel = projectionsFile->GetListOfKeys()->Last()->GetName();
    TDirectory* dirOfLastVoxel = (TDirectory*)projectionsFile->Get(nameOfLastVoxel);
    TString nameOfLastS_dc = dirOfLastVoxel->GetListOfKeys()->Last()->GetName();
    NbinsProjections = ((TH1F*)dirOfLastVoxel->Get(nameOfLastS_dc))->GetNbinsX();
    delete dirOfLastVoxel;

    checkValidity();
    if (!isCalculationValid){
        return;
    }

    // iterate through all voxels v
    nextVoxel->Reset();
    while ((keyVoxel = (TKey*)nextVoxel->Next())){
        TString nameOfVoxel = keyVoxel->GetName();

        TString description;
        description.Form("Probabilities to contribute one count to N_dcb from voxel %s", nameOfVoxel.Data());

        TString internalName;
        internalName.Form("p_dcb_%s", nameOfVoxel.Data());
        TH3F* p_dcb = new TH3F(internalName, description,
                               numberOfDetectors, 0, numberOfDetectors,
                               numberOfDetectors, 0, numberOfDetectors,
                               NbinsProjections, 0, NbinsProjections);

        internalName.Form("p_dcb_%s_prime", nameOfVoxel.Data());
        TH3F* p_dcbPrime = new TH3F(internalName, description,
                                    numberOfDetectors, 0, numberOfDetectors,
                                    numberOfDetectors, 0, numberOfDetectors,
                                    NbinsProjections, 0, NbinsProjections);

        TDirectory* dirOfVoxel = (TDirectory*)projectionsFile->Get(nameOfVoxel);
        TIter nextS_dc(dirOfVoxel->GetListOfKeys());
        while ((keyS_dcVoxel = (TKey*)nextS_dc())){
            // iterate through all spectra in voxel v

            TString nameOfS_dc = keyS_dcVoxel->GetName();
            TH1F* S_dc = (TH1F*)dirOfVoxel->Get(nameOfS_dc);
            S_dc->Scale(1.0 / S_dc->Integral());  // faster convergence

            Int_t detectorD;
            Int_t detectorC;
            getDetectorIndices(nameOfS_dc(3, 4), detectorD, detectorC);
            for (Int_t bin = 1; bin <= NbinsProjections; ++bin){

                Double_t N_dcbBinContent = N_dcb->GetBinContent(detectorD + 1, detectorC + 1, bin);
                if (N_dcbBinContent == 0){

                    p_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         0);

                    p_dcbPrime->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                              0);
                } else{

                    Double_t p_dcbBinContent = S_dc->GetBinContent(bin);
                    Double_t p_dcbPrimeBinContent = p_dcbBinContent * N_dcbBinContent;

                    p_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         p_dcbBinContent);

                    p_dcbPrime->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                              p_dcbPrimeBinContent);
                }
            }

            delete S_dc;
        }

        p_dcbv->AddLast(p_dcb);
        p_dcbvPrime->AddLast(p_dcbPrime);
        p_dcbvSum.push_back(p_dcb->Integral());

        delete dirOfVoxel;
    }

    // projections histogram used for the calculation
    projections = new TH3F("projections", "Projections",
                           numberOfDetectors, 0, numberOfDetectors,
                           numberOfDetectors, 0, numberOfDetectors,
                           NbinsProjections, 0, NbinsProjections);
}

void ReconstructionMLEM::createA_v(){
    // create the activity distribution A (=A_v)
    // all voxels of the image space contain an activity A_v

    // the image space is based on the last folder title in the projections
    // title is expected to have the form "(X,Y,Z)"
    Int_t indexX;
    Int_t indexY;
    Int_t indexZ;

    // get volume
    TString titleOfLastVoxel;
    titleOfLastVoxel = projectionsFile->GetListOfKeys()->Last()->GetTitle();
    getImageSpaceIndices(titleOfLastVoxel, indexX, indexY, indexZ);
    fillImageSpaceIndices(indexX, indexY, indexZ);

    if (imageVolume.size() != 6){
        imageVolume = {0, static_cast<Double_t>(indexX + 1),
                       0, static_cast<Double_t>(indexY + 1),
                       0, static_cast<Double_t>(indexZ + 1)};
    }

    // create the 3d histogram
    A_v = new TH3F("A_v", "3D Activity Distribution",
                   indexX + 1, imageVolume[0], imageVolume[1],
                   indexY + 1, imageVolume[2], imageVolume[3],
                   indexZ + 1, imageVolume[4], imageVolume[5]);

    makeA_vHomogeneous();

    // Plotting options
    A_v->GetXaxis()->SetTitle("#it{x} / mm");
    A_v->GetXaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * indexX), "X");

    A_v->GetYaxis()->SetTitle("#it{y} / mm");
    A_v->GetXaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * indexY), "Y");

    A_v->GetZaxis()->SetTitle("#it{z} / mm");
    A_v->GetZaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * indexZ), "Z");
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
    // calculates the projection

    projections->Reset();
    for (UInt_t v = 0; v < imageIndices.size(); ++v){
        std::array<Int_t, 3> coordinate = imageIndices.at(v);

        TH3F* p_dcb = (TH3F*)p_dcbv->At(v);
        Double_t activityInVoxel = A_v->GetBinContent(coordinate[0], coordinate[1], coordinate[2]);
        projections->Add(p_dcb, activityInVoxel);
    }
}

void ReconstructionMLEM::backprojection(){
    // create TH3 for activity correction factors

    for (UInt_t v = 0; v < imageIndices.size(); ++v){
        std::array<Int_t, 3> coordinate = imageIndices.at(v);
        Double_t correctionFactor = 0.0;

        Double_t activityInVoxel = A_v->GetBinContent(coordinate[0], coordinate[1], coordinate[2]);
        if (activityInVoxel >= activityThreshold){
            TH3F* p_dcbPrime = (TH3F*)p_dcbvPrime->At(v);

            for (Int_t d = 1; d <= numberOfDetectors; ++d){
                for (Int_t c = 1; c <= numberOfDetectors; ++c){
                    for (Int_t bin = 1; bin <= NbinsProjections; ++bin){

                        Double_t N_dcbPrimeContent = projections->GetBinContent(d, c, bin);
                        if (N_dcbPrimeContent != 0){

                            Double_t p_dcbPrimeContent = p_dcbPrime->GetBinContent(d, c, bin);
                            if (p_dcbPrimeContent != 0){

                                correctionFactor += p_dcbPrimeContent / N_dcbPrimeContent;
                            }
                        }
                    }
                }
            }

            correctionFactor = correctionFactor / p_dcbvSum[v];
            correctionFactor = std::pow(correctionFactor, accelerator);  // accelerate the convergence

        } else{
            correctionFactor = 1.0;
        }

        A_v->SetBinContent(coordinate[0], coordinate[1], coordinate[2],
                           correctionFactor * activityInVoxel);
    }
}

// ##### ROOT #####
void MLEM_Old_Modified(){
    TBenchmark b;
    b.Start("total");

    // Specify the location of the measurement file
    TString pathToMeasurements = "../folder/subfolder/*.root";

    // Specify the location of the projections file
    TString pathToProjection = "../folder/subfolder/*.root";

    // create an instance of the reconstruction class
    ReconstructionMLEM* reco = new ReconstructionMLEM(pathToMeasurements, pathToProjection);
    reco->setAccelerator(1.5);
    reco->setActivityThreshold(0);
    reco->setImageVolume({-30, 30, -30, 30, 10, 20});

    Int_t numberOfIterations = 100;
    Double_t stoppingCriterion = 0.999;
    reco->start(numberOfIterations, stoppingCriterion, kTRUE);

    b.Stop("total");
    std::cout << "\nTotal Time:\t\t" << b.GetRealTime("total") << " seconds\n";
}

// ##### COMPILE FILE #####
#ifndef __CINT__
int main(){
    MLEM_Old_Modified();

    return 0;
}
#endif
