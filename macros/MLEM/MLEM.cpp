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


const Bool_t normalizeSpectra = kTRUE;

// ##### UTILITIES #####
namespace Utilities {

    Bool_t checkForSameNumberOfBins(const Int_t NbinsMeasurements, const Int_t NbinsSystemMatrix){
        // image reconstruction is not possible unless the number of bins are identical

        if (NbinsMeasurements == NbinsSystemMatrix){
            return kTRUE;
        } else {
            std::cout << "Binning pattern is not consistent\n";
            return kFALSE;
        }
    }

    void getDetectorIndices(const TString nameOfSpectrum, Int_t& d, Int_t& c){
        // extract index of both the detectors from the name of the spectrum

        d = ((TString)nameOfSpectrum(0, 2)).Atoi();
        c = ((TString)nameOfSpectrum(2, 2)).Atoi();
    }

    void getImageSpaceIndices(const TString titleOfVoxel, Int_t &x, Int_t &y, Int_t &z){
        // extract location in image space from the title (=name of directory) of the voxel v

        // prepare title
        TString title = titleOfVoxel.Copy();
        title.Remove(0, 1);
        title.Replace(title.Length() - 1, 1, ',');

        // actual indices / coordinates
        std::vector<Int_t> location = { 0, 0, 0 };
        for (Int_t i = 0; i < 3; ++i){
            Int_t positionOfSeparator = title.First(',');
            Int_t lengthOfString = title.Last(',');

            location[i] =  ((TString)title(0, positionOfSeparator)).Atoi();
            title = title(positionOfSeparator + 1, lengthOfString - positionOfSeparator + 1);
        }

        // assign indices / coordinates
        x = location[0];
        y = location[1];
        z = location[2];
    }
}

// ##### MEASUREMENTS #####
class Measurements{
public:
    Measurements(const TString pathToMeasurements);
    ~Measurements();

    void createN_dcb();

    Int_t numberOfDetectors;
    Int_t numberOfBins;

    TH3F* N_dcb = nullptr;  // contains all detected events in a 3D histogram

private:
    void getNumbers();
    void fillN_dcb(const Bool_t normalized);

    TFile* measurementsFile = nullptr;
    TIter* nextS_dcMeasurements = nullptr;
    TKey* keyS_dcMeasurements = nullptr;
};

Measurements::Measurements(const TString pathToMeasurements) : numberOfDetectors(0), numberOfBins(0){
    // open the measurement *.root file

    measurementsFile = new TFile(pathToMeasurements, "READ");
    nextS_dcMeasurements = new TIter(measurementsFile->GetListOfKeys());
    getNumbers();
}

Measurements::~Measurements(){
    delete N_dcb;

    delete keyS_dcMeasurements;
    delete nextS_dcMeasurements;
    delete measurementsFile;
}

void Measurements::createN_dcb(){
    // create empty 3d histogram of N_dcb
    // x = detector d, y = detector c, z = bins in spectrum
    // histogram is based on the last spectrum

    N_dcb = new TH3F("N_dcb", "Number of events in bin b in S_dc",
                     numberOfDetectors, 0, numberOfDetectors,
                     numberOfDetectors, 0, numberOfDetectors,
                     numberOfBins, 0, numberOfBins);

    fillN_dcb(normalizeSpectra);

    TCanvas* canvasN = new TCanvas("N_dcbCanvas", "Measurements", 600, 600);
    N_dcb->Draw("BOX2Z");
    canvasN->Update();
}

void Measurements::getNumbers(){
    // extract number of detectors and bins

    // title of spectra is expected to be of "ddcc" format
    Int_t detD;
    Int_t detC;
    TString nameOfLastS_dc;
    nameOfLastS_dc = measurementsFile->GetListOfKeys()->Last()->GetName();
    Utilities::getDetectorIndices(nameOfLastS_dc, detD, detC);

    // get number of detectors
    numberOfDetectors = std::max(detD, detC) + 1;

    // get number of bins
    TH1F* lastS_dc = (TH1F*)measurementsFile->Get(nameOfLastS_dc);
    numberOfBins = lastS_dc->GetNbinsX();
    delete lastS_dc;
}

void Measurements::fillN_dcb(const Bool_t normalized){
    // fill N_dcb with the number of events extracted from the spectra

    nextS_dcMeasurements->Reset();

    if (normalized){
        while ((keyS_dcMeasurements = (TKey*)nextS_dcMeasurements->Next())){
            // iterate through every measurement spectrum

            Int_t detectorD;
            Int_t detectorC;
            TString nameOfS_dc = keyS_dcMeasurements->GetName();
            Utilities::getDetectorIndices(nameOfS_dc, detectorD, detectorC);

            // get spectrum
            TH1F* S_dc = (TH1F*)measurementsFile->Get(nameOfS_dc);

            Double_t integral = S_dc->Integral();
            if (integral != 0){

                S_dc->Scale(1.0 / integral);
                for (Int_t bin = 1; bin <= numberOfBins; ++bin){
                    N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         S_dc->GetBinContent(bin));
                }

            } else{

                for (Int_t bin = 1; bin <= numberOfBins; ++bin){
                    N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         0);
                }
            }

            delete S_dc;
        }
    }

    else {
        while ((keyS_dcMeasurements = (TKey*)nextS_dcMeasurements->Next())){
            // iterate through every measurement spectrum

            Int_t detectorD;
            Int_t detectorC;
            TString nameOfS_dc = keyS_dcMeasurements->GetName();
            Utilities::getDetectorIndices(nameOfS_dc, detectorD, detectorC);

            // get spectrum
            TH1F* S_dc = (TH1F*)measurementsFile->Get(nameOfS_dc);

            if (S_dc->Integral() != 0){
                for (Int_t bin = 1; bin <= numberOfBins; ++bin){
                    N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                         S_dc->GetBinContent(bin));
                }
            }

            delete S_dc;
        }
    }
}

// ##### SYSTEM MATRIX #####
class SystemMatrix{
public:
    SystemMatrix(const TString pathToProjections);
    ~SystemMatrix();

    void createSystemMatrix(const Bool_t normalized, const TH3F* N_dcb, const Int_t nDet);
    void createP_dcbvPrime(const TH3F* N_dcb);

    Int_t numberOfBins;
    std::array<Int_t, 3> size;

    TList* systemMatrix = nullptr;
    TList* p_dcbvPrime = nullptr;  // systemMatrix x N_dcb
    std::vector<Double_t> sensitivities;

private:
    void getNumbers();

    TFile* systemMatrixFile = nullptr;
    TIter* nextVoxel = nullptr;
    TKey* keyVoxel = nullptr;
    TKey* keyS_dcVoxel = nullptr;
};

SystemMatrix::SystemMatrix(const TString pathToProjections) : numberOfBins(0){
    // open the projections *.root file

    systemMatrixFile = new TFile(pathToProjections, "READ");
    nextVoxel = new TIter(systemMatrixFile->GetListOfKeys());
    getNumbers();

    TString titleOfLastVoxel = systemMatrixFile->GetListOfKeys()->Last()->GetTitle();
    Int_t sizeX, sizeY, sizeZ;
    Utilities::getImageSpaceIndices(titleOfLastVoxel, sizeX, sizeY, sizeZ);
    size = { sizeX, sizeY, sizeZ };
}

SystemMatrix::~SystemMatrix(){
    p_dcbvPrime->Delete();
    systemMatrix->Delete();

    delete p_dcbvPrime;
    delete systemMatrix;
    delete keyS_dcVoxel;
    delete keyVoxel;
    delete nextVoxel;
    delete systemMatrixFile;
}

void SystemMatrix::createSystemMatrix(const Bool_t normalized, const TH3F *N_dcb, const Int_t nDet){
    // create list of 3d histograms containing the probabilities for each voxel = system matrix

    systemMatrix = new TList();

    // iterate through all voxels v
    nextVoxel->Reset();
    while ((keyVoxel = (TKey*)nextVoxel->Next())){
        TString nameOfVoxel = keyVoxel->GetName();

        TString description;
        description.Form("Probabilities to contribute one count to N_dcb from voxel %s", nameOfVoxel.Data());

        TString internalName;
        internalName.Form("p_dcb_%s", nameOfVoxel.Data());
        TH3F* p_dcb = new TH3F(internalName, description,
                               nDet, 0, nDet,
                               nDet, 0, nDet,
                               numberOfBins, 0, numberOfBins);

        TDirectory* dirOfVoxel = (TDirectory*)systemMatrixFile->Get(nameOfVoxel);
        TIter nextS_dc(dirOfVoxel->GetListOfKeys());
        while ((keyS_dcVoxel = (TKey*)nextS_dc())){
            // iterate through all spectra in voxel v

            TString nameOfS_dc = keyS_dcVoxel->GetName();
            TH1F* S_dc = (TH1F*)dirOfVoxel->Get(nameOfS_dc);

            if (normalized){
                S_dc->Scale(1.0 / S_dc->Integral());
            }

            Int_t detectorD;
            Int_t detectorC;
            Utilities::getDetectorIndices(nameOfS_dc(3, 4), detectorD, detectorC);
            for (Int_t bin = 1; bin <= numberOfBins; ++bin){

                Double_t N_dcbBinContent = N_dcb->GetBinContent(detectorD + 1, detectorC + 1, bin);
                if (N_dcbBinContent != 0){
                    p_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin, S_dc->GetBinContent(bin));
                }
            }

            delete S_dc;
        }

        systemMatrix->AddLast(p_dcb);
        sensitivities.push_back(p_dcb->Integral());
        delete dirOfVoxel;
    }
}

void SystemMatrix::createP_dcbvPrime(const TH3F *N_dcb){
    // Multiply p_dcbv * N_dcb to reduce computation time

    p_dcbvPrime = new TList();

    TIter next(systemMatrix);
    TH3F* p_dcb;
    while ((p_dcb = (TH3F*)next())){

        TString name;
        name.Form("%s_prime", p_dcb->GetName());

        TH3F* p_dcbPrime = (TH3F*)p_dcb->Clone(name);
        p_dcbPrime->Multiply(N_dcb);

        p_dcbvPrime->AddLast(p_dcbPrime);
    }
}

void SystemMatrix::getNumbers(){
    // get number of bins

    TString nameOfLastVoxel = systemMatrixFile->GetListOfKeys()->Last()->GetName();
    TDirectory* dirOfLastVoxel = (TDirectory*)systemMatrixFile->Get(nameOfLastVoxel);
    TString nameOfLastS_dc = dirOfLastVoxel->GetListOfKeys()->Last()->GetName();
    numberOfBins = ((TH1F*)dirOfLastVoxel->Get(nameOfLastS_dc))->GetNbinsX();
    delete dirOfLastVoxel;
}

// ##### IMAGE SPACE / ACTIVITY #####
class ImageSpace{
public:
    ImageSpace(const std::vector<Double_t> volume, const std::array<Int_t, 3> size);
    ~ImageSpace();

    void createA_v();
    void makeA_vHomogeneous();

    // Image space info
    std::vector<Double_t> imageVolume;                      // in mm, xMin, xMax, yMin, yMax, zMin, zMax
    std::vector<std::array<Int_t, 3> > imageIndices;        // coordinates of voxel v

    // Activity distribution
    Int_t numberOfBinsX;
    Int_t numberOfBinsY;
    Int_t numberOfBinsZ;

    Int_t numberOfVoxels;

    TH3F* A_v = nullptr;

private:
    void setImageVolume(const std::vector<Double_t> volume){imageVolume = volume;}
    void setImageIndices();
};

ImageSpace::ImageSpace(const std::vector<Double_t> volume, const std::array<Int_t, 3> size){
    setImageVolume(volume);

    numberOfBinsX = size[0] + 1;
    numberOfBinsY = size[1] + 1;
    numberOfBinsZ = size[2] + 1;

    setImageIndices();
    numberOfVoxels = imageIndices.size();

    createA_v();
}

ImageSpace::~ImageSpace(){
    delete A_v;
}

void ImageSpace::createA_v(){
    // create the activity distribution contained in the image space

    A_v = new TH3F("A_v", "3D Activity Distribution",
                   numberOfBinsX, imageVolume[0], imageVolume[1],
                   numberOfBinsY, imageVolume[2], imageVolume[3],
                   numberOfBinsZ, imageVolume[4], imageVolume[5]);

    // Plotting options
    A_v->GetXaxis()->SetTitle("#it{x} / mm");
    A_v->GetXaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * (numberOfBinsX - 1)), "X");

    A_v->GetYaxis()->SetTitle("#it{y} / mm");
    A_v->GetXaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * (numberOfBinsY - 1)), "Y");

    A_v->GetZaxis()->SetTitle("#it{z} / mm");
    A_v->GetZaxis()->SetTitleOffset(1.5);
    A_v->SetNdivisions(std::max(5.0, 0.5 * (numberOfBinsZ - 1)), "Z");
}

void ImageSpace::makeA_vHomogeneous(){
   // set all cells to 1

    for (UInt_t v = 0; v < imageIndices.size(); ++v){
        // iterate through every voxel

        std::array<Int_t, 3> coordinate = imageIndices.at(v);
        A_v->SetBinContent(coordinate[0], coordinate[1], coordinate[2], 1.0);
    }
}

void ImageSpace::setImageIndices(){
    // assign coordinates to voxel number v

    Int_t x = numberOfBinsX - 1;
    Int_t y = numberOfBinsY - 1;
    Int_t z = numberOfBinsZ - 1;

    for (Int_t zIndex = 0; zIndex <= z; ++zIndex){
        for (Int_t yIndex = 0; yIndex <= y; ++yIndex){
            for (Int_t xIndex = 0; xIndex <= x; ++xIndex){

                std::array<Int_t, 3> coordinate = {xIndex + 1, yIndex + 1, zIndex + 1};
                imageIndices.push_back(coordinate);
            }
        }
    }
}

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

// ##### ROOT #####
void MLEM(){
    TBenchmark b;
    b.Start("total");

    // Specify the location of the measurement file
    TString pathToMeasurements = "../folder/subfolder/*.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos1.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos24.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos48.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos6.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourcePos0+24+26+36.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourceCross.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase49/Bins50/SourceDiagonal.root";

//    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins50/SourceSquare.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins50/Source0.root";
    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins50/Source20.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins50/SourceH.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins700/SourceHZDR.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins700/SourceH.root";

    // Specify the location of the projections file
    TString pathToProjection = "../folder/subfolder/*.root";
    pathToProjection = "../../data/SystemMatrix/Original/SPCIBase49_120.root";
//    pathToProjection = "../../data/SystemMatrix/Bins50/SPCIBase441.root";
//    pathToProjection = "../../data/SystemMatrix/Original/SPCIBase441.root";

    // create an instance of the reconstruction class
    ReconstructionMLEM* reco = new ReconstructionMLEM(pathToMeasurements,
                                                      pathToProjection,
                                                      {-30, 30, -30, 30, 5, 15});  // {-50, 50, -50, 50, 5, 10}
    reco->setAccelerator(1.5);
    reco->start(50, 0.999);

    b.Stop("total");
    std::cout << "\nTotal Time:\t\t" << b.GetRealTime("total") << " seconds\n";
}

// ##### COMPILE FILE #####
#ifndef __CINT__
int main(){
    MLEM();

    return 0;
}
#endif
