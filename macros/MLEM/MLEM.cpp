// ROOT Macro using the Maximum Likelihood Expectation Maximization algorithm for SPCI-Reconstruction

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-05-06


#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2.h>
#include <TH3.h>
#include <TBenchmark.h>
#include <TCanvas.h>


// ##### RECONSTRUCTION #####
class ReconstructionMLEM{
public:
    ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections);
    ~ReconstructionMLEM();
    void start(Int_t maxNumberOfIterations, Bool_t doPlotting);

    // Activity distribution
    TH3F* A_v;

private:
    // ##### UTILITIES #####
    void openFile(TString fileType, TString pathToFile);
    void getDetectorIndices(TString nameOfSpectrum, Int_t& d, Int_t& c);
    void getImageSpaceIndices(TString titleOfVoxel, Int_t& x, Int_t& y, Int_t& z);
    void fillN_dcb();
    void makeA_vHomogeneous();
    void checkValidity();

    // ##### PREPARATION FUNCTIONS #####
    void createN_dcb();
    void createP_dcbv();
    void createA_v();

    // ##### CALCULATION FUNCTIONS #####
    void calculate();
    void projection(TH3F& p);
    void backprojection(TH3F& p, TH3F& bp);

    // ##### MEMBERS #####

    TCanvas* canvas;
    Bool_t isCalculationValid;

    // Measurement Data
    TFile* measurementsFile;
    TH3F* N_dcb;

    Int_t numberOfDetectors;
    Int_t NbinsMeasurements;

    TIter* nextS_dcMeasurements;
    TKey* keyS_dcMeasurements;

    // Projection Data
    // Probabilities
    TFile* projectionsFile;
    TList* p_dcbv;
    TList* p_dcbvPrime;  // p_dcbv * N_dcb
    std::vector<Double_t> p_dcbvSum;

    Int_t NbinsProjections;

    TIter* nextVoxel;
    TKey* keyVoxel;
    TKey* keyS_dcVoxel;
};

// ##### PUBLIC FUNCTIONS #####
ReconstructionMLEM::ReconstructionMLEM(TString pathToMeasurements, TString pathToProjections){
    // prepare all data to start the reconstruction

    // open files
    this->openFile("MEASUREMENT", pathToMeasurements);
    this->nextS_dcMeasurements = new TIter(this->measurementsFile->GetListOfKeys());

    this->openFile("PROJECTION", pathToProjections);
    this->nextVoxel = new TIter(this->projectionsFile->GetListOfKeys());

    TBenchmark b;
    // prepare the measurement data N_dcb
    // N_dcb = number of events in detector pair (d, c) separated by bins
    b.Start("N_dcb");
    this->createN_dcb();
    b.Stop("N_dcb");
    std::cout << "\nN_dcb Creation Time: " << b.GetRealTime("N_dcb") << " seconds\n";

    // prepare the probabilities p_dcbv
    // p_dcbv = voxel v specific probabilitiy to contribute one count to N_dcb;
    // there are v 3d histograms p_dcb in the list p_dcbv
    b.Start("p_dcb");
    this->createP_dcbv();
    b.Stop("p_dcb");
    std::cout << "\np_dcb Creation Time: " << b.GetRealTime("p_dcb") << " seconds\n";

    // prepare the homogeneous activity distribution A_v
    // A_v = activity in voxel v of the image space
    b.Start("A_v");
    this->createA_v();
    b.Stop("A_v");
    std::cout << "\nA_v Creation Time: " << b.GetRealTime("A_v") << " seconds\n";

    // show the image reconstruction process
    this->canvas = new TCanvas("activity", "Reconstruction", 10, 10, 1000, 450);
    this->canvas->Divide(2, 1);
}

ReconstructionMLEM::~ReconstructionMLEM(){
    delete this->A_v;

    delete this->keyS_dcMeasurements;
    delete this->nextS_dcMeasurements;
    delete this->N_dcb;
    delete this->measurementsFile;

    delete this->keyS_dcVoxel;
    delete this->keyVoxel;
    delete this->nextVoxel;
    p_dcbv->Delete();
    delete this->p_dcbv;
    delete this->projectionsFile;
}

void ReconstructionMLEM::start(Int_t maxNumberOfIterations, Bool_t doPlotting){
    // execute the calculation

    this->checkValidity();
    if (!this->isCalculationValid){
        return;
    }

    TBenchmark b;
    b.Start("stats");
    for (Int_t i = 0; i < maxNumberOfIterations ; ++i){

        this->calculate();

        if (doPlotting && (i % 10 == 0)){
            this->canvas->cd(1);
            this->A_v->Draw("BOX2Z");

            this->canvas->cd(2);
            TH2D* projection = (TH2D*)this->A_v->Project3D("yx");
            projection->Draw("COLZ");

            this->canvas->Update();
            delete projection;
        }
    }
    b.Stop("stats");

    // Inform user
    std::cout << "Image reconstruction done. Steps: " << maxNumberOfIterations << "\n";
    std::cout << "\nCalculation Time: " << b.GetRealTime("stats") << " seconds\n";

    this->canvas->cd(1);
    this->A_v->Draw("BOX2Z");

    this->canvas->cd(2);
    TH2D* projection = (TH2D*)this->A_v->Project3D("yx");
    projection->Draw("COLZ");
    this->canvas->Update();
}

// ##### PRIVATE FUNCTIONS #####
// ##### UTILITIES #####
void ReconstructionMLEM::openFile(TString fileType, TString pathToFile){
    // open file of specified type

    if (fileType == "MEASUREMENT"){
        this->measurementsFile = new TFile(pathToFile, "READ");
        return;
    }

    if (fileType == "PROJECTION"){
        this->projectionsFile = new TFile(pathToFile, "READ");
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
    std::vector<Double_t> location = { 0, 0, 0 };
    for (Int_t i = 0; i < 3; ++i){
        Int_t positionOfSeparator = title.First(',');
        Int_t lengthOfString = title.Last(',');

        location[i] =  ((TString)title(0, positionOfSeparator)).Atof();
        title = title(positionOfSeparator + 1, lengthOfString - positionOfSeparator + 1);
    }

    x = location[0];
    y = location[1];
    z = location[2];
}

void ReconstructionMLEM::fillN_dcb(){
    // fill N_dcb with number of events extracted from the spectra

    this->nextS_dcMeasurements->Reset();
    while ((this->keyS_dcMeasurements = (TKey*)this->nextS_dcMeasurements->Next())){

        Int_t detectorD;
        Int_t detectorC;
        TString nameOfS_dc = this->keyS_dcMeasurements->GetName();
        this->getDetectorIndices(nameOfS_dc, detectorD, detectorC);

        // get spectrum
        TH1F* S_dc = (TH1F*)this->measurementsFile->Get(nameOfS_dc);
        S_dc->Scale(1.0 / S_dc->Integral());  // faster convergence

        // fill N_dcb
        for (Int_t bin = 1; bin <= this->NbinsMeasurements; ++bin){
            this->N_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                       S_dc->GetBinContent(bin));
        }

        delete S_dc;
    }
}

void ReconstructionMLEM::makeA_vHomogeneous(){
    // set all bin contents in A_v = 1

    // iterate through all voxels
    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){
        Int_t indexX;
        Int_t indexY;
        Int_t indexZ;

        TString titleOfVoxel = this->keyVoxel->GetTitle();
        this->getImageSpaceIndices(titleOfVoxel, indexX, indexY, indexZ);
        this->A_v->SetBinContent(indexX + 1, indexY + 1, indexZ + 1,
                                 1.0);
    }
}

void ReconstructionMLEM::checkValidity(){
    // calculation is only possible if binning pattern in the measurement and
    // projection spectra is identical

    if (this->NbinsMeasurements == this->NbinsProjections){
        this->isCalculationValid = kTRUE;
    } else{
        std::cout << "Calculation is not valid.\nBinning pattern is not consistent\n";
        this->isCalculationValid = kFALSE;
    }
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

    nameOfLastSpectrum = this->measurementsFile->GetListOfKeys()->Last()->GetName();
    this->getDetectorIndices(nameOfLastSpectrum, detectorD, detectorC);

    // get number of detectors
    this->numberOfDetectors = std::max(detectorD, detectorC) + 1;

    // get number of bins
    this->NbinsMeasurements = ((TH1F*)this->measurementsFile->Get(nameOfLastSpectrum))->GetNbinsX();

    // create the 3d histogram
    this->N_dcb = new TH3F("N_dcb", "Number of events in bin b in S_dc",
                           this->numberOfDetectors, 0, this->numberOfDetectors,
                           this->numberOfDetectors, 0, this->numberOfDetectors,
                           this->NbinsMeasurements, 0, this->NbinsMeasurements);
    this->fillN_dcb();
}

void ReconstructionMLEM::createP_dcbv(){
    // create list of 3d histograms containing the probabilities for each voxel

    this->p_dcbv = new TList();
    this->p_dcbvPrime = new TList();

    // get number of bins
    TString nameOfLastVoxel = this->projectionsFile->GetListOfKeys()->Last()->GetName();
    TDirectory* dirOfLastVoxel = (TDirectory*)this->projectionsFile->Get(nameOfLastVoxel);
    TString nameOfLastS_dc = dirOfLastVoxel->GetListOfKeys()->Last()->GetName();
    this->NbinsProjections = ((TH1F*)dirOfLastVoxel->Get(nameOfLastS_dc))->GetNbinsX();
    delete dirOfLastVoxel;

    // iterate through all voxels v
    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){
        TString nameOfVoxel = this->keyVoxel->GetName();

        TString internalName;
        internalName.Form("p_dcb_%s", nameOfVoxel.Data());

        TString description;
        description.Form("Probabilities to contribute one count to N_dcb from voxel %s", nameOfVoxel.Data());

        TH3F* p_dcb = new TH3F(internalName, description,
                               this->numberOfDetectors, 0, this->numberOfDetectors,
                               this->numberOfDetectors, 0, this->numberOfDetectors,
                               this->NbinsProjections, 0, this->NbinsProjections);

        TDirectory* dirOfVoxel = (TDirectory*)this->projectionsFile->Get(nameOfVoxel);
        TIter nextS_dc(dirOfVoxel->GetListOfKeys());
        while ((this->keyS_dcVoxel = (TKey*)nextS_dc())){
            // iterate through all spectra in voxel v
            TString nameOfS_dc = this->keyS_dcVoxel->GetName();
            TH1F* S_dc = (TH1F*)dirOfVoxel->Get(nameOfS_dc);
            S_dc->Scale(1.0 / S_dc->Integral());  // faster convergence

            Int_t detectorD;
            Int_t detectorC;
            this->getDetectorIndices(nameOfS_dc(3, 4), detectorD, detectorC);
            for (Int_t bin = 1; bin <= this->NbinsProjections; ++bin){
                p_dcb->SetBinContent(detectorD + 1, detectorC + 1, bin,
                                     S_dc->GetBinContent(bin));
            }

            delete S_dc;
        }

        this->p_dcbv->AddLast(p_dcb);
        this->p_dcbvSum.push_back(p_dcb->Integral());

        TH3F* p_dcbPrime = (TH3F*)p_dcb->Clone("Prime");
        p_dcbPrime->Multiply(this->N_dcb);
        this->p_dcbvPrime->AddLast(p_dcbPrime);

        delete dirOfVoxel;
    }
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
    titleOfLastVoxel = this->projectionsFile->GetListOfKeys()->Last()->GetTitle();
    this->getImageSpaceIndices(titleOfLastVoxel, indexX, indexY, indexZ);

    // create the 3d histogram
    this->A_v = new TH3F("A_v", "Activity distribution",
                         indexX + 1, 0, indexX + 1,
                         indexY + 1, 0, indexY + 1,
                         indexZ + 1, 0, indexZ + 1);
    this->makeA_vHomogeneous();
}

// ##### CALCULATION FUNCTIONS #####
void ReconstructionMLEM::calculate(){
    // execute the maximum likelihood expectation maximization algorithm
    // to calculate the activity distribution A (=A_v)

    // calculate the projection
    TH3F* projections = new TH3F("projections", "Projections",
                                 this->numberOfDetectors, 0, this->numberOfDetectors,
                                 this->numberOfDetectors, 0, this->numberOfDetectors,
                                 this->NbinsProjections, 0, this->NbinsProjections);
    this->projection(*projections);

    // calculate the backprojection
    TH3F* backprojections = new TH3F("backprojections", "Backprojections",
                                     this->A_v->GetNbinsX(), 0, this->A_v->GetNbinsX(),
                                     this->A_v->GetNbinsY(), 0, this->A_v->GetNbinsY(),
                                     this->A_v->GetNbinsZ(), 0, this->A_v->GetNbinsZ());
    this->backprojection(*projections, *backprojections);

    // calculate new activity
    this->A_v->Multiply(backprojections);

//    // Chi²-Test
//    TH3F* Chi2 = (TH3F*)this->N_dcb->Clone("Chi");
//    Chi2->Add(projections, -1.0);
//    Chi2->Multiply(Chi2);
//    Chi2->Divide(projections);
//    Double_t testVariable = Chi2->Integral();
//    std::cout << testVariable << "\n";

    delete backprojections;
    delete projections;
}

void ReconstructionMLEM::projection(TH3F& p){

    Int_t voxel = 0;
    for (Int_t z = 1; z <= this->A_v->GetNbinsZ(); ++z){
        for (Int_t y = 1; y <= this->A_v->GetNbinsY(); ++y){
            for (Int_t x = 1; x <= this->A_v->GetNbinsX(); ++x){

                TH3F* p_dcb = (TH3F*)this->p_dcbv->At(voxel);
                Double_t activityInVoxel = this->A_v->GetBinContent(x, y, z);
                p.Add(p_dcb, activityInVoxel);

                ++voxel;
            }
        }
    }
}

void ReconstructionMLEM::backprojection(TH3F& p, TH3F& bp){
    // create TH3 for activity correction factors

    Int_t voxel = 0;
    Double_t lowerThreshold = 0.0001;  // heavily dependent of the real activity distribution
                                       // needs adjustment

    this->nextVoxel->Reset();
    while ((this->keyVoxel = (TKey*)this->nextVoxel->Next())){

        Int_t indexX;
        Int_t indexY;
        Int_t indexZ;
        this->getImageSpaceIndices(this->keyVoxel->GetTitle(), indexX, indexY, indexZ);

        Double_t correctionFactor = 0.0;
        Double_t activityInVoxel = this->A_v->GetBinContent(indexX + 1, indexY + 1, indexZ + 1);

        // only iterate over voxels with enough activity content --> save computing time
        if (activityInVoxel >= lowerThreshold) {
            TH3F* p_dcbPrime = (TH3F*)this->p_dcbvPrime->At(voxel);

            for (Int_t d = 1; d <= this->numberOfDetectors; ++d){
                for (Int_t c = 1; c <= this->numberOfDetectors; ++c){
                    for (Int_t bin = 1; bin <= this->NbinsProjections; ++bin){

                        Double_t N_dcbPrime = p.GetBinContent(d, c, bin);
                        if (N_dcbPrime != 0) {
                            correctionFactor += p_dcbPrime->GetBinContent(d, c, bin) / N_dcbPrime;
                        }
                    }
                }
            }

            correctionFactor = correctionFactor / this->p_dcbvSum[voxel];
            correctionFactor = std::pow(correctionFactor, 1.3);  // accelerate the convergence

        } else{
            correctionFactor = 1.0;
        }

        bp.SetBinContent(indexX + 1, indexY + 1, indexZ + 1,
                        correctionFactor);

        ++voxel;
    }
}

// ##### ROOT #####
void MLEM(){
    TBenchmark bench;
    bench.Start("total");

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

    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins50/SourceSquare.root";
//    pathToMeasurements = "../../data/Measurements/SPCIBase441/Bins50/Source0.root";

    // Specify the location of the projections file
    TString pathToProjection = "../folder/subfolder/*.root";
//    pathToProjection = "../../data/SystemMatrix/Bins50/SPCIBase49.root";
    pathToProjection = "../../data/SystemMatrix/Bins50/SPCIBase441.root";

    Int_t numberOfIterations = 100;

    // create an instance of the reconstruction class
    ReconstructionMLEM* reco = new ReconstructionMLEM(pathToMeasurements, pathToProjection);
    reco->start(numberOfIterations, kFALSE);

    bench.Stop("total");
    std::cout << "\nTotal Time: " << bench.GetRealTime("total") << " seconds\n";
}

// ##### COMPILE FILE #####
#ifndef __CINT__
int main(){
    MLEM();
    return 0;
}
#endif
