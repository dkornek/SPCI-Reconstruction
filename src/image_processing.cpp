// image_processing.cpp

#include "image_processing.h"

#include <iostream>
//#include <TStyle.h>


// PUBLIC FUNCTIONS
ActivityDistribution::ActivityDistribution(TH3F *a, TString fileName) : activity(a){
    // open file to save reconstructed activity

    this->saveData = new TFile(fileName, "RECREATE");
}

ActivityDistribution::~ActivityDistribution(){
    std::cout << "\nSaved reconstructed activity to file.\n";

    if (this->saveData->IsOpen()){
        delete this->saveData;
    }
}

void ActivityDistribution::save3DHistogram(){
    // save the original TH3F activity in the image space

    if (this->saveData->IsOpen()){
        this->activity->SetDrawOption("BOX");

        this->activity->SetNameTitle("A_v_3D", "3D Activity Distribution");

        this->activity->SetXTitle("x / Voxel");
        this->activity->SetNdivisions(this->activity->GetNbinsX(), "X");

        this->activity->SetYTitle("y / Voxel");
        this->activity->SetNdivisions(this->activity->GetNbinsY(), "Y");

        this->activity->SetZTitle("z / Voxel");
        this->activity->SetNdivisions(this->activity->GetNbinsZ(), "Z");

        this->activity->Write();
    }
}

void ActivityDistribution::save2DProjection(){
    // save the 2D projection of the 3D activity distribution

    if (this->saveData->IsOpen()){
        TH2D* projection = (TH2D*)activity->Project3D("yx");

        projection->SetDrawOption("COLZ");

        projection->SetNameTitle("A_v_2D", "2D Projection of Activity Distribution");

        projection->SetXTitle("x / Pixel");
        projection->SetNdivisions(projection->GetNbinsX(), "X");

        projection->SetYTitle("y / Pixel");
        projection->SetNdivisions(projection->GetNbinsY(), "Y");

        projection->Write();
        delete projection;
    }
}

void ActivityDistribution::save2DSlices(){
    // save the z-Slices of the 3D activity distribtion

    if (this->saveData->IsOpen()){
        Int_t detectorD = this->activity->GetNbinsX();
        Int_t detectorC = this->activity->GetNbinsY();
        Int_t numberOfSlices = this->activity->GetNbinsZ();

        TDirectory* subDir = this->saveData->mkdir("Slices", "yx Slices of Activity Distribution");
        subDir->cd();

        TH2F* slice;
        for (Int_t nSlice = 1; nSlice <= numberOfSlices; ++nSlice){
            // iterate through every slice
            TString internalName;
            TString description;

            internalName.Form("Slice_%03i", nSlice);
            description.Form("Slice_%03i of the 3D Activity Distribution", nSlice);
            slice = new TH2F(internalName, description,
                             detectorD, 0, detectorD,
                             detectorC, 0, detectorC);

            for (Int_t d = 1; d <= detectorD; ++d){
                for (Int_t c = 1; c <= detectorC; ++c){
                    // reconstruct the slice

                    slice->SetBinContent(d, c, this->activity->GetBinContent(d, c, nSlice));
                }
            }

            slice->SetDrawOption("COLZ");

            slice->SetXTitle("x / Pixel");
            slice->SetNdivisions(detectorD, "X");

            slice->SetYTitle("y / Pixel");
            slice->SetNdivisions(detectorC, "Y");

            slice->Write();
        }

        delete slice;
    }
}
