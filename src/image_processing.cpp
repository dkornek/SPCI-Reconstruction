// image_processing.cpp

#include "image_processing.h"

#include <iostream>
#include <TStyle.h>


// PUBLIC FUNCTIONS
ActivityDistribution::ActivityDistribution(TH3F *a, TString fileName) : activity(a){
    // open file to save reconstructed activity

    this->saveData = new TFile(fileName, "RECREATE");
}

ActivityDistribution::~ActivityDistribution(){
    std::cout << "Saved reconstructed activity to file.\n";
    delete this->projection;
    delete this->saveData;
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
        this->projection = (TH2D*)activity->Project3D("yx");
        this->projection->SetNameTitle("A_v_2D", "2D Projection of Activity Distribution");
        this->projection->SetXTitle("x / Pixel");
        this->projection->SetYTitle("y / Pixel");

        this->projection->Write();
    }
}

//void mlem(){

////    TCanvas* alpha = new TCanvas("a", "3D Activity distribution", 600, 400);
////    reco->A_v->Draw("BOX");
////    alpha->Update();

////    TCanvas* beta = new TCanvas("b", "2D Projection", 600, 400);
////    TH2D* projection = (TH2D*)reco->A_v->Project3D("yx");
////    projection->Draw("COLZ");
////    beta->Update();

////    TCanvas* gammma = new TCanvas("c", "2D Slice", 600, 400);
////    TH2F* slice = new TH2F("2d_slice", "2D Activity distribution",
////                           7, 0, 7,
////                           7, 0, 7);
////    for (Int_t d = 1; d <= 7; ++d){
////        for (Int_t c = 1; c <= 7; ++c){
////            slice->SetBinContent(d, c, reco->A_v->GetBinContent(d, c, 1));
////        }
////    }

////    slice->Draw("COLZ");
////    gamma->Update();

////    delete reco;
//}
