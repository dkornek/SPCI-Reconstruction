// imagespace.cpp

#include "imagespace.h"

// ##### IMAGE SPACE / ACTIVITY #####
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

    A_v = new TH3F("A_v", "3D Emission Density",
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
