// imagespace.h

#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <TH3.h>

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
