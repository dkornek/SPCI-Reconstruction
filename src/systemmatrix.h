// systemmatrix.h

#pragma once
#include <TH3.h>
#include <TKey.h>
#include <TFile.h>

// ##### SYSTEM MATRIX #####
class SystemMatrix{
public:
    SystemMatrix(const TString pathToProjections);
    ~SystemMatrix();

    void createSystemMatrix(const Int_t nDet);  // OE
    void createSystemMatrix(const Bool_t normalized, const TH3F* N_dcb, const Int_t nDet);  // MLEM
    void createP_dcbvPrime(const TH3F* N_dcb);

    Int_t numberOfBins;
    std::array<Int_t, 3> size;

    TList* systemMatrix = nullptr;
    TList* p_dcbvPrime = nullptr;  // systemMatrix x N_dcb
    std::vector<Double_t> sensitivities;

    // vector mode for system matrix: 1) voxel 2) d 3) c 4) bin
    std::vector<std::vector<std::vector<std::vector<Double_t> > > > systemMatrixVector;

private:
    void getNumbers();

    TFile* systemMatrixFile = nullptr;
    TIter* nextVoxel = nullptr;
    TKey* keyVoxel = nullptr;
    TKey* keyS_dcVoxel = nullptr;
};
