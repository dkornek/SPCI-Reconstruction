// systemmatrix.cpp

#include "systemmatrix.h"
#include "utilities.h"

// ##### SYSTEM MATRIX #####
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
    // p_dcbvPrime->Clear();
    // systemMatrix->Delete();

    delete p_dcbvPrime;
    delete systemMatrix;
    delete keyS_dcVoxel;
    delete keyVoxel;
    delete nextVoxel;
    delete systemMatrixFile;
}

void SystemMatrix::createSystemMatrix(const Int_t nDet){
    // (OE MODE) create 4d vector containing the probabilities for each voxel = system matrix

    // iterate through all voxels v
    nextVoxel->Reset();
    while ((keyVoxel = (TKey*)nextVoxel->Next())){
        std::vector<std::vector<std::vector<Double_t> > > p_dcb;

        for (Int_t d = 0; d < nDet; ++d){
            std::vector<std::vector<Double_t> > p_cb;

            for (Int_t c = 0; c < nDet; ++c){
                std::vector<Double_t> p_b;
                p_cb.push_back(p_b);
            }

            p_dcb.push_back(p_cb);
        }

        Double_t sensitivity = 0;
        TString nameOfVoxel = keyVoxel->GetName();
        TDirectory* dirOfVoxel = (TDirectory*)systemMatrixFile->Get(nameOfVoxel);
        TIter nextS_dc(dirOfVoxel->GetListOfKeys());
        while ((keyS_dcVoxel = (TKey*)nextS_dc())){
            // iterate through all spectra in voxel v

            TString nameOfS_dc = keyS_dcVoxel->GetName();
            TH1F* S_dc = (TH1F*)dirOfVoxel->Get(nameOfS_dc);

            // Divide counts by numbers of EMISSIONS (usually unkown) in the voxel to maintain absolute counts
            S_dc->Scale(1.0 / 5000000);  // Tonis Simulation 20181212_5x5mm2.root with 441 positions, solid angle correction not considered

            Int_t detectorD;
            Int_t detectorC;
            Utilities::getDetectorIndices(nameOfS_dc(3, 4), detectorD, detectorC);
            for (Int_t bin = 1; bin <= numberOfBins; ++bin){
                Double_t binContent = S_dc->GetBinContent(bin);
                p_dcb.at(detectorD).at(detectorC).push_back(binContent);

                sensitivity += binContent;
            }

            delete S_dc;
        }

        delete dirOfVoxel;
        sensitivities.push_back(sensitivity);
        systemMatrixVector.push_back(p_dcb);
    }
}

void SystemMatrix::createSystemMatrix(const Bool_t normalized, const TH3F *N_dcb, const Int_t nDet){
    // (ML-EM MODE) create list of 3d histograms containing the probabilities for each voxel = system matrix

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

            // Divide counts by numbers of EMISSIONS (usually unkown) in the voxel to maintain absolute counts
            S_dc->Scale(1.0 / 5000000);  // Tonis Simulation 20181212_5x5mm2.root with 441 positions, solid angle correction not considered

            if (normalized){
                // faster converge because distant detectors, which have more precise position information, get more weight
                // only working for MLEM

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
    // Multiply p_dcbv * N_dcb once to reduce computation time during iterations

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
