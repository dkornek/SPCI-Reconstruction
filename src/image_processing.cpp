// image_processing.cpp

#include "image_processing.h"

#include <iostream>
#include <TCanvas.h>


// PUBLIC FUNCTIONS
ActivityDistribution::ActivityDistribution(TString fileName){
    // open file to save reconstructed activity

    saveData = new TFile(fileName, "RECREATE");
}

ActivityDistribution::~ActivityDistribution(){
    std::cout << "\nSaved reconstructed activity to file.\n";

    if (saveData->IsOpen()){
        delete saveData;
    }
}

void ActivityDistribution::save3DHistogram(TH3F *activity3D){
    // save the original TH3F activity in the image space

    if (saveData->IsOpen()){

        saveData->cd();

        TCanvas* a3DCanvas = new TCanvas("A_v", "3D Activity Distribution");
        activity3D->Draw("BOX2");
        a3DCanvas->Write();
        delete a3DCanvas;
    }
}

void ActivityDistribution::saveChiStatistics(TGraph *chi){
    // save the Chi²-Statistics

    if (saveData->IsOpen()){

        saveData->cd();

        TCanvas* chiCanvas = new TCanvas("chi", "#chi^{2}-Statistics");
        chiCanvas->SetLogy();
        chiCanvas->SetLogx();
        chi->Draw("AL*");
        chiCanvas->Write();
        delete chiCanvas;
    }
}

void ActivityDistribution::save2DSlices(TH3F *activity3D){
    // save the z-axis slices of the 3D activity distribution

    if (saveData->IsOpen()){

        TDirectory* subDir = saveData->mkdir("Slices", "Z-axis Slices");
        subDir->cd();

        Int_t detectorD = activity3D->GetNbinsX();
        Int_t detectorC = activity3D->GetNbinsY();
        Int_t numberOfSlices = activity3D->GetNbinsZ();

        for (Int_t nSlice = 1; nSlice <= numberOfSlices; ++nSlice){
            TString internalName;
            internalName.Form("A_v3D_Slice_%02i", nSlice);

            TString description;
            description.Form("Slice %02i of final activity distribution", nSlice);

            TH2F* slice = new TH2F(internalName, description,
                                   detectorD, 0, detectorD,
                                   detectorC, 0, detectorC);
            slice->GetXaxis()->SetTitle("#it{x} / mm");
            slice->GetYaxis()->SetTitle("#it{y} / mm");
            slice->GetZaxis()->SetTitle("Relative Activity");
            slice->GetZaxis()->SetTitleOffset(Float_t(1.4));
            slice->SetStats(kFALSE);

            for (Int_t d = 1; d <= detectorD; ++d){
                for (Int_t c = 1; c <= detectorC; ++c){
                    // reconstruct the slice

                    slice->SetBinContent(d, c, activity3D->GetBinContent(d, c, nSlice));
                }
            }

            TCanvas* a2DCanvas = new TCanvas(internalName, description);
            a2DCanvas->SetRightMargin(0.2);
            slice->Draw("COLZ");
            a2DCanvas->Write();
            delete slice;
            delete a2DCanvas;
        }
    }
}

void ActivityDistribution::save2DSteps(TList *steps){
    // save the iteration steps of the reconstruction process

    if (saveData->IsOpen()){

        TDirectory* subDir = saveData->mkdir("Steps", "Reconstruction Steps");
        subDir->cd();

        TIter next(steps);
        TH2F* step;
        while (step = (TH2F*)next()){

            TCanvas* stepCanvas = new TCanvas(step->GetName(), step->GetTitle());
            stepCanvas->SetRightMargin(0.2);
            step->Draw("COLZ");
            stepCanvas->Write();
            delete stepCanvas;
        }

        TDirectory* subDir2 = saveData->mkdir("StepsTH2", "Raw Steps");
        subDir2->cd();

        steps->Write();
    }
}

void ActivityDistribution::saveAll(TH3F *a3D, TList *s, TGraph *chi){
    // save all data

    save3DHistogram(a3D);
    saveChiStatistics(chi);
    save2DSlices(a3D);
    save2DSteps(s);
}
