#ifndef IMAGERECONSTRUCTION_H
#define IMAGERECONSTRUCTION_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TStyle.h"
#include "THistPainter.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TBenchmark.h"


using namespace std;

char *StrHist = new char[50];


class ImageReconstruction {

public:
    ImageReconstruction(UInt_t dim,
                        UInt_t voxel, Double_t voxel_len,
                        UInt_t det_n, UInt_t det_m,
                        UInt_t bins, UInt_t rebin = 1);
    ~ImageReconstruction() {delete mVoxelIndex; delete mA_v; delete sum_v; delete mN_dcb;}

    void        SetVoxelID();

    void        SetN_dcbToZero(vector<vector<vector<Double_t> > >*);
    void        SetProbN_dcb(UInt_t voxel);
    void        SetMeasN_dcb();

    void        ResetA_v_alike(vector<vector<vector<Double_t> > >*);
    void        CalculateA_v(UInt_t iterations = 1);
    Double_t    CalculateCorrFactorA_v(UInt_t voxel);
    void        CalculateSum_v();
        
    void        SetInputFileName(TString StrFileName)  { glStrFileIn     = StrFileName;}
    void        SetOutputFileName(TString StrFileName) { glStrFileOut    = StrFileName;}
    void        SetBaseFileName(TString StrFileName)   { glStrBaseFile   = StrFileName;}
        
    TString     GetInputFileName()  {return glStrFileIn;}
    TString     GetOutputFileName() {return glStrFileOut;}
    TString     GetBaseFileName()   {return glStrBaseFile;}

    TCanvas*    ShowA_v();
    TCanvas*    ShowN_dcbReconstructed();
		
private:
    TCanvas*    CanvasA_v;
    TCanvas*    CanvasN_dcbReconstructed;
    TH2D*       hist2d;

    map<UInt_t, map<UInt_t, TH1D*>>     histReconstructed;  //~ TH1* histReconstructed[16][16];
    map<UInt_t, map<UInt_t, TH1*>>      hist_voxel;  //~ TH1* hist_voxel[16][16];

    UInt_t      mDim;
    UInt_t      mNx;
    UInt_t      mNy;
    UInt_t      mNz;
    UInt_t      mVoxel;
    Double_t    mVoxelLen; // in mm
    UInt_t      mDet_n;
    UInt_t      mDet_m;
    UInt_t      mNDet;
    UInt_t      mBins;
    UInt_t      mRebin;
    Double_t    mIntegral2d;
    Double_t    mEBinsMax;
    Double_t    mEBinsMin;

    TString     glStrBaseFile, glStrFileIn, glStrFileOut;

    vector<vector<vector<Int_t> > >*                mVoxelIndex;
    vector<vector<vector<Double_t> > >*             mA_v;
    vector<vector<vector<Double_t> > >*             sum_v;
    vector<vector<vector<Double_t> > >*             mN_dcb;
    vector<vector<vector<vector<Double_t> > >* >    mVoxelN_dcb;
};

ImageReconstruction::ImageReconstruction(UInt_t dim,
                                         UInt_t voxel, Double_t voxel_len,
                                         UInt_t det_n, UInt_t det_m,
                                         UInt_t bins, UInt_t rebin) :
    CanvasA_v(NULL), CanvasN_dcbReconstructed(NULL),
    hist2d(NULL),
    mDim(dim),
    mVoxel(voxel), mVoxelLen(voxel_len),
    mDet_n(det_n), mDet_m(det_m), mNDet(det_n * det_m),
    mBins(bins/rebin), mRebin(rebin),
    mIntegral2d(0), mEBinsMax(0) {

    TBenchmark bench;

    // Initilize image dimensions
    mNz = pow(mVoxel, 1./mDim);
    mNy = 1;
    mNx = 1;

    if (mDim > 1)
        mNy = pow(mVoxel, 1./mDim);

    if (mDim > 2)
        mNx = pow(mVoxel, 1./mDim);

    // Initilize probabilities
    bench.Start("p_dcb");
    for (UInt_t iVoxel = 0; iVoxel < voxel; iVoxel++) {
        // Initialize with nullptr
        mN_dcb = NULL;

        //Initialize with free address
        mN_dcb = new vector<vector<vector<Double_t> > >();

        //Initilize matrix elements with 0
        SetN_dcbToZero(mN_dcb);

        //Initilize matrix with values
        SetProbN_dcb(iVoxel);

        //Add address to voxel vector
        mVoxelN_dcb.push_back(mN_dcb);
    }
    mN_dcb = NULL;
    bench.Stop("p_dcb");
    cout << "\np_dcb Creation Time: " << bench.GetRealTime("p_dcb") << " seconds\n";
	
    // Initilize mVoxelIndex matrix elements with voxel ID
    mVoxelIndex = NULL;
    mVoxelIndex = new vector<vector<vector<Int_t> > >();
    SetVoxelID();

    //Initilize A_v matrix elements with 1
    bench.Start("A_v");
    mA_v = NULL;
    mA_v = new vector<vector<vector<Double_t> > >();
    ResetA_v_alike(mA_v);
    bench.Stop("A_v");
    cout << "\nA_v Creation Time: " << bench.GetRealTime("A_v") << " seconds\n";

    // Initilize N_dcb matrix
    mN_dcb = new vector<vector<vector<Double_t> > >();
    SetN_dcbToZero(mN_dcb);

    //Initilize sum_v matrix with 0
    sum_v = NULL;
    sum_v = new vector<vector<vector<Double_t> > >();
    SetN_dcbToZero(sum_v);
}

void ImageReconstruction::SetVoxelID() {
    ///   z
    ///   /2/4/
    ///  /1/3/ y
    /// |1|3|
    /// |5|7|
    /// x

    UInt_t ID = 0;
    for (UInt_t iNx = 0; iNx < mNx; iNx++) {
        mVoxelIndex->push_back(vector<vector<Int_t> >());

        for (UInt_t iNy = 0; iNy < mNy; iNy++) {
            mVoxelIndex->at(iNx).push_back(vector<Int_t>());

            for (UInt_t iNz = 0; iNz < mNz; iNz++) {
                mVoxelIndex->at(iNx).at(iNy).push_back(ID);
                ID++;
            }
        }
    }
}

void ImageReconstruction::ResetA_v_alike(vector<vector<vector<Double_t> > >* matrix) {
    // creates homogeneous activity distribution (cell = 1)

    matrix->clear();

    for (UInt_t iNx = 0; iNx < mNx; iNx++) {
        matrix->push_back(vector<vector<Double_t> >());

        for (UInt_t iNy = 0; iNy < mNy; iNy++) {
            matrix->at(iNx).push_back(vector<Double_t>());

            for (UInt_t iNz = 0; iNz < mNz; iNz++) {
                matrix->at(iNx).at(iNy).push_back(1);
            }
        }
    }
}

void ImageReconstruction::SetN_dcbToZero(vector<vector<vector<Double_t> > >* matrix) {
    // Delete whole matrix, then build it up again to desired size and set every cell = 0

    matrix->clear();

    for (UInt_t iN = 0; iN < mNDet; iN++) {
        matrix->push_back(vector<vector<Double_t> >());

        for (UInt_t iM = 0; iM < mNDet; iM++) {
            matrix->at(iN).push_back(vector<Double_t>());

            for (UInt_t iBins = 0; iBins < mBins; iBins++) {
                matrix->at(iN).at(iM).push_back(0);
            }
        }
    }
}

void ImageReconstruction::SetProbN_dcb(UInt_t voxel) {

    // TString StrFileIn = GetBaseFileName();
    TFile* file_voxel = new TFile("../../data/SystemMatrix/Bins50/SPCIBase441.root", "READ"); //path to database

    TString StrVoxelID;
	
    StrVoxelID.Form("%03i", voxel);
    TDirectory *pDir = (TDirectory*)file_voxel->Get(StrVoxelID);
	
    for (UInt_t iN = 0; iN < mNDet; iN++){

        for (UInt_t iM = iN + 1; iM < mNDet; iM++){

            sprintf(StrHist, "%03i%02i%02i", voxel, iN, iM);

            TH1* hist_voxel = (TH1*)pDir->Get(StrHist); //format for data base

            hist_voxel->Scale(1./hist_voxel->Integral());
            hist_voxel->Rebin(mRebin);

            if (mBins != (UInt_t)hist_voxel->GetNbinsX()) {
                cout << "ERROR in 'ImageReconstruction' line " << __LINE__ << ": 'mBins' does not match data." << endl;
            }

            for (UInt_t iBins = 0; iBins < mBins; iBins++) {

                Double_t BinContent = hist_voxel->GetBinContent(iBins + 1);
                if ( BinContent != BinContent ) {
                    //only true for nan
                    break;

                } else {
                    mN_dcb->at(iN).at(iM).at(iBins) = BinContent;
                }
            }

            delete hist_voxel;
        }
    }

    file_voxel->Close();
}

void ImageReconstruction::SetMeasN_dcb() {
    
    TString StrFileIn = GetInputFileName();
    TFile* file_voxel = new TFile(StrFileIn, "READ");  //path to measurement

    ///////////////////////////////////
    for (UInt_t iN = 0; iN < mNDet; iN++) {

        for (UInt_t iM = iN + 1; iM < mNDet; iM++) {

            sprintf(StrHist,"%02i%02i", iN, iM);
            hist_voxel[iN][iM] = (TH1*)file_voxel->Get(StrHist);
            hist_voxel[iN][iM]->Scale(1./hist_voxel[iN][iM]->Integral());
            hist_voxel[iN][iM]->Rebin(mRebin);
            hist_voxel[iN][iM]->SetLineColor(2);

            /////////////////////////////////
            if (mBins != (UInt_t)hist_voxel[iN][iM]->GetNbinsX()) {
                cout << "ERROR in 'ImageReconstruction' line " << __LINE__ << ": 'mBins' does not match data." << endl;
            }

            for (UInt_t iBins = 0; iBins < mBins; iBins++) {

                Double_t BinContent = hist_voxel[iN][iM]->GetBinContent(iBins + 1);
                    if ( BinContent != BinContent ) {
                        //only true for nan
                        for (UInt_t iNanBins = 0; iNanBins <= mBins + 1; iNanBins++) {
                            hist_voxel[iN][iM]->SetBinContent(iNanBins, 0);
                        }

                        break;

                    } else{
                        mN_dcb->at(iN).at(iM).at(iBins) = BinContent;
                    }
            }

            if (mEBinsMax == 0) {
                mEBinsMax = hist_voxel[iN][iM]->GetXaxis()->GetXmax();
                mEBinsMin = hist_voxel[iN][iM]->GetXaxis()->GetXmin();
            }
        }
    }

    file_voxel->Close();
}

void ImageReconstruction::CalculateA_v(UInt_t iterations){
    TBenchmark bench;

    UInt_t voxel;

    //init measurements
    bench.Start("N_dcb");
    SetMeasN_dcb();
    bench.Stop("N_dcb");
    cout << "\nN_dcb Creation Time: " << bench.GetRealTime("N_dcb") << " seconds\n";

    ShowA_v();
//    TString StrFileOut = GetOutputFileName();
//    TFile *pFile = new TFile(StrFileOut, "RECREATE");
//    pFile->cd();

//    TString StrHistName;
//    StrHistName.Form("%s_%i", hist2d->GetName(), 0);
//    hist2d->Write(StrHistName);
    
    bench.Start("calc");
    for (UInt_t it = 0; it < iterations; it++) {

        CalculateSum_v();

        for (UInt_t iNx = 0; iNx < mNx; iNx++){
            for (UInt_t iNy = 0; iNy < mNy; iNy++){
                for (UInt_t iNz = 0; iNz < mNz; iNz++) {
                    voxel = mVoxelIndex->at(iNx)[iNy][iNz];
                    mA_v->at(iNx).at(iNy).at(iNz) *= CalculateCorrFactorA_v(voxel);
                }
            }
        }

//        if (it % 10 == 0) {
//            cout << "Iteration: " << it << "\n";
//            ShowA_v();
//            // ShowN_dcbReconstructed();
//        }
		
//        hist2d->Reset();
//        for (UInt_t iNy = 0; iNy < mNy; iNy++){
//            for (UInt_t iNz = 0; iNz < mNz; iNz++) {
//                hist2d->SetBinContent(iNz + 1, iNy + 1, mA_v->at(0)[iNy][iNz]);
//                mIntegral2d += mA_v->at(0)[iNy][iNz];
//            }
//        }

//        pFile->cd();
//        StrHistName.Form("%s_%i", hist2d->GetName(), it + 1);
//        hist2d->Write(StrHistName);
    }
    bench.Stop("calc");
    cout << "Image reconstruction done. Steps: " << iterations << "\n";
    cout << "\nCalculation time: " << bench.GetRealTime("calc") << " seconds\n";

    hist2d->Reset();
    for (UInt_t iNy = 0; iNy < mNy; iNy++){
        for (UInt_t iNz = 0; iNz < mNz; iNz++) {
            hist2d->SetBinContent(iNz + 1, iNy + 1, mA_v->at(0)[iNy][iNz]);
            mIntegral2d += mA_v->at(0)[iNy][iNz];
        }
    }

    ShowA_v();
	
//    pFile->Save();
//    pFile->Close();
}

Double_t ImageReconstruction::CalculateCorrFactorA_v(UInt_t voxel){
    Double_t CorrFactor = 0;
    Double_t N_dcb;
    Double_t p_dcbv;

    for (UInt_t iN = 0; iN < mNDet; iN++){
        for (UInt_t iM = iN + 1; iM < mNDet; iM++){
            for (UInt_t iBins = 0; iBins < mBins; iBins++) {
                p_dcbv = mVoxelN_dcb[voxel]->at(iN)[iM][iBins];
                N_dcb  = mN_dcb->at(iN)[iM][iBins];

                if (sum_v->at(iN)[iM][iBins] != 0) {
                    CorrFactor += p_dcbv * N_dcb / sum_v->at(iN)[iM][iBins];
                }

            }
        }
    }

    return CorrFactor;
}

void ImageReconstruction::CalculateSum_v() {
    // Reconstruct matrix
    SetN_dcbToZero(sum_v);
    for (UInt_t iN = 0; iN < mNDet; iN++){
        for (UInt_t iM = 0; iM < mNDet; iM++){
            for (UInt_t iBins = 0; iBins < mBins; iBins++) {

                for (UInt_t iNx = 0; iNx < mNx; iNx++) {
                    for (UInt_t iNy = 0; iNy < mNy; iNy++){
                        for (UInt_t iNz = 0; iNz < mNz; iNz++) {

                            UInt_t voxel = mVoxelIndex->at(iNx)[iNy][iNz];
                            Double_t p_dcbv = mVoxelN_dcb[voxel]->at(iN)[iM][iBins];
                            Double_t A_k_v  = mA_v->at(iNx)[iNy][iNz];
                            sum_v->at(iN)[iM][iBins] += p_dcbv * A_k_v;
                        }
                    }
                }
            }
        }
    }
}

TCanvas* ImageReconstruction::ShowA_v() {
    mIntegral2d = 0;

    if (CanvasA_v == NULL) {
        CanvasA_v = new TCanvas("CanvasA_v", "CanvasA_v");
    }

    if (mDim <= 2) {
        if (hist2d == NULL) {
            hist2d = new TH2D("Hist2dA_v", "Hist2dA_v",
                              mNz, -(mNz / 2.0) * mVoxelLen, (mNz / 2.0) * mVoxelLen,
                              mNy, -(mNy / 2.0) * mVoxelLen, (mNy / 2.0) * mVoxelLen);
            gStyle->SetOptStat("ksiourmen");
            hist2d->GetXaxis()->SetTitle("#it{x} / mm");
            hist2d->GetYaxis()->SetTitle("#it{y} / mm");

        } else {
            hist2d->Reset();
        }

        for (UInt_t iNy = 0; iNy < mNy; iNy++){
            for (UInt_t iNz = 0; iNz < mNz; iNz++) {
                //~ hist2d->SetBinContent(iNz+1, iNy+1, mA_v->at(0)[iNy][iNz]);
                hist2d->SetBinContent(iNz + 1, iNy + 1, mA_v->at(0)[iNy][iNz]);
                mIntegral2d += mA_v->at(0)[iNy][iNz];
            }
        }
    }

    CanvasA_v->SetGrid();
    CanvasA_v->cd();
    hist2d->Scale(1./mIntegral2d); // the sum of all probability distributions on the graph is 1.

    // gStyle->SetNumberContours(50);
    // hist2d->Draw("CONT4Z");
    hist2d->Draw("COLZ");

    gPad->Update();
    CanvasA_v->Update();
    TPaletteAxis *palette = (TPaletteAxis*) hist2d->GetListOfFunctions()->FindObject("palette");
//    TColor *color = gROOT->GetColor(palette->GetValueColor(0)); // e.g. zero value
//    color->SetRGB(255, 255, 255); // all positions with value 0 are white

    return CanvasA_v;
}

TCanvas* ImageReconstruction::ShowN_dcbReconstructed() {

    //~ Double_t sum = 0;
    if (CanvasN_dcbReconstructed == NULL) {
        CanvasN_dcbReconstructed = new TCanvas("CanvasN_dcbReconstructed", "CanvasN_dcbReconstructed",
                                               600, 600);
        CanvasN_dcbReconstructed->Divide(mNDet, mNDet);
    }
	
    for (UInt_t iN = 0; iN < mNDet; iN++){
        for (UInt_t iM = iN + 1; iM < mNDet; iM++){
            if (histReconstructed[iN][iM] == NULL) {

                histReconstructed[iN][iM] = new TH1D(("HistN_dcbReconstructed" + to_string(iN + mNDet * iM)).c_str(),
                                                     ("HistN_dcbReconstructed" + to_string(iN + mNDet * iM)).c_str(),
                                                     mBins, mEBinsMin, mEBinsMax);
                gStyle->SetOptStat(0);

                histReconstructed[iN][iM]->GetXaxis()->SetTitle("Energy");
                histReconstructed[iN][iM]->GetYaxis()->SetTitle("Frequency");
                histReconstructed[iN][iM]->SetLineColor(4);

            } else {
                histReconstructed[iN][iM]->Reset();
            }
				
            if (mDim <= 2) {
                for (UInt_t iNy = 0; iNy < mNy; iNy++){
                    for (UInt_t iNz = 0; iNz < mNz; iNz++){
                        for (UInt_t iBins = 0; iBins < mBins; iBins++) {
                            UInt_t voxel = mVoxelIndex->at(0)[iNy][iNz];
                            Double_t p_dcbv = mVoxelN_dcb[voxel]->at(iN)[iM][iBins];
                            //~ histReconstructed[iN][iM]->AddBinContent(iBins+1,  sum_v->at(iN)[iM][iBins]);
                            histReconstructed[iN][iM]->AddBinContent(iBins+1, p_dcbv*mA_v->at(0)[iNy][iNz]);
                        }

                        if (histReconstructed[iN][iM]->Integral() != 0) {
                            histReconstructed[iN][iM]->Scale(1./histReconstructed[iN][iM]->Integral());
                        }

                    }
                }
            }

            CanvasN_dcbReconstructed->cd(iN + mNDet * iM + 1);
//            hist_voxel[iN][iM]->Draw("hist");
            histReconstructed[iN][iM]->Draw("histsame");
        }
    }

    CanvasN_dcbReconstructed->Update();
    return CanvasN_dcbReconstructed;
}

#endif
