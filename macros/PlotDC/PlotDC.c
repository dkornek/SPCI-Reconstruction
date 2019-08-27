// Plot projection probabilities of the image space of a fixed detector pair in the measurement space

void plotProbabilities(TString pair){
    TString pathToSystemMatrix = "../../data/SystemMatrix/SPCIBase441_B700.root";
    TFile* input = new TFile(pathToSystemMatrix, "READ");

    if (!input->IsOpen()){
        std::cout << "File not found!\nShut down macro.\n";
    }

    TString dc = pair;
    TString name;
    // name.Form("Sensitivity of Image Space for %s", dc.Data());
    name.Form("");

    TH2F* sensitivityOfDC = new TH2F("sens", name,
                                     21, -52.5, 52.5,
                                     21, -52.5, 52.5);

    sensitivityOfDC->SetStats(kFALSE);
    sensitivityOfDC->SetContour(99);

    sensitivityOfDC->GetXaxis()->SetTitle("#it{x} / mm");
    sensitivityOfDC->GetXaxis()->SetTitleOffset(1);
    sensitivityOfDC->GetXaxis()->SetNdivisions(5);
    sensitivityOfDC->GetXaxis()->SetLabelSize(0.07);
    sensitivityOfDC->GetXaxis()->SetTitleSize(0.07);

    sensitivityOfDC->GetYaxis()->SetTitle("#it{y} / mm");
    sensitivityOfDC->GetYaxis()->SetTitleOffset(1);
    sensitivityOfDC->GetYaxis()->SetNdivisions(5);
    sensitivityOfDC->GetYaxis()->SetLabelSize(0.07);
    sensitivityOfDC->GetYaxis()->SetTitleSize(0.07);

    sensitivityOfDC->GetZaxis()->SetTitle("Projection Probability / %");
    sensitivityOfDC->GetZaxis()->SetTitleOffset(1.6);
    sensitivityOfDC->GetZaxis()->SetNdivisions(5);
    sensitivityOfDC->GetZaxis()->SetLabelSize(0.07);
    sensitivityOfDC->GetZaxis()->SetTitleSize(0.07);

    Int_t yVoxel = 0;
    TIter nextVoxel(input->GetListOfKeys());
    TKey* keyVoxel;
    while ((keyVoxel = (TKey*)nextVoxel())){

        TString nameOfVoxel = keyVoxel->GetName();
        TDirectory* directoryOfVoxel = (TDirectory*)input->Get(nameOfVoxel);

        Int_t globalVoxel = nameOfVoxel.Atoi();
        Int_t xVoxel = globalVoxel % 21 + 1;

        if ((globalVoxel % 21) == 0){
            ++yVoxel;
        }

        TIter nextEDS(directoryOfVoxel->GetListOfKeys());
        TKey* keyEDS;
        while ((keyEDS = (TKey*)nextEDS())){

            TString nameOfEDS = keyEDS->GetName();
            TString nameOfEDSReduced = (TString)nameOfEDS(3, 4);

            if (nameOfEDSReduced == dc){
                TH1F* EDSInVoxel = (TH1F*)directoryOfVoxel->Get(nameOfEDS);
                EDSInVoxel->Scale(1.0 / 5000000);

                Double_t binContent = 0;
                for (Int_t bin = 1;; ++bin){

                    if (EDSInVoxel->IsBinOverflow(bin)){
                        break;
                    }
                    binContent += EDSInVoxel->GetBinContent(bin);
                }

                sensitivityOfDC->SetBinContent(xVoxel, yVoxel, binContent);
            }
        }
        delete keyEDS;
        delete directoryOfVoxel;
    }
    delete keyVoxel;

    gStyle->SetPalette(55, 0);
    TGaxis::SetMaxDigits(4);

    TCanvas* c = new TCanvas("c", "Plot", 450, 410);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.26);
    c->SetBottomMargin(0.15);

    sensitivityOfDC->Draw("COLZ2");
    c->Update();

    name.Form("%s.png", dc.Data());
    c->Print(name);

    delete sensitivityOfDC;
    delete c;
    input->Close();
}

void PlotDC(){

    std::vector<TString> pairs = {"0506", "0509", "0015"};

    for (UInt_t p = 0; p < pairs.size(); ++p){
        plotProbabilities(pairs.at(p));
    }
}

#ifndef __CINT__
int main(){
    PlotDC();

    return 0;
}
#endif
