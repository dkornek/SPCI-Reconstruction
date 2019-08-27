// Plot projection probabilities of the image space of a fixed detector element in the measurement space

void plotProbabilities(TString pair, Int_t bin){
    TString pathToSystemMatrix = "../../data/SystemMatrix/SPCIBase441_B700.root";
    TFile* input = new TFile(pathToSystemMatrix, "READ");

    if (!input->IsOpen()){
        cout << "File not found!\nShut down macro.\n";
    }

    TString dc = pair;
    Int_t b = bin;
    TString name;
    name.Form("");

    TH2F* sensitivityOfDCB = new TH2F("sens", name,
                                      21, -52.5, 52.5,
                                      21, -52.5, 52.5);

    sensitivityOfDCB->SetStats(kFALSE);
    sensitivityOfDCB->SetContour(99);

    sensitivityOfDCB->GetXaxis()->SetTitle("#it{x} / mm");
    sensitivityOfDCB->GetXaxis()->SetTitleOffset(1);
    sensitivityOfDCB->GetXaxis()->SetNdivisions(5);
    sensitivityOfDCB->GetXaxis()->SetLabelSize(0.07);
    sensitivityOfDCB->GetXaxis()->SetTitleSize(0.07);

    sensitivityOfDCB->GetYaxis()->SetTitle("#it{y} / mm");
    sensitivityOfDCB->GetYaxis()->SetTitleOffset(1);
    sensitivityOfDCB->GetYaxis()->SetNdivisions(5);
    sensitivityOfDCB->GetYaxis()->SetLabelSize(0.07);
    sensitivityOfDCB->GetYaxis()->SetTitleSize(0.07);

    sensitivityOfDCB->GetZaxis()->SetTitle("Projection Probability / %");
    sensitivityOfDCB->GetZaxis()->SetTitleOffset(1.6);
    sensitivityOfDCB->GetZaxis()->SetNdivisions(5);
    sensitivityOfDCB->GetZaxis()->SetLabelSize(0.07);
    sensitivityOfDCB->GetZaxis()->SetTitleSize(0.07);

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
                sensitivityOfDCB->SetBinContent(xVoxel, yVoxel,
                                    EDSInVoxel->GetBinContent(b));

                delete EDSInVoxel;
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

    sensitivityOfDCB->Draw("COLZ2");
    c->Update();

    name.Form("%s_%03i.png", dc.Data(), b);
    c->SaveAs(name);

    delete sensitivityOfDCB;
    delete c;
    input->Close();
}

void PlotDCB(){

    std::vector<TString> pairs = {"0506", "0509", "0015"};
    std::vector<Int_t> bins = {175, 350, 525};

    for (Int_t p = 0; p < pairs.size(); ++p){
        for (Int_t b = 0; b < bins.size(); ++b){
            plotProbabilities(pairs.at(p), bins.at(b));
        }
    }
}
