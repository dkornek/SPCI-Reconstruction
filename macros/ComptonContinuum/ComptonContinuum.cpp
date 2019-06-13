// Show the underlying Compton continuum of a difference spectrum

void ComptonContinuum(){

    TString pathToMeasurement = "../../data/Measurements/SPCIBase441/Bins700/SourceH.root";
    TFile* input = new TFile(pathToMeasurement, "READ");

    TList* keysList = input->GetListOfKeys();
    TH1F* diffSpectrum = (TH1F*)input->Get(keysList->First()->GetName());
    TAxis* xAxisDiff = diffSpectrum->GetXaxis();

    TH1F* continuumE1 = new TH1F("cont", "Compton Continuum", 700, 0, 700);
    TAxis* xAxisCont = continuumE1->GetXaxis();

    Int_t numberOfBins = diffSpectrum->GetNbinsX();
    for (Int_t b = 1; b <= numberOfBins; ++b){
        Double_t energyLow = xAxisDiff->GetBinLowEdge(b);
        Double_t energyUp = xAxisDiff->GetBinUpEdge(b);
        Double_t content = diffSpectrum->GetBinContent(b);

        Double_t energyCenter = xAxisDiff->GetBinCenter(b);
        Double_t e1 = 0.5 * (662 + energyCenter);
        Double_t e2 = 0.5 * (662 - energyCenter);

        Double_t binInCont = xAxisCont->FindBin(e2);
        continuumE1->SetBinContent(binInCont, content);

        cout << energyLow << " - " << energyUp << "\t" << content << "\t" << e1 << " - " << e2 << "\n";
    }

    TCanvas* c = new TCanvas("diffSpec", "Difference Spectrum", 500, 500);
    continuumE1->Draw();
    c->Update();
}
