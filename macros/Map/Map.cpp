// ROOT Macro to create a map for the defined coordinate systems

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-05-10

#include "/home/MED/koeglerto/Programs/Root-Work/StyleSheets/StyleSheet.C"

void Map(){
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle();

    Int_t numberOfVoxelsX = 3; //20;
    Int_t numberOfVoxelsY = 3; //20;

    Double_t rangeXmin = -16; //-60;
    Double_t rangeXmax = 16; //60;
    Double_t rangeYmin = -16; //-60;
    Double_t rangeYmax = 16; //60;

    // prepare the layout
    TLatex text;
    text.SetTextColor(4);
    text.SetTextSize(0.08); //(0.03);
    text.SetTextAlign(33);

    TMarker marker;
    marker.SetMarkerColor(4);
    marker.SetMarkerStyle(20);
    marker.SetMarkerSize(0.8);
    TString nameOfMarker;
    Int_t numberOfMarker = 0;
    Double_t markerX;  // position of marker in x axis
    Double_t markerY;  // position of marker in y axis
    Double_t sizeX = 8.15; //5;  // size (mm) of marker in x axis
    Double_t sizeY = 8.15; //5;  // size (mm) of marker in y axis

    TCanvas *canvas = new TCanvas("base", "BaseGrid", 250, 250);
    TH2F *histogram = new TH2F("BaseGridMap", "SCI Base Grid",
                               700, rangeXmin, rangeXmax,
                               700, rangeYmin, rangeYmax);
    histogram->GetXaxis()->SetTitle("#it{x} / mm");
    histogram->GetYaxis()->SetTitle("#it{y} / mm");
    histogram->Draw();

    // draw the markers
    for (Int_t j = 0; j <= numberOfVoxelsY; ++j){
        for (Int_t i = 0; i <= numberOfVoxelsX; ++i){

            // draw the marker
            markerX = i * sizeX - rangeXmax + 4.075; // + 10;
            markerY = j * sizeY - rangeYmax + 4.075; // + 10;
            marker.DrawMarker(markerX, markerY);

            if (i % 1 == 0 && j % 1 == 0){
                nameOfMarker.Form("%i", numberOfMarker);
                text.DrawLatex(markerX - 0.5, markerY - 0.5, nameOfMarker);
            }

            ++numberOfMarker;
        }
    }

    canvas->Update();
}
