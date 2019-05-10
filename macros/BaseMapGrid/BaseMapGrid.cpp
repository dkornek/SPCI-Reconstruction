// ROOT Macro to create a base map grid for the corresponding GEANT4 simulation

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-05-10

#include "/home/MED/koeglerto/Programs/Root-Work/StyleSheets/StyleSheet.C"

void BaseMapGrid()
{
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle();

    Int_t numberOfVoxelsX = 20;
    Int_t numberOfVoxelsY = 20;

    Double_t rangeXmin = -60.0;
    Double_t rangeXmax = 60.0;
    Double_t rangeYmin = -60.0;
    Double_t rangeYmax= 60.0;

    // prepare the layout
    TLatex text;
    text.SetTextColor(4);
    text.SetTextSize(0.02);
    text.SetTextAlign(33);

    TMarker marker;
    marker.SetMarkerColor(4);
    marker.SetMarkerStyle(20);
    marker.SetMarkerSize(0.75);
    TString nameOfMarker;
    Int_t numberOfMarker = 0;
    Double_t markerX;  // position of marker in x axis
    Double_t markerY;  // position of marker in y axis
    Double_t sizeX = 5.0;  // size (mm) of marker in x axis
    Double_t sizeY = 5.0;  // size (mm) of marker in y axis

    TCanvas *canvas = new TCanvas("base", "BaseGrid", 800, 600);
    TH2F *histogram = new TH2F("BaseGridMap", "SCI Base Grid",
                               700, rangeXmin, rangeXmax,
                               700, rangeYmin, rangeYmax);
    histogram->GetXaxis()->SetTitle("#it{x} / mm");
    histogram->GetYaxis()->SetTitle("#it{y} / mm");
    histogram->Draw();

    // draw the markers
    for (Int_t j = numberOfVoxelsY; j >= 0; --j){
        for (Int_t i = 0; i <= numberOfVoxelsX; ++i){

            // draw the marker
            markerX = i * sizeX - rangeXmax + 10;
            markerY = j * sizeY - rangeYmax + 10;
            marker.DrawMarker(markerX, markerY);

            if (i % 2 == 0){
                nameOfMarker.Form("%i", numberOfMarker);
                text.DrawLatex(markerX - 0.5, markerY - 0.5, nameOfMarker);
            }

            ++numberOfMarker;
        }
    }

    canvas->Update();
