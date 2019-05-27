// navigation.cpp

#include "navigation.h"

#include "mlem.h"
#include "image_processing.h"

#include <boost/lexical_cast.hpp>
#include <TFile.h>


int promptChoice(std::string promptString = "SELECT: "){
    // prompt user for choosing an option

    std::string input;
    int choice;

    for (;;){
        std::cout << promptString;
        std::cin >> input;  // get user input

        try{
            choice = boost::lexical_cast<int>(input);
            break;  // success

        } catch(const boost::bad_lexical_cast &e){
            std::cerr << e.what() << "\n";  // print error message
        }
    }

    return choice;
}

double promptParameter(std::string promptString, double lowerLimit, double upperLimit){
    // prompt user for choosing any needed parameter

    std::string input;
    double choice;

    for (;;){
        std::cout << promptString;
        std::cin >> input;  // get user input

        try{
            choice = boost::lexical_cast<double>(input);

            if ((choice < lowerLimit) || (choice >= upperLimit)){
                continue;

            } else{
                break;  // success
            }

        } catch(const boost::bad_lexical_cast &e){
            std::cerr << e.what() << "\n";  // print error message
        }
    }

    return choice;
}

bool checkPath(TString path){
    // check validity of path to .root file

    TFile* file = new TFile(path, "READ");

    if (file->IsOpen()){
        std::cout << "File opened successfully.\n";
        file->Close();
        return true;

    } else{
        std::cout << "Invalid file path.\n";
        return false;
    }
}

TString promptPath(std::string promptString){
    // prompt user for a file path

    TString input;

    for (;;){
        std::cout << promptString;
        std::cin >> input;  // get user input

        if (checkPath(input)){
            break;  // success
        }
    }

    return input;
}

TString promptFileName(){
    // prompt user for a file name

    TString input;
    std::cout << "TYPE PATH TO SAVED FILE: ";
    std::cin >> input;
    return input;
}

// ##### MAIN MENU #####
MainMenu::MainMenu(){
    // display main menu

    std::string title("\nMAIN MENU\n");
    std::string option1("1 - Quit\n");
    std::string option2("2 - Image Reconstruction\n");

    messageText = title + option1 + option2;
}

TemplateMenu* MainMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user for next menu

    TemplateMenu* nextMenu;
    nextMenu = nullptr;

    switch (promptChoice()){
    case 1:
        std::cout << "\nApplication shut down.\n";
        isQuitOptionSelected = true;
        break;

    case 2:
        nextMenu = new ReconstructionMenu();
        break;

    default:
        break;
    }

    return nextMenu;
}

// ##### MEASUREMENTS MENU #####
ReconstructionMenu::ReconstructionMenu(){
    // display image reconstruction menu

    std::string title("\nIMAGE RECONSTRUCTION MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Choose *.root files\n");

    messageText = title + option1 + option2;
}

TemplateMenu* ReconstructionMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user to specify the files

    TemplateMenu* nextMenu;
    nextMenu = nullptr;

    TString pathToMeasurements;
    TString pathToProjections;

    switch (promptChoice()) {
    case 1:
        nextMenu = new MainMenu();
        break;

    case 2:
//        pathToMeasurements = promptPath("TYPE PATH TO MEASUREMENTS FILE: ");
        pathToMeasurements = "../data/Measurements/SPCIBase49/Bins50/SourcePos0+24+26+36.root";

//        pathToProjections = promptPath("TYPE PATH TO PROJECTIONS FILE: ");
        pathToProjections = "../data/SystemMatrix/Bins50/SPCIBase49.root";

        nextMenu = new AlgorithmMenu(pathToMeasurements, pathToProjections);
        break;

    default:
        break;
    }

    isQuitOptionSelected = false;
    return nextMenu;
}

// ##### ALGORITHM MENU #####
AlgorithmMenu::AlgorithmMenu(TString filePathToMeasurements, TString filePathToProjections) :
    pathToMeasurements(filePathToMeasurements), pathToProjections(filePathToProjections){
    // display algorithm menu

    std::string title("\nALGORITHM MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Maximum Likelihood Expectation Maximization\n");
    std::string option3("3 - Origin Ensemble\n");

    messageText = title + option1 + option2 + option3;
}

TemplateMenu* AlgorithmMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user to choose the algorithm

    TemplateMenu* nextMenu;
    nextMenu = nullptr;

    switch (promptChoice()) {
    case 1:
        nextMenu = new ReconstructionMenu();
        break;

    case 2:
        nextMenu = new MLEMMenu(pathToMeasurements, pathToProjections);
        break;

    case 3:
        std::cout << "Not yet implemented!\n";
        nextMenu = new MainMenu();
        break;

    default:
        break;
    }

    isQuitOptionSelected = false;
    return nextMenu;
}

// ##### MLEM MENU #####
MLEMMenu::MLEMMenu(TString filePathToMeasurements, TString filePathToProjections) :
    pathToMeasurements(filePathToMeasurements), pathToProjections(filePathToProjections){
    // display MLEM menu

    std::string title("\nMLEM MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Choose calculation criteria\n");

    messageText = title + option1 + option2;
}

TemplateMenu* MLEMMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user to enter calculation criteria

    TemplateMenu* nextMenu;
    nextMenu = nullptr;

    TString saveDataName;
    int maxNumberOfIterations;
    double stopCriterion;
    double accelerator;
    ReconstructionMLEM* reco;
    ActivityDistribution* activity;

    switch (promptChoice()) {
    case 1:
        nextMenu = new AlgorithmMenu(pathToMeasurements, pathToProjections);
        break;

    case 2:
        saveDataName = promptFileName();
//        saveDataName = "../test/Activity.root";
        maxNumberOfIterations = promptChoice("CHOOSE MAXIMUM NUMBER OF ITERATIONS: ");
//        maxNumberOfIterations = 300;
        stopCriterion = promptParameter("CHOOSE A STOPPING CRITERION ( 0 < ... < 1, 0 = NO ABORTION): ", 0.0, 1.0);
//        stopCriterion = 0;
        accelerator = promptParameter("CHOOSE AN OVER-RELAXATION PARAMETER ( 1 < ... < 2, 1 = NO OVER-RELAXATION: ", 1.0, 2.0);
//        accelerator = 1.5;

        reco = new ReconstructionMLEM(pathToMeasurements, pathToProjections);
        reco->setAccelerator(accelerator);
        reco->setActivityThreshold(0);                      // make function
        reco->setImageVolume({-50, 50, -50, 50, 5, 10});    // make function
        reco->start(maxNumberOfIterations, stopCriterion);

        activity = new ActivityDistribution(saveDataName);
        activity->saveAll(reco->A_v, reco->A_vProject3DSteps, reco->plotChi);

        nextMenu = new MainMenu();
        break;

    default:
        break;
    }

    delete activity;
    delete reco;
    isQuitOptionSelected = false;
    return nextMenu;
}
