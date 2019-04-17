// navigation.cpp

#include "navigation.h"
#include "mlem.h"
#include <boost/lexical_cast.hpp>
#include <TFile.h>


int promptChoice(std::string iPrompt = "SELECT: "){
    // prompt user for choosing an option

    std::string input;
    int iChoice;

    for (;;){
        std::cout << iPrompt;
        std::cin >> input;  // get user input

        try{
            iChoice = boost::lexical_cast<int>(input);
            break;  // success
        } catch(const boost::bad_lexical_cast &e){
            std::cerr << e.what() << "\n";  // print error message
        }
    }

    return iChoice;
}

bool checkPath(TString iPath){
    // check validity of path to .root file

    TFile iFile(iPath, "READ");
    if (iFile.IsOpen()){
        std::cout << "File opened successfully.\n";
        iFile.Close();
        return true;
    } else{
        std::cout << "Invalid file path.\n";
        return false;
    }
}

TString promptPath(){
    // prompt user for a file path

    TString input;
    std::cout << "TYPE PATH TO FILE: ";
    std::cin >> input;

    return input;
}

// ##### MAIN MENU #####
MainMenu::MainMenu(){
    // display main menu

    std::string title("\nMAIN MENU\n");
    std::string option1("1 - Image Reconstruction\n");
    std::string option2("2 - Quit\n");

    messageText = title + option1 + option2;
}

TemplateMenu* MainMenu::getNextMenu(bool &iIsQuitOptionSelected){
    // prompt user for next menu

    TemplateMenu* iNextMenu;
    iNextMenu = nullptr;

    switch (promptChoice()){
    case 1:
        iNextMenu = new ReconstructionMenu;
        break;
    case 2:
        std::cout << "\nApplication shut down.\n";
        iIsQuitOptionSelected = true;
        break;
    default:
        break;
    }

    return iNextMenu;
}

// ##### MEASUREMENTS MENU #####
ReconstructionMenu::ReconstructionMenu(){
    // display image reconstruction menu

    std::string title("\nIMAGE RECONSTRUCTION MENU\n");
    std::string option1("1 - Choose *.root file of measurements\n");
    std::string option2("2 - Back\n");

    messageText = title + option1 + option2;
}

TemplateMenu* ReconstructionMenu::getNextMenu(bool &iIsQuitOptionSelected){
    // prompt user to locate a measurements file or to go one step back

    TemplateMenu* iNextMenu;
    iNextMenu = nullptr;

    TString iPath;
    bool isPathValid = false;

    switch (promptChoice()) {
    case 1:
        while (!isPathValid){
            // prompt user until path is valid
            // iPath = promptPath();
            iPath = "../test/Position_x1.5_y-1.5.root";
            isPathValid = checkPath(iPath);
        }

        iNextMenu = new AlgorithmMenu(iPath);
        break;
    case 2:
        iNextMenu = new MainMenu;
        break;
    default:
        break;
    }

    iIsQuitOptionSelected = false;
    return iNextMenu;
}

// ##### ALGORITHM MENU #####
AlgorithmMenu::AlgorithmMenu(TString iPathToMeasurements){
    // display algorithm menu

    std::string title("\nALGORITHM MENU\n");
    std::string option1("1 - Maximum Likelihood Expectation Maximization\n");
    std::string option2("2 - Origin Ensemble\n");
    std::string option3("3 - Back\n");

    messageText = title + option1 + option2 + option3;

    pathToMeasurements = iPathToMeasurements;
}

TemplateMenu* AlgorithmMenu::getNextMenu(bool &iIsQuitOptionSelected){
    // prompt user to choose the algorithm

    TemplateMenu* iNextMenu;
    iNextMenu = nullptr;

    switch (promptChoice()) {
    case 1:
        // MLEM Menue
        iNextMenu = new MLEMMenu(pathToMeasurements);
        break;
    case 2:
        // OE Menu
        std::cout << "Not yet implemented!\n";
        break;
    case 3:
        iNextMenu = new ReconstructionMenu;
    default:
        break;
    }

    iIsQuitOptionSelected = false;
    return iNextMenu;
}

// ##### MLEM MENU #####
MLEMMenu::MLEMMenu(TString iPathToMeasurements){
    // display MLEM menu

    std::string title("\nMLEM MENU\n");
    std::string option1("1 - Choose *.root file of projection matrix\n");
    std::string option2("2 - Back\n");

    messageText = title + option1 + option2;

    pathToMeasurements = iPathToMeasurements;
}

TemplateMenu* MLEMMenu::getNextMenu(bool &iIsQuitOptionSelected){
    // prompt user to select the file of the projections

    TemplateMenu* iNextMenu;
    iNextMenu = nullptr;

    TString iPath;
    bool isPathValid = false;

    switch (promptChoice()) {
    case 1:
        while (!isPathValid){
            // prompt user until path is valid
            // iPath = promptPath();
            iPath = "../test/SPCIBase49.root";
            isPathValid = checkPath(iPath);
        }

        iNextMenu = new MLEMRecoMenu(pathToMeasurements, iPath);
        break;
    case 2:
        iNextMenu = new AlgorithmMenu(pathToMeasurements);
        break;
    default:
        break;
    }

    iIsQuitOptionSelected = false;
    return iNextMenu;
}

MLEMRecoMenu::MLEMRecoMenu(TString iPathToMeasurements, TString iPathToProjections){
    // display MLEM Reconstruction Menu

    std::string title("\nMLEM ALGORITHM MENU\n");
    std::string option1("1 - Set number of iterations\n");
    std::string option2("2 - Back\n");

    messageText = title + option1 + option2;

    pathToMeasurements = iPathToMeasurements;
    pathToProjections = iPathToProjections;
}

TemplateMenu* MLEMRecoMenu::getNextMenu(bool &iIsQuitOptionSelected){
    // prompt user to set the iterations

    TemplateMenu* iNextMenu;
    iNextMenu = nullptr;

    int iNumberOfIterations;

    switch (promptChoice()) {
    case 1:
        iNumberOfIterations = promptChoice("NUMBER OF ITERATIONS: ");
        std::cout << "Calculating ...\n";

        // start reco, save image to disk, then show it

        iNextMenu = new MainMenu();
        break;
    case 2:
        iNextMenu = new MLEMMenu(pathToMeasurements);
    default:
        break;
    }

    iIsQuitOptionSelected = false;
    return iNextMenu;
}
