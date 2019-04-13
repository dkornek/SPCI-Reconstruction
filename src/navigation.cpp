// navigation.cpp

#include "navigation.h"

#include <boost/lexical_cast.hpp>


// ##### TEMPLATE MENU #####
void TemplateMenu::promptChoice(int& iChoice){
    // prompt user for choosing an option

    std::string input;
    for (;;){
        std::cout << "SELECT: ";
        std::cin >> input;  // get user input

        try{
            iChoice = boost::lexical_cast<int>(input);
            break;  // success
        } catch(const boost::bad_lexical_cast &e){
            std::cerr << e.what() << "\n";  // print error message
        }
    }
}

std::string TemplateMenu::promptPath(){
    // prompt user for a file path

    std::string input;
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
    iNextMenu = 0;

    int iChoice;
    this->promptChoice(iChoice);

    switch (iChoice){
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
    iNextMenu = 0;

    std::string iPath;

    int iChoice;
    this->promptChoice(iChoice);

    switch (iChoice) {
    case 1:
        iPath = this->promptPath();
        std::cout << "You said: " << iPath << "\n";

        // if path found -> algorithm menu

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


//int Navigation::checkInput(int choice, int lowerValue, int upperValue){
//    // check input for validity
//    // if not valid, call getInput()

//    if (choice >= lowerValue && choice <= upperValue){
//        return choice;
//    } else{
//        std::cout << "SELECT: ";
//        return this->checkInput(this->promptInput(), lowerValue, upperValue);
//    }


//void Navigation::algorithmMenu(){
//    // display algorithm menu

//    std::cout << "ALGORITHM MENU\n";
//    std::cout << "1 - Maximum Likelihood Expectation Maximization\n";
//    std::cout << "2 - Origin Ensemble\n";
//    std::cout << "SELECT: ";
//}



//int projectionMatrixPrompt(){
//    // prompt user to set the projection matrix
//    // the projection matrix will be used for image reconstruction in ML-EM

//    std::cout << "Set the projection matrix for image reconstruction:" << std::endl;
//    std::cout << "1: Create projection matrix" << std::endl;
//    std::cout << "2: Load projection matrix" << std::endl;
//    std::cout << "Please choose option (1-2):";

//    return 0;
//}
