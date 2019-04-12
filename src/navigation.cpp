// navigation.cpp

#include "navigation.h"

#include <iostream>
#include <boost/lexical_cast.hpp>


int Navigation::promptInput(){
    // prompt user for choosing an option

    std::string input;
    int choice;

    for (;;){
        std::cin >> input;

        try{
            choice = boost::lexical_cast<int>(input);
            break;  // success
        } catch(const boost::bad_lexical_cast &e){
            std::cerr << e.what() << '\n';  // print error message
            std::cout << "MUST BE INTEGER, SELECT: ";
        }
    }

    return choice;
}

std::string Navigation::promptPath(){
    // prompt user for a file path

    std::string filePath;
    std::cout << "TYPE PATH: ";
    std::cin >> filePath;
    return filePath;
}

int Navigation::checkInput(int choice, int lowerValue, int upperValue){
    // check input for validity
    // if not valid, call getInput()

    if (choice >= lowerValue && choice <= upperValue){
        return choice;
    } else{
        std::cout << "SELECT: ";
        return this->checkInput(this->promptInput(), lowerValue, upperValue);
    }
}

void Navigation::mainMenu(){
    // display main menu

    std::cout << "MAIN MENU\n";
    std::cout << "1 - Image Reconstruction\n";
    std::cout << "2 - Quit\n";
    std::cout << "SELECT: ";
}

void Navigation::measurementsMenu(){
    // display image reconstruction menu

    std::cout << "IMAGE RECONSTRUCTION MENU\n";
    std::cout << "1 - Choose *.root file\n";
    std::cout << "2 - Back\n";
    std::cout << "SELECT: ";
}

void Navigation::algorithmMenu(){
    // display algorithm menu

    std::cout << "ALGORITHM MENU\n";
    std::cout << "1 - Maximum Likelihood Expectation Maximization\n";
    std::cout << "2 - Origin Ensemble\n";
    std::cout << "SELECT: ";
}




//int reconstructionPrompt(){
//    // prompt user to set the reconstruction algorithm

//    int answer;
//    answer = 0;

//    std::cout << "Mode:" << std::endl;


//    return answer;
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
