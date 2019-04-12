/*
This project executes the image reconstruction for SPCI measurement data.
To do so, two different algorithms can be used:
- maximum likelihood expectation maximization (ML-EM),
- origin ensemble (OE).

Author: Dominik Kornek <dominik.kornek@gmail.com>
*/


/*
TO-DO:

1) create projection matrix:
    1.1 read measurement data (histograms)
    1.2 get N_dbc (Number of events of detector pair d/c in bin b)
    1.3 model real activity distribution provided by the measurements A_v
    1.3 calculate matrix elements p_dcbv

2) implement ML-EM algorithm:
    2.1 create a-priori activity distribution (e.g. homogeneous distribution)
    2.2 implement formula

3) implement OE algorithm:
    3.1 ?

4) analysis:
    4.1 create plots for every iteration
    4.2 update graph of maximization (for ML-EM (log?))
    4.3 generate file with image properties (e.g. SNR, resolution ...)
*/

#include <iostream>
#include <string>

#include "navigation.h"


int main(){

    // start application
    Navigation menu;
    int choice;
    std::string measurementsPath;

    for (;;){
        // choose option from main menu
        menu.mainMenu();
        choice = menu.promptInput();
        choice = menu.checkInput(choice, 1, 2);

        if (choice == 2){
            // quit application
            std::cout << "Application shut down.\n";
            break;

        } else{
            // choose option from image reconstruction menu
            menu.measurementsMenu();
            choice = menu.promptInput();
            choice = menu.checkInput(choice, 1, 2);

            if (choice == 2){
                // go one step back
                continue;

            } else{
                // prompt for path of measurement data
                measurementsPath = menu.promptPath();
                std::cout << "Location: " << measurementsPath << std::endl;

                // check for valid path
                // if (invalid path){
                //     try again
                // } else{
                //     // choose option from algorithm menu
                // }
            }
        }
    }

    return 0;
}
