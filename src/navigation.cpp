// navigation.cpp

#include "navigation.h"
#include "user_prompts.h"
#include "mlem.h"
#include "oe.h"

TString pathToSystemMatrix;
TString pathToMeasurements;

// ##### MAIN MENU #####
MainMenu::MainMenu(){
    // display main menu

    std::string title("\nMAIN MENU\n");
    std::string option1("1 - Quit\n");
    std::string option2("2 - Start Image Reconstruction\n");

    messageText = title + option1 + option2;
}

TemplateMenu* MainMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user for shut down or image reconstruction

    TemplateMenu* nextMenu = nullptr;
    switch (promptChoice()){
        case 1:
            std::cout << "\nApplication shut down.\n";
            isQuitOptionSelected = true;
            break;

        case 2:
            nextMenu = new SystemMatrixMenu();
            break;

        default:
            break;
    }
    return nextMenu;
}

// ##### RECONSTRUCTION MENU #####
SystemMatrixMenu::SystemMatrixMenu(){
    // display system matrix menu

    std::string title("\nSYSTEM MATRIX MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Choose System Matrix (*.root)\n");

    messageText = title + option1 + option2;
}

TemplateMenu* SystemMatrixMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user for system matrix

    TemplateMenu* nextMenu = nullptr;
    switch (promptChoice()){
        case 1:
            nextMenu = new MainMenu();
            break;

        case 2:
            pathToSystemMatrix = promptPath("TYPE PATH TO SYSTEM MATRIX: ");
            nextMenu = new MeasurementsMenu();
            break;

        default:
            break;
    }

    isQuitOptionSelected = false;
    return nextMenu;
}

MeasurementsMenu::MeasurementsMenu(){
    // display measurements menu

    std::string title("\nMEASUREMENTS MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Choose Measurements (*.root)\n");

    messageText = title + option1 + option2;
}

TemplateMenu* MeasurementsMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user for measurements file

    TemplateMenu* nextMenu = nullptr;
    switch (promptChoice()){
        case 1:
            nextMenu = new SystemMatrixMenu();
            break;

        case 2:
            pathToMeasurements = promptPath("TYPE PATH TO MEASUREMENTS: ");
            nextMenu = new AlgorithMenu();
            break;

        default:
            break;
    }

    isQuitOptionSelected = false;
    return nextMenu;
}

AlgorithMenu::AlgorithMenu(){
    // display algorithms menu

    std::string title("\nALGORITHM MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Maximum-Likelihood Expectation-Maximization\n");
    std::string option3("3 - Origin Ensemble\n");

    messageText = title + option1 + option2 + option3;
}

TemplateMenu* AlgorithMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user for algorithm

    TemplateMenu* nextMenu = nullptr;
    switch (promptChoice()){
        case 1:
            nextMenu = new MeasurementsMenu();
            break;

        case 2:
            nextMenu = new MLEMMenu();
            break;

        case 3:
            nextMenu = new OEMenu();

        default:
            break;
    }

    isQuitOptionSelected = false;
    return nextMenu;
}

MLEMMenu::MLEMMenu(){
    // display ML-EM menu

    std::string title("\nMAXIMUM-LIKELIHOOD EXPECTATION-MAXIMIZATION MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Use standard settings\n");
    std::string option3("3 - Advanced\n");

    messageText = title + option1 + option2 + option3;
}

TemplateMenu* MLEMMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user for measurements file

    TBenchmark b;
    int iterations;
    double accelerator;
    ReconstructionMLEM* reco = nullptr;

    TemplateMenu* nextMenu = nullptr;
    switch (promptChoice()){
        case 1:
            nextMenu = new AlgorithMenu();
            break;

        case 2:
            b.Start("totalMLEM");

            reco = new ReconstructionMLEM(pathToMeasurements,
                                          pathToSystemMatrix,
                                          {-52.5, 52.5, -52.5, 52.5, 5, 10});
            accelerator = 1.9;
            reco->setAccelerator(accelerator);

            iterations = promptChoice("NUMBER OF ITERATIONS: ");
            reco->start(iterations);

            b.Stop("totalMLEM");
            std::cout << "\nTotal Time:\t\t" << b.GetRealTime("totalMLEM") << " seconds\n";

            nextMenu = new MainMenu();
            break;

        case 3:
            b.Start("totalMLEM2");

            reco = new ReconstructionMLEM(pathToMeasurements,
                                          pathToSystemMatrix,
                                          {-52.5, 52.5, -52.5, 52.5, 5, 10});

            accelerator = promptParameter("SET THE EXPONENT ( 1 < ... < 2, 1 = NO ACCELERATION: ", 1.0, 2.0);
            reco->setAccelerator(accelerator);

            iterations = promptChoice("NUMBER OF ITERATIONS: ");
            reco->start(iterations);

            b.Stop("totalMLEM2");
            std::cout << "\nTotal Time:\t\t" << b.GetRealTime("totalMLEM2") << " seconds\n";

            nextMenu = new MainMenu();
            break;

        default:
            break;
    }

    isQuitOptionSelected = false;
    delete reco;
    return nextMenu;
}

OEMenu::OEMenu(){
    // display OE menu

    std::string title("\nORIGIN ENSEMBLE MENU\n");
    std::string option1("1 - Back\n");
    std::string option2("2 - Choose settings\n");

    messageText = title + option1 + option2;
}

TemplateMenu* OEMenu::getNextMenu(bool& isQuitOptionSelected){
    // prompt user for measurements file

    TBenchmark b;
    int states, samples;
    ReconstructionOE* reco = nullptr;

    TemplateMenu* nextMenu = nullptr;
    switch (promptChoice()){
        case 1:
            nextMenu = new AlgorithMenu();
            break;

        case 2:
            b.Start("totalOE");

            reco = new ReconstructionOE(pathToMeasurements,
                                        pathToSystemMatrix,
                                        {-52.5, 52.5, -52.5, 52.5, 5, 10});

            states = promptChoice("NUMBER OF STATES TO REACH EQUILIBRIUM: ");
            samples = promptChoice("NUMBER OF SAMPLES IN EQUILIBRIUM: ");
            reco->start(states, samples);

            b.Stop("totalOE");
            std::cout << "\nTotal Time:\t\t" << b.GetRealTime("totalOE") << " seconds\n";

            nextMenu = new MainMenu();
            break;

        default:
            break;
    }

    isQuitOptionSelected = false;
    delete reco;
    return nextMenu;
}
