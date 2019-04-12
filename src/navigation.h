// navigation.h
// main app menu navigation

#pragma once
#include <string>

class Navigation{

public:
    int promptInput();
    std::string promptPath();
    int checkInput(int choice, int lowerValue, int upperValue);

    void mainMenu();
    void measurementsMenu();
    void algorithmMenu();
};
