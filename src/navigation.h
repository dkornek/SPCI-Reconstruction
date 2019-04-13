// navigation.h
// main app menu navigation

#pragma once
#include <iostream>
#include <string>

//class Navigation{

//public:
//    std::string promptPath();
//    int checkInput(int choice, int lowerValue, int upperValue);

//    void algorithmMenu();
//};

class TemplateMenu{
public:
    TemplateMenu(){ messageText = "TemplateMenu"; }
    virtual ~TemplateMenu(){}
    virtual TemplateMenu* getNextMenu(bool& iIsQuitOptionSelected) = 0;
    virtual void show(){ std::cout << messageText; }
    void promptChoice(int& iChoice);
    std::string promptPath();

protected:
    std::string messageText;
};

class MainMenu : public TemplateMenu{
public:
    MainMenu();
    TemplateMenu* getNextMenu(bool& iIsQuitOptionSelected);
};

class ReconstructionMenu : public TemplateMenu{
public:
    ReconstructionMenu();
    TemplateMenu* getNextMenu(bool &iIsQuitOptionSelected);
};
