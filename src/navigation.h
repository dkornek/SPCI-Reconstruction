// navigation.h
// main app menu navigation

#pragma once
#include <iostream>
#include <string>
#include <TString.h>


class TemplateMenu{
public:
    TemplateMenu(){ messageText = "TemplateMenu"; }
    virtual ~TemplateMenu(){}
    virtual TemplateMenu* getNextMenu(bool& isQuitOptionSelected) = 0;
    virtual void show() const { std::cout << messageText; }

protected:
    std::string messageText;
};

class MainMenu : public TemplateMenu{
public:
    MainMenu();
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
};

class ReconstructionMenu : public TemplateMenu{
public:
    ReconstructionMenu();
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
};

class AlgorithmMenu : public TemplateMenu{
public:
    AlgorithmMenu(TString filePathToMeasurements, TString filePathToProjections);
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
    TString pathToMeasurements;
    TString pathToProjections;
};

class MLEMMenu : public TemplateMenu{
public:
    MLEMMenu(TString filePathToMeasurements, TString filePathToProjections);
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
    TString pathToMeasurements;
    TString pathToProjections;
};
