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
    virtual TemplateMenu* getNextMenu(bool& iIsQuitOptionSelected) = 0;
    virtual void show() const { std::cout << messageText; }

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

class AlgorithmMenu : public TemplateMenu{
public:
    AlgorithmMenu(TString iPathToMeasurements);
    TemplateMenu* getNextMenu(bool &iIsQuitOptionSelected);
    TString pathToMeasurements;
};

class MLEMMenu : public TemplateMenu{
public:
    MLEMMenu(TString iPathToMeasurements);
    TemplateMenu* getNextMenu(bool &iIsQuitOptionSelected);
    TString pathToMeasurements;
};

class MLEMRecoMenu : public TemplateMenu{
public:
    MLEMRecoMenu(TString iPathToMeasurements, TString iPathToProjections);
    TemplateMenu* getNextMenu(bool &iIsQuitOptionSelected);
    TString pathToMeasurements, pathToProjections;
};
