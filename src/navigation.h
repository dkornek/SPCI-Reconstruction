// navigation.h
// main app menu navigation

#pragma once
#include <iostream>
#include <string>

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

class SystemMatrixMenu : public TemplateMenu{
public:
    SystemMatrixMenu();
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
};

class MeasurementsMenu : public TemplateMenu{
public:
    MeasurementsMenu();
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
};

class AlgorithMenu : public TemplateMenu{
public:
    AlgorithMenu();
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
};

class MLEMMenu : public TemplateMenu{
public:
    MLEMMenu();
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
};

class OEMenu : public TemplateMenu{
public:
    OEMenu();
    TemplateMenu* getNextMenu(bool& isQuitOptionSelected);
};
