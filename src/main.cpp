/*
This project executes the image reconstruction for SPCI measurement data.
To do so, two different algorithms can be used:
- maximum likelihood expectation maximization (ML-EM),
- origin ensemble (OE).
*/

#include "navigation.h"

int main(){
    // start application
    TemplateMenu* currentMenu;
    currentMenu = new MainMenu;

    bool isQuitOptionSelected;
    isQuitOptionSelected = false;  // shut down application if true

    while (!isQuitOptionSelected){
        currentMenu->show();
        TemplateMenu* nextMenu = currentMenu->getNextMenu(isQuitOptionSelected);

        if (nextMenu){
            delete currentMenu;
            currentMenu = nextMenu;
        }
    }

    return 0;
}
