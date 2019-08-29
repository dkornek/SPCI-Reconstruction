// user_prompts.h

#pragma once
#include <boost/lexical_cast.hpp>
#include <TFile.h>
#include <TString.h>

int promptChoice(std::string promptString = "SELECT: "){
   // prompt user for choosing an option

   std::string input;
   int choice;

   for (;;){
       std::cout << promptString;
       std::cin >> input;  // get user input

       try{
           choice = boost::lexical_cast<int>(input);
           break;  // success

       } catch(const boost::bad_lexical_cast &e){
           std::cerr << e.what() << "\n";  // print error message
       }
   }

   return choice;
}

bool checkPath(TString path){
   // check validity of path to .root file

   TFile* file = new TFile(path, "READ");

   if (file->IsOpen()){
       std::cout << "File opened successfully.\n";
       file->Close();
       return true;

   } else{
       std::cout << "Invalid file path.\n";
       return false;
   }
}

TString promptPath(std::string promptString){
   // prompt user for a file path

   TString input;

   for (;;){
       std::cout << promptString;
       std::cin >> input;  // get user input

       if (checkPath(input)){
           break;  // success
       }
   }

   return input;
}

double promptParameter(std::string promptString, double lowerLimit, double upperLimit){
    // prompt user for choosing any needed parameter

    std::string input;
    double choice;

    for (;;){
       std::cout << promptString;
       std::cin >> input;  // get user input

       try{
           choice = boost::lexical_cast<double>(input);

           if ((choice < lowerLimit) || (choice >= upperLimit)){
               continue;

           } else{
               break;  // success
           }

       } catch(const boost::bad_lexical_cast &e){
           std::cerr << e.what() << "\n";  // print error message
       }
    }

    return choice;
}
