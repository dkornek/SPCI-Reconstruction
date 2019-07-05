// state.h

#pragma once
#include <array>
#include <vector>
#include <random>
#include <algorithm>

#include <TH3.h>
#include <TList.h>

// ##### STATE CHARACTERIZATION #####

class State{
public:
    State(const TH3F* N_dcb);
    ~State(){}

    void generateRandomOrigins(const Int_t numberOfVoxels,
                               const std::vector<std::vector<std::vector<std::vector<Double_t> > > > systemMatrix);
    void MCMCNextState(const std::vector<std::vector<std::vector<std::vector<Double_t> > > > systemMatrix,
                       const std::vector<Double_t> sensitivities);

    // ##### MEMBERS #####
    std::vector<Int_t> events;                  // events n from 1 ... N
    std::vector<Int_t> origins;                 // origins (=voxels) v of event n
    std::vector<std::array<Int_t, 3>> bin;      // bin (d/c/bin) of event n

    std::vector<Int_t> countsInVoxel;           // counts C_sv in voxel v
};

State::State(const TH3F *N_dcb){
    // fill the events- and bin-vector in pseudo-list-mode format

    Int_t detD = N_dcb->GetNbinsX();
    Int_t detC = N_dcb->GetNbinsY();
    Int_t nBins = N_dcb->GetNbinsZ();

    Int_t numberOfEvents = 1;
    for (Int_t d = 1; d <= detD; ++d){
        for (Int_t c = 1; c <= detC; ++c){
            for (Int_t b = 1; b <= nBins; ++b){

                std::array<Int_t, 3> dcb = { d, c, b };
                Int_t numberOfEventsInBin = Int_t(N_dcb->GetBinContent(d, c, b));
                for (Int_t n = 0; n < numberOfEventsInBin; ++n){

                    events.push_back(numberOfEvents);
                    bin.push_back(dcb);
                    ++numberOfEvents;
                }
            }
        }
    }
}

void State::generateRandomOrigins(const Int_t numberOfVoxels,
                                  const std::vector<std::vector<std::vector<std::vector<Double_t> > > > systemMatrix){
    // generate random origins for each event
    // function is used to generate inital state s_0 for OE algorithm

    // fill countsInVoxel-Vector
    for (Int_t v = 0; v < numberOfVoxels; ++v){
        countsInVoxel.push_back(0);
    }

    // generate random seed
    std::random_device seed;

    // Mersenne-Twister random number engine
    std::mt19937 generate(seed());

    // uniform distribution
    std::uniform_real_distribution<Double_t> uniDisVox(0, numberOfVoxels);

    // accelerate by drawing random numbers in advance
    std::vector<Int_t> randomOrigins;
    for (UInt_t n = 0; n < events.size(); ++n){
        Int_t origin = Int_t(uniDisVox(generate));
        randomOrigins.push_back(origin);
    }

    for (UInt_t n = 0; n < events.size(); ++n){
        // draw random origin

        Int_t origin = randomOrigins[n];
        while (true){
            // system matrix element must be > 0

            std::array<Int_t, 3> b = bin[n];
            Double_t probability = systemMatrix.at(origin).at(b[0] - 1).at(b[1] - 1).at(b[2] - 1);

            if (probability > 0){
                break;
            }

            origin = Int_t(uniDisVox(generate));
        }

        origins.push_back(origin);
        ++countsInVoxel[origin];
    }
}

void State::MCMCNextState(const std::vector<std::vector<std::vector<std::vector<Double_t> > > > systemMatrix,
                          const std::vector<Double_t> sensitivities){
    // generate new state in Marcov Chain Monte Carlo manner

    // ##### PREPARATION #####
    // generate random seed
    std::random_device seed;

    // Mersenne-Twister random number engine
    std::mt19937 generate(seed());

    // uniform distribution for events
    UInt_t numberOfEvents = events.size();
    std::uniform_real_distribution<Double_t> uniDisEv(0, numberOfEvents);

    // uniform distribution for voxels
    UInt_t numberOfVoxels = countsInVoxel.size();
    std::uniform_real_distribution<Double_t> uniDisVox(0, numberOfVoxels);

    // uniform distribution for transition of state
    std::uniform_real_distribution<Double_t> uniDisTrans(0, 1);

    // ##### GENERATE RANDOM NUMBERS #####
    // accelerate by drawing random numbers in advance
    std::vector<Int_t> randomEvents;
    std::vector<Int_t> randomOrigins;
    std::vector<Double_t> transitionChances;
    for (UInt_t n = 0; n < numberOfEvents; ++n){
        Int_t randomEvent = Int_t(uniDisEv(generate));
        randomEvents.push_back(randomEvent);

        Int_t randomOrigin = Int_t(uniDisVox(generate));
        randomOrigins.push_back(randomOrigin);

        Double_t transitionChance = uniDisTrans(generate);
        transitionChances.push_back(transitionChance);
    }

    // ##### NEXT STATE #####
    Double_t successfulTransitions = 0;
    for (UInt_t n = 0; n < numberOfEvents; ++n){

        // randomly select event n
        Int_t randomEvent = randomEvents[n];

        // origin v of random event n
        Int_t originFrom = origins[randomEvent];

        // randomly select new origin v' for event n
        Int_t originTo = randomOrigins[n];

        if (originTo == originFrom){
            // if the origins are the same, continue the loop
            continue;
        }

        // calculate transition probability
        std::array<Int_t, 3> b = bin.at(randomEvent);
        Double_t p_dcbvFrom = systemMatrix.at(originFrom).at(b[0] - 1).at(b[1] - 1).at(b[2] - 1);
        Double_t C_svFrom = countsInVoxel[originFrom];
        Double_t sensitivityFrom = sensitivities.at(originFrom);

        Double_t p_dcbvTo = systemMatrix.at(originTo).at(b[0] - 1).at(b[1] - 1).at(b[2] - 1);
        Double_t C_svTo = countsInVoxel[originTo];
        Double_t sensitivityTo = sensitivities.at(originTo);

        // Sitek, Statistical Computing in Nuclear Imaging, 2015, Equation (6.32)
         Double_t ratio = (p_dcbvTo * (C_svTo + 1) * sensitivityFrom) / (p_dcbvFrom * C_svFrom * sensitivityTo);

        // Sitek, Statistical Computing in Nuclear Imaging, 2015, Equation (6.39)
        // Double_t ratio = ((p_dcbvTo * sensitivityFrom) / (p_dcbvFrom * sensitivityTo)) * std::pow((C_svTo + 1) / C_svFrom, 0.5) ;

        // Sitek, 2011, Eq. (7)
        // Double_t ratio = (p_dcbvTo / p_dcbvFrom) * (sensitivityFrom / sensitivityTo) * ((C_svTo + 1) / (C_svFrom - 1))
        //         * std::pow((C_svFrom - 1) / C_svFrom, C_svFrom)
        //         * std::pow((C_svTo + 1) / C_svTo, C_svTo);

        Double_t transitionProbability = std::min(1.0, ratio);

        // draw random chance of transition
        Double_t chance = transitionChances[n];

        // move origin to new origin if successful
        if (chance <= transitionProbability){

            origins.at(randomEvent) = originTo;
            --countsInVoxel[originFrom];
            ++countsInVoxel[originTo];
            ++successfulTransitions;
        }
    }

    std::cout << successfulTransitions/numberOfEvents << "\n";
}
