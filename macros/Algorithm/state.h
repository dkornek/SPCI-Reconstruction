// state.h

#pragma once
#include <array>
#include <vector>
#include <random>

#include <TH3.h>
#include <TList.h>

// ##### STATE CHARACTERIZATION #####

class State{
public:
    State(const TH3F* N_dcb);
    ~State(){}

    void generateRandomOrigins(const Int_t numberOfVoxels, const TList* systemMatrix);
    void MCMCNextState(const TList* systemMatrix);

    // ##### MEMBERS #####
    std::vector<Int_t> events;                  // events n from 1 ... N
    std::vector<Int_t> origins;                 // origins (=voxels) v of event n
    std::vector<std::array<Int_t, 3>> bin;      // bin (d/c/bin) of event n

    std::vector<Int_t> countsInVoxel;           // counts C_sv in voxel v

private:
    Double_t calculateTransitionProbability(const TList* systemMatrix);
};

State::State(const TH3F *N_dcb){
    // fill the events- and bin-vector

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

void State::generateRandomOrigins(const Int_t numberOfVoxels, const TList *systemMatrix){
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

    for (UInt_t n = 0; n < events.size(); ++n){

        Int_t origin;
        while (true){
            // system matrix element must be > 0

            origin = Int_t(uniDisVox(generate));
            TH3F* p_dcb = (TH3F*)systemMatrix->At(origin);
            Double_t probability = p_dcb->GetBinContent(bin.at(n)[0], bin.at(n)[1], bin.at(n)[2]);

            if (probability > 0){
                break;
            }
        }

        origins.push_back(origin);
        ++countsInVoxel[origin];
    }
}

void State::MCMCNextState(const TList *systemMatrix){
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

    // ##### NEXT STATE #####
    for (UInt_t n = 0; n < numberOfEvents; ++n){

        // randomly select event n
        Int_t randomEvent = Int_t(uniDisEv(generate));

        // origin v of random event n
        Int_t originOfRandomEvent = origins.at(randomEvent);

        // randomly select new origin v' for event n
        Int_t nextOrigin = Int_t(uniDisVox(generate));

        if (nextOrigin == originOfRandomEvent){
            // if the origins are the same, continue the loop
            continue;
        }

        // calculate transition probability
        Double_t transitionProbability = calculateTransitionProbability(systemMatrix);

        // draw random chance of transition
        Double_t chance = uniDisTrans(generate);

        // move origin to new origin if successful
        if (chance <= transitionProbability){

            origins.at(randomEvent) = nextOrigin;
            --countsInVoxel[originOfRandomEvent];
            ++countsInVoxel[nextOrigin];
        }
    }
}

Double_t State::calculateTransitionProbability(const TList *systemMatrix){

    Double_t probability;

    // stats of previous state
    Double_t p_dcbvOld = 1;
    Double_t C_svOld = 1;
    Double_t p_vOld = 1;  // (sensitivity)

    // stats of next state
    Double_t p_dcbvNew = 1;
    Double C_svNew = 1;
    Double_t p_vNew = 1;

    // numerator

    // denominator

    // ratio

//    return probability;
    return 0.0;
}
