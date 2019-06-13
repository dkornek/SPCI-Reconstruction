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

    void generateRandomOrigins(const Int_t numberOfVoxels, const TList* systemMatrix);
    void MCMCNextState(const TList* systemMatrix, const std::vector<Double_t> sensitivies);

    // ##### MEMBERS #####
    std::vector<Int_t> events;                  // events n from 1 ... N
    std::vector<Int_t> origins;                 // origins (=voxels) v of event n
    std::vector<std::array<Int_t, 3>> bin;      // bin (d/c/bin) of event n

    std::vector<Int_t> countsInVoxel;           // counts C_sv in voxel v

private:
    Double_t calculateTransitionProbability(const TList* systemMatrix, const std::vector<Double_t> sensitivies,
                                            const Int_t numberOfEvent,
                                            const Int_t originOld, const Int_t originNew);
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

void State::MCMCNextState(const TList *systemMatrix, const std::vector<Double_t> sensitivies){
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
    Double_t successfulTransitions = 0;
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
        Double_t transitionProbability = calculateTransitionProbability(systemMatrix, sensitivies,
                                                                        randomEvent,
                                                                        originOfRandomEvent, nextOrigin);

        // draw random chance of transition
        Double_t chance = uniDisTrans(generate);

        // move origin to new origin if successful
        if (chance <= transitionProbability){

            origins.at(randomEvent) = nextOrigin;
            --countsInVoxel[originOfRandomEvent];
            ++countsInVoxel[nextOrigin];
            ++successfulTransitions;
        }
    }

    std::cout << successfulTransitions/numberOfEvents << "\n";
}

Double_t State::calculateTransitionProbability(const TList *systemMatrix, const std::vector<Double_t> sensitivies,
                                               const Int_t numberOfEvent,
                                               const Int_t originOld, const Int_t originNew){
    // calculate transition probability (see Eq. (7) in A. SITEK, 2011)

    // PREVIOUS STATE
    // system matrix probability
    TH3F* p_dcbOld = (TH3F*)systemMatrix->At(originOld);
    Double_t p_dcbvOld = p_dcbOld->GetBinContent(bin.at(numberOfEvent)[0],
                                                 bin.at(numberOfEvent)[1],
                                                 bin.at(numberOfEvent)[2]);

    // counts in voxel
    Double_t C_svOld = countsInVoxel.at(originOld);

    // sensitivity
    Double_t p_vOld = sensitivies.at(originOld);

    // NEXT STATE
    // system matrix probability
    TH3F* p_dcbNew = (TH3F*)systemMatrix->At(originNew);
    Double_t p_dcbvNew = p_dcbNew->GetBinContent(bin.at(numberOfEvent)[0],
                                                 bin.at(numberOfEvent)[1],
                                                 bin.at(numberOfEvent)[2]);

    // counts in voxel
    Double_t C_svNew = countsInVoxel.at(originNew);

    // sensitivity
    Double_t p_vNew = sensitivies.at(originNew);

    // CALCULATE THE PROBABILITY
    Double_t lambda = (p_dcbvNew / p_dcbvOld) * (p_vOld / p_vNew);

    Double_t kappa = (C_svNew + 1) / (C_svOld - 1);

    // countsRatioOld/-New are nearly 1 and can thus be neglected (losing some precision)
//    Double_t countsRatioOld = (C_svOld - 1) / C_svOld;
//    countsRatioOld = std::pow(countsRatioOld, C_svOld);

//    Double_t countsRatioNew = (C_svNew + 1) / C_svNew;
//    countsRatioNew = std::pow(countsRatioNew, C_svNew);

//    Double_t probability = lambda * kappa * countsRatioNew * countsRatioOld;

    Double_t probability = lambda * kappa;
    probability = std::min(1.0, probability);

    return probability;
}
