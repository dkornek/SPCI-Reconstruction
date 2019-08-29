// state.h

#pragma once
#include <array>
#include <vector>
#include <TH3.h>

// ##### STATE CHARACTERIZATION #####
class State{
public:
    State(const TH3F* N_dcb);
    ~State(){}

    std::vector<Int_t> generateRandomOrigins(const Int_t numberOfVoxels,
                                             const std::vector<std::vector<std::vector<std::vector<Double_t> > > > systemMatrix);
    std::vector<Int_t> MCMCNextState(const std::vector<std::vector<std::vector<std::vector<Double_t> > > > systemMatrix,
                                     const std::vector<Double_t> sensitivities,
                                     Double_t& relTransitions);

    // ##### MEMBERS #####
    std::vector<Int_t> events;                  // events n from 1 ... N
    std::vector<Int_t> origins;                 // origins (=voxels) v of event n
    std::vector<std::array<Int_t, 3>> bin;      // detector element/bin (d/c/b) of event n

    std::vector<Int_t> countsInVoxel;           // counts C_sv in voxel v
};
