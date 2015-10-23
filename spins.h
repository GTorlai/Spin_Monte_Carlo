#ifndef SPINS_H
#define SPINS_H

// spins.h
// Class that contains a vector of lattice Ising spins

#include <vector>
#include "MersenneTwister.h"

using namespace std;

class Spins {
    
    public:
    
        int N;  //Number of spins


        //Vector containing the spins values
        vector<int> spin;

        //Functions
        Spins(int N_);
        Spins();
        void resize(int N_);
        void flip(int index);
        void print();
        void randomize(MTRand & random);
    
};

#endif
