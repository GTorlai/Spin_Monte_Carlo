#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include <vector>
#include <iostream>

using namespace std;

class Hypercube {
    
	public:
        int L; //linear size
        int D; //dimension
        int N; //total number of sites

        //the lattice is a vector of vectors: no double counting
        vector<vector<int> > Neighbors; //neighbors with positive offsets
        vector<vector<int> > Negatives; //neighbors with negative offsets
        vector<vector<int> > Coordinates;

        //public functions
        Hypercube(int L_, int D_);
        void print();
         
    private:
        int myPow(int, int);

};

#endif
