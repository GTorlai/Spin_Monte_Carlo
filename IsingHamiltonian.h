#ifndef ISINGHAMILTONIAN_H
#define ISINGHAMILTONIAN_H

#include "hypercube.cpp"
#include "spins.cpp"

//IsingHamiltonian.h
//A class that simulates the classical Ising model

#include <string>

class IsingHamiltonian {

public:

	int N;	//Number of lattice sites		
	int D;  //Dimension

	double Energy;	//Total energy of the system

	vector<double> J;	//Vector of interaction couplings

	vector<vector<int> > NearestNeighbors;

	vector<vector<int> > bonds;

	//Functions
	IsingHamiltonian(Spins & sigma, Hypercube & cube, MTRand & random);

	void GetEnergy(Spins & sigma);
	void RandomizeInteractions(double p, MTRand & random);
	void LocalUpdate(Spins & sigma, double & T, MTRand & random);
	void print();
};

#endif
