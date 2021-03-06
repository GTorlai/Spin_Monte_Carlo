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
    long int Magn;
    
    double Z;
    int K;
    int M;

    vector<double> beta;

    vector<int> LocalMagn;
    vector<vector<int> > SpinSpinCorr;

    vector<double> J;	//Vector of interaction couplings
    
    vector<vector<int> > NearestNeighbors;

    vector<vector<int> > bonds;
    
    //Functions
    IsingHamiltonian(Spins & sigma, Hypercube & cube, MTRand & random);

    void GetEnergy(Spins & sigma);
    void GetMagnetization(Spins & sigma);

    void RandomizeInteractions(double p, MTRand & random);
    void Update(Spins & sigma,double dE, int site);
    double aisSweep(Spins & sigma,MTRand & random,double T,double T0);

    void AnnealedImportanceSampling(Spins & sigma, MTRand & random,double T,double T0);
    void SimulatedAnnealing(Spins & sigma, double T_f, MTRand & random);
    void LocalUpdate(Spins & sigma, double & T, MTRand & random);
    void GrowCluster(Spins & sigma, double & T, MTRand & ran, int & site);
    void GaugeUpdate(Spins & sigma, double & T, MTRand & ran);
    
    

    void print();
};

#endif
