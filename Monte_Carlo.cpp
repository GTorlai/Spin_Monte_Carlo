#include "IsingHamiltonian.cpp"
#include "param.h"


int main() {

	Params par;

	par.print();

    MTRand random;

	Hypercube cube(3,2);
   
	Spins sigma;

	IsingHamiltonian ising(sigma,cube,random);	
	//ising.RandomizeInteractions(0.5,random);	
	ising.print();
	ising.GetEnergy(sigma);

}
