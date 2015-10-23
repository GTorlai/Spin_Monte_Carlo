#include "IsingHamiltonian.cpp"
#include "params.h"
#include "measurements.cpp"

int main() {

	//Contruct parameter class
	Params par;

	//Contruct random number generator
    MTRand random(par.SEED);

	//Contruct hypercube lattice
	Hypercube cube(par.nX,par.Dim);
   
	//Contruct spin class
	Spins sigma;

	//Construct Ising model
	IsingHamiltonian ising(sigma,cube,random);	

	//Compute energy and magnetization
    ising.GetEnergy(sigma);
    ising.GetMagnetization(sigma);
    
	//Construct measurements class
    Measurements measure(sigma.N,par);
    
	//Create file for data
    //measure.createFileName(par.nX,par.Dim,"Ising_ferromagnet",par.MCS);
    measure.createFileName(par.nX,par.Dim,"Ising_random",par.MCS);
 
	ofstream fileData(measure.fileName.c_str());
    
	int counter = -1; 	//Progress counter
	double T;
    //for(T = par.Tlow; T<par.Thigh; T += par.Tstep) {
    for(double p = 0.05; p < 0.15; p += 0.01) {  
  		
		counter++;
		cout << "..." << counter/(par.Thigh-par.Tlow)/par.Tstep << "% progress " << endl;
		
		T = 2/log((1-p)/p);

		//Reset the observables values
		measure.reset();

		for(int r=0; r<par.ROD; r++) {
			
			//For the random bond Ising model
			ising.RandomizeInteractions(p,random);
			//ising.print();
			sigma.randomize(random);
			//sigma.print();
			ising.GetEnergy(sigma);
			ising.GetMagnetization(sigma);

			//Equilibrate the system
			for(int k=0; k<par.EQL; k++) {
				ising.LocalUpdate(sigma,T,random);
	        }//k

			//Run Markov chain Monte Carlo
	        for(int k=0; k<par.MCS; k++) {
	            ising.LocalUpdate(sigma,T,random);
	            //ising.GaugeUpdate(sigma,T,random);
				measure.record(ising.Energy,ising.Magn,sigma);
	        }//k
		
		}//ROD

		//Print measurements on file
        measure.output(T,fileData);
    
	}//T
    
	fileData.close();

    return 0;
}//main
