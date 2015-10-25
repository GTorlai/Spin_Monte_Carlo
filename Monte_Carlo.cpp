#include "IsingHamiltonian.cpp"
#include "params.h"
#include "measurements.cpp"

int main(int argc, char *argv[]) {

    double p;
    double T;
    char* model;
    
    for (int i=1;i<argc;i++) {
        if (strcmp(argv[i],"-p") == 0) {
            p = double(atof(argv[i+1]));
        }//if
        if (strcmp(argv[i],"-T") == 0) {
            T = double(atof(argv[i+1]));
        }//if
        else if (strcmp(argv[i], "-model") == 0) {
            model = argv[i+1];
        }//else if
    }//i
    
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
    measure.createFileName(par,model);
 
	ofstream fileData(measure.fileName.c_str());
    
	//int counter = -1; 	//Progress counter
    //for(T = par.Tlow; T<par.Thigh; T += par.Tstep) {
    //for(p = 0.05; p < 0.15; p += 0.01) {
  		
		//counter++;
		//cout << "..." << counter/(par.Thigh-par.Tlow)/par.Tstep << "% progress " << endl;
		
		cout << "Temperature: :" << T << endl;
		//cout << "Disorder strength: " << p << endl;
		//T = 1.0/log((1-p)/p);

		//Reset the observables values
		measure.reset();

		for(int r=0; r<par.ROD; r++) {
			
			//RANDOM BOND ISING MODEL ONLY
			//ising.RandomizeInteractions(p,random);
			//sigma.randomize(random);
			//ising.GetEnergy(sigma);
			//ising.GetMagnetization(sigma);


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
		
		//Calculate correlation lengths
		measure.GetCorrelationLength(par.nX,cube.Coordinates);
		
		//Print measurements on file
        measure.output(T,fileData);
		sigma.print();
    
	//}//T-p
    
	fileData.close();

    return 0;
}//main
