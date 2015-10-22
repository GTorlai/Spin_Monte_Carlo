#include "IsingHamiltonian.cpp"
#include "params.h"
#include "measurements.cpp"

int main() {

	Params par;

    MTRand random(par.SEED);

	Hypercube cube(par.nX,par.Dim);
   
	Spins sigma;

	IsingHamiltonian ising(sigma,cube,random);	

    ising.GetEnergy(sigma);
    ising.GetMagnetization(sigma);
    
    Measurements measure(sigma.N,par);
    
    measure.createFileName(par.nX,"Ising_ferromagnet",par.MCS);
    ofstream fileData(measure.fileName.c_str());
    
    for(double T = par.Tlow; T<par.Thigh; T += par.Tstep) {
        for(int k=0; k<par.EQL; k++) {
            ising.LocalUpdate(sigma,T,random);
        }//k
        
        measure.reset();
        
        for(int k=0; k<par.MCS; k++) {
            ising.LocalUpdate(sigma,T,random);
            measure.record(ising.Energy,ising.Magn,sigma);
        }//k
        measure.output(T,fileData);
        sigma.print();
    }//T
    fileData.close();

    return 0;
}//main
