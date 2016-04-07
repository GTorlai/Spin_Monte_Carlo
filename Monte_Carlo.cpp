#include "IsingHamiltonian.cpp"
//#include "ToricCodeError.cpp"
#include "params.h"
#include "measurements.cpp"
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char *argv[]) {

    //Contruct parameter class
    Params par;
    
    int simNum = 0;
    double p = 0.0;
    double T = 0.0;
    char* model;
    int seed_add = 0;
    
    for (int i=1;i<argc;i++) {
        if (strcmp(argv[i],"--n") == 0) {        
            simNum = atoi(argv[i+1]);
        }//if
        if (strcmp(argv[i],"--seed") == 0) {
            seed_add = atoi(argv[i+1]);
        }//if
        if (strcmp(argv[i],"--p") == 0) {
            p = double(atof(argv[i+1]));
        }//if
        if (strcmp(argv[i],"--T") == 0) {
            T = double(atof(argv[i+1]));
        }//if
        if (strcmp(argv[i],"--L") == 0) {
            par.nX = atoi(argv[i+1]);
        }//if
        else if (strcmp(argv[i], "--model") == 0) {
            model = argv[i+1];
        }//else if
    }//i

    //Contruct random number generator
    MTRand random(par.SEED+seed_add);

    //Contruct hypercube lattice
    Hypercube cube(par.nX,par.Dim);
    
    //Contruct spin class
    Spins sigma;

    //Construct Ising model
    IsingHamiltonian ising(sigma,cube,random);
    //ising.print();
    
    
    //Compute energy and magnetization
    ising.GetEnergy(sigma);
    ising.GetMagnetization(sigma);
   
    //Construct measurements class
    Measurements measure(sigma.N,par);
    
    //Create file for data
    measure.createFileName(par,model,simNum);       
    ofstream fileData(measure.fileName.c_str());
    
    double T0 = 20.0;
    //ising.AnnealedImportanceSampling(sigma,random,T,T0);


    measure.printHeaders(fileData);
    for(int k=0; k<par.EQL; k++) {
        ising.LocalUpdate(sigma,T,random);
    }//k

    //for(T = par.Tlow; T<par.Thigh; T += par.Tstep) {  //ISING FERROMAGNET
    //for(p = 0.05; p < 0.15; p += 0.01) {                //ISING RANDOM
                
        cout << "Temperature: :" << T << endl;
        //cout << "Disorder strength: " << p << endl;

        //T = 2.0/log((1-p)/p);                       //ISING RANDOM

        //Reset the observables values
        //measure.reset();
        //for(int r=0; r<par.ROD; r++) {
            //ising.RandomizeInteractions(p,random);  //ISING RANDOM
            //sigma.randomize(random);                //ISING RANDOM
                        
            //ising.GetEnergy(sigma);                 //ISING RANDOM
            //ising.GetMagnetization(sigma);          //ISING RANDOM

            //ising.SimulatedAnnealing(sigma,T,random);
            
            //Equilibrate the system
            //for(int k=0; k<par.EQL; k++) {
            //    ising.LocalUpdate(sigma,T,random);
            //}//k

            //Run Markov chain Monte Carlo
            int counter =0;
            for(int k=0; k<par.MCS; k++) {
                
                if (counter == 100) {
                    ising.GaugeUpdate(sigma,T,random);
                    counter = 0;
                }
                ising.LocalUpdate(sigma,T,random);
                counter++;

                ising.GaugeUpdate(sigma,T,random);
                measure.record(ising.Energy,ising.Magn,sigma);
            }//k
            
            //sigma.filePrint(fileData);       
            //sigma.print();
        //}//ROD
                
        //Calculate correlation lengths
        //measure.GetCorrelationLength(par.nX,cube.Coordinates);
                
        //Print measurements on file
        measure.output(T,"Ising_Ferromagnet",fileData);    //ISING FERROMAGNET
        //measure.output(p,"Ising_Random",fileData);   //ISING RANDOM

    
    //}//T-p
    
    fileData.close();
    
    return 0;

}//main
