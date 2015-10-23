#include "measurements.h"

//Constructor
Measurements::Measurements(const int & N_, const Params & par){

    N = N_;
	MCS = par.MCS; //Include realization of disorder
	ROD = par.ROD;

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_Mag = 0.0;
    TOT_Mag2 = 0.0;
	TwoPointCorr = 0.0;
}

//Reset
void Measurements::reset(){

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_Mag = 0.0;
    TOT_Mag2 = 0.0;
	TwoPointCorr = 0.0;
}

//Update the measurements
void Measurements::record(double & energy, long int & magn, Spins & sigma){

    TOT_energy += energy;
    TOT_energy2 += energy * energy;
    TOT_Mag += 1.0*abs(magn);
    TOT_Mag2 += 1.0*magn*magn;
    TwoPointCorr += 1.0*sigma.spin[0]*sigma.spin[N/2];

}//update


//Write average on file
void Measurements::output(const double & T, ofstream & file){
    
    file<< T <<" ";
    file<< TOT_energy/(1.0*MCS * N*ROD) <<" ";
    file<< TOT_energy2/(1.0*MCS * N * N*ROD) <<" ";
    double Cv = TOT_energy2/(1.0*MCS*ROD) - TOT_energy*TOT_energy/(1.0*MCS*MCS*ROD*ROD);
    file<< Cv/(T*T*1.0*N) <<" ";
    file<< TOT_Mag/(1.0*MCS * N*ROD) <<" ";
    file<< TOT_Mag2/(1.0*MCS * N*N*ROD) <<" ";
    double susc = TOT_Mag2/(1.0*MCS*ROD) - TOT_Mag*TOT_Mag/(1.0*MCS*MCS*ROD*ROD);
    file<< susc/(T*1.0*N) <<" ";
	file<< TwoPointCorr/(1.0*MCS*ROD);
    file << endl;
    
}//output

//Create name-file
void Measurements::createFileName(int L_, int D_, const char* model, int MCS_) {
    
    fileName.clear();
    stringstream str;
    
    str << "MC_" << D_ << "D_" << model;
    str << "_L" << L_;
    str << "_MCS" << MCS_/1000;
    str << "k.dat";
    
    fileName = str.str();
}
