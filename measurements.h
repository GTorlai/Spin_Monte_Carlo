#ifndef MEASURE_H
#define MEASURE_H

// measure.h: a class that performs statistical measurements of estimators
// Energy, Specific Heat, Magnetization, Susceptibility

#include "spins.h"
#include "params.h"
#include <sstream>
#include <fstream>
#include <string>

class Measurements {
private:
    
    int N;
    int MCS;
    int ROD;
    
public:

    double TOT_energy;   		//energy
    double TOT_energy2;  		//energy^2
    double TOT_Mag;    			//magnetization s
    double TOT_Mag2;    		//magnetization squared
	double TwoPointCorr;	//Spin-Spin correlation at opposite side

    string fileName;

    Measurements(const int & N_, const Params & par);
    
    void createFileName(int L_, int D_, const char* model, int MCS_);
    void reset();
    void record(double & energy, long int & magn, Spins & sigma);
    void output(const double & T, ofstream & file);
 
};

#endif
