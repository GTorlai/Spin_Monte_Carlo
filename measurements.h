#ifndef MEASURE_H
#define MEASURE_H

// measure.h: a class that performs statistical measurements of estimators
// Energy, Specific Heat, Magnetization, Susceptibility

#include "spins.h"
#include "params.h"
#include <sstream>
#include <fstream>
#include <string>

#define PI 3.14159265

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
	double TOT_Mag4;			//Fourth power of magnetization

	double TwoPointCorr;		//Spin-Spin correlation at opposite side
	double CorrLengthA;			//Correlation length using Sandvik def a
	double CorrLengthB;			//Correlation length using Sandvik def b
	double CorrLengthC;			//Connected correlation 
	double CorrLengthD;			//Connected correlation squared

	vector<vector<double> > SpinSpinCorr;	//Spin SPin correlation matrix
    vector<double> LocalMagn;

    string fileName;

    Measurements(const int & N_, const Params & par);
    
    void createFileName(Params & par, const char* model, double simPar);
    void reset();
    void printHeaders(ofstream & file);
	void record(double & energy, long int & magn, Spins & sigma);
    void GetCorrelationLength(const int & L,vector<vector<int> > &coordinate);
	void output(const double & T, ofstream & file);
 
};

#endif
