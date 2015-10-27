#ifndef MEASURE_H
#define MEASURE_H

// measure.h: a class that performs statistical measurements of estimators

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
		double CorrLength;			//Correlation length using Sandvik def a
		//double CorrLength2;		//Correlation length using Sandvik def b
		double CorrLengthConnected;	//Connected correlation squared
		
		vector<vector<double> > SpinSpinCorr;	//Spin SPin correlation matrix
    	vector<double> LocalMagn;

    	string fileName;

    	Measurements(const int & N_, const Params & par);
    
    	void createFileName(Params & par, const char* model, int simNum);
    	void reset();
    	void printHeaders(ofstream & file);
		void record(double & energy, long int & magn, Spins & sigma);
    	void GetCorrelationLength(const int & L,vector<vector<int> > &coordinate);
		void output(const double simPar, const char* model, ofstream & file);
};

#endif
