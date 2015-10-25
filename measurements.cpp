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
	TOT_Mag4 = 0.0;
	TwoPointCorr = 0.0;

	LocalMagn.assign(N,0.0);
	
	for(int i=0; i<N; i++) {
		SpinSpinCorr.push_back(LocalMagn);
	}//i

}//Constructor

//Reset
void Measurements::reset(){

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_Mag = 0.0;
    TOT_Mag2 = 0.0;
	TOT_Mag4 = 0.0;
	TwoPointCorr = 0.0;

	for(int i=0; i<N; i++) {
		LocalMagn[i] = 0.0;
		for(int j=0; j<N; j++) {
			SpinSpinCorr[i][j] = 0.0;
		}//j
	}//i

}//reset

//Update the measurements
void Measurements::record(double & energy, long int & magn, Spins & sigma){

    TOT_energy += energy;
    TOT_energy2 += energy * energy;
    TOT_Mag += 1.0*abs(magn);
    TOT_Mag2 += 1.0*magn*magn;
	TOT_Mag4 += 1.0*magn*magn*magn*magn;

    TwoPointCorr += 1.0*sigma.spin[0]*sigma.spin[N/2];
 
 	for(int i=0; i<N; i++) {
		LocalMagn[i] += sigma.spin[i];
		for(int j=0; j<N; j++) {
			SpinSpinCorr[i][j] += sigma.spin[i]*sigma.spin[j];
		}//j
	}//i

}//update

//Computer different correlation lengths
void Measurements::GetCorrelationLength(const int & L, vector<vector<int> > &coordinate) {	

	double q1 = 2*PI/L;
	double q2 = 4*PI/L;	
	double suscept0 =0.0;
	double suscept1 =0.0;

	vector<vector<double> > ConnectedCorr;
	ConnectedCorr.resize(N,vector<double>(N));
	
	//Get Average of correlation and local magnetization
	for(int i=0; i<N; i++) {
		LocalMagn[i] /= 1.0*MCS*ROD;
		for(int j=0; j<N; j++) {
			ConnectedCorr[i][j] = SpinSpinCorr[i][j]/(1.0*MCS*ROD);
		}//i
	}//j

	//CORRELATION LENGTH A
	for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
			suscept0 += ConnectedCorr[i][j];
			suscept1 += ConnectedCorr[i][j]*cos(q1*(coordinate[i][0]-coordinate[j][0])); 	
		}//i
	}//j

	CorrLengthA = (1.0/q1)*sqrt(suscept0/suscept1 - 1.0)/L;
	
	//CORRELATION LENGTH B
	suscept0 = 0.0;
	suscept1 = 0.0;

	for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            suscept0 += ConnectedCorr[i][j]*cos(q2*(coordinate[i][0]-coordinate[j][0]));	
			suscept1 += ConnectedCorr[i][j]*cos(q1*(coordinate[i][0]-coordinate[j][0]));//*(ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j]);
        }//i
    }//j

	CorrLengthB = (1.0/q1)*sqrt((suscept1/suscept0 - 1.0)/(4.0 - suscept1/suscept0))/L;

	//CORRELATION LENGTH C
	suscept0 = 0.0;
	suscept1 = 0.0;

	for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            suscept0 += (ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j]);//*(ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j]);	
			suscept1 += (ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j])*cos(q1*(coordinate[i][0]-coordinate[j][0]));//*(ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j]);
        }//i
    }//j

	CorrLengthC = (1.0/(2.0*sin(q1/2.0)))*sqrt(suscept0/suscept1 - 1.0)/L;

	//CORRELATION LENGTH D
	suscept0 = 0.0;
	suscept1 = 0.0;

	for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            suscept0 += (ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j])*(ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j]);	
			suscept1 += (ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j])*cos(q1*(coordinate[i][0]-coordinate[j][0]))*(ConnectedCorr[i][j]-LocalMagn[i]*LocalMagn[j]);
        }//i
    }//j

	CorrLengthD = (1.0/(2.0*sin(q1/2.0)))*sqrt(suscept0/suscept1 - 1.0)/L;

}//GetCorrelationLength


//Write average on file
void Measurements::output(const double & T, ofstream & file){
    
    file<< T <<" ";
    file<< TOT_energy/(1.0*MCS * N*ROD) <<" ";
    //file<< TOT_energy2/(1.0*MCS * N * N*ROD) <<" ";
    double Cv = TOT_energy2/(1.0*MCS*ROD) - TOT_energy*TOT_energy/(1.0*MCS*MCS*ROD*ROD);
    file<< Cv/(T*T*1.0*N) <<" ";
    file<< TOT_Mag/(1.0*MCS * N*ROD) <<" ";
    //file<< TOT_Mag2/(1.0*MCS * N*N*ROD) <<" ";
    double susc = TOT_Mag2/(1.0*MCS*ROD) - TOT_Mag*TOT_Mag/(1.0*MCS*MCS*ROD*ROD);
    file<< susc/(T*1.0*N) <<" ";
	file<< TwoPointCorr/(1.0*MCS*ROD)<<" ";
    double AndersonOrder = 0.0;
 	for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
			AndersonOrder += SpinSpinCorr[i][j];
		}//j
	}//i
	file<< AndersonOrder/(1.0*MCS*ROD*N*N) <<" ";
	file<< CorrLengthA << " ";
	file<< CorrLengthB << " ";
	file<< CorrLengthC << " ";
	file<< CorrLengthD << " ";
	double BinderCumulant = 1.5*(1.0-(1.0/3.0)*(TOT_Mag4/(TOT_Mag2*TOT_Mag2))*ROD*MCS);
	file<< BinderCumulant;
	file << endl;
    
}//output

//Create name-file
void Measurements::createFileName(Params & par,  const char* model) {
    
    fileName.clear();
    stringstream str;
    
    str << "MC_" << par.Dim << "D_" << model;
    str << "_L" << par.nX;
	if(par.ROD > 1) 
		str << "_ROD" << par.ROD;
    str << "_MCS" << par.MCS/1000;
    str << "k.dat";
    
    fileName = str.str();
}
