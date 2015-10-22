#ifndef PARAM_H
#define PARAM_H

#include <fstream>
#include <iostream>
using namespace std;

//Simple class to read in the simulation parameters from a file
class Params
{
    public:
        int nX_;      //linear size of lattice
        int Dim_;      //dimension of lattice
        double Tlow_;      //Temperature limit 1
        double Thigh_;      //Temperature limit 2
        double Tstep_;      //Temperature step
        int EQL_;     //the number of equilibration steps
        int MCS_;     //the number of Monte Carlo steps
        int nBin_;    //number of production bins
        long SEED_;   //the random number seed

        Params();
        void print();


}; //PARAMS

Params::Params(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("params.dat");

    pfin >> nX_;
    pfin >> Dim_;
    pfin >> Tlow_;
    pfin >> Thigh_;
    pfin >> Tstep_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> nBin_;
    pfin >> SEED_;
    pfin.close();

}//constructor


void Params::print(){

    cout<<"Linear size "<<nX_<<endl;
    cout<<"Dimension "<<Dim_<<endl;
    cout<<"# Equil steps "<<EQL_<<endl;
    cout<<"# MC steps "<<MCS_<<endl;
    cout<<"# data bins "<<nBin_<<endl;
    cout<<"RNG seed "<<SEED_<<endl;

}

#endif
