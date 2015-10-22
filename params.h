#ifndef PARAM_H
#define PARAM_H

#include <fstream>
#include <iostream>
using namespace std;

//Simple class to read in the simulation parameters from a file
class Params
{
    public:
        int nX;      //linear size of lattice
        int Dim;      //dimension of lattice
        double Tlow;      //Temperature limit 1
        double Thigh;      //Temperature limit 2
        double Tstep;      //Temperature step
        int EQL;     //the number of equilibration steps
        int MCS;     //the number of Monte Carlo steps
        int nBin;    //number of production bins
        long SEED;   //the random number seed

        Params();
        void print();


}; //PARAMS

Params::Params(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("params.dat");

    pfin >> nX;
    pfin >> Dim;
    pfin >> Tlow;
    pfin >> Thigh;
    pfin >> Tstep;
    pfin >> EQL;
    pfin >> MCS;
    pfin >> nBin;
    pfin >> SEED;
    pfin.close();

}//constructor


void Params::print(){

    cout<<"Linear size "<<nX<<endl;
    cout<<"Dimension "<<Dim<<endl;
    cout<<"# Equil steps "<<EQL<<endl;
    cout<<"# MC steps "<<MCS<<endl;
    cout<<"# data bins "<<nBin<<endl;
    cout<<"RNG seed "<<SEED<<endl;

}

#endif
