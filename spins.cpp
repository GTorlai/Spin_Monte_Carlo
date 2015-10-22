#include "spins.h"

//Constructor 1
Spins::Spins(){

    spin.clear(); 

}

//Constructor 2
Spins::Spins(int N_){

    N = N_;

    spin.resize(N,1); //assign every spin as 1

}

//takes the total number of lattice sites
void Spins::resize(int N_){

    N = N_;

    spin.resize(N,1); //assign every spin as 1

}

//Randomize the spin values
void Spins::randomize(){

    MTRand irand(129345); //random number 

    int ising_spin;
    
    for (int i = 0; i<spin.size(); i++){
        ising_spin = 2*irand.randInt(1)-1;

        spin.at(i) = ising_spin;
    }//i


}//randomize

//Perform a single-spin flip
void Spins::flip(int index){

    spin.at(index) *= -1;

}//flip


//Print function
void Spins::print(){

    for (int i=0;i<spin.size();i++){
        cout<<(spin[i]+1)/2<<" ";
    }//i
    cout<<endl;

}//print
