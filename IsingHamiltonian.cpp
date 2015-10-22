#include "IsingHamiltonian.h"
#include <vector>

//Constructor
IsingHamiltonian::IsingHamiltonian(Spins & sigma, Hypercube & cube, MTRand & random) {

    N = cube.N;
	D = cube.D;

	sigma.resize(N);
	sigma.randomize();

	J.assign(N*D,1.0);

	//Build the nearest neighbors connections
	NearestNeighbors.resize(N,vector<int>(2*D));
	
	for(int i=0; i<NearestNeighbors.size(); i++) {
		for(int j=0; j<D; j++) {
			NearestNeighbors[i][j] = cube.Neighbors[i][j];
			NearestNeighbors[i][j+D] = cube.Negatives[i][j];
			//NearestNeighbors[cube.Neighbors[i][j]][j+cube.D] = i;			
		}//j
	}//i

	//The one cells connected to each zero-cell.  Used for cluster updates
	bonds.resize(N,vector<int>(2*D));
    	
    int back_site;
	
	for (int i=0; i<N; i++) {
		for (int j=0; j<D; j++){
			bonds[i][j]= D*i + j;
			back_site = cube.Negatives[i][j];
			bonds[i][j+D] = D*back_site + j;
        }//j
	}//i

}//constructor

//Initialize the interaction parameter according to model definition
void IsingHamiltonian::RandomizeInteractions(double p, MTRand & random) {
	
	for(int i=0; i<J.size(); i++) {
		if(random.rand() > p)
			J[i] = 1.0;
		else
			J[i] = -1.0;
	}//i

}//RandomizeInteractions


//Computer the energy in the system
void IsingHamiltonian::GetEnergy(Spins & sigma) {
	
	Energy=0.0;
	for(int i=0; i<sigma.N; i++) {
		for(int j=0; j<NearestNeighbors[i].size(); j++) {
			Energy += - sigma.spin[i]*sigma.spin[NearestNeighbors[i][j]]*J[bonds[i][j]];
		}//j
	}//i

	Energy /= 2.0;
    cout << "Initial Energy: " << Energy << endl;

}//GetEnergy

//Computer the magnetization in the system
void IsingHamiltonian::GetMagnetization(Spins & sigma) {
    
    Magn=0;
    for(int i=0; i<sigma.N; i++) {
        Magn += sigma.spin[i];
    }//i
    
    cout << "Initial Magnetization: " << Magn << endl;

}//GetMagnetization


//Perform Metropolis Update
void IsingHamiltonian::LocalUpdate(Spins & sigma, double & T, MTRand & random) {

	int site;
	double deltaE;
	double ran_num;

	for(int k=0; k<N; k++) {
        site = random.randInt(N-1);
		deltaE = 0.0;
		
		for(int i=0; i<NearestNeighbors[site].size(); i++) {
			deltaE += 2*sigma.spin[site]*sigma.spin[NearestNeighbors[site][i]]*J[bonds[site][i]];
		}//i
		//cout << "Energy difference: " << deltaE << endl;

		//Metropoli Algorithm

		if(deltaE<0) {
			sigma.flip(site);
			Energy += deltaE;
            Magn += 2*sigma.spin[site];
		}//if
		
		else {
			ran_num = random.rand();
			//cout << "Metropolis weight: " << exp(-deltaE/T) << "\t Random:" << ran_num << endl;
		   if (exp(-deltaE/T) > ran_num) {
		   		sigma.flip(site);
				Energy += deltaE;
                Magn += 2*sigma.spin[site];
           }//if
		   
		   //Otherwise reject
		   //else cout << " REJECT " << endl;
		}//else
	}//k

}//LocalUpdate


//Print Function
void IsingHamiltonian::print() {

	cout << "...printing nearest neighbors:" << endl << endl;
	
	for(int i=0; i<NearestNeighbors.size(); i++) {
		cout << "Site # ";
	    PRINT_RED(i);
		cout << "  -> ";
		
		for(int j=0; j<2*D; j++) {
		
			PRINT_GREEN(NearestNeighbors[i][j]);
		   	cout << " ";
		}//j
        
		cout << endl;
	}//i

	cout << endl;
	cout << "...printing bonds:" << endl << endl;

	for(int i=0; i<bonds.size(); i++) {
		cout << "Site # ";
	    PRINT_RED(i);
		cout << "  -> ";

		for(int j=0; j<bonds[i].size(); j++) {
			PRINT_GREEN(bonds[i][j]);
		    cout << "  ";
		}//j
		cout << endl;
	}//i
	
	cout << endl;
	cout << "...printing interactions:" << endl << endl;

	for(int i=0; i<J.size(); i++) {
		cout << "Bond # ";
	    PRINT_RED(i);
		cout << "  -> ";
		cout << "J = ";
	   	PRINT_GREEN(J[i]);
	    cout << endl;

	}//i
    
	cout << endl;

}//print
