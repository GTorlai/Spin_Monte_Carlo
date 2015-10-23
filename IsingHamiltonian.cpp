#include "IsingHamiltonian.h"
#include <vector>

//Constructor
IsingHamiltonian::IsingHamiltonian(Spins & sigma, Hypercube & cube, MTRand & random) {

    N = cube.N;
	D = cube.D;

	sigma.resize(N);
	sigma.randomize(random);

	J.assign(N*D,1.0);
	
	//Initialize the Spin Spin Correlation matrix
	vector<int> temp;
	temp.assign(N,0);
	for(int i=0; i<N; i++) {
		SpinSpinCorr.push_back(temp);
	}//i

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

//Initialize the interactions coupling at random between 1 and -1 according
//to the probability p
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

}//GetEnergy

//Computer the magnetization in the system
void IsingHamiltonian::GetMagnetization(Spins & sigma) {
    
    Magn=0;
    for(int i=0; i<sigma.N; i++) {
        Magn += sigma.spin[i];
    }//i
    
}//GetMagnetization

//Update Energy, Magnetization and Correlation after a spin flip
void IsingHamiltonian::Update(Spins & sigma, double dE, int site) {
	
	//Update Energy
	Energy += dE;

	//Update Total Magnetization
	Magn += 2*sigma.spin[site];

}
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
			Update(sigma,deltaE,site);
			//Energy += deltaE;
            //Magn += 2*sigma.spin[site];
		}//if
		
		else {
			ran_num = random.rand();
			//cout << "Metropolis weight: " << exp(-deltaE/T) << "\t Random:" << ran_num << endl;
		   if (exp(-deltaE/T) > ran_num) {
		   		sigma.flip(site);
				Update(sigma,deltaE,site);
				//Energy += deltaE;
                //Magn += 2*sigma.spin[site];
           }//if
		   
		   //Otherwise reject
		   //else cout << " REJECT " << endl;
		}//else
	}//k

}//LocalUpdate

void IsingHamiltonian::GrowCluster(Spins & sigma, double & T, MTRand & random, int & site) {
    
    int newSite;
    int State = sigma.spin[site];
    double deltaE;
    
    deltaE = 0.0;
    for(int i=0; i<NearestNeighbors[site].size(); i++) {
        deltaE += 2*sigma.spin[site]*sigma.spin[NearestNeighbors[site][i]]*J[bonds[site][i]];
    }//i
    
    sigma.flip(site);
    Update(sigma,deltaE,site);
	//Energy += deltaE;
    //Magn += 2*sigma.spin[site];
    
    double ran_num;
    
    for(int i=0; i<NearestNeighbors[site].size(); i++) {
        newSite = NearestNeighbors[site][i];
        
        if (sigma.spin[newSite] == State) {
            ran_num=random.rand();
            //if(m_rand < 1-exp(2*J[OnesConnectedToZero[site][i]]/T)) {
            if(ran_num < 1-exp(2/T)) {
                GrowCluster(sigma,T,random,newSite);
            }
        }
    }
    
}
//Calculate the Wolff cluster update
void IsingHamiltonian::GaugeUpdate(Spins & sigma, double & T, MTRand & random){
    
    int site; //random site where the cluster starts
    
    site = random.randInt(N-1);
    
    GrowCluster(sigma,T,random,site);
    
}

//Print Function
void IsingHamiltonian::print() {

/*	cout << "...printing nearest neighbors:" << endl << endl;
	
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
*/
	for(int i=0; i<J.size(); i++) {
		//cout << "Bond # ";
	    //PRINT_RED(i);
		//cout << "  -> ";
		//cout << "J = ";
	   	PRINT_GREEN(J[i]);
		cout << " ";
	    //cout << endl;

	}//i
    
	cout << endl;

}//print
