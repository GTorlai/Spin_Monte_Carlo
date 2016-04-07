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
    
    //Annealed Importance Sampling Parameters
    
    K = 1000;
    M = 10;

    for (int i=0; i< K/10; i++) {
        beta.push_back(0.5*i/(0.1*K));
    }
    for (int i=0; i< 4*K/10; i++) {
        beta.push_back(0.5+0.4*i/(0.4*K));
    }
     for (int i=0; i< K/2; i++) {
        beta.push_back(0.9+0.1*i/(0.5*K));
    }
    beta.push_back(1.0);

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

//Performed simulated annealing to desired temperature
void IsingHamiltonian::SimulatedAnnealing(Spins & sigma, double T_f, MTRand & random) {
	
    double Temp = 10.0;	
    int TSteps = 50;
    double deltaT = 1.0*(Temp-T_f)/TSteps;
    int annealSteps = 4000;

    for(int k=0; k<TSteps; k++) {
		
	for(int i=0; i<annealSteps; i++) {
	    LocalUpdate(sigma,Temp,random);
	}//i
		
	Temp -= deltaT;
	//cout << "Current temperature: " << Temp << " towards " << T_f << endl;
    }//k

}//SimulatedAnnealing

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
    //for(int i=0; i<NearestNeighbors[site].size(); i++) {
    //    deltaE += 2*sigma.spin[site]*sigma.spin[NearestNeighbors[site][i]]*J[bonds[site][i]];
    //}//i
    
    sigma.flip(site);
    //Update(sigma,deltaE,site);
    //Energy += deltaE;
    //Magn += 2*sigma.spin[site];
    
    double ran_num;
    
    for(int i=0; i<NearestNeighbors[site].size(); i++) {
        newSite = NearestNeighbors[site][i];
        
        if (sigma.spin[newSite] == State) {
            ran_num=random.rand();
            //if(m_rand < 1-exp(2*J[OnesConnectedToZero[site][i]]/T)) {
            if(ran_num < 1-exp(-2.0/T)) {
                GrowCluster(sigma,T,random,newSite);
            }
        }
    }
    
}
//Calculate the Wolff cluster update
void IsingHamiltonian::GaugeUpdate(Spins & sigma, double & T, MTRand & random){
    
    int site; //random site where the cluster starts
    
    for (int i=0; i<N; i++) {
        site = random.randInt(N-1);
        GrowCluster(sigma,T,random,site);
    }
}


//Sweep of AIS
double IsingHamiltonian::aisSweep(Spins & sigma,MTRand & random,double T,double T0) {

    double W = 0.0;
    double p_new;
    double p_old;
    double T_old;
    double T_new;

    for (int k=0; k<K; k++) {

        T_old = (T0*T)/((1-beta[k-1])*T + beta[k-1]*T0);
        T_new = (T0*T)/((1-beta[k])*T + beta[k]*T0);

        LocalUpdate(sigma, T_new,random);    
        GetEnergy(sigma);

        p_new = exp(-1.0*Energy/T_new);
        p_old = exp(-1.0*Energy/T_old);

        W += log(p_new) - log(p_old);
    }

    return W;

}

void IsingHamiltonian::AnnealedImportanceSampling(Spins & sigma, MTRand & random, double T, double T0){

    
    double w = 0.0;
    int eq = 10000;
    double sweep;

    for (int i=0;i<M;i++) {

        sigma.randomize(random);
        
        for (int i=0; i<eq;i++) {
            LocalUpdate(sigma,T0,random);
        }
        
        sweep = aisSweep(sigma,random,T,T0);
        w += sweep;
   }
    
    w /= M;

    cout << "Ratio: " << w << endl << endl;



}

//Print Function
void IsingHamiltonian::print() {

    //cout << "...printing nearest neighbors:" << endl << endl;
    
    for(int i=0; i<NearestNeighbors.size(); i++) {
    	cout << "Site # ";
        PRINT_RED(i);
    	cout << "  -> ";
        //cout << i << " ";
    	for(int j=0; j<2*D; j++) {
    	
    	    PRINT_GREEN(NearestNeighbors[i][j]);
            //cout << NearestNeighbors[i][j];
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
