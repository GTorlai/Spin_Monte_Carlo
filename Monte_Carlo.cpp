
#include "hypercube.cpp"
#include "spins.cpp"

int main() {

    
	Hypercube cube(4,2);
	cube.print();
   
	Spins sigma(20);
	sigma.randomize();
	sigma.print();
   	
    
}
