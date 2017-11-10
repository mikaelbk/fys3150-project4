/*
Program to solve the two-dimensional Ising model
The coupling constant J = 1
Boltzmann's constant = 1, temperature has thus dimension energy
Metropolis sampling is used. Periodic boundary conditions.
*/
#include <iostream>
#include <armadillo>
#include <random>
using namespace std;
using namespace arma;

// inline function for periodic boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) { 
	return (i+limit+add) % (limit);
}

int main(int argc, char const *argv[])
{
	int L = 20; // lattice size is L²
	double T = 1.0;

	// Initialize the seed and call the Mersienne algo
	std::random_device rd;
	std::mt19937_64 gen(rd());
	// Set up the uniform distribution for x \in [[0, 1]
	std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
	// Initialize the lattice spin values
	mat SpinMatrix = zeros<mat>(L,L);

	return 0;
}

void alignSpin(int L, mat& lattice, double& energy, double& magMom){
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			lattice(i,j) = 1.0;
		}
	}
	magMom = L*L;

	for(int x =0; x < L; x++) {
		for (int y= 0; y < L; y++){
			Energy -=  (double) SpinMatrix(x,y)*
			(SpinMatrix(PeriodicBoundary(x,L,-1),y) +
			SpinMatrix(x,PeriodicBoundary(y,L,-1)));
		}
	}
}


void alignSpin(int L, mat& lattice, double& energy, double& magMom){
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			lattice(i,j) = 1.0;
			magMom += lattice(i,j);
		}
	}
}