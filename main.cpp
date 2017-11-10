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
inline int periodicBoundary(int i, int limit, int add) { 
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
	std::uniform_real_distribution<double> rand(0.0,1.0);
	
	mat lattice = zeros<mat>(L,L); // the lattice containing L² spins

	return 0;
}

void alignSpin(int L, mat& lattice, double& energy, double& magMom){
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			lattice(i,j) = 1.0;
		}
	}
	magMom += L*L;
	energy -= 2*L*L;
}


void randomizeSpin(int L, mat& lattice, double& energy, double& magMom){
	// cointoss function
	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_int_distribution<int> cointoss(0,1);

	// randomization loop
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			lattice(i,j) = cointoss(gen)*2-1;
			magMom += lattice(i,j);
		}
	}

	// set energy
	for(int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++){
			energy -=  (double) lattice(i,j) *
			(lattice(periodicBoundary(i,L,-1),j) +
			lattice(i,periodicBoundary(j,L,-1)));
		}
	}
}