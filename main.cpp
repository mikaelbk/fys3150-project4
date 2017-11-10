
#include <iostream>
#include <armadillo>
#include <random>
#include <math.h>
using namespace std;
using namespace arma;

// inline function for periodic boundary conditions
inline int periodicBoundary(int i, int limit, int add) { 
	return (i+limit+add) % (limit);
}

int main(int argc, char const *argv[])
{
	// thrown message with bad usage of program call
	if (argc != 7) {
		cout << "\tError: Invalid number of arguments. Program requires the following call:" << endl;
		cout << "\t " << argv[0] << " outFileName L mcc ti tf dt" << endl;
		cout << "\tWhere the arguments are:\n\toutFileName:\tName that will be given to the output file\n\tL:\tSquare root of the number of lattice points for the LxL lattice\n\tmcc:\tNumber of Montecarlo cycles. If less than 20 then is set to 10^mcc\n\tti:\tstarting temperature of experiment\n\ttf:\tfinal temperature of experiment\n\tdt:\ttemperature step between measurements" << endl;
		exit(1);
	}

	// declaration of input arguments
	string outFileName = argv[1];
	int L = atoi(argv[2]);
	int mcc = atoi(argv[3]);
	if(mcc < 20){mcc = int(pow(10,mcc));};
	double ti = atof(argv[4]);
	double tf = atof(argv[5]);
	double dt = atof(argv[6]);

	cout << outFileName << endl;
	cout << L << endl;
	cout << mcc << endl;
	cout << ti << endl;
	cout << tf << endl;
	cout << dt << endl;

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