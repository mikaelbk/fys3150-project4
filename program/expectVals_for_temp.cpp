
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
using namespace std;
using namespace arma;
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodicBoundary(int i, int limit, int add) { 
	return (i+limit+add) % (limit);
}

// function declaration
void metropolis(int, int, double, vec&, bool);
void alignSpin(int, mat&, double&, double&);
void randomizeSpin(int, mat&, double&, double&);
void WriteResultstoFile(int, int, double, vec);

int main(int argc, char const *argv[])
{
	// thrown message with bad usage of program call
	if (argc != 6) {
		cout << "\tError: Invalid number of arguments. Program requires the following call:" << endl;
		cout << "\t " << argv[0] << " L mcc ti tf nt" << endl;
		cout << "\tWhere the arguments are:\n\tL:\tSquare root of the number of lattice points for the LxL lattice\n\tmcc:\tNumber of Montecarlo cycles. If less than 20 then is set to 10^mcc\n\tti:\tstarting temperature of experiment\n\ttf:\tfinal temperature of experiment\n\tnt:\tnumber of measurements between ti and tf" << endl;
		exit(1);
	}

	// declaration of input arguments
	//string filename = argv[1];	// filename of the data .txt file
	int L = atoi(argv[1]);			// LxL size of lattice
	int mcc = atoi(argv[2]);		// number of MC cycles
	if(mcc < 20){mcc = int(pow(10,mcc));};
	double ti = atof(argv[3]);		// initial temperature
	double tf = atof(argv[4]);		// final temperature
	int nt = atoi(argv[5]);		// temperature step
	double dt;
	if(nt == 1){dt = 0;}else{dt = (tf-ti)/double(nt-1);}
	double energy = 0;	// variable for the energy
	double magMom = 0;	// variable for the magnetic momentum
	vec expectVals = zeros<mat>(5);	// list that contains various expectation values
	double T;
	//string filename = "test";
	string fileout = to_string(L) + "x" + to_string(L) + "_" + to_string(int(0.3+log10(mcc)));

	ofile.open(fileout);

	
	// loop that runs the 'experiment' and takes measurements at different temperatures
	for (int i = 0; i < nt; ++i)
	{
		T = ti + i*dt;
		//printf("T = %.8f\n", T);
		cout << "T = " << T << endl; 
		metropolis(L,mcc,T,expectVals,false);
		WriteResultstoFile(L, mcc, T, expectVals);
		expectVals = zeros<mat>(5);
	}
	ofile.close();

	return 0;
}

void metropolis(int L, int mcc, double T, vec& expectVals, bool rand){
	// RNG
	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

	// initializations
	double energy = 0.;
	double magMom = 0.;

	mat lattice = zeros<mat>(L,L);
	if(rand){randomizeSpin(L, lattice, energy, magMom);} else{alignSpin(L, lattice, energy, magMom);}

	// initialize array for expectation ExpectationValues
	// setup array for possible energy changes
	vec deVector = zeros<mat>(17); 
	for( int de =-8; de <= 8; de+=4) deVector(de+8) = exp(-de/T);
	// Start Monte Carlo cycles
	for (int cycle = 1; cycle <= mcc; cycle++){
		// The sweep over the lattice, looping over all spin sites
		for(int x =0; x < L; x++) {
			for (int y= 0; y < L; y++){
				int ix = (int) (RandomNumberGenerator(gen)*(double)L);
				int iy = (int) (RandomNumberGenerator(gen)*(double)L);
				int deltaE =  2*lattice(ix,iy)*
				(lattice(ix,periodicBoundary(iy,L,-1))+
				lattice(periodicBoundary(ix,L,-1),iy) +
				lattice(ix,periodicBoundary(iy,L,1)) +
				lattice(periodicBoundary(ix,L,1),iy));
				if ( RandomNumberGenerator(gen) <= deVector(deltaE+8) ) {
					lattice(ix,iy) *= -1.0;  // flip one spin and accept new spin config
					magMom += (double) 2*lattice(ix,iy);
					energy += (double) deltaE;
				}
			}
		}
		if(cycle > 0.5*mcc){
			expectVals(0) += energy;    expectVals(1) += energy*energy;
			expectVals(2) += magMom;    
			expectVals(3) += magMom*magMom; 
			expectVals(4) += fabs(magMom);
		}
	}
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

void WriteResultstoFile(int NSpins, int mcc, double temperature, vec ExpectationValues)
{
	double norm = 1.0/((double) (mcc/2));  // divided by  number of cycles
	double E_ExpectationValues = ExpectationValues(0)*norm;
	double E2_ExpectationValues = ExpectationValues(1)*norm;
	double M_ExpectationValues = ExpectationValues(2)*norm;
	double M2_ExpectationValues = ExpectationValues(3)*norm;
	double Mabs_ExpectationValues = ExpectationValues(4)*norm;

	// all expectation values are per spin, divide by 1/NSpins/NSpins
	double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
	double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setw(15) << setprecision(8) << temperature;
	ofile << setw(15) << setprecision(8) << E_ExpectationValues/NSpins/NSpins;
	ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
	ofile << setw(15) << setprecision(8) << M_ExpectationValues/NSpins/NSpins;
	ofile << setw(15) << setprecision(8) << Mvariance/temperature;
	ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins << endl;
} // end output function