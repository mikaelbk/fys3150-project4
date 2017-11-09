/*
Program to solve the two-dimensional Ising model
The coupling constant J = 1
Boltzmann's constant = 1, temperature has thus dimension energy
Metropolis sampling is used. Periodic boundary conditions.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

// inline function for periodic boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) { 
	return (i+limit+add) % (limit);
}

int main(int argc, char const *argv[])
{
	cout << "hello world\n";
	return 0;
}

