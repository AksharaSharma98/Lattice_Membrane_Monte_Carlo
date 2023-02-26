// system header files

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>

// program header files

//#include "lipid.h"
#include "membrane.h"
#include "parameters.h"
#include "evolution.h"
#include "Lattice_Membrane_Monte_Carlo.h"

#include "math_functions.h"

using namespace std;

// global variables
double kT = 1.0;
int n_sp = 3;
int n = 40;
std::string species[2][3] = { {"DPPC", "DIPC", "CHOL"},
							  {"DPPC", "DIPC", "CHOL"} };
int pop[2][3] = { {640, 480, 480},
				  {640, 480, 480} };

// initialize forcefield
parameters forcefield;


int main()
{
	// timekeeper
	clock_t clkStart = clock();

	// system initialization
	membrane upper(n, species[0], pop[0], n_sp);
	membrane lower(n, species[1], pop[1], n_sp); printf("Completed system setup\n");
	// system evolution
	printf("Starting MC evolution\n");
	evolve_mc(upper, lower, 100000, 1, 10, 10000);

	// testing (temporary)
	
	// wrap-up

	clock_t clkFinish = clock();
	printf("Runtime = %Lf seconds\n", (double)((clkFinish - clkStart)/CLOCKS_PER_SEC));

	return 0;
}
