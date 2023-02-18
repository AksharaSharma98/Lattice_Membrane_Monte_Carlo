// system header files

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <map>

// program header files

#include "lipid.h"
#include "membrane.h"
#include "parameters.h"
#include "evolution.h"
#include "math_functions.h"
#include "output.h"
#include "qol_functions.h"
#include "Lattice_Membrane_Monte_Carlo.h"
#include "periodicboundary.h"


using namespace std;

// global variables
double kT = 1.0;

// initialize forcefield
parameters forcefield;


int main()
{
	// timekeeper
	clock_t clkStart = clock();

	// system variable definitions
	int n = 8;
	std::string species[2][2] = { {"DPPC", "DIPC"}, 
								  {"DPPC", "DIPC"} };
	int pop[2][2] = { {32, 32}, 
					  {32, 32} };
	double tail[2][2] = { { 0.5, 0.2 }, 
						  { 0.5, 0.2 } };

	// system initialization
	membrane upper(n, species[0], pop[0], tail[0], 2);
	membrane lower(n, species[1], pop[1], tail[1], 2);
	
	// system evolution
	evolve_mc(upper, lower, 1000000, 1, 1000);

	// testing

	// wrap-up

	clock_t clkFinish = clock();
	printf("Runtime = %Lf seconds\n", (long double)((clkFinish - clkStart)/CLOCKS_PER_SEC));

	return 0;
}
