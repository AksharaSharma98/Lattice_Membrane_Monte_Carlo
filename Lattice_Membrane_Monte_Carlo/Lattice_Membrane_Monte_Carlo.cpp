﻿// system header files

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
#include "system.h"
#include "initialize.h"
#include "evolution.h"
#include "Lattice_Membrane_Monte_Carlo.h"

#include "math_functions.h"

using namespace std;

// global variables
double e = 1.0;
double kT = 0.9;    // 0.802 is 280 K, 0.9 is 314 K

// initialize forcefield
Parameters forcefield;

// initialize system
System sys;


int main()
{
	// timekeeper
	clock_t clkStart = clock();

	// system initialization
	int steps, energy_out_freq, config_out_freq, restart_out_freq;
	std::pair<membrane, membrane> bilayer = initialize(steps, energy_out_freq, config_out_freq, restart_out_freq);
	printf("Completed system setup\n");

	// system evolution
	evolve_mc_farago(bilayer.first, bilayer.second, steps, energy_out_freq, config_out_freq, restart_out_freq);

	// testing (temporary)
	
	// wrap-up
	clock_t clkFinish = clock();
	printf("Runtime = %Lf seconds\n", (double)((clkFinish - clkStart)/CLOCKS_PER_SEC));

	return 0;
}
