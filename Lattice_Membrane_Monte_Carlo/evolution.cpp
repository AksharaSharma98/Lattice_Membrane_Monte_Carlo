#include <math.h>
#include <string>
#include <assert.h>

#include "lipid.h"
#include "membrane.h"
#include "evolution.h"
#include "parameters.h"
#include "energy.h"
#include "math_functions.h"
#include "output.h"
#include "qol_functions.h"
#include "Lattice_Membrane_Monte_Carlo.h"

// evolves two leaflets via metropolis monte carlo exchange moves.
void evolve_mc(membrane& upper, membrane& lower, int steps, int energy_output_freq, int config_output_freq) {

	// output initialization
	FILE* config_file; FILE* energy_file;
	config_file = fopen("config.txt", "w+");
	energy_file = fopen("energy.txt", "w+");
	//fprintf(config_file, "header");                         // insert function to write header

	double energy = system_energy(upper, lower);
	int c1[2] = {0,0}, c2[2] = {0,0};

	assert(steps >= 1 && "Invalid number of steps");
	for (int t = 0; t < steps; t++) {

		// pick non-degenerate lipids to exchange, attempt a monte-carlo exchange move
		lipid_picker(upper, c1, c2);
		energy += monte_carlo_move(upper, lower, c1, c2);

		// pick non-degenerate lipids to exchange, attempt a monte-carlo exchange move
		lipid_picker(lower, c1, c2);
		energy += monte_carlo_move(lower, upper, c1, c2);

		// output configuration at specified frequency
		if (t % energy_output_freq == 0) {
			write_energy(energy_file, energy);
		}
		if (t % config_output_freq == 0) {
			write_config_int(config_file, upper, lower);
		}
	}

	// clean-up
	fclose(config_file);
	fclose(energy_file);
}

// pick two non-degenerate grid points to swap
void lipid_picker(membrane leaflet, int c1[2], int c2[2]) {

	int n = leaflet.getgrid().size();

	bool picker = false;

	// pick first point
	for (int i = 0; i < 2; i++) {
		c1[i] = rand_int(0, n-1);
	}
	std::string l1 = leaflet.getlipid(c1[0], c1[1]).getspecies();

	// pick second point to be non-degenerate
	while (picker == false) {
		for (int i = 0; i < 2; i++) {
			c2[i] = rand_int(0, n-1);
		}
		std::string l2 = leaflet.getlipid(c2[0], c2[1]).getspecies();

		if (c1[0] != c2[0] && c1[1] != c2[1] && l1 != l2) {
			picker = true;
		}
	}
}


double monte_carlo_move(membrane& current, membrane& opposing, int* a, int* b) {

	// calculate initial energy
	double Ei = local_energy(current, opposing, a) + local_energy(current, opposing, b);

	// swap lipids specified for the move
	current.swap(a, b);

	// calculate final energy, energy difference
	double Ej = local_energy(current, opposing, a) + local_energy(current, opposing, b);
	double delE = Ej - Ei;

	// accept/reject attempted move
	bool accept = metropolis_accept(delE);

	// return system to initial state if move is rejected
	if (accept == false) {
		current.swap(a, b);
		return 0.0;
	}
	else {
		return delE;
	}
}


bool metropolis_accept(double delE) {
	
	double p = rand01();
	double accept = exp(-delE/kT);

	if (p<=accept) {
		return true;
	}
	else {
		return false;
	}
}