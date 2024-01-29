#include <math.h>
#include <string>
#include <assert.h>
#include <vector>

#include "lipid.h"
#include "membrane.h"
#include "evolution.h"
#include "parameters.h"
#include "mc_moves.h"
#include "energy.h"
#include "math_functions.h"
#include "output.h"
#include "qol_functions.h"
#include "Lattice_Membrane_Monte_Carlo.h"


// evolves two leaflets via metropolis monte carlo exchange moves.
void evolve_mc_farago(membrane& upper, membrane& lower, int steps, int energy_output_freq, int config_output_freq) {

	// output initialization
	FILE* config_file; FILE* tailconfig_file; FILE* energy_file;
	config_file = fopen("config.txt", "w+");
	//tailconfig_file = fopen("tailconfig.txt", "w+");
	energy_file = fopen("energy.txt", "w+");
	//fprintf(config_file, "header");                         // insert function to write header
	
	double energy = (system_energy_farago(upper) + system_energy_farago(lower)) / e; 
	printf("System energy = %lf\n", energy);
	
	assert(steps >= 1 && "Invalid number of steps");
	printf("Start of system evolution\nCompletion: ");

	for (int t = 0; t < steps; t++) {
		printf("Energy before moves = %lf\n",energy);
		// pick non-degenerate lipids to exchange, attempt a chain of exchange moves
		energy += multi_swap(upper, 10) / e;
		energy += multi_swap(lower, 10) / e;
		
		// pick non-degenerate DPPC lipids, attempt a chain of state swap moves
		//energy += state_swap(upper, 10) / e;
		//energy += state_swap(lower, 10) / e;
		
		// pick non-degenerate lipids to exchange, attempt a patch-swap move
		//energy += patch_swap(upper) / e;
		//energy += patch_swap(lower) / e;

		// output configuration at specified frequency
		if (t % energy_output_freq == 0) {
			write_energy(energy_file, energy);
		}
		if (t % config_output_freq == 0) {
			write_config_int(config_file, upper, lower);
			//write_tailconfig(tailconfig_file, upper, lower);
		}
		if (t % (steps / 10) == 0) {
			printf("%d%% ", t / (steps / 100));
		}
	}
	printf("100%%\nSimulation complete\n");

	// clean-up
	fclose(config_file);
	//fclose(tailconfig_file);
	fclose(energy_file);
}




// Old functions

/*
// evolves two leaflets via metropolis monte carlo exchange moves.
void evolve_mc(membrane& upper, membrane& lower, int steps, int tailorder_update_freq, int energy_output_freq, int config_output_freq) {

	// output initialization
	FILE* config_file; FILE* tailconfig_file; FILE* energy_file;
	config_file = fopen("config.txt", "w+");
	tailconfig_file = fopen("tailconfig.txt", "w+");
	energy_file = fopen("energy.txt", "w+");
	//fprintf(config_file, "header");                         // insert function to write header
	
	double energy = system_energy(upper, lower) / e; printf("System energy = %lf\n", energy);
	int c1[2] = {0,0}, c2[2] = {0,0};
	
	assert(steps >= 1 && "Invalid number of steps");
	printf("Start of system evolution\nCompletion: ");
	for (int t = 0; t < steps; t++) {

		// pick non-degenerate lipids to exchange, attempt a monte-carlo exchange move
		lipid_picker(upper, c1, c2);
		energy += single_swap(upper, lower, c1, c2)/e;// , swap_sites_upper);
		energy += update_tailorder_mc(upper, lower, c1)/e;
		energy += update_tailorder_mc(upper, lower, c2)/e;

		// pick non-degenerate lipids to exchange, attempt a monte-carlo exchange move
		lipid_picker(lower, c1, c2);
		energy += single_swap(lower, upper, c1, c2)/e;// , swap_sites_lower);
		energy += update_tailorder_mc(lower, upper, c1)/e;
		energy += update_tailorder_mc(lower, upper, c2)/e;

		// output configuration at specified frequency
		if (t % energy_output_freq == 0) {
			write_energy(energy_file, energy);
		}
		if (t % config_output_freq == 0) {
			write_config_int(config_file, upper, lower);
			write_tailconfig(tailconfig_file, upper, lower);
		}
		if (t % (steps / 10) == 0) {
			printf("%d%% ", t / (steps / 100));
		}
	}
	printf("100%%\nSimulation complete\n");

	// clean-up
	fclose(config_file);
	fclose(tailconfig_file);
	fclose(energy_file);
}


double update_tailorder_mc(membrane& current, membrane& opposing, int* a) {

	lipid l = current.getlipid(a[0], a[1]);

	// calculate initial energy
	double si = l.gettail_order();
	double Ei = local_enthalpy(current, opposing, a) + local_planeentropy_env(current, a) + local_interentropy_env(current, opposing, a);

	// sample a new tail order value
	double sf = sample_tailorder(l.getspecies(), si);
	current.tail_update(a, sf);

	// calculate final energy, energy difference
	double Ef = local_enthalpy(current, opposing, a) + local_planeentropy_env(current, a) + local_interentropy_env(current, opposing, a);
	double delE = Ef - Ei;

	// accept/reject attempted tail order update
	bool accept = metropolis_accept(delE);

	// return tail order to initial state if move is rejected
	if (accept == false) {
		current.tail_update(a, si);
		return 0.0;
	}
	else {
		return delE;
	}

	return 0.0;
}
*/