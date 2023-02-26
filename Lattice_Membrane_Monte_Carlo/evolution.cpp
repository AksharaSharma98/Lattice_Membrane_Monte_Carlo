#include <math.h>
#include <string>
#include <assert.h>
#include <vector>

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
void evolve_mc(membrane& upper, membrane& lower, int steps, int tailorder_update_freq, int energy_output_freq, int config_output_freq) {

	// output initialization
	FILE* config_file; FILE* tailconfig_file; FILE* energy_file;
	config_file = fopen("config.txt", "w+");
	tailconfig_file = fopen("tailconfig.txt", "w+");
	energy_file = fopen("energy.txt", "w+");
	//fprintf(config_file, "header");                         // insert function to write header
	
	double energy = system_energy(upper, lower) / e; printf("System energy = %lf\n", energy);
	int c1[2] = {0,0}, c2[2] = {0,0};
	//std::vector<std::vector<int>> swap_sites_upper(2, std::vector<int>(0, 0));
	//std::vector<std::vector<int>> swap_sites_lower(2, std::vector<int>(0, 0));
	
	assert(steps >= 1 && "Invalid number of steps");
	printf("Start of system evolution\nCompletion: ");
	for (int t = 0; t < steps; t++) {

		// pick non-degenerate lipids to exchange, attempt a monte-carlo exchange move
		lipid_picker(upper, c1, c2);
		energy += monte_carlo_move(upper, lower, c1, c2)/e;// , swap_sites_upper);
		energy += update_tailorder_mc(upper, lower, c1)/e;
		energy += update_tailorder_mc(upper, lower, c2)/e;

		// pick non-degenerate lipids to exchange, attempt a monte-carlo exchange move
		lipid_picker(lower, c1, c2);
		energy += monte_carlo_move(lower, upper, c1, c2)/e;// , swap_sites_lower);
		energy += update_tailorder_mc(lower, upper, c1)/e;
		energy += update_tailorder_mc(lower, upper, c2)/e;

		// update system tail orders at specified frequency
		/*if (t % tailorder_update_freq == 0) {
			update_tailorder(upper, lower, swap_sites_upper);
			update_tailorder(lower, upper, swap_sites_lower);
		}*/

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


// pick two non-degenerate grid points to swap.
// use when lipids of the same species are different
void lipid_picker(membrane& leaflet, int c1[2], int c2[2]) {

	int n = leaflet.getgrid().size();

	bool picker = false;

	// pick first point
	for (int i = 0; i < 2; i++) {
		c1[i] = rand_int(0, n-1);
	}

	// pick second point to be non-degenerate
	while (picker == false) {
		for (int i = 0; i < 2; i++) {
			c2[i] = rand_int(0, n-1);
		}
		
		if (c1[0] != c2[0] && c1[1] != c2[1]) {
			picker = true;
		}
	}
}


// pick two non-degenerate grid points to swap.
// use when all lipids of the same species are identical
void lipid_picker_species(membrane& leaflet, int c1[2], int c2[2]) {

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

		if (l1 != l2) {
			picker = true;
		}
	}
}


double monte_carlo_move(membrane& current, membrane& opposing, int* a, int* b) {//, std::vector<std::vector<int>>& swap_sites) {
	
	// calculate initial energy
	double Ei = local_enthalpy(current, opposing, a) + local_enthalpy(current, opposing, b);
	Ei += local_planeentropy_env(current, a) + local_planeentropy_env(current, b);
	Ei += local_interentropy_env(current, opposing, a) + local_interentropy_env(current, opposing, b);

	// swap lipids specified for the move
	current.swap(a, b);

	// calculate final energy, energy difference
	double Ef = local_enthalpy(current, opposing, a) + local_enthalpy(current, opposing, b);
	Ef += local_planeentropy_env(current, a) + local_planeentropy_env(current, b);
	Ef += local_interentropy_env(current, opposing, a) + local_interentropy_env(current, opposing, b);
	double delE = Ef - Ei;

	// accept/reject attempted move
	bool accept = metropolis_accept(delE);

	// return system to initial state if move is rejected
	if (accept == false) {
		current.swap(a, b);
		return 0.0;
	}
	else {
		//swap_sites[0].push_back(a[0]);
		//swap_sites[1].push_back(a[1]);
		//swap_sites[0].push_back(b[0]);
		//swap_sites[1].push_back(b[1]);
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

void update_system_tailorder(membrane& current, membrane& opposing, std::vector<std::vector<int>>& swap_sites) {

	int n = current.getgrid().size();
	std::vector<std::vector<int>> seq(2, std::vector<int>(0, 0));
	update_sequence(n, swap_sites, seq);

}

void update_sequence(int n, std::vector<std::vector<int>> swap_sites, std::vector<std::vector<int>>& seq) {

}