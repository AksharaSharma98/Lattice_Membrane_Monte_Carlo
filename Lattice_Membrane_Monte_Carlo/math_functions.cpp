#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <vector>
#include <string>

#include "membrane.h"
#include "math_functions.h"
#include "Lattice_Membrane_Monte_Carlo.h"

// set RNG seed
std::mt19937 mt(1790);
std::uniform_real_distribution<double> uniform(0.0, 1.0);


// returns a random number between (including) 0 and 1
double rand01() {

	return uniform(mt);
}


// returns a random integer between (including) lower and 
// upper bounds.
// Example: If you want to sample the numbers 0, 1, and 2
// randomly, call rand_int(0,2)
int rand_int(int lower_bound, int upper_bound) {

	std::uniform_int_distribution<int> uniform_int(lower_bound, upper_bound);

	return uniform_int(mt);
}


// samples the swap patch size distribution
int sample_swap_size() {
	std::vector<int> sizes = sys.get_swap_sizes();
	std::vector<double> weights = sys.get_swap_weights();

	std::piecewise_constant_distribution<double> dist(sizes.begin(), sizes.end(), weights.begin());

	return dist(mt);
}


// samples the average S_CD distribution for a lipid species
double sample_tailorder(std::string species, double s_old) {

	double s = 0.0;
	std::vector<double> bins = forcefield.tailorder_bins(species);
	std::vector<double> weights = forcefield.tailorder_weights(species);

	std::piecewise_constant_distribution<double> dist(bins.begin(), bins.end(), weights.begin());
	bool picker = false;
	while (picker == false) {
		s = dist(mt);
		if (s != s_old) {
			picker = true;
		}
	}

	return s;
}


void mole_fraction_plane(membrane& current, int* x, std::vector<double>& phi, int nbs[][2]) {

	std::string lc;
	int n_sp = sys.get_num_species();

	std::vector<double> count(n_sp, 0.0);
	for (int i = -1; i < 6; i++) {
		if (i == -1) {
			lc = current.getlipid(x[0], x[1]).getspecies();
		}
		else {
			lc = current.getlipid(nbs[i][0], nbs[i][1]).getspecies();
		}

		for (int j = 0; j < n_sp; j++) {
			if (lc == sys.get_species(0,j)) {
				count[j] += 1.0;
			}
		}

		for (int j = 0; j < n_sp; j++) {
			phi[j] = count[j] / 7.0;
		}
	}
}


void mole_fraction_total(membrane& current, membrane& opposing, int* x, std::vector<double>& phi, int nbs[][2]) {

	std::string lc, lo;
	int n_sp = sys.get_num_species();

	std::vector<double> count(n_sp, 0.0);
	for (int i = -1; i < 6; i++) {
		if (i == -1) {
			lc = current.getlipid(x[0], x[1]).getspecies();
			lo = opposing.getlipid(x[0], x[1]).getspecies();
		}
		else {
			lc = current.getlipid(nbs[i][0], nbs[i][1]).getspecies();
			lo = opposing.getlipid(nbs[i][0], nbs[i][1]).getspecies();
		}
		
		for (int j = 0; j < n_sp; j++) {
			if (lc == sys.get_species(0, j)) {
				count[j] += 1.0;
			}
			if (lo == sys.get_species(0, j)) {
				count[j] += 1.0;
			}
		}

		for (int j = 0; j < n_sp; j++) {
			phi[j] = count[j] / 14.0;
		}
	}
}