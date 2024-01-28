#include <math.h>
#include <string>
#include <assert.h>
#include <vector>
#include <set>

#include "lipid.h"
#include "membrane.h"
#include "parameters.h"
#include "mc_moves.h"
#include "energy.h"
#include "periodicboundary.h"
#include "math_functions.h"
#include "qol_functions.h"
#include "Lattice_Membrane_Monte_Carlo.h"


double single_swap(membrane& current) {

	// pick non-degenerate lipids to exchange
	std::vector<int> a{ 0, 0 }, b{ 0, 0 };
	//lipid_picker_species(current, a, b);              ///// NEED TO MODIFY BEFORE USING

	// calculate initial energy
	double Ei = local_energy_farago(current, a) + local_energy_farago(current, b);

	// swap lipids specified for the move
	current.swap(a, b);

	// calculate final energy, energy difference
	double Ef = local_energy_farago(current, a) + local_energy_farago(current, b);
	double delE = Ef - Ei;

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


double multi_swap(membrane& current, int batch_size) {

	// pick and store a sequence of pairs of non-degenerate lipids to exchange
	std::vector<int> a{ 0, 0 }, b{ 0, 0 };
	std::vector<std::vector<int> > sites1, sites2;

	for (int i = 0; i < batch_size; i++) {
		lipid_picker_species(current, sites1, sites2);
	}
	
	// find and store unique interaction pairs involved in the swap chain
	std::map<std::set<int>, std::vector<int>> s;
	unique_interaction_pairs(current, sites1, sites2, s);
	
	// calculate initial energy by looping over unique sites
	double Ei = 0.0;
	for (const auto& [key, value] : s) {
		std::vector<int> pair_coords = value;
		Ei += pair_energy_farago(current, pair_coords);
	}
	
	// loop over stored swap sites and execute swaps
	for (int i = 0; i < batch_size; i++) {
		current.swap(sites1[i], sites2[i]);
	}
	
	// calculate final energy by looping over unique sites, energy difference
	double Ef = 0.0;
	for (const auto& [key, value] : s) {
		std::vector<int> pair_coords = value;
		Ef += pair_energy_farago(current, pair_coords);
	}
	double delE = Ef - Ei;
	
	// accept/reject attempted move
	bool accept = metropolis_accept(delE);

	// return system to initial state if move is rejected
	if (accept == false) {
		// loop over stored swap sites in reverse and execute swaps
		for (int i = batch_size-1; i >= 0; i--) {
			current.swap(sites1[i], sites2[i]);
		}
		return 0.0;
	}
	else {
		return delE;
	}
}


double state_swap(membrane& current, int batch_size) {

	// pick and store a sequence of pairs of non-degenerate lipids to exchange
	std::vector<int> a{ 0, 0 };
	std::vector<std::vector<int> > sites;
	
	for (int i = 0; i < batch_size; i++) {
		lipid_picker_DPPC(current, sites);
	}

	// find and store unique interaction pairs involved in the swap chain
	std::map<std::set<int>, std::vector<int> > s;
	unique_interaction_pairs(current, sites, s);

	// calculate initial energy by looping over unique sites
	double Ei = 0.0;
	for (const auto& [key, value] : s) {
		std::vector<int> pair_coords = value;
		Ei += pair_energy_farago(current, pair_coords);
	}
	for (int i = 0; i < batch_size; i++) {
		Ei += entropy_farago(current, sites[i]);
	}
	// loop over stored sites and execute state-swaps
	for (int i = 0; i < batch_size; i++) {
		current.swap_DPPC_state(sites[i]);
	}

	// calculate final energy by looping over unique sites, energy difference
	double Ef = 0.0;
	for (const auto& [key, value] : s) {
		std::vector<int> pair_coords = value;
		Ef += pair_energy_farago(current, pair_coords);
	}
	for (int i = 0; i < batch_size; i++) {
		Ef += entropy_farago(current, sites[i]);
	}
	double delE = Ef - Ei;

	// accept/reject attempted move
	bool accept = metropolis_accept(delE);
	
	// return system to initial state if move is rejected
	if (accept == false) {
		// loop over stored swap sites in reverse and execute swaps
		for (int i = batch_size - 1; i >= 0; i--) {
			current.swap_DPPC_state(sites[i]);
		}
		return 0.0;
	}
	else {
		return delE;
	}
}


double patch_swap(membrane& current) {

	// sample a patch size from input distribution
	int size = sample_swap_size();
	assert(size > 1 && "Invalid patch size used. Patches must be larger than 1 lipid\n");
	
	// pick two valid patch centers
	std::vector<int> a{ 0, 0 }, b{ 0, 0 };
	patch_center_picker(current, a, b, size);
	
	// find and store patch boundaries
	std::vector<std::vector<int> > boundary1, boundary2;
	std::vector<int> bounds1, bounds2; // contains {x1, x2, y1, y2} bounds of the respective patches
	patch_boundaries(current, a, size, boundary1, bounds1);
	patch_boundaries(current, b, size, boundary2, bounds2);

	// find and store unique interaction pairs involved in the swap chain
	std::map<std::set<int>, std::vector<int> > s;
	unique_interaction_pairs(current, boundary1, boundary2, s);

	// calculate initial energy along patch boundaries
	double Ei = 0.0;
	for (const auto& [key, value] : s) {
		std::vector<int> pair_coords = value;
		Ei += pair_energy_farago(current, pair_coords);
	}
	
	// swap patches
	current.patch_swap(bounds1, bounds2, size);

	// calculate final energy, energy difference
	double Ef = 0.0;
	for (const auto& [key, value] : s) {
		std::vector<int> pair_coords = value;
		Ef += pair_energy_farago(current, pair_coords);
	}
	double delE = Ef - Ei;
	
	// accept/reject attempted move
	bool accept = metropolis_accept(delE);

	// return system to initial state if move is rejected
	if (accept == false) {
		current.patch_swap(bounds1, bounds2, size);
		return 0.0;
	}
	else {
		return delE;
	}
}


// pick two non-degenerate grid points to swap.
// use when lipids of the same species are different
void lipid_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2) {

	int n = leaflet.getgrid().size();

	bool picker = false;

	// pick first point
	for (int i = 0; i < 2; i++) {
		c1[i] = rand_int(0, n - 1);
	}

	// pick second point to be non-degenerate
	while (picker == false) {
		for (int i = 0; i < 2; i++) {
			c2[i] = rand_int(0, n - 1);
		}

		if (c1[0] != c2[0] && c1[1] != c2[1]) {
			picker = true;
		}
	}
}


// pick two non-degenerate grid points to swap.
// use when all lipids of the same species are identical
void lipid_picker_species(membrane& leaflet, std::vector<std::vector<int> > &sites1, std::vector<std::vector<int> > &sites2) {

	int n = leaflet.getgrid().size();

	bool picker = false;
	std::vector<int> c1 = { 0,0 };
	std::vector<int> c2 = { 0,0 };

	// pick second point to be non-degenerate
	while (picker == false) {
		// pick 2 points
		for (int i = 0; i < 2; i++) {
			c1[i] = rand_int(0, n - 1);
			c2[i] = rand_int(0, n - 1);
		}
		// ensure they're not redundant
		std::string l1 = leaflet.getlipid(c1[0], c1[1]).getspecies();
		std::string l2 = leaflet.getlipid(c2[0], c2[1]).getspecies();
		if (l1 != l2) {
			picker = true;
		}

		// ensure they're a unique pair of points
		for (int i = 0; i < sites1.size(); i++) {
			if (c1[0] == sites1[i][0] && c1[1] == sites1[i][1] && c2[0] == sites2[i][0] && c2[1] == sites2[i][1]) {
				picker = false;
				break;
			}
			if (c1[0] == sites2[i][0] && c1[1] == sites2[i][1] && c2[0] == sites1[i][0] && c2[1] == sites1[i][1]) {
				picker = false;
				break;
			}
		}

	}
	if (picker == true) {
		sites1.push_back(c1);
		sites2.push_back(c2);
	}
}


void lipid_picker_DPPC(membrane& leaflet, std::vector<std::vector<int> >& sites) {

	int n = leaflet.getgrid().size();

	bool picker = false;
	std::vector<int> c = { 0,0 };

	// pick a grid point that is a DPPC lipid in either of its states
	while (picker == false) {
		for (int i = 0; i < 2; i++) {
			c[i] = rand_int(0, n - 1);
		}

		std::string l = leaflet.getlipid(c[0], c[1]).getspecies();
		if (l == "DPPCd" || l == "DPPCo") {
			picker = true;
		}

		for (int i = 0; i < sites.size(); i++) {
			if (c[0] == sites[i][0] && c[1] == sites[i][1]) {
				picker = false;
				break;
			}
		}
	}
	if (picker == true) {
		sites.push_back(c);
	}
}


bool metropolis_accept(double delE) {

	double p = rand01();
	double accept = exp(-delE / kT);

	if (p <= accept) {
		return true;
	}
	else {
		return false;
	}
}


void unique_interaction_pairs(membrane& current, const std::vector<std::vector<int> > &sites1, const std::vector<std::vector<int> > &sites2, std::map<std::set<int>, std::vector<int> > &unique) {

	int x, y;
	int n = current.getgrid().size();
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	// append sites into one list
	std::vector<std::vector<int> > sites = sites1;
	sites.insert(sites.end(), sites2.begin(), sites2.end());

	for (int i = 0; i < sites.size(); i++) {
		
		x = n * sites[i][0] + sites[i][1];

		// get list of neighbours with PBC correction for ith site
		periodic_neighbours(current, sites[i][0], sites[i][1], nbs);

		for (int j = 0; j < 6; j++) {
			y = n * nbs[j][0] + nbs[j][1];
			std::set<int> key = {x, y};
			unique.insert({ key , {sites[i][0], sites[i][1], nbs[j][0], nbs[j][1]} });
		}
	}
}


void unique_interaction_pairs(membrane& current, const std::vector<std::vector<int> >& sites, std::map<std::set<int>, std::vector<int> >& unique) {

	int x, y;
	int n = current.getgrid().size();
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	for (int i = 0; i < sites.size(); i++) {

		x = n * sites[i][0] + sites[i][1];

		// get list of neighbours with PBC correction for ith site
		periodic_neighbours(current, sites[i][0], sites[i][1], nbs);

		for (int j = 0; j < 6; j++) {
			y = n * nbs[j][0] + nbs[j][1];
			std::set<int> key = { x, y };
			unique.insert({ key , {sites[i][0], sites[i][1], nbs[j][0], nbs[j][1]} });
		}
	}
}


void patch_boundaries(membrane& current, const std::vector<int> &a, int patch_size, std::vector<std::vector<int> >& boundary, std::vector<int>& bounds) {
	
	int n = current.getgrid().size();
	
	int shift = patch_size / 2;
	//printf("shift = %d\n", shift);
	if (patch_size % 2 == 0) {
		int x = (a[0] + 1 - shift) % n, y = (a[1] + 1 - shift) % n;
		if (x < 0) { x = n + x; }
		if (y < 0) { y = n + y; }
		bounds.push_back(x);
		bounds.push_back((a[0] + shift) % n);
		bounds.push_back(y);
		bounds.push_back((a[1] + shift) % n);
	}
	else {
		int x = (a[0] - shift) % n, y = (a[1] - shift) % n;
		if (x < 0) { x = n + x; }
		if (y < 0) { y = n + y; }
		bounds.push_back(x);
		bounds.push_back((a[0] + shift) % n);
		bounds.push_back(y);
		bounds.push_back((a[1] + shift) % n);
	}

	for (int i = 0; i < patch_size-1; i++) {
		int x = (bounds[1] - i) % n, y = (bounds[3] - i) % n;
		if (x < 0) { x = n + x; }
		if (y < 0) { y = n + y; }
		boundary.push_back(std::vector<int>{ bounds[0], (bounds[2] + i) % n });
		boundary.push_back(std::vector<int>{ (bounds[0] + i) % n, bounds[3] });
		boundary.push_back(std::vector<int>{ bounds[1], y });
		boundary.push_back(std::vector<int>{ x, bounds[2] });
	}
}


void patch_center_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int patch_size) {

	int n = leaflet.getgrid().size();

	// pick first point
	for (int i = 0; i < 2; i++) {
		c1[i] = rand_int(0, n - 1);
	}

	// calculate invalid area to avoid picking overlapping patches
	int invalid_bounds[4] = {0, 0, 0, 0};
	if (patch_size%2 == 0){
		int x = (c1[0] + 1 - patch_size) % n, y = (c1[1] + 1 - patch_size) % n;
		if (x < 0) { x = n + x; }
		if (y < 0) { y = n + y; }
		invalid_bounds[0] = x;
		invalid_bounds[1] = (c1[0] + patch_size) % n;
		invalid_bounds[2] = y;
		invalid_bounds[3] = (c1[1] + patch_size) % n;
	}
	else {
		int x = (c1[0] - patch_size) % n, y = (c1[1] - patch_size) % n;
		if (x < 0) { x = n + x; }
		if (y < 0) { y = n + y; }
		invalid_bounds[0] = x;
		invalid_bounds[1] = (c1[0] + patch_size) % n;
		invalid_bounds[2] = y;
		invalid_bounds[3] = (c1[1] + patch_size) % n;
	}

	// pick second point such that patches don't overlap
	bool invalid_x = true, invalid_y = true;
	while (invalid_x == true && invalid_y == true) {
		for (int i = 0; i < 2; i++) {
			c2[i] = rand_int(0, n - 1);
		}

		if (invalid_bounds[0] < invalid_bounds[1]) {
			if (c2[0] < invalid_bounds[0] || c2[0] > invalid_bounds[1]) {
				invalid_x = false;
			}
		}
		else {
			if (c2[0] < invalid_bounds[0] && c2[0] > invalid_bounds[1]) {
				invalid_x = false;
			}
		}

		if (invalid_bounds[2] < invalid_bounds[3]) {
			if (c2[1] < invalid_bounds[2] || c2[1] > invalid_bounds[3]) {
				invalid_y = false;
			}
		}
		else {
			if (c2[1] < invalid_bounds[2] && c2[1] > invalid_bounds[3]) {
				invalid_y = false;
			}
		}
	}
}





// Old functions

/*double single_swap(membrane& current, membrane& opposing, int* a, int* b) {

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
		return delE;
	}
}*/