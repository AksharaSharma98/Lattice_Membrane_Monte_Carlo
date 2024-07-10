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


double multi_swap(membrane& current, int batch_size, std::vector<double>& acceptance) {

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
	for (std::map<std::set<int>, std::vector<int> >::iterator it = s.begin(); it != s.end(); it++) {
		Ei += pair_energy_farago(current, it->second);
	}

	// loop over stored swap sites and execute swaps
	for (int i = 0; i < batch_size; i++) {
		current.swap(sites1[i], sites2[i]);
	}

	// calculate final energy by looping over unique sites, energy difference
	double Ef = 0.0;
	for (std::map<std::set<int>, std::vector<int> >::iterator it = s.begin(); it != s.end(); it++) {
		Ef += pair_energy_farago(current, it->second);
	}
	double delE = Ef - Ei;

	// accept/reject attempted move
	bool accept = metropolis_accept(delE);

	// return system to initial state if move is rejected
	if (accept == false) {
		// loop over stored swap sites in reverse and execute swaps
		for (int i = batch_size - 1; i >= 0; i--) {
			current.swap(sites1[i], sites2[i]);
		}
		// log acceptance ratio
		acceptance[1] += 1.0;

		return 0.0;
	}
	else {
		// log acceptance ratio
		acceptance[0] += 1.0; acceptance[1] += 1.0;

		return delE;
	}
}


double state_swap(membrane& current, int batch_size, std::vector<double>& acceptance) {

	if (sys.get_population(current.getleafletindex(), 1) == 0 && sys.get_population(current.getleafletindex(), 2) == 0) {
		return 0.0;
	}

	// pick and store a sequence of pairs of non-degenerate lipids to exchange
	std::vector<int> a{ 0, 0 };
	std::vector<std::vector<int> > sites;

	for (int i = 0; i < batch_size; i++) {
		site_picker_DPPC(current, sites);
	}

	// find and store unique interaction pairs involved in the swap chain
	std::map<std::set<int>, std::vector<int> > s;
	unique_interaction_pairs(current, sites, s);

	// calculate initial energy by looping over unique sites
	double Ei = 0.0;
	for (std::map<std::set<int>, std::vector<int> >::iterator it = s.begin(); it != s.end(); it++) {
		Ei += pair_energy_farago(current, it->second);
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
	for (std::map<std::set<int>, std::vector<int> >::iterator it = s.begin(); it != s.end(); it++) {
		Ef += pair_energy_farago(current, it->second);
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
		// log acceptance ratio
		acceptance[1] += 1.0;

		return 0.0;
	}
	else {
		// log acceptance ratio
		acceptance[0] += 1.0; acceptance[1] += 1.0;

		return delE;
	}
}


double patch_swap(membrane& current, std::map<int, std::vector<double> >& patch_accept) {

	// pick two valid patches
	std::vector<int> a{ 0, 0 }, b{ 0, 0 };
	int size;
	patch_picker(current, a, b, size);

	// find and store patch sites and square boundaries
	std::vector<std::vector<int> > patch1, patch2;
	std::vector<int> bounds1, bounds2; // contains {x1, x2, y1, y2} bounds of the respective patches
	patch_sites(current, a, size, patch1, bounds1);
	patch_sites(current, b, size, patch2, bounds2);

	// find and store optimal indent-matching sections for patches
	std::vector<std::vector<int> > indent_sections;
	indentation_sections(current, a, b, indent_sections);

	// find and store unique interaction pairs involved in the swap chain
	std::map<std::set<int>, std::vector<int> > s;
	unique_interaction_pairs(current, patch1, patch2, s);

	// calculate initial energy along patch boundaries
	double Ei = 0.0;
	for (std::map<std::set<int>, std::vector<int> >::iterator it = s.begin(); it != s.end(); it++) {
		Ei += pair_energy_farago(current, it->second);
	}

	// swap patches
	current.patch_swap(bounds1, bounds2, size);

	// calculate final energy, energy difference
	double Ef = 0.0;
	for (std::map<std::set<int>, std::vector<int> >::iterator it = s.begin(); it != s.end(); it++) {
		Ef += pair_energy_farago(current, it->second);
	}
	double delE = Ef - Ei;

	// accept/reject attempted move
	bool accept = metropolis_accept(delE);

	// return system to initial state if move is rejected
	if (accept == false) {
		current.patch_swap(bounds1, bounds2, size);
		// log acceptance ratio
		patch_accept[size][1] += 1.0;

		return 0.0;
	}
	else {
		// log acceptance ratio
		patch_accept[size][0] += 1.0; patch_accept[size][1] += 1.0;

		return delE;
	}
}


void site_picker_DPPC(membrane& leaflet, std::vector<std::vector<int> >& sites) {

	int n = sys.get_grid_size();

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


void patch_sites(membrane& current, const std::vector<int>& a, int patch_size, std::vector<std::vector<int> >& patch, std::vector<int>& bounds) {

	int n = sys.get_grid_size();

	int shift = patch_size / 2;

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

	for (int i = 0; i < patch_size; i++) {
		int x = (bounds[0] + i) % n;
		for (int j = 0; j < patch_size; j++) {
			int y = (bounds[2] + j) % n;
			patch.push_back(std::vector<int>{x, y});
		}
	}
}


void patch_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int& patch_size) {

	int n = sys.get_grid_size();
	// pick swappable patches with mathcing indentation counts
	int size; bool match = false;

	while (match!=true){
		// sample a patch size from input distribution
		size = sample_swap_size();
		assert(size > 1 && "Invalid patch size used. Patches must be larger than 1 lattice site\n");

		int pick_attempts = 0;
		while (pick_attempts < 0.5*(n*n) && match!=true) {
			// pick two patch centers
			patch_center_picker(leaflet, c1, c2, size);
			match = indentation_match(leaflet, c1, c2, size);
			if (match == false) {
				pick_attempts += 1;
			}
		}
	}
}


void patch_center_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int patch_size) {

	int n = sys.get_grid_size();

	// pick first point
	for (int i = 0; i < 2; i++) {
		c1[i] = rand_int(0, n - 1);
	}

	// calculate invalid area to avoid picking overlapping patches
	int invalid_bounds[4] = { 0, 0, 0, 0 };
	if (patch_size % 2 == 0) {
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


bool indentation_match(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int patch_size) {

	return true;
}


void indentation_sections(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, std::vector<std::vector<int> >& indent_sections) {

}