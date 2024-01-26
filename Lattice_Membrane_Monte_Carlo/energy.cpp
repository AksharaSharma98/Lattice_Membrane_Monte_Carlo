#include <math.h>
#include <string>
#include <tuple>
#include <assert.h>

#include "energy.h"
#include "lipid.h"
#include "membrane.h"
#include "parameters.h"
#include "system.h"
#include "periodicboundary.h"
#include "math_functions.h"
#include "Lattice_Membrane_Monte_Carlo.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>


double system_energy_farago(membrane& current) {

	int n = current.getgrid().size();
	double energy = 0.0, local_energy = 0.0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::vector x = { i,j };
			lipid lc = current.getlipid(x[0], x[1]);
			if (lc.getspecies() == "DPPCd") { energy -= (forcefield.getplane_entropy_const() * kT); }
			local_energy += local_energy_farago(current, x);
		}
	}

	return energy + (0.5 * local_energy);
	// *** Very Important ***
	// verify if *0.5 is still correct after any Hamiltonian modificaitons
}


double local_energy_farago(membrane& current, std::vector<int>& x) {

	double local_energy = 0.0;
	
	// get list of neighbours with PBC correction
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };
	periodic_neighbours(current, x[0], x[1], nbs);
	
	lipid lc = current.getlipid(x[0], x[1]);
	std::string c = lc.getspecies();

	for (int k = 0; k < 6; k++) {
		lipid l = current.getlipid(nbs[k][0], nbs[k][1]);
		std::string cn = l.getspecies();
		if (c == "DPPCo" || cn == "DPPCo") {
			local_energy -= forcefield.getplane_pair_energy(std::make_pair(c, cn));
		}
	}
	
	return local_energy;
}


double pair_energy_farago(membrane& current, std::vector<int>& x) {
	
	double pair_energy = 0.0;

	lipid l1 = current.getlipid(x[0], x[1]);
	lipid l2 = current.getlipid(x[2], x[3]);
	std::string s1 = l1.getspecies();
	std::string s2 = l2.getspecies();

	if (s1 == "DPPCo" || s2 == "DPPCo") {
		pair_energy -= forcefield.getplane_pair_energy(std::make_pair(s1, s2));
	}

	return pair_energy;
}


double entropy_farago(membrane& current, std::vector<int>& x) {
	lipid lc = current.getlipid(x[0], x[1]);
	if (lc.getspecies() == "DPPCd") { 
		return -(forcefield.getplane_entropy_const() * kT); 
	}
	else {
		return 0.0;
	}
}


double system_energy(membrane& upper, membrane& lower) {

	int n = upper.getgrid().size();
	double enthalpy = 0.0, entropy = 0.0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int x[2] = { i,j };
			enthalpy += local_enthalpy(upper, lower, x) + local_enthalpy(lower, upper, x);
			entropy += local_planeentropy_lipid(upper, x) + local_planeentropy_lipid(lower, x);
			entropy += 2 * local_interentropy_lipid(lower, upper, x);
		}
	}
	
	return 0.5*enthalpy + entropy;
	// *** Very Important ***
	// verify if *0.5 is still correct after any Hamiltonian modificaitons
}


double local_enthalpy(membrane& current, membrane& opposing, int* x) {

	std::string c, o, cn;
	double sc, so, scn;
	double loc_e_plane = 0.0, loc_e_inter = 0.0, epsc = 0.0, epso = 0.0, epsi = 0.0;
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	periodic_neighbours(current, x[0], x[1], nbs);

	lipid lc = current.getlipid(x[0], x[1]);
	lipid lo = opposing.getlipid(x[0], x[1]);
	c = lc.getspecies();
	o = lo.getspecies();
	sc = lc.gettail_order();
	so = lo.gettail_order();
	epsi = forcefield.getinter_pair_energy(std::make_pair(c, o));
	for (int k = 0; k < 6; k++) {
		lipid l = current.getlipid(nbs[k][0], nbs[k][1]);
		cn = l.getspecies();
		scn = l.gettail_order();
		epsc = forcefield.getplane_pair_energy(std::make_pair(c, cn));
		loc_e_plane += epsc * sc * scn;                                         // e_i = eps_ij*s_i*s_j
	}
	loc_e_inter += epsi* sc* so;                                              // e_p = eps'_ii*s_i*s'_i
	
	return loc_e_plane + loc_e_inter;
}

double local_planeentropy_lipid(membrane& current, int* x) {

	double temp = 0.0;
	double epss = forcefield.getplane_entropy_const();
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	periodic_neighbours(current, x[0], x[1], nbs);

	std::vector<double> phi(sys.get_num_species(), 0.0);
	mole_fraction_plane(current, x, phi, nbs);
	for (int i = 0; i < sys.get_num_species(); i++) {
		if (phi[i] == 0.0) {
			temp += 0.0;
		}
		else {
			temp += phi[i] * log(phi[i]);
		}
	}

	return epss * kT * temp;
}

double local_planeentropy_env(membrane& current, int* x) {

	double loc_e_plane = 0.0;
	double epss = forcefield.getplane_entropy_const();
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	periodic_neighbours(current, x[0], x[1], nbs);

	for (int i = -1; i < 6; i++) {
		if (i == -1) {
			loc_e_plane += local_planeentropy_lipid(current, x);
		}
		else {
			loc_e_plane += local_planeentropy_lipid(current, nbs[i]);
		}
	}

	return loc_e_plane;
}

double local_interentropy_lipid(membrane& current, membrane& opposing, int* x) {

	double temp = 0.0;
	double epss = forcefield.getinter_entropy_const();
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	periodic_neighbours(current, x[0], x[1], nbs);

	std::vector<double> phi(sys.get_num_species(), 0.0);
	mole_fraction_total(current, opposing, x, phi, nbs);
	for (int i = 0; i < sys.get_num_species(); i++) {
		if (phi[i] == 0.0) {
			temp += 0.0;
		}
		else {
			temp += phi[i] * log(phi[i]);
		}
	}

	return epss*kT*temp;
}

double local_interentropy_env(membrane& current, membrane& opposing, int* x) {

	double loc_e_plane = 0.0;
	double epss = forcefield.getinter_entropy_const();
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	periodic_neighbours(current, x[0], x[1], nbs);

	for (int i = -1; i < 6; i++) {
		if (i == -1) {
			loc_e_plane += 2*local_interentropy_lipid(current, opposing, x);
		}
		else {
			loc_e_plane += 2*local_interentropy_lipid(current, opposing, nbs[i]);
		}
	}

	return loc_e_plane;
}