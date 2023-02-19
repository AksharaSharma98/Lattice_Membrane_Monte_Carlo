#include <math.h>
#include <string>
#include <tuple>
#include <assert.h>

#include "energy.h"
#include "lipid.h"
#include "membrane.h"
#include "parameters.h"
#include "periodicboundary.h"
#include "Lattice_Membrane_Monte_Carlo.h"


double system_energy(membrane upper, membrane lower) {

	int n = upper.getgrid().size();
	double energy = 0.0;
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int x[2] = { i,j };
			energy += local_energy(upper, lower, x);
			energy += local_energy(lower, upper, x);
		}
	}

	return energy/2.0;
	// *** Very Important ***
	// verify if /2.0 is still correct after any Hamiltonian modificaitons
}


double local_energy(membrane current, membrane opposing, int* x) {

	std::string c, o, cn, on;
	double sc, so, scn;
	double loc_e_plane = 0.0, loc_e_inter = 0.0, epsc = 0.0, epso = 0.0, epsi = 0.0;
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	periodic_neighbours(current, x[0], x[1], nbs);

	c = current.getlipid(x[0], x[1]).getspecies();
	o = opposing.getlipid(x[0], x[1]).getspecies();
	sc = current.getlipid(x[0], x[1]).gettail_order();
	so = opposing.getlipid(x[0], x[1]).gettail_order();
	epsi = forcefield.getinter_pair_energy(std::make_pair(c, o));
	for (int k = 0; k < 6; k++) {
		cn = current.getlipid(nbs[k][0], nbs[k][1]).getspecies();
		on = opposing.getlipid(nbs[k][0], nbs[k][1]).getspecies();
		scn = current.getlipid(nbs[k][0], nbs[k][1]).gettail_order();
		epsc = forcefield.getplane_pair_energy(std::make_pair(c, cn));
		loc_e_plane += epsc * sc * scn;                                         // e_i = eps_ij*s_i*s_j
	}
	loc_e_inter += epsi * sc * so;                                          // e_p = eps'_ii*s_i*s'_i

	return loc_e_plane +loc_e_inter;
}