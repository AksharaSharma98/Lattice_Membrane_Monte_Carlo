#include <math.h>
#include <string>
#include <assert.h>

#include "energy.h"
#include "lipid.h"
#include "membrane.h"
#include "parameters.h"
#include "periodicboundary.h"
#include "Lattice_Membrane_Monte_Carlo.h"


double system_energy(membrane upper, membrane lower) {

	int n = upper.getgrid().size();
	double energy = 0.0, epsu = 0.0, epsl = 0.0;
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			periodic_neighbours(upper, i, j, nbs);
			std::string u = upper.getlipid(i, j).getspecies();
			std::string l = lower.getlipid(i, j).getspecies();
			for (int k = 0; k < 6; k++) {
				std::string un = upper.getlipid(nbs[k][0], nbs[k][1]).getspecies();
				std::string ln = lower.getlipid(nbs[k][0], nbs[k][1]).getspecies();
				epsu = forcefield.getpair_energy(std::make_pair(u, un));
				epsl = forcefield.getpair_energy(std::make_pair(l, ln));
				energy += epsu + epsl;
			}
			// need to add leaflet coupling terms later
		}
	}

	return energy/2.0; // verify if /2.0 is still correct after any Hamiltonian modificaitons
}


double local_energy(membrane current, membrane opposing, int* a, int* b) {

	double loc_energy = 0.0, epsc = 0.0, epso = 0.0;
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };

	periodic_neighbours(current, a[0], a[1], nbs);
	std::string c = current.getlipid(a[0], a[1]).getspecies();
	std::string o = opposing.getlipid(a[0], a[1]).getspecies();
	for (int k = 0; k < 6; k++) {
		std::string cn = current.getlipid(nbs[k][0], nbs[k][1]).getspecies();
		std::string on = opposing.getlipid(nbs[k][0], nbs[k][1]).getspecies();
		epsc = forcefield.getpair_energy(std::make_pair(c, cn));
		epso = forcefield.getpair_energy(std::make_pair(o, on));
		loc_energy += epsc + epso;
	}

	periodic_neighbours(current, b[0], b[1], nbs);
	c = current.getlipid(b[0], b[1]).getspecies();
	o = opposing.getlipid(b[0], b[1]).getspecies();
	for (int k = 0; k < 6; k++) {
		std::string cn = current.getlipid(nbs[k][0], nbs[k][1]).getspecies();
		std::string on = opposing.getlipid(nbs[k][0], nbs[k][1]).getspecies();
		epsc = forcefield.getpair_energy(std::make_pair(c, cn));
		epso = forcefield.getpair_energy(std::make_pair(o, on));
		loc_energy += epsc + epso;
	}

	return loc_energy;
}