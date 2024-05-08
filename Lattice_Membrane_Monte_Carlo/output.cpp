#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "lipid.h"
#include "membrane.h"
#include "output.h"
#include "Lattice_Membrane_Monte_Carlo.h"


void write_config_int(FILE* config_file, membrane& upper, membrane& lower) {

	int n = upper.getgrid().size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			lipid l = upper.getlipid(i, j);
			fprintf(config_file, "%d ", sys.get_output_type(l.getspecies()));
		}
		fprintf(config_file, " ");
		for (int j = 0; j < n; j++) {
			lipid l = lower.getlipid(i, j);
			fprintf(config_file, "%d ", sys.get_output_type(l.getspecies()));
		}
		fprintf(config_file, "\n");
	}
}


void write_config_bin(std::ofstream& config_bin_file, membrane& upper, membrane& lower) {

	int n = upper.getgrid().size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			lipid l = upper.getlipid(i, j);
			short int site = sys.get_output_type(l.getspecies());
			config_bin_file.write(reinterpret_cast<char*> (&site), sizeof(site));
		}
		for (int j = 0; j < n; j++) {
			lipid l = lower.getlipid(i, j);
			short int site = sys.get_output_type(l.getspecies());
			config_bin_file.write(reinterpret_cast<char*> (&site), sizeof(site));
		}
	}
}


void write_tailconfig(FILE* tailconfig_file, membrane& upper, membrane& lower) {

	int n = upper.getgrid().size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fprintf(tailconfig_file, "%.3f ", upper.getlipid(i, j).gettail_order());
		}
		fprintf(tailconfig_file, " ");
		for (int j = 0; j < n; j++) {
			fprintf(tailconfig_file, "%.3f ", lower.getlipid(i, j).gettail_order());
		}
		fprintf(tailconfig_file, "\n");
	}
}


void write_config_species(FILE* config_file, membrane& upper, membrane& lower) {

	int n = upper.getgrid().size();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			lipid l = upper.getlipid(i, j);
			fprintf(config_file, "%s ", l.getspecies().c_str());
		}
		fprintf(config_file, " ");
		for (int j = 0; j < n; j++) {
			lipid l = lower.getlipid(i, j);
			fprintf(config_file, "%s ", l.getspecies().c_str());
		}
		fprintf(config_file, "\n");
	}
}


void write_energy(FILE* energy_file, double energy) {

	fprintf(energy_file, "%lf\n", energy);

}


void write_header(FILE* energy_file) {
	// this is very crude and needs proper rewriting
	fprintf(energy_file, "# ");
	fprintf(energy_file, "%lf ", forcefield.getplane_pair_energy(std::make_pair("DPPCo", "DPPCo")));
	fprintf(energy_file, "%lf ", forcefield.getplane_pair_energy(std::make_pair("DPPCo", "CHOL")));
	fprintf(energy_file, "%lf ", forcefield.getplane_pair_energy(std::make_pair("DPPCo", "DOPC")));
	fprintf(energy_file, "\n");
}


void print_leaflet_species(membrane& leaflet) {

	int n = leaflet.getgrid().size();

	printf("\n");
	for (int i = 0; i < n; i++) {
		printf("\n");
		if (i % 2 == 0) {
			printf("   ");
		}
		for (int j = 0; j < n; j++) {
			lipid l = leaflet.getlipid(i, j);
			printf("%s  ", l.getspecies().c_str());
		}
		printf("\n");
	}
	printf("\n");
}