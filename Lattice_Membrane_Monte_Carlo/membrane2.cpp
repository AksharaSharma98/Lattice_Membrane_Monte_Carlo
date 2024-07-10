#include <string>
#include <vector>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <iostream>

#include "lipid2.h"
#include "math_functions.h"
#include "periodicboundary.h"
#include "qol_functions.h"
#include "membrane2.h"
#include "Lattice_Membrane_Monte_Carlo.h"


// default constructor

membrane2::membrane2(int leaflet) {
	
	Grid_size = sys.get_grid_size();
	Leaflet_index = leaflet;

	assert(Grid_size % 2 == 0 && "Lattice size must be even");
	assert(Grid_size >= 1 && "Invalid grid size");

	int num_species = sys.get_num_species();

	// initialize Grid of lipids, Lipid_lists
	for (int i = 0; i < Grid_size; i++) {
		Grid.push_back(std::vector<int>(Grid_size, -1));
	}
	for (int i = 0; i < num_species; i++) {
		Lipid_lists.push_back(lipid_list());
	}
	
	// initialize lattice matrix

	std::ifstream restart_file;
	restart_file.open("restart.bin", std::ios::in | std::ios::binary);

	// initialize from restart file if it exists
	if (restart_file) {

		// read leaflet config using grid size as guide, assuming data is in short int form

		if(leaflet==1){ // skip over upper leaflet config
			restart_file.seekg(Grid_size * Grid_size * sizeof(short int));
		}
		// populate lattice matrix
		for (int i = 0; i < Grid_size; i++) {
			for (int j = 0; j < Grid_size; j++) {
				short int temp;
				restart_file.read(reinterpret_cast<char*> (&temp), sizeof(short int));
				Grid[i].push_back((int)temp);
			}
		}
		if (leaflet == 0) { // skip over lower leaflet config
			restart_file.seekg(2 * Grid_size * Grid_size * sizeof(short int));
		}

		// read connectivity info header for appropriate leaflet
		int current_pos = 2 * Grid_size * Grid_size * sizeof(short int);
		
		int num_two_site_upper, num_one_site_upper, num_two_site_lower, num_one_site_lower;
		restart_file.read(reinterpret_cast<char*> (&num_two_site_upper), sizeof(int));
		restart_file.read(reinterpret_cast<char*> (&num_one_site_upper), sizeof(int));
		restart_file.read(reinterpret_cast<char*> (&num_two_site_lower), sizeof(int));
		restart_file.read(reinterpret_cast<char*> (&num_one_site_lower), sizeof(int));
		std::vector<int> num_two_site = { num_two_site_upper, num_two_site_lower };
		std::vector<int> num_one_site = { num_one_site_upper, num_one_site_lower };

		// read connectivity data
		current_pos = current_pos + 4 * sizeof(int);
		
		restart_file.seekg(current_pos + leaflet * (2 * num_two_site_upper + num_one_site_upper) * sizeof(int));
		for (int i = 0; i < num_two_site[leaflet]; i++) {
			int index1, index2;
			restart_file.read(reinterpret_cast<char*> (&index1), sizeof(int));
			restart_file.read(reinterpret_cast<char*> (&index2), sizeof(int));

			// initialize lipid object
			std::vector<int> coords = coords_from_index(index1, Grid_size);
			std::vector<int> coords2 = coords_from_index(index2, Grid_size);
			coords.insert(coords.end(), coords2.begin(), coords2.end());
			int sp = Grid[coords[0]][coords[1]];
			std::vector<double> s = { sample_tailorder(sys.get_species(leaflet, sp), -1.0) , sample_tailorder(sys.get_species(leaflet, sp), -1.0) };
			std::vector<int> indices = { index1, index2 };
			lipid2 a(sys.get_species(leaflet, sp), s, coords, indices);

			// insert lipid into corresponding lipid list
			Lipid_lists[sp].push_back(a);
		}
		for (int i = 0; i < num_one_site[leaflet]; i++) {
			int index;
			restart_file.read(reinterpret_cast<char*> (&index), sizeof(int));

			// initialize lipid object
			std::vector<int> coords = coords_from_index(index, Grid_size);
			int sp = Grid[coords[0]][coords[1]];
			std::vector<double> s = { sample_tailorder(sys.get_species(leaflet, sp), -1.0) };
			std::vector<int> indices = { index };
			lipid2 a(sys.get_species(leaflet, sp), s, coords, indices);

			// insert lipid into corresponding lipid list
			Lipid_lists[sp].push_back(a);
		}
		
	}

	// initialize by randomly distributing lipids in absence of restart file
	else {
		// create lists of edges and lattice sites
		std::vector<std::vector<int> > edges = create_edge_list(Grid_size);
		std::vector<int> sites;
		for (int i = 0; i < Grid_size; i++) {
			sites.push_back(i);
		}
		// source this from parameter file instead -------------------------------------TEMPORARY//
		std::vector<int> num_lipid_sites = { 1, 2, 2, 1, 2 };

		// create list of lipids that have to be placed on the lattice
		std::vector<std::vector<int> > lipids;
		for (int i = 0; i < num_species; i++) {
			std::vector<int> temp;
			for (int k = 0; k < num_lipid_sites[i]; k++) {
				temp.push_back(i);
			}
			int pop = sys.get_population(leaflet, i);
			for (int j = 0; j < pop; j++) {
				lipids.push_back(temp);
			}
		}

		while (sites.size() != 0) {
			// randomly select a lipid from the total list to place on the grid
			int l = rand_int(0, lipids.size() - 1);
			// pop lipid from list
			std::vector<int> temp_lipid = lipids[l];
			lipids.erase(lipids.begin() + l);

			// place on first available space from edge/site lists based on number of sites needed by lipid
			if (temp_lipid.size() == 2) {
				// initialize lipid object
				std::vector<double> s = { sample_tailorder(sys.get_species(leaflet, l), -1.0) , sample_tailorder(sys.get_species(leaflet, l), -1.0) };
				std::vector<int> coords = coords_from_index(edges[0][0], Grid_size);
				std::vector<int> coords2 = coords_from_index(edges[0][1], Grid_size);
				coords.insert(coords.end(), coords2.begin(), coords2.end());
				std::vector<int> indices = { edges[0][0], edges[0][1] };
				lipid2 a(sys.get_species(leaflet, temp_lipid[0]), s, coords, indices);

				// insert lipid into corresponding lipid list, populate lattice matrix
				Lipid_lists[temp_lipid[0]].push_back(a);
				Grid[coords[0]][coords[1]] = temp_lipid[0];
				Grid[coords[2]][coords[3]] = temp_lipid[0];

				// pop occupied and obstructed edges, indices from respective lists
				int index1 = edges[0][0]; int index2 = edges[0][1];
				sites.erase(std::find(sites.begin(), sites.end(), index1));
				sites.erase(std::find(sites.begin(), sites.end(), index2));
				std::vector<int> temp;
				for (int i = 0; i < edges.size(); i++) {
					if (edges[i][0] == index1 || edges[i][1] == index1) {
						temp.push_back(i);
					}
					if (edges[i][0] == index2 || edges[i][1] == index2) {
						temp.push_back(i);
					}
				}
				for (int i = temp.size() - 1; i >= 0; i--) {
					edges.erase(edges.begin() + temp[i]);
				}
			}
			else {
				// initialize lipid object
				std::vector<double> s = { sample_tailorder(sys.get_species(leaflet, l), -1.0) };
				std::vector<int> coords = coords_from_index(edges[0][0], Grid_size);
				std::vector<int> indices = { edges[0][0] };
				lipid2 a(sys.get_species(leaflet, temp_lipid[0]), s, coords, indices);
				
				// insert lipid into corresponding lipid list, populate lattice matrix
				Lipid_lists[temp_lipid[0]].push_back(a);
				Grid[coords[0]][coords[1]] = temp_lipid[0];

				// pop occupied and obstructed edges, indices from respective lists
				int index = edges[0][0];
				sites.erase(std::find(sites.begin(), sites.end(), index));
				std::vector<int> temp;
				for (int i = 0; i < edges.size(); i++) {
					if (edges[i][0] == index || edges[i][1] == index) {
						temp.push_back(i);
					}
				}
				for (int i = temp.size() - 1; i >= 0; i--) {
					edges.erase(edges.begin() + temp[i]);
				}
			}
		}
	}

}


// default destructor

membrane2::~membrane2() {

}


// accessor functions

int membrane2::grid_size() {
	return Grid_size;
}


int membrane2::leaflet_index() {
	return Leaflet_index;
}


membrane2::lipid_list& membrane2::species_list(int species) {
	return Lipid_lists[species];
}


std::vector<std::vector<int> > membrane2::get_lattice() {
	return Grid;
}


// modifier functions

void membrane2::swap(std::vector<std::vector<int> > list_indices) {
	int x1 = list_indices[0][0]; int y1 = list_indices[0][1];
	int x2 = list_indices[1][0]; int y2 = list_indices[1][1];

	lipid2 temp = Lipid_lists[x1][y1];
	Lipid_lists[x1][y1].update_coords(Lipid_lists[x2][y2].coordinates());
	Lipid_lists[x1][y1].update_indices(Lipid_lists[x2][y2].indices());
	Lipid_lists[x2][y2].update_coords(temp.coordinates());
	Lipid_lists[x2][y2].update_indices(temp.indices());
}


void membrane2::swap_DPPC_state(std::vector<int> list_indices) {
	int x = list_indices[0]; int y = list_indices[1];
	if (Lipid_lists[x][y].species() == "DPPCd") {
		Lipid_lists[x][y].update_species("DPPCo");
	}
	else {
		Lipid_lists[x][y].update_species("DPPCd");
	}
}