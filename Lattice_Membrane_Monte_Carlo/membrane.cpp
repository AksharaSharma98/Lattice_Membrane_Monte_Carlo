#include <string>
#include <vector>
#include <assert.h>

#include "lipid.h"
#include "math_functions.h"
#include "membrane.h"
#include "Lattice_Membrane_Monte_Carlo.h"


// default constructor

membrane::membrane (int leaflet)
{
	size = sys.get_grid_size();
	
	int num_species = sys.get_num_species();
	
	assert(size%2 == 0 && "Grid size must be even");
	assert(size >= 1 && "Invalid grid size");

	// initialize grid of lipids

	std::vector<std::vector<int> > matrix(size, std::vector<int>(size, -1));
	std::vector<int> count;
	for (int i = 0; i < num_species; i++) {
		count.push_back(0);
	}

	bool full = false;
	int x, y, l, total = 0;
	while (full == false) {
		bool valid_position = false;
		while (valid_position == false) {
			x = rand_int(0, size - 1);
			y = rand_int(0, size - 1);
			if (matrix[x][y] == -1) {
				valid_position = true;
			}
		}
		bool valid_species = false;
		while (valid_species == false) {
			l = rand_int(0, num_species - 1);
			if (count[l] < sys.get_population(leaflet,l)) {
				valid_species = true;
			}
		}
		if (valid_position == true && valid_species == true) {
			matrix[x][y] = l;
			count[l] += 1;
			total += 1;
		}
		if (total == size*size) {
			full = true;
		}
	}

	for (int i = 0; i < size; i++) {
		grid.push_back(Array());
	}
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			l = matrix[i][j];
			int position[2] = { i,j };
			double s = sample_tailorder(sys.get_species(leaflet,l), -1.0);
			lipid a(sys.get_species(leaflet,l), s, position);
			grid[i].push_back(a);  // vector.push_back(value) appends value at end of vector
			count[l] += 1;
		}
	}
}


// default destructor

membrane::~membrane() {

}


// accessor functions

int membrane::getsize() {
	return size;
}


membrane::Grid& membrane::getgrid() {
	return grid;
}


lipid& membrane::getlipid(int i, int j) {
	return grid[i][j];
}


// modifier functions

void membrane::swap(std::vector<int> &a, std::vector<int> &b) {
	lipid l = grid[a[0]][a[1]];
	grid[a[0]][a[1]] = grid[b[0]][b[1]];
	grid[b[0]][b[1]] = l;
}


void membrane::swap_DPPC_state(std::vector<int>& a) {
	lipid l = grid[a[0]][a[1]];
	int position[2] = { a[0], a[1] };
	
	if (l.getspecies() == "DPPCd") {
		lipid l1("DPPCo", 0.0, position);
		grid[a[0]][a[1]] = l1;
	}
	else {
		lipid l1("DPPCd", 0.0, position);
		grid[a[0]][a[1]] = l1;
	}
}


void membrane::patch_swap(std::vector<int>& bounds1, std::vector<int>& bounds2, int patch_size) {
	
	// iterate within patch bounds and swap lipid by lipid
	for (int i = 0; i < patch_size; i++) {
		int x1 = (bounds1[0] + i) % size;
		int x2 = (bounds2[0] + i) % size;

		for (int j = 0; j < patch_size; j++) {
			int y1 = (bounds1[2] + j) % size;
			int y2 = (bounds2[2] + j) % size;
			
			lipid l = grid[x1][y1];
			grid[x1][y1] = grid[x2][y2];
			grid[x2][y2] = l;
		}
	}

	// For overlapping patches, incomplete: Needs resolution scheme for overlapping region(s)
	/*// Create a deep copy of both patches
	Grid patch1, patch2;

	for (int i = 0; i < patch_size; i++) {
		int x1 = (bounds1[0] + i) % size;
		int x2 = (bounds2[0] + i) % size;
		patch1.push_back(Array());
		patch2.push_back(Array());
		for (int j = 0; j < patch_size; j++) {
			int y1 = (bounds1[2] + j) % size;
			int y2 = (bounds2[2] + j) % size;

			patch1[i].push_back(grid[x1][y1]);
			patch2[i].push_back(grid[x2][y2]);
		}
	}
	
	// swap the patches using the copies
	for (int i = 0; i < patch_size; i++) {
		int x1 = (bounds1[0] + i) % size;
		int x2 = (bounds2[0] + i) % size;
		
		for (int j = 0; j < patch_size; j++) {
			int y1 = (bounds1[2] + j) % size;
			int y2 = (bounds2[2] + j) % size;
			
			grid[x1][y1] = patch2[i][j];
			grid[x2][y2] = patch1[i][j];
		}
	}*/

}


void membrane::tail_update(int* a, double s_new) {
	grid[a[0]][a[1]].update_tailorder(s_new);
}