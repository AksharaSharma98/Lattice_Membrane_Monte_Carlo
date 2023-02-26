#include <string>
#include <vector>
#include <assert.h>

#include "lipid.h"
#include "math_functions.h"
#include "membrane.h"
#include "Lattice_Membrane_Monte_Carlo.h"


// default constructor

membrane::membrane (int n, std::string* sp, int* pop)
{
	assert(n%2 == 0 && "Grid size must be even");
	assert(n >= 1 && "Invalid grid size");
	size = n;

	// initialize grid of lipids

	std::vector<std::vector<int>> matrix(n, std::vector<int>(n, -1));
	std::vector<int> count;
	for (int i = 0; i < n_sp; i++) {
		count.push_back(0);
	}

	bool full = false;
	int x, y, l, total = 0;
	while (full == false) {
		bool valid_position = false;
		while (valid_position == false) {
			x = rand_int(0, n - 1);
			y = rand_int(0, n - 1);
			if (matrix[x][y] == -1) {
				valid_position = true;
			}
		}
		bool valid_species = false;
		while (valid_species == false) {
			l = rand_int(0, n_sp - 1);
			if (count[l] < pop[l]) {
				valid_species = true;
			}
		}
		if (valid_position == true && valid_species == true) {
			matrix[x][y] = l;
			count[l] += 1;
			total += 1;
		}
		if (total == n * n) {
			full = true;
		}
	}

	for (int i = 0; i < n; i++) {
		grid.push_back(Array());
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			l = matrix[i][j];
			int position[2] = { i,j };
			double s = sample_tailorder(sp[l], -1.0);
			lipid a(sp[l], s, position);
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

void membrane::swap(int* a, int* b) {
	lipid l = grid[a[0]][a[1]];
	grid[a[0]][a[1]] = grid[b[0]][b[1]];
	grid[b[0]][b[1]] = l;
}


void membrane::tail_update(int* a, double s_new) {
	grid[a[0]][a[1]].update_tailorder(s_new);
}