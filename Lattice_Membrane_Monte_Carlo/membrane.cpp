#include <string>
#include <vector>
#include <assert.h>

#include "lipid.h"
#include "math_functions.h"
#include "membrane.h"


// default constructor

membrane::membrane (int n, std::string* sp, int* pop, double* tail, int sp_num)
{
	assert(n%2 == 0 && "Grid size must be even");
	assert(n >= 1 && "Invalid grid size");
	size = n;

	// initialize grid of lipids

	std::vector<int> count;
	for (int i = 0; i < sp_num;i++) {
		count.push_back(0);
	}

	for (int i = 0; i < n; i++) {
		grid.push_back(Array());
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int t = 0, l = 0;
			while (t == 0) {
				l = rand_int(0, sp_num-1);
				if (count[l] < pop[l]) {
					t = 1;
				}
			}
			int position[2] = { i,j };
			lipid a(sp[l], tail[l], position);
			grid[i].push_back(a);
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