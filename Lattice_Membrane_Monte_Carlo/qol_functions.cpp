#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "lipid.h"
#include "membrane.h"
#include "periodicboundary.h"
#include "qol_functions.h"


std::vector<int> coords_from_index(int index, int grid_size) {
	return std::vector<int> {index / grid_size, index% grid_size};
}


std::vector<std::vector<int> > create_edge_list(int grid_size) {
	std::vector<std::vector<int> > edge_list;
	int nbs[6][2] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };
	for (int i = 0; i < grid_size; i++) {
		for (int j = 0; j < grid_size; j++) {
			int current_index = i * grid_size + j;
			periodic_neighbours(grid_size, i, j, nbs);
			int neighbour_indices[6] = {0,0,0,0,0,0};
			for (int k = 0; k < 6; k++) {
				neighbour_indices[k] = nbs[k][0] * grid_size + nbs[k][1];
				if (current_index < neighbour_indices[k]) {
					edge_list.push_back({ current_index, neighbour_indices[k] });
				}
			}
		}
	}

	return edge_list;
}