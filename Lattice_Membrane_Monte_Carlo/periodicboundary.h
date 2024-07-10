#ifndef periodic_bound_h
#define periodic_bound_h


void periodic_neighbours(membrane& leaflet, int i, int j, int nbs[][2]);

void periodic_neighbours(int grid_size, int i, int j, int nbs[][2]);

#endif