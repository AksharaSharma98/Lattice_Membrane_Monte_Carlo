#ifndef evolution_h
#define evolution_h


void evolve_mc(membrane& upper, membrane& lower, int steps, int energy_output_freq, int config_output_freq);

void lipid_picker(membrane leaflet, int c1[2], int c2[2]);

double monte_carlo_move(membrane& current, membrane& opposing, int* a, int* b);

bool metropolis_accept(double delE);

#endif
