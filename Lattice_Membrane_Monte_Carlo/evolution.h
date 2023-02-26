#ifndef evolution_h
#define evolution_h


void evolve_mc(membrane& upper, membrane& lower, int steps, int tailorder_update_freq, int energy_output_freq, int config_output_freq);

void lipid_picker(membrane& leaflet, int c1[2], int c2[2]);

void lipid_picker_species(membrane& leaflet, int c1[2], int c2[2]);

double monte_carlo_move(membrane& current, membrane& opposing, int* a, int* b);// , std::vector<std::vector<int>>& swap_sites);

bool metropolis_accept(double delE);

double update_tailorder_mc(membrane& current, membrane& opposing, int* a);

void update_system_tailorder(membrane& current, membrane& opposing, std::vector<std::vector<int>>& swap_sites);

void update_sequence(int n, std::vector<std::vector<int>> swap_sites, std::vector<std::vector<int>>& seq);

#endif
