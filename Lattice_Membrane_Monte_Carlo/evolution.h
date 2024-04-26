#ifndef evolution_h
#define evolution_h


void evolve_mc_farago(membrane& upper, membrane& lower, int steps, int energy_output_freq, int config_output_freq, int restart_output_freq);

//void evolve_mc(membrane& upper, membrane& lower, int steps, int tailorder_update_freq, int energy_output_freq, int config_output_freq);

//double update_tailorder_mc(membrane& current, membrane& opposing, int* a);

#endif
