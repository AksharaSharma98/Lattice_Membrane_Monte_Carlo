#ifndef energy_h
#define energy_h

#include "membrane.h"


double system_energy_farago(membrane& current);

double local_energy_farago(membrane& current, std::vector<int>& x);

double pair_energy_farago(membrane& current, std::vector<int>& x);

double entropy_farago(membrane& current, std::vector<int>& x);

double system_energy(membrane& upper, membrane& lower);

double local_enthalpy(membrane& current, membrane& opposing, int* x);

double local_planeentropy_lipid(membrane& current, int* x);

double local_planeentropy_env(membrane& current, int* x);

double local_interentropy_lipid(membrane& current, membrane& opposing, int* x);

double local_interentropy_env(membrane& current, membrane& opposing, int* x);

#endif
