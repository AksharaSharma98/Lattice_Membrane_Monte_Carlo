#ifndef energy_h
#define energy_h

#include "membrane.h"


double system_energy(membrane& upper, membrane& lower);

double local_enthalpy(membrane& current, membrane& opposing, int* x);

double local_entropy_lipid(membrane& current, membrane& opposing, int* x);

double local_entropy_env(membrane& current, membrane& opposing, int* x);

#endif
