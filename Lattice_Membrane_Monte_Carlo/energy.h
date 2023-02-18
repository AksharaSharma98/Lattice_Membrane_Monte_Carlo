#ifndef energy_h
#define energy_h

#include "membrane.h"


double system_energy(membrane upper, membrane lower);

double local_energy(membrane current, membrane opposing, int* a, int* b);

#endif
