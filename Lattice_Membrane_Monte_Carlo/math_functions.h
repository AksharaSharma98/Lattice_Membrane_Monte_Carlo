#ifndef math_h
#define math_h

#include "membrane.h"

double rand01();

int rand_int(int lower_bound, int upper_bound);

double sample_tailorder(std::string species, double s_old);

void mole_fraction(membrane& current, membrane& opposing, int* x, std::vector<double>& phi, int nbs[][2]);

#endif
