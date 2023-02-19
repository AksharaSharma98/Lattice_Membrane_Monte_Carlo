#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>

#include "math_functions.h"

// set RNG seed
std::mt19937 mt(1790);
std::uniform_real_distribution<double> uniform(0.0, 1.0);


double rand01()
{
	return uniform(mt);
}

int rand_int(int lower_bound, int upper_bound)
{
	std::uniform_int_distribution<int> uniform_int(lower_bound, upper_bound);
	return uniform_int(mt);
}