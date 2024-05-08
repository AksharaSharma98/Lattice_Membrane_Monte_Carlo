#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "lipid.h"
#include "membrane.h"
#include "initialize.h"
#include "Lattice_Membrane_Monte_Carlo.h"


std::pair<membrane, membrane> initialize() {
	membrane upper(0);
	membrane lower(1);
	return std::make_pair(upper, lower);
}