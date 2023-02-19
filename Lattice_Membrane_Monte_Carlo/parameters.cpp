#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <iostream>

#include "parameters.h"


// default constructor

parameters::parameters() {

	// Note: Both possible combinations of a pair-key must be added unless they are the same
	plane_pair_energy = { {std::make_pair("DPPC","DPPC"), -2.0},
						  {std::make_pair("DIPC","DIPC"), -1.0},
						  {std::make_pair("DPPC","DIPC"), -1.2},
						  {std::make_pair("DIPC","DPPC"), -1.2} };

	inter_pair_energy = { {std::make_pair("DPPC","DPPC"), -1.0},
						  {std::make_pair("DIPC","DIPC"), -2.0},
						  {std::make_pair("DPPC","DIPC"), -1.5},
						  {std::make_pair("DIPC","DPPC"), -1.5} };

	output_type = { {"DPPC", 0},
				    {"DIPC", 1} };
}

// default destructor

parameters::~parameters() {

}


// accessor functions

double parameters::getplane_pair_energy(std::pair <std::string, std::string> species_pair) {
	return plane_pair_energy[species_pair];
}

double parameters::getinter_pair_energy(std::pair <std::string, std::string> species_pair) {
	return inter_pair_energy[species_pair];
}

int parameters::getoutput_type(std::string species) {
	return output_type[species];
}