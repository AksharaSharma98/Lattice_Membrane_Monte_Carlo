#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <iostream>

#include "parameters.h"


// default constructor

parameters::parameters() {
	pair_energy = { {std::make_pair("DPPC","DPPC"), -2.0},
					{std::make_pair("DIPC","DIPC"), -1.0},
					{std::make_pair("DIPC","DPPC"), -1.2} };
	output_type = { {"DPPC", 0},
				    {"DIPC", 1}};
}

// default destructor

parameters::~parameters() {

}


// accessor functions

double parameters::getpair_energy(std::pair <std::string, std::string> species_pair) {
	return pair_energy[species_pair];
}

int parameters::getoutput_type(std::string species) {
	return output_type[species];
}