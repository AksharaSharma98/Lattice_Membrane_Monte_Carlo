#include <string>
#include <vector>
#include <map>

#include "parameters.h"


// default constructor

parameters::parameters() {

	// Note: Both possible combinations of a pair-key must be added unless they are the same
	plane_pair_energy = { {std::make_pair("DPPC","DPPC"), -2.0},
						  {std::make_pair("DIPC","DIPC"), -0.8},
						  {std::make_pair("CHOL","CHOL"), -0.5},
						  {std::make_pair("DPPC","DIPC"), -1.2},
						  {std::make_pair("DIPC","DPPC"), -1.2},
						  {std::make_pair("DPPC","CHOL"), -1.0},
						  {std::make_pair("CHOL","DPPC"), -1.0},
						  {std::make_pair("CHOL","DIPC"), -1.0},
						  {std::make_pair("DIPC","CHOL"), -1.0} };

	inter_pair_energy = { {std::make_pair("DPPC","DPPC"), -1.0},
						  {std::make_pair("DIPC","DIPC"), -2.0},
						  {std::make_pair("CHOL","CHOL"), 0.0},
						  {std::make_pair("DPPC","DIPC"), -1.5},
						  {std::make_pair("DIPC","DPPC"), -1.5},
						  {std::make_pair("DPPC","CHOL"), 0.0},
						  {std::make_pair("CHOL","DPPC"), 0.0},
						  {std::make_pair("CHOL","DIPC"), 0.0},
						  {std::make_pair("DIPC","CHOL"), 0.0} };

	plane_entropy_const = 1.0;

	inter_entropy_const = 5.0;

	tail_order_bins = { {"DPPC",{-0.5,-0.35,-0.2,-0.05,0.1,0.25,0.4,0.55,0.7,0.85,1.0}},
						{"DIPC",{-0.5,-0.35,-0.2,-0.05,0.1,0.25,0.4,0.55,0.7,0.85,1.0}},
						{"CHOL",{0.25,0.3}} };

	tail_order_weights = { {"DPPC",{0.0,0.0,0.0,0.0,0.17,0.29,0.36,0.15,0.03,0.0}},
						   {"DIPC",{0.0,0.0,0.07,0.25,0.41,0.20,0.07,0.0,0.0,0.0}},
						   {"CHOL",{1.0}} };

	output_type = { {"DPPC", 0},
					{"CHOL", 1},
				    {"DIPC", 2} };
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

double parameters::getplane_entropy_const() {
	return plane_entropy_const;
}

double parameters::getinter_entropy_const() {
	return inter_entropy_const;
}

std::vector<double> parameters::tailorder_bins(std::string species) {
	return tail_order_bins[species];
}

std::vector<double> parameters::tailorder_weights(std::string species) {
	return tail_order_weights[species];
}

int parameters::getoutput_type(std::string species) {
	return output_type[species];
}