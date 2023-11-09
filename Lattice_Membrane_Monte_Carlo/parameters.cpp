#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>   // temporary, for dev only

#include "parameters.h"


// parameters class
// 
// default constructor

parameters::parameters() {

	read_pair_energy("plane", plane_pair_energy);
	read_pair_energy("inter", inter_pair_energy);

	// map output to test if it was built correctly
	/*for (auto it = plane_pair_energy.cbegin(); it != plane_pair_energy.cend(); ++it)
	{
		printf("(%s,%s) -> %lf\n",it->first.first.c_str(), it->first.second.c_str(), it->second);
	}*/

	plane_entropy_const = 0.06;

	inter_entropy_const = 0.04;

	tail_order_bins = { {"DPPC",{-0.5,-0.35,-0.2,-0.05,0.1,0.25,0.4,0.55,0.7,0.85,1.0}},
						{"DIPC",{-0.5,-0.35,-0.2,-0.05,0.1,0.25,0.4,0.55,0.7,0.85,1.0}},
						{"CHOL",{0.4,0.45}} };

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


// IO functions for parameter class

void read_pair_energy(const std::string& domain, std::map<std::pair<std::string, std::string>, double>& pair_energy) {

	std::ifstream file("Parameters.txt");
	bool read_complete = false;
	std::string section_marker = "#" + domain;
	
	if (file.is_open()) {
		std::string line;

		while (std::getline(file, line) && !read_complete) {
			
			if (line.find(section_marker) != std::string::npos) {
				
				while (std::getline(file, line)) {
					if (line.empty())
					{
						read_complete = true;
						break;
					}
					else if (line.rfind("#", 0) != 0)
					{
						std::string lipid1, lipid2;
						double interaction_param;
						std::istringstream line_string(line);

						line_string >> lipid1 >> lipid2 >> interaction_param;

						// Note: Both possible combinations of a pair-key must be added for ease of use
						pair_energy.insert({ std::make_pair(lipid1,lipid2), interaction_param });
						pair_energy.insert({ std::make_pair(lipid2,lipid1), interaction_param });
					}
				}
			}
		}
		file.close();
	}
	/*for (auto it = pair_energy.cbegin(); it != pair_energy.cend(); ++it)
	{
		printf("(%s,%s) -> %lf\n", it->first.first.c_str(), it->first.second.c_str(), it->second);
	}
	printf("\n");*/
}