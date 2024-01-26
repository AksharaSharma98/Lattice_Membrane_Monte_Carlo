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

Parameters::Parameters() {

	read_pair_energy("Plane", plane_pair_energy);
	read_pair_energy("Inter", inter_pair_energy);

	// map output to test if it was built correctly
	/*for (auto it = plane_pair_energy.cbegin(); it != plane_pair_energy.cend(); ++it)
	{
		printf("(%s,%s) -> %lf\n",it->first.first.c_str(), it->first.second.c_str(), it->second);
	}*/

	plane_entropy_const = 3.9; //0.06;

	inter_entropy_const = 0.04;

	tail_order_bins = { {"DPPC",{-0.5,-0.35,-0.2,-0.05,0.1,0.25,0.4,0.55,0.7,0.85,1.0}},
						{"DIPC",{-0.5,-0.35,-0.2,-0.05,0.1,0.25,0.4,0.55,0.7,0.85,1.0}},
						{"CHOL",{0.4,0.45}} };

	tail_order_weights = { {"DPPC",{0.0,0.0,0.0,0.0,0.17,0.29,0.36,0.15,0.03,0.0}},
						   {"DIPC",{0.0,0.0,0.07,0.25,0.41,0.20,0.07,0.0,0.0,0.0}},
						   {"CHOL",{1.0}} };
}

// default destructor

Parameters::~Parameters() {

}


// accessor functions

double Parameters::getplane_pair_energy(std::pair <std::string, std::string> species_pair) {
	return plane_pair_energy[species_pair];
}

double Parameters::getinter_pair_energy(std::pair <std::string, std::string> species_pair) {
	return inter_pair_energy[species_pair];
}

double Parameters::getplane_entropy_const() {
	return plane_entropy_const;
}

double Parameters::getinter_entropy_const() {
	return inter_entropy_const;
}

std::vector<double> Parameters::tailorder_bins(std::string species) {
	return tail_order_bins[species];
}

std::vector<double> Parameters::tailorder_weights(std::string species) {
	return tail_order_weights[species];
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
}