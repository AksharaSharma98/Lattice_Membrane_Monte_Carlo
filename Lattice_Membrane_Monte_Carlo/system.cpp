#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>   // temporary, for dev only
#include <assert.h>

#include "system.h"


// system class
// 
// default constructor

System::System() {

	read_system_parameters(species, population, output_type, swap_sizes, swap_weights, grid_size);

	num_species = species[0].size();

	// testing
	/*for (int i = 0; i < species.size(); i++) {
		for (int j = 0; j < species[i].size(); j++) {
			std::cout << species[i][j] << " ";
		}
		printf("\n");
	}
	for (int i = 0; i < swap_sizes.size(); i++) {
		std::cout << swap_sizes[i] << " ";
	}
	printf("\n");
	for (int i = 0; i < swap_weights.size(); i++) {
		std::cout << swap_weights[i] << " ";
	}
	printf("\n");*/

}

System::~System() {

}

// accessor functions

int System::get_num_species() {
	return num_species;
}

int System::get_grid_size() {
	return grid_size;
}

std::string System::get_species(int leaflet, int i) {
	return species[leaflet][i];
}

int System::get_population(int leaflet, int i) {
	return population[leaflet][i];
}

std::vector<int> System::get_swap_sizes() {
	return swap_sizes;
}

std::vector<double> System::get_swap_weights() {
	return swap_weights;
}

int System::get_output_type(std::string species) {
	return output_type[species];
}

// IO functions for system class

void read_system_parameters(std::vector<std::vector<std::string>>& species, std::vector<std::vector<int>>& population, std::map<std::string, int>& output_type, std::vector<int>& swap_sizes, std::vector<double>& swap_weights, int& grid_size) {
	std::ifstream file("Parameters.txt");
	std::string section_marker = "#System";

	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {

			if (line.find(section_marker) != std::string::npos) {
				int count = 0;
				while (std::getline(file, line) && count <=5) {
					if (line.empty())
					{
						break;
					}
					else if (line.rfind("#", 0) != 0)
					{
						count++;
						std::istringstream line_string(line);
						
						if (count == 1) {
							line_string >> grid_size;
						}
						else if (count == 2){
							initialize_lipid_species(species, line_string);
						}
						else if(count == 3) {
							initialize_lipid_populations(population, line_string);
						}
						else if (count == 4) {
							initialize_output_type(species, output_type, line_string);
						}
						else if (count == 5) {
							initialize_swap_sizes(swap_sizes, line_string);
						}
						else if (count == 6) {
							initialize_swap_weights(swap_weights, line_string);
						}
					}
				}
			}
		}
		file.close();
	}
}


void initialize_lipid_species(std::vector<std::vector<std::string>>& species, std::istringstream& line) {

	std::string lipid;
	std::vector<std::string> lipids;
	while (line >> lipid) {
		lipids.push_back(lipid);
	}
	species.push_back(lipids);
	species.push_back(lipids);
}


void initialize_lipid_populations(std::vector<std::vector<int>>& population, std::istringstream& line) {
	int pop;
	std::vector<int> pops;
	while (line >> pop) {
		pops.push_back(pop);
	}
	population.push_back(pops);
	population.push_back(pops);
}


void initialize_output_type(std::vector<std::vector<std::string>>& species, std::map<std::string, int>& output_type, std::istringstream& line) {
	int otype;
	int count = 0;
	while (line >> otype) {
		output_type.insert({species[0][count], otype});
		count++;
	}
}


void initialize_swap_sizes(std::vector<int>& swap_sizes, std::istringstream& line) {
	int size;
	while (line >> size) {
		swap_sizes.push_back(size);
	}
}


void initialize_swap_weights(std::vector<double>& swap_weights, std::istringstream& line) {
	int weight;
	while (line >> weight) {
		swap_weights.push_back(weight);
	}
	// normalize weights
	double total = 0.0;
	for (int i = 0; i < swap_weights.size(); i++) {
		total += swap_weights[i];
	}

	assert(total != 0.0 && "No non-zero weights specified for any swap patch size\n");
	double norm = 1.0 / total;
	for (int i = 0; i < swap_weights.size(); i++) {
		swap_weights[i] = swap_weights[i]*norm;
	}
}

