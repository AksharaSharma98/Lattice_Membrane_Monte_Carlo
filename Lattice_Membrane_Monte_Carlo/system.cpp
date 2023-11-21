#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>   // temporary, for dev only

#include "system.h"


// system class
// 
// default constructor

System::System() {

	read_system_parameters(species, population);

	// testing
	/*for (int i = 0; i < population.size(); i++) {
		for (int j = 0; j < population[i].size(); j++) {
			std::cout << population[i][j] << " ";
		}
		printf("\n");
	}*/

	num_species = species[0].size();
	grid_size = 100;

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

// IO functions for system class

void read_system_parameters(std::vector<std::vector<std::string>>& species, std::vector<std::vector<int>>& population) {
	std::ifstream file("Parameters.txt");
	std::string section_marker = "#system";

	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {

			if (line.find(section_marker) != std::string::npos) {
				int count = 0;
				while (std::getline(file, line) && count <=2) {
					if (line.empty())
					{
						break;
					}
					else if (line.rfind("#", 0) != 0)
					{
						count++;
						std::istringstream line_string(line);

						if (count==1){
							std::string lipid;
							std::vector<std::string> lipids;
							while (line_string >> lipid) { 
								lipids.push_back(lipid);
							}
							species.push_back(lipids);
							species.push_back(lipids);
						}

						else if(count == 2) {
							int pop;
							std::vector<int> pops;
							while (line_string >> pop) {
								pops.push_back(pop);
							}
							population.push_back(pops);
							population.push_back(pops);
						}
					}
				}
			}
		}
		file.close();
	}
}