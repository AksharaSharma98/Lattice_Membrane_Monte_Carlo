#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

#include "lipid.h"
#include "membrane.h"
#include "initialize.h"
#include "Lattice_Membrane_Monte_Carlo.h"


std::pair<membrane, membrane> initialize(int& steps, int& energy_output_freq, int& config_output_freq, int& restart_output_freq) {

    // read simulation parameters
	std::ifstream file("Parameters.txt");
	std::string section_marker = "#Simulation";

	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {

			if (line.find(section_marker) != std::string::npos) {
				int count = 0;
				while (std::getline(file, line) && count <= 3) {
					if (line.empty())
					{
						break;
					}
					else if (line.rfind("#", 0) != 0)
					{
						count++;
						std::istringstream line_string(line);
						std::string keyword;
						int temp;
						line_string >> keyword >> temp;

						if (keyword == "steps:") {
							steps = temp;
						}
						else if (keyword == "write_energy:") {
							energy_output_freq = temp;
						}
						else if (keyword == "write_config:") {
							config_output_freq = temp;
						}
						else if (keyword == "write_restart:") {
							restart_output_freq = temp;
						}
					}
				}
			}
		}
	}

	membrane upper(0);
	membrane lower(1);
	return std::make_pair(upper, lower);
}