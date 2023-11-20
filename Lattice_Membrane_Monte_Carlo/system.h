#ifndef System_h
#define System_h

#include <string>
#include <vector>
#include <map>


class System {
public:
	// default constructor, destructor
	System();
	~System();

	// accessor functions
	std::string get_species(int leaflet, int i);

	int get_population(int leaflet, int i);

protected:
	int num_species;
	int grid_size;

	// Vector of vectors containing the names of species in each leaflet
	std::vector<std::vector<std::string>> species;

	// Vector of vectors containing the corresponding population of each species in each leaflet
	std::vector<std::vector<int>> population;

};

// file input for species names and corresponding populations
void read_system_parameters(std::vector<std::vector<std::string>>& species, std::vector<std::vector<int>>& population);

#endif