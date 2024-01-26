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

	int get_num_species();

	int get_grid_size();

	std::string get_species(int leaflet, int i);

	int get_population(int leaflet, int i);

	std::vector<int> get_swap_sizes();

	std::vector<double> get_swap_weights();

	int get_output_type(std::string species);

protected:
	int num_species;
	int grid_size;

	// Vector of vectors containing the names of species in each leaflet
	std::vector<std::vector<std::string> > species;
	// Vector of vectors containing the corresponding population of each species in each leaflet
	std::vector<std::vector<int> > population;

	// Vector containing the different swap patch sizes to be used for MC swapping
	std::vector<int> swap_sizes;
	// Vector containing the different weights for the patch sizes
	std::vector<double> swap_weights;

	// contains integer-type config output values for keys of string-type species names
	std::map<std::string, int> output_type;
};

// file input for species names and corresponding populations
void read_system_parameters(std::vector<std::vector<std::string> >& species, std::vector<std::vector<int> >& population, std::map<std::string, int>& output_type, std::vector<int>& swap_sizes, std::vector<double>& swap_weights, int& grid_size);

// initializes the 2D lipid species vector with species names for each leaflet symmetrically
void initialize_lipid_species(std::vector<std::vector<std::string> >& species, std::istringstream& line);

// initializes the 2D lipid population vector with species populations for each leaflet symmetrically
void initialize_lipid_populations(std::vector<std::vector<int> >& population, std::istringstream& line);

// initializes the map containing species name and output integer pairs
void initialize_output_type(std::vector<std::vector<std::string> >& species, std::map<std::string, int>& output_type, std::istringstream& line);

// initializes the swap patch sizes vector
void initialize_swap_sizes(std::vector<int>& swap_sizes, std::istringstream& line);

// initializes the weights for the corresponding swap patch sizes
void initialize_swap_weights(std::vector<double>& swap_weights, std::istringstream& line);

#endif