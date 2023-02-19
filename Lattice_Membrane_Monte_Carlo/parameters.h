#ifndef parameters_h
#define parameters_h

#include <string>
#include <map>


class parameters {
public:
	// default constructor, destructor
	parameters();
	~parameters();

	// accessor functions
	double getplane_pair_energy(std::pair <std::string, std::string> species_pair);
	double getinter_pair_energy(std::pair <std::string, std::string> species_pair);
	int getoutput_type(std::string species);
	
	// modifier functions

protected:
	// contains pairwise in_plane interaction energy values for keys of 2 species name strings
	std::map<std::pair<std::string, std::string>, double> plane_pair_energy;
	// contains pairwise interleaflet interaction energy values for keys of 2 species name strings
	std::map<std::pair<std::string, std::string>, double> inter_pair_energy;
	// contains integer-type config output values for keys of string-type species names
	std::map<std::string, int> output_type;
};

#endif
