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
	double getpair_energy(std::pair <std::string, std::string> species_pair);
	int getoutput_type(std::string species);
	
	// modifier functions

protected:
	// contains pairwise interaction energy values for keys of 2 species type strings
	std::map<std::pair<std::string, std::string>, double> pair_energy;
	std::map<std::string, int> output_type;
};

#endif
