#ifndef parameters_h
#define parameters_h

#include <string>
#include <vector>
#include <map>


class parameters {
public:
	// default constructor, destructor

	parameters();
	~parameters();

	// accessor functions

	double getplane_pair_energy(std::pair <std::string, std::string> species_pair);

	double getinter_pair_energy(std::pair <std::string, std::string> species_pair);

	double getplane_entropy_const();

	double getinter_entropy_const();

	std::vector<double> tailorder_bins(std::string species);

	std::vector<double> tailorder_weights(std::string species);

	int getoutput_type(std::string species);
	
	// modifier functions

protected:
	// map that contains pairwise in_plane interaction energy values for a species-pair string key
	std::map<std::pair<std::string, std::string>, double> plane_pair_energy;
	// map that contains pairwise interleaflet interaction energy values for a species-pair string key
	std::map<std::pair<std::string, std::string>, double> inter_pair_energy;

	double plane_entropy_const;
	double inter_entropy_const;

	// map that contains tail-order distribution values for each species string key
	std::map<std::string, std::vector<double>> tail_order_bins;
	// map that contains tail-order distribution weights for each species string key
	std::map<std::string, std::vector<double>> tail_order_weights;

	// contains integer-type config output values for keys of string-type species names
	std::map<std::string, int> output_type;
};

// file input for pair energy parameters
void read_pair_energy(const std::string& domain, std::map<std::pair<std::string, std::string>, double>& pair_energy);

#endif
