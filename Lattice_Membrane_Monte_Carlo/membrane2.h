#ifndef membrane2_h
#define membrane2_h

#include <string>
#include <vector>

#include "lipid2.h"


class membrane2 {
	// typedef for the 2D vector of lipids
	typedef std::vector<lipid2> lipid_list;
	typedef std::vector<lipid_list> lipid_lists;

public:
	// default constructor, destructor
	membrane2(int leaflet);
	~membrane2();

	// accessor functions
	int grid_size();
	int leaflet_index();
	lipid_list& species_list(int species);
	std::vector<std::vector<int> > get_lattice();

	// modifier functions
	void swap(std::vector<std::vector<int> > list_indices);
	void swap_DPPC_state(std::vector<int> list_indices);

protected:
	int Grid_size;
	int Leaflet_index;
	lipid_lists Lipid_lists;
	std::vector<std::vector<int> > Grid;

};


#endif