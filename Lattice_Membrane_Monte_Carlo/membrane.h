#ifndef membrane_h
#define membrane_h

#include <string>
#include <vector>

#include "lipid.h"


class membrane {
	// typedef for the 2D vector of lipids
	typedef std::vector<lipid> Array;
	typedef std::vector<Array> Grid;

public:
	// default constructor, destructor
	membrane(int leaflet);
	~membrane();

	// accessor functions
	int getsize();
	Grid& getgrid();
	lipid& getlipid(int i, int j);
	

	// modifier functions
	void swap(std::vector<int> &a, std::vector<int> &b);
	void swap_DPPC_state(std::vector<int>& a);
	void patch_swap(std::vector<int>& bounds1, std::vector<int>& bounds2, int patch_size);
	void tail_update(int* a, double s_new);

protected:

	int size;
	Grid grid;

};

#endif
