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
	membrane(int n, std::string* sp, int* pop, int sp_num);
	~membrane();

	// accessor functions
	int getsize();
	Grid& getgrid();
	lipid& getlipid(int i, int j);
	

	// modifier functions
	void swap(int* a, int* b);
	void tail_update(int* a, double s_new);

protected:

	int size;
	Grid grid;

};

#endif
