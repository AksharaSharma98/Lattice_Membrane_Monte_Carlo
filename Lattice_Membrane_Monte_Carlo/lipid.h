#ifndef lipid_h
#define lipid_h

#include <string>


class lipid {
public:
	// default constructor, destructor
	lipid(std::string sp, double tail_order, int* position);
	~lipid();

	// accessor functions

	std::string getspecies();
	double gettail_order();
	int* getposition();

	// modifier functions

	void update_tailorder(double s_new);

protected:
	std::string species;
	double s;
	int coords[2];

};

#endif
