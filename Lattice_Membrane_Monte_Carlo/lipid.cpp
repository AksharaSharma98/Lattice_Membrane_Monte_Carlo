#include <string>

#include "lipid.h"


// default constructor

lipid::lipid(std::string sp, double tail_order, int* position)
{
	species = sp;
	s = tail_order;
	for (int i = 0; i < 2; i++) {
		coords[i] = position[i];
	}
}

// default destructor

lipid::~lipid() {

}


// accessor functions

std::string lipid::getspecies() {
	return species;
}


double lipid::gettail_order() {
	return s;
}


int* lipid::getposition() {
	return coords;
}


// modifier functions

void lipid::updatetail_order(double s_new) {
	s = s_new;
}