#include <string>
#include <assert.h>

#include "lipid2.h"


// default constructor

lipid2::lipid2(std::string sp, std::vector<double> tail_order, std::vector<int> coordinates, std::vector<int> indices)
{
	Species = sp;

	assert(coordinates.size() <=4 && "Too many 2D coordinates provided for lipid initialization\n"
									 "The program only supports upto 2-site lipids\n");
	for (int i = 0; i < coordinates.size(); i += 2) {
		std::vector<int> coord = { coordinates[i], coordinates[i+1] };
		Coords.push_back(coord);
	}

	for (int i = 0; i < indices.size(); i++) {
		Indices.push_back(indices[i]);
		Tail_order.push_back(tail_order[i]);
	}
}


// default destructor

lipid2::~lipid2() {

}


// accessor functions

std::string lipid2::species() {
	return Species;
}


std::vector<double> lipid2::tail_order() {
	return Tail_order;
}


std::vector<std::vector<int> > lipid2::coordinates() {
	return Coords;
}


std::vector<int> lipid2::indices() {
	return Indices;
}


// modifier functions

void lipid2::update_species(std::string species) {
	Species = species;
}


void lipid2::update_coords(std::vector<std::vector<int> > coords) {
	for (int i = 0; i < coords.size(); i++) {
		Coords[i] = coords[i];
	}
}


void lipid2::update_indices(std::vector<int> indices) {
	for (int i = 0; i < indices.size(); i++) {
		Indices[i] = indices[i];
	}
}


void lipid2::update_tail_order(std::vector<double> tail_order) {
	assert(tail_order.size() == Tail_order.size() && "Incorrect number of tail order values provided for update\n");
	for (int i = 0; i < tail_order.size(); i++) {
		Tail_order[i] = tail_order[i];
	}
}

