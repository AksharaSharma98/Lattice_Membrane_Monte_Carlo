#ifndef lipid2_h
#define lipid2_h

#include <string>
#include <vector>


class lipid2 {
public:
	// default constructor, destructor
	lipid2(std::string sp, std::vector<double> tail_order, std::vector<int> coordinates, std::vector<int> indices);
	~lipid2();

	// accessor functions
	std::string species();
	std::vector<double> tail_order();
	std::vector<std::vector<int> > coordinates();
	std::vector<int> indices();

	// modifier functions
	void update_species(std::string species);
	void update_coords(std::vector<std::vector<int> > coords);
	void update_indices(std::vector<int> indices);
	void update_tail_order(std::vector<double> tail_order);

protected:
	std::string Species;
	std::vector<double> Tail_order;
	std::vector<std::vector<int> > Coords;
	std::vector<int> Indices;
};

#endif