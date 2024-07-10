#ifndef moves_h
#define moves_h

#include <set>


double single_swap(membrane& current);

double multi_swap(membrane& current, int batch_size, std::vector<double>& acceptance);

double state_swap(membrane& current, int batch_size, std::vector<double>& acceptance);

double patch_swap(membrane& current, std::map<int, std::vector<double> >& patch_accept);

void site_picker_DPPC(membrane& leaflet, std::vector<std::vector<int> >& sites);

bool metropolis_accept(double delE);

void patch_sites(membrane& current, const std::vector<int>& a, int patch_size, std::vector<std::vector<int> >& patch, std::vector<int>& bounds);

void patch_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int& patch_size);

void patch_center_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int patch_size);

bool indentation_match(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int patch_size);

void indentation_sections(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, std::vector<std::vector<int> >& indent_sections);


#endif