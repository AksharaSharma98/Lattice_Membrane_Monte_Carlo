#ifndef moves_h
#define moves_h

#include <set>


double single_swap(membrane& current);

double multi_swap(membrane& current, int batch_size);

double state_swap(membrane& current, int batch_size);

double patch_swap(membrane& current);

void lipid_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2);

void lipid_picker_species(membrane& leaflet, std::vector<std::vector<int>>& sites1, std::vector<std::vector<int>>& sites2);

void lipid_picker_DPPC(membrane& leaflet, std::vector<std::vector<int>>& sites);

bool metropolis_accept(double delE);

void unique_interaction_pairs(membrane& current, const std::vector<std::vector<int>>& sites1, const std::vector<std::vector<int>>& sites2, std::map<std::set<int>, std::vector<int>> &unique);

void unique_interaction_pairs(membrane& current, const std::vector<std::vector<int>>& sites, std::map<std::set<int>, std::vector<int>>& unique);

void patch_boundaries(membrane& current, const std::vector<int>& a, int patch_size, std::vector<std::vector<int>> &boundary, std::vector<int> &bounds);

void patch_center_picker(membrane& leaflet, std::vector<int>& c1, std::vector<int>& c2, int patch_size);

//double mc_move_entr(membrane& current, membrane& opposing, int* a, int* b);

#endif