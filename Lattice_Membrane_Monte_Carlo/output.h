#ifndef output_h
#define output_h


void write_config_int(FILE* config_file, membrane& upper, membrane& lower);

void write_config_bin(std::ofstream& config_bin_file, membrane& upper, membrane& lower);

void write_tailconfig(FILE* config_file, membrane& upper, membrane& lower);

void write_config_species(FILE* config_file, membrane& upper, membrane& lower);

void write_energy(FILE* energy_file, double energy);

void write_header(FILE* energy_file);

void print_leaflet_species(membrane& leaflet);

#endif
