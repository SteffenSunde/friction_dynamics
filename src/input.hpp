#ifndef INPUT_HPP
#define INPUT_HPP

//#include "data_structures.hpp"

#include "yaml-cpp/yaml.h"
#include <string>


struct SimulationData {

};


SimulationData parse_simulation_data(std::string const input_file);

void run_on_input_file(std::string const input_file);

#endif