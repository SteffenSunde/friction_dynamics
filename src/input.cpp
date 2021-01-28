#include "input.hpp"
#include "multi_rate_sine.hpp"

#include <iostream>
#include <vector>

SimulationData parse_simulation_data(std::string const input_file)
{
    return SimulationData();
}


void run_on_input_file(std::string const input_file)
{
    YAML::Node input = YAML::LoadFile(input_file);

    if(input["model"].as<std::string>() == "HertzRateSine") {
        int num_blocks = input["num_blocks"].as<int>();
        double peak_pressure = input["peak_pressure"].as<double>();
        int free_blocks = input["free_blocks"].as<int>();
        HertzRateSine model(num_blocks, peak_pressure, free_blocks);

        if(input["natural_frequencies"]) {
            std::vector<double> eigenvals = model.calc_natural_frequencies();
            double min_freq = 1e10;
            double max_freq = 0;
            for(auto i = 0; i < eigenvals.size(); ++i) {
                //std::cout << eigenvals[i] << "\n";
                min_freq = std::min(min_freq, eigenvals[i]);
                max_freq = std::max(max_freq, eigenvals[i]);
            }
            std::cout << "Highest natural frequency: " << max_freq << " Hz\n";
            std::cout << "Lowest natural frequency: " << min_freq << " Hz\n";
        }
    }
}