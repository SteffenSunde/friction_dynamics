#include "single_rate_sine.hpp"
#include "multi_rate_sine.hpp"
#include "input.hpp"

#include <boost/program_options.hpp>

#include <stdio.h>
#include <iostream>
#include <vector>
#include <chrono>

namespace po = boost::program_options;
using Clock = std::chrono::system_clock; 

int main(int argc, const char* argv[]) 
{
    auto start_time = Clock::now();

    try {
        po::options_description desc{"Allowed options"};
        desc.add_options()
            ("help,h", "Display help message")
            ("input,i", po::value<std::string>(), "Run simulator on given YAML input file")
            ("integrators", "List available integrators")
            ("todo", "Display most important TODOS.")
            ("SingleRateSine", "Single block with rate-dependent friction and Sine driver.")
            ("HertzRateSine", "N (100) blocks with rate-and state dependent friction and sine driver.")
        ;

        po::variables_map vm;
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
        po::store(parsed, vm);
        po::notify(vm);   
        std::vector<std::string> unrecognized = po::collect_unrecognized(parsed.options, po::include_positional);

        // for (auto el: vm) {
        //     std::cout << el.first << "; ";
        // }

        if (vm.empty()) {
            std::cout << "Spring-block simulator. See --help for more info.\n";
        } else if (vm.count("help")) {
            std::cout << desc << "\n";
        } else if (unrecognized.size() > 0) {
            std::cout << "Error: Did not understand the following parameters\n";
            for (auto c: unrecognized) std::cout << c << " ";
            std::cout << ". \n Please see the --help command";
        } else if( vm.count("integrators")) {
            std::cout << "None\n";
        } else if (vm.count("SingleRateSine")) {
            std::cout << "Running single block with Rate-dependent friction and sine-driver.\n";
            if (!vm["input"].empty()) { // Run on input file
                single_rate_sine(vm["input"].as<std::string>());
            //} else if (!vm["input"].empty()) { // Run with standard parameters
            //     single_rate_sine());
            } else {
                // double frequency = vm["frequency"].as<double>();
                // int transient_periods = vm["transients"].as<int>();
                // int num_intersections = vm["intersections"].as<int>();
                single_poincare_chaos_finder(20.0, 110.0);
            }
        } else if (vm.count("HertzRateSine")) {
            std::cout << "HertzRateSine\n";
            //std::cout << "Running HertzRateSine model with standard parameters\n";
            //hertz_rate_sine();
            // for (int freq=10; freq < 21; ++freq) {
            //     hertz_rate_sine_shear((double)freq);
            // }
            
            hertz_rate_sine_slip();
            //calculate_multi_poincare_sections();
        } else if (vm.count("input")) {
            std::string const& file = vm["input"].as<std::string>();
            std::cout << "Running on input-file " << file << "\n";
            run_on_input_file(file);
        }
    } catch (const po::error& ex) {
        std::cerr << ex.what() << "\n";
    }

    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - start_time).count();
    printf("Elapsed: %.3f seconds", elapsed);

    return 0;
} 