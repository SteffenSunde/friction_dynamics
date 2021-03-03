#include "single_rate_sine.hpp"
#include "multi_rate_sine.hpp"
#include "input.hpp"

#include <boost/program_options.hpp>

#include <stdio.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <random>

namespace po = boost::program_options;
using Clock = std::chrono::system_clock; 

int main(int argc, const char* argv[]) 
{
    auto start_time = Clock::now();

    try {
        po::options_description global{"Allowed options"};
        global.add_options()
            ("help,h", "Display help message")
            ("input,i", po::value<std::string>(), "Run simulator on given YAML input file")
            ("integrators", "List available integrators")
            ("todo", "Display most important TODOS.")
            ("SingleRateSine", "Single block with rate-dependent friction and Sine driver.")
            ("HertzRateSine", "N (100) blocks with rate-and state dependent friction and sine driver.")
            ("SingleRateHistory", po::value<double>()->implicit_value(15.0), "Calculate steady state for single DOF velocity-weakening friction at given frequency")
            ("SingleRatePoincare", po::value<double>()->implicit_value(15.0), "Calculate Poincare maps for single DOF velocity-weakening friction at given frequency")
            //("SingleRatePoincare", po::value<std::vector<std::string> >()->default_value({}), "Poincare map for single oscillator")
        ;

        // po::positional_options_description pos;
        // pos.add("SingleRatePoincare", 1)
        //     .add("subargs", -1);

        po::variables_map vm;
        po::parsed_options parsed = po::command_line_parser(argc, argv)
            .options(global)
            .allow_unregistered()
            .run();
        po::store(parsed, vm);
        po::notify(vm);   
        std::vector<std::string> unrecognized = po::collect_unrecognized(parsed.options, po::include_positional);

        // for (auto el: vm) {
        //     std::cout << el.first << "; ";
        // }

        if (vm.empty()) {
            std::cout << "Spring-block simulator. See --help for more info.\n";
        } else if (vm.count("help")) {
            std::cout << global << "\n";
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
                // std::mt19937_64 generator(100);  // TODO Which seed?
                // std::uniform_real_distribution<long double> dis(0.0, 1.0);
                // #pragma omp parallel for
                // for(int x=0; x < (int)1e6; x++){
                    single_poincare_chaos_finder(2, 30.0);
                // }
            }
        } else if (vm.count("SingleRateHistory")) {
            // po::options_description sub_opts("SingleRateHistory options");
            // sub_opts.add_options()("f", po::value<double>()->implicit_value(15.0));

            // po::store(po::command_line_parser(opts).options(sub_opts)).run(), vm);
            double frequency = vm["SingleRateHistory"].as<double>();
            printf("Calculating steady-state for a single velocity weakening oscillator (Frequency: %.2f)\n", frequency);
            single_rate_sine_history(frequency);
        } else if(vm.count("SingleRatePoincare")) {
            // po::options_description singleratepoincare_desc("SingleRatePoincare options");
            // singleratepoincare_desc.add_options()
            //     ("frequency", po::value<double>()->default_value(15.0), "Frequency of system");
            // std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            // opts.erase(opts.begin());  // New options include positional command which is not needed.
            // po::store(po::command_line_parser(opts).options(singleratepoincare_desc).run(), vm);  // Reparse
            // //std::cout << vm["frequency"].as<std::string>();
            // double frequency = vm["frequency"].as<double>();
            //printf("Calculating Poincaré map for single block with f: %.2f\n", frequency);
            double frequency = vm["SingleRatePoincare"].as<double>();
            single_rate_sine_poincare(frequency);
        } else if (vm.count("HertzRateSine")) {
            double const frequency = vm["HertzRateSine"].as<double>();
            hertz_rate_sine_slip(frequency);
            //calculate_multi_poincare_sections();
        } else if (vm.count("input")) {
            std::string const& file = vm["input"].as<std::string>();
            std::cout << "Running on input-file " << file << "\n";
            //run_on_input_file(file);
        }
    } catch (const po::error& ex) {
        std::cerr << ex.what() << "\n";
    }

    double elapsed = (double)std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - start_time).count();
    printf("Elapsed: %.3f seconds", elapsed);

    return 0;
} 