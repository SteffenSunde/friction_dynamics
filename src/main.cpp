#include "single_rate_sine.hpp"
#include "multi_rate_sine.hpp"
#include "input.hpp"

#include <boost/program_options.hpp>

#include <stdio.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <random>

#include "omp.h"

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
            ("SingleChaosFinder", "Chaos finder.")
            ("SingleRateSine", "Single block with rate-dependent friction and Sine driver.")
            ("HertzRateSine", po::value<double>()->implicit_value(15.0), "N (100) blocks with rate-and state dependent friction and sine driver.")
            ("SingleRateHistory", po::value<double>()->implicit_value(15.0), "Calculate steady state for single DOF velocity-weakening friction at given frequency")
            ("SingleRatePoincare", po::value<double>()->implicit_value(15.0), "Calculate Poincare maps for single DOF velocity-weakening friction at given frequency")
            ("HertzEvolve", "Calculate evolving Hertzian fretting contact.")
            ("xi", po::value<double>()->implicit_value(0.05), "Damping ratio of lowest natural frequency")
            ("delta", po::value<double>()->implicit_value(1.0), "Friction slope")
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
                    
                // }
            }   
        } else if (vm.count("SingleChaosFinder")) {
            single_poincare_chaos_finder(10, 20.0);
        } else if (vm.count("SingleRateHistory")) {
            // po::options_description sub_opts("SingleRateHistory options");
            // sub_opts.add_options()("f", po::value<double>()->implicit_value(15.0));

            // po::store(po::command_line_parser(opts).options(sub_opts)).run(), vm);
            // std::cout << vm["delta"].as<std::string>() << "\n";
            // std::cout << vm["xi"].as<std::string>() << "\n";
            // double delta = vm["delta"].as<double>();
            // double xi = vm["xi"].as<double>();
            double delta = 5;
            //double xi = 0.05;
            //double frequency = vm["SingleRateHistory"].as<double>();
            printf("Calculating steady-state for a single velocity weakening oscillator\n");
            //for (double xi: {0.01, 0.05, 0.075, 0.1}) {
            //std::vector<long double> deltas = {0.5, 1.0, 2.0, 5.0};
            //std::vector<long double> xis = {0.001, 0.01, 0.05, 0.1};
            //std::vector<long double> xis = {0.0001};
            //#pragma omp parallel for
            //for (int i=0; i<xis.size(); ++i) {
                auto xi = vm["SingleRateHistory"].as<double>();
                double frequency = 5;
                single_rate_sine_history(frequency, xi, delta);
            //}
        } else if(vm.count("SingleRatePoincare")) {
            // po::options_description singleratepoincare_desc("SingleRatePoincare options");
            // singleratepoincare_desc.add_options()
            //     ("frequency", po::value<double>()->default_value(15.0), "Frequency of system");
            // std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            // opts.erase(opts.begin());  // New options include positional command which is not needed.
            // po::store(po::command_line_parser(opts).options(singleratepoincare_desc).run(), vm);  // Reparse
            // //std::cout << vm["frequency"].as<std::string>();
            double frequency = vm["frequency"].as<double>();
            printf("Calculating Poincaré maps for single block with f: %.2f\n", frequency);
            for (long double damping_ratio : {0.0001, 0.1}) {
                single_rate_sine_poincare(frequency, damping_ratio);
            }
        } else if (vm.count("HertzRateSine")) {
            //double const frequency = vm["HertzRateSine"].as<double>();

            std::vector<double> pressures = {15, 150};
            #pragma omp parallel for
            for(int i=0; i < pressures.size(); ++i) {
                double const frequency = 15;
                double const xi = 0.05;
                double const delta = 1.0;
                double const pressure = pressures[i];
                hertz_rate_sine_slip(frequency, delta, xi, pressure);
            }

            // std::for_each(std::execution::par_unseq, 
            // frequencies.begin(), 
            // frequencies.end(), [delta](double frequency) {
            //     hertz_rate_sine_slip(frequency, delta);
            // });
        } else if (vm.count("input")) {
            std::string const& file = vm["input"].as<std::string>();
            std::cout << "Running on input-file " << file << "\n";
        } else if (vm.count("HertzEvolve")) {

            std::vector<double> deltas = {0.5, 5.0};
            #pragma omp parallel for
            for(int i=0; i < deltas.size(); ++i) {
                double const frequency = 15.0;
                double const delta = deltas[i];
                double const damping_ratio = 0.05;
                double const evolve_rate = 0.01;
                
                hertz_evolve(frequency, delta, damping_ratio, evolve_rate);
            }
        }
    } catch (const po::error& ex) {
        std::cerr << ex.what() << "\n";
    }

    double elapsed = (double)std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - start_time).count();
    printf("Elapsed: %.3f seconds", elapsed);

    return 0;
} 