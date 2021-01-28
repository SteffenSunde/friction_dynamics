#include "data_structures.hpp"
#include "single_rate_sine.hpp"
#include "integrators.hpp"
#include "output.hpp"

#include <iostream>
#include <random>
#include <memory>

#include "omp.h"

auto SingleRateSine::slope(Vec3 const& state) const -> Vec3 {
    Vec3 dxdt;
    double const belt_velocity = 2.0*M_PI*f*d*std::cos(2.0*M_PI*f*state[2]);
    double const belt_acceleration = -std::pow(2.0*M_PI*f, 2.0)*d*std::sin(2.0*M_PI*f*state[2]);
    double const relative_velocity = state[1] - belt_velocity;
    double const external_force = -k*state[0] - c*state[1];
    if (std::abs(relative_velocity) < eps) {
        double const friction_limit = cof_static*p;
        if (std::abs(external_force + m*belt_acceleration) <= friction_limit) {
            dxdt[0] = belt_velocity;
            dxdt[1] = belt_acceleration;
        } else {
            dxdt[0] = state[1];
            dxdt[1] = 1.0/m * (external_force - friction_limit*sgn(external_force));
        }
    } else {
        dxdt[0] = state[1];
        dxdt[1] = 1.0/m * (external_force - friction(relative_velocity)*cof_kinetic*p*sgn(relative_velocity));
    }

    dxdt[2] = 1.0;
    return dxdt;
}

void single_rate_sine()
{
    std::cout << "Not implemented yet!\n";
}

void single_rate_sine(
    double frequency, 
    int transient_periods, 
    int num_intersections) 
{
    double const time_step = 1e-5;
    //int const transient_periods = 200;
    int const num_initial_positions = 100;
    //int const num_intersections = 100;
    //double const frequency = 19.5;
    double const perturbance = 1e-14;
    std::string const meta = "_b005";
    std::string const output_file = "poincare_"+ std::to_string(frequency) +"_"+std::to_string(transient_periods)+"_"+std::to_string(num_intersections)+meta+".txt";

    SingleRateSine system(frequency);
    system.c = 0.005*system.k;

    std::cout << "Running SingleRateSine model with frequency " << std::to_string(frequency) << " Hz driver.\n"
              << "recording " << num_intersections << " intersections at whole periods. \n"
              << std::to_string(transient_periods) << " transient periods are discarded."
              << "Damping: " << system.c << ", initial_positions: " << num_initial_positions
              << "Perturbance: " << perturbance << "\n";
    
    auto intersections = calculate_poincare_sections(
        system, time_step, num_initial_positions, num_intersections, transient_periods, perturbance);

    printf("Storing results in file %s \n", output_file.c_str());
    writeToCSVfile(output_file, intersections);
}


void single_rate_sine(std::string const& input_file)
{
    std::cout << "Not implemented yet!!\n";
}


auto calculate_poincare_sections(
    SingleRateSine const& system, 
    double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    double const& perturbance) -> Eigen::MatrixXd
{
    double const period_time = 1.0/system.f;
    int period_steps = (int)(period_time/time_step) + 1;
    int total_points = num_initial_positions * num_intersections;
    double const initial_velocity = 2.0*M_PI*system.f*system.d;

    std::mt19937_64 generator(100);  // TODO Which seed?
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    Eigen::MatrixXd intersections(total_points, 2);

    int const num_threads = omp_get_max_threads();

    #pragma omp parallel shared(intersections)
    {   // TODO: Distribute better? (e.g. if one in last thread)
        int const taks_per_thread = static_cast<int>(num_initial_positions/omp_get_num_threads());
        int const thread_id = omp_get_thread_num();
        int const start = thread_id * taks_per_thread;

        int stop;
        if(thread_id != num_threads-1) {
            stop = start + taks_per_thread;
        } else {
            stop = num_initial_positions;
        }  

        for (int i=start; i<stop; ++i) {
            double const perturb = dis(generator) * perturbance;
            Vec3 state(perturb, initial_velocity, 0.0);

            for (int j=0; j < transient_periods; ++j) {
                double cycle_time = 0.0;
                for (int k=0; k < period_steps; ++k) {
                    double const dt = std::min(time_step, period_time-cycle_time);
                    state += step_rk4(system, state, dt);
                    cycle_time += dt;
                }
            }

            for (int j=0; j < num_intersections; ++j) {
                double cycle_time = 0.0;
                for (int k=0; k < period_steps; ++k) {
                    double const dt = std::min(time_step, period_time - cycle_time);
                    state += step_rk4(system, state, dt);
                    cycle_time += dt;
                } // TODO: Is it necessary to check direction of intersection?? No?
                intersections(i*num_intersections + j, 0) = state[0];
                intersections(i*num_intersections + j, 1) = state[1];
            }
        }
    }

    return intersections;
}