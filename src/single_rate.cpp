#include "data_structures.hpp"
#include "single.hpp"
#include "integrators.hpp"
#include "output.hpp"

#include <iostream>
#include <random>
#include <memory>

#include "omp.h"

auto SingleRate::slope(Vec3 const& state) const -> Vec3 {
    Vec3 dxdt;
    double const relative_velocity = state[1] - v0;
    double const external_force = -k*state[0] - beta*k*state[1] - alpha*m*state[1];
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


void SingleRateSine::stiffness_damping(double const frequency, double const ratio) {
    /*
    Set Rayleigh coefficients assuming a purely stiffness-proportional damping.

    Frequency is in Hz.
    Ratio is ratio of damping to the critical amount for the given frequency.

    xi = 0.5 * (omega * beta)
    beta = 2*xi/omega = 2*xi/(2*pi*f) = xi /(pi*f)

    */
    alpha = 0.0;
    beta = ratio/(M_PI*frequency);
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
    double const damping_ratio = 0.05;
    std::string const meta = "_b005";
    std::string const output_file = "poincare_"+ std::to_string(frequency) +"_"+std::to_string(transient_periods)+"_"+std::to_string(num_intersections)+meta+".txt";

    SingleRateSine system(frequency);
    system.stiffness_damping(frequency, damping_ratio);

    std::cout << "Running SingleRateSine model with frequency " << std::to_string(frequency) << " Hz driver.\n"
              << "recording " << num_intersections << " intersections at whole periods. \n"
              << std::to_string(transient_periods) << " transient periods are discarded."
              << "Damping ratio : " << damping_ratio << ", initial_positions: " << num_initial_positions
              << "Perturbance: " << perturbance << "\n";
    
    auto intersections = calculate_poincare_sections_perturbed(
        system, time_step, num_initial_positions, num_intersections, transient_periods, perturbance);

    printf("Storing results in file %s \n", output_file.c_str());
    writeToCSVfile(output_file, intersections);
}


void single_rate_sine(std::string const& input_file)
{
    std::cout << "Not implemented yet!!\n";
}


auto calculate_poincare_sections_perturbed (
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


auto calculate_poincare_sections(
    SingleRateSine const& system, 
    double const& time_step, 
    int const& num_intersections, 
    int const& transient_periods
) -> Eigen::MatrixXd
{
    double const period_time = 1.0/system.f;
    int period_steps = (int)(period_time/time_step) + 1;
    double const initial_velocity = 2.0*M_PI*system.f*system.d;

    Eigen::MatrixXd intersections(num_intersections, 2);

    Vec3 state(0.0, initial_velocity, 0.0);

    for (int i=0; i < transient_periods; ++i) {
        double cycle_time = 0.0;
        for (int k=0; k < period_steps; ++k) {
            double const dt = std::min(time_step, period_time-cycle_time);
            state += step_rk4(system, state, dt);
            cycle_time += dt;
        }
    }

    for (int i=0; i < num_intersections; ++i) {
        double cycle_time = 0.0;
        for (int k=0; k < period_steps; ++k) {
            double const dt = std::min(time_step, period_time - cycle_time);
            state += step_rk4(system, state, dt);
            cycle_time += dt;
        } // TODO: Is it necessary to check direction of intersection?? No?
        intersections(i, 0) = state[0];
        intersections(i, 1) = state[1];
    }

    return intersections;
}


void single_chaos_finder(double const frequency_start, double const frequency_stop)
{
    /*
    Calculates Poincaré maps for a large series of random frequencies within the given range.
    Saves the resulting intersections only if variations in the displacement and velocity
    are within a certain criteria. Loops "forever"!

    TODO: 
    - Make parallel
    - Implement starting perturbed initial conditions

    */
    double const time_step = 1e-5;
    int const transient_periods = 1000;
    int const num_intersections = (int)1e4;
    int const num_initial_positions = 1;  // This is used for "perturbed maps"
    double const perturbance = 0.0;
    double const damping_ratio = 0.05;

    std::mt19937_64 generator(101);  // TODO Which seed?
    std::uniform_real_distribution<double> distribution(frequency_start, frequency_stop);

    printf("Searching for interesting Poincaré maps between frequencies %.1f and %.1f Hz...\n", frequency_start, frequency_stop);

    for (int i = 0; i < (int)1e6; ++i) {
        double frequency = distribution(generator);
        std::string const output_file = "poincare_"+ std::to_string(frequency)+".txt";

        SingleRateSine system(frequency);
        system.stiffness_damping(frequency, damping_ratio);
        
        auto intersections = calculate_poincare_sections(
            system, time_step, num_intersections, transient_periods);

        double ratio_displacement = (intersections.col(0).maxCoeff()-intersections.col(0).minCoeff())/(intersections.col(0).minCoeff()+1e-10);
        double ratio_velocity = (intersections.col(1).maxCoeff()-intersections.col(1).minCoeff())/(intersections.col(1).minCoeff()+1e-10);

        if (std::abs(ratio_displacement) > 0.1 && std::abs(ratio_velocity) > 0.1) {
            printf("Chaos may have been found for frequency %f! Storing results in file %s \n", frequency, output_file.c_str());
            writeToCSVfile(output_file, intersections);
        } else {
            printf("Nothing interesting found for frequency %f. Moving on...\n\n", frequency);
        }
    }    
}