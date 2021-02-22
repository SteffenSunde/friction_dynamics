#include "data_structures.hpp"
#include "single_rate_sine.hpp"
#include "integrators.hpp"
#include "output.hpp"

#include <iostream>
#include <random>
#include <memory>

#include "omp.h"

auto SingleRateSine::velocity_at_time(long double const& time) const -> long double 
{
    return 2.0*M_PI*f*d*std::cos(2.0*M_PI*f*time);
}

auto SingleRateSine::initial_state() const -> Vec3 
{
    return Vec3{0.0, velocity_at_time(0), 0.0};
}

auto SingleRateSine::slope(Vec3 const& state) const -> Vec3 {
    Vec3 dxdt;
    long double const belt_velocity = 2.0*M_PI*f*d*std::cos(2.0*M_PI*f*state[2]);
    long double const belt_acceleration = -std::pow(2.0*M_PI*f, 2.0)*d*std::sin(2.0*M_PI*f*state[2]);
    long double const relative_velocity = state[1] - belt_velocity;
    long double const external_force = -k*state[0] - beta*k*state[1] - alpha*m*state[1];
    if (std::abs(relative_velocity) < eps) {
        long double const friction_limit = cof_static*p;
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

    dxdt[2] = (long double)1.0;
    return dxdt;
}


void SingleRateSine::stiffness_damping(long double const frequency, long double const ratio) {
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
    long double frequency, 
    int transient_periods, 
    int num_intersections) 
{
    long double const time_step = 1e-5;
    //int const transient_periods = 200;
    int const num_initial_positions = 100;
    //int const num_intersections = 100;
    //double const frequency = 19.5;
    long double const perturbance = 1e-14;
    long double const damping_ratio = 0.05;
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
    long double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    long double const& perturbance) -> Eigen::MatrixXd
{
    long double const period_time = 1.0/system.f;
    int period_steps = (int)(period_time/time_step) + 1;
    int total_points = num_initial_positions * num_intersections;
    long double const initial_velocity = 2.0*M_PI*system.f*system.d;

    std::mt19937_64 generator(100);  // TODO Which seed?
    std::uniform_real_distribution<long double> dis(0.0, 1.0);

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
    long double const& time_step, 
    int const& num_intersections, 
    int const& transient_periods
) -> Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>
{
    /*

    TODO: Investigate time drift. Increases with increasing frequency.
    */

    bool half_period = false;  // Flag to set if half-periods are to be intersected. TODO: Check 

    long double const period_time = (long double)1.0/system.f;
    int period_steps = (int)(period_time/time_step)+1;
    long double const initial_velocity = system.velocity_at_time(0.0);

    Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> intersections(num_intersections, 3);

    Vec3 state(0.0, initial_velocity, 0.0);

    long double acc_drift = 0.0;
    long double time = 0.0;
    for (int i=0; i < transient_periods; ++i) {
        long double cycle_time = 0.0;
        for (int k=0; k < period_steps; ++k) {
            long double const dt = std::min(time_step, period_time-cycle_time);
            cycle_time += dt;
            // time += dt;
            //auto dstate =   // WTF why is time in system drifting but not diff in cycle_time??
            //if (std::abs(dstate[2]-dt) > 1e-9) throw std::string("ERROR: Drift");
            state += step_rk4(system, state, dt);  // Maybe it is fmod which is wrong. See below.
            // state[2] = cycle_time;
            // time += dt;
            // state[2] = time;
        }
        acc_drift += period_time - cycle_time;
    }
    
    printf("Transient drift: %e\n", state[2]/period_time-(long double)transient_periods);

    if (false) {
        for (int i=0; i < (int)(transient_periods/2); ++i) {
            long double cycle_time = 0.0;
            for (int k=0; k < (int)(period_steps/2); ++k) {
                long double const dt = std::min(time_step, period_time-cycle_time);
                cycle_time += dt;
                state += step_rk4(system, state, dt);
                // state[2] = cycle_time;
            }
        }
    }

    for (int i=0; i < num_intersections; ++i) {
        long double cycle_time = 0.0;
        long double dt = time_step;
        for (int k=0; k < period_steps; ++k) {
            dt = std::min(time_step, period_time - cycle_time);
            state += step_rk4(system, state, dt);
            cycle_time += dt;
            // time += dt;
            // state[2] = time;
        }
        acc_drift += period_time - period_time;
        intersections.row(i) = state;
        // intersections(i, 0) = state[0];
        // intersections(i, 1) = state[1];
    }
    printf("period_time - cycle_time: %e\n", acc_drift);
    printf("Calculated %d intersections with frequency %.2f and time step %e. Accumulated time drift: %e\n",
        num_intersections, system.f, time_step, state[2]/period_time - (long double)(num_intersections+transient_periods));  //std::fmod(state[2], period_time)

    return intersections;
}


void single_poincare_chaos_finder(long double const frequency_start, long double const frequency_stop)
{
    /*
    Calculates Poincaré maps for a large series of random frequencies within the given range.
    Saves the resulting intersections only if variations in the displacement and velocity
    are within a certain criteria. Loops "forever"!

    TODO: 
    - Make parallel
    - Implement starting perturbed initial conditions

    */
    long double const time_step = 1e-6;
    int const transient_periods = 1000;
    int const num_intersections = (int)1e5;
    int const num_initial_positions = 1;  // This is used for "perturbed maps"
    long double const perturbance = 0.0;
    long double const damping_ratio = 0.05;

    std::mt19937_64 generator(101);  // TODO Which seed?
    std::uniform_real_distribution<long double> distribution(frequency_start, frequency_stop);

    printf("Searching for interesting Poincaré maps between frequencies %.1f and %.1f Hz...\n", frequency_start, frequency_stop);

    for (int i = 0; i < (int)1e6; ++i) {
        long double frequency = distribution(generator);
        std::string const output_file = "poincare_"+ std::to_string(frequency)+".txt";

        SingleRateSine system(frequency);
        system.stiffness_damping(frequency, damping_ratio);
        
        auto intersections = calculate_poincare_sections(
            system, time_step, num_intersections, transient_periods);

        long double ratio_displacement = (intersections.col(0).maxCoeff()-intersections.col(0).minCoeff())/(intersections.col(0).minCoeff()+1e-10);
        long double ratio_velocity = (intersections.col(1).maxCoeff()-intersections.col(1).minCoeff())/(intersections.col(1).minCoeff()+1e-10);

        if (std::abs(ratio_displacement) > 0.1 && std::abs(ratio_velocity) > 0.1) {
            printf("Chaos may have been found for frequency %f! Storing results in file %s \n", frequency, output_file.c_str());
            writeToCSVfile(output_file, intersections);
        } else {
            printf("Nothing interesting found for frequency %f. Moving on...\n", frequency);
        }
    }    
}

void 
single_rate_sine_history(long double const frequency) 
{
    // Physical parameters
    //double const frequency = 19.0;
    long double const displacement = 0.01;

    // Integration parameters
    long double const time_step = 1e-6;
    int const transient_periods = 1000;
    long double const simulation_time = 20.0/frequency;
    int const write_frequency = 1;

    // Setup equation system
    SingleRateSine system(frequency);
    system.d = displacement;
    system.k = 1e4;
    system.m = 1;
    system.p = 100;
    system.cof_static = 1.0;
    system.cof_kinetic = 0.5;
    system.delta = 3.0;
    system.stiffness_damping(frequency, 0.05);
    

    // Storage info
    std::string const output_file = "history_single.csv";
    std::string const file_header = "f: " + std::to_string(frequency);

    long double const period_time = 1.0/frequency;
    int const period_steps = (int)std::floor(period_time/time_step) + 1;
    int const num_time_steps = (int)std::floor(simulation_time/time_step);
    int const num_saves = (int)std::floor((long double)num_time_steps/(long double)write_frequency);

    Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> output(num_saves, 3);  // [Position, Velocity, Time]
    auto state = system.initial_state();

    for (int i=0; i < transient_periods; ++i) {
        long double cycle_time = 0.0;
        for (int k=0; k < period_steps; ++k) {
            long double const dt = std::min(time_step, period_time-cycle_time);
            state += step_rk4(system, state, dt);
            cycle_time += dt;
        }
    }
    
    int const simulation_steps = (int)(simulation_time/time_step) + 1;
    int j = 0;
    for (int i = 0; i < simulation_steps; ++i) {
        state += step_rk4(system, state, time_step);
        if(i % write_frequency == 0) {
            output.row(j) = state;
            j += 1;
        }
    }

    std::cout << "Storing results in file " << output_file << ".\n";
    writeToCSVfile(output_file, output, file_header);
}


void 
single_rate_sine_poincare(long double const frequency) 
{


    // Physical parameters
    //long double const frequency = 30.0;
    long double const displacement = 0.01;

    // Integration parameters
    long double const time_step = 1e-5;
    int const transient_periods = 1000;
    long double const simulation_time = 1.0;
    int const write_frequency = 1;
    int num_intersections = (int)1e6;

    // Setup equation system
    SingleRateSine system(frequency);
    system.d = displacement;
    system.k = 1e4;
    system.m = 1;
    system.p = 100;
    system.cof_static = 1.0;
    system.cof_kinetic = 0.5;
    system.delta = 2.0;
    system.stiffness_damping(frequency, 0.05);
    
    // Storage info
    std::string const output_file = "poincare_single.csv";
    std::string const header = "f: " + std::to_string(frequency);

    long double const period_time = 1.0/frequency;
    int const period_steps = (int)(period_time/time_step) + 1;
    int const num_time_steps = (int)(simulation_time/time_step);
    int const num_saves = (int)((long double)num_time_steps/(long double)write_frequency);
    
    printf("Calculating %d intersections of first return map for single dof (k=%.2f, p=%.2f etc.)\n", num_intersections, system.k, system.p);

    auto intersections = calculate_poincare_sections(
        system, time_step, num_intersections, transient_periods);

    std::cout << "Integration finished! Storing results in file " << output_file << ".\n";
    writeToCSVfile(output_file, intersections);
}