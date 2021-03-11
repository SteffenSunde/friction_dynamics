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

auto SingleRateSine::natural_frequency() const -> double
{
    return sqrt(k/m)/2.0/M_PI;
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
    //beta = ratio/(M_PI*frequency);
    beta = 4.0*M_PI*ratio/frequency;
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
    for (int i=0; i < transient_periods; ++i) {
        long double cycle_time = 0.0;
        // for (int k=0; k < period_steps; ++k) {
        while (cycle_time < period_time) {
            long double const dt = std::min(time_step, period_time-cycle_time);
            cycle_time += dt;
            //if (std::abs(dstate[2]-dt) > 1e-9) throw std::string("ERROR: Drift");
            state += step_rk4(system, state, dt);
        }
        acc_drift += period_time - cycle_time;
    }
    
    printf("Transient drift: %e\n", state[2]/period_time-(long double)transient_periods);

    if (half_period) {
        long double cycle_time = 0.0;
        for (int k=0; k < (int)(period_steps/2)+1; ++k) {
            long double const dt = std::min(time_step, (period_time/2.0)-cycle_time);
            cycle_time += dt;
            state += step_rk4(system, state, dt);
        }
    }

    for (int i=0; i < num_intersections; ++i) {
        long double cycle_time = 0.0;
        //for (int k=0; k < period_steps; ++k) {
        while(cycle_time < period_time) {
            long double const dt = std::min(time_step, period_time - cycle_time);
            state += step_rk4(system, state, dt);
            cycle_time += dt;
        }
        acc_drift += period_time - period_time;
        intersections.row(i) = state;
        // intersections(i, 0) = state[0];
        // intersections(i, 1) = state[1];
    }
    double time_drift = state[2] - (long double)(num_intersections + transient_periods)*period_time;
    if (half_period) { time_drift -= (long double)0.5*period_time; }
    printf("period_time - cycle_time: %e\n", acc_drift);
    printf("Calculated %d intersections with frequency %.2f and time step %e. Accumulated time drift: %e\n",
        num_intersections, system.f, time_step, time_drift);  //std::fmod(state[2], period_time)

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
    long double const time_step = 1e-7;
    int const transient_periods = 5000;
    int const num_intersections = (int)1e6;
    long double const damping_ratio = 0.005;

    std::mt19937_64 generator(102);  // TODO Which seed?
    //std::mt19937_64 generator2(101);
    std::uniform_real_distribution<long double> distribution(frequency_start, frequency_stop);
    std::uniform_real_distribution<long double> pressure_distribution(1250.0, 1500.0);

    printf("Searching for interesting Poincaré maps between frequencies %.1f and %.1f Hz...\n", frequency_start, frequency_stop);

    //#pragma omp parallel for
    for (int i = 0; i < (int)1e6; ++i) {
        long double frequency = distribution(generator);  // TODO: mt19937_64 thread safe? Make thread id part of seed?
        std::string const output_file = "chaos/del05/poincare_"+ std::to_string(frequency)+".txt";

            SingleRateSine system(frequency);
            system.d = 0.01;
            system.k = 1e5;
            system.m = 0.05;
            system.p = pressure_distribution(generator);
            system.cof_static = 0.75;
            system.cof_kinetic = 0.5;
            system.delta = 0.50;
            system.stiffness_damping(system.natural_frequency(), damping_ratio);

        std::string const file_header = "f:" + std::to_string(frequency)   // TODO Also store seed?
                                + ",xi:" + std::to_string(damping_ratio)
                                + ",k:" + std::to_string(system.k)
                                + ",d:" + std::to_string(system.d)
                                + ",m:" + std::to_string(system.m)
                                + ",del:" + std::to_string(system.delta)
                                + ",p:" + std::to_string(system.p)
                                + ",cs:" + std::to_string(system.cof_static)
                                + ",cd:" + std::to_string(system.cof_kinetic)
                                + ",dt:" + std::to_string(time_step)
                                + ",trans:" + std::to_string(transient_periods);
        
        auto intersections = calculate_poincare_sections(
            system, time_step, num_intersections, transient_periods);

        long double const ratio_displacement = (intersections.col(0).maxCoeff()-intersections.col(0).minCoeff())/system.d;
        long double const ratio_velocity = (intersections.col(1).maxCoeff()-intersections.col(1).minCoeff())/system.d;

        if (std::abs(ratio_displacement) > 0.1 && std::abs(ratio_velocity) > 0.1) {
            printf("Chaos may have been found for frequency %f! Storing results in file %s \n\n", frequency, output_file.c_str());
            writeToCSVfile(output_file, intersections, file_header);
        } else {
            printf("Nothing interesting found for frequency %f. Moving on...\n\n", frequency);
        }
    }    
}

 
auto single_rate_sine_history(
    long double const frequency,    
    long double const damping_ratio,
    long double const delta
) -> void
{
    // Physical parameters
    long double const displacement = 0.01;

    // Integration parameters
    long double const time_step = 1e-7;
    int const transient_periods = 4000;
    long double const simulation_time = 10/frequency;
    int const write_frequency = 100;
    //long double const damping_ratio = 0.05;

    // Setup equation system
    SingleRateSine system(frequency);
    system.d = displacement;
    system.k = 1e5;
    system.m = 0.05;
    system.p = 1200.0;
    system.cof_static = 0.75;
    system.cof_kinetic = 0.5;
    system.delta = delta;
    double const natural_frequency = system.natural_frequency();
    system.stiffness_damping(natural_frequency, damping_ratio);
    
    // Storage info
    std::string const output_file = "chaos/history_single_xi"+std::to_string(damping_ratio)+"_trans"+std::to_string(transient_periods) +".csv";
    std::string const file_header = "f:" + std::to_string(frequency)
                            + ",xi:" + std::to_string(damping_ratio)
                            + ",k:" + std::to_string(system.k)
                            + ",disp:" + std::to_string(system.d)
                            + ",m:" + std::to_string(system.m)
                            + ",del:" + std::to_string(system.delta)
                            + ",p:" + std::to_string(system.p)
                            + ",cs:" + std::to_string(system.cof_static)
                            + ",cd:" + std::to_string(system.cof_kinetic)
                            + ",dt:" + std::to_string(time_step)
                            + ",trans:" + std::to_string(transient_periods);

    long double const period_time = 1.0/frequency;
    int const period_steps = (int)(period_time/time_step) + 1;
    int const num_time_steps = (int)(simulation_time/time_step) + 1;
    int const num_saves = (int)(num_time_steps/write_frequency) + 1;

    Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> output(num_saves, 3);  // [Position, Velocity, Time]
    auto state = system.initial_state();

    std::cout << "System: "<< file_header << "\n";

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


auto single_rate_sine_poincare(
    long double const frequency,
    long double const damping_ratio
) -> void
{
    // Physical parameters
    long double const displacement = 0.01;
    //long double const damping_ratio = 0.0001;

    // Integration parameters
    long double const time_step = 1e-7;
    int const transient_periods = 4000;
    int num_intersections = (int)1e4;

    // Setup equation system
    SingleRateSine system(frequency);
    system.d = displacement;
    system.k = 1e5;
    system.m = 0.05;
    system.p = 1200; //1255.301837;
    system.cof_static = 0.75;
    system.cof_kinetic = 0.5;
    system.delta = 1.0;
    double const natural_frequency = system.natural_frequency();
    system.stiffness_damping(natural_frequency, damping_ratio);
    
    // Storage info
    std::string const output_file = "chaos/poincare_single_f" +std::to_string(frequency) +
                                    "xi_"+ std::to_string(damping_ratio)+".csv";
    std::string const file_header = "f:" + std::to_string(frequency)
                        + ",xi:" + std::to_string(damping_ratio)
                        + ",k:" + std::to_string(system.k)
                        + ",disp:" + std::to_string(system.d)
                        + ",m:" + std::to_string(system.m)
                        + ",del:" + std::to_string(system.delta)
                        + ",p:" + std::to_string(system.p)
                        + ",cs:" + std::to_string(system.cof_static)
                        + ",cd:" + std::to_string(system.cof_kinetic)
                        + ",trans:" + std::to_string(transient_periods);

    long double const period_time = (long double)1.0/frequency;
    int const period_steps = (int)(period_time/time_step) + 1;

    std::cout << "System: " << file_header << "\n";

    auto intersections = calculate_poincare_sections(
        system, time_step, num_intersections, transient_periods);

    std::cout << "Integration finished! Storing results in file " << output_file << ".\n";
    writeToCSVfile(output_file, intersections, file_header);
}