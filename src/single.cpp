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
    double const external_force = -k*state[0] - beta*k*state[1];
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
        //dxdt[1] = 1.0/m * (external_force - friction(relative_velocity)*cof_kinetic*p*sgn(relative_velocity));
        dxdt[1] = 1.0/m * (external_force - friction(relative_velocity)*p*sgn(relative_velocity));
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
    std::cout << "Running Single DOF with rate-dependent friction and sine driver\n";
    double const time_step = 5e-6;
    //int const transient_periods = 200;
    int const num_initial_positions = 1;
    //int const num_intersections = 100;
    //double const frequency = 19.5;
    double const perturbance = 0.0;
    double const damping_ratio = 0.005;
    double const cof_static = 1.0;
    double const cof_dynamic = 0.5;
    double const delta = 1.0; // Friction delta... TODO
    std::string const output_file = "poincare_"+ std::to_string(frequency) + "Hz.txt"; //\
        + "_" + std::to_string(transient_periods) \
        + "_" + std::to_string(num_intersections) \
        + "_xi " +std::to_string(damping_ratio) \
        + "_delta"+std::to_string(delta) \
        +".txt";

    SingleRateSine system(frequency);
    system.stiffness_damping(frequency, damping_ratio);
    system.set_friction(cof_static, cof_dynamic, delta);
    system.k = 1e4;
    system.m = 1.0;
    system.p = 100;

    printf("(k: %f, m: %f, p: %f, xi: %f)\n", system.k, system.m, system.p, damping_ratio);

    auto intersections = calculate_poincare_sections(
        system, time_step, num_initial_positions, num_intersections, transient_periods, perturbance);

    double r1 = (intersections.col(0).maxCoeff() - intersections.col(0).minCoeff())/(intersections.col(0).minCoeff()+1e-10);
    double r2 = (intersections.col(1).maxCoeff() - intersections.col(1).minCoeff())/(intersections.col(1).minCoeff()+1e-10);
    if (std::abs(r1) > 0.1 && std::abs(r2) > 0.1) {
        std::cout << "Chaos may have been found! Storing results in file " << output_file << ".\n";
        writeToCSVfile(output_file, intersections);
    } else {
        std::cout << "Nothing interesting found for frequency " << frequency << " Hz, moving on...\n";
    }

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
    //period_steps *= 4;
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
            double const perturbance = dis(generator) * perturbance;
            Vec3 state(perturbance, initial_velocity, 0.0);

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