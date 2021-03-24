#ifndef SINGLE_RATE_SINE_HPP
#define SINGLE_RATE_SINE_HPP

#include "ode.hpp"
#include "utilities.hpp"

#include <string>
#include <memory>
#include <iostream>


struct SingleRateSine {
    /*
    A single oscillator with rate-dependent (velocity-weakening)
    friction law and harmonic (sine) belt.
    */

    // Physical properties
    long double p = 2000.0;            // Pressure
    long double m = 0.01;               // Mass
    long double k = 1e4;               // Stiffness
    //long //double c = 500;             // Damping (Depricated in favour of Rayleigh coefficients)
    long double alpha = 0.0;           // Rayleigh mass-proportional damping
    long double beta = 0.0;            // Rayleigh stiffness-proportional damping
    //long // Driver settings
    long double d = 0.01;              // Displacement amplitude for belt
    long double f = 20.0;              // Frequency of belt

    // Friction settings
    long double eps = 1e-4;            // Limit relative velocity for stick
    long double cof_static = 0.75;      // Static coefficient of friction
    long double cof_kinetic = 0.50;    // Kinetic coefficient of friction
    long double delta = 1.0;

    SingleRateSine(long double const& frequency) : f(frequency) {}

    long double friction(long double const& v_rel) const {
        //return 1.0/(1.0+delta*std::abs(v_rel));
        return cof_kinetic + (cof_static - cof_kinetic)/(1.0 + delta*std::abs(v_rel));
    }

    void stiffness_damping(long double const frequency, long double const ratio);
    
    auto natural_frequency() const -> double;
    auto initial_state() const -> Vec3;
    auto velocity_at_time(long double const& time) const -> long double;
    auto slope(Vec3 const& state) const -> Vec3;
};


auto calculate_poincare_sections(
    SingleRateSine const& system, 
    long double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    long double const& perturbance
) -> Eigen::MatrixXd;


void single_rate_sine();


void single_rate_sine(
    long double frequency,
    int transient_periods, 
    int num_intersections
);


void single_rate_sine(std::string const& input_file);


auto calculate_poincare_sections_perturbed(
    SingleRateSine const& system, 
    long double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    long double const& perturbance
) -> Eigen::MatrixXd;


auto calculate_poincare_sections(
    SingleRateSine const& system, 
    long double const& time_step, 
    int const& num_intersections, 
    int const& transient_periods
) -> Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;

void single_poincare_chaos_finder(long double const frequency_start, long double const frequency_stop);

auto single_rate_sine_history(
    long double const frequency, 
    long double const pressure,
    long double const ratio, 
    long double const delta
) -> void;

void single_rate_sine_poincare(long double const frequency, long double const damping_ratio);

#endif