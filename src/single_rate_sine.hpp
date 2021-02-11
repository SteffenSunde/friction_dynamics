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
    double p = 2000.0;            // Pressure
    double m = 1.0;               // Mass
    double k = 1e5;               // Stiffness
    //double c = 500;             // Damping (Depricated in favour of Rayleigh coefficients)
    double alpha = 0.0;           // Rayleigh mass-proportional damping
    double beta = 0.0;            // Rayleigh stiffness-proportional damping
    // Driver settings
    double d = 0.01;              // Displacement amplitude for belt
    double f = 20.0;              // Frequency of belt

    // Friction settings
    double eps = 1e-4;            // Width of "stick window" on relative velocity
    double cof_static = 0.7;      // Static coefficient of friction
    double cof_kinetic = 0.55;    // Kinetic coefficient of friction

    SingleRateSine(double const& frequency) : f(frequency) {}

    inline double friction(double const& v_rel) const {
        return 1.0/(1.0+std::abs(v_rel));
    }

    inline void stiffness_damping(double const frequency, double const ratio);

    auto slope(Vec3 const& state) const -> Vec3;
};


auto calculate_poincare_sections(
    SingleRateSine const& system, 
    double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    double const& perturbance
) -> Eigen::MatrixXd;


void single_rate_sine();


void single_rate_sine(
    double frequency,
    int transient_periods, 
    int num_intersections
);


void single_rate_sine(std::string const& input_file);


auto calculate_poincare_sections_perturbed(
    SingleRateSine const& system, 
    double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    double const& perturbance
) -> Eigen::MatrixXd;


auto calculate_poincare_sections(
    SingleRateSine const& system, 
    double const& time_step, 
    int const& num_intersections, 
    int const& transient_periods
) -> Eigen::MatrixXd;

void single_poincare_chaos_finder(double const frequency_start, double const frequency_stop);

#endif