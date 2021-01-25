#ifndef SINGLE_RATE_SINE_HPP
#define SINGLE_RATE_SINE_HPP

#include "ode.hpp"
#include "utilities.hpp"

#include <string>
#include <memory>
#include <iostream>


struct SingleRateSine {
    // Physical properties
    double const p = 2000.0;
    double const m = 1.0;
    double const k = 1e5;
    double const c = 0.1;

    // Driver settings
    double const d = 0.01;
    double const f = 20.0;

    // Friction settings
    double const eps = 1e-4;
    double const cof_static = 0.75;
    double const cof_kinetic = 0.55;

    SingleRateSine(double frequency) : f(frequency) {}

    double friction(double const& v_rel) const {
        return 1.0/(1.0+std::abs(v_rel));
    }

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


auto calculate_poincare_sections(
    SingleRateSine const& system, 
    double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    double const& perturbance
) -> Eigen::MatrixXd;


#endif