#ifndef SINGLE_RATE_SINE_HPP
#define SINGLE_RATE_SINE_HPP

#include "ode.hpp"
#include "utilities.hpp"

#include <string>
#include <memory>
#include <iostream>


struct SingleConstant {
    // Physical properties
    double p = 200.0;
    double m = 1.0;
    double k = 1e5;
    double c = 0.0;

    // Friction settings
    double eps = 1e-4;
    double cof_static = 1.0;
    double cof_kinetic = 0.5;
    double delta = 1.0;

    // Rayleigh parameters
    double alpha = 0.0;
    double beta = 0.0;

    SingleRateSine(double frequency) : f(frequency) {}

    double friction(double const& v_rel) const {
        //return 1.0/(1.0+delta*std::abs(v_rel));
        return cof_kinetic + (cof_static-cof_kinetic)*std::exp(-delta*std::abs(v_rel));
    }

    void stiffness_damping(double const frequency, double const ratio) { 
        beta = ratio/(M_PI*frequency);
    }

    void set_friction(double _cof_static, double _cof_kinetic, double _delta) {
        cof_static = _cof_static;
        cof_kinetic = _cof_kinetic;
        delta = _delta;
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