#ifndef SINGLE_RATE_SINE_HPP
#define SINGLE_RATE_SINE_HPP

#include "ode.hpp"
#include "utilities.hpp"

#include <string>
#include <memory>
#include <iostream>


struct SingleRate {
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
    double v0 = 1.0;              // Constant velocity of belt

    // Friction settings
    double eps = 1e-4;            // Width of "stick window" on relative velocity
    double cof_static = 0.7;      // Static coefficient of friction
    double cof_kinetic = 0.55;    // Kinetic coefficient of friction

    //SingleRate(double const& frequency) : f(frequency) {}

    inline double friction(double const& v_rel) const {
        return 1.0/(1.0+std::abs(v_rel));
    }

    inline void stiffness_damping(double const frequency, double const ratio);

    auto slope(Vec3 const& state) const -> Vec3;
};

#endif