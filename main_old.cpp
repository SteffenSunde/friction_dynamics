#include <iostream>
#include <fstream>
#include <random>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include <boost/program_options.hpp>
#include <Eigen/Dense>
#include "omp.h"

namespace po = boost::program_options;

auto calculate_poincare_sections(
    SingleBlockRateSine const& system, 
    double const& time_step, 
    int const& num_initial_positions, 
    int const& num_intersections, 
    int const& transient_periods, 
    double const& perturbance) -> Eigen::MatrixXd;

auto step_rk4(SingleBlockRateSine const& system, State const& x, double const& dt) -> State
{
    State const k1 = system.slope(x);
    State const k2 = system.slope(x + 0.5*dt*k1);
    State const k3 = system.slope(x + 0.5*dt*k2);
    State const k4 = system.slope(x + dt*k3);

    return 1.0/6.0 * dt * (k1 + 2*k2 + 2*k3 + k4);
}


