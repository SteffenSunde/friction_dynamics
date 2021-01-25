#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include "ode.hpp"
#include "single_rate_sine.hpp"

#include <Eigen/Dense>


template<typename ODE, typename State>
auto step_rk4(
    ODE const& system, 
    State const& x, 
    double const& dt) -> State
{
    State const k1 = system.slope(x);
    State const k2 = system.slope(x + 0.5*dt*k1);
    State const k3 = system.slope(x + 0.5*dt*k2);
    State const k4 = system.slope(x + dt*k3);

    return 1.0/6.0 * dt * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

#endif