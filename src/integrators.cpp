#include "integrators.hpp"
#include "ode.hpp"
#include "single_rate_sine.hpp"

// template<typename Derived, typename ODE>
// auto step_rk4(
//     ODE const& system, 
//     Eigen::EigenBase<Derived> const& x, 
//     double const& dt) -> Eigen::EigenBase<Derived>
// {
//     auto const k1 = system->slope(x);
//     auto const k2 = system->slope(x + 0.5*dt*k1);
//     auto const k3 = system->slope(x + 0.5*dt*k2);
//     auto const k4 = system->slope(x + dt*k3);

//     return 1.0/6.0 * dt * (k1 + 2*k2 + 2*k3 + k4);
// }

// // Explicit template instantiation
// template auto step_rk4(SingleRateSine const& system, Vec3 const& x, double const& dt) -> Vec3;
