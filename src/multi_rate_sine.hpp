#ifndef MULTI_RATE_SINE_HPP
#define MULTI_RATE_SINE_HPP

#include "data_structures.hpp"
#include "output.hpp"

#include <vector>

struct HertzRateSine {
    /*
    A N-block system with rate-dependent (velocity-weakening) friction function.
    State-dependent friction is possible through linear scaling related to slip.
    */

    // Physical properties
    int const N = 100;              // Number of blocks
    double const p = 2000.0;        // Peak pressure (at center block)
    int const num_free_blocks = 5;  // Number of free blocks in each end

    double m = 1.0;                 // Block mass
    double k0 = 5000;               // Stiffness to ground
    double k = 1e3;                 // Stiffness between blocks
    double c0 = 50;                 // Damping
    double c = 10;                  // Damping between blocks
    Eigen::VectorXd pressure;       // TODO: unique_ptr?

    // Belt settings
    double d = 0.01;                // Displacement amplitude
    double f = 10.0;                // Displacement frequency

    // Friction settings
    double eps = 1e-4;              // Friction "stick belt" on relative velocity
    double cof_static = 0.75;       // Static CoF
    double cof_kinetic = 0.55;      // Kinetic CoF
    double evolve_cof = 0.0;        // Coefficient for linear CoF evolution (TODO)

    HertzRateSine(int _N, double _p, int _num_free_blocks);

    double friction(double const& v_rel) const {
        return 1.0/(1.0+std::abs(v_rel));
    }

    auto slope(Vec const& state) const -> Vec;
    double scale_friction(double const& vrel) const {return 1.0/std::abs(1.0 + vrel);}  // TODO
    double velocity_at_time(double const& time) const;
    double position_at_time(double const& time) const;
    double calc_displacement_amplitude(Vec const& state) const;
    auto get_initial() const -> Vec; 
    auto shear_force(Vec const& state, int block) const -> double;
    auto resultant_shear_force(Vec const& state) const -> double;
    void set_roughness(double const& height, double const wavelength);
    std::vector<double> calc_natural_frequencies() const;
};

auto calculate_hertz_rate_sine(
    HertzRateSine const& system, 
    double const& time_step, 
    double const& transient_time,
    double const& simulation_time,
    int const& write_frequency) -> Eigen::MatrixXd;

void hertz_rate_sine();
void hertz_rate_sine_shear(double const frequency);
void hertz_rate_sine_slip();

void calculate_multi_poincare_sections();

#endif