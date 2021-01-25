#ifndef MULTI_RATE_SINE_HPP
#define MULTI_RATE_SINE_HPP

#include "data_structures.hpp"

struct HertzRateSine {
    // Physical properties
    int const N = 100;               // Number of blocks
    int const num_free_blocks = 5;
    double p = 200.0;        // Peak pressure (at center block)
    double m = 1.0;           // Block mass
    double k0 = 5000;          // Stiffness to ground
    double k = 1e3;           // Stiffness between blocks
    double c0 = 0.0;          // Damping
    double c = 0.0;           // Damping between blocks
    Eigen::VectorXd pressure;  // TODO: unique_ptr?

    // Driver settings
    double d = 0.01;
    double f = 10.0;

    // Friction settings
    double eps = 1e-4;
    double cof_static = 0.75;
    double cof_kinetic = 0.55;
    double evolve_cof = 0.0;

    HertzRateSine(int N) {  // TODO: Add free blocks
        pressure.resize(N);
        pressure.setZero();
        if (N < num_free_blocks + 3) {
            for(int i=0; i < N; ++i) {
                pressure(i) = p;
            }
        } else {
            int const pressured_blocks = N - 2*num_free_blocks;
            for(int i=0; i < pressured_blocks; ++i) {
                double const pos = -0.5 + 1.0/(double)pressured_blocks * (double)i;
                double const pres = p*std::sqrt(0.5 - 2.0*pos*pos);
                pressure(num_free_blocks + i) = pres;  
            }
        }
    }

    double friction(double const& v_rel) const {
        return 1.0/(1.0+std::abs(v_rel));
    }

    auto slope(Vec const& state) const -> Vec;
    double scale_friction(double const& vrel) const {return 1.0/std::abs(1.0 + vrel);}  // TODO
    double velocity_at_time(double const& time) const;
    double position_at_time(double const& time) const;
    double calc_displacement_amplitude(Vec const& state) const;
    auto get_initial() const -> Vec; 
    auto shear_force(Vec const& state, int block) -> double;
    auto resultant_shear_force(Vec const& state) -> double;
};

auto calculate_hertz_rate_sine(
    HertzRateSine const& system, 
    double const& time_step, 
    double const& transient_time,
    double const& simulation_time,
    int const& write_frequency) -> Eigen::MatrixXd;

void hertz_rate_sine();
void hertz_rate_sine_shear();

void calculate_multi_poincare_sections();

#endif