#include "multi_rate_sine.hpp"
#include <iostream>
#include <string>
#include "output.hpp"
#include "utilities.hpp"
#include "integrators.hpp"

#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <random>

// Main and only(?) constructor
HertzRateSine::HertzRateSine(int _N, double _p, int _num_free_blocks) 
    : N(_N), p(_p), num_free_blocks(_num_free_blocks)
{
    pressure.resize(N);
    pressure.setZero();
    if (N < num_free_blocks + 3) {
        for(int i=0; i < N; ++i) {
            pressure(i) = p;
        }
    } else {
        int const pressured_blocks = N - 2*num_free_blocks;
        for(int i=0; i < pressured_blocks; ++i) {
            double const pos = -1.0 + 2.0/(double)pressured_blocks * (double)i;
            double const pres = p*std::sqrt(1.0 - pos*pos);
            pressure(num_free_blocks + i) = pres;
        }
    }
    //writeToCSVfile("pressure.txt", pressure);  // For debug
}

void set_roughness(double const& height, double& const wavelength)
{
    std::cerr << "Not implemented!";
}

void HertzRateSine::set_roughness(double const& scale)
{
    std::random_device random_device;
    std::mt19937_64 generator(100);
    std::normal_distribution<double> distribution(0, scale*p);
    for(int i=num_free_blocks; i<N-num_free_blocks; ++i) {
        pressure(i) = pressure(i) + distribution(generator);
    }
}


auto HertzRateSine::slope(Vec const& x) const -> Vec 
{
    /*
    TODO: Remove redundant code and clean up (See friction calculation)
    */
    Vec dxdt(3*N+1);

    double const& time = x(3*N);
    double const belt_position = d*std::sin(2.0*M_PI*f*time);
    double const belt_velocity = 2.0*M_PI*f*d*std::cos(2.0*M_PI*f*time);
    double const belt_acceleration = -std::pow(2.0*M_PI*f, 2.0)*d*std::sin(2.0*M_PI*f*time);

    double external_force = -k0*x(0) - (alpha*m + beta*k0)*x(1);
    if (N > 1) {
        external_force += k*(x(3)-x(0)) + (beta*k)*(x(4) - x(1));  // Todo: Check if (N > 1)* .. is faster
    }
    
    double relative_velocity = x(1) - belt_velocity;
    if (std::abs(relative_velocity) < eps) {
        double const friction_limit = (cof_static + x(2))*pressure(0);
        double const stick_force = std::abs(external_force + m*belt_acceleration);
        if (stick_force <= friction_limit) {
            dxdt(0) = belt_velocity;
            dxdt(1) = belt_acceleration;
        } else {
            dxdt(0) = x(1);
            dxdt(1) = 1.0/m*(external_force - friction_limit*sgn(external_force));
        }
    } else {
        dxdt(0) = x(1);
        double const friction_force =  friction(relative_velocity, x(2))*pressure(0)*sgn(relative_velocity);
        dxdt(1) = 1.0/m*(external_force - friction_force);
    }

    dxdt(2) = evolve_cof*std::abs(shear_force(x, 0)*(x(0) - belt_position));  // TODO Rather based on shear work?

    // Loop over interiour blocks
    if (N > 1) {
        for (int i=1; i < N-1; ++i) { 
            external_force = k*(x(3*i+3) - x(3*i)) + beta*k*(x(3*i+4) - x(3*i+1))
                            -k*(x(3*i) - x(3*i-3)) - beta*k*(x(3*i+1) - x(3*i-2))
                            -k0*x(3*i) - (alpha*m + beta*k0)*x(3*i+1);
            relative_velocity = x(3*i+1) - belt_velocity;
            
            if (std::abs(relative_velocity) < eps) {
                double const friction_limit = (cof_static + x(3*i+2))*pressure(i);
                double const stick_force = std::abs(external_force + m*belt_acceleration);
                if(stick_force <= friction_limit) {
                    dxdt(3*i) = belt_velocity;
                    dxdt(3*i+1) = belt_acceleration;
                } else {
                    dxdt(3*i) = x(3*i+1);
                    dxdt(3*i+1) = 1.0/m*(external_force - friction_limit*sgn(external_force));
                }
            } else {
                dxdt(3*i) = x(3*i+1);
                double const friction_force = friction(relative_velocity, x(3*i+2))*pressure(i)*sgn(relative_velocity);
                dxdt(3*i+1) = 1.0/m*(external_force - friction_force);
            }

            dxdt(3*i+2) = evolve_cof*std::abs(shear_force(x, i)*(x(3*i) - belt_position));
        }

        // Right-end block  // TODO: Check damping coefficients
        external_force = - k*(x(3*N-3)-x(3*N-6)) - beta*k*(x(3*N-2)-x(3*N-5))
                         - k0*x(3*N-3) - (alpha*m + beta*k0)*x(3*N-2);
        relative_velocity = x(3*N-2)- belt_velocity;
        
        if (std::abs(relative_velocity) < eps) {
            double const friction_limit = (cof_static + x(3*N-1)) * pressure(N-1);
            double const stick_force = std::abs(external_force + m*belt_acceleration);
            if(stick_force <= friction_limit) {
                dxdt(3*N-3) = belt_velocity; 
                dxdt(3*N-2) = belt_acceleration;
            } else {
                dxdt(3*N-3) = x(3*N-2);
                dxdt(3*N-2) = 1.0/m * (external_force - friction_limit*sgn(external_force));
            }
        } else {
            dxdt(3*N-3) = x(3*N-2);
            double const friction_force = friction(relative_velocity, x(3*N-1))*pressure(N-1)*sgn(relative_velocity);
            dxdt(3*N-2) = 1.0/m*(external_force - friction_force);
        }
        dxdt(3*N-1) = evolve_cof*std::abs(shear_force(x, N)*(x(3*N-3) - belt_position));
    }

    dxdt(3*N) = 1;  // dt/dt == 1

    return dxdt;
}


auto HertzRateSine::get_initial() const -> Vec
{
    Vec state(3*N+1);
    for (int i = 0; i < N; ++i) {
        state(3*i) = 0.0;
        state(3*i+1) = velocity_at_time(0.0);
        state(3*i+2) = 0.0;
    }
    state(3*N) = 0.0;

    return state;
}


auto HertzRateSine::velocity_at_time(double const& time) const -> double
{
    return 2.0*M_PI*f*d*std::cos(2.0*M_PI*f*time);
}


auto HertzRateSine::position_at_time(double const& time) const -> double
{
    return d*std::sin(2.0*M_PI*f*time);
}


auto HertzRateSine::calc_displacement_amplitude(Vec const& state) const -> double
{    
    double displacement = 0.0;
    for (int i=0; i < N; ++i) {
        displacement += state(3*i);
    }
    return displacement;
}

auto HertzRateSine::friction_force(Vec const& x, int block) const -> double
{
    /*
    TODO: Irradicate behaviour when natural frequency band is wide! Why?
    TODO: Rename to friction_force
    */
    double const& time = x(3*N);
    double const belt_velocity = velocity_at_time(time);
    double const belt_acceleration = -std::pow(2.0*M_PI*f, 2.0)*d*std::sin(2.0*M_PI*f*time);
    double const relative_velocity = x(3*block+1) - belt_velocity; 

    double external_force = 0;
    if (block == 0) {  // Left end
        external_force = -k0*x(0) - (alpha*m + beta*k0)*x(1);
        if (N > 1) {
            external_force += k*(x(3)-x(0)) + (beta*k)*(x(4) - x(1));
        }
    } else if (block == N-1) {  // Right end
        external_force = - k*(x(3*N-3)-x(3*N-6)) - beta*k*(x(3*N-2)-x(3*N-5))
                         - k0*x(3*N-3) - (alpha*m + beta*k0)*x(3*N-2);
    } else {  // In-between
        external_force = k*(x(3*block+3) - x(3*block)) + beta*k*(x(3*block+4) - x(3*block+1))
                    -k*(x(3*block) - x(3*block-3)) - beta*k*(x(3*block+1) - x(3*block-2))
                    -k0*x(3*block) - (alpha*m + beta*k0)*x(3*block+1);        
    }
    
    if (std::abs(relative_velocity) < eps) {
        double const friction_limit = (cof_static + x(3*block+2))*pressure(block);
        double const stick_force = std::abs(external_force + m*belt_acceleration);
        if(stick_force <= friction_limit) {
            return -stick_force*sgn(external_force);
        } else {
            return -friction_limit*sgn(external_force);
        }
    } else {
        double const kinetic_friction = (cof_kinetic + x(3*block+2))*scale_friction(relative_velocity)*pressure(block)*sgn(relative_velocity);
        return kinetic_friction;
    }
}


auto HertzRateSine::resultant_friction_force(Vec const& state) const -> double
{
    double force = 0.0;
    for (auto i=0; i < N; ++i) {
        force += friction_force(state, i);
    }
    return force;
}


auto HertzRateSine::shear_force(Vec const& x, int block) const -> double
{
    /*
    TODO: Rename to friction_force
    */
    double const& time = x(3*N);
    double const belt_position = position_at_time(time);
    double const belt_velocity = velocity_at_time(time);
    double const belt_acceleration = -std::pow(2.0*M_PI*f, 2.0)*d*std::sin(2.0*M_PI*f*time);
    double const relative_velocity = x(3*block+1) - belt_velocity; 

    double external_force = 0;
    if (block == 0) {  // Left end
        external_force = -k0*x(0) - (alpha*m + beta*k0)*x(1);
        if (N > 1) {
            external_force += k*(x(3)-x(0)) + (beta*k)*(x(4) - x(1));
        }
    } else if (block == N-1) {  // Right end
        external_force = - k*(x(3*N-3)-x(3*N-6)) - beta*k*(x(3*N-2)-x(3*N-5))
                         - k0*x(3*N-3) - (alpha*m + beta*k0)*x(3*N-2);
    } else {  // In-between
        external_force = k*(x(3*block+3) - x(3*block)) + beta*k*(x(3*block+4) - x(3*block+1))
                    -k*(x(3*block) - x(3*block-3)) - beta*k*(x(3*block+1) - x(3*block-2))
                    -k0*x(3*block) - (alpha*m + beta*k0)*x(3*block+1);        
    }
    
    if (std::abs(relative_velocity) < eps) {
        double const friction_limit = (cof_static + x(3*block+2))*pressure(block);
        double const stick_force = external_force + m*belt_acceleration;
        if(std::abs(stick_force) <= friction_limit) {
            return -stick_force; //*sgn(external_force + m*belt_acceleration);
        } else {
            return -friction_limit*sgn(external_force);
        }
    } else {
        double const kinetic_friction = friction(relative_velocity, x(3*block+2))*pressure(block)*sgn(relative_velocity);
        return -kinetic_friction;
    }
}


auto HertzRateSine::resultant_shear_force(Vec const& state) const -> double
{
    double force = 0.0;
    for (auto i=0; i < N; ++i) {
        force += shear_force(state, i);
    }
    return force;
}


std::vector<double> HertzRateSine::calc_natural_frequencies() const
{
    /*
    Calculates the natural frequencies of the chain. 
    Returns frequencies in Hz.
    Warning: Expensive using dense matrices
    TODO: Solve sparse system with Spectra??
    */
    // // Set up equation system
    // using Triplet = Eigen::Triplet<double>;
    // Eigen::SparseMatrix<double> K;
    // std::vector<Triplet> tripletlist;
    // tripletlist.reserve(3*(N-2) + 4);
    // tripletlist.emplace_back(Triplet(0, 0, k0+k));
    // tripletlist.emplace_back(Triplet(0, 1, -k));
    // for (auto i=1; i < N-1; ++i) {
    //     tripletlist.emplace_back(Triplet(i, i-1, -k));
    //     tripletlist.emplace_back(Triplet(i, i, 2*k+k0));
    //     tripletlist.emplace_back(Triplet(i, i+1, -k));
    // }
    // tripletlist.emplace_back(Triplet(3*N, 3*N-1, -k));
    // tripletlist.emplace_back(Triplet(3*N, 3*N, k+k0));
    // K.setFromTriplets(tripletlist.begin(), tripletlist.end());

    // Eigen::SparseMatrix<double> M;
    // tripletlist.clear();
    // for(auto i=0; i < N; ++i) {
    //     tripletlist.emplace_back(Triplet(i,i, m));
    // }

    // Eigen::EigenSolver<Eigen::SparseMatrix<double>> solver;

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);
    K(0,0) = k+k0;
    K(0,1) = -k;
    M(0,0) = m;
    for (auto i=1; i < N-1; ++i) {
        K(i,i-1) = -k;
        K(i,i) = 2*k+k0;
        K(i, i+1) = -k;

        M(i,i) = m;
    }
    K(N-1,N-2) = -k;
    K(N-1,N-1) = k+k0;
    M(N-1,N-1) = m;

    // Use Eigen dense solver
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver;

    solver.compute(K, M);
    auto const alphas = solver.alphas();
    auto const betas = solver.betas();
    std::vector<double> eigenvalues;
    for (auto i=0; i<N; ++i) {
        eigenvalues.emplace_back(2.0*M_PI*std::sqrt(std::real(alphas(i)/betas(i)))); 
    }

    return eigenvalues;
}


void HertzRateSine::damping_ratio(double const ratio_lowest, double const ratio_highest)
{
    /*
    Calculates and applies Rayleigh damping term based on damping ratio
    for lowest natural mode and highest natural mode.
    NOTE: that other physical parameters must be set first (mass and stiffnes).

    */

   double const lowest_frequency = lowest_natural_frequency();
   double const highest_frequency = highest_natural_frequenc();

   std::cerr << "Not implemented yet!";

}


double HertzRateSine::lowest_natural_frequency() const
{
    /*
    Returns the lowest natural frequency of the system in radians
    Note: Must be called after other physical properties of the system
    is set.

    This mode corresponds to all block masses oscillating in-phase.

    TODO: Check with calc_natural_frequencies() 
    */
    return std::sqrt(k0/m);
}


double HertzRateSine::highest_natural_frequenc() const 
{
    /*
    Returns the highest natural frequency of the system in radians
    Note: Must be called after other physical properties of the system
    is set.

    This mode corresponds to all block masses oscillating in opposite phase.

    TODO: Check with calc_natural_frequencies() 
    */
    return std::sqrt((4*k + k0)/m);
}


void HertzRateSine::stiffness_damping(double const frequency, double const ratio)
{
    /*
    Set stiffness-proportional damping to be a ratio of critical damping for the given
    frequency [Hz]. For ratio < 1.0 frequency is underdamped, 1 == is critically damped
    and > 1.0 overdamped. Note that other physical quantities must already be set (mass and
    stiffness).

    C = alpha * M + beta * K
    xi = 0.5 * beta * omega
    beta = 2xi/omega = 2 xi /(2*pi f) = xi(pi * f)

    TODO: Verify coeffients

    */
    this->alpha = 0.0;
    this->beta = ratio / (M_PI*frequency);
}


void HertzRateSine::rayleigh_coefficients(double const xi1, double const xi2)
{
    /*
    Set Rayleigh damping coefficients according to damping ratio for first natural frequency
    and last natural frequency.
    */
    auto const w1 = lowest_natural_frequency();
    auto const w2 = highest_natural_frequenc();
    double const denom = std::pow(w2, 2.0) - std::pow(w1, 2.0);
    this->alpha = (2.0*w1*w2)/denom*(w2*xi1-w1*xi2);
    this->beta = 2.0/denom*(w2*xi2 - w1*xi1);
}


/*
Hard-coded utility functions follow. To be written into builder-pattern!
*/


auto calculate_hertz_rate_sine(
    HertzRateSine const& system, 
    double const& time_step, 
    double const& transient_time,
    double const& simulation_time,
    int const& write_frequency) -> Eigen::MatrixXd
{
    /*
    Runs a N-block Hertzian simulator and returns the state space for all blocks
    for a given write frequency. 
    */
    int const num_saves = (int)(simulation_time/time_step/write_frequency)+1;
    auto x = system.get_initial();
    Eigen::MatrixXd storage = Eigen::MatrixXd::Zero(num_saves, 3*system.N+1);  // TODO: Store only a selection of dofs

    double const transient_steps = (int)(transient_time/time_step);
    for(int i=0; i < transient_steps; ++i) {
        x += step_rk4(system, x, time_step); // TODO: Total time not accurate.
    }

    int const simulation_steps = (int)(simulation_time/time_step);
    int j = 0;
    for (int i = 0; i < simulation_steps; ++i) {
        x += step_rk4(system, x, time_step);
        if(i % write_frequency == 0) {
            storage.row(j) = x;
            j += 1;
        }
    }

    return storage;
}


void hertz_rate_sine() 
{
    /*
    Calculates and stores full N-block history at a given write frequency. 
    High memory usage for long histories!
    */
    double const time_step = 1e-5;
    double const frequency = 10.0;
    int const num_blocks = 100;
    int const num_free_blocks = 5;  // On each side 
    double transient_time = 1.0;
    double const simulation_time = 1.0;
    double const peak_pressure = 2000.0;
    std::string const output_file = "history_" + std::to_string(frequency) + "Hz_" + std::to_string(num_blocks) + "blocks.csv";

    HertzRateSine system(num_blocks, peak_pressure, num_free_blocks);
    system.f = frequency;
    system.evolve_cof = 0.0;
    auto result = calculate_hertz_rate_sine(system, time_step, transient_time, simulation_time, 1);

    std::cout << "Storing results in file " << output_file << ".\n";
    writeToCSVfile(output_file, result);
}


void hertz_rate_sine_slip(
    double const frequency, 
    double const delta,
    double const damping_ratio,
    double const pressure)
{
    /*
    Calculates and stores the slip history for left edge, contact edge and contact center
    */
    double const time_step = 1e-7;
    //double const frequency = 19.0;
    int const num_blocks = 100;
    int const num_free_blocks = 5;  // On each side
    //double const pressure = 150; 
    double transient_time = 500/frequency;
    double const simulation_time = 2;
    int const write_frequency = 100;
    double const stiffness = 1e5;
    double const mass = 0.05/num_blocks;
    double const cof_static = 0.75;
    double const cof_kinetic = 0.5;
    //double const delta = 1.0;
    //double const stiffness_damping_ratio = 0.05;

    HertzRateSine system(num_blocks, pressure, num_free_blocks);  // TODO make into builder pattern
    system.f = frequency;
    system.k0 = stiffness/num_blocks;
    system.k = stiffness;
    system.m = mass;
    system.cof_static = cof_static;
    system.cof_kinetic = cof_kinetic;
    system.delta = delta;
    
    double const first_natural_frequency = system.lowest_natural_frequency() / (2*M_PI);
    //system.stiffness_damping(first_natural_frequency, damping_ratio);
    system.rayleigh_coefficients(0.01, 0.05);
    
    auto state = system.get_initial();
    
    //writeToCSVfile("mdof/P1200/pressure_profile", system.pressure);

    std::string const output_history = "mdof/mass_damped/mdof_slip_history_f" + std::to_string(frequency) 
                                   + "_P" + std::to_string(pressure)
                                   + "_xi" + std::to_string(damping_ratio)
                                   + "_del"+ std::to_string(delta) + ".csv";

    std::string const output_loop = "mdof/mass_damped/mdof_fretting_loop_f" + std::to_string(frequency)
                                     + "_P" + std::to_string(pressure)
                                     + "_xi" + std::to_string(damping_ratio)
                                     + "_del"+ std::to_string(delta) +".csv";

    std::string const file_header = "f:" + std::to_string(frequency)
                    + ",xi:" + std::to_string(damping_ratio)
                    + ",k:" + std::to_string(system.k)
                    + ",k0:" + std::to_string(system.k0)
                    + ",disp:" + std::to_string(system.d)
                    + ",m:" + std::to_string(system.m)
                    + ",del:" + std::to_string(system.delta)
                    + ",p:" + std::to_string(system.p)
                    + ",cs:" + std::to_string(system.cof_static)
                    + ",cd:" + std::to_string(system.cof_kinetic)
                    + ",dt:" + std::to_string(time_step)
                    + ",trans_time:" + std::to_string(transient_time);

    int const num_saves = (int)(simulation_time/time_step/write_frequency);
    Eigen::MatrixXd storage = Eigen::MatrixXd::Zero(num_saves, 5);
    Eigen::MatrixXd fretting_map = Eigen::MatrixXd::Zero(num_saves, 2);

    std::cout << "Calculating slip values for end, edge and center for frequency " << frequency << " Hz\n";
    std::cout << "Settings: " << file_header << "\n";

    double const lowest_frequency = system.lowest_natural_frequency();
    double const highest_frequency = system.highest_natural_frequenc();

    std::cout << "Lowest natural frequency: " << std::to_string(lowest_frequency/(2*M_PI)) << " Hz\n"
              << "Highest natural frequency:" << std::to_string(highest_frequency/(2*M_PI)) << "Hz\n";
    
    //std::string const header = "f: " + std::to_string(system.f) + ",p:" + std::to_string(pressure) +"\n";

    int const transient_steps = (int)(transient_time/time_step);
    for(int i=0; i < transient_steps; ++i) {
        state += step_rk4(system, state, time_step);
    }

    int const simulation_steps = (int)(simulation_time/time_step);
    int j = 0;
    for (int i = 0; i < simulation_steps; ++i) {
        state += step_rk4(system, state, time_step);
        if(i % write_frequency == 0) {
            double const& time = state(3*num_blocks);
            double const pad_pos = system.position_at_time(time);
            storage(j, 0) = time;
            storage(j, 1) = state(0) - pad_pos;                                   // Trailing end
            storage(j, 2) = state(3*(system.num_free_blocks + 1)) - pad_pos;      // Trailing edge
            storage(j, 3) = state(3*num_blocks/2) - pad_pos;                      // Contact center   
            storage(j, 4) = system.calc_displacement_amplitude(state) - pad_pos;  // Sum displacement
            fretting_map(j,0) = system.resultant_shear_force(state);
            fretting_map(j,1) = pad_pos; //system.calc_displacement_amplitude(state);
            j += 1;
        }
    }

    std::cout << "Storing slip histories for "<< std::to_string(num_saves) << " steps in file " << output_history << ".\n";
    writeToCSVfile(output_history, storage, file_header);
    writeToCSVfile(output_loop, fretting_map, file_header);
}


void hertz_rate_sine_shear(double const frequency)
{
    /*
    Calculates the resulting shear amplitude
    for a given belt frequency.
    TODO: Shear or friction amplitude??
    */
    double const time_step = 1e-5;
    int const num_blocks = 100;
    int const num_free_blocks = 5;  // On each side 
    double const pressure = 2000.0;
    double transient_time = 100;
    double const simulation_time = 10;
    int const write_frequency = 10;
    std::string const output_file = "shear_" + std::to_string(frequency) + "Hz_"+std::to_string(num_blocks) + "blocks_c005.csv";

    HertzRateSine system(num_blocks, pressure, num_free_blocks);
    system.f = frequency;

    std::cout << "Calculating resulting shear stress " << frequency << " Hz.";

    int const num_saves = (int)(simulation_time/time_step/write_frequency)+1;
    Eigen::MatrixXd storage = Eigen::MatrixXd::Zero(num_saves, 2);
    auto state = system.get_initial();

    int const transient_steps = (int)(transient_time/time_step);
    for(int i=0; i < transient_steps; ++i) {
        state += step_rk4(system, state, time_step); // TODO: Total time not accurate.
    }

    int const simulation_steps = (int)(simulation_time/time_step);
    int j = 0;
    for (int i = 0; i < simulation_steps; ++i) {
        state += step_rk4(system, state, time_step);
        if(i % write_frequency == 0) {
            double const& time = state(3*num_blocks);
            double const pad_pos = system.position_at_time(time);
            storage(j, 0) = time;
            storage(j, 1) = system.resultant_shear_force(state);
            j += 1;
        }
    }

    std::cout << "Storing total shear force for "<< std::to_string(num_saves) << " steps in file " << output_file << ".\n";
    writeToCSVfile(output_file, storage);
}


void hertz_rate_sine_roughness(double const frequency, double const ratio, double const wave)
{
    /*
    Adds a random "roughness" to the pressure profile on a Hertzian contact.
    Calculates 

    */
    double const time_step = 1e-5;
    //double const frequency = 15.0;
    int const num_blocks = 100;
    int const num_free_blocks = 5;  // On each side 
    double const pressure = 2000.0;
    double transient_time = 100;
    double const simulation_time = 10;
    int const write_frequency = 10;
    std::string const output_file = "shear_" + std::to_string(frequency) + "Hz_"+std::to_string(num_blocks) + "blocks.csv";

    HertzRateSine system(num_blocks, pressure, num_free_blocks);  // TODO make factory
    system.f = frequency;

    std::cout << "Calculating resulting shear stress " << frequency << " Hz.";

    int const num_saves = (int)(simulation_time/time_step/write_frequency)+1;
    Eigen::MatrixXd storage = Eigen::MatrixXd::Zero(num_saves, 2);
    auto state = system.get_initial();

    int const transient_steps = (int)(transient_time/time_step);
    for(int i=0; i < transient_steps; ++i) {
        state += step_rk4(system, state, time_step); // TODO: Total time not accurate.
    }

    int const simulation_steps = (int)(simulation_time/time_step);
    int j = 0;
    for (int i = 0; i < simulation_steps; ++i) {
        state += step_rk4(system, state, time_step);
        if(i % write_frequency == 0) {
            double const& time = state(3*num_blocks);
            double const pad_pos = system.position_at_time(time);
            storage(j, 0) = time;
            storage(j, 1) = system.resultant_shear_force(state);
            j += 1;
        }
    }

    std::cout << "Storing total shear force for "<< std::to_string(num_saves) << " steps in file " << output_file << ".\n";
    writeToCSVfile(output_file, storage);
}


void hertz_evolve(
    double const frequency, 
    double const pressure,
    double const delta, 
    double const damping_ratio,
    double const evolve_rate)
{
    /*
    Calculates evolving fretting loops at given snapshots.
    */
    double const time_step = 1e-7;
    int num_steps = 5;
    int num_loops = 5;
    double const time_difference = 100;
    double transient_time = 50/frequency;
    double const simulation_time = 10;
    int const write_frequency = 10;

    int const num_blocks = 100;
    int const num_free_blocks = 5;  // On each side
    //double const pressure = 15; 
    double const stiffness = 1e5;
    double const mass = 0.05/100;  // Per block
    double const cof_static = 0.75;
    double const cof_kinetic = 0.50;
    //double const delta = 1.0;
    double const stiffness_damping_ratio = 0.05;
    double const roughness = 1.0;

    HertzRateSine system(num_blocks, pressure, num_free_blocks);  // TODO make into builder pattern
    system.f = frequency;
    system.k0 = stiffness/num_blocks;
    system.k = stiffness;  // TODO: Increase stiffness with N??
    system.m = mass;
    system.cof_static = cof_static;
    system.cof_kinetic = cof_kinetic;
    system.delta = delta;
    system.evolve_cof = evolve_rate;
    
    double const first_natural_frequency = system.lowest_natural_frequency() / (2*M_PI);
    //system.stiffness_damping(first_natural_frequency, stiffness_damping_ratio);
    //system.set_roughness(roughness);
    system.rayleigh_coefficients(0.01, 0.05);

    //writeToCSVfile("mdof/evolve/pressure_profile_"+std::to_string(roughness)+".csv", system.pressure);
    
    auto state = system.get_initial();

    // std::string const output_history = "mdof/evolve/multi_slip_f" + std::to_string(frequency)
    //                               + "_del" + std::to_string(system.delta)
    //                               + "_xi" + std::to_string(stiffness_damping_ratio)+".csv";

    std::string const output_loop = "mdof/evolve/fretting_loop_f" + std::to_string(frequency)
                                    + "_p" + std::to_string(pressure)
                                    + "_del" + std::to_string(system.delta)
                                    + "_xi" + std::to_string(stiffness_damping_ratio)+".csv";

    std::string const file_header = "f:" + std::to_string(frequency)
                    + ",xi:" + std::to_string(stiffness_damping_ratio)
                    + ",k:" + std::to_string(system.k)
                    + ",k0:" + std::to_string(system.k0)
                    + ",disp:" + std::to_string(system.d)
                    + ",m:" + std::to_string(system.m)
                    + ",del:" + std::to_string(system.delta)
                    + ",p:" + std::to_string(system.p)
                    + ",cs:" + std::to_string(system.cof_static)
                    + ",cd:" + std::to_string(system.cof_kinetic)
                    + ",dt:" + std::to_string(time_step)
                    + ",trans_time:" + std::to_string(transient_time)
                    + ",num_loops:" + std::to_string(num_loops);
    
    int const period_steps = (int)(1.0/(time_step * frequency)) + 1;
    int const num_saves = (int)(num_loops*num_steps*period_steps/write_frequency)+0.5;

    //Eigen::MatrixXd storage = Eigen::MatrixXd::Zero(num_saves, 5);
    Eigen::MatrixXd fretting_map = Eigen::MatrixXd::Zero(num_saves, 5);

    std::cout << "Calculating evolving fretting loop \n";
    std::cout << "Settings: " << file_header << "\n";

    double const lowest_frequency = system.lowest_natural_frequency();
    double const highest_frequency = system.highest_natural_frequenc();

    std::cout << "Lowest natural frequency: " << std::to_string(lowest_frequency/(2*M_PI)) << " Hz\n"
              << "Highest natural frequency:" << std::to_string(highest_frequency/(2*M_PI)) << "Hz\n";
    
    //std::string const header = "f: " + std::to_string(system.f) + ",p:" + std::to_string(pressure) +"\n";

    int const transient_steps = (int)(transient_time/time_step);
    for(int i=0; i < transient_steps; ++i) {
        state += step_rk4(system, state, time_step);
    }

    int j = 0;
    for (int i = 0; i < num_steps; ++i) {
        double wait_time = 0.0;
        while (wait_time < time_difference) {
            state += step_rk4(system, state, time_step);
            wait_time += time_step;
        }
        
        for (int k=0; k < num_loops*period_steps; ++k) {
            if(j == fretting_map.rows()-1) break;
            if(k % write_frequency == 0) {
                double const& time = state(3*num_blocks);
                double const pad_pos = system.position_at_time(time);
                fretting_map(j,0) = time;
                fretting_map(j,1) = system.resultant_shear_force(state);
                fretting_map(j,2) = state(3*(system.num_free_blocks + 1)+2);      // Friction at edge¨
                fretting_map(j,3) = state(3*(int)(num_blocks/4)+2);
                fretting_map(j,4) = state(3*num_blocks/2+2);                      // Friction at center
                j += 1;
            }

            state += step_rk4(system, state, time_step);          
        }

    }
    std::cout << "j: " << j << ", fretting_map.rows():" << fretting_map.rows() << "\n";

    std::cout << "Storing fretting loop for steps in file " << output_loop << ".\n";
    //writeToCSVfile(output_history, storage, file_header);
    writeToCSVfile(output_loop, fretting_map, file_header);
}


void calculate_multi_poincare_sections()
{   
    int const num_blocks = 100;
    double const pressure = 2000;
    int const num_free_blocks = 5;
    HertzRateSine system(num_blocks, pressure, num_free_blocks);
    system.f = 13;

    double const time_step = 1e-5;
    int const num_initial_positions = 5;
    int const num_intersections = (int)2e6;
    int const transient_periods = 100;
    double const& perturbance = 1e-14;
    double const period_time = 1.0/system.f;
    int period_steps = (int)(period_time/time_step) + 1;
    int total_points = num_initial_positions * num_intersections;
    double const initial_velocity = 2.0*M_PI*system.f*system.d;
    std::string const output_file = "poincare_multi_"+ std::to_string(system.f) +"_"+std::to_string(transient_periods)+"_"+std::to_string(num_intersections)+".txt";

    std::mt19937_64 generator(100);  // TODO Which seed?
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    Eigen::MatrixXd intersections(total_points, 2);

    int const num_threads = omp_get_max_threads();

    #pragma omp parallel shared(intersections)
    {   // TODO: Distribute better? (e.g. if one in last thread)
        int const taks_per_thread = static_cast<int>(num_initial_positions/omp_get_num_threads());
        int const thread_id = omp_get_thread_num();
        int const start = thread_id * taks_per_thread;

        int stop;
        if(thread_id != num_threads-1) {
            stop = start + taks_per_thread;
        } else {
            stop = num_initial_positions;
        }  

        for (int i=start; i<stop; ++i) {
            double const perturbance = dis(generator) * perturbance;
            Vec state(3*system.N+1);
            for (int i=0; i < system.N; ++i) {
                state(3*i) = perturbance;
                state(3*i+1) = system.velocity_at_time(0.0);
                state(3*i+2) = 0;
            }

            for (int j=0; j < transient_periods; ++j) {
                double cycle_time = 0.0;
                for (int k=0; k < period_steps; ++k) {
                    double const dt = std::min(time_step, period_time-cycle_time);
                    state += step_rk4(system, state, dt);
                    cycle_time += dt;
                }
            }

            for (int j=0; j < num_intersections; ++j) {
                double cycle_time = 0.0;
                for (int k=0; k < period_steps; ++k) {
                    double const dt = std::min(time_step, period_time - cycle_time);
                    state += step_rk4(system, state, dt);
                    cycle_time += dt;
                } // TODO: Is it necessary to check direction of intersection?? No?
                intersections(i*num_intersections + j, 0) = state(3*(system.num_free_blocks + 1));
                intersections(i*num_intersections + j, 1) = state(3*(system.num_free_blocks + 1)+1);
            }
        }
    }

    std::cout << "Storing results in file " << output_file << ".\n";
    writeToCSVfile(output_file, intersections);
}

