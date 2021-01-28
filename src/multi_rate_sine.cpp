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
            double const pos = -0.5 + 1.0/(double)pressured_blocks * (double)i;
            double const pres = p*std::sqrt(0.5 - 2.0*pos*pos);
            pressure(num_free_blocks + i) = pres;  
        }
    }
    //writeToCSVfile("pressure.txt", pressure);  // For debug
}


auto HertzRateSine::slope(Vec const& x) const -> Vec 
{
    Vec dxdt(3*N+1);

    double const& time = x(3*N);
    double const belt_position = d*std::sin(2.0*M_PI*f*time);
    double const belt_velocity = 2.0*M_PI*f*d*std::cos(2.0*M_PI*f*time);
    double const belt_acceleration = -std::pow(2.0*M_PI*f, 2.0)*d*std::sin(2.0*M_PI*f*time);

    double external_force = -k0*x(0) - c0*x(1);
    if (N > 1) {
        external_force += k*(x(3)-x(0)) + c*(x(4) - x(1));
    }
    
    //auto friction = [&](double v_rel) { return mu_d + (mu_s-mu_d)*std::exp(-std::abs(v_rel)/0.5); };
    double relative_velocity = x(1) - belt_velocity;
    if (std::abs(relative_velocity) < eps) {  // TODO More should be required for block to stick?? no
        double friction_limit = (cof_static + x(2))*pressure(0);  // TODO: log function?
        double stick_force = std::abs(external_force + m*belt_acceleration);
        if (stick_force <= friction_limit) {
            dxdt(0)= belt_velocity;
            dxdt(1) = belt_acceleration;
        } else {
            dxdt(0) = x(1);
            dxdt(1) = 1.0/m*(external_force - friction_limit*sgn(external_force));  // Should acceleration be checked here?
        }
    } else {
        dxdt(0) = x(1);
        double kinetic_friction =  (cof_kinetic + x(2))*scale_friction(relative_velocity)*pressure(0)*sgn(relative_velocity);
        dxdt(1) = 1.0/m*(external_force - kinetic_friction);
    }

    dxdt(2) = evolve_cof*std::abs(x(0) - belt_position); // TODO: Friction work rather?? Ruiz??

    // Loop over interiour blocks
    if (N > 1) { // TODO Check state coefficients! Step through for N=3
        for (int i=1; i < N-1; ++i) { 
            external_force = k*(x(3*i+3) - x(3*i)) + c*(x(3*i+4) - x(3*i+1))
                            -k*(x(3*i) - x(3*i-3)) - c*(x(3*i+1) - x(3*i-2))
                            -k0*x(3*i) - c0*x(3*i+1);
            relative_velocity = x(3*i+1) - belt_velocity;
            
            if (std::abs(relative_velocity) < eps) {
                double const friction_limit = (cof_static + x(3*i+2))*pressure(i);
                double const stick_force = std::abs(external_force + m*belt_acceleration);
                if(stick_force <= friction_limit) {
                    dxdt(3*i) = belt_velocity;
                    dxdt(3*i+1) = belt_acceleration;
                } else {
                    dxdt(3*i) = x(3*i+1);
                    dxdt(3*i+1) = 1.0/m*(external_force - friction_limit*sgn(external_force));  // *sgn(v_rel)TODO Check. Should belt acc be accounted for?
                }
            } else {
                dxdt(3*i) = x(3*i+1);
                double const kinetic_friction = (cof_kinetic + x(3*i+2))*scale_friction(relative_velocity)*pressure(i)*sgn(relative_velocity);
                dxdt(3*i+1) = 1.0/m*(external_force - kinetic_friction);
            }

            dxdt(3*i+2) = evolve_cof*std::abs(x(3*i) - belt_position);
        }

        // Right-end block
        external_force = -k*(x(3*N-3)-x(3*N-6))-c*(x(3*N-2)-x(3*N-5))-k*x(3*N-3)-c*x(3*N-2);
        relative_velocity = x(3*N-2)- belt_velocity;
        
        if (std::abs(relative_velocity) < eps) {
            double const friction_limit = (cof_static + x(3*N-1)) * pressure(N-1);
            double const stick_force = abs(external_force + m*belt_acceleration);  // TODO Check sign
            if(stick_force <= friction_limit) {
                dxdt(3*N-3) = belt_velocity; 
                dxdt(3*N-2) = belt_acceleration;
            } else {
                dxdt(3*N-3) = x(3*N-2);
                dxdt(3*N-2) = 1.0/m * (external_force - friction_limit*sgn(external_force));  // *sgn(v_rel) TODO sgn(shear)
            }
        } else {
            dxdt(3*N-3) = x(3*N-2);
            double const kinetic_friction = (cof_kinetic + x(3*N-1))*scale_friction(relative_velocity)*pressure(N-1)*sgn(relative_velocity);
            dxdt(3*N-2) = 1.0/m*(external_force - kinetic_friction);
        }
        dxdt(3*N-1) = evolve_cof*std::abs(x(3*N-3) - belt_position);
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


auto HertzRateSine::shear_force(Vec const& state, int block) const -> double
{
    /*
    TODO: Check if correct. Must check dynamic equilbrium in stick condition??
    */
    double const& time = state(3*N);
    double const pad_velocity = velocity_at_time(time);
    double const relative_velocity = pad_velocity - state(3*block+1);
    if (std::abs(relative_velocity) < eps) {
        return (cof_static + state(3*block+2))*pressure(block)*sgn(pad_velocity);
    } else {
        return (cof_kinetic + state(3*block+2))*relative_velocity*pressure(block)*friction(relative_velocity)*sgn(pad_velocity);
    }
    // double shear_force = 0.0;
    // if(block == 0) {
    //     double external_force = -k0*state(0) - c0*state(1) 
    //         + k*(state(3) - state(0)) 
    //         + c * (state(4) - state(1));
    //     double const relative_velocity = pad_velocity - state(1);

    //     if (std::abs(relative_velocity) < eps) { // Stick
    //         // double const friction_limit = (cof_static + state(3))*pressure(0);
    //         // if (external_force > friction_limit) {  // Breaking
    //         //     shear_force = 
    //         // } else {
    //         //     shear_force = 0.0;
    //         // }
    //         shear_force = (cof_static + state(2))*pressure(0);
    //     } else {  // Slipping
    //         shear_force = (cof_kinetic + state(2))*relative_velocity*pressure(0)*friction(relative_velocity);
    //     }
    // } else if (block == N) {

    // }
    
    //return shear_force;
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

    double const simulation_steps = (int)(simulation_time/time_step);
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


void hertz_rate_sine_slip()
{
    /*
    Calculates and stores the slip history for left edge, contact edge and contact center
    */
    double const time_step = 1e-5;
    double const frequency = 13.0;
    int const num_blocks = 100;
    int const num_free_blocks = 5;  // On each side
    double const pressure = 2000; 
    double transient_time = 100.0;
    double const simulation_time = 10;
    int const write_frequency = 10;
    std::string const output_file = "slip_" + std::to_string(frequency) + "Hz_"+std::to_string(num_blocks) + "blocks_b001.csv";

    HertzRateSine system(num_blocks, pressure, num_free_blocks);  // TODO make into builder pattern
    system.f = frequency;
    auto state = system.get_initial();

    int const num_saves = (int)(simulation_time/time_step/write_frequency)+1;
    Eigen::MatrixXd storage = Eigen::MatrixXd::Zero(num_saves, 5);

    std::cout << "Calculating slip values for end, edge and center for frequency " << frequency << " Hz\n";

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
            j += 1;
        }
    }

    std::cout << "Storing total shear for "<< std::to_string(num_saves) << " steps in file " << output_file << ".\n";
    writeToCSVfile(output_file, storage);
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

