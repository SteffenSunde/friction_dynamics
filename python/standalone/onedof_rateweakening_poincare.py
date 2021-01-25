import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
import time
from numba import jit, int32, float64
import random

# Static variables :(
pressure = 200.0              # Pressure on block
mass = 1                      # Mass [kg]
stiffness = 100000.0          # Stiffness [N/mm]
damping = 0.1                 # Viscous damping coefficient [-]
displacement = 0.01           # Displacement amplitude [mm]
frequency = 18.0              # Frequency [Hz]
eps = 1e-4                    # Width of "stick-belt" [mm/s]
cof_static = 0.75             # Static coefficient of friction [-]
cof_kinetic = 0.4             # Kinetic coefficient of friction [-]


def main():
    # Integration parameters
    time_step = 0.00001
    transient_periods = 100
    num_intersections = 10000
   
    print(f"Caulculation Poincaré map with {num_intersections} intersections at frequency {frequency} Hz. {transient_periods} periods are discared, so could take some time...")

    start_time = time.perf_counter()
    intersections = calculate_poincare(time_step, num_intersections, transient_periods)
    elapsed_time = time.perf_counter() - start_time

    print("Integration finished in {:.2f} seconds".format(elapsed_time))
    
    result_file = "poincare_{}_{}.txt".format(num_intersections, frequency)

    print(f"Saving results in file {result_file}")
    #np.savetxt(result_file, intersections)

    # Plot
    fig, ax = plt.subplots()
    ax.scatter(intersections[:, 0], intersections[:,1], c='black', s=0.25)

    plt.grid()
    ax.set_title("Poincaré portrait at {:.2f} Hz, m/s".format(frequency))
    plt.show()

    return


@jit
def fric(vrel):
    return 1/(1+abs(vrel))

@jit
def slope(x):
    # Calculates the tangent of the system at current state x
    # x is state of system: [position, velocity, time]
    dxdt = np.array([x[1], 0.0, 0.0], dtype=np.float64)

    # TODO: Velocity, and position must be integrated if frequency is changed throughout analysis
    belt_velocity = 2*pi*frequency*displacement*cos(2*pi*frequency*x[2])
    belt_acceleration = -(2*pi*frequency)**2 * displacement * sin(2*pi*frequency*x[2])
    relative_velocity = x[1] - belt_velocity
    external_force = -stiffness*x[0] - damping*x[1]
    if abs(relative_velocity) < eps:
        friction_limit = cof_static * pressure
        if(abs(external_force + mass * belt_acceleration) <= friction_limit):  # Stick
            dxdt[0] = belt_velocity
            dxdt[1] = belt_acceleration
        else:  # Transition from stick to slip
            dxdt[1] = 1.0/mass * (external_force - friction_limit*np.sign(external_force))
    else:  # Slip
        dxdt[1] = 1.0/mass * (external_force - fric(relative_velocity)*cof_kinetic * pressure * np.sign(relative_velocity))
    
    dxdt[2] = 1  # dt/dt == 1

    return dxdt


@jit(nopython=True)
def step_rk4(x, h):
    # Returns the slope of the system at x using classical 4th order Runge-Kutta
    k1 = slope(x)
    k2 = slope(x + 0.5*h*k1)
    k3 = slope(x + 0.5*h*k2)
    k4 = slope(x + h*k3)

    return 1/6*h*(k1 + 2*k2 + 2*k3 + k4)


@jit
def calculate_poincare(time_step, num_intersections, transient_periods):
    period_time = 1.0/frequency
    period_steps = int(period_time/time_step) + 1
    num_initial_positions = 100
    num_points = num_initial_positions*num_intersections
    intersections = np.zeros((num_points, 2), dtype=np.float64)
    
    for i in range(num_initial_positions):
        perturbance = random.random()*1e-10
        x = np.array([perturbance, 2*pi*frequency*displacement, .0], dtype=np.float64)
        
        for _ in range(transient_periods):
            cycle_time = 0.0
            for _ in range(period_steps):
                h = min(time_step, period_time - cycle_time)
                x += step_rk4(x, h)
                cycle_time += h

        for j in range(num_intersections):
            cycle_time = 0.0
            for _ in range(period_steps):
                h = min(time_step, period_time - cycle_time)
                x += step_rk4(x, h)
                cycle_time += h
            intersections[i*num_intersections + j, 0] = x[0]
            intersections[i*num_intersections + j, 1] = x[1]

    return intersections

if __name__ == "__main__":
    main()