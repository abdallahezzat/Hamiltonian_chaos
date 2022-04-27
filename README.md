# Hamiltonian chaos
Exploring the behaviour of stars orbits

## The program

This code plots 3D orbits and Poincaré maps of _N_ stars in a galaxy assuming a logarithmic potential with parameters _b_ and _c_.

## Paramters and initial conditions
The number of stars _N_, timestep _dt_, total time _T_, as well as the parameters _b_ and _c_ are specified in **data.txt**.

The initial conditions for the stars (positions and velocities) can be manually changed in **line 36** of **main.py**.

## The integrator
The integrators used for now are _Euler_ and _4th-order Runge-Kutta_ integrators. You can choose the method in **line 58** of **main.py** by typing _euler_ or _rk4_.

## Usage

Running the file main.py will produce three plots for the stars trajectories, Poincaré maps and energy.
