# Vessel Motion Simulation with Water Ingress

## Overview

This C++ program simulates the motion of a hollow vessel with water gradually entering through a hole at the bottom. The simulation solves a coupled system of differential equations using numerical methods to model the vessel's displacement, velocity, and water level over time.

## Problem Description

The system models a hollow vessel with the following characteristics:
- Mass `m` with infinite length and negligible wall thickness
- Water enters through a hole at the bottom at a flow rate proportional to the water level difference
- The vessel experiences weight, buoyancy, and drag forces
- Motion is governed by coupled first and second-order differential equations

### Governing Equations

**Water Level Equation:**
```
dh/dt = (K/(ρA)) * (ȳ + y(t) - h(t))
```

**Vessel Motion Equation:**
```
dy/dt = (ρA/(m + ρAh(t))) * [2gh(t) - gy(t) - (1/2)Cₐ(dy/dt)|dy/dt|]
```

Where:
- `h(t)` = water level inside the vessel
- `y(t)` = vessel displacement from equilibrium
- `ȳ = m/(ρA)` = equilibrium position of empty vessel
- `K` = flow rate proportionality constant
- `ρ` = water density
- `A` = cross-sectional area
- `Cₐ` = drag coefficient
- `g` = gravitational acceleration

## Project Structure

```
 src/
    main.cpp         # Main program with Vessel class and numerical methods
    parameters.txt   # Input parameters file
    freq.py         # Frequency analysis script
    vis.py          # Visualization script
 .gitignore          # Git ignore rules
 README.md           # This file
 output.txt          # Generated simulation results
```

## Features

### Vessel Class
- Encapsulates vessel properties (mass, area, flow coefficient, initial conditions)
- Implements getter/setter methods for all parameters
- Contains numerical solution methods

### Numerical Methods
1. **Euler Method**: First-order explicit method
2. **Runge-Kutta Method**: Fourth-order method for higher accuracy

### Input/Output
- Reads simulation parameters from `parameters.txt`
- Supports multiple vessel configurations
- Outputs time-series data to `output.txt`
- Handles comment lines (starting with `#`) in input files

## Usage

### Input File Format (`parameters.txt`)

```
# Simulation parameters
scheme T dt g rho Cd
# Vessel parameters (one per line)
m A K y0 v0 h0
```

**Parameters:**
- `scheme`: Numerical method (0 = Euler, 1 = Runge-Kutta)
- `T`: Total simulation time
- `dt`: Time step
- `g`: Gravitational acceleration
- `rho`: Water density
- `Cd`: Drag coefficient
- `m`: Vessel mass
- `A`: Cross-sectional area
- `K`: Flow rate coefficient
- `y0`: Initial displacement
- `v0`: Initial velocity
- `h0`: Initial water level

### Compilation and Execution

```bash
# Compile
g++ -o vessel_simulation src/main.cpp

# Run
./vessel_simulation
```

### Output Format

The program generates `output.txt` with columns:
```
Time | Vessel1_y | Vessel1_v | Vessel1_h | Vessel2_y | Vessel2_v | Vessel2_h | ...
```

## Mathematical Implementation

### Force Analysis
- **Weight**: `Fw = mg + ρgAh(t)`
- **Buoyancy**: `Fb = -ρgA(ȳ + y(t) - h(t))`
- **Drag**: `Fd = -(1/2)ρACₐ(dy/dt)|dy/dt|`

### Numerical Integration
The program converts the second-order ODE into a system of first-order ODEs:
- `dy/dt = v`
- `dv/dt = f_v(y, v, h, t)`
- `dh/dt = f_h(y, h, t)`

Both Euler and Runge-Kutta methods solve this coupled system simultaneously.

## Dependencies

- C++ Standard Library
- `<vector>`, `<fstream>`, `<iostream>`, `<string>`, `<sstream>`, `<cmath>`, `<iomanip>`

## Error Handling

- File I/O validation
- Parameter parsing error detection
- Graceful handling of malformed input files

## Example Usage

```cpp
// Create vessel with parameters
Vessel vessel(mass, area, flow_coeff, init_y, init_v, init_h);

// Solve using Runge-Kutta method
auto results = vessel.rungeKuttaMethod(total_time, time_step, gravity, density, drag_coeff);

// Access results: results[0] = y(t), results[1] = v(t), results[2] = h(t)
```

## Notes

- The simulation assumes downward direction as positive
- Water flow rate is independent of vessel motion
- Vessel walls have negligible thickness
- The system exhibits complex dynamics due to coupling between vessel motion and water level

## Author

Part of C++ Programming Coursework - Numerical Solution of Differential Equations
