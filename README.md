# Computational Chaos: Controlling Chaos in Networks of Oscillators

This repository contains MATLAB implementations for studying and controlling chaos in networks of coupled oscillators, based on research in nonlinear dynamics and chaos control theory.

## Overview

The project implements various algorithms and methods for:
- Finding Unstable Periodic Orbits (UPOs) in chaotic systems
- Implementing the OGY (Ott-Grebogi-Yorke) method for chaos control
- Analyzing network dynamics of coupled oscillators
- Visualizing phase space trajectories and Poincaré sections

## Project Structure

### Core MATLAB Files

- **`device.m`** - Main implementation of the OGY chaos control method for network oscillators
- **`Get_upo.m`** - Algorithm for finding Unstable Periodic Orbit candidates using Poincaré sections
- **`ogy.m`** - Implementation of the OGY method for chaos control
- **`Refine_upo.m`** - Refinement algorithm for improving UPO approximations
- **`refine_upo2.m`** - Alternative UPO refinement method
- **`ploting_orbit.m`** - Visualization utilities for plotting orbits and trajectories
- **`plot_re_w.m`** - Plotting functions for real parts of oscillator states
- **`sample.m`** - Sample code and examples

### Data Files

- **`C3_upo.mat`** - Pre-computed UPO data for 3x3 oscillator network
- **`C2_refined_upo.mat`** - Refined UPO data for 2x2 oscillator network
- **`stabilized_upo.mat`** - Stabilized UPO states (large file: 187MB)
- **`upo_t=46.mat`** - UPO data at specific time point
- **`refined_upo_state.mat`** - Refined UPO state data
- **`upo_states.mat`** - UPO state information

### Visualization Files

- **`stabilized_upo.png`** - Visualization of stabilized UPO
- **`stabilized_upo_ps.png`** - Phase space visualization
- **`md_plus.png`** - Plus pattern visualization
- **`md_cross.png`** - Cross pattern visualization
- **`untitled1.png`** - Additional visualization

### Documentation

- **`Comp_with_chaos.pdf`** - Main project documentation
- **`Controlling chaos in a network of oscillators - PhysRevE.48.945.pdf`** - Reference paper

## Mathematical Background

### Network Dynamics

The system models a network of N×N coupled oscillators described by:

```
Ẋ = f(X, Y, W, I)
Ẏ = g(X, Y, W, I)
```

Where:
- X, Y are the oscillator state matrices
- W represents coupling weights
- I is the input pattern matrix
- f, g are nonlinear functions defining the oscillator dynamics

### OGY Method

The Ott-Grebogi-Yorke method is implemented for controlling chaos by:
1. Finding unstable periodic orbits in the chaotic attractor
2. Computing the stable and unstable manifolds
3. Applying small perturbations to stabilize the system around the UPO

### Poincaré Sections

The code uses Poincaré sections to:
- Detect when trajectories cross a specific hyperplane
- Find periodic orbit candidates
- Analyze the structure of the chaotic attractor

## Usage

### Prerequisites

- MATLAB (R2018b or later recommended)
- MATLAB's ODE solver toolbox
- Sufficient memory for large data files (stabilized_upo.mat is 187MB)

### Basic Usage

1. **Find UPOs**: Run `Get_upo.m` to find unstable periodic orbit candidates
2. **Refine UPOs**: Use `Refine_upo.m` or `refine_upo2.m` to improve UPO approximations
3. **Control Chaos**: Execute `ogy.m` to implement chaos control
4. **Visualize Results**: Use plotting functions to analyze results

### Example Workflow

```matlab
% 1. Initialize system and find UPOs
Get_upo;

% 2. Refine the found UPOs
Refine_upo;

% 3. Apply chaos control
ogy;

% 4. Visualize results
ploting_orbit;
```

## Key Parameters

- **N**: Network size (default: 9×9 oscillators)
- **D**: Coupling strength (default: 1.3)
- **α, β**: Oscillator parameters (default: α=-10, β=2)
- **Threshold**: Distance threshold for UPO detection (default: 0.1)

## Research Applications

This code is designed for research in:
- Nonlinear dynamics and chaos theory
- Network synchronization
- Chaos control in coupled oscillator systems
- Pattern formation in complex networks

## References

- Ott, E., Grebogi, C., & Yorke, J. A. (1990). Controlling chaos. Physical Review Letters, 64(11), 1196.
- "Controlling chaos in a network of oscillators" - Physical Review E, 48(2), 945

## Contributing

This is an academic research project. For questions or contributions, please refer to the main documentation in `Comp_with_chaos.pdf`.

## License

This project is for academic research purposes. Please cite the original research papers when using this code in publications.

---
