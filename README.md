# A Benchmark Simulator for Advanced Control of Ethanol Steam Reforming
This repository contains the benchmark simulator for a staged-separation membrane reactor (SSMR) designed for pure hydrogen production. The simulator implements a nonlinear dynamic model of a two-stage process: ethanol steam reforming (ESR) and hydrogen separation. It is intended as a benchmark for testing control and estimation strategies.

### Installation

1. Download the ZIP folder from the GitHub repository 
2. Extract the contents to a folder of your choice
3. Launch MATLAB R2020a or later (recommended). Additional toolboxes are not necessary
4. Set the extracted folder as the current folder in MATLAB
5. Open the main simulation script: SSMR_simulation.m 

### Instructions

1. Specify simulation conditions in SSMR_simulator/SSMR_simulation.m
   * Select the normal operating conditions: Mode 1, 2 or 3 (see Table 2 in the paper)
   * Select the disturbance scenario (see subsection 4.2 in the paper)
   * Specify the time at which the disturbance is to be applied (if any)
   * Select the initial conditions (see more details in the supplementary material)
   * Select the overall simulation time
   * Select the set-point profile (see Fig. 4 in the paper)
   * Select the simulation type (open loop or with control)
2. Run SSMR_simulator/SSMR_simulation.m
3. The simulation results will be displayed in plots that appear on the screen
4. To modify or implement custom controllers, open control.m

### Run case studies

1.

2.
