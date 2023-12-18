# WINd Turbine EmulatoR (WINTER)

This repository contains the work done for the master thesis with title: **Development of a simulator for control strategy testing in wind turbine generators**

The scope of the thesis is twofold: developing a Wind Turbine (WT) emulator and use it to test some control laws. <br>
The target Wind Turbine (WT) is the DTU 10 MW. 

## Capabilities of the emulator
The simulator is intended to study static and dynamic aspects of a WT. In particular, one can test the standard MPPT control below rated wind speed (WS) and the pitch regulation above rated ones. The MPPT is firstly implemented for maximizing the power extracted from the resource, then the generator is included in the model and a new control law maximizing the generator electrical power is included, and finally the Interactive Multiple Model (IMM) is used for takeing into account the variabilities in the physical parameters. <br>
Different input WSs may be applied as input, in particular constant, ramp, and synthetic WS defined starting from a pre-defined Power Spectral Density (PSD).  

## Structure of the repository
*report*: contains the text used for the report and rhe presentation <br>
*references*: contains the reference used for the report <br>
*src/simulation_5*: contains the code used in the simulations<br>
*src/simulation_5/aerodynamic_functions*: Blade Element Momentum (BEM) algorithm for computing the aerodynamic of the blades. No effects of mooring system, nor hydrodynamic implications are taken into account. <br>
*src/simulation_5/aerodynamic_functions_frozen_wake*: *Frozen Wake BEM* for computing the variable gains necessary for controlling the blade pitching by means of the gain scheduling approach. <br>
*src/simulation_5/controllers*: functions necessary for tuning the voltage controller for the generator machine. <br>
*src/simulation_5/lookup*: lookup tables used for the extracting the static data during the dynamic simulations are reported in 
*src/simulation_5/plot*: scripts for plotting the results
*src/simulation_5/run*: scripts for running the Simulink simulations
*src/simulation_5/simulink*: Simulink files
*src/simulation_5/utilities*: function used for performing post-process analysis and other things
*src/simulation_5/wind_series*: function used for generating the wind series.

## Workflow
### 1. Define the parameters
The parameters that have to be chosen in *parameters.m* are:
1. Simulink model (Control based on the MPPT on the rotor side, MPPT on the generator electrical power, use of the IMM)
2. Wind speed (constant value, ramp, synthetic) time for starting/stopping the WS and value of the WS
3. Simulation time
4. Choice if use or not the IMM

### 2. Simulate
After the definition of the parameters one can start the simulation by running the *main.m*.

### 3. Post process and plots
After the simulation, the post processing of the data and the plots are done. One can choose the type of plots that have to be done.