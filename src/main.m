clear 
close all
clc

%% Load parameters
parameters
R = 9.5; % rotor diameter (m)
rho = 1.225; % air density (kg/m^3)
sincronous_velocity = 5; % (rad/s)
MG_coefficient = 2.55e5; % generator curve coefficient
I = 15e3; % rotor iniertia [kgm^2]
delta_t = 0.005; % discretization time (s)
% load the look-up table for the power coefficient
cP_vs_lambda = load('Cp_Lambda.txt');
% load the wind time history
wind_speed = load('usim.dat');


omega_r = 0;
%% Simulink simulation
mdl = 'winter_simulink'; 
open_system(mdl);
in = Simulink.SimulationInput(mdl);
out = sim(in, 'ShowProgress','on'); % performing the simulation