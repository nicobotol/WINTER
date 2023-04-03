% This script computes the inertia moments of the blade and the rotors 
% wrt rotations around the rotor axis. 
% Axial positions and distributed masses are loaded from tables

clear 
close all
clc

%% DTU 10 MW
addpath("airfoil_data");
data_10 = readmatrix('inertia_data_DTU10MW.txt'); % read data 
r_10 = data_10(:, 1); % radial position [m]
rho_10 = data_10(:, 2); % linear mass [kg/m]

[I_blade_10, I_rotor_10, m_blade_10, m_rotor_10] = ...
  rotor_inertia(r_10, rho_10, 'DTU 10 MW'); % results of the DTU 10 MW

%% NREL 5 MW
addpath("airfoil_data");
data_5 = readmatrix('inertia_data_NREL5MW.txt'); % read data 
r_5 = data_5(:, 1); % radial position [m]
rho_5 = data_5(:, 5); % linear mass [kg/m]

[I_blade_5, I_rotor_5, m_blade_5, m_rotor_5] = ...
  rotor_inertia(r_5, rho_5, 'NREL 5 MW'); % results of the DTU 10 MW
