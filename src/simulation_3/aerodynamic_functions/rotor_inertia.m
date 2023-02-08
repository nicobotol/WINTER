clear 
close all
clc

addpath("airfoil_data");
data = readmatrix('inertia_data.txt'); % read data 
r = data(:, 1); % radial position [m]
rho = data(:, 2); % linear mass [kg/m]

I_blade = 0;  % initialize blade inertia
for i = 2:length(r)
  mi = rho(i)*(r(i) - r(i - 1)); % mass of the section [kg]
  I_blade = I_blade + mi*r(i)^2; % inertia of the section [kgm^2]
end

I_rotor = 3*I_blade;
fprintf('The rotor inertia is I = %6.4e [kgm^2]\n', I_rotor);

