%% unsteady cP as function of pitch angle and Tip Speed ratio (TSR)
% This file is aimed to create a mesh of unsteady cP and cT as function of the pitch angle and the tip speed ratio given the corresponding static maps.

clear 
close all
clc

% add path for the functions for the aerodynamic
addpath("..\")
addpath("..\aerodynamic_functions")
addpath("..\aerodynamic_functions\airfoil_data")


% load the parameters
parameters

lookup_cP_3v_u = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); % matrix to store cP
lookup_cT_3v_u = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); % matrix to store cT
% lookup_cP(theta, lambda)

for i = 1:velocity_item_3v
  V0 = velocity_vector_3v(i);
  % T = 1/2*rho*rotor.A*lookup_cT_3v(:, :, i).*V0.^2;         % [N] thrust
  % lambda_s_1 = V0 + sqrt(V0.^2 - 2*T/(rho*pi*rotor.R^2));  % [m/s] speed at the streamtube section
  % lambda_s_2 = V0 - sqrt(V0.^2 - 2*T/(rho*pi*rotor.R^2)); % [m/s] speed at the streamtube section

  lambda_s1 = V0/2*(1 + sqrt(1 - lookup_cT(:, :, i)));
  lambda_s1 = V0/2*(1 - sqrt(1 - lookup_cT(:, :, i)));
end

%%
% Save the results 
% save('lookup_cP_3v_u_theta_lambda.mat', 'lookup_cP_3v_u');
% save('lookup_cT_3v_u_theta_lambda.mat', 'lookup_cT_3v_u');