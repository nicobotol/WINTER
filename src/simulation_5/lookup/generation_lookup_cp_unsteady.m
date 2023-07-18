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
lambda_s1 = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); 
lambda_s2 = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); 


for i = 1:velocity_item_3v
  V0 = velocity_vector_3v(i);
  lambda_s1(:, :, i) = V0/2*(1 + sqrt(1 - lookup_cT_3v(:, :, i)));
  lambda_s2(:, :, i) = V0/2*(1 - sqrt(1 - lookup_cT_3v(:, :, i)));
  pause
end

%%
% Save the results 
% save('lookup_cP_3v_u_theta_lambda.mat', 'lookup_cP_3v_u');
% save('lookup_cT_3v_u_theta_lambda.mat', 'lookup_cT_3v_u');