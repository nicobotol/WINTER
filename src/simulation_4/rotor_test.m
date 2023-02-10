%% Test the rotor performance

clear 
close all
clc

%% Load parameters
addpath("lookup"); % add the lookup tables path
parameters;

%% Simulink simulation
% Run a simulation with wind blowing at constant speed, for all the
% windspeeds between cut-in and cut-out
mdl = 'rotor_simulink';                         % model's name
open_system(mdl);                               % open the model
set_param(mdl, 'StopTime', num2str(stop_time)); % set simulation time
in = Simulink.SimulationInput(mdl);             % set simulation parameters

K_opt =  0.5*(lambda_opt)^(-3)*rho*pi*rotor.R^5*cp_max; % constant for the generator

%% Single simulation
wind_speed = 15;
WS = wind_speed;
out = sim(in, 'ShowProgress','on');             % run the simulation
omega_r_store = mean(out.omega_r.Data(end - 20:end));
pitch_store = mean(out.pitch.Data(end - 20:end));
P_r_store = mean(out.P_r.Data(end - 20:end));

%% Multiple simulation

% WS = V0_cut_in:2:V0_cut_out;                    % range of ws to test [m/s]
% WS = V0_cut_in:1:11;
% WS = 12:1:V0_cut_out; 
% WS_length = length(WS);
% omega_r_store = zeros(1, WS_length);
% pitch_store = zeros(1, WS_length);
% P_r_store = zeros(1, WS_length);
% Pr_ref = zeros(1, WS_length);
% for i = 1:WS_length
%   wind_speed = WS(i);
%   out = sim(in, 'ShowProgress','on');             % run the simulation
%   omega_r_store(i) = mean(out.omega_r.Data(end - 20:end));
%   pitch_store(i) = mean(out.pitch.Data(end - 80:end));
%   P_r_store(i) = mean(out.P_r.Data(end - 20:end));  
% end

%% Plot the results
plots