%% WINTER simulation

clear 
close all
clc

%% Load parameters
addpath("lookup"); % add the lookup tables path
parameters;

%% Load PMSM transfer functions

PMSM_TF;

%% Simulink simulation
% Run a simulation with wind blowing at constant speed, for all the
% windspeeds between cut-in and cut-out
mdl = 'winter_simulink';                        % model's name
open_system(mdl);                               % open the model
set_param(mdl, 'StopTime', num2str(stop_time)); % set simulation time
in = Simulink.SimulationInput(mdl);             % set simulation parameters

% WS = V0_cut_in:1:V0_cut_out;                  % range of ws to test [m/s]
WS = 10;
WS_length = length(WS);
omega_r_store = zeros(1, WS_length);
pitch_store = zeros(1, WS_length);
P_r_store = zeros(1, WS_length);
for i = 1:WS_length
  wind_speed = WS(i);
  out = sim(in, 'ShowProgress','on');             % run the simulation
  omega_r_store(i) = mean(out.omega_r.Data(end - 80:end));
  pitch_store(i) = mean(out.pitch.Data(end - 80:end));
  P_r_store(i) = mean(out.P_r.Data(end - 80:end));
end

%% Plot the results
plots