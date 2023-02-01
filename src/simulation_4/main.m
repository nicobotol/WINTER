%% WINTER simulation

clear 
close all
clc

%% Load parameters
addpath("lookup"); % add the lookup tables
parameters;

%% Simulink simulation
mdl = 'winter_simulink';                        % model's name
open_system(mdl);                               % open the model
set_param(mdl, 'StopTime', num2str(stop_time)); % set simulation time
in = Simulink.SimulationInput(mdl);             % set simulation parameters
% out = sim(in, 'ShowProgress','on');             % run the simulation

WS = V0_cut_in:1:V0_cut_out;
WS_length = length(WS);
omega_r_store = zeros(1, WS_length);
pitch_store = zeros(1, WS_length);
for i = 1:WS_length
  wind_speed = WS(i);
  out = sim(in, 'ShowProgress','on');             % run the simulation
  omega_r_store(i) = mean(out.omega_r.Data(end - 80:end));
  pitch_store(i) = mean(out.pitch.Data(end - 80:end));
end

%% Plot the results
plots