%% WINTER simulation
clear 
close all
clc

%% Load parameters
parameters;
% parameters_sim2;

%% Load PMSM transfer functions and design the controller
[Yiq, Gc, Riq, GR] = PMSM_TF_pid(generator.design, generator.bode_plot);

%% Simulink simulation
open_system(simulation.mdl);                    % open the model
in = Simulink.SimulationInput(simulation.mdl);  % set simulation parameters
set_param(strcat(simulation.mdl,'/Gain1'),'SampleTime','simulation.time_step_H');
set_param(strcat(simulation.mdl,'/Gain'),'SampleTime','simulation.time_step_H');
set_param(strcat(simulation.mdl,'/Somma'),'SampleTime','simulation.time_step_L');
set_param(strcat(simulation.mdl,'/Step'),'SampleTime','simulation.time_step_H');
% set_param(strcat(simulation.mdl,'/sum1'),'SampleTime','simulation.time_step_H');
% set_param(strcat(simulation.mdl,'/sum2'),'SampleTime','simulation.time_step_H');
out_store = cell(wind.WS_len);

tic
for i = 1:wind.WS_len
  stop_time = simulation.stop_time(i); % set the stop time 
  wind_speed = zeros(stop_time/wind.sample_t, 2);
  
  switch simulation.type
    case 1  % constant wind speed  
      wind_speed = run_single_constant_speed(wind.mean(i), wind_speed, ...
        stop_time);
    case 2  % wind speed ramp
      wind_speed = run_ramp(wind.ramp_WS_start, wind.ramp_WS_stop, ...
        wind.ramp_time_start(i), wind.ramp_time_stop(i), wind_speed, ...
        stop_time);
    case 3  % generated time series
      wind_speed = run_generated_wind_series(wind.mean(i), ...
        wind.turbulence(i), wind_speed, stop_time);
  end

  % Run the simulation
  out = sim(in, 'ShowProgress','on'); 
  out_store{i} = out; % store the results of the simulation
end
toc

%% Plot the results
% plots