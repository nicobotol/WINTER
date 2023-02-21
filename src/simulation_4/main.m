%% WINTER simulation
clear 
close all
clc

%% Load parameters
parameters;

%% Load PMSM transfer functions
[Yiq, Gc, Riq] = PMSM_TF(generator.design, generator.bode_plot);

%% Simulink simulation
% Run a simulation with wind blowing at constant speed, for all the
% windspeeds between cut-in and cut-out
mdl = 'winter_simulink_no_delay';                        % model's name
open_system(mdl);                               % open the model
in = Simulink.SimulationInput(mdl);             % set simulation parameters
wind_speed = zeros(simulation.stop_time/wind.sample_t, 2);
out_store = cell(wind.WS_len);

for i = 1:wind.WS_len
  switch simulation.type
    case 1  % constant wind speed  
      wind_speed = run_single_constant_speed(wind.mean(i));
    case 2  % wind speed ramp
      wind_speed = run_ramp(wind.ramp_WS_start, wind.ramp_WS_stop, ...
        wind.ramp_time_start(i), wind.ramp_time_stop(i));
    case 3  % generated time series
      wind_speed = run_generated_wind_series(wind.mean(i), ...
        wind.turbulence(i));
  end

  % Run the simulation
  out = sim(in, 'ShowProgress','on'); 
  out_store{i} = out; % store the results of the simulation
end

%% Plot the results
% plots