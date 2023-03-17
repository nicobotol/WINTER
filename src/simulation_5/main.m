%% WINTER simulation
clear 
close all
clc

%% Load parameters
parameters;
% parameters_sim2;

rng(simulation.seed); % set the seed for the random number generator

%% Load PMSM transfer functions and design the controller
[Yiq, Gc, Riq, GR] = PMSM_TF_pid(generator.design, generator.bode_plot);

%% Simulink simulation
open_system(simulation.mdl);                    % open the model
in = Simulink.SimulationInput(simulation.mdl);  % set simulation parameters

set_simulink_parameters(simulation.mdl, simulation.type);

tic
for i = 1:wind.WS_len
  stop_time = simulation.stop_time(i); % set the stop time 
  wind_speed = zeros(stop_time/wind.sample_t, 2);
  
  switch simulation.type
    case 1  % constant wind speed  
      wind_speed = run_single_constant_speed(wind.mean(i), wind_speed, ...
        stop_time);
    case {2, 6}  % wind speed ramp
      wind_speed = run_ramp(wind.ramp_WS_start, wind.ramp_WS_stop, ...
        wind.ramp_time_start(i), wind.ramp_time_stop(i), wind_speed, ...
        stop_time);
    case {3, 5} % generated time series
      wind_speed = run_generated_wind_series(wind.mean(i), ...
        wind.turbulence(i), wind_speed, stop_time);
    case 4 % generator step response
      [stop_time] = run_generator_step();
  end

  % Run the simulation
  out = sim(in, 'ShowProgress','on'); 
  out_store{i} = out; % store the results of the simulation

%   % Uncomment all the blocks
%   uncomment_all(simulation.mdl)
end
toc

%% Check if aby if the simulations ended with error
simulation_fails = uint8(0); 
for i=1:wind.WS_len
  if out_store{1}.SimulationMetadata.ModelInfo.StopTime < simulation.stop_time(i);
    simulation_fails = 1;
  end
end
%% Plot the results
if simulation_fails == 0
  plots
else
  fprintf('Smulation ended with error\n');
end