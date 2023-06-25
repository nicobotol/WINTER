%% WINTER simulation
clear 
close all
clc

%% Load parameters
parameters;
% parameters_sim2;

rng(simulation.seed); % set the seed for the random number generator

%% Tuning the controllers
tuning_controllers

%% Run the simualtions
tic
for i=1:wind.WS_len
  if i==1
    simulation.mdl = 'winter_simulink_with_PC';
  else 
    simulation.mdl = 'winter_simulink_with_PC_generator_control';
  end

  open_system(simulation.mdl);                    % open the model
  in = Simulink.SimulationInput(simulation.mdl);  % set simulation parameters
  set_simulink_parameters(simulation.mdl, simulation.type);
  stop_time = simulation.stop_time(i); % set the stop time 
  wind_speed = zeros(stop_time/wind.sample_t, 2);
  wind_speed = run_ramp(wind.ramp_WS_start, wind.ramp_WS_stop, ...
    wind.ramp_time_start(i), wind.ramp_time_stop(i), wind_speed, ...
    stop_time);
  if i==1
    [rotor, generator, blade, T_R0] = initial_conditions(rho,lambda_vector, pitch_vector, lookup_cP,rotor, blade, ...
      generator, gearbox, wind_speed(1, 2), rated_values, lookup_Pitch);
  else
    [rotor, generator, blade, T_R0] = initial_conditions_GE(rho, ...
      lambda_vector, pitch_vector, lookup_cP, rotor, blade, ...
      generator, gearbox,wind_speed(1, 2), V0_rated, rated_values_P_GE, lookup_pitch_P_GE);
  end
  out = sim(in, 'ShowProgress','on'); 
  out_store{i} = out;
end 

toc

%% Check if any of the simulations ended with error
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
  fprintf('Simulation ended with error\n');
end

%% Post processing
if simulation.type ~= 4
  RMS_errors = post_process(out_store, wind, omega_rated, generator, simulation);
end