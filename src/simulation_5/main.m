%% WINTER simulation
clear 
close all
clc

%% Load parameters
parameters;

rng(simulation.seed); % set the seed for the random number generator

%% Tuning the controllers
tuning_controllers

%% Simulink simulation
open_system(simulation.mdl);                    % open the model
in = Simulink.SimulationInput(simulation.mdl);  % set simulation parameters

% set_simulink_parameters(simulation.mdl, simulation.type);

tic
for i = 1:wind.WS_len
  stop_time = simulation.stop_time(i); % set the stop time 
  wind_speed = zeros(stop_time/wind.sample_t, 2);

  switch simulation.type
    case 1  % constant wind speed  
      wind_speed = run_single_constant_speed(wind.mean(i), wind_speed, stop_time);
    case {2, 6}  % wind speed ramp
      wind_speed = run_ramp(wind.ramp_WS_start(i), wind.ramp_WS_stop(i),      wind.ramp_time_start(i), wind.ramp_time_stop(i), wind_speed,      stop_time);
    case {3, 5} % generated time series
      wind_speed = run_generated_wind_series(wind.mean(i), wind.turbulence(i), wind_speed, stop_time, simulation.seed);
    case 4 % generator step response
      [in, stop_time, out_store] = run_generator_step();
    case 7 % with or without gain schdeuling
      wind_speed = run_generated_wind_series(wind.mean(i),wind.turbulence(i), wind_speed, stop_time, simulation.seed);
      [blade.kp_schedule, blade.ki_schedule] = run_blade_gains(blade, i);
    case 8 % with gain schdeuling or stall regulated
      wind_speed = run_ramp(wind.ramp_WS_start(i), wind.ramp_WS_stop(i),      wind.ramp_time_start(i), wind.ramp_time_stop(i), wind_speed,      stop_time);
      % wind_speed = run_generated_wind_series(wind.mean(i), wind.turbulence(i), wind_speed, stop_time, simulation.seed);
      [blade.kp_schedule, blade.ki_schedule] = run_stall_regulation(blade, i);
    case 9 % with or without gain schdeuling
      wind_speed = run_generated_wind_series(wind.mean(i), wind.turbulence(i), wind_speed, stop_time, simulation.seed);
      [blade] = run_pitch_dynamics(blade, i);
    case 10 % test with K_opt and K_opt_GE
      wind_speed = run_single_constant_speed(wind.mean(i), wind_speed, stop_time);
      % wind_speed = run_ramp(wind.ramp_WS_start(i), wind.ramp_WS_stop(i), wind.ramp_time_start(i), wind.ramp_time_stop(i), wind_speed, stop_time);
      [rotor, generator, blade, T_R0, omega_rated_GE] = run_Kopt_KoptGE(rho, lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, wind_speed(i, 2), rated_values, rated_values_P_GE_no_B, lookup_Pitch, lookup_pitch_P_GE, simulation, i);
    case 11 % wind speed ramp and change the value of K_opt
      wind_speed = run_ramp(wind.ramp_WS_start(i), wind.ramp_WS_stop(i),      wind.ramp_time_start(i), wind.ramp_time_stop(i), wind_speed,      stop_time);
      [generator] = run_K_opt_sensitivity(generator, i);

  end

  % Set the initial conditions
  if simulation.model ~= 3 && simulation.type ~= 10 
    [rotor, generator, blade, T_R0, x_est_initial, P_est_initial] = initial_conditions(rho,lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, wind_speed(1, 2), rated_values, lookup_Pitch, IMM);
  elseif simulation.model == 3 && simulation.type ~= 10
    [rotor, generator, blade, T_R0] = initial_conditions_GE(rho,lambda_vector, pitch_vector, lookup_cP, rotor, blade, generator, gearbox, wind_speed(1, 2), V0_rated, rated_values_P_GE, lookup_pitch_P_GE);
  end

  % Run the simulation
  out_store{i} = sim(in, 'ShowProgress','on'); % store the results of the simulation

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
    RMS_errors = post_process(out_store, wind, omega_rated, generator, simulation, lookup_static_values, lookup_Pitch, lookup_P_GE);
  end