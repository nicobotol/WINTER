function [wind_speed] = run_single_constant_speed(WS)
% This function runs the simulation for a constant windspeed provided as
% input
% WS -> wind speed [m/s]

parameters;

% Define the wind speed vector
wind_speed        = zeros(simulation.stop_time/wind.sample_t, 2);
wind_speed(:, 1)  = [wind.sample_t:wind.sample_t:...
  simulation.stop_time]'; % [s]
wind_speed(:, 2)  = WS;

end