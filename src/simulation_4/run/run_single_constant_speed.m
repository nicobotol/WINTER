function [wind_speed] = run_single_constant_speed(WS, wind_speed, ...
  stop_time)
% This function runs the simulation for a constant windspeed provided as
% input
% WS -> wind speed [m/s]

parameters;

% Define the wind speed vector
wind_speed(:, 1)  = [wind.sample_t:wind.sample_t:stop_time]'; % [s]
wind_speed(:, 2)  = WS;

end