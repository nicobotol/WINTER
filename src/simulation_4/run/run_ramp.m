function [wind_speed] = run_ramp(WS_start, WS_stop, time_start, ...
  time_stop, wind_speed, stop_time)
% This function runs the simulation for a wind speed ramp
% WS_start -> wind speed at the start of the ramp [m/s]
% WS_stop -> wind speed at the stop of the ramp [m/s]
% time_start -> time when the ramp starts [s]
% time_stop -> time when the ramp stops [s]
% ts -> wind serie sample time [s]

parameters

% Check if the time instants are correct
if time_start <= time_stop && time_stop <= stop_time

else
  error('Incorrect time steps for the ramp')
end


% Define the wind speed vector
wind_speed(:, 1)  = [wind.sample_t:wind.sample_t:stop_time]'; % [s]
N_start = time_start/wind.sample_t;
N_stop = time_stop/wind.sample_t;

wind_speed(1:N_start, 2)  = WS_start;
wind_speed(N_start+1:N_stop, 2) = WS_start + (WS_stop - WS_start)/...
  (N_stop - N_start).*([0:1:N_stop-N_start-1]');
wind_speed(N_stop+1:end, 2) = WS_stop;

end