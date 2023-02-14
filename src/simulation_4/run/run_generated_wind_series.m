function [wind_speed] = run_generated_wind_series(wind_mean, ...
  wind_turbulence)
% This function generates s wind series based on the parameters setted in
% parameters.m and runs the simulation for it

parameters
wind_speed = zeros(simulation.stop_time/wind.sample_t, 2);

% l parameter based on the height
if wind.height <= 30 
  l = 20;
else
  l = 600;
end

delta_ts = wind.sample_t; % sampling time [s]
N = simulation.stop_time/delta_ts;           % number of samples

I = wind_turbulence/wind_mean;  % turbulence intensity

wind_speed(:, 1) = delta_ts*[1:1:N]; % [s] time  
cos_v = zeros(N, 1);  % vector of cosines
p_sum = 0;            % partial sum
for n = 1:N/2
  f_n = n/simulation.stop_time;                           % frequency [Hz]
  PSD = I^2*wind_mean*l/(1 + 1.5*f_n*l/wind_mean)^(5/3);  % PSD
  phi_n = rand(1)*2*pi;                                   % phase [rad]
  cos_v = cos(2*pi*f_n.*wind_speed(:, 1) - phi_n);        % cosines
  p_sum = p_sum + sqrt(2*PSD/simulation.stop_time)*cos_v; % partial sum
end

wind_speed(:, 2) = wind_mean + p_sum;  % add the mean to the WS [m/s]

% Rescale the std
wind_speed(:, 2) = wind_turbulence*wind_speed(:, 2) + ...
  (1 - wind_turbulence)*wind_mean;

% close all
% figure()
% plot(wind_speed(:, 1), wind_speed(:, 2))
% xline(mean(wind_speed(:, 2)))
% grid on
end