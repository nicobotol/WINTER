clear 
close all
clc

%% Load parmeters
addpath('C:\Users\Niccolò\Documents\UNIVERSITA\TESI_MAGISTRALE\src\simulation_4');
parameters

rng(simulation.seed);

[u, t] = wind_series(wind.mean, 0.5, wind.sample_f, ...
  wind.height, simulation.stop_time);
% [u1, ~] = wind_series(wind.mean, wind.turbulence, wind.sample_f, ...
%   wind.height, simulation.stop_time);
% [u2, ~] = wind_series(wind.mean, wind.turbulence, wind.sample_f, ...
%   wind.height, simulation.stop_time);

%%

figure()
plot(t, u)
% hold on
% plot(t, u1)
% plot(t, u2)
yline(wind.mean, '--r')
grid on

