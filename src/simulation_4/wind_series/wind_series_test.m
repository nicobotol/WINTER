clear 
close all
clc

%% Load parmeters
parameters

[u, t] = wind_series(wind.mean, wind.turbulence, wind.sample_f, ...
  wind.height, simulation.stop_time);
[u1, ~] = wind_series(wind.mean, wind.turbulence, wind.sample_f, ...
  wind.height, simulation.stop_time);
[u2, ~] = wind_series(wind.mean, wind.turbulence, wind.sample_f, ...
  wind.height, simulation.stop_time);

%%
close all
figure()
plot(t, u)
hold on
plot(t, u1)
plot(t, u2)
yline(wind.mean, '--r')
grid on

std(u)
std(u1)