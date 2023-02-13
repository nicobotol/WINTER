%% Plot of the results
% This file is used for print the results of the simulation

%% Static simulation
% Comparison between the results from the DTU report and the static
% simulation ones
fig_omega = figure('Position', get(0, 'Screensize'));
plot(WS, omega_r_store*30/pi, 'o')
hold on
plot(reference(:, 1), reference(:, 3))
yline(omega_rated*30/pi, 'g')
xlabel('Wind speed [m/s]')
ylabel('\omega_r [rpm]')
legend('simulation', 'reference', 'location', 'best')
grid on

fig_power = figure('Position', get(0, 'Screensize'));
plot(WS, P_r_store, 'o')
hold on
plot(reference(:, 1), reference(:, 6)*1e3)
xlabel('Wind speed [m/s]')
ylabel('P [W]')
legend('simulation', 'reference', 'location', 'best')
title('Mechanical power')
grid on


fig_pitch = figure('Position', get(0, 'Screensize'));
plot(WS, rad2deg(pitch_store), 'o')
hold on
plot(reference(:, 1), reference(:, 2))
xlabel('Wind speed [m/s]')
ylabel('\theta [°]')
legend('simulation', 'reference', 'location', 'best')
grid on