%% Plot of the results
% This file is used for print the results of the simulation

%% Static simulation
% Comparison between the results from the DTU report and the static
% simulation ones
figure(1)
plot(WS, omega_r_store*30/pi, 'o')
hold on
plot(reference(:, 1), reference(:, 3))
xlabel('Wind speed [m/s]')
ylabel('\omega_r [rpm]')
legend('simulation', 'reference', 'location', 'best')
grid on

figure(2)
plot(WS, rad2deg(pitch_store), 'o')
hold on
plot(reference(:, 1), reference(:, 2))
xlabel('Wind speed [m/s]')
ylabel('\theta [°]')
legend('simulation', 'reference', 'location', 'best')
grid on