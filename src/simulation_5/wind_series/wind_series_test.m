clear 
close all
clc

%% Load parmeters
parameters

rng(simulation.seed);

V10 = 10;           % mean wind speed [m/s]
V10_std = 0.1*V10;  % std deviaition [m/s]
T = 100;            % simulation time [s]
fs = 50;            % sampling frequency [Hz]
height = 150;

[u, t, PSD] = wind_series(V10, V10_std, fs, height, T);
[u1, ~, ~] = wind_series(V10, V10_std, fs, height, T);
[u2, ~, ~] = wind_series(V10, V10_std, fs, height, T);

figure()
plot(t, u)
xlabel('Time [s]')
ylabel('WS [m/s]')
hold on
plot(t, u1)
plot(t, u2)
yline(V10, '--r')
yline(V10 + 3*V10_std, '--r')
yline(V10 - 3*V10_std, '--r')
grid on

%% Analysis of the generated wind series
fprintf('Mean of WS: %7.4f\n', mean(u));
fprintf('Standard deviation of WS: %7.4f\n', std(u));

% Compute the PSD of the wind speed from the generated series
N = length(u);
udft = fft(u);
udft = udft(1:N/2+1);
psdx = (1/(fs*N)) * abs(udft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(u):fs/2;

figure()
semilogx([1:1:N]/T, PSD, 'DisplayName', 'PSD from definition', ...
  'LineWidth', line_width);
hold on
semilogx(freq, psdx, '--', 'DisplayName', 'PSD from data', ...
  'LineWidth', line_width);
grid on
ylabel('PSD [$\frac{m^2}{s}$]', 'Interpreter','latex')
xlabel('Frequency [Hz]',  'Interpreter','latex')
title('Wind PSD')
legend('Location', 'northeast')

psd_integral = trapz([1:1:N]/T, PSD);
fprintf('PSD integral: %6.3f\n', psd_integral);
