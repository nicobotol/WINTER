clear 
close all
clc

font_size = 25;             % fontsize for plots
line_width = 2;             % line width for plots

set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  font_size)
set(0,'DefaultLegendFontSize', font_size)

%% Load parmeters
parameters

rng(simulation.seed);

V10 = 10;           % mean wind speed [m/s]
V10_std = 0.1*V10;  % std deviaition [m/s]
T = 300;            % simulation time [s]
fs = 50;            % sampling frequency [Hz]
height = 150;

[u, t, PSD] = wind_series(V10, V10_std, fs, height, T);
% [u1, ~, ~] = wind_series(V10, V10_std, fs, height, T);
% [u2, ~, ~] = wind_series(V10, V10_std, fs, height, T);

fig1 = figure('Color','w');
plot(t, u, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('WS [m/s]')
% hold on
% plot(t, u1)
% plot(t, u2)
ylim([V10-3.5*V10_std, V10+3.5*V10_std]);
yline(V10, 'r', '$V_{10}$', 'LineWidth', line_width , 'Interpreter','latex', 'FontSize',font_size)
yline(V10 + 3*V10_std, '--r', '$V_{10} + 3\sigma_{V_{10}}$', 'LineWidth', line_width, 'Interpreter','latex', 'FontSize',font_size)
yline(V10 - 3*V10_std, '--r', '$V_{10} - 3\sigma_{V_{10}}$', 'LineWidth', line_width, 'Interpreter','latex', 'FontSize',font_size, 'LabelVerticalAlignment','bottom')
title('Wind speed simulation')
grid on
path = 'C:\Users\Niccolò\Documents\UNIVERSITA\TESI_MAGISTRALE\report\images\';
export_fig('fig1',[path, 'wind_generation.svg']);

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

fig2 = figure('Color','w');
semilogx([1:1:N]/T, PSD, 'DisplayName', 'PSD from definition', ...
  'LineWidth', line_width);
hold on
semilogx(freq, psdx, '--', 'DisplayName', 'PSD from data', ...
  'LineWidth', line_width);
grid on
ylabel('PSD [$\frac{m^2}{s}$]', 'Interpreter','latex')
xlabel('Frequency [Hz]',  'Interpreter','latex')
title('Wind speed PSD')
legend('Location', 'northeast')
export_fig('fig2',[path, 'wind_generation_PSD.svg']);

psd_integral = trapz([1:1:N]/T, PSD);
fprintf('PSD integral: %6.3f\n', psd_integral);
