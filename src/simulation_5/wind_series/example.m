% Test the wind seris generation
close all
clear
clc

parameters;

V10 = 10;
V10_std = 1;
fs = 50;
h = 150;
T = 300;
[u, t, PSD] = wind_series(V10, V10_std, fs, h, T);

%% Plot the time serie
fig = figure('Color','w');
hold on
plot(t, u, 'LineWidth', line_width);
hold on
yline(V10, 'r', 'LineWidth',1.2*line_width);
yline(V10+3*V10_std , '--r', '$V_{10}+3\sigma_{V_{10}}$', 'LineWidth',1.2*line_width,'FontSize', font_size, 'Interpreter','latex');
yline(V10-3*V10_std , '--r', '$V_{10}-3\sigma_{V_{10}}$', 'LineWidth',1.2*line_width,'FontSize', font_size, 'Interpreter','latex');
xlabel('Time [s]','interpreter','latex')
ylabel('Wind speed [m/s]','interpreter','latex')
title('Wind speed series generation', 'Interpreter','latex')
grid on
box on
ylim([6, 14])
set(gca, 'FontSize', font_size)

fig_name = strcat(path_images,'\', date_fig, 'wind_generation','.eps');
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'wind_generation', '.eps'), path_images);
end
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 

%% Plot the PSD
N = length(u);
udft = fft(u);
udft = udft(1:N/2+1);
psdx = (1/(fs*N)) * abs(udft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/N:fs/2;

fig = figure('Color', 'w');
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

fig_name = strcat(path_images,'\', date_fig, 'wind_generation_PSD','.eps');
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'wind_generation_PSD', '.eps'), path_images);
end
