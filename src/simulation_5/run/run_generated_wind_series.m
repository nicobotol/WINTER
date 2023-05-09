function [wind_speed] = run_generated_wind_series(wind_mean, ...
  wind_turbulence, wind_speed, stop_time, seed)
% This function generates s wind series based on the parameters setted in
% parameters.m and runs the simulation for it

parameters
fs = wind.sample_f; % sampling frequency [Hz]
T = stop_time;      % total time [s]
rng(seed)
% Generate the winds
[wind_speed(:, 2), wind_speed(:, 1), PSD] = wind_series(wind_mean, wind_turbulence, fs, wind.height, stop_time);

% Compute the PSD of the wind speed from the generated series
if simulation.print_figure == 1
  N = length(wind_speed(:, 2));
  udft = fft(wind_speed(:, 2));
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
  
  plot_name = 'wind_EPS';
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);

end

end