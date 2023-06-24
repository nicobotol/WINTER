function [E_R, E_G, E_GE] = ...
  power_check(out_cell, I, B, n_series, line_width, date_fig, plot_name,...
  simulation, font_size, path_images)
% This fucntion checks the sum of the power of the system:
% the generator power must be equal to the sum of the genrator, a term proportional to the inertia and one proportional to the damping
% I -> equivalent inertia (generator side) [kg*m^2]
% B -> equivalent damping (generator side) [Ns/m]
% n_series -> of wind number of series

P_R = cell(1, n_series);
P_G = cell(1, n_series);
time = cell(1, n_series);
P_inertia = cell(1, n_series);
P_damping = cell(1, n_series);
P_sum = cell(1, n_series);
P_GE = cell(1, n_series);

for i = 1:n_series
  P_R{i} = out_cell{i}.P_R.Data;        % rotor power [W]
  P_G{i} = out_cell{i}.P_G.Data;        % generator input power [W]
  P_GE{i} = out_cell{i}.P_GE.Data;      % generator output power [W]
  omega_R = out_cell{i}.omega_R.Data;   % rotor speed [rad/s]
  time{i} = out_cell{i}.tout;           % time [s]
  domega_R = diff(omega_R)./diff(time{i}); % rotor speed derivative [rad/s^2]
  P_inertia{i} = I.*omega_R(1:end - 1).*domega_R;  % inertia power [W]
  P_damping{i} = B.*omega_R(1:end - 1).^2;         % damping power [W]
  P_sum{i} = P_G{i}(1:end-1) + P_damping{i};          % sum of the power [W]

  % compute the energy
  E_R{i} = trapz(time{i}, P_R{i});         % rotor energy [J]
  E_G{i} = trapz(time{i}, P_G{i});         % generator energy [J]
  E_GE{i} = trapz(time{i}, P_GE{i});       % generator energy [J]

end

fig = figure('Color', 'w');
for i = 1:n_series
  yyaxis left
  plot(time{i}, P_R{i}/1e6, 'DisplayName', '$P_{rotor}$', 'LineWidth', ...
    line_width);
  hold on
  plot(time{i}, P_G{i}/1e6, '--', 'DisplayName', '$P_{generator}$', ...
    'LineWidth', line_width);
  plot(time{i}(1:end - 1), P_sum{i}/1e6, ':', 'DisplayName', ...
    '$P_{generator} + P_{damping}$', 'LineWidth', line_width);
  ylabel('Power [MW]');

  yyaxis right
  plot(time{i}(10:end - 1), P_inertia{i}(10:end)/1e6, '-', 'DisplayName', ...
    '$P_{inertia}$', 'LineWidth', line_width);
  plot(time{i}(10:end - 1), P_damping{i}(10:end)/1e6, '--', ...
    'DisplayName', '$P_{damping}$', 'LineWidth', line_width);
  xlabel('Time [s]');
  ylabel('Power lost and stored [MW]');
  legend('Location', 'south', 'NumColumns', 3);
  title('Comparison of the powers');
  grid on
  set(gca, 'FontSize', font_size)
end
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date_fig, plot_name,'.eps');
  export_figure(fig, strcat(date_fig, plot_name, '.eps'), path_images);
%   export_IEEE_figure(strcat(date_fig, plot_name, '.eps'), path_images); 
end
  
end