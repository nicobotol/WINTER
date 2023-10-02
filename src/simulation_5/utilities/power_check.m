function [E_R, E_G, E_GE] = power_check(out_cell, I, B, n_series, line_width, date_fig, plot_name, simulation, font_size, path_images)
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
   % Time from where start to print
   t_start = max(0, out_cell{i}.P_R.Time(end) - simulation.plot_time(i)); % [s]
   [~, s_start] = min(abs(out_cell{i}.P_R.Time - t_start)); % sample from where to start
  s_start = max(1, s_start);
  
  P_R{i} = out_cell{i}.P_R.Data;        % rotor power [W]
  P_G{i} = out_cell{i}.P_G.Data;        % generator input power [W]
  P_GE{i} = out_cell{i}.P_GE.Data;      % generator output power [W]
  P_GInductance{i} = out_cell{i}.P_GInductance.Data; % generator inductance power [W]
  P_GJoule{i} = out_cell{i}.P_GJoule.Data; % generator joule power [W]
  omega_R = out_cell{i}.omega_R.Data;   % rotor speed [rad/s]
  time{i} = out_cell{i}.tout;           % time [s]
  domega_R = min(1,max(-1,[0; diff(omega_R)./diff(time{i})])); % rotor speed derivative [rad/s^2]
  P_inertia{i} = out_cell{i}.P_inertia.Data; %I.*omega_R.*domega_R;  % inertia power [W]
  P_damping{i} = out_cell{i}.P_damping.Data;%B.*omega_R.^2;         % damping power [W]
  P_sum{i} = P_G{i} + P_damping{i};          % sum of the power [W]

  % compute the energy
  E_R{i} = trapz(time{i}, P_R{i});         % rotor energy [J]
  E_G{i} = trapz(time{i}, P_G{i});         % generator energy [J]
  E_GE{i} = trapz(time{i}, P_GE{i});       % generator energy [J]

end

fig = figure('Color', 'w');
for i = 1:n_series
  yyaxis left
  plot(time{i}(s_start:end), P_R{i}(s_start:end)/1e6, 'DisplayName', '$P_{R}$', 'LineWidth',line_width);
  hold on
  plot(time{i}(s_start:end), P_G{i}(s_start:end)/1e6, '--', 'DisplayName', '$P_{G}$','LineWidth', line_width);
  plot(time{i}(s_start:end), P_sum{i}(s_start:end)/1e6, ':', 'DisplayName','$P_{G} + P_{D}$', 'LineWidth', line_width);
  plot(time{i}(s_start:end), (P_G{i}(s_start:end) + P_damping{i}(s_start:end) + P_inertia{i}(s_start:end))/1e6, '--', 'DisplayName', '$P_{G} + P_{D} + P_{I}$', 'LineWidth', line_width, 'Color', color(3));
  ylabel('Power [MW]');

  yyaxis right
  plot(time{i}(s_start:end), P_inertia{i}(s_start:end)/1e6, '-', 'DisplayName','$P_{I}$', 'LineWidth', line_width);
  plot(time{i}(s_start:end), P_damping{i}(s_start:end)/1e6, '--','DisplayName', '$P_{D}$', 'LineWidth', line_width);
  xlabel('Time [s]');
  ylabel('Power lost and stored [MW]');
  legend('Location', 'southeast', 'NumColumns', 3);
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
