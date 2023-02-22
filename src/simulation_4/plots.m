%% Plot of the results
% This file is used for print the results of the simulation
close all
%% Wind time series
fig_wind_TS = figure('Position', get(0, 'Screensize'));
leg = cell(1, wind.WS_len);
hold on
for i = 1:wind.WS_len
  plot(out_store{i}.wind.Time, out_store{i}.wind.Data, 'LineWidth', ...
    line_width);
  leg{i} = ['Mean: ', num2str(wind.mean(i)), '; Std: ',...
    num2str(wind.turbulence(i))];
end
legend(leg, 'Location', 'best', 'FontSize', font_size);
xlabel('Time [s]')
ylabel('Wind speed [m/s]')
title('Wind speed time serie')
grid on
set(gca, 'FontSize', font_size)
clear leg

%% Static simulation
% Comparison between the results from the DTU report and the static
% simulation ones

% compute the reference rotational speed
omega_ref = zeros(length(reference(:,1)),1);
for i=1:length(reference(:,1))
  if (reference(i,1) < rated_values(1))
    omega_ref(i) = reference(i, 1).*rated_values(4)/rotor.R;
  else
    omega_ref(i) = omega_rated;
  end
end
omega_ref = omega_ref*30/pi;

fig_omega = figure('Position', get(0, 'Screensize'));
leg = cell(1, wind.WS_len + 2);
hold on
for i = 1:wind.WS_len
  plot(wind.mean(i), ...
    mean(out_store{i}.omega_r.Data(end-simulation.plot_step:end))*30/pi,...
    'o')
  leg{i} = ['Simulation', num2str(i)];
end
plot(reference(:, 1), reference(:, 3)); % DTU reference
plot(reference(:, 1), omega_ref, 'g'); % reference with my TSR
leg{wind.WS_len + 1} =  ['DTU 10MW ref.'];
leg{wind.WS_len + 2} =  ['Computed ref.'];
xlabel('Wind speed [m/s]')
ylabel('\omega_r [rpm]')
legend(leg, 'location', 'best', 'FontSize', font_size);
set(gca, 'FontSize', font_size)
grid on

fig_power = figure('Position', get(0, 'Screensize'));
hold on
for i = 1:wind.WS_len
  plot(wind.mean(i), ...
    mean(out_store{i}.P_r.Data(end-simulation.plot_step:end)), 'o')
end
plot(reference(:, 1), reference(:, 6)*1e3)
xlabel('Wind speed [m/s]')
ylabel('P [W]')
legend(leg, 'location', 'best', 'FontSize', font_size);
title('Mechanical power')
set(gca, 'FontSize', font_size)
grid on


fig_pitch = figure('Position', get(0, 'Screensize'));
hold on
for i = 1:wind.WS_len
  plot(wind.mean(i), ...
   rad2deg(mean(out_store{i}.pitch.Data(end-simulation.plot_step:end))),...
   'o')
end
plot(reference(:, 1), reference(:, 2))
xlabel('Wind speed [m/s]')
ylabel('\theta [°]')
legend(leg, 'location', 'best', 'FontSize', font_size);
set(gca, 'FontSize', font_size)
grid on

%% Dynamic simulation
fig_pitch_dynamic = figure('Position', get(0, 'Screensize'));
hold on
leg(end) = []; % remove last element from the legend name
for i = 1:wind.WS_len
  plot(out_store{i}.pitch.Time, rad2deg(out_store{i}.pitch.Data), ...
    'LineWidth', line_width)
end
xlabel('Wind speed [m/s]')
ylabel('\theta [°]')
legend(leg, 'location', 'best', 'FontSize', font_size)
set(gca, 'FontSize', font_size)
grid on

fig_power_dynamic = figure('Position', get(0, 'Screensize'));
hold on
leg = cell(2*wind.WS_len, 1); % remove last element from the legend name
for i = 1:wind.WS_len
 plot(out_store{i}.P_r.Time, out_store{i}.P_r.Data, ...
    'LineWidth', line_width, 'Color', colors_vect(i, :))
 plot(out_store{i}.P_g.Time, out_store{i}.P_g.Data, ...
  'LineStyle', '--', 'LineWidth', line_width, 'Color', colors_vect(i, :))
  leg{2*i - 1} = ['Aero series ', num2str(i)];
  leg{2*i} = ['Generator series ', num2str(i)];
end
xlabel('Wind speed [m/s]')
ylabel('P [W]')
legend(leg, 'location', 'best', 'FontSize', font_size)
set(gca, 'FontSize', font_size)
grid on