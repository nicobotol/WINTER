%% Plot of the results
% This file is used for print the results of the simulation
close all
clc

date = string(datetime('now','TimeZone','local','Format', ...
        'y_MM_d_HH_mm_ss')); % save the date to identify the figures

%% Wind time series
fig_wind_TS = figure('Position', get(0, 'Screensize'));
leg = cell(1, wind.WS_len);
hold on
for i = 1:wind.WS_len
  plot(out_store{i}.wind.Time, out_store{i}.wind.Data, 'LineWidth', ...
    line_width);
  if simulation.type == 1       % constant
    leg{i} = ['Mean: ', num2str(wind.mean(i)), '; Std: 0'];
  elseif simulation.type == 2   % ramp
    leg{i} = ['Sim. ', num2str(i)];
  elseif simulation.type == 3   % generated wind series
    leg{i} = ['Mean: ', num2str(wind.mean(i)), '; Std: ',...
    num2str(wind.turbulence(i))];
  end
end
legend(leg, 'Location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
xlabel('Time [s]','interpreter','latex')
ylabel('Wind speed [m/s]','interpreter','latex')
title('Wind speed time serie')
grid on
set(gca, 'FontSize', font_size)
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_wind_TS','.png');
  saveas(fig_wind_TS, fig_name, 'png');
end

%% Static simulation
% Comparison between the results from the DTU report and the static
% simulation ones

if simulation.type == 1 || simulation.type == 3 % static or turbulence
clear leg

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
    'o', 'LineWidth', line_width)
  leg{i} = ['Sim.', num2str(i)];
end
plot(reference(:,1),reference(:,3),'LineWidth',line_width); % DTU reference
plot(reference(:,1),omega_ref, 'Color', colors_vect(i + 2,:), ...
  'LineWidth',line_width);%ref. with my TSR
leg{wind.WS_len + 1} =  ['DTU 10MW ref.'];
leg{wind.WS_len + 2} =  ['Computed ref.'];
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('$\omega_{r}$ [rpm]','interpreter','latex')
title('Rotor rotational speed')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_omega','.png');
  saveas(fig_omega, fig_name,'png');
end

fig_power = figure('Position', get(0, 'Screensize'));
hold on
for i = 1:wind.WS_len
  plot(wind.mean(i), ...
    mean(out_store{i}.P_g.Data(end-simulation.plot_step:end)), 'o', ...
    'LineWidth', line_width)
end
plot(reference(:, 1), reference(:, 6)*1e3, 'LineWidth', line_width)
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('P [W]','interpreter','latex')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
title('Generator power')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_power','.png');
  saveas(fig_power, fig_name,'png');
end


fig_pitch = figure('Position', get(0, 'Screensize'));
hold on
for i = 1:wind.WS_len
  plot(wind.mean(i), ...
   rad2deg(mean(out_store{i}.pitch.Data(end-simulation.plot_step:end))),...
   'o', 'LineWidth', line_width)
end
plot(reference(:, 1), reference(:, 2), 'LineWidth', line_width)
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('$\theta \ [deg.]$','interpreter','latex')
title('Collective blade pitch angle')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_pitch','.png');
  saveas(fig_pitch, fig_name,'png');
end

leg(end) = []; % remove last element from the legend name
else
end

%% Dynamic simulation
fig_pitch_dynamic = figure('Position', get(0, 'Screensize'));
hold on
for i = 1:wind.WS_len
  plot(out_store{i}.pitch.Time, rad2deg(out_store{i}.pitch.Data), ...
    'LineWidth', line_width)
end
xlabel('Time [s]','interpreter','latex')
ylabel('$\theta \ [deg.]$ ','interpreter','latex')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex')
set(gca, 'FontSize', font_size)
title('Collective blade pitch angle')
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_pitch_dynamic','.png');
  saveas(fig_pitch_dynamic, fig_name,'png');
end

fig_omega_dynamic = figure('Position', get(0, 'Screensize'));
hold on
for i = 1:wind.WS_len
 plot(out_store{i}.omega_r.Time, out_store{i}.omega_r.Data, ...
    'LineWidth', line_width, 'Color', colors_vect(i, :))
end
xlabel('Time [s]','interpreter','latex')
ylabel('$\omega_{r}$ [rad/s]','interpreter','latex')
title('Rotor rotational speed')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_omega_dynamic', '.png');
  saveas(fig_omega_dynamic, fig_name,'png');
end

fig_power_dynamic = figure('Position', get(0, 'Screensize'));
hold on
leg = cell(2*wind.WS_len + 1, 1);
for i = 1:wind.WS_len
 plot(out_store{i}.P_r.Time, out_store{i}.P_r.Data, ...
    'LineWidth', line_width*0.6, 'Color', colors_vect(i, :))
 plot(out_store{i}.P_g.Time, out_store{i}.P_g.Data, ...
  'LineStyle', '--', 'LineWidth', line_width*1.2, ...
  'Color', colors_vect(i, :))
  leg{2*i - 1} = ['Aero. sim. ', num2str(i)];
  leg{2*i} = ['Gen. sim. ', num2str(i)];
end
yline(rotor.P_rated, 'LineStyle', '-.', 'LineWidth', line_width,...
  'Color', colors_vect(i + 1, :));
leg{2*wind.WS_len + 1} = ['Rated'];
xlabel('Time [s]', 'interpreter','latex')
ylabel('P [W]', 'interpreter','latex')
title('Rotor and generator powers')
legend(leg, 'location', 'best', 'FontSize', font_size, ...
  'interpreter','latex')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_power_dynamic','.png');
  saveas(fig_power_dynamic, fig_name,'png');
end

