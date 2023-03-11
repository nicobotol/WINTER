%% Plot of the results
% This file is used for print the results of the simulation
close all

date = string(datetime('now','TimeZone','local','Format', ...
        'y_MM_d_HH_mm_ss')); % save the date to identify the figures

%% Wind time series
if simulation.type ~= 4
fig_wind_TS = figure('Position', get(0, 'Screensize'), 'Color','w');
leg = cell(1, wind.WS_len);
hold on
for i = 1:wind.WS_len
  plot(out_store{i}.wind.Time, out_store{i}.wind.Data, 'LineWidth', ...
    line_width);
  leg{i} = ['Sim. ', num2str(i)];
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
  export_fig('fig_wind_TS', fig_name);
end
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

fig_omega = figure('Position', get(0, 'Screensize'), 'Color','w');
leg = cell(1, wind.WS_len + 2);
hold on
for i = 1:wind.WS_len
  plot(wind.mean(i), ...
    mean(out_store{i}.omega_R.Data(end-simulation.plot_step:end))*30/pi,...
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
  export_fig('fig_omega', fig_name);
end

fig_power = figure('Position', get(0, 'Screensize'), 'Color','w');
hold on
for i = 1:wind.WS_len
  plot(wind.mean(i), ...
    mean(out_store{i}.P_G.Data(end-simulation.plot_step:end)), 'o', ...
    'LineWidth', line_width)
end
plot(reference(:, 1), reference(:, 6)*1e-3, 'LineWidth', line_width)
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('P [MW]','interpreter','latex')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
title('Generator power')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_power','.png');
  export_fig('fig_power', fig_name);
end

fig_pitch = figure('Position', get(0, 'Screensize'), 'Color','w');
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
  export_fig('fig_pitch', fig_name);
end

leg(end) = []; % remove last element from the legend name
else
end

%% Dynamic simulation
if simulation.type ~= 4 

fig_pitch_dynamic = figure('Position', get(0, 'Screensize'), 'Color','w');
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
  export_fig('fig_pitch_dynamic', fig_name);
end

fig_omega_dynamic = figure('Position', get(0, 'Screensize'), 'Color','w');
hold on
for i = 1:wind.WS_len
 plot(out_store{i}.omega_R.Time, out_store{i}.omega_R.Data, ...
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
  export_fig('fig_omega_dynamic', fig_name);
end

fig_power_dynamic = figure('Position', get(0, 'Screensize'), 'Color','w');
hold on
leg = cell(2*wind.WS_len + 1, 1);
for i = 1:wind.WS_len
 plot(out_store{i}.P_R.Time, out_store{i}.P_R.Data/1e6, ...
    'LineWidth', line_width*0.6, 'Color', colors_vect(i, :))
 plot(out_store{i}.P_G.Time, out_store{i}.P_G.Data/1e6, ...
  'LineStyle', '--', 'LineWidth', line_width*1.2, ...
  'Color', colors_vect(i, :))
  leg{2*i - 1} = ['Aero. sim. ', num2str(i)];
  leg{2*i} = ['Gen. sim. ', num2str(i)];
end
yline(rotor.P_rated/1e6, 'LineStyle', '-.', 'LineWidth', line_width,...
  'Color', colors_vect(i + 1, :));
leg{2*wind.WS_len + 1} = ['Rated'];
xlabel('Time [s]', 'interpreter','latex')
ylabel('P [MW]', 'interpreter','latex')
title('Rotor and generator powers')
legend(leg, 'location', 'best', 'FontSize', font_size, ...
  'interpreter','latex')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_power_dynamic','.png');
  export_fig('fig_power_dynamic', fig_name);
end

fig_torque_dynamic = figure('Position', get(0, 'Screensize'), 'Color','w');
hold on
leg = cell(2*wind.WS_len, 1);
for i = 1:wind.WS_len
 plot(out_store{i}.T_G_reference.Time, ...
   out_store{i}.T_G_reference.Data/1e6,'LineWidth',line_width*0.6, ...
   'LineStyle', '--','Color',colors_vect(i, :))
 plot(out_store{i}.T_G.Time, out_store{i}.T_G.Data/1e6, ...
  'LineWidth', line_width*1.2,'Color', colors_vect(i, :))
  leg{2*i - 1} = ['Ref. sim. ', num2str(i)];
  leg{2*i} = ['Sim. ', num2str(i)];
end
xlabel('Time [s]', 'interpreter','latex')
ylabel('T [MNm]', 'interpreter','latex')
title('Generator torque')
legend(leg, 'location', 'best', 'FontSize', font_size, ...
  'interpreter','latex')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_torque_dynamic','.png');
  export_fig('fig_torque_dynamic', fig_name);
end

end

%% Paramertrization plot
% print only the steady state time horizon (i.e. the last 
% simulation.plot_time seconds of the simulation)
if simulation.type == 5 || simulation.type == 6

leg = cell(wind.WS_len + 1, 1);     % struct for the legend

% Elements to be plotted
% low sample time
print_index_L = zeros(wind.WS_len, 2); % transform time into indeces
print_index_L(:, 1) = (simulation.stop_time - simulation.plot_time)/ ...
  simulation.time_step_L;
print_index_L(:, 2) = simulation.stop_time/simulation.time_step_L;
print_index_L = ceil(print_index_L);
% high sample time
print_index_H = zeros(wind.WS_len, 2); % transform time into indeces
print_index_H(:, 1) = (simulation.stop_time - simulation.plot_time)/ ...
  simulation.time_step_H;
print_index_H(:, 2) = simulation.stop_time/simulation.time_step_H;
print_index_H = ceil(print_index_H);

wind_resampled_L = cell(wind.WS_len, 1);
wind_resampled_H = cell(wind.WS_len, 1);
% Interpolate the time for the two simulations
for i=1:wind.WS_len
  % interpolate the simulation with smaller time
  wind_resampled_L{i} = interp1(out_store{i}.wind.Time, ...
  out_store{i}.wind.Data, out_store{i}.P_G.Time);
  % interpolate the simulation with higher time
  wind_resampled_H{i} = interp1(out_store{i}.wind.Time, ...
  out_store{i}.wind.Data, out_store{i}.omega_R.Time);
end

fig_power_param = figure('Position', get(0, 'Screensize'),'Color','w');
for i = 1:wind.WS_len
  plot(wind_resampled_L{i}(print_index_L(i, 1):print_index_L(i, 2)), ...
    out_store{i}.P_G.Data(print_index_L(i, 1):print_index_L(i, 2))/1e6, ...
    'LineWidth', line_width, 'Color', colors_vect(i, :));
  hold all;
  leg{i} = ['Sim. ', num2str(i)];
end
plot(reference(:, 1), reference(:, 6)*1e-3, '--', 'LineWidth', line_width)
leg{wind.WS_len + 1} =  ['DTU 10MW ref.'];
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('P [MW]','interpreter','latex')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
title('Generator power')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_power_param','.png');
  export_fig('fig_power_param', fig_name);
end

fig_gen_torque_param=figure('Position', get(0, 'Screensize'),'Color','w');
for i = 1:wind.WS_len
  % plot
  plot(wind_resampled_L{i}(print_index_L(i, 1):print_index_L(i, 2)), ...
    out_store{i}.T_G.Data(print_index_L(i, 1):print_index_L(i, 2))/1e6, ...
    'LineWidth', line_width, 'Color', colors_vect(i, :));
  hold all;
end
plot(reference(:, 1), reference(:, 7)/1e6, '--', 'LineWidth', line_width)
leg{wind.WS_len + 1} =  ['DTU 10MW ref.'];
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('T [MNm]','interpreter','latex')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
title('Generator torque')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_gen_torque_param','.png');
  export_fig('fig_gen_torque_param', fig_name);
end

fig_pitch_param = figure('Position', get(0, 'Screensize'),'Color','w');
for i = 1:wind.WS_len
  % plot
  plot(wind_resampled_L{i}(print_index_L(i, 1):print_index_L(i, 2)), ...
    out_store{i}.pitch.Data(print_index_L(i, 1):print_index_L(i, 2)) ...
    *180/pi, 'LineWidth', line_width, 'Color', colors_vect(i, :));
  hold all;
end
plot(reference(:, 1), reference(:, 2), '--', 'LineWidth', line_width)
leg{wind.WS_len + 1} =  ['DTU 10MW ref.'];
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('$\theta$ [deg.]','interpreter','latex')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
title('Pitch angle')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_pitch_param','.png');
  export_fig('fig_pitch_param', fig_name);
end

fig_omega_param = figure('Position', get(0, 'Screensize'),'Color','w');
for i = 1:wind.WS_len
  % plot
  plot(wind_resampled_L{i}(print_index_L(i, 1):print_index_L(i, 2)), ...
    out_store{i}.omega_R.Data(print_index_L(i, 1):print_index_L(i, 2)) ...
    *30/pi, 'LineWidth', line_width, 'Color', colors_vect(i, :));
  hold all;
end
plot(reference(:, 1), reference(:, 3), '--', 'LineWidth', line_width)
leg{wind.WS_len + 1} =  ['DTU 10MW ref.'];
xlabel('Wind speed [m/s]','interpreter','latex')
ylabel('$\omega$ [rpm]','interpreter','latex')
legend(leg, 'location', 'best', 'FontSize', font_size,...
  'interpreter','latex');
title('Rotor rotational speed')
set(gca, 'FontSize', font_size)
grid on
if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_omega_param','.png');
  export_fig('fig_omega_param', fig_name);
end
  
end

%% Generator step response
if simulation.type == 4

  N_start = ceil(step_t_start/simulation.time_step_L); % item where start 
  N_end = ceil(step_t_stop/simulation.time_step_L);    % item where stop 

  fig_gen_step = figure('Position', get(0, 'Screensize'),'Color','w');
  plot(out.T_G_reference.Time, out.T_G_reference.Data, '--',...
    'LineWidth', line_width);
  hold on
  plot(out.T_G.Time, out.T_G.Data, 'LineWidth', line_width*0.75);
  legend('Reference', 'Response', 'location', 'best', 'FontSize', ...
    font_size, 'interpreter','latex');
  grid on;
  title('Generator torque step response')
  xlabel('Time [s]', 'interpreter','latex')
  ylabel('T [Nm]', 'interpreter','latex')
  set(gca, 'FontSize', font_size)
  rectangle('Position', [step_t_start step_y_min ...
    step_step_t_stop-step_t_start  step_y_max-step_y_min])

  axes('position',[2.8/5 0.4/1.2 .30 .30])
  plot(out.T_G_reference.Time(N_start:N_end), ...
    out.T_G_reference.Data(N_start:N_end), '--','LineWidth', 1.5, ...
    'MarkerSize',6, 'Color', colors_vect(1,:)); %axis tight
  hold on
  plot(out.T_G.Time(N_start:N_end), ...
    out.T_G.Data(N_start:N_end), 'LineWidth', 1.5, ...
    'MarkerSize',6, 'Color', colors_vect(2,:)); %axis tight
  xlim([step_t_start step_t_stop])
  ylim([step_y_min step_y_max]);
  grid on
  title('Overshoot of the response')
  set(gca, 'FontSize', font_size)

  if simulation.print_figure == 1
  fig_name = strcat(path_images,'\', date, 'fig_gen_step','.png');
  export_fig('fig_gen_step', fig_name);
  end

end
