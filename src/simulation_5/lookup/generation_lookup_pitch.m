%% This code computes the pitching angle 
% This file is aimed to create a look up table for the pitch angle of the
% turbine. The parameters are loaded in parameters.m and the cP lookup
% table is used

clear 
close all
clc

% load the physical parameters from the file "parameters.m"
addpath("..\")
addpath("..\aerodynamic_functions")
addpath("..\aerodynamic_functions\airfoil_data")
parameters

%% 
P = zeros(p_item, 1);                   % vector of power
stall = zeros(velocity_item, 1);        % pitch angle for stall [rad]
feather = zeros(velocity_item, 1);      % pitch angle for feathering [rad]
stall_deg = zeros(velocity_item, 1);    % stall angle [°]
feather_deg = zeros(velocity_item, 1);  % feathering angle [°]
cp = zeros(p_item, 1);

for v = 1:velocity_item % loop over different velocities
  
  V0 = velocity_vector(v);
  if V0 <= (V0_rated)
    stall(v) = 0;
    feather(v) = 0;
  else
    lambda = omega_rated*rotor.R/V0;

    % adapt the space where to look for solutions
    p_vector = linspace(stall_lim, feather_lim, p_item); 
    
    cP(:) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, ...
      p_vector); % interpolate the look-up table
    P(:) = 0.5*rotor.A*V0^3.*cP*rho; % compute the power
  
    [~, max_pos] = max(P); % find the peak of the power curve
  
    % split the power and pitch vectors before and after the power peak
    P_part1 = P(1:max_pos);                       % before the peak 
    P_part2 = P(max_pos:end);                     % after the peak
    Theta_p_part1 = p_vector(1 : max_pos);
    Theta_p_part2 = p_vector(max_pos : end);
    
    % interpolate the power curve to find the pitch corresponding to the
    % rated power
    stall(v) = interp1(P_part1, Theta_p_part1, rotor.P_rated);    % [rad]
    feather(v) = interp1(P_part2, Theta_p_part2, rotor.P_rated);  % [rad]
    feathering_deg = rad2deg(feather);                            % [deg]
    stall_deg = rad2deg(stall);                                   % [deg]
  
    % adapt the stall and feather limits
    stall_lim = max(stall(v) - 8*pi/180, pitch_range(1));
    feather_lim = min(feather(v) + 8*pi/180, pitch_range(2));
  end
end

%% Generate the data for the plots
v_plot = [4, 9, 12, 18, 25];
v_plot_len = length(v_plot);
P_plot = zeros(p_item, v_plot_len);
cp = zeros(p_item, 1);
p_vector = linspace(-15, 25, p_item)*pi/180;
for v=1:v_plot_len
  V0 = v_plot(v);
  if V0 <= V0_rated
    lambda = rated_values(4);
  else
    lambda = omega_rated*rotor.R/V0;    
  end

  cP(:) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, ...
      p_vector); % interpolate the look-up table
  P_plot(:, v) = 0.5*rotor.A*V0^3.*cP*rho; % compute the power

end

%% Save the results
lookup_Pitch = zeros(3, velocity_item);
lookup_Pitch(1, :) = velocity_vector;
lookup_Pitch(2, :) = stall;
lookup_Pitch(3, :) = feather;
save('lookup_pitch.mat', 'lookup_Pitch');

%% Plot the results

% Plot the pitch angles
fig_pitch_vs_V0 = figure('Color','w');
plot(velocity_vector', feathering_deg, 'LineWidth', line_width)
hold on
plot(velocity_vector', stall_deg, 'LineWidth', line_width)
plot(reference(:, 1), reference(:, 2), 'ro', 'LineWidth', line_width)
legend('Feather', 'Stall', 'Validation','Location', 'southwest' ...
  )
grid on
box on
xlabel('Wind speed [m/s]')
ylabel('Pitch [deg.]')
title('Pitch as function of windspeed')
set(gca, 'FontSize', font_size)
export_figure(fig_pitch_vs_V0, '\fig_pitch_vs_V0.eps', path_images);

%% Plot the power curve for some velocities
fig_power_vs_pitch = figure('Color','w');
hold on;
leg = cell(v_plot_len + 1, 1);
for i = 1:v_plot_len
  plot(p_vector'*180/pi, P_plot(:, i)/1e6, 'LineWidth', line_width);
  leg{i} = ['$V_0$ = ', num2str(v_plot(i)), ' [m/s]'];
end
yline(rotor.P_rated/1e6, '--r', 'LineWidth', line_width); 
leg{end} = 'Rated power';
legend(leg, 'location', 'southwest', 'FontSize', font_size, ...
  'Interpreter','latex','NumColumns',3);
xlabel('Pitch [deg.]')
ylabel('Power [MW]')
% ylim([-1e6, 6e7]);
grid on
box on
title('Mechanical power as function of pitch angle')
set(gca, 'FontSize', font_size)
export_figure(fig_power_vs_pitch, '\fig_power_vs_pitch.eps', path_images);