%% This code computes the pitching angle 
% This file is aimed to create a look up table for the pitch angle of the
% turbine. The parameters are loaded in parameters.m and the cP lookup
% table is used

clear 
close all
clc

% load the physical parameters from the file "parameters.m"
parameters

P = zeros(p_item, 1);               % vector of power
stall = zeros(velocity_item, 1);    % pitch angle for stall [rad]
feather = zeros(velocity_item, 1);  % pitch angle for feathering [rad]
stall_deg = zeros(velocity_item, 1);    % stall angle [°]
feather_deg = zeros(velocity_item, 1);  % feathering angle [°]

for v = 1:velocity_item % loop over different velocities
  
  V0 = velocity_vector(v);
  if V0 <= (V0_rated)
    stall(v) = 0;
    feather(v) = 0;
  else
    lambda = omega_rated*rotor.R/V0;

    % adapt the space where to look for solutions
    p_vector = linspace(stall_lim, feather_lim, p_item); 
    
    for p = 1:p_item % loop over different pitch
      Theta_p = p_vector(p);
      cP = interp2(lambda_vector, pitch_vector, cP_store, lambda, ...
        Theta_p); % interpolate the look-up table
      P(p) = 0.5*rotor.A*V0^3*cP*rho; % compute the power
    end
  
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

%% Save the results
lookup_Pitch = zeros(3, velocity_item);
lookup_Pitch(1, :) = velocity_vector;
lookup_Pitch(2, :) = stall;
lookup_Pitch(3, :) = feather;
save('lookup_pitch.mat', 'lookup_Pitch');

%% Plot the results

fig_pitch = figure();
plot(velocity_vector', feathering_deg, 'LineWidth', line_width)
hold on
plot(velocity_vector', stall_deg, 'LineWidth', line_width)
plot(velocity_reference, pitch_reference, 'ro', 'LineWidth', line_width)
legend('Feather', 'Stall', 'Validation','Location', 'southwest')
grid on
xlabel('Wind speed [m/s]')
ylabel('Pitch angle [°]')
title('Pitch as function of windspeed')
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)