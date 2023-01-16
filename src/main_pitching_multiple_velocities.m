%% This code computes the pitching angle 
% Inputs parameters for the code to work, to be setted in parameters.m are
% the wind speed V0, the rotational speed omega, and the power where to
% pitch at

clear 
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS IN parameters.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load the physical parameters from the file "parameters.m"
parameters

% load the airfoil data from the different files
[aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);

% load the blade data from "bladedat.txt"
[r_vector, c_vector, beta_vector, thick_vector] = ...
load_blade_data(blade_filename);

r_item = size(r_vector, 2); % number of cross sections along the blade
r_item_no_tip = r_item - 1; % numer of cross section to take into account

% pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item); % 
% vector of pitch equally distributed in the range 
P = zeros(pitch_item, 1); % vector of power where to store the power

stall = zeros(velocity_item, 1);
feather = zeros(velocity_item, 1);

% chose Theta_p
for v=1:velocity_item
  V0 = velocity_pitch(v);
  lambda = omega * R / V0;

  % adapt the space where to look fpr solutions
  pitch_vector = linspace(stall_lim - 8*pi/180, feather_lim + 6*pi/180, pitch_item); % 
  
  for p = 1:pitch_item % loop over different pitch
    Theta_p = pitch_vector(p);
  
    [cp_partial, cT_partial,~,~] = cP_cT_partial(r_item_no_tip, ...
      r_vector, beta_vector, thick_vector, c_vector, B, a_guess, ...
      a_prime_guess, R, lambda, Theta_p, aoa_mat, cl_mat, cd_mat, ...
      thick_prof, fake_zero, rho,V0,omega, i_max);
  
    P(p) = 0.5*omega*B*rho*V0^2*trapezoidal_integral( ...
      r_vector(1:r_item_no_tip), cp_partial);
    
  end

  [~, max_pos] = max(P); % find the element correspondignto the peak of the
  % power curve

  % split the power vector in two parts: the first before the max power and
  % the second after
  P_part1 = P(1:max_pos);
  P_part2 = P(max_pos:end);
  % split the corresponding theta
  Theta_p_part1 = pitch_vector(1 : max_pos);
  Theta_p_part2 = pitch_vector(max_pos : end);
  
  stall(v) = interp1(P_part1, Theta_p_part1, P_rated); % pitch angle for 
  % stalling
  feather(v) = interp1(P_part2, Theta_p_part2, P_rated); % pitch angle for 
  % feathering
  
  feathering_deg = rad2deg(feather); % convert hte anagle in degrees
  stall_deg = rad2deg(stall);
  
  disp('%%%%%%%%')
  disp(strcat('At the velocity of ', num2str(V0)))
  disp(strcat('The angle for feathering is:', num2str(feather(v)), ...
    ' [rad] = ', num2str(feathering_deg(v)), ' [°]'));
  disp(strcat('The angle for stalling is:', num2str(stall(v)), ' [rad] = ', ...
    num2str(stall_deg(v)), ' [°]'));

  % adapt the stall and feather limits
  stall_lim = stall(v);
  feather_lim = feather(v);
end

%% PLOT THE RESULTS
feathering_deg(1) = 0; % cheat a little bit to have better graphs
stall_deg(1) = 0;

fig_pitch = figure('Position', get(0, 'Screensize'));
plot(velocity_pitch', feathering_deg, 'LineWidth', line_size)
hold on
plot(velocity_pitch', stall_deg, 'LineWidth', line_size)
plot([0 min(velocity_pitch)], [0 0], 'LineWidth', line_size )
legend('Feather', 'Stall', 'Below rated', 'Location', 'southwest')
grid on
xlabel('Wind speed [m/s]')
ylabel('Pitch angle [°]')
title('Pitch as function of windspeed')
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_pitch, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ' ...
  'ANNO\WIND_ENERGY\exercise_pitching\fig_pitch.png'],'png');

