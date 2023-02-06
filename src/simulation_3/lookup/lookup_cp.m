%% cP as function of pitch angle and Tip Speed ratio (TSR)
% This file is aimed to create a mesh of cP and cT as function of the pitch 
% angle and the tip speed ratio.
% Furthermore the optimum TSR, pitch angle are found.
% Finally the rated wind speed and rotational speed are found 

clear 
close all
clc

% add path for the functions for the aerodynamic
addpath("..\")

% load the parameters
parameters

lookup_cP = zeros(pitch_item, lambda_item); % matrix to store cP
lookup_cT = zeros(pitch_item, lambda_item); % matrix to store cT
for l = 1:lambda_item % loop over the TSR
  lambda = lambda_vector(l); 
  for p = 1:pitch_item  % loop over the pitch 
    theta = pitch_vector(p);
    cP_partial = zeros(1, r_item_no_tip);
    cT_partial = zeros(1, r_item_no_tip);

    % find the local cP and cT along the blade
    [cP_partial, cT_partial,~,~] = cP_cT_partial(r_item_no_tip, ...
      r_vector, beta_vector, thick_vector, c_vector, rotor.blades, ...
      a_guess, a_prime_guess, rotor.R, lambda, theta, aoa_mat, cl_mat, ...
      cd_mat, thick_prof, fake_zero, rho, 0, 0, i_max);

    % integrate the local cP and cT along the blade
    cP = lambda*rotor.blades/(rotor.A*rotor.R)*...
      trapezoidal_integral(r_vector(1:r_item_no_tip), cP_partial); 
    cT = rotor.blades/rotor.A*...
      trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial);

    % store the results
    lookup_cP(p, l) = cP;
    lookup_cT(p, l) = cT;
  end
end

%% Optimium point
max_tmp = max(lookup_cP);       
[~, lambda_pos] = max(max_tmp);                     
[cP_max, theta_pos] = max(lookup_cP(:, lambda_pos)); % max cP
lambda_opt = lambda_vector(lambda_pos);             % TSR for cP_max
theta_opt = pitch_vector(theta_pos);                % pitch for cP_max

%% Plot
% contour plot for cP
contour_plot_cP = figure('Position', get(0, 'Screensize'));
[C, h] = contourf(lambda_vector, rad2deg(pitch_vector), lookup_cP, ...
  'ShowText', 'on'); % To display value on the plot use ,'ShowText','on'
clabel(C,h,'FontSize',font_size*0.8)
hold on
plot(lambda_opt, theta_opt, 'r.', 'MarkerSize',30)
text(7.9, 0.4, num2str(cP_max), 'Color','r', 'FontSize', font_size)
hold off
colorbar()
xlabel('\lambda')
ylabel('\theta_p (°)')
title('Contour plot of c_P')
ax = gca;
ax.FontSize = font_size;

% contour plot for cT
contour_plot_cT = figure('Position', get(0, 'Screensize'));
[C, h] = contourf(lambda_vector, rad2deg(pitch_vector), lookup_cT, ...
  'ShowText', 'on');
clabel(C,h,'FontSize',font_size*0.8)
colorbar()
xlabel('\lambda')
ylabel('\theta_p (°)')
title('Contour plot of c_T')
ax = gca;
ax.FontSize = font_size;

%% Rated wind speed and rated velocity
V0_rated = (rotor.P_rated/(0.5*rho*cP_max*rotor.A))^(1/3);  % [m/s]
omega_rated = V0_rated*lambda_opt/rotor.R;                  % [rad/s]
omega_rated_rpm = 30/pi*omega_rated;                        % [rpm]
rated_values(1) = V0_rated;
rated_values(2) = omega_rated;
rated_values(3) = omega_rated_rpm;
rated_values(4) = lambda_opt;
rated_values(5) = cP_max;

% % Save the results 
save('lookup_cP_theta_lambda.mat', 'lookup_cP');
save('lookup_cT_theta_lambda.mat', 'lookup_cT');
save('rated_values.mat', 'rated_values');