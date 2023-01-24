%% cP as function of pitch angle and tip speed ratio
% This file is aimed to create a mesh of cP and cT as function of the pitch 
% angle and the tip speed ratio.
% Furthermore the optimum TSR, pitch angle are found.
% Finally the rated wind speed and rotational speed are found 

clear 
close all
clc

parameters

% distribute lambda and TSR in their ranges
lambda_vector = linspace(lambda_range(1), lambda_range(2), lambda_item); 
pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item);

cP_store = zeros(pitch_item, lambda_item); % matrix to store cP
cT_store = zeros(pitch_item, lambda_item); % matrix to store cT
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
    cP_store(p, l) = cP;
    cT_store(p, l) = cT;
  end
end

%% Optimium point
max_tmp = max(cP_store);       
[~, lambda_pos] = max(max_tmp);                     
[cP_max, theta_pos] = max(cP_store(:, lambda_pos)); % max cP
lambda_opt = lambda_vector(lambda_pos);             % TSR for cP_max
theta_opt = pitch_vector(theta_pos);                % pitch for cP_max

%% Plot
% contour plot for cP
contour_plot_cP = figure('Position', get(0, 'Screensize'));
[C, h] = contourf(lambda_vector, rad2deg(pitch_vector), cP_store, ...
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
[C, h] = contourf(lambda_vector, rad2deg(pitch_vector), cT_store, ...
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

%% Save the results 
save('mesh_cP_theta_lambda.mat', 'cP_store');
save('mesh_cT_theta_lambda.mat', 'cT_store');
save('rated_values.mat', 'rated_values');