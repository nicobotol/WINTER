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
addpath("aerodynamic_functions");
addpath("aerodynamic_functions\airfoil_data")

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
[cP_max, theta_pos] = max(lookup_cP(:, lambda_pos));  % max cP
lambda_opt = lambda_vector(lambda_pos);               % TSR for cP_max
theta_opt = pitch_vector(theta_pos);                  % pitch for cP_max

% Remove the outlyer
while cP_max > 1
  cP_tmp = mean([lookup_cP(theta_pos - 1, lambda_pos), lookup_cP(theta_pos + 1, lambda_pos), lookup_cP(theta_pos, lambda_pos - 1), lookup_cP(theta_pos, lambda_pos + 1)]);
  lookup_cP(theta_pos, lambda_pos) = cP_tmp;
  max_tmp = max(lookup_cP);       
  [~, lambda_pos] = max(max_tmp);                     
  [cP_max, theta_pos] = max(lookup_cP(:, lambda_pos));  % max cP
  lambda_opt = lambda_vector(lambda_pos);               % TSR for cP_max
  theta_opt = pitch_vector(theta_pos);                  % pitch for cP_max
end

%% Plot
% contour plot for cP
contour_plot_cP = figure('Position', get(0, 'Screensize'), 'Color', 'w');
[~, lambda_pos_5] = min(abs(lambda_vector - 5));
[~, lambda_pos_10] = min(abs(lambda_vector - 10));
[~, theta_pos_2] = min(abs(pitch_vector + 2*pi/180));
[~, theta_pos_5] = min(abs(pitch_vector - 5*pi/180));
[C, h] = contourf(lambda_vector(lambda_pos_5:lambda_pos_10), rad2deg(pitch_vector(theta_pos_2:theta_pos_5)), lookup_cP(theta_pos_2:theta_pos_5, lambda_pos_5:lambda_pos_10),'ShowText', 'on'); % To display value on the plot use ,'ShowText','on'
clabel(C,h,'FontSize',font_size*0.8, 'Interpreter', 'latex')
hold on
plot(lambda_opt, theta_opt, 'r.', 'MarkerSize',30)
text(7.9, 0.4, num2str(cP_max, '%.3f'), 'Color','r', 'FontSize', font_size, 'Interpreter', 'latex')
hold off
colorbar()
xlabel('$\lambda$', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
title('Contour plot of $c_P$', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cP, 'contour_plot_cP.eps', path_images);

% contour plot for cT
contour_plot_cT = figure('Position', get(0, 'Screensize'), 'Color', 'w');
[C, h] = contourf(lambda_vector(lambda_pos_5:lambda_pos_10), rad2deg(pitch_vector(theta_pos_2:theta_pos_5)), lookup_cT(theta_pos_2:theta_pos_5, lambda_pos_5:lambda_pos_10), 'ShowText', 'on');
clabel(C,h,'FontSize',font_size*0.8, 'Interpreter', 'latex')
colorbar()
xlabel('$\lambda$', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
title('Contour plot of $c_T$', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cT, 'contour_plot_cT.eps', path_images);

%% Rated wind speed and rated velocity
V0_rated = (rotor.P_rated/(0.5*rho*cP_max*rotor.A))^(1/3);  % [m/s]
omega_rated = V0_rated*lambda_opt/rotor.R;                  % [rad/s]
omega_rated_rpm = 30/pi*omega_rated;                        % [rpm]
rated_values(1) = V0_rated;
rated_values(2) = omega_rated;
rated_values(3) = omega_rated_rpm;
rated_values(4) = lambda_opt;
rated_values(5) = cP_max;
rated_values 


%% Comparison with the NREL 5MW results
NREL_data = readmatrix('Cp_vs_TSR_vs_pitch.txt', "CommentStyle", "%");
NREL_data_theta = NREL_data(1,2:end);
NREL_data_lambda = NREL_data(2:end,1);
[~, t_min_pos] = min(abs(NREL_data_theta+1.8));
[~, t_max_pos] = min(abs(NREL_data_theta-5));
[~, l_min_pos] = min(abs(NREL_data_lambda-5));
[~, l_max_pos] = min(abs(NREL_data_lambda-10));
NREL_data_theta = NREL_data(1,t_min_pos:t_max_pos);
NREL_data_lambda = NREL_data(l_min_pos:l_max_pos,1);
max_NREL =  max(NREL_data(2:end, 2:end), [], 'all'); % max cp value
NREL_contour_plot_cP = figure('Position', get(0, 'Screensize'), 'Color', 'w');
[C, h] = contourf(NREL_data_lambda, NREL_data_theta, NREL_data(l_min_pos:l_max_pos, t_min_pos:t_max_pos)','ShowText', 'on'); % To display value on the plot use ,'ShowText','on'
clabel(C,h,'FontSize',font_size*0.8, 'Interpreter', 'latex')
% hold on
% plot(lambda_opt, theta_opt, 'r.', 'MarkerSize',30)
% text(7.9, 0.4, num2str(cP_max, '%.3f'), 'Color','r', 'FontSize', font_size, 'Interpreter', 'latex')
% hold off
colorbar()
xlabel('$\lambda$', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
title('Contour plot of $c_P$ from reference data', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cP, 'contour_plot_cP.eps', path_images);

lambda_vec = 5:1:10; % tip speed ratio
pitch = -2:1:5;    % [deg]
NREL_data_theta = NREL_data(1,2:end);
NREL_data_lambda = NREL_data(2:end,1);
figure();
hold on;
plot([NaN, NaN], 'o', 'Color', 'k', 'DisplayName', 'Ref.')
plot([NaN, NaN], '-', 'Color', 'k', 'DisplayName', 'Comp.')
for i = 1:length(pitch)
  cp_comp =  interp2(lambda_vector, pitch_vector, lookup_cP, lambda_vec, pitch(i)*3.14/180);
  cp_NREL = interp2(NREL_data_lambda, NREL_data_theta, NREL_data(2:end, 2:end)', lambda_vec, pitch(i));
  plot(lambda_vec, cp_comp, 'Color', color(i), 'HandleVisibility', 'off', 'LineWidth', line_width);
  plot(lambda_vec, cp_NREL, 'o', 'Color', color(i), 'DisplayName', [num2str(pitch(i)), '$^\circ$' ], 'LineWidth', line_width);
end
legend('Location', 'EastOutside')
title('Comparison between reference and computed $c_P$ for NREL 5 MW')
xlabel('$\lambda$ [-]')
ylabel('$c_P [-]$')
grid on
box on

tmp(1, :) = kron(lambda_vector, ones(1, pitch_item));
tmp(2, :) = kron(ones(1, lambda_item), pitch_vector);
tmp(3, :) = reshape(lookup_cP, 1, lambda_item*pitch_item);
idx = tmp(3, :) > 0;
figure()
plot3(tmp(1, idx), tmp(2, idx)*180/3.14, tmp(3, idx), 'o')
xlabel('$\lambda$')
ylabel('$\theta [^\circ]$')
zlabel('$c_P$')

%%
% Save the results 
save('lookup\lookup_cP_theta_lambda.mat', 'lookup_cP');
save('lookup\lookup_cT_theta_lambda.mat', 'lookup_cT');
save('rated_values.mat', 'rated_values');