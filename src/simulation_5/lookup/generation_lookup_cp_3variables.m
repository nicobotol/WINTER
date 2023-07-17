%% cP as function of pitch angle, V0 and omega
% This file is aimed to create a mesh of cP and cT as function of the pitch 
% angle, the wind velocity and the rotational speed

clear 
close all
clc

% add path for the functions for the aerodynamic
addpath("..\")
addpath("..\aerodynamic_functions")
addpath("..\aerodynamic_functions\airfoil_data")


% load the parameters
parameters

lookup_cP_3v = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); % matrix to store cP
lookup_cT_3v = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); % matrix to store cT
A = zeros(omega_item_3v, velocity_item_3v, pitch_item_3v); % matrix to store the feasibility
A_store = []; % matrix to store the feasibility
V0_max = max(velocity_vector_3v);
theta_max = max(pitch_vector_3v);
omega_max = max(omega_vector_3v);
contour_plot_A = figure('Position', get(0, 'Screensize'), 'Color', 'w');
hold on
for k = 1:omega_item_3v
  lhs = [];
  omega = omega_vector_3v(k);
  % rhs = V0_max*(omega*theta_max - omega_max*pitch_vector_3v)/(omega_max*theta_max);
  rhs = pitch_vector_3v*V0_max/theta_max*ones(1, velocity_item_3v)';
  lhs = velocity_vector_3v;
  A(k, :, :) = lhs.' + rhs < 0.1;

  A_tmp = reshape(A(k, :, :), 1, velocity_item_3v*pitch_item_3v);
  V0_tmp = kron(ones(1, pitch_item_3v), velocity_vector_3v);
  theta_tmp = kron(pitch_vector_3v, ones(1, velocity_item_3v));
  omega_tmp = omega_vector_3v(k)*ones(1, velocity_item_3v*pitch_item_3v);
  A_store = [A_store, [omega_tmp; theta_tmp; V0_tmp; A_tmp]];
end

keep = A_store(4, :) == 1;
scatter3(A_store(1, keep), A_store(2, keep), A_store(3, keep), 40, A_store(4, keep), 'filled')
hold off
colorbar()
xlabel('$\omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
zlabel('$V_0$ [m/s]', 'Interpreter', 'latex')
title('3 varaibles $c_P$ map', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cP, 'contour_plot_cP_3v.eps', path_images);

store = [];
for i = 1:omega_item_3v
  lambda = omega_vector_3v(i)*rotor.R./velocity_vector_3v;
  cP = interp2(lambda_vector, pitch_vector, lookup_cP, lambda', pitch_vector_3v); % interpolate the look-up table
  cT = interp2(lambda_vector, pitch_vector, lookup_cT, lambda', pitch_vector_3v); % interpolate the look-up table
  cP(isnan(cP)) = 0; % set to zero the NaN values
  cT(isnan(cT)) = 0; % set to zero the NaN values
  cP(cP < 0) = 0; % set to zero the negative values
  cT(cT < 0) = 0; % set to zero the negative values

  % store the results
  lookup_cP_3v(i, A(i, :, :)>0) = cP( A(i, :, :)>0);
  lookup_cT_3v(i, A(i, :, :)>0) = cT( A(i, :, :)>0);

  cP_tmp = reshape(cP', 1, velocity_item_3v*pitch_item_3v);
  cT_tmp = reshape(cT', 1, velocity_item_3v*pitch_item_3v);
  V0_tmp = kron(ones(1, pitch_item_3v), velocity_vector_3v);
  theta_tmp = kron(pitch_vector_3v, ones(1, velocity_item_3v));
  omega_tmp = omega_vector_3v(i)*ones(1, velocity_item_3v*pitch_item_3v);
  store = [store, [omega_tmp; theta_tmp; V0_tmp; cP_tmp; cT_tmp]];
end

%% Plot
% contour plot for cP
keep = A_store(4, :) > 0;
contour_plot_cP = figure('Position', get(0, 'Screensize'), 'Color', 'w');
scatter3(store(1, keep), store(2, keep), store(3, keep), 40, store(4, keep), 'filled')
hold off
colorbar()
xlabel('$\omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
zlabel('$V_0$ [m/s]', 'Interpreter', 'latex')
title('3 varaibles $c_P$ map', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cP, 'contour_plot_cP_3v.eps', path_images);

% contour plot for cT
contour_plot_cT = figure('Position', get(0, 'Screensize'), 'Color', 'w');
scatter3(store(1, keep), store(2, keep), store(3, keep), 40, store(5, keep), 'filled')
hold off
colorbar()
xlabel('$\omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
zlabel('$V_0$ [m/s]', 'Interpreter', 'latex')
title('3 varaibles $c_T$ map', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cT, 'contour_plot_cT_3v.eps', path_images);

%%
% Save the results 
% save('lookup_cP_3v_theta_lambda.mat', 'lookup_cP_3v');
% save('lookup_cT_3v_theta_lambda.mat', 'lookup_cT_3v');