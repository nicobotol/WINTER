%% unstedy cP as function of pitch angle, V0 and omega
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

%% Extract some configurations that are feasible. This is done to avoid to have cT grather than 1. This conditions happens for some configurations that are not feasible such as having high wind speed but also low rotational speed and no pitch
lookup_cP_3v = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); % matrix to store cP
lookup_cT_3v = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); % matrix to store cT
A = zeros(omega_item_3v, velocity_item_3v, pitch_item_3v); % matrix to store the feasibility
A_store = []; % matrix to store the feasibility
A_cell = cell(omega_item_3v, 1);
V0_max = max(velocity_vector_3v);
theta_max = max(pitch_vector_3v)*1;
omega_max = max(omega_vector_3v);
contour_plot_A = figure('Position', get(0, 'Screensize'), 'Color', 'w');
hold on
a = [0 0 0];
b = [omega_max, theta_max/2, 0];
c = [omega_max/2, theta_max/3, V0_max];
for i = 1:omega_item_3v
  lhs = [];
  omega = omega_vector_3v(i);
  rhs = plane(a, b, c, omega, pitch_vector_3v);
  lhs = velocity_vector_3v;
  A(i, :, :) = lhs.' - rhs < 0.1;
  A_cell{i} = lhs.' - rhs < 0.1;
  A_tmp = reshape(A(i, :, :), 1, velocity_item_3v*pitch_item_3v);
  V0_tmp = kron(ones(1, pitch_item_3v), velocity_vector_3v);
  theta_tmp = kron(pitch_vector_3v, ones(1, velocity_item_3v));
  omega_tmp = omega_vector_3v(i)*ones(1, velocity_item_3v*pitch_item_3v);
  A_store = [A_store, [omega_tmp; theta_tmp; V0_tmp; A_tmp]];
end

keep = A_store(4, :) == 1;
% scatter3(A_store(1, keep), A_store(2, keep), A_store(3, keep), 40, A_store(4, keep), 'filled')
scatter3(A_store(1, :), A_store(2, :), A_store(3, :), 40, A_store(4, :), 'filled')
hold off
colorbar()
xlabel('$\omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
zlabel('$V_0$ [m/s]', 'Interpreter', 'latex')
title('3 varaibles $c_P$ map', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cP, 'contour_plot_cP_3v.eps', path_images);

%% Compute the cp and cT for the configurations found previously
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
  % lookup_cP_3v(i, A_cell{i} > 0) = cP(A_cell{i} > 0);
  % lookup_cT_3v(i, A_cell{i} > 0) = cT(A_cell{i} > 0);
  lookup_cP_3v(i,:,:) = cP;
  lookup_cT_3v(i,:,:) = cT;
  
  cP_tmp = reshape(cP', 1, velocity_item_3v*pitch_item_3v);
  cT_tmp = reshape(cT', 1, velocity_item_3v*pitch_item_3v);
  V0_tmp = kron(ones(1, pitch_item_3v), velocity_vector_3v);
  theta_tmp = kron(pitch_vector_3v, ones(1, velocity_item_3v));
  omega_tmp = omega_vector_3v(i)*ones(1, velocity_item_3v*pitch_item_3v);
  store = [store, [omega_tmp; theta_tmp; V0_tmp; cP_tmp; cT_tmp]];
end

%% Genration cP and cT unsteady map
lookup_cP_3v_u_1 = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); 
lookup_cP_3v_u_2 = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); 
lookup_cT_3v_u_1 = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); 
lookup_cT_3v_u_2 = zeros(omega_item_3v, pitch_item_3v, velocity_item_3v); 
u_store = [];

for i = 1:omega_item_3v
  V0 = velocity_vector_3v;
  cT_reshape = reshape(lookup_cT_3v(i,:,:), pitch_item_3v, velocity_item_3v);
  cP_reshape = reshape(lookup_cP_3v(i,:,:), pitch_item_3v, velocity_item_3v);
  lambda_s1(i,:,:) = V0/2.*(1 + sqrt(1 - cT_reshape));
  lambda_s2(i,:,:) = V0/2.*(1 - sqrt(1 - cT_reshape));
  
  lambda_s1_reshape = reshape(lambda_s1(i,:,:), pitch_item_3v, velocity_item_3v);
  lambda_s2_reshape = reshape(lambda_s2(i,:,:), pitch_item_3v, velocity_item_3v);

  lookup_cP_3v_u_1(i,:,:) = cP_reshape./(1-lambda_s1_reshape./V0).^3;
  lookup_cP_3v_u_2(i,:,:) = cP_reshape./(1-lambda_s2_reshape./V0).^3;
  lookup_cT_3v_u_1(i,:,:) = cT_reshape./(1-lambda_s1_reshape./V0).^2;
  lookup_cT_3v_u_2(i,:,:) = cT_reshape./(1-lambda_s2_reshape./V0).^2;

  l_cP_rs_1 = reshape(lookup_cP_3v_u_1(i,:,:), pitch_item_3v, velocity_item_3v);
  l_cP_rs_2 = reshape(lookup_cP_3v_u_2(i,:,:), pitch_item_3v, velocity_item_3v);
  l_cT_rs_1 = reshape(lookup_cT_3v_u_1(i,:,:), pitch_item_3v, velocity_item_3v);
  l_cT_rs_2 = reshape(lookup_cT_3v_u_2(i,:,:), pitch_item_3v, velocity_item_3v);
  

  cP_1_tmp = reshape(l_cP_rs_1', 1, velocity_item_3v*pitch_item_3v);
  cP_2_tmp = reshape(l_cP_rs_2', 1, velocity_item_3v*pitch_item_3v);
  cT_1_tmp = reshape(l_cT_rs_1', 1, velocity_item_3v*pitch_item_3v);
  cT_2_tmp = reshape(l_cT_rs_2', 1, velocity_item_3v*pitch_item_3v);
  omega_tmp = omega_vector_3v(i)*ones(1, velocity_item_3v*pitch_item_3v);
  V0_tmp = kron(ones(1, pitch_item_3v), velocity_vector_3v);
  theta_tmp = kron(pitch_vector_3v, ones(1, velocity_item_3v));
  
  u_store = [u_store, [omega_tmp; theta_tmp; V0_tmp; cP_1_tmp; cP_2_tmp; cT_1_tmp; cT_2_tmp]];
end

% Remove from the u_store mat the values corresponding to lmbda_1 because they are not phesible 
keep = A_store(4, :) > 0;
lookup_cP_u = zeros(4, sum(keep == 1));
lookup_cT_u = zeros(4, sum(keep == 1));
lookup_cP_u(1:3, :) = u_store(1:3, keep);
lookup_cP_u(4, :)   = u_store(5, keep);
lookup_cT_u(1:3, :) = u_store(1:3, keep);
lookup_cT_u(4, :)   = u_store(7, keep);

%% Plot
% contour plot for cP
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

% contour plot for cP_u
contour_plot_cP_u = figure('Position', get(0, 'Screensize'), 'Color', 'w');
scatter3(u_store(1, keep), u_store(2, keep), u_store(3, keep), 40, u_store(5, keep), 'filled')
hold off
colorbar()
xlabel('$\omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
zlabel('$V_0$ [m/s]', 'Interpreter', 'latex')
title('3 varaibles $c_P^u$ map', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cT_u, 'contour_plot_cT_3v_u.eps', path_images);

contour_plot_cT_u = figure('Position', get(0, 'Screensize'), 'Color', 'w');
scatter3(u_store(1, keep), u_store(2, keep), u_store(3, keep), 40, u_store(7, keep), 'filled')
hold off
colorbar()
xlabel('$\omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\theta_p [^\circ]$', 'Interpreter', 'latex')
zlabel('$V_0$ [m/s]', 'Interpreter', 'latex')
title('3 varaibles $c_T^u$ map', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = font_size;
% export_figure(contour_plot_cT_u, 'contour_plot_cT_3v_u.eps', path_images);

%%
% Save the results
save('lookup_cP_3v_omega_theta_V0.mat', 'lookup_cP_3v');
save('lookup_cT_3v_omega_theta_V0.mat', 'lookup_cT_3v');
save('lookup_cP_3v_omega_theta_V0_u.mat', 'lookup_cP_u');
save('lookup_cT_3v_omega_theta_V0_u.mat', 'lookup_cT_u');