%% This code computes the gain scheduling for the blade pitch controller
% The gains for the controller are computed based on a frozen wake BEM, see
% 'NREL 5 MW report' and 'Aerodynamics of Wind Turbines' by M.O.L. Hansen
clear 
close all
clc

warning ('off','all');

addpath("..\..\")
% addpath("..\lookup")

% load parameters
parameters

% Defien the rated parameters and the angle ranges used for the evaluations
theta_offset = 8*pi/180;  % angle range [rad]
delta_theta = 0.1*pi/180; % angular discretization [rad]
theta_vector = -5*pi/180:delta_theta:theta_offset;
[min_tv, pos_min_tv] = min(abs(theta_vector));
if min_tv > fake_zero % Check if 0 [rad] is included
  theta_vector(end + 1) = 0;
  theta_vector = sort(theta_vector);
  [min_tv, pos_min_tv] = min(abs(theta_vector));
end
len_theta = length(theta_vector);
cP = zeros(len_theta, 1);

P = zeros(len_theta, 1);            % power for each pitch angle [W]
dPdtheta = zeros(len_theta - 1, 1); % deriv. of the power wrt pitch [W/rad] 

%% Deriative of the power wrt the pitcha angle, with "frozen wake"
% Determine the values of the induction factors along the blade, in the
% rated condition
lambda = lambda_opt;  % optimal TSR
theta = 0;            % pitch angle [rad]
V0 = V0_rated;        % rated wind speed [m/s]
[~, ~, ~, ~, a_store, a_prime_store] = ...
  cP_cT_partial(r_item_no_tip, r_vector, beta_vector, thick_vector, ...
  c_vector, rotor.blades, a_guess, a_prime_guess, rotor.R, lambda, ...
  theta, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0, ...
  omega_rated, i_max);

% Perturb the pitch angle and use the already computed induced factors to 
% find the corresponding power  
for t=1:len_theta
  theta = theta_vector(t);  % [rad]
  [cp_partial, ~, ~, ~] = cP_cT_partial_FW( a_store, a_prime_store,...
    r_item_no_tip, r_vector, beta_vector, thick_vector, c_vector, ...
    rotor.R, lambda, theta, aoa_mat, cl_mat, cd_mat, thick_prof, rho, ...
    V0, omega_rated);

  cP(t) = lambda*rotor.blades/(rotor.A*rotor.R)*...
      trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial); % [-]
end
P(:) = 0.5*rotor.A*cP*V0^3*rho; % power [W]

% Compute the derivative of the power wrt the pitch
dPdtheta = diff(P)/delta_theta; % dP/dTheta [W/rad]

%% Gains for the controller
% Angle where the derivative is doubled
dPdt0 = dPdtheta(pos_min_tv);   % dP/dTheta(theta = 0)
[~, pos] = min(abs(dPdtheta - 2*dPdt0));
KK = theta_vector(pos);

theta_v = pi/180*[0:1:25];      % [rad]
theta_v2 = pi/180*[-2.5:1:27.5];    % angle for evalaute overfitting [rad]
GK = 1./(1 + theta_v/KK);       % gain reduction
kp = 2*I_eq*omega_rated*blade.zeta_phi*blade.omega_phin/...
  (-gearbox.ratio*dPdt0)*GK;    % proportional gain
ki = I_eq*omega_rated*blade.omega_phin^2/(-gearbox.ratio*dPdt0)*GK; % [1/s]

% Interpolate the results to have a polynomial fitting
index = 9; % degree of the polynomial fitting
A = zeros(length(theta_v), index);
for i=1:index
  A(:, i) = theta_v.^(index-i);
end
coeffs_kp = (A'*A)\A'*kp';                % poly coeff. prop. term  
coeffs_ki = (A'*A)\A'*ki';                % poly coeff. integral term 
kp_interp = polyval(coeffs_kp, theta_v2); % fit of the proportional terms 
ki_interp = polyval(coeffs_ki, theta_v2); % fit of the integral terms

% Gains from Olimpo's
kp_gain = polyval(blade.kp_schedule, theta_v);
ki_gain = polyval(blade.ki_schedule, theta_v);

%% Plots
% Power and its own derivative
figure()
yyaxis left
plot(theta_vector, P/1e6,'LineWidth', line_width);
ylabel('P [MW]', 'Interpreter','latex','LineWidth', line_width)
yyaxis right
plot(theta_vector(1:end-1), dPdtheta/1e6,'LineWidth', line_width);
hold on
plot(min_tv, dPdt0/1e6, 'x','LineWidth', line_width, 'MarkerSize', 10)
plot(KK, dPdtheta(pos)/1e6, 'o','LineWidth', line_width, 'MarkerSize', 10)
ylabel('$\frac{dP}{d \theta}$ [MW/rad]', 'Interpreter','latex', ...
  'LineWidth', line_width)
xlabel('$\theta$ [rad]')
grid on

% Gains
fig = figure('Position', get(0, 'Screensize'), 'Color','w');
yyaxis left
plot(theta_v*180/pi, kp_gain,'LineWidth', line_width, 'DisplayName', 'Poly');
hold on
plot(theta_v*180/pi, kp,'LineWidth', line_width, 'DisplayName', 'Cmp');
plot(theta_v2*180/pi, kp_interp,'LineWidth', line_width, 'DisplayName', 'Int');
ylabel('$k_{p}$ [-]', 'Interpreter','latex')
yyaxis right
plot(theta_v*180/pi, ki_gain,'LineWidth', line_width, 'DisplayName', 'Poly');
hold on
plot(theta_v*180/pi, ki,'LineWidth', line_width, 'DisplayName', 'Cmp');
plot(theta_v2*180/pi, ki_interp,'LineWidth', line_width, 'DisplayName', 'Int');
ylabel('$k_{i} [s^{-1}]$', 'Interpreter','latex')
xlabel('Pitch angle [deg.]')
title('Gain scheduling coefficients', 'Interpreter','latex')
grid on
set(gca, 'FontSize', font_size)
legend();

% figure()
% plot(theta_v, kp);
% hold on
% plot(theta_v, ki);

% figure()
% for v=1:length(V0_vect)
%   plot3(V0_vect(v)*ones(len_theta, 1), theta_vector{v}, P{v}/1e6);
%   hold on
% end
% xlabel('$V_0$ [m/s]')
% ylabel('$\theta$ [rad]')
% zlabel('P [MW]')
% grid on