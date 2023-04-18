%% This code computes the gain scheduling for the blade pitch controller
% The gains for the controller are computed based on a frozen wake BEM, see
% 'NREL 5 MW report' and 'Aerodynamics of Wind Turbines' by M.O.L. Hansen
clear 
close all
clc

warning ('off','all');

addpath("..\..\")
addpath("..\..\lookup")

% load parameters
parameters_NREL5MW

% wind speeds
V0_vect = gs.V0_vect;         % WS [m/s]
V0_vect(end + 1) = V0_rated;
V0_vect = sort(V0_vect);
V0_len = length(V0_vect);

% parameters for the angle generation
theta_offset = gs.theta_offset;
delta_theta = gs.delta_theta;

% Initializations
dPdtheta = cell(1, V0_len);     % derivative of the power wrt to pitch 
dTdtheta = cell(1, V0_len);     % derivative of the torque wrt to pitch
theta_vector = cell(1, V0_len); % vector of pitch for each WS
theta = zeros(V0_len, 1);       % pitch for rated power [rad]
P = cell(1, V0_len);            % power for each wind speed for each pitch
T = cell(1, V0_len);            % torque for each wind speed for each pitch
theta_pos = zeros(V0_len, 1);   % position of rated pitch in theta_vector 
dPdt = zeros(V0_len, 1);        % power derivative at rated pitch angle
dTdt = zeros(V0_len, 1);        % torque derivative at rated pitch angle

for v=1:V0_len
  V0 = V0_vect(v);                                             % [m/s]
  lambda = omega_rated*rotor.R/V0;                             % TSR
  theta(v) = interp1(lookup_Pitch(1,:), lookup_Pitch(3,:), V0);%pitch [rad]

  %% Induction factors in nominal conditions "standard BEM"
  [cp_part, ~, ~, ~, a_store, a_prime_store] = ...
    cP_cT_partial(r_item_no_tip, r_vector, beta_vector, thick_vector, ...
    c_vector, rotor.blades, a_guess, a_prime_guess, rotor.R, lambda, ...
    theta(v), aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0, ...
    omega_rated, i_max);
  
  %% Induction factors with "frozen wake BEM" 
  theta_vector{v}(:) = theta(v)-theta_offset:delta_theta: ...
    theta(v)+theta_offset; % [rad]
  [~, theta_pos(v)] = min(abs(theta_vector{v}(:) - theta(v))); 
  len_theta = length(theta_vector{v});
  cP = zeros(len_theta, 1);
  for t=1:len_theta                       % loop over the different angles
    [cp_partial, ~, pt_partial, ~] = cP_cT_partial_FW( a_store, ...
      a_prime_store, r_item_no_tip, r_vector, beta_vector, ...
      thick_vector, c_vector, rotor.R, lambda, theta_vector{v}(t), ...
      aoa_mat, cl_mat, cd_mat, thick_prof, rho, V0, omega_rated);
  
    T{v}(t) = rotor.blades*trapezoidal_integral( ...
      r_vector(1:r_item_no_tip), pt_partial.*r_vector( ...
      1:r_item_no_tip));                  % [Nm]
    cP(t) = lambda*rotor.blades/(rotor.A*rotor.R)*...
        trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial); % [-]
  end
  P{v}(:) = 0.5*rho*rotor.A*V0.^3.*cP;

  %% Torque and power sensitivities
  % Differentiation of the power wrt the pitch angle
  dPdtheta{v} = diff(P{v}(:))./diff(theta_vector{v}(:));  % [W/rad]
  dPdt(v) = dPdtheta{v}(theta_pos(v));                    % [W/rad]
  % Differentiation of the torque wrt the pitch angle
  dTdtheta{v} = diff(T{v})./diff(theta_vector{v});        % [Nm/rad]
  dTdt(v) = dTdtheta{v}(theta_pos(v));                    % [Nm/rad]
end

%% Interpolations
% Interpolation of the power derivative wrt the pitch
% 1st order interpolation
A = ones(V0_len, 2);
A(:, 1) = theta;                  % [rad]
coeff_dPdt = (A'*A)\A'*dPdt;      % pseudo inversion
% 2nd order interpolation
A2 = ones(V0_len, 3);
A2(:, 1) = theta.^2;              % [rad^2]
A2(:, 2) = theta;                 % [rad]
coeff_dPdt2 = (A2'*A2)\A2'*dPdt;  % pseudo inversion

% Interpolation of the torque derivative wrt the pitch
% 1st order interpolation
coeff_dTdt = (A'*A)\A'*dTdt;      % pseudo inversion
% 2nd order interpolation
coeff_dTdt2 = (A2'*A2)\A2'*dTdt;  % pseudo inversion


% KK value for the gain scheduling
% solve(coeff(1)*t + coeff(2) == 2*coeff(2), t);
KK = coeff_dPdt(2)/coeff_dPdt(1); % [rad]
KK2 = (-coeff_dPdt2(2) - sqrt(coeff_dPdt2(2)^2 + 4*coeff_dPdt2(1)*coeff_dPdt2(3)))/2/coeff_dPdt2(1);

% Gain schdeling
theta_v = gs.theta_v;                         % [rad]
GK = 1./(1 + theta_v/KK2);                     % gain reduction
kp = 2*I_eq*omega_rated*gs.zeta_phi*gs.omega_phin/...
  (-1/gearbox.ratio*coeff_dPdt2(3))*GK;          % proportional gain [-]
ki = I_eq*omega_rated*gs.omega_phin^2/...
  (-1/gearbox.ratio*coeff_dPdt2(3))*GK;          % integral gain [1/s]


% Polinomial interpolation between angle and gain
index_kp = 7; % index of the interpolation
AA = ones(length(theta_v), index_kp);
for i=1:index_kp
  AA(:, i) = (theta_v).^(index_kp - i);
end
coeff_kp = (AA'*AA)\AA'*kp';

index_ki = 6; % index of the interpolation
AA = ones(length(theta_v), index_ki);
for i=1:index_ki
  AA(:, i) = (theta_v).^(index_ki - i);
end
coeff_ki = (AA'*AA)\AA'*ki';

%% Plots
close all
% poly interp wrt the wind speed
% A = ones(V0_len, 2);
% A(:, 1) = V0_vect;
% coeff_V0 = (A'*A)\A'*dPdt;

% Power derivative vs WS
% figure()
% plot(V0_vect, dPdt, '-o', 'DisplayName', 'Computed')
% hold on
% % plot(V0_vect, lookup_dPdtheta, '-o', 'DisplayName','Rated')
% plot(V0_vect, polyval(coeff_V0, V0_vect), 'DisplayName', 'Interpolation')
% xlabel('[m/s]')
% ylabel('$\frac{dP}{d\theta}$ [W/rad]')
% title('Derivative of the power w.r.t. the pitch')

% Power derivative vs pitch
dPdtheta_eval = polyval(coeff_dPdt, theta);
dPdtheta2_eval = polyval(coeff_dPdt2, theta);
fig_dPdtheta = figure( 'Color','w'); clf;
plot(theta*180/pi, dPdt/1e6, '-o', 'DisplayName', 'BEM FW', 'LineWidth', line_width)
hold on
plot(lookup_Pitch(3, :)*180/pi, lookup_dPdtheta, '-x', 'DisplayName','Reference', 'LineWidth', line_width)
plot(theta*180/pi, dPdtheta_eval/1e6, 'DisplayName', 'Interp. $1^{st}$','LineWidth', line_width)
plot(theta*180/pi, dPdtheta2_eval/1e6, 'DisplayName', 'Interp. $2^{nd}$','LineWidth', line_width)
xlabel('Pitch [deg]')
ylabel('$\frac{dP}{d\theta}$ [MW/rad]')
legend()
grid on
title('Derivative of the power w.r.t. the pitch')
%export_fig(fig_dPdtheta, [path_images, 'dPdtheta_eval.svg'])

% Torque deriative vs pitch
dTdtheta_eval = polyval(coeff_dTdt, theta);   % [Nm/rad]
dTdtheta_eval2 = polyval(coeff_dTdt2, theta); % [Nm/rad]
DTU10MW_ref = -2.8*(1 + (theta*180/pi)/164.13 + (theta*180/pi).^2/702.09); % [MNm/deg]
figure(2); clf;
plot(theta*180/pi, dTdt/1e6, '-o', 'DisplayName', 'BEM FW', 'LineWidth', line_width)
hold on
plot(theta*180/pi, dTdtheta_eval/1e6, 'DisplayName', 'Interp. $1^{st}$', 'LineWidth', line_width)
plot(theta*180/pi, dTdtheta_eval2/1e6, 'DisplayName', 'Interp. $2^{nd}$', 'LineWidth', line_width)
% plot(theta*180/pi, DTU10MW_ref*pi/180, 'DisplayName','DTU ref', 'LineWidth', line_width)
xlabel('Pitch [deg]')
ylabel('$\frac{dT}{d\theta}$ [MNm/rad]')
legend()
grid on
title('Derivative of the torque w.r.t. the pitch')

% % cP_ref = interp1(lookup_cP(1, :), lookup_cP(2, :), V0_vect);
% figure()
% for i=1:V0_len
%   plot(V0_vect(i), cP_rated{i}, '-o', 'DisplayName', 'Computed')
%   hold on
% end
% % plot(V0_vect, cP_ref, 'DisplayName', 'Reference')
% ylabel('$c_p$ [-]')
% xlabel('wind speed [m/s]')
% title('Power coefficient')

% Power derivative for each WS
figure(3); clf;
for i=1:V0_len
  plot(theta_vector{i}(1:end-1)*180/pi, dPdtheta{i}/1e6, 'DisplayName',num2str(V0_vect(i)))
  hold on
end
plot(theta*180/pi, dPdtheta_eval/1e6, 'DisplayName', 'Interp. $1^{st}$', 'LineWidth', line_width)
xlabel('Pitch [deg]')
ylabel('$\frac{dP}{d\theta}$ [MW/rad]')
title('Derivative of the power w.r.t. the pitch')
grid on


% Torque deriviative for each WS
% figure()
% for i=1:V0_len
%   plot(theta_vector{i}(1:end-1)*180/pi, dTdtheta{i}, 'DisplayName',num2str(V0_vect(i)))
%   hold on
% end
% title('Derivative of the torque wrt the pitch')
%%
% Gains vs pitch
fig_gain_pitch = figure('Color', 'w'); clf;
yyaxis left
plot(theta_v*180/pi, GK, 'DisplayName','GK', 'LineWidth',line_width)
ylabel('Gain reduction GK')
yyaxis right
ylim([0, 0.025])
plot(theta_v*180/pi, kp, '-','DisplayName','kp [s]', 'Color', colors_vect(3,:), 'LineWidth',line_width)
hold on
plot(theta_ref, kp_ref, '--','DisplayName','ref. kp [s]', 'Color', colors_vect(3,:), 'LineWidth',line_width)
plot(theta_v*180/pi, ki, '-', 'DisplayName','ki [-]', 'Color', colors_vect(2,:), 'LineWidth',line_width)
plot(theta_ref, ki_ref, '--', 'DisplayName','ref. ki [-]', 'Color', colors_vect(2,:), 'LineWidth',line_width)
plot(25, 0.025, 'o')
legend()
grid on
xlabel('[deg]')
ylabel('kp, ki gains')
title('Gain reduction and scheduling')
%%
% export_fig(fig_gain_pitch, [path_images, 'dPdtheta_eval.svg'])
%%
% Gains vs pitch
theta_expanded = pi/180*[-2:1:30]; % expand the range in order to see overfitting
ki_expanded = polyval(coeff_ki, theta_expanded);  % [-]
kp_expanded = polyval(coeff_kp, theta_expanded);  % [s]
figure(5); clf;
yyaxis left
plot(theta_v*180/pi, kp, 'DisplayName','BEM FW', 'LineWidth',line_width, 'Color',colors_vect(1,:))
hold on
plot(theta_expanded*180/pi, kp_expanded,'--', 'DisplayName','BEM FW interp.', 'LineWidth',line_width, 'Color',colors_vect(1,:))
plot(theta_v*180/pi, polyval(blade.kp_schedule_report, theta_v), 'o-.', 'LineWidth',line_width, 'DisplayName','ref.')
ylabel('kp [s]')
yyaxis right
plot(theta_v*180/pi, ki, 'DisplayName','BEM FW', 'LineWidth',line_width, 'Color',colors_vect(2,:))
hold on
plot(theta_expanded*180/pi, ki_expanded,'--', 'DisplayName','BEM FW interp.', 'LineWidth',line_width, 'Color',colors_vect(2,:))
plot(theta_v*180/pi, polyval(blade.ki_schedule_report, theta_v), 'o-.', 'LineWidth',line_width, 'DisplayName','ref.')
ylabel('ki [-]')
legend()
grid on
xlabel('Pitch [deg]')
title('Gain scheduling for pitch control')
