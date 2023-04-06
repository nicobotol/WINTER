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
V0_vect = 12:0.5:25;
V0_vect(end + 1) = V0_rated;
V0_vect = sort(V0_vect);
V0_len = length(V0_vect);

% parameters for the angle generation
theta_offset = 5*pi/180;
delta_theta = 0.1*pi/180;

% Initializations
P = cell(1, V0_len); %  power for each wind speed for each pitch 
dPdtheta = cell(1, V0_len); % derivative of the power wrt to pitch 
theta_vector = cell(1, V0_len);
theta_pos = zeros(V0_len, 1);
P_rated =  zeros(V0_len, 1);
dPdt =  zeros(V0_len, 1);
cP_rated =  zeros(V0_len, 1);
theta =  zeros(V0_len, 1); % angle to have the rated power [rad]
M =  cell(1, V0_len);
Torque = zeros(V0_len, 1);
dMdt = zeros(V0_len, 1);
dMdtheta  =  cell(1, V0_len);
min_t = zeros(V0_len, 1);
cP_th = zeros(V0_len, 1);

for v=1:V0_len
  V0 = V0_vect(v);
  lambda = omega_rated*rotor.R/V0; 
  theta(v) = interp1(lookup_Pitch(1,:), lookup_Pitch(3,:), V0); % angle corresponding to the pitch one [rad]

  % Determine the values of the induction factors along the blade
  [cp_part, ~, ~, ~, a_store, a_prime_store] = ...
    cP_cT_partial(r_item_no_tip, r_vector, beta_vector, thick_vector, ...
    c_vector, rotor.blades, a_guess, a_prime_guess, rotor.R, lambda, theta(v), ...
    aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0, ...
    omega_rated, i_max);
  
  % Power coefficient and aerodynamical power produced
%   cP_th(v) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, theta(v)); % cP from the lookup
  cP_rated(v) = lambda*rotor.blades/(rotor.A*rotor.R)*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_part);
  P_rated(v) = 0.5*rho*rotor.A*V0^3*cP_rated(v); % power produced at rated pitch [W]

  % Induction faction with "frozen wake" 
  theta_vector{v}(:) = theta(v)-theta_offset:delta_theta:theta(v)+theta_offset; % vector to use for the derivative
%   theta_vector{v}(end + 1) = theta(v); % add the rated pitch angle
%   theta_vector{v}(:) = sort(theta_vector{v}(:));
  [~, theta_pos(v)] = min(abs(theta_vector{v}(:) - theta(v))); % position of the angle corresponding to the one of pitch
  
  len_theta = length(theta_vector{v});
  cP = zeros(len_theta, 1);
  
  for t=1:len_theta % loop over the different angle
    [cp_partial, ~, pt_partial, ~] = cP_cT_partial_FW( a_store, a_prime_store,...
      r_item_no_tip, r_vector, beta_vector, thick_vector, c_vector, rotor.R, ...
      lambda, theta_vector{v}(t), aoa_mat, cl_mat, cd_mat, thick_prof, rho, V0, omega_rated);
  
    M{v}(t) = rotor.blades*trapezoidal_integral(r_vector(1:r_item_no_tip),...
      pt_partial.*r_vector(1:r_item_no_tip)); % torque for the range of pitch angles
    cP(t) = lambda*rotor.blades/(rotor.A*rotor.R)*...
        trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial); 
  end
  Torque(v) = M{v}(theta_pos(v)); % torque at rated pitch but extracted from the vector of multiple pitches
  Power = Torque*omega_rated; % power at rated pitch but extracted from the vector of multiple pitches
  P{v}(:) = 0.5*rotor.A*cP*V0^3*rho;

  % Differentiation of the power wrt the pitch angle
  dPdtheta{v} = diff(P{v}(:))./diff(theta_vector{v}(:)); % 0.031; % 2*pi/180
  dPdt(v) = dPdtheta{v}(theta_pos(v));

  % Differentiation of the torque wrt the pitch angle
  dMdtheta{v} = diff(M{v})./diff(theta_vector{v}); %*0.031
  dMdt(v) = dMdtheta{v}(theta_pos(v));
end

% power produced at each velocity
P_store = 0.5*rho*V0_vect.^3*rotor.A*cP_rated;

% Interpolation of the power derivative wrt the pitch
A = ones(V0_len, 2);
A(:, 1) = theta*180/pi; % [deg]
coeff_dPdtheta = (A'*A)\A'*dPdt;

A2 = ones(V0_len, 3);
A2(:, 1) = (theta*180/pi).^2; % [deg]
A2(:, 2) = theta*180/pi; % [deg]
coeff_dPdtheta2 = (A2'*A2)\A2'*dPdt;

% Interpolation of the torque derivative wrt the pitch
coeff_dMdtheta = (A'*A)\A'*dMdt;

coeff_dMdtheta2 = (A2'*A2)\A2'*dMdt;


% Find the KK value for the gain scheduling
% solve(coeff(1)*t + coeff(2) == 2*coeff(2), t);
KK_deg = coeff_dPdtheta(2)/coeff_dPdtheta(1); % [deg]
%KK = KK_deg*180/pi;                     % [rad]

% Gain schdeling
I_eq_HS = rotor.I*gearbox.ratio^2 + generator.I;  % inertia at the High Speed side
theta_v = [0:1:25];                               % [deg]
GK = 1./(1 + theta_v/KK_deg);                     % gain reduction
kp = 2*I_eq_HS*omega_rated*30/pi*blade.zeta_phi*blade.omega_phin/...
  (-gearbox.ratio*coeff_dPdtheta(2))*GK;             % proportional gain
ki = I_eq_HS*omega_rated*30/pi*blade.omega_phin^2/...
  (-gearbox.ratio*coeff_dPdtheta(2))*GK;             % integral gain [1/s]

% Interpolating gain scheduling
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
A = ones(V0_len, 2);
A(:, 1) = V0_vect;
coeff_V0 = (A'*A)\A'*dPdt;

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
dPdtheta_eval = polyval(coeff_dPdtheta, theta*180/pi);
figure()
plot(theta*180/pi, dPdt, '-o', 'DisplayName', 'Computed')
hold on
% plot(theta*180/pi, lookup_dPdtheta, '-o', 'DisplayName','Rated')
plot(theta*180/pi, dPdtheta_eval, 'DisplayName', 'Interpolation')
xlabel('[deg]')
ylabel('[W/rad]')
legend()
title('Derivative of the power w.r.t. the pitch')

% Torque deriative vs pitch
dMdtheta_eval = polyval(coeff_dMdtheta, theta*180/pi);
dMdtheta_eval2 = polyval(coeff_dMdtheta2, theta*180/pi);
DTU10MW_ref = -2.8*(1 + theta_v/164.13 + theta_v.^2/702.09);
figure()
plot(theta*180/pi, dMdt/1e6, '-o', 'DisplayName', 'Computed')
hold on
plot(theta*180/pi, dMdtheta_eval/1e6, 'DisplayName', 'Interpolation 1')
plot(theta*180/pi, dMdtheta_eval2/1e6, 'DisplayName', 'Interpolation 2')
plot(theta_v, DTU10MW_ref, 'DisplayName','DTU ref')
xlabel('[deg]')
ylabel('[MNm/rad]')
legend()
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
figure()
for i=1:V0_len
  plot(theta_vector{i}(1:end-1)*180/pi, dPdtheta{i}, 'DisplayName',num2str(V0_vect(i)))
  hold on
end
title('Derivative of the power wrt the pitch')

% Torque deriviative for each WS
% figure()
% for i=1:V0_len
%   plot(theta_vector{i}(1:end-1)*180/pi, dMdtheta{i}, 'DisplayName',num2str(V0_vect(i)))
%   hold on
% end
% title('Derivative of the torque wrt the pitch')

% Gains vs pitch
figure()
yyaxis left
plot(theta_v, GK, 'DisplayName','GK', 'LineWidth',line_width)
ylabel('GK')
yyaxis right
plot(theta_v, kp, 'DisplayName','kp', 'LineWidth',line_width)
hold on
plot(theta_v, ki, 'DisplayName','ki', 'LineWidth',line_width)
legend()
grid on
xlabel('[deg]')
ylabel('kp ki gains')

% Gains vs pitch
theta_expanded = -5:1:30;
ki_expanded = polyval(coeff_ki, theta_expanded);
kp_expanded = polyval(coeff_kp, theta_expanded);
figure()
yyaxis left
plot(theta_v, kp, 'DisplayName','kp', 'LineWidth',line_width, 'Color',colors_vect(1,:))
hold on
plot(theta_expanded, kp_expanded,'--', 'DisplayName','kp interp', 'LineWidth',line_width, 'Color',colors_vect(1,:))
% plot(theta_v, polyval(blade.kp_schedule_report, theta_v*pi/180), '.-', 'DisplayName','kp ref')
ylabel('kp gains')
yyaxis right
plot(theta_v, ki, 'DisplayName','ki', 'LineWidth',line_width, 'Color',colors_vect(2,:))
hold on
plot(theta_expanded, ki_expanded,'--', 'DisplayName','ki interp', 'LineWidth',line_width, 'Color',colors_vect(2,:))
% plot(theta_v, polyval(blade.ki_schedule_report, theta_v*pi/180), '.-', 'DisplayName','ki ref')
ylabel(' ki gains')
legend()
grid on
xlabel('[deg]')

% rated torque and power vs V0
figure()
yyaxis left
plot(V0_vect, omega_rated*Torque/1e6);
ylabel('Power [MW]')
yyaxis right
plot(V0_vect, Torque/1e6);
xlabel('V0 [m/s]')
ylabel('Torque [MNm]')
title('Generator torque and power')

figure()
for v=1:V0_len
  plot(theta_vector{v}, P{v}/1e6);
  plot(theta_vector{v}(theta_pos(v)), omega_rated*Torque/1e6, 'o')
  hold on
end
xlabel('$\theta$ [rad]')
ylabel('P [MW]')
grid on

figure()
plot(V0_vect, P_rated/1e6)
hold on
plot(V0_vect, Power/1e6, 'o')
title('Power at rated pitch')

% Comparison between my and rated cP
figure()
plot(V0_vect, cP_rated, 'o', 'DisplayName', 'Computed here')
hold on
plot(V0_vect, cP_th)
title('cP at rated pitch')

% Comparison between computed and lookup pitch
figure()
plot(V0_vect, theta, 'o')
hold on
plot(V0_vect, interp1(lookup_Pitch(1,:), lookup_Pitch(3,:), V0_vect))
title('Pitch comparison')

blade_schedule_gains = cell(2, 1);
blade_schedule_gains{1} = coeff_kp;
blade_schedule_gains{2} = coeff_ki;
save 'C:\Users\Niccolò\Documents\UNIVERSITA\TESI_MAGISTRALE\src\simulation_5\lookup\blade_schedule_gains.mat' blade_schedule_gains;


% figure()
% for v=1:V0_len
%   plot3(V0_vect(v)*ones(len_theta, 1), theta_vector{v}, P{v}/1e6);
%   hold on
% end
% xlabel('$V_0$ [m/s]')
% ylabel('$\theta$ [rad]')
% zlabel('P [MW]')
% grid on