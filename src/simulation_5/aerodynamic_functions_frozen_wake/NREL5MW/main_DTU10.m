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

V0_vect = 12:0.5:25;
V0_vect(end + 1) = V0_rated;
V0_vect = sort(V0_vect);
theta_offset = 5*pi/180;
delta_theta = 0.1*pi/180;
P = cell(1, length(V0_vect)); %  power for each wind speed for each pitch 
dPdtheta = cell(1, length(V0_vect)); % derivative of the power wrt to pitch 
theta_vector = cell(1, length(V0_vect));
theta_pos = cell(1, length(V0_vect));
P_tot =  cell(1, length(V0_vect));
dPdt =  zeros(length(V0_vect), 1);
cP_store =  cell(1, length(V0_vect));
theta_store =  zeros(length(V0_vect), 1);
M =  cell(1, length(V0_vect));
Torque = zeros(length(V0_vect), 1);
dMdt = zeros(length(V0_vect), 1);
dMdtheta  =  cell(1, length(V0_vect));

for v=1:length(V0_vect)
  V0 = V0_vect(v);
  lambda = omega_rated*rotor.R/V0; 
  theta = interp1(lookup_Pitch(1,:), lookup_Pitch(3,:), V0); % angle corresponding to the pitch one
  theta_store(v) = theta; % [rad]

  % Determine the values of the induction factors along the blade
  [cp_part, ~, ~, ~, a_store, a_prime_store] = ...
    cP_cT_partial(r_item_no_tip, r_vector, beta_vector, thick_vector, ...
    c_vector, rotor.blades, a_guess, a_prime_guess, rotor.R, lambda, theta, ...
    aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0, ...
    omega_rated, i_max);
  
  % Power coefficient and aerodynamical power produced
  cP_store{v} = lambda*rotor.blades/(rotor.A*rotor.R)*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_part);
  P_tot{v} = 0.5*rho*V0^3*rotor.A*cP_store{v};
        
  % Induction faction with "frozen wake" 
  theta_vector{v} = theta-theta_offset:delta_theta:theta+theta_offset;
  [~, theta_pos{v}] = min(abs(theta_vector{v} - theta)); % position of the angle corresponding to the one of pitch
  len_theta = length(theta_vector{v});
  cP = zeros(len_theta, 1);
  Mr = zeros(len_theta, 1);
  
  for t=1:len_theta
    theta = theta_vector{v}(t);
    [cp_partial, ~, pt_partial, ~] = cP_cT_partial_FW( a_store, a_prime_store,...
      r_item_no_tip, r_vector, beta_vector, thick_vector, c_vector, rotor.R, ...
      lambda, theta, aoa_mat, cl_mat, cd_mat, thick_prof, rho, V0, omega_rated);
  
    Mr(t) = rotor.blades*trapezoidal_integral(r_vector(1:r_item_no_tip), pt_partial.*r_vector(1:r_item_no_tip)); % torque
    cP(t) = lambda*rotor.blades/(rotor.A*rotor.R)*...
        trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial); 
    Torque(v) = Mr(theta_pos{v});
  end
  P{v} = zeros(len_theta, 1);
  P{v}(:) = 0.5*rotor.A*cP*V0^3*rho;
  M{v}(:) = Mr;

  % Differentiation of the power wrt the pitch angle
  dPdtheta{v} = diff(P{v})/delta_theta;
  dPdt(v) = dPdtheta{v}(theta_pos{v});

  % Differentiation of the torque wrt the pitch angle
  dMdtheta{v} = diff(M{v})/delta_theta;
  dMdt(v) = dMdtheta{v}(theta_pos{v});
end

% Interpolation of the power derivative wrt the pitch
A = ones(length(V0_vect), 2);
A(:, 1) = theta_store*180/pi; % [deg]
coeff_theta = (A'*A)\A'*dPdt;

% Find the KK value for the gain scheduling
% solve(coeff(1)*t + coeff(2) == 2*coeff(2), t);
KK_deg = coeff_theta(2)/coeff_theta(1); % [deg]
%KK = KK_deg*180/pi;                     % [rad]

% Gain schdeling
I_eq_HS = rotor.I*gearbox.ratio^2 + generator.I;  % inertia at the High Speed side
theta_v = [0:1:25];                               % [deg]
GK = 1./(1 + theta_v/KK_deg);                     % gain reduction
kp = 2*I_eq_HS*omega_rated*30/pi*blade.zeta_phi*blade.omega_phin/...
  (-gearbox.ratio*coeff_theta(2))*GK;             % proportional gain
ki = I_eq_HS*omega_rated*30/pi*blade.omega_phin^2/...
  (-gearbox.ratio*coeff_theta(2))*GK;             % integral gain [1/s]

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
A = ones(length(V0_vect), 2);
A(:, 1) = V0_vect;
coeff_V0 = (A'*A)\A'*dPdt;

figure()
plot(V0_vect, dPdt, '-o', 'DisplayName', 'Computed')
hold on
% plot(V0_vect, lookup_dPdtheta, '-o', 'DisplayName','Rated')
plot(V0_vect, polyval(coeff_V0, V0_vect), 'DisplayName', 'Interpolation')
xlabel('[m/s]')
ylabel('$\frac{dP}{d\theta}$ [W/rad]')
title('Derivative of the power w.r.t. the pitch')


pitch_eval = polyval(coeff_theta, theta_store*180/pi);
figure()
plot(theta_store*180/pi, dPdt, '-o', 'DisplayName', 'Computed')
hold on
% plot(theta_store*180/pi, lookup_dPdtheta, '-o', 'DisplayName','Rated')
plot(theta_store*180/pi, pitch_eval, 'DisplayName', 'Interpolation')
xlabel('[deg]')
ylabel('[W/rad]')
legend()
title('Derivative of the power w.r.t. the pitch')

figure()
plot(theta_store*180/pi, dMdt/1e6, '-o', 'DisplayName', 'Computed')
hold on
% plot(theta_store*180/pi, pitch_eval, 'DisplayName', 'Interpolation')
xlabel('[deg]')
ylabel('[MNm/rad]')
legend()
title('Derivative of the torque w.r.t. the pitch')

% % cP_ref = interp1(lookup_cP(1, :), lookup_cP(2, :), V0_vect);
% figure()
% for i=1:length(V0_vect)
%   plot(V0_vect(i), cP_store{i}, '-o', 'DisplayName', 'Computed')
%   hold on
% end
% % plot(V0_vect, cP_ref, 'DisplayName', 'Reference')
% ylabel('$c_p$ [-]')
% xlabel('wind speed [m/s]')
% title('Power coefficient')

figure()
for i=1:length(V0_vect)
  plot(theta_vector{i}(1:end-1)*180/pi, dPdtheta{i}, 'DisplayName',num2str(V0_vect(i)))
  hold on
end
title('Derivative of the power wrt the pitch')

figure()
for i=1:length(V0_vect)
  plot(theta_vector{i}(1:end-1)*180/pi, dMdtheta{i}, 'DisplayName',num2str(V0_vect(i)))
  hold on
end
title('Derivative of the torque wrt the pitch')

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

theta_expanded = -5:1:30;
ki_expanded = polyval(coeff_ki, theta_expanded);
kp_expanded = polyval(coeff_kp, theta_expanded);
figure()
yyaxis left
plot(theta_v, kp, 'DisplayName','kp', 'LineWidth',line_width, 'Color',colors_vect(1,:))
hold on
plot(theta_expanded, kp_expanded,'--', 'DisplayName','kp interp', 'LineWidth',line_width, 'Color',colors_vect(1,:))
% plot(theta_v, polyval(blade.kp_schedule_report, theta_v*pi/180)*180/pi, '.-', 'DisplayName','kp ref')
ylabel('kp gains')
yyaxis right
plot(theta_v, ki, 'DisplayName','ki', 'LineWidth',line_width, 'Color',colors_vect(2,:))
hold on
plot(theta_expanded, ki_expanded,'--', 'DisplayName','ki interp', 'LineWidth',line_width, 'Color',colors_vect(2,:))
% plot(theta_v, polyval(blade.ki_schedule_report, theta_v*pi/180)*180/pi, '.-', 'DisplayName','ki ref')
ylabel(' ki gains')
legend()
grid on
xlabel('[deg]')

figure()
plot(V0_vect, Torque);

% blade_schedule_gains = cell(2, 1);
% blade_schedule_gains{1} = coeff_kp;
% blade_schedule_gains{2} = coeff_ki;
% save 'C:\Users\Niccolò\Documents\UNIVERSITA\TESI_MAGISTRALE\src\simulation_5\lookup\blade_schedule_gains.mat' blade_schedule_gains;
% figure()
% for v=1:length(V0_vect)
%   plot3(V0_vect(v)*ones(len_theta, 1), theta_vector{v}, P{v}/1e6);
%   hold on
% end
% xlabel('$V_0$ [m/s]')
% ylabel('$\theta$ [rad]')
% zlabel('P [MW]')
% grid on