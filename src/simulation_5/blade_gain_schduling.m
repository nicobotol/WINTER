%% This functions shows the polynomials for gain schduling

close all
clear 
clc

parameters

zeta_phi = 0.65;
omega_phi_eta = 0.66;

dtheta = 0.5;
theta = [-5:dtheta:25];
[~, theta_0_pos] = min(abs(theta - 0));
theta_rad = theta*pi/180;

cp_vect = zeros(length(theta_rad), 1);
cp_vect = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_opt, theta_rad);
P = 0.5*rotor.A*V0_rated^3*rho.*cp_vect;
dPdtheta = diff(P)/dtheta;
figure()
plot(theta(1:end-1), dPdtheta);
grid on

dPdtheta_rat = dPdtheta(theta_0_pos);
[~, pos_KK]  = min(abs(dPdtheta - 2*dPdtheta_rat));
KK = 0.00001;%theta_rad(pos_KK);
kp = 2*I_eq*omega_rated*zeta_phi*omega_phi_eta/(-gearbox.ratio*dPdtheta_rat);
ki = I_eq*omega_rated*omega_phi_eta^2/(-gearbox.ratio*dPdtheta_rat);

dtheta = 0.5;
theta = [dtheta:dtheta:25];
theta_rad = theta*pi/180;

GK = 1./(1 + theta_rad/KK);
Ki = GK*ki;
Kp = GK*kp;

kp_gain = polyval(blade.kp_schedule, theta_rad);
ki_gain = polyval(blade.ki_schedule, theta_rad);

fig = figure('Position', get(0, 'Screensize'), 'Color','w');
yyaxis left
plot(theta, Kp,'LineWidth', line_width);
ylabel('$k_{p}$ [-]', 'Interpreter','latex')
yyaxis right
plot(theta, Ki,'LineWidth', line_width);
ylabel('$k_{i} [s^{-1}]$', 'Interpreter','latex')
xlabel('Pitch angle [deg.]')
title('Gain scheduling coefficients', 'Interpreter','latex')
grid on
set(gca, 'FontSize', font_size)

fig = figure('Position', get(0, 'Screensize'), 'Color','w');
yyaxis left
plot(theta, kp_gain,'LineWidth', line_width);
ylabel('$k_{p}$ [-]', 'Interpreter','latex')
yyaxis right
plot(theta, ki_gain,'LineWidth', line_width);
ylabel('$k_{i} [s^{-1}]$', 'Interpreter','latex')
xlabel('Pitch angle [deg.]')
title('Gain scheduling coefficients', 'Interpreter','latex')
grid on
set(gca, 'FontSize', font_size)


fig_name = strcat(path_images,'\fig_gain_scheduling','.png');
export_fig('fig', fig_name);