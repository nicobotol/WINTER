% this function is necessary to identify the lower and upper value where  to distribute the values of Kopt

close all;
clc;
clear;

parameters;


mu_omega = 0.1;
mu_rho = 0.01;
mu_R = 2;
mu_V0_rated = 1;
mu_theta = 1*pi/180; 

N = 1e6;

theta_d = random('Normal', 0, mu_theta, N, 1);
rho_d = random('Normal', rho, mu_rho, N, 1);
R_d = random('Normal', rotor.R, mu_R, N, 1);
omega_d = random('Normal', omega_rated_GE, mu_omega, N, 1);
V0_d = random('Normal', V0_rated, mu_V0_rated, N, 1);

lambda_d = omega_d.*R_d./V0_d;
cp_d = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_d, theta_d);

K = 0.5*cp_d.*rho_d.*pi.*R_d.^2.*(V0_d./omega_d).^3;
q_10 = quantile(K, 0.1);
q_90 = quantile(K, 0.9);

figure()
histogram(K, 'Normalization', 'pdf', 'BinMethod', 'sturges')
xlabel('K')
ylabel('pdf')
grid()
xline(generator.K_opt_GE, 'r--', 'LineWidth', 2)
xline(q_10, 'b--', 'LineWidth', 2)
xline(q_90, 'b--', 'LineWidth', 2)


figure()
cdfplot(K)
yline(0.1, 'r--', 'LineWidth', 2)
yline(0.9, 'r--', 'LineWidth', 2)
