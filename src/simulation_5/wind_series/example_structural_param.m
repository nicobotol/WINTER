clear; close all; clc;

parameters

fs = 1/simulation.time_step_H;
h = 150;
T = 300;

[u_R, t, PSD_R] = wind_series(rotor.R , IMM.sigma_R, 1/simulation.time_step_H, h, T);
[u_rho, t, PSD_rho] = wind_series(rho, IMM.sigma_rho, 1/simulation.time_step_H, h, T);
[u_theta, t, PSD_theta] = wind_series(0, IMM.sigma_theta, 1/simulation.time_step_H, h, T);
[u_omega, t, PSD_omega] = wind_series(omega_rated, IMM.sigma_omega, 1/simulation.time_step_H, h, T);


u_R(u_R > rotor.R) = rotor.R;
figure();
plot(t, u_R);
yline(rotor.R - 3*IMM.sigma_R, '--r')

figure();
plot(t, u_omega);
yline(omega_rated + 3*std(u_omega), '--r')
yline(omega_rated - 3*std(u_omega), '--r')