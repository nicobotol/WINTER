%% This functions shows the polynomials for gain schduling

close all
clear 
clc

parameters

dtheta = 0.5;
theta = [-5:dtheta:25];
[~, theta_0_pos] = min(abs(theta - 0));
theta_rad = theta*pi/180;

kp_gain = polyval(blade.kp_schedule, theta_rad);
ki_gain = polyval(blade.ki_schedule, theta_rad);

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


% fig_name = strcat(path_images,'\fig_gain_scheduling','.png');
% export_fig('fig', fig_name);