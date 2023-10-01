%% This functions shows the polynomials for gain schduling

close all
clear 
clc

parameters

dtheta = 0.5;
theta = [-5:dtheta:25];
[~, theta_0_pos] = min(abs(theta - 0));
theta_rad = theta*pi/180;

kp_gain = polyval(blade.kp_schedule_report, theta_rad);
ki_gain = polyval(blade.ki_schedule_report, theta_rad);

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
box on


if simulation.print_figure == 1
  fig_name = strcat('C:\Users\Niccolò\Documents\UNIVERSITA\TESI_MAGISTRALE\report\images\vectorial\fig_gain_scheduling','.eps');
  export_fig(fig, fig_name);
end