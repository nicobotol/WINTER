% This function plots the different types of regions 

clear; close all; clc;
parameters;
A = rotor.A;

V0_1 = [0:0.1:4];
V0_2 = [4:0.1:11];
V0_3 = [11:0.1:11.4];
V0_4 = [11.4 25];
V0_5 = [25 35];
V = [V0_1 V0_2 V0_3 V0_4 V0_5];
P_1 = zeros(1, length(V0_1));
P_2 = 1/2*rho*A*V0_2.^3*cp_max;
P_3 = 1/2*rho*A*V0_3.^3*cp_max;
P_4 = P_3(end)*[1 1];
P_5 = [0 0];
P = [P_1 P_2 P_3 P_4 P_5];
Vw = [0:0.1:10];
Pw = 1/2*rho*A*Vw.^3;
line_width = 1.2*line_width;

fig = figure('Color','w');
grid minor;hold on;box on;
plot(V, P, 'LineWidth', line_width, 'DisplayName', 'Mech. power');
plot(Vw, Pw, 'LineWidth', line_width, 'DisplayName', 'Wind power');
xline(V0_1(end), ':', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off', 'FontSize', font_size);
xline(V0_2(end), ':', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off');
xline(V0_3(end),  ':', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off');
xline(V0_4(end), ':', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off', 'FontSize', font_size);
xlim([0 35])
xticks([4 25]);
xticklabels({'Cut in', 'Cut out'});
text(2, 14e6, 'I', 'FontSize', font_size);
text(7.5, 14e6, 'II', 'FontSize', font_size);
text(10.6, 14e6, 'III', 'FontSize', font_size);
text(18.2, 14e6, 'IV', 'FontSize', font_size);
text(30, 14e6, 'V', 'FontSize', font_size);
yticks([P_2(1) P_4(end)]);
yticklabels({'$P_{min}$', '$P_{max}$'});
xlabel('$V_0$ [m/s]'); ylabel('P [W]');
set(gca, 'FontSize', font_size); 
legend('Location', 'SouthEast');
title('Power curve');

if simulation.print_figure == 1
  export_figure(fig, '\operating_reagions.eps', path_images);
end