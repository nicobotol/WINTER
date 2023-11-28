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
Vbetz = [0:0.1:15];
Pbetz = 16/27*1/2*rho*A*Vbetz.^3;
line_width = 1.2*line_width;


% for v=1:length(V)
%   V0 = V(v);
%   if V0 <= V0_rated
%     lambda = rated_values(4);
%   else
%     lambda = omega_rated*rotor.R/V0;    
%   end

%   cP(:) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda, 0); % interpolate the look-up table
%   P_plot(:, v) = 0.5*rotor.A*V0^3.*cP*rho; % compute the power

% end

fig = figure('Color','w');
grid minor;hold on;box on;
xline(V0_1(end), '-', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off', 'FontSize', font_size);
xline(V0_2(end), '-', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off');
xline(V0_3(end), '-', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off');
xline(V0_4(end), '-', 'LineWidth', line_width, 'Color', color(5), 'HandleVisibility','off', 'FontSize', font_size);
plot(V, P, 'LineWidth', 1.2*line_width, 'DisplayName', 'Mech. power');
plot(Vw, Pw, 'LineWidth', 1.2*line_width, 'DisplayName', 'Wind power');
plot(Vbetz, Pbetz, 'LineWidth', 1.2*line_width, 'DisplayName', 'Betz limit');
x_min = 0;
x_max = 35;
y_min = 0;
y_max = max(Pw);
xlim([x_min x_max]);
ylim([y_min y_max]);
xBox_b = [x_min, V0_3(end), V0_3(end), x_min];
yBox_b = [y_min, y_min, y_max, y_max];
patch(xBox_b, yBox_b, 'black', 'FaceColor', color(1), 'FaceAlpha', 0.1, 'DisplayName', 'Below rated');
xBox_a = [V0_3(end), V0_5(end), V0_5(end), V0_3(end)];
yBox_a = [y_min, y_min, y_max, y_max];
patch(xBox_a, yBox_a, 'black', 'FaceColor', color(3), 'FaceAlpha', 0.3, 'DisplayName', 'Above rated');
xticks([4 11.4 25]);
xticklabels({'Cut in', 'Rated', 'Cut out'});
text(2, 13e6, 'I', 'FontSize', font_size);
text(7.5, 13e6, 'II', 'FontSize', font_size);
text(10.6, 13e6, 'III', 'FontSize', font_size);
text(18.2, 13e6, 'IV', 'FontSize', font_size);
text(30, 13e6, 'V', 'FontSize', font_size);
yticks([P_2(1) P_4(end)]);
yticklabels({'$P_{min}$', '$P_{max}$'});
xlabel('$V_0$ [m/s]'); ylabel('P [W]');
set(gca, 'FontSize', font_size); 
legend('Location', 'SouthEast');
title('Power curve');

if simulation.print_figure == 1
  export_figure(fig, '\operating_reagions.eps', path_images);
end