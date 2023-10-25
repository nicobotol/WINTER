clear;
close all;

parameters;

N = 1e6;
rho_10 = 0.1*rho;
sigma_rho = (rho_10/3);
rho_d = random('Normal', rho, sigma_rho, N, 1);
q_10_d = quantile(rho_d, 0.01);
q_90_d = quantile(rho_d, 0.99);

figure();hold on;grid on;box on;
plot(NaN, NaN, 'k--', 'DisplayName', '10 percentile', 'LineWidth', line_width)
plot(NaN, NaN, 'k-.', 'DisplayName', '90 percentile', 'LineWidth', line_width)
histogram(rho_d, 'Normalization', 'pdf', 'BinMethod', 'sturges', 'FaceColor', color(1), 'DisplayName', 'Undeform Dist.')
xline([q_10_d], 'b--', 'LineWidth', line_width, 'FontSize', font_size, 'Handlevisibility', 'off')
xline([q_90_d], 'b-.', 'LineWidth', line_width, 'FontSize', font_size, 'Handlevisibility', 'off')
xline(rho, 'k-', 'FontSize', font_size, 'LineWidth', line_width, 'DisplayName', 'Undeform mean')
xlabel('R [m]')
ylabel('Probability [-]')
legend()
set(gca, 'FontSize', font_size, 'LineWidth', line_width)
