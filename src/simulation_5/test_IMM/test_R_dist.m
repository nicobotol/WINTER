clear;
close all;

parameters;

N = 1e6;
sigma_R = (4/3)^2;
R_d = random('Normal', rotor.R, sigma_R, 3*N, 1);
R_d_s = R_d(R_d<=rotor.R);
R_d_s = R_d_s(1:N);
q_10_d = quantile(R_d, 0.05);
q_90_d = quantile(R_d, 0.95);
q_10_d_s = quantile(R_d_s, 0.1);
q_90_d_s = quantile(R_d_s, 0.9);

figure();hold on;grid on;box on;
plot(NaN, NaN, 'k--', 'DisplayName', '10 percentile', 'LineWidth', line_width)
plot(NaN, NaN, 'k-.', 'DisplayName', '90 percentile', 'LineWidth', line_width)
histogram(R_d, 'Normalization', 'pdf', 'BinMethod', 'sturges', 'FaceColor', color(1), 'DisplayName', 'Undeform Dist.')
histogram(R_d_s, 'Normalization', 'pdf', 'BinMethod', 'sturges', 'FaceColor', color(2), 'DisplayName', 'deform Dist.')
xline([q_10_d], 'b--', 'LineWidth', line_width, 'FontSize', font_size, 'Handlevisibility', 'off')
xline([q_90_d], 'b-.', 'LineWidth', line_width, 'FontSize', font_size, 'Handlevisibility', 'off')
xline([q_10_d_s], 'r--', 'LineWidth', line_width, 'FontSize', font_size, 'HandleVisibility', 'off')
xline([q_90_d_s], 'r-.', 'LineWidth', line_width, 'FontSize', font_size, 'HandleVisibility', 'off')
xline(rotor.R, 'k-', 'FontSize', font_size, 'LineWidth', line_width, 'DisplayName', 'Undeform mean')
xlabel('R [m]')
ylabel('Probability [-]')
legend()
set(gca, 'FontSize', font_size, 'LineWidth', line_width)
