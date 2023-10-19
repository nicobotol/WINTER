% this function is necessary to identify the lower and upper value where  to distribute the values of Kopt

close all;
clc;
clear;

parameters;
N = 1e6; % sampes to generate

%                                        _                      _      _ 
%   _ __   __ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___   ___| |_ __| |
%  | '_ \ / _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __| / __| __/ _` |
%  | |_) | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \ \__ \ || (_| |
%  | .__/ \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/ |___/\__\__,_|
%  |_|                                                                   

mu_omega = omega_rated*0.05/3; % fixed as 5% of the nominal value
mu_rho = rho*0.05/3; % fixed as 5% of the nominal value
mu_R = (4/3); % assuming a  deflection of 4 meters
mu_V0_rated = V0_rated*0.15/3; % fixed as 15% of the nominal value 
mu_theta = 1*pi/180/3; % assuming 1 deg of uncertainty

%   ____  _     _        _ _           _   _                 
%  |  _ \(_)___| |_ _ __(_) |__  _   _| |_(_) ___  _ __  ___ 
%  | | | | / __| __| '__| | '_ \| | | | __| |/ _ \| '_ \/ __|
%  | |_| | \__ \ |_| |  | | |_) | |_| | |_| | (_) | | | \__ \
%  |____/|_|___/\__|_|  |_|_.__/ \__,_|\__|_|\___/|_| |_|___/
                                                           
theta_d = random('Normal', 0, mu_theta, N, 1); % pitch angle
rho_d = random('Normal', rho, mu_rho, N, 1); % air density
R_d = random('Normal', rotor.R, mu_R, 3*N, 1); % undeformed rotor radius
R_d_s = R_d(R_d<=rotor.R); % take only the radius smaller than the undeformed
R_d_s = R_d_s(1:N); % remove extra values to have N samples
omega_d = random('Normal', omega_rated_GE, mu_omega, N, 1);
V0_d = random('Normal', V0_rated, mu_V0_rated, N, 1);
lambda_d = omega_d.*R_d_s./V0_d; % tip speed ratio at rated conditions
cp_d = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_d, theta_d);

%   __  __             _          ____           _        
%  |  \/  | ___  _ __ | |_ ___   / ___|__ _ _ __| | ___   
%  | |\/| |/ _ \| '_ \| __/ _ \ | |   / _` | '__| |/ _ \  
%  | |  | | (_) | | | | ||  __/ | |__| (_| | |  | | (_) | 
%  |_|  |_|\___/|_| |_|\__\___|  \____\__,_|_|  |_|\___/  
                                                        
                                                    
K = 0.5*cp_d.*rho_d.*pi.*R_d_s.^2.*(V0_d./omega_d).^3;
p1 = 1; % lower perentile
p2 = 99; % higher perentile
q_10 = quantile(K, p1/100);
q_90 = quantile(K, p2/100);

fig = figure('Color', 'w'); grid on; box on; hold on;
histogram(K, 'Normalization', 'pdf', 'BinMethod', 'sturges', 'DisplayName', 'Disribution')
xlabel('K [$Nms^2$]')
ylabel('pdf [-]')
xline(generator.K_opt_GE, 'r--', 'LineWidth', line_width, 'DisplayName', 'Nominal $K_{opt,GE}$')
x1 = xline(q_10, 'b--', {num2str(q_10, '%.2e')}, 'LineWidth', line_width, 'DisplayName', [num2str(p1), ' percentile'], 'FontSize', font_size);
x1.LabelHorizontalAlignment = 'left';
xline(q_90, 'b-.', {num2str(q_90, '%.2e')}, 'LineWidth', line_width, 'DisplayName',  [num2str(p2), ' percentile'], 'FontSize', font_size)
legend()
set(gca, 'FontSize', font_size)
title('$K_{opt,GE}$ distribution')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'K_GE_distribution.eps'), path_images);
end

%   ____        _ _     _        _ _           _   _             
%  |  _ \    __| (_)___| |_ _ __(_) |__  _   _| |_(_) ___  _ __  
%  | |_) |  / _` | / __| __| '__| | '_ \| | | | __| |/ _ \| '_ \ 
%  |  _ <  | (_| | \__ \ |_| |  | | |_) | |_| | |_| | (_) | | | |
%  |_| \_\  \__,_|_|___/\__|_|  |_|_.__/ \__,_|\__|_|\___/|_| |_|
                            

p3 = 1; % lower perentile
p4 = 99; % higher perentile
q_10_d = quantile(R_d, p3/100);
q_90_d = quantile(R_d, p4/100);
q_10_d_s = quantile(R_d_s, p3/100);
q_90_d_s = quantile(R_d_s, p4/100);
fig = figure('Color', 'w');hold on;grid on;box on;
plot(NaN, NaN, 'k--', 'DisplayName',  [num2str(p3), ' percentile'], 'LineWidth', line_width)
plot(NaN, NaN, 'k-.', 'DisplayName',  [num2str(p4), ' percentile'], 'LineWidth', line_width)
histogram(R_d, 'Normalization', 'pdf', 'BinMethod', 'sturges', 'FaceColor', color(1), 'DisplayName', 'Undeform Dist.')
histogram(R_d_s, 'Normalization', 'pdf', 'BinMethod', 'sturges', 'FaceColor', color(2), 'DisplayName', 'deform Dist.')
xline([q_10_d], 'b--', 'LineWidth', line_width, 'FontSize', font_size, 'Handlevisibility', 'off')
xline([q_90_d], 'b-.', 'LineWidth', line_width, 'FontSize', font_size, 'Handlevisibility', 'off')
xline([q_10_d_s], 'r--', 'LineWidth', line_width, 'FontSize', font_size, 'HandleVisibility', 'off')
xline([q_90_d_s], 'r-.', 'LineWidth', line_width, 'FontSize', font_size, 'HandleVisibility', 'off')
xline(rotor.R, '-', 'Color', color(3), 'FontSize', font_size, 'LineWidth', line_width, 'DisplayName', 'Undeform mean')
xlabel('R [m]')
ylabel('pdf [-]')
legend()
title('R distribution')
set(gca, 'FontSize', font_size)
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'R_distribution.eps'), path_images);
end
