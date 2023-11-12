clear 
close all
clc

parameters;
Rs = generator.Rs; % [Ohm] stator resistance
A = rotor.A; % [m^2] rotor swept area
B = B_eq; % [kgm^2] transmission equiuvalent damping
p = generator.p; % [-] number of pole pairs
Lambda = generator.Lambda; % [Wb] generator flux linkage
eta = generator.eta; % [-] generator efficiency
R = rotor.R; % [m]
                                                                          
physical(1) = generator.Rs; % [Ohm] stator resistance
physical(2) = rotor.A; % [m^2] rotor swept area
physical(3) = B_eq; % [kgm^2] transmission equiuvalent damping
physical(4) = generator.p; % [-] number of pole pairs
physical(5) = generator.Lambda; % [Wb] generator flux linkage
physical(6) = generator.eta; % [-] generator efficiency
physical(7) = rotor.R; % [m]
physical(8) = rho; % [kg/m^3]

% parameters without the damping
physical_no_B = physical;
physical_no_B(3) = 0;

V0_b = [4:0.05:V0_rated+0.001]'; % [m/s] wind speed 
V0_a = [V0_rated:0.5:25]';         % [m/s] wind speed above rated one  
V0_v = unique(sort([V0_b; V0_a])); % [m/s] windspeed
omega_rotor = lambda_opt.*V0_b/R; % [rad/s] omega based on the maximization of the rotor power

% Maximize the power output by selecting the proper combination of TSR and pitch angle
% The variable x containts x = (TSR, pitch angle)
lb = [7.5, -5*pi/180]; % variable lower bound
ub = [10, 5*pi/180]; % variable upper bound
P = zeros(length(V0_b),1);
min_v = zeros(length(V0_b),2); % minimization values (lambda, theta)
x0 = [7.5; 0]; % initial guess value

[lambda_GE_mean, lambda_GE, theta_opt, theta, feather, omega_rated, omega_GE, cp_GE_mean, cp_GE, P_R, P_G, P_GE, P_electro, P_joule]= fun_compute_P(physical, lambda_vector, pitch_vector, lookup_cP, lookup_Pitch, V0_b, V0_a, V0_v, V0_rated, x0, lb, ub, line_width);
[lambda_GE_mean_no_B, lambda_GE_no_B, theta_opt_no_B, theta_no_B, feather_no_B, omega_rated_no_B, omega_GE_no_B, cp_GE_mean_no_B, cp_GE_no_B, P_R_no_B, P_G_no_B, P_GE_no_B, P_electro_no_B, P_joule_no_B]= fun_compute_P(physical_no_B, lambda_vector, pitch_vector, lookup_cP, lookup_Pitch, V0_b, V0_a, V0_v, V0_rated, x0, lb, ub, line_width);

%% Plots
% Below rated
fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, omega_GE, 'LineWidth',line_width, 'DisplayName','Gen.')
plot(V0_b, omega_GE_no_B, 'LineWidth',line_width, 'DisplayName','Gen. B=0 $[\frac{kgm^2}{s}]$')
plot(V0_b, omega_rotor,'--', 'LineWidth',line_width, 'DisplayName','Rotor', 'Color', colors_vect(5,:))
title('Rotor rotational speed')
xlabel('$V_0$ [m/s]')
ylabel('$\omega$ [rad/s]')
legend('location', 'northwest')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_omega_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, theta*180/pi, 'LineWidth',line_width, 'DisplayName','Gen.')
plot(V0_b, theta_no_B*180/pi, 'LineWidth',line_width, 'DisplayName','Gen. B=0$[\frac{kgm^2}{s}]$')
yline(0, '--', 'Color', colors_vect(5,:),'LineWidth',line_width, 'FontSize', font_size, 'DisplayName','Rotor')
title('Pitch angle')
xlabel('$V_0$ [m/s]')
ylabel('$\theta$ [deg]')
legend('location','east')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_pitch_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, lambda_GE, 'LineWidth',line_width, 'DisplayName','Gen.','Color',colors_vect(1,:))
plot(V0_b, lambda_GE_no_B, 'LineWidth',line_width, 'DisplayName','Gen. B=0$[\frac{kgm^2}{s}]$','Color',colors_vect(2,:))
yl_1 = yline(lambda_GE_mean, '--', 'Mean', 'LineWidth',line_width, 'HandleVisibility', 'off', 'Color',colors_vect(1,:), 'FontSize', font_size);
yl_1.LabelHorizontalAlignment = 'left';
yl_2 = yline(lambda_GE_mean_no_B, '--', 'Mean', 'LineWidth',line_width, 'HandleVisibility', 'off', 'Color',colors_vect(2,:), 'FontSize', font_size);
yl_2.LabelHorizontalAlignment = 'left';
yl_3 = yline(lambda_opt, '--', 'Rotor','LineWidth',line_width, 'HandleVisibility', 'off','Color',colors_vect(5,:), 'FontSize', font_size);
yl_3.LabelHorizontalAlignment = 'left';
title('Tip speed ratio')
xlabel('$V_0$ [m/s]')
ylabel('$\lambda$ [-]')
legend('location', 'southeast')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_lambda_GE', '.eps'), path_images); 
end

fig = figure('Color', 'w');
hold on; box on; grid on;
plot(V0_b, cp_GE, 'LineWidth',line_width, 'DisplayName','Gen.')
plot(V0_b, cp_GE_no_B, 'LineWidth',line_width, 'DisplayName','Gen. B=0$[\frac{kgm^2}{s}]$')
y1 = yline(cp_GE_mean, '--', 'Mean', 'LineWidth',line_width, 'HandleVisibility','off','Color',colors_vect(1,:), 'FontSize', font_size);
y2 = yline(cp_GE_mean_no_B, '--', 'Mean', 'LineWidth',line_width,'HandleVisibility','off','Color',colors_vect(2,:), 'FontSize', font_size);
y1.LabelHorizontalAlignment = 'left';
y1.LabelVerticalAlignment = 'top';
y2.LabelHorizontalAlignment = 'left';
y2.LabelVerticalAlignment = 'bottom';
yline(cp_max, '--', 'LineWidth',line_width, 'DisplayName','Rotor','Color',colors_vect(5,:))
title('Power coefficient')
xlabel('$V_0$ [m/s]')
ylabel('$c_P$ [-]')
legend('location', 'east')
if simulation.print_figure == 1
  export_figure(fig, strcat(date_fig, 'fig_cP_GE', '.eps'), path_images); 
end


% Pitch angle as function of V0
fig_new_pitch_map = figure('Color','w');
hold on
plot(V0_v, feather*180/pi, 'LineWidth', line_width, 'DisplayName', 'Limit on the generator' )
plot(lookup_Pitch(1,:), lookup_Pitch(3,:)*180/pi, 'LineWidth', line_width, 'DisplayName','Limit on the rotor')
xlabel('$V_{0}$ [m/s]')
ylabel('Pitch angle [deg]')
title('Comparison of the angle for the pitching to feather')
legend('location', 'northwest')
grid on
box on
if simulation.print_figure == 1
  export_figure(fig_new_pitch_map, '\fig_new_pitch_map.eps', path_images);
end


%   ____                    _             _                
%  / ___|  __ ___   _____  | | ___   ___ | | ___   _ _ __  
%  \___ \ / _` \ \ / / _ \ | |/ _ \ / _ \| |/ / | | | '_ \ 
%   ___) | (_| |\ V /  __/ | | (_) | (_) |   <| |_| | |_) |
%  |____/ \__,_| \_/ \___| |_|\___/ \___/|_|\_\\__,_| .__/ 
%                                                   |_|    
clear("lookup_P_GE")
clear("rated_values_P_GE")
clear("lookup_pitch_P_GE")
 
rated_values_P_GE(1) = lambda_GE_mean;
rated_values_P_GE(2) = cp_GE_mean;
rated_values_P_GE(3) = omega_rated;
rated_values_P_GE(4) = theta_opt;
save('lookup\rated_values_P_GE.mat', "rated_values_P_GE");

rated_values_P_GE_no_B(1) = lambda_GE_mean_no_B;
rated_values_P_GE_no_B(2) = cp_GE_mean_no_B;
rated_values_P_GE_no_B(3) = omega_rated_no_B;
rated_values_P_GE_no_B(4) = theta_opt_no_B;
save('lookup\rated_values_P_GE_no_B.mat', "rated_values_P_GE_no_B");
 
% pitch to feather
lookup_pitch_P_GE = zeros(2, length(V0_v));
lookup_pitch_P_GE(1, :) = V0_v;
lookup_pitch_P_GE(2, :) = feather;
save('lookup\lookup_pitch_P_GE.mat', "lookup_pitch_P_GE");

% generator electrical power
lookup_P_GE(1, :) = V0_v;
lookup_P_GE(2, :) = P_GE;
lookup_P_GE(3, :) = P_GE_no_B;
save('lookup\lookup_P_GE.mat', "lookup_P_GE");


%    __                                                _           ____  
%   / _|_   _ _ __      ___ ___  _ __ ___  _ __  _   _| |_ ___    |  _ \ 
%  | |_| | | | '_ \    / __/ _ \| '_ ` _ \| '_ \| | | | __/ _ \   | |_) |
%  |  _| |_| | | | |  | (_| (_) | | | | | | |_) | |_| | ||  __/   |  __/ 
%  |_|  \__,_|_| |_|___\___\___/|_| |_| |_| .__/ \__,_|\__\___|___|_|    
%                 |_____|                 |_|                |_____|     
function [lambda_opt, lambda, theta_opt, theta, feather, omega_rated, omega, cp_mean, cp, P_R, P_G, P_GE, P_electro, P_joule] = fun_compute_P(physical, lambda_vector, pitch_vector, lookup_cP, lookup_Pitch, V0_b, V0_a, V0_v, V0_rated, x0, lb, ub, line_width)
  %   ____                                _                
  %  |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
  %  | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
  %  |  __/ (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
  %  |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
                                                        
  Rs = physical(1);  % [Ohm] stator resistance
  A = physical(2);  % [m^2] rotor swept area
  B_eq  = physical(3);  % [kgm^2] transmission equiuvalent damping
  p = physical(4);  % [-] number of pole pairs
  Lambda = physical(5);  % [Wb] generator flux linkage
  eta = physical(6);  % [-] generator efficiency
  R = physical(7);  % [m]
  rho  = physical(8);  % [kg/m^3]
  len_V0b = length(V0_b);
  len_V0a = length(V0_a);

  %   ____       _                           _           _ 
  %  | __ )  ___| | _____      __  _ __ __ _| |_ ___  __| |
  %  |  _ \ / _ \ |/ _ \ \ /\ / / | '__/ _` | __/ _ \/ _` |
  %  | |_) |  __/ | (_) \ V  V /  | | | (_| | ||  __/ (_| |
  %  |____/ \___|_|\___/ \_/\_/   |_|  \__,_|\__\___|\__,_|
                                                    
  options = optimoptions(@fmincon,'Display','off');
  for i=1:length(V0_b)
    [min_v(i, :), P(i, :)] = fmincon(@(x)compute_P_GE(x, physical, lambda_vector, pitch_vector, lookup_cP, V0_b(i)), x0, [], [], [], [], lb, ub, [], options);
  end
  lambda = min_v(:, 1);
  omega = lambda.*V0_b/R; % [rad/s] rotational speed
  theta = min_v(:, 2); % [rad] pitch angle

  % scalar values
  lambda_opt = mean(min_v(:, 1)); % optimal TSR based on generator maximization
  theta_opt = mean(min_v(:, 2)); % pitch coming from minimization
  omega_rated = lambda_opt*V0_rated/R; % [rad/s] rated velocity
  
  % vector values
  P_GE_rated = -P(end); % [W] rated generator electrical power 
  P_GE_b = -P;
  
  %      _    _                                _           _ 
  %     / \  | |__   _____   _____   _ __ __ _| |_ ___  __| |
  %    / _ \ | '_ \ / _ \ \ / / _ \ | '__/ _` | __/ _ \/ _` |
  %   / ___ \| |_) | (_) \ V /  __/ | | | (_| | ||  __/ (_| |
  %  /_/   \_\_.__/ \___/ \_/ \___| |_|  \__,_|\__\___|\__,_|
                                                          
  theta_p = [-1:0.1:25]*pi/180;       % [rad] range of angle where to search for the optimum 
  lambda_a = zeros(length(V0_a), 1);  % [-] TSR
  cp_a = zeros(length(theta_p), length(V0_a));
  Pg_out_a = zeros(length(theta_p), length(V0_a)); % [W] power above the rated
  P_GE_a = zeros(length(V0_a), 1);

  for i=1:length(V0_a)
    lambda_a = omega_rated*R/V0_a(i); % TSR
    cp_a(:, i) = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_a, theta_p','cubic',0);
    x = [lambda_a*ones(1,length(theta_p)); theta_p];
    Pg_out_a(:, i) = -compute_P_GE(x, physical, lambda_vector, pitch_vector, lookup_cP, V0_a(i));
  end

  [~, max_p_a] = min(abs(Pg_out_a - P_GE_rated), [], 1); % extract the maximum of the power
  for i=1:length(V0_a)
    P_GE_a(i) = Pg_out_a(max_p_a(i), i); % [W] maximum power 
  end
  feather_a = theta_p(max_p_a)'; % [rad] angle for feather pitching 

  %    ___                      _ _ 
  %   / _ \__   _____ _ __ __ _| | |
  %  | | | \ \ / / _ \ '__/ _` | | |
  %  | |_| |\ V /  __/ | | (_| | | |
  %   \___/  \_/ \___|_|  \__,_|_|_|
                                
  lambda_v = [lambda_opt*ones(len_V0b,1); omega_rated*R./V0_a]; % [-] TSR
  omega_v = [lambda_opt.*V0_b/R; omega_rated*ones(len_V0a,1)]; % [rad/s] rotational speed
  feather = [theta; feather_a]; % [rad]
  cP_R = interp2(lambda_vector, pitch_vector, lookup_cP, lambda_v, feather,'cubic',0); % [-] feather power coefficient
  P_R = 0.5*A*rho.*cP_R.*V0_v.^3; % [W] power to the rotor side
  cp = P_R(1:len_V0b)./(0.5*A*rho.*V0_b.^3); % [-] power coefficient
  cp_mean = mean(cp); % mean power coefficient
  P_G = P_R - B_eq*omega_v.^2; % [W] mechancial power to the generator
  P_GE = [P_GE_b; P_GE_a];
  iq = 2/3*P_G./(Lambda.*omega_v*p);  % [A] generator current
  uq =  omega_v*p*Lambda - Rs*iq;     % [V] generator voltage
  P_electro = 1.5*uq.*iq;           % [W] electrical power
  P_joule = 1.5*Rs*iq.^2;           % [W] Joule losses

end